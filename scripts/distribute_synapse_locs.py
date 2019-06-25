import os, sys, gc, logging, string, time, itertools
from mpi4py import MPI
import click
from collections import defaultdict
import numpy as np
import dentate
from dentate import cells, neuron_utils, synapses, utils
from dentate.env import Env
from dentate.neuron_utils import configure_hoc_env
from dentate.utils import *
from neuroh5.io import NeuroH5TreeGen, append_cell_attributes, read_population_ranges
import h5py

#sys_excepthook = sys.excepthook
#def mpi_excepthook(type, value, traceback):
#    sys_excepthook(type, value, traceback)
#    if MPI.COMM_WORLD.size > 1:
#        MPI.COMM_WORLD.Abort(1)
#sys.excepthook = mpi_excepthook


def update_syn_stats(env, syn_stats_dict, syn_dict):

    syn_type_excitatory = env.Synapse_Types['excitatory']
    syn_type_inhibitory = env.Synapse_Types['inhibitory']

    this_syn_stats_dict = { 'section': defaultdict(lambda: { 'excitatory': 0, 'inhibitory': 0 }), \
                            'layer': defaultdict(lambda: { 'excitatory': 0, 'inhibitory': 0 }), \
                            'swc_type': defaultdict(lambda: { 'excitatory': 0, 'inhibitory': 0 }), \
                            'total': { 'excitatory': 0, 'inhibitory': 0 } }

    for (syn_id,syn_sec,syn_type,swc_type,syn_layer) in \
        zip(syn_dict['syn_ids'],
                       syn_dict['syn_secs'],
                       syn_dict['syn_types'],
                       syn_dict['swc_types'],
                       syn_dict['syn_layers']):
        
        if syn_type == syn_type_excitatory:
            syn_type_str = 'excitatory'
        elif syn_type == syn_type_inhibitory:
            syn_type_str = 'inhibitory'
        else:
            raise ValueError('Unknown synapse type %s' % str(syn_type))

        syn_stats_dict['section'][syn_sec][syn_type_str] += 1
        syn_stats_dict['layer'][syn_layer][syn_type_str] += 1
        syn_stats_dict['swc_type'][swc_type][syn_type_str] += 1
        syn_stats_dict['total'][syn_type_str] += 1

        this_syn_stats_dict['section'][syn_sec][syn_type_str] += 1
        this_syn_stats_dict['layer'][syn_layer][syn_type_str] += 1
        this_syn_stats_dict['swc_type'][swc_type][syn_type_str] += 1
        this_syn_stats_dict['total'][syn_type_str] += 1

    return this_syn_stats_dict


def global_syn_summary(comm, syn_stats, global_count, root):
    res = []
    for population in syn_stats:
        pop_syn_stats = syn_stats[population]
        for part in ['layer', 'swc_type']:
            syn_stats_dict = pop_syn_stats[part]
            for part_name in syn_stats_dict:
                for syn_type in syn_stats_dict[part_name]:
                    global_syn_count = comm.gather(syn_stats_dict[part_name][syn_type], root=root)
                    if comm.rank == root:
                        res.append("%s %s %s: mean %s synapses per cell: %f" % (population, part, part_name, syn_type, old_div(np.sum(global_syn_count),global_count)))
        total_syn_stats_dict = pop_syn_stats['total']
        for syn_type in total_syn_stats_dict:
            global_syn_count = comm.gather(total_syn_stats_dict[syn_type], root=root)
            if comm.rank == root:
                res.append("%s: mean %s synapses per cell: %f" % (population, syn_type, old_div(np.sum(global_syn_count),global_count)))
        
    return str.join('\n', res)

def local_syn_summary(syn_stats_dict):
    res = []
    for part_name in ['layer','swc_type']:
        for part_type in syn_stats_dict[part_name]:
            syn_count_dict = syn_stats_dict[part_name][part_type]
            for syn_type, syn_count in list(syn_count_dict.items()):
                res.append("%s %i: %s synapses: %i" % (part_name, part_type, syn_type, syn_count))
    return str.join('\n', res)


def check_syns(gid, morph_dict, syn_stats_dict, seg_density_per_sec, layer_set_dict, swc_set_dict, env, logger):

    layer_stats = syn_stats_dict['layer']
    swc_stats = syn_stats_dict['swc_type']

    warning_flag = False
    for syn_type, layer_set in list(layer_set_dict.items()):
        for layer in layer_set:
            if layer in layer_stats:
                if layer_stats[layer][syn_type] <= 0:
                    warning_flag = True
            else:
                warning_flag = True
    if warning_flag:
        logger.warning('Rank %d: incomplete synapse layer set for cell %d: %s' % (env.comm.Get_rank(), gid, str(layer_stats)))
        logger.info('layer_set_dict: %s' % str(layer_set_dict))
        logger.info('gid %d: seg_density_per_sec: %s' % (gid, str(seg_density_per_sec)))
        logger.info('gid %d: morph_dict: %s' % (gid, str(morph_dict)))
    for syn_type, swc_set in viewitems(swc_set_dict):
        for swc_type in swc_set:
            if swc_type in swc_stats:
                if swc_stats[swc_type][syn_type] <= 0:
                    warning_flag = True
            else:
                warning_flag = True
    if warning_flag:
        logger.warning('Rank %d: incomplete synapse swc type set for cell %d: %s' % (env.comm.Get_rank(), gid, str(swc_stats)))
        logger.info('swc_set_dict: %s' % str(swc_set_dict.items))
        logger.info('gid %d: seg_density_per_sec: %s' % (gid, str(seg_density_per_sec)))
        logger.info('gid %d: morph_dict: %s' % (gid, str(morph_dict)))
                

            
        
@click.command()
@click.option("--config", required=True, type=str)
@click.option("--config-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), default='config')
@click.option("--template-path", type=str)
@click.option("--output-path", type=click.Path(exists=False, file_okay=True, dir_okay=False))
@click.option("--forest-path", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--populations", '-i', required=True, multiple=True, type=str)
@click.option("--distribution", type=str, default='uniform')
@click.option("--io-size", type=int, default=-1)
@click.option("--chunk-size", type=int, default=1000)
@click.option("--value-chunk-size", type=int, default=1000)
@click.option("--cache-size", type=int, default=10000)
@click.option("--verbose", "-v", is_flag=True)
@click.option("--dry-run", is_flag=True)
def main(config, config_prefix, template_path, output_path, forest_path, populations, distribution, io_size, chunk_size, value_chunk_size,
         cache_size, verbose, dry_run):
    """

    :param config:
    :param config_prefix:
    :param template_path:
    :param forest_path:
    :param populations:
    :param distribution:
    :param io_size:
    :param chunk_size:
    :param value_chunk_size:
    :param cache_size:
    """

    utils.config_logging(verbose)
    logger = utils.get_script_logger(os.path.basename(__file__))
        
    comm = MPI.COMM_WORLD
    rank = comm.rank
    
    if rank == 0:
        logger.info('%i ranks have been allocated' % comm.size)

    comm.barrier()
    env = Env(comm=MPI.COMM_WORLD, config_file=config, config_prefix=config_prefix, template_paths=template_path)

    configure_hoc_env(env)
    
    if io_size == -1:
        io_size = comm.size

    if output_path is None:
        output_path = forest_path

    if not dry_run:
        if rank==0:
            if not os.path.isfile(output_path):
                input_file  = h5py.File(forest_path,'r')
                output_file = h5py.File(output_path,'w')
                input_file.copy('/H5Types',output_file)
                input_file.close()
                output_file.close()
        comm.barrier()
        
    (pop_ranges, _) = read_population_ranges(forest_path, comm=comm)
    start_time = time.time()
    syn_stats = {}
    for population in populations:
        logger.info('Rank %i population: %s' % (rank, population))
        (population_start, _) = pop_ranges[population]
        template_class = env.load_cell_template(population)
        
        density_dict = env.celltypes[population]['synapses']['density']
        layer_set_dict = defaultdict(set)
        swc_set_dict = defaultdict(set)
        for sec_name, sec_dict in viewitems(density_dict):
            for syn_type, syn_dict in viewitems(sec_dict):
                swc_set_dict[syn_type].add(env.SWC_Types[sec_name])
                for layer_name in syn_dict:
                    if layer_name != 'default':
                        layer = env.layers[layer_name]
                        layer_set_dict[syn_type].add(layer)
        
        syn_stats_dict = { 'section': defaultdict(lambda: { 'excitatory': 0, 'inhibitory': 0 }), \
                           'layer': defaultdict(lambda: { 'excitatory': 0, 'inhibitory': 0 }), \
                           'swc_type': defaultdict(lambda: { 'excitatory': 0, 'inhibitory': 0 }), \
                           'total': { 'excitatory': 0, 'inhibitory': 0 } }

        count = 0
        for gid, morph_dict in NeuroH5TreeGen(forest_path, population, io_size=io_size, comm=comm, topology=True):
            local_time = time.time()
            synapse_dict = {}
            if gid is not None:
                logger.info('Rank %i gid: %i' % (rank, gid))
                cell = cells.make_neurotree_cell(template_class, neurotree_dict=morph_dict, gid=gid)
                cell_sec_dict = {'apical': (cell.apical, None), 'basal': (cell.basal, None), 'soma': (cell.soma, None), 'ais': (cell.ais, None)}
                cell_secidx_dict = {'apical': cell.apicalidx, 'basal': cell.basalidx, 'soma': cell.somaidx, 'ais': cell.aisidx}

                random_seed = env.modelConfig['Random Seeds']['Synapse Locations'] + gid
                if distribution == 'uniform':
                    syn_dict, seg_density_per_sec = synapses.distribute_uniform_synapses(random_seed, env.Synapse_Types, env.SWC_Types, env.layers,
                                                                                         density_dict, morph_dict,
                                                                                         cell_sec_dict, cell_secidx_dict)
                                                                    
                    
                elif distribution == 'poisson':
                    syn_dict, seg_density_per_sec = synapses.distribute_poisson_synapses(random_seed, env.Synapse_Types, env.SWC_Types, env.layers,
                                                                                         density_dict, morph_dict,
                                                                                         cell_sec_dict, cell_secidx_dict)
                else:
                    raise Exception('Unknown distribution type: %s' % distribution)

                synapse_dict[gid] = syn_dict
                this_syn_stats = update_syn_stats (env, syn_stats_dict, syn_dict)
                check_syns(gid, morph_dict, this_syn_stats, seg_density_per_sec, layer_set_dict, swc_set_dict, env, logger)
                
                del cell
                num_syns = len(synapse_dict[gid]['syn_ids'])
                logger.info('Rank %i took %i s to compute %d synapse locations for %s gid: %i' % (rank, time.time() - local_time, num_syns, population, gid))
                logger.info('%s gid %i synapses: %s' % (population, gid, local_syn_summary(this_syn_stats)))
                count += 1
            else:
                logger.info('Rank %i gid is None' % rank)
            if not dry_run:
                append_cell_attributes(output_path, population, synapse_dict,
                                       namespace='Synapse Attributes', comm=comm, io_size=io_size, chunk_size=chunk_size,
                                       value_chunk_size=value_chunk_size, cache_size=cache_size)
            syn_stats[population] = syn_stats_dict
            del synapse_dict
            gc.collect()

        global_count = comm.gather(count, root=0)

        if count > 0:
            color = 1
        else:
            color = 0
            
        comm0 = comm.Split(color, 0)
        if color == 1:
            summary = global_syn_summary(comm0, syn_stats, np.sum(global_count), root=0)
            if rank == 0:
                logger.info('target: %s, %i ranks took %i s to compute synapse locations for %i cells' % (population, comm.size,time.time() - start_time,np.sum(global_count)))
                logger.info(summary)
        comm0.Free()
        comm.barrier()
            
    MPI.Finalize()


if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == os.path.basename(__file__), sys.argv)+1):])
