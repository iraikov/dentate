
import os, sys, click
import dentate
from dentate import env, plot, utils, cells, neuron_utils
from dentate.neuron_utils import h, configure_hoc_env
from dentate.env import Env
from dentate.cells import make_biophys_cell, init_biophysics, report_topology
from mpi4py import MPI

script_name = os.path.basename(__file__)

@click.command()
@click.option("--config-file", '-c', required=True, type=str, help='model configuration file name')
@click.option("--population", '-p', required=True, type=str, help='target population')
@click.option("--gid", '-g', required=True, type=int, help='target cell gid')
@click.option("--template-paths", type=str, default="templates",
              help='colon-separated list of paths to directories containing hoc cell templates')
@click.option("--dataset-prefix", type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='path to directory containing required neuroh5 data files')
@click.option("--config-prefix", type=click.Path(exists=True, file_okay=False, dir_okay=True),
              default='config', help='path to directory containing network and cell mechanism config files')
@click.option("--data-file", required=False, type=click.Path(exists=True, file_okay=True, dir_okay=False),
              help='input data file (overrides file provided in configuration)')
@click.option("--load-synapses", "-s", type=bool, default=False, is_flag=True)
@click.option("--syn-sources", multiple=True, type=str, help='synapse filter for source populations')
@click.option("--syn-source-threshold", type=float, help='only show synapses for sources in top n percentile')
@click.option("--syn-types", multiple=True, type=str, help='synapse filter for synapse types')
@click.option("--font-size", type=float, default=14)
@click.option("--bgcolor", type=(float,float,float), default=(0.,0.,0.))
@click.option("--colormap", type=str, default='coolwarm')
@click.option("--verbose", "-v", type=bool, default=False, is_flag=True)
def main(config_file, population, gid, template_paths, dataset_prefix, config_prefix, data_file, load_synapses, syn_types, syn_sources, syn_source_threshold, font_size, bgcolor, colormap, verbose):

    utils.config_logging(verbose)
    logger = utils.get_script_logger(script_name)

    if dataset_prefix is None and data_file is None:
        raise RuntimeError('Either --dataset-prefix or --data-file must be provided.')

    params = dict(locals())
    env = Env(**params)
    configure_hoc_env(env)

    if env.data_file_path is None:
        env.data_file_path = data_file
        env.data_file_path = data_file
        env.connectivity_file_path = data_file
        env.load_celltypes()

    ## Determine if a mechanism configuration file exists for this cell type
    if 'mech_file_path' in env.celltypes[population]:
        mech_file_path = env.celltypes[population]['mech_file_path']
    else:
        mech_file_path = None

    ## Determine if correct_for_spines flag has been specified for this cell type
    synapse_config = env.celltypes[population]['synapses']
    if 'correct_for_spines' in synapse_config:
        correct_for_spines_flag = synapse_config['correct_for_spines']
    else:
        correct_for_spines_flag = False

    logger.info('loading cell %i' % gid)

    load_weights = False
    biophys_cell = make_biophys_cell(env, population, gid, 
                                     load_synapses=load_synapses,
                                     load_weights=load_weights, 
                                     load_edges=load_synapses,
                                     mech_file_path=mech_file_path)
    
    init_biophysics(biophys_cell, reset_cable=True, 
                    correct_cm=correct_for_spines_flag,
                    correct_g_pas=correct_for_spines_flag, env=env)
    report_topology(biophys_cell, env)
    
    
    if len(syn_types) == 0:
        syn_types = None
    else:
        syn_types = list(syn_types)
    if len(syn_sources) == 0:
        syn_sources = None
    else:
        syn_sources = list(syn_sources)
        
    plot.plot_biophys_cell_tree (env, biophys_cell, saveFig=True,
                                 syn_source_threshold=syn_source_threshold,
                                 synapse_filters={'syn_types': syn_types, 'sources': syn_sources},
                                 bgcolor=bgcolor, colormap=colormap)
    

if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == script_name, sys.argv)+1):])
