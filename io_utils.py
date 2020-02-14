import os, itertools, pprint
from collections import defaultdict
from mpi4py import MPI
import h5py
import numpy as np
import dentate
from dentate.utils import Struct, range, str, viewitems, basestring, Iterable, compose_iter, get_module_logger
from neuroh5.io import write_cell_attributes, append_cell_attributes, append_cell_trees, write_graph, read_cell_attribute_selection, read_tree_selection, read_graph_selection
from neuron import h


# This logger will inherit its settings from the root logger, created in dentate.env
logger = get_module_logger(__name__)


grp_h5types = 'H5Types'
grp_projections = 'Projections'
grp_populations = 'Populations'

path_population_labels = '/%s/Population labels' % grp_h5types
path_population_range = '/%s/Population range' % grp_h5types

grp_population_projections = 'Population projections'
grp_valid_population_projections = 'Valid population projections'
path_population_projections = '/%s/Population projections' % grp_h5types

# Default I/O configuration
default_io_options = Struct(io_size=-1, chunk_size=1000, value_chunk_size=1000, cache_size=50, write_size=10000)


def list_concat(a, b, datatype):
    return a+b

mpi_op_concat = MPI.Op.Create(list_concat, commute=True)


def h5_get_group(h, groupname):
    if groupname in h:
        g = h[groupname]
    else:
        g = h.create_group(groupname)
    return g


def h5_get_dataset(g, dsetname, **kwargs):
    if dsetname in g:
        dset = g[dsetname]
    else:
        dset = g.create_dataset(dsetname, (0,), **kwargs)
    return dset


def h5_concat_dataset(dset, data):
    dsize = dset.shape[0]
    newshape = (dsize + len(data),)
    dset.resize(newshape)
    dset[dsize:] = data
    return dset


def make_h5types(env, output_path, gap_junctions=False):
    populations = []
    for pop_name, pop_idx in viewitems(env.Populations):
        layer_counts = env.geometry['Cell Layer Counts'][pop_name]
        pop_count = 0
        for layer_name, layer_count in viewitems(layer_counts):
            pop_count += layer_count
        populations.append((pop_name, pop_idx, pop_count))
    populations.sort(key=lambda x: x[1])

    projections = []
    if gap_junctions:
        for (post, pre), connection_dict in viewitems(env.gapjunctions):
            projections.append((env.Populations[pre], env.Populations[post]))
    else:
        for post, connection_dict in viewitems(env.connection_config):
            for pre, _ in viewitems(connection_dict):
                projections.append((env.Populations[pre], env.Populations[post]))

    # create an HDF5 enumerated type for the population label
    mapping = {name: idx for name, idx in viewitems(env.Populations)}
    dt_population_labels = h5py.special_dtype(enum=(np.uint16, mapping))

    with h5py.File(output_path, "a") as h5:

        h5[path_population_labels] = dt_population_labels

        dt_populations = np.dtype([("Start", np.uint64), ("Count", np.uint32),
                                   ("Population", h5[path_population_labels].dtype)])
        h5[path_population_range] = dt_populations

        # create an HDF5 compound type for population ranges
        dt = h5[path_population_range].dtype

        g = h5_get_group(h5, grp_h5types)

        dset = h5_get_dataset(g, grp_populations, maxshape=(len(populations),), dtype=dt)
        dset.resize((len(populations),))
        a = np.zeros(len(populations), dtype=dt)

        start = 0
        for name, idx, count in populations:
            a[idx]["Start"] = start
            a[idx]["Count"] = count
            a[idx]["Population"] = idx
            start += count

        dset[:] = a

        dt_projections = np.dtype([("Source", h5[path_population_labels].dtype),
                                   ("Destination", h5[path_population_labels].dtype)])

        h5[path_population_projections] = dt_projections

        dt = h5[path_population_projections]
        dset = h5_get_dataset(g, grp_valid_population_projections,
                              maxshape=(len(projections),), dtype=dt)
        dset.resize((len(projections),))
        a = np.zeros(len(projections), dtype=dt)
        idx = 0
        for i, prj in enumerate(projections):
            src, dst = prj
            a[i]["Source"] = int(src)
            a[i]["Destination"] = int(dst)

        dset[:] = a


def mkout(env, results_filename):
    """
    Creates simulation results file and adds H5Types group compatible with NeuroH5.

    :param env:
    :param results_filename:
    :return:
    """
    dataset_path = os.path.join(env.dataset_prefix, env.datasetName)
    data_file_path = os.path.join(dataset_path, env.model_config['Cell Data'])
    data_file = h5py.File(data_file_path, 'r')
    results_file = h5py.File(results_filename)
    if 'H5Types' not in results_file:
        data_file.copy('/H5Types', results_file)
    data_file.close()
    results_file.close()


def spikeout(env, output_path, t_start=None, clear_data=False):
    """
    Writes spike times to specified NeuroH5 output file.

    :param env:
    :param output_path:
    :param clear_data: 
    :return:
    """

    t_vec = np.array(env.t_vec, dtype=np.float32)
    id_vec = np.array(env.id_vec, dtype=np.uint32)

    binlst = []
    typelst = sorted(env.celltypes.keys())
    binvect = np.asarray([env.celltypes[k]['start'] for k in typelst ])
    sort_idx = np.argsort(binvect, axis=0)
    pop_names = [typelst[i] for i in sort_idx]
    bins = binvect[sort_idx][1:]
    inds = np.digitize(id_vec, bins)

    if env.results_namespace_id is None:
        namespace_id = "Spike Events"
    else:
        namespace_id = "Spike Events %s" % str(env.results_namespace_id)

    equilibration_duration = float(env.stimulus_config['Equilibration Duration'])
    for i, pop_name in enumerate(pop_names):
        spkdict = {}
        sinds = np.where(inds == i)
        if len(sinds) > 0:
            ids = id_vec[sinds]
            ts = t_vec[sinds]
            for j in range(0, len(ids)):
                gid = ids[j]
                t = ts[j]
                if (t_start is None) or (t >= t_start):
                    if gid in spkdict:
                        spkdict[gid]['t'].append(t)
                    else:
                        spkdict[gid] = {'t': [t]}
            for gid in spkdict:
                spiketrain = np.array(spkdict[gid]['t'], dtype=np.float32)
                spiketrain -= equilibration_duration
                if gid in env.spike_onset_delay:
                    spiketrain -= env.spike_onset_delay[gid]
                spkdict[gid]['t'] = spiketrain
        append_cell_attributes(output_path, pop_name, spkdict, namespace=namespace_id, comm=env.comm, io_size=env.io_size)
        del (spkdict)

    if clear_data:
        env.t_vec.resize(0)
        env.id_vec.resize(0)

    env.comm.barrier()
    if env.comm.Get_rank() == 0:
        logger.info("*** Output spike results to file %s" % output_path)

def recsout(env, output_path, t_start=None, clear_data=False):
    """
    Writes intracellular state traces to specified NeuroH5 output file.

    :param env:
    :param output_path:
    :param clear_data:
    :param reduce_data:
    :return:
    """
    t_rec = env.t_rec
    equilibration_duration = float(env.stimulus_config['Equilibration Duration'])
    reduce_data = env.recording_profile.get('reduce', None)

    for pop_name in sorted(env.celltypes.keys()):
        local_rec_types = list(env.recs_dict[pop_name].keys())
        rec_types = sorted(set(env.comm.allreduce(local_rec_types, op=mpi_op_concat)))
        for rec_type in rec_types:
            recs = env.recs_dict[pop_name][rec_type]
            attr_dict = defaultdict(lambda: {})
            for rec in recs:
                gid = rec['gid']
                data_vec = np.array(rec['vec'], copy=clear_data, dtype=np.float32)
                time_vec = np.array(t_rec, copy=clear_data, dtype=np.float32)[:-1]
                time_vec -= equilibration_duration
                if t_start is not None:
                    time_inds = np.where(time_vec >= t_start)[0]
                    time_vec = time_vec[time_inds]
                    data_vec = data_vec[time_inds]
                label = rec['label']
                if label in attr_dict[gid]:
                    if reduce_data is None:
                        raise RuntimeError('recsout: duplicate recorder labels and no reduce strategy specified')
                    elif reduce_data is True:
                        attr_dict[gid][label] += data_vec
                    else:
                        raise RuntimeError('recsout: unsupported reduce strategy specified')
                else:
                    attr_dict[gid][label] = data_vec
                    attr_dict[gid]['t'] = time_vec
                distance = rec.get('distance', None)
                if distance is not None:
                    attr_dict[gid]['distance'] = np.asarray([distance], dtype=np.float32)
                section = rec.get('section', None)
                if section is not None:
                    attr_dict[gid]['section'] = np.asarray([section], dtype=np.int16)
                loc = rec.get('loc', None)
                if loc is not None:
                    attr_dict[gid]['loc'] = np.asarray([loc], dtype=np.float32)
                if clear_data:
                    rec['vec'].resize(0)
            if env.results_namespace_id is None:
                namespace_id = "Intracellular %s" % (rec_type)
            else:
                namespace_id = "Intracellular %s %s" % (rec_type, str(env.results_namespace_id))
            append_cell_attributes(output_path, pop_name, attr_dict, namespace=namespace_id,
                                   comm=env.comm, io_size=env.io_size)
    if clear_data:
        env.t_rec.resize(0)

    env.comm.barrier()
    if env.comm.Get_rank() == 0:
        logger.info("*** Output intracellular state results to file %s" % output_path)



def lfpout(env, output_path):
    """
    Writes local field potential voltage traces to specified HDF5 output file.

    :param env:
    :param output_path:
    :param clear_data:
    :return:
    """

    for lfp in list(env.lfp.values()):

        if env.results_namespace_id is None:
            namespace_id = "Local Field Potential %s" % str(lfp.label)
        else:
            namespace_id = "Local Field Potential %s %s" % (str(lfp.label), str(env.results_namespace_id))
        import h5py
        output = h5py.File(output_path, 'a')

        grp = output.create_group(namespace_id)

        grp['t'] = np.asarray(lfp.t, dtype=np.float32)
        grp['v'] = np.asarray(lfp.meanlfp, dtype=np.float32)

        output.close()

    if env.comm.Get_rank() == 0:
        logger.info("*** Output LFP results to file %s" % output_path)


def get_h5py_attr(attrs, key):
    """
    str values are stored as bytes in h5py container attrs dictionaries. This function enables py2/py3 compatibility by
    always returning them to str type upon read. Values should be converted during write with the companion function
    set_h5py_str_attr.
    :param attrs: :class:'h5py._hl.attrs.AttributeManager'
    :param key: str
    :return: val with type converted if str or array of str
    """
    if key not in attrs:
        raise KeyError('get_h5py_attr: invalid key: %s' % key)
    val = attrs[key]
    if isinstance(val, basestring):
        val = np.string_(val).astype(str)
    elif isinstance(val, Iterable) and len(val) > 0:
        if isinstance(val[0], basestring):
            val = np.array(val, dtype='str')
    return val


def set_h5py_attr(attrs, key, val):
    """
    str values are stored as bytes in h5py container attrs dictionaries. This function enables py2/py3 compatibility by
    always converting them to np.string_ upon write. Values should be converted back to str during read with the
    companion function get_h5py_str_attr.
    :param attrs: :class:'h5py._hl.attrs.AttributeManager'
    :param key: str
    :param val: type converted if str or array of str
    """
    if isinstance(val, basestring):
        val = np.string_(val)
    elif isinstance(val, Iterable) and len(val) > 0:
        if isinstance(val[0], basestring):
            val = np.array(val, dtype='S')
    attrs[key] = val


def write_cell_selection(env, write_selection_file_path, populations=None, write_kwds={}):
    """
    Writes out the data necessary to instantiate the selected cells.

    :param env: an instance of the `dentate.Env` class
    """

    if 'comm' not in write_kwds:
        write_kwds['comm'] = env.comm
    if 'io_size' not in write_kwds:
        write_kwds['io_size'] = env.io_size

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    dataset_path = env.dataset_path
    data_file_path = env.data_file_path

    if populations is None:
        pop_names = sorted(env.cell_selection.keys())
    else:
        pop_names = populations
    
    for pop_name in pop_names:

        gid_range = [gid for gid in env.cell_selection[pop_name] if gid % nhosts == rank]

        trees_output_dict = {}
        coords_output_dict = {}
        num_cells = 0
        if (pop_name in env.cell_attribute_info) and ('Trees' in env.cell_attribute_info[pop_name]):
            if rank == 0:
                logger.info("*** Reading trees for population %s" % pop_name)

            cell_tree_iter, _ = read_tree_selection(data_file_path, pop_name, selection=gid_range, \
                                                 topology=False, comm=env.comm)
            if rank == 0:
                logger.info("*** Done reading trees for population %s" % pop_name)

            for i, (gid, tree) in enumerate(cell_tree_iter):
                trees_output_dict[gid] = tree
                num_cells += 1

        elif (pop_name in env.cell_attribute_info) and ('Coordinates' in env.cell_attribute_info[pop_name]):
            if rank == 0:
                logger.info("*** Reading coordinates for population %s" % pop_name)

            cell_attributes_iter = read_cell_attribute_selection(data_file_path, pop_name, selection=gid_range, \
                                                                 namespace='Coordinates', comm=env.comm)

            if rank == 0:
                logger.info("*** Done reading coordinates for population %s" % pop_name)

            for i, (gid, coords) in enumerate(cell_attributes_iter):
                coords_output_dict[gid] = coords
                num_cells += 1

            
        if rank == 0:
            logger.info("*** Writing cell selection for population %s to file %s" % (pop_name, write_selection_file_path))
        append_cell_trees(write_selection_file_path, pop_name, trees_output_dict, create_index=True, **write_kwds)
        write_cell_attributes(write_selection_file_path, pop_name, coords_output_dict, namespace='Coordinates', **write_kwds)



def write_connection_selection(env, write_selection_file_path, populations=None, write_kwds={}):
    """
    Loads NeuroH5 connectivity file, and writes the corresponding
    synapse and network connection mechanisms for the selected postsynaptic cells.

    :param env: an instance of the `dentate.Env` class
    """
    
    if 'comm' not in write_kwds:
        write_kwds['comm'] = env.comm
    if 'io_size' not in write_kwds:
        write_kwds['io_size'] = env.io_size

    
    connectivity_file_path = env.connectivity_file_path
    forest_file_path = env.forest_file_path
    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())
    syn_attrs = env.synapse_attributes

    if populations is None:
        pop_names = sorted(env.cell_selection.keys())
    else:
        pop_names = populations

    input_sources = {pop_name: set([]) for pop_name in env.celltypes}

    for (postsyn_name, presyn_names) in sorted(viewitems(env.projection_dict)):

        if rank == 0:
            logger.info('*** Writing connection selection of population %s' % (postsyn_name))

        if postsyn_name not in pop_names:
            continue

        gid_range = [gid for gid in env.cell_selection[postsyn_name] if gid % nhosts == rank]

        synapse_config = env.celltypes[postsyn_name]['synapses']

        if 'weights' in synapse_config:
            has_weights = synapse_config['weights']
        else:
            has_weights = False

        weights_namespaces = []
        if 'weights' in synapse_config:
            has_weights = synapse_config['weights']
            if has_weights:
                if 'weights namespace' in synapse_config:
                    weights_namespaces.append(synapse_config['weights namespace'])
                elif 'weights namespaces' in synapse_config:
                    weights_namespaces.extend(synapse_config['weights namespaces'])
                else:
                    weights_namespaces.append('Weights')
        else:
            has_weights = False

        if rank == 0:
            logger.info('*** Reading synaptic attributes of population %s' % (postsyn_name))

        syn_attributes_iter = read_cell_attribute_selection(forest_file_path, postsyn_name, selection=gid_range,
                                                            namespace='Synapse Attributes', comm=env.comm)
        
        syn_attributes_output_dict = dict(list(syn_attributes_iter))
        write_cell_attributes(write_selection_file_path, postsyn_name, syn_attributes_output_dict, namespace='Synapse Attributes', **write_kwds)
        del syn_attributes_output_dict
        del syn_attributes_iter
        
        if has_weights:
            for weights_namespace in sorted(weights_namespaces):
                weight_attributes_iter = read_cell_attribute_selection(forest_file_path, postsyn_name,
                                                                       selection=gid_range,
                                                                       namespace=weights_namespace, comm=env.comm)
                weight_attributes_output_dict = dict(list(weight_attributes_iter))
                write_cell_attributes(write_selection_file_path, postsyn_name, weight_attributes_output_dict, namespace=weights_namespace, **write_kwds)
                del weight_attributes_output_dict
                del weight_attributes_iter

                
        logger.info('*** Rank %i: reading connectivity selection from file %s for postsynaptic population: %s: selection: %s' % (rank, connectivity_file_path, postsyn_name, str(gid_range)))

        (graph, attr_info) = read_graph_selection(connectivity_file_path, selection=gid_range, \
                                                  projections=[ (presyn_name, postsyn_name) for presyn_name in sorted(presyn_names) ], \
                                                  comm=env.comm, namespaces=['Synapses', 'Connections'])
        

        for presyn_name in sorted(presyn_names):
            gid_dict = {}
            edge_count = 0
            node_count = 0
            if postsyn_name in graph:

                if postsyn_name in attr_info and presyn_name in attr_info[postsyn_name]:
                    edge_attr_info = attr_info[postsyn_name][presyn_name]
                else:
                    raise RuntimeError('write_connection_selection: missing edge attributes for projection %s -> %s' % \
                                       (presyn_name, postsyn_name))
                
                if 'Synapses' in edge_attr_info and \
                        'syn_id' in edge_attr_info['Synapses'] and \
                        'Connections' in edge_attr_info and \
                        'distance' in edge_attr_info['Connections']:
                    syn_id_attr_index = edge_attr_info['Synapses']['syn_id']
                    distance_attr_index = edge_attr_info['Connections']['distance']
                else:
                    raise RuntimeError('write_connection_selection: missing edge attributes for projection %s -> %s' % \
                                           (presyn_name, postsyn_name))
            
                edge_iter = compose_iter(lambda edgeset: input_sources[presyn_name].update(edgeset[1][0]), \
                                         graph[postsyn_name][presyn_name])
                for (postsyn_gid, edges) in edge_iter:

                    presyn_gids, edge_attrs = edges
                    edge_syn_ids = edge_attrs['Synapses'][syn_id_attr_index]
                    edge_dists = edge_attrs['Connections'][distance_attr_index]
                    
                    gid_dict[postsyn_gid] = (presyn_gids,
                                             {'Synapses': {'syn_id': edge_syn_ids},
                                              'Connections': {'distance': edge_dists} })
                    edge_count += len(presyn_gids)
                    node_count += 1

            logger.info('*** Rank %d: Writing projection %s -> %s selection: %d nodes, %d edges' % (rank, presyn_name, postsyn_name, node_count, edge_count))
            write_graph(write_selection_file_path, \
                        src_pop_name=presyn_name, dst_pop_name=postsyn_name, \
                        edges=gid_dict, comm=env.comm, io_size=env.io_size)
        env.comm.barrier()

    return input_sources

                    
def write_input_cell_selection(env, input_sources, write_selection_file_path, populations=None, write_kwds={}):
    """
    Writes out predefined spike trains when only a subset of the network is instantiated.

    :param env: an instance of the `dentate.Env` class
    :param input_sources: a dictionary of the form { pop_name, gid_sources }
    """

    if 'comm' not in write_kwds:
        write_kwds['comm'] = env.comm
    if 'io_size' not in write_kwds:
        write_kwds['io_size'] = env.io_size

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    dataset_path = env.dataset_path
    input_file_path = env.data_file_path

    if populations is None:
        pop_names = sorted(env.celltypes.keys())
    else:
        pop_names = populations

    for pop_name, gid_range in sorted(viewitems(input_sources)):

        if pop_name not in pop_names:
            continue

        spikes_output_dict = {}

        if (env.cell_selection is not None) and (pop_name in env.cell_selection):
            local_gid_range = gid_range.difference(set(env.cell_selection[pop_name]))
        else:
            local_gid_range = gid_range

        gid_ranges = env.comm.allgather(local_gid_range)
        this_gid_range = set([])
        for gid_range in gid_ranges:
            for gid in gid_range:
                if gid % nhosts == rank:
                    this_gid_range.add(gid)


        has_spike_train = False
        spike_input_source_loc = []
        if (env.spike_input_attribute_info is not None) and (env.spike_input_ns is not None):
            if (pop_name in env.spike_input_attribute_info) and \
                    (env.spike_input_ns in env.spike_input_attribute_info[pop_name]):
                has_spike_train = True
                spike_input_source_loc.append((env.spike_input_path, env.spike_input_ns))
        if (env.cell_attribute_info is not None) and (env.spike_input_ns is not None):
            if (pop_name in env.cell_attribute_info) and \
                    (env.spike_input_ns in env.cell_attribute_info[pop_name]):
                has_spike_train = True
                spike_input_source_loc.append((input_file_path,env.spike_input_ns))

        if rank == 0:
            logger.info('*** Reading spike trains for population %s: %d cells: has_spike_train = %s' % (pop_name, len(this_gid_range), str(has_spike_train)))

        if has_spike_train:

            vecstim_attr_set = set(['t'])
            if env.spike_input_attr is not None:
                vecstim_attr_set.add(env.spike_input_attr)
            if 'spike train' in env.celltypes[pop_name]:
                vecstim_attr_set.add(env.celltypes[pop_name]['spike train']['attribute'])
                    
            cell_spikes_iters = [ read_cell_attribute_selection(input_path, pop_name, \
                                                                list(this_gid_range), \
                                                                namespace=input_ns, \
                                                                mask=vecstim_attr_set, \
                                                                comm=env.comm) for (input_path, input_ns) in spike_input_source_loc ]

                
            for cell_spikes_iter in cell_spikes_iters:
                spikes_output_dict.update(dict(list(cell_spikes_iter)))

        if rank == 0:
            logger.info('*** Writing spike trains for population %s: %s' % (pop_name, str(spikes_output_dict)))

                
        write_cell_attributes(write_selection_file_path, pop_name, spikes_output_dict,  \
                              namespace=env.spike_input_ns, **write_kwds)
