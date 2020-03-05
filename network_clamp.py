"""
Routines for Network Clamp simulation.
"""
import os, sys, copy, uuid, pprint
from collections import defaultdict
from mpi4py import MPI
import numpy as np
import click
from dentate import io_utils, spikedata, synapses, stimulus, cell_clamp
from dentate.cells import h, make_input_cell, register_cell, record_cell
from dentate.env import Env
from dentate.neuron_utils import h, configure_hoc_env, make_rec
from dentate.utils import is_interactive, Context, Closure, list_find, list_index, range, str, viewitems, zip_longest, get_module_logger
from dentate.cell_clamp import init_biophys_cell
from neuroh5.io import read_cell_attribute_selection, read_cell_attribute_info

# This logger will inherit its settings from the root logger, created in dentate.env
logger = get_module_logger(__name__)

context = Context()

def mpi_excepthook(type, value, traceback):
    """

    :param type:
    :param value:
    :param traceback:
    :return:
    """
    sys_excepthook(type, value, traceback)
    sys.stderr.flush()
    if MPI.COMM_WORLD.size > 1:
        MPI.COMM_WORLD.Abort(1)

sys_excepthook = sys.excepthook
sys.excepthook = mpi_excepthook


def distgfs_reduce_fun(xs):
    return xs[0]

def generate_weights(env, weight_source_rules, this_syn_attrs):
    """
    Generates synaptic weights according to the rules specified in the
    Weight Generator section of network clamp configuration.
    """
    weights_dict = {}

    if len(weight_source_rules) > 0:

        for presyn_id, weight_rule in viewitems(weight_source_rules):
            source_syn_dict = defaultdict(list)

            for syn_id, syn in viewitems(this_syn_attrs):
                this_presyn_id = syn.source.population
                this_presyn_gid = syn.source.gid
                if this_presyn_id == presyn_id:
                    source_syn_dict[this_presyn_gid].append(syn_id)

            if weight_rule['class'] == 'Sparse':
                weights_name = weight_rule['name']
                rule_params = weight_rule['params']
                fraction = rule_params['fraction']
                seed_offset = int(env.model_config['Random Seeds']['Sparse Weights'])
                seed = int(seed_offset + 1)
                weights_dict[presyn_id] = \
                    synapses.generate_sparse_weights(weights_name, fraction, seed, source_syn_dict)
            elif weight_rule['class'] == 'Log-Normal':
                weights_name = weight_rule['name']
                rule_params = weight_rule['params']
                mu = rule_params['mu']
                sigma = rule_params['sigma']
                clip = None
                if 'clip' in rule_params:
                    clip = rule_params['clip']
                seed_offset = int(env.model_config['Random Seeds']['GC Log-Normal Weights 1'])
                seed = int(seed_offset + 1)
                weights_dict[presyn_id] = \
                    synapses.generate_log_normal_weights(weights_name, mu, sigma, seed, source_syn_dict, clip=clip)
            elif weight_rule['class'] == 'Normal':
                weights_name = weight_rule['name']
                rule_params = weight_rule['params']
                mu = rule_params['mu']
                sigma = rule_params['sigma']
                seed_offset = int(env.model_config['Random Seeds']['GC Normal Weights'])
                seed = int(seed_offset + 1)
                weights_dict[presyn_id] = \
                    synapses.generate_normal_weights(weights_name, mu, sigma, seed, source_syn_dict)
            else:
                raise RuntimeError('network_clamp.generate_weights: unknown weight generator rule class %s' % \
                                   weight_rule['class'])

    return weights_dict



def init(env, pop_name, gid, spike_events_path, generate_inputs_pops=set([]), generate_weights_pops=set([]),
         spike_events_namespace='Spike Events', t_var='t', t_min=None, t_max=None, write_cell=False, plot_cell=False):
    """
    Instantiates a cell and all its synapses and connections and loads
    or generates spike times for all synaptic connections.

    :param env: an instance of env.Env
    :param pop_name: population name
    :param gid: list of cell gid
    :param spike_events_path:

    """
    if env.results_file_path is not None:
        io_utils.mkout(env, env.results_file_path)

    if env.cell_selection is None:
        env.cell_selection = {}
    selection = env.cell_selection.get(pop_name, [])
    env.cell_selection[pop_name] = [gid] + [selection]
    
    ## If specified, presynaptic spikes that only fall within this time range
    ## will be loaded or generated
    if t_max is None:
        t_range = None
    else:
        if t_min is None:
            t_range = [0.0, t_max]
        else:
            t_range = [t_min, t_max]

    ## Attribute namespace that contains recorded spike events
    namespace_id = spike_events_namespace

    ## Determine presynaptic populations that connect to this cell type
    presyn_names = env.projection_dict[pop_name]

    ## Load cell gid and its synaptic attributes and connection data
    cell = init_biophys_cell(env, pop_name, gid, write_cell=write_cell)

    ## Load spike times of presynaptic cells
    spkdata = spikedata.read_spike_events(spike_events_path, \
                                          presyn_names, \
                                          spike_events_namespace, \
                                          spike_train_attr_name=t_var, \
                                          time_range=t_range)
    spkindlst = spkdata['spkindlst']
    spktlst = spkdata['spktlst']
    spkpoplst = spkdata['spkpoplst']

    ## Organize spike times by index of presynaptic population and gid
    input_source_dict = {}
    weight_source_dict = {}
    for presyn_name in presyn_names:
        presyn_index = int(env.Populations[presyn_name])
        spk_pop_index = list_index(presyn_name, spkpoplst)
        if spk_pop_index is None:
            logger.warning("No spikes found for population %s in file %s" % (presyn_name, spike_events_path))
            continue
        spk_inds = spkindlst[spk_pop_index]
        spk_ts = spktlst[spk_pop_index]
        spk_ts += float(env.stimulus_config['Equilibration Duration'])

        
        if presyn_name in generate_inputs_pops:
            if (presyn_name in env.netclamp_config.input_generators):
                spike_generator = env.netclamp_config.input_generators[presyn_name]
            else:
                raise RuntimeError('network_clamp.init: no input generator specified for population %s' % presyn_name)
        else:
            spike_generator = None

        input_source_dict[presyn_index] = {'gen': spike_generator,
                                           'spiketrains': {'gid': spk_inds, 't': spk_ts, }}

        if presyn_name in generate_weights_pops:
            if (presyn_name in env.netclamp_config.weight_generators[pop_name]):
                weight_rule = env.netclamp_config.weight_generators[pop_name][presyn_name]
            else:
                raise RuntimeError(
                    'network_clamp.init: no weights generator rule specified for population %s' % presyn_name)
        else:
            weight_rule = None

        if weight_rule is not None:
            weight_source_dict[presyn_index] = weight_rule

    min_delay = float('inf')
    syn_attrs = env.synapse_attributes
    this_syn_attrs = syn_attrs[gid]
    weight_params = defaultdict(dict)
    for syn_id, syn in viewitems(this_syn_attrs):
        presyn_id = syn.source.population
        presyn_gid = syn.source.gid
        delay = syn.source.delay
        if presyn_id in input_source_dict:
            ## Load presynaptic spike times into the VecStim for each synapse;
            ## if spike_generator_dict contains an entry for the respective presynaptic population,
            ## then use the given generator to generate spikes.
            if not (presyn_gid in env.gidset):
                cell = make_input_cell(env, presyn_gid, presyn_id, input_source_dict)
                register_cell(env, presyn_id, presyn_gid, cell)

    source_weight_params = generate_weights(env, weight_source_dict, this_syn_attrs)

    for presyn_id, weight_params in viewitems(source_weight_params):
        weights_syn_ids = weight_params['syn_id']
        for syn_name in (syn_name for syn_name in weight_params if syn_name != 'syn_id'):
            weights_values = weight_params[syn_name]
            syn_attrs.add_mech_attrs_from_iter(gid, syn_name, \
                                               zip_longest(weights_syn_ids, \
                                                           [{'weight': x} for x in weights_values]))
    synapses.config_biophys_cell_syns(env, gid, pop_name, insert=True, insert_netcons=True, verbose=True)
    record_cell(env, pop_name, gid)

    biophys_cell = env.biophys_cells[pop_name][gid]

    if plot_cell:
        import dentate.plot
        from dentate.plot import plot_synaptic_attribute_distribution
        syn_attrs = env.synapse_attributes
        biophys_cell = env.biophys_cells[pop_name][gid]
        syn_name = 'AMPA'
        syn_mech_name = syn_attrs.syn_mech_names[syn_name]
        for param_name in ['weight', 'g_unit']:
            param_label = '%s; %s; %s' % (syn_name, syn_mech_name, param_name)
            plot_synaptic_attribute_distribution(biophys_cell, env, syn_name, param_name, filters=None, from_mech_attrs=True,
                                                 from_target_attrs=True, param_label=param_label,
                                                 export='syn_params_%d.h5' % gid, description='network_clamp', show=False,
                                                 svg_title="Synaptic parameters for gid %d" % (gid),
                                                 output_dir=env.results_path)
        
        
    cell = env.pc.gid2cell(gid)
    for sec in list(cell.all):
        h.psection(sec=sec)

        
    env.pc.set_maxstep(10)
    h.stdinit()

    if is_interactive:
        context.update(locals())


def run(env, cvode=False):
    """
    Runs network clamp simulation. Assumes that procedure `init` has been
    called with the network configuration provided by the `env`
    argument.

    :param env: instance of env.Env
    :param cvode: whether to use adaptive integration
    """

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    if env.recording_profile is not None:
        rec_dt = env.recording_profile.get('dt', 0.1) 
        env.t_rec.record(h._ref_t, rec_dt)
    env.t_vec.resize(0)
    env.id_vec.resize(0)

    st_comptime = env.pc.step_time()

    h.cvode_active(1 if cvode else 0)
    
    h.t = 0.0
    h.dt = env.dt
    h.tstop = env.tstop
    if 'Equilibration Duration' in env.stimulus_config:
        h.tstop += float(env.stimulus_config['Equilibration Duration'])
    h.finitialize(env.v_init)

    if rank == 0:
        logger.info("*** Running simulation with dt = %.03f and tstop = %.02f" % (h.dt, h.tstop))

    env.pc.barrier()
    env.pc.psolve(h.tstop)

    if rank == 0:
        logger.info("*** Simulation completed")
    env.pc.barrier()

    comptime = env.pc.step_time() - st_comptime
    avgcomp = env.pc.allreduce(comptime, 1) / nhosts
    maxcomp = env.pc.allreduce(comptime, 2)

    if rank == 0:
        logger.info("Host %i  ran simulation in %g seconds" % (rank, comptime))

    env.pc.runworker()
    env.pc.done()

    return spikedata.get_env_spike_dict(env, include_artificial=None)


def run_with(env, param_dict, cvode=False):
    """
    Runs network clamp simulation with the specified parameters for the given gid(s).
    Assumes that procedure `init` has been called with
    the network configuration provided by the `env` argument.

    :param env: instance of env.Env
    :param param_dict: dictionary { gid: params }
    :param cvode: whether to use adaptive integration
    """

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    for pop_name, gid_param_dict in viewitems(param_dict):
        biophys_cell_dict = env.biophys_cells[pop_name]

        synapse_config = env.celltypes[pop_name]['synapses']
        weights_dict = synapse_config.get('weights', {})
        param_expr_dict = {}
        if 'expr' in weights_dict:
            param_expr_dict['weight'] = weights_dict['expr']

        for gid, params_tuples in viewitems(gid_param_dict):
            for update_operator, destination, source, sec_type, syn_name, param_path, param_value in params_tuples:
                if isinstance(param_path, tuple):
                    p, s = param_path
                    if p in param_expr_dict:
                        param_expr = param_expr_dict[p]
                        param_expr[s] = param_value
            biophys_cell = biophys_cell_dict[gid]
            for update_operator, destination, source, sec_type, syn_name, param_path, param_value in params_tuples:
                if isinstance(param_path, tuple):
                    p, s = param_path
                else:
                    p, s = param_path, None

                sources = None
                if isinstance(source, tuple):
                    sources = list(source)
                else:
                    if source is not None:
                        sources = [source]
                synapses.modify_syn_param(biophys_cell, env, sec_type, syn_name,
                                          param_name=p, value=param_value,
                                          filters={'sources': sources} if sources is not None else None,
                                          update_operator=update_operator,
                                          origin='soma', update_targets=True)
            cell = env.pc.gid2cell(gid)

    if env.recording_profile is not None:
        rec_dt = env.recording_profile.get('dt', 0.1) 
        env.t_rec.record(h._ref_t, rec_dt)

    env.t_vec.resize(0)
    env.id_vec.resize(0)

    st_comptime = env.pc.step_time()

    h.cvode_active(1 if cvode else 0)

    h.t = 0.0
    h.tstop = env.tstop
    if 'Equilibration Duration' in env.stimulus_config:
        h.tstop += float(env.stimulus_config['Equilibration Duration'])

    h.dt = env.dt
    h.finitialize(env.v_init)

    if rank == 0:
        logger.info("*** Running simulation with dt = %.03f and tstop = %.02f" % (h.dt, h.tstop))
        logger.info("*** Parameters: %s" % pprint.pformat(param_dict))


    env.pc.barrier()
    env.pc.psolve(h.tstop)

    if rank == 0:
        logger.info("*** Simulation completed")
    env.pc.barrier()

    comptime = env.pc.step_time() - st_comptime
    avgcomp = env.pc.allreduce(comptime, 1) / nhosts
    maxcomp = env.pc.allreduce(comptime, 2)

    if rank == 0:
        logger.info("Host %i  ran simulation in %g seconds" % (rank, comptime))

    env.pc.runworker()
    env.pc.done()
    
    return spikedata.get_env_spike_dict(env, include_artificial=None)



def optimize_params(env, pop_name, param_type, param_config_name):
                        
    param_bounds = {}
    param_names = []
    param_initial_dict = {}
    param_range_tuples = []

    synapse_config = env.celltypes[pop_name]['synapses']
    weights_dict = synapse_config.get('weights', {})
    param_expr_dict = {}
    if 'expr' in weights_dict:
        param_expr_dict['weight'] = weights_dict['expr']

    if param_type == 'synaptic':
        if (pop_name in env.netclamp_config.optimize_parameters['synaptic']['Parameter ranges'][param_config_name]):
            opt_params = env.netclamp_config.optimize_parameters['synaptic']['Parameter ranges'][param_config_name]
            param_ranges = opt_params[pop_name]
        else:
            raise RuntimeError(
                "network_clamp.optimize_params: population %s does not have optimization configuration" % pop_name)
        keyfun = lambda kv: str(kv[0])
        for source, source_dict in sorted(viewitems(param_ranges), key=keyfun):
            for sec_type, sec_type_dict in sorted(viewitems(source_dict), key=keyfun):
                for syn_name, syn_mech_dict in sorted(viewitems(sec_type_dict), key=keyfun):
                    for param_fst, param_rst in sorted(viewitems(syn_mech_dict), key=keyfun):
                        if isinstance(param_rst, dict):
                            update_operator = lambda gid, syn_id, old, new: old
                            for const_name, const_range in sorted(viewitems(param_rst)):
                                param_path = (param_fst, const_name)
                                param_range_tuples.append((update_operator, pop_name, source, sec_type, syn_name, param_path, const_range))
                                param_key = '%s.%s.%s.%s.%s.%s' % (pop_name, str(source), sec_type, syn_name, param_fst, const_name)
                                param_initial_value = (const_range[1] - const_range[0]) / 2.0
                                param_initial_dict[param_key] = param_initial_value
                                param_bounds[param_key] = const_range
                                param_names.append(param_key)
                        else:
                            update_operator = None
                            param_name = param_fst
                            param_range = param_rst
                            param_range_tuples.append((update_operator, pop_name, source, sec_type, syn_name, param_name, param_range))
                            param_key = '%s.%s.%s.%s.%s' % (pop_name, source, sec_type, syn_name, param_name)
                            param_initial_value = (param_range[1] - param_range[0]) / 2.0
                            param_initial_dict[param_key] = param_initial_value
                            param_bounds[param_key] = param_range
                            param_names.append(param_key)
        
    else:
        raise RuntimeError("network_clamp.optimize_params: unknown parameter type %s" % param_type)

    return param_bounds, param_names, param_initial_dict, param_range_tuples


def init_state_objfun(config_file, population, gid, generate_inputs, generate_weights, t_max, t_min, tstop, opt_iter, template_paths, dataset_prefix, config_prefix, results_path, spike_events_path, spike_events_namespace, spike_events_t, param_type, param_config_name, recording_profile, state_variable, target_value, **kwargs):

    params = dict(locals())
    env = Env(**params)
    env.results_file_path = None
    configure_hoc_env(env)
    init(env, population, gid, spike_events_path, 
         generate_inputs_pops=set(generate_inputs), 
         generate_weights_pops=set(generate_weights), 
         spike_events_namespace=spike_events_namespace, 
         t_var=spike_events_t, t_min=t_min, t_max=t_max)

    time_step = env.stimulus_config['Temporal Resolution']
    equilibration_duration = float(env.stimulus_config['Equilibration Duration'])
    
    param_bounds, param_names, param_initial_dict, param_range_tuples = \
      optimize_params(env, population, param_type, param_config_name)
    
    def from_param_dict(params_dict):
        result = []
        for param_pattern, (update_operator, population, source, sec_type, syn_name, param_name, param_range) in zip(param_names, param_range_tuples):
            result.append((update_operator, population, source, sec_type, syn_name, param_name, params_dict[param_pattern]))
        return result

    def gid_state_value(spkdict, gid, t_offset, t_rec, state_recs):
        time_vec = np.asarray(t_rec.to_python(), dtype=np.float32) - t_offset
        t_inds = np.where(time_vec > 0.)[0]
        state_values = []
        for rec in state_recs:
            vec = np.asarray(rec['vec'].to_python(), dtype=np.float32)
            data = vec[t_inds]
            state_values.append(np.mean(data))
        m = np.mean(np.asarray(state_values))
        logger.info('state value objective: mean value of %s of gid %i is %.2f' % (state_variable, gid, m))
        return m

    recording_profile = { 'label': 'network_clamp.state.%s' % state_variable,
                          'dt': 0.1,
                          'section quantity': {
                              state_variable: { 'swc types': ['soma'] }
                            }
                        }
    env.recording_profile = recording_profile
    state_recs = record_cell(env, population, gid, recording_profile=recording_profile)

    f = lambda **v: (-abs(gid_state_value(run_with(env, {population: {gid: from_param_dict(v)}}), gid,
                            equilibration_duration, env.t_rec, state_recs) - target_value))
    return f


def init_rate_objfun(config_file, population, gid, generate_inputs, generate_weights, t_max, t_min, tstop, opt_iter, template_paths, dataset_prefix, config_prefix, results_path, spike_events_path, spike_events_namespace, spike_events_t, param_type, param_config_name, recording_profile, target_rate, **kwargs):

    params = dict(locals())
    env = Env(**params)
    env.results_file_path = None
    configure_hoc_env(env)
    init(env, population, gid, spike_events_path, 
         generate_inputs_pops=set(generate_inputs), 
         generate_weights_pops=set(generate_weights), 
         spike_events_namespace=spike_events_namespace, 
         t_var=spike_events_t, t_min=t_min, t_max=t_max)

    time_step = env.stimulus_config['Temporal Resolution']
    equilibration_duration = float(env.stimulus_config['Equilibration Duration'])

    param_bounds, param_names, param_initial_dict, param_range_tuples = \
      optimize_params(env, population, param_type, param_config_name)
    
    def from_param_dict(params_dict):
        result = []
        for param_pattern, (update_operator, population, source, sec_type, syn_name, param_name, param_range) in zip(param_names, param_range_tuples):
            result.append((update_operator, population, source, sec_type, syn_name, param_name, params_dict[param_pattern]))
        return result

    def gid_firing_rate(spkdict, gid):
        if gid in spkdict[pop_name]:
            spkdict1 = {gid: spkdict[pop_name][gid]}
        else:
            spkdict1 = {gid: np.asarray([], dtype=np.float32)}
        rate_dict = spikedata.spike_rates(spkdict1)
        if gid in spkdict[pop_name]:
            logger.info('firing rate objective: spikes times of gid %i: %s' % (gid, pprint.pformat(spkdict[pop_name][gid])))
        logger.info('firing rate objective: rate of gid %i is %.2f' % (gid, rate_dict[gid]))
        return rate_dict[gid]['rate']

    logger.info("firing rate objective: target rate: %.02f" % target_rate)

    f = lambda **v: (-abs(gid_firing_rate(run_with(env, {population: {gid: from_param_dict(v)}}), gid) -
                         target_rate))

    return f


def init_rate_dist_objfun(config_file, population, gid, generate_inputs, generate_weights, t_max, t_min, tstop, opt_iter, 
             template_paths, dataset_prefix, config_prefix, results_path, spike_events_path, spike_events_namespace, spike_events_t,
             param_type, param_config_name, recording_profile, target_rate_map_path, target_rate_map_namespace, target_rate_map_arena, target_rate_map_trajectory,
             **kwargs):

    params = dict(locals())
    env = Env(**params)
    env.results_file_path = None
    configure_hoc_env(env)
    init(env, population, gid, spike_events_path, 
         generate_inputs_pops=set(generate_inputs), 
         generate_weights_pops=set(generate_weights), 
         spike_events_namespace=spike_events_namespace, 
         t_var=spike_events_t, t_min=t_min, t_max=t_max)

    time_step = env.stimulus_config['Temporal Resolution']
    equilibration_duration = float(env.stimulus_config['Equilibration Duration'])

    input_namespace = '%s %s %s' % (target_rate_map_namespace, target_rate_map_arena, target_rate_map_trajectory)
    it = read_cell_attribute_selection(target_rate_map_path, population, namespace=input_namespace,
                                        selection=[gid], mask=set(['Trajectory Rate Map']))
    trj_rate_map = dict(it)[gid]['Trajectory Rate Map']

    trj_x, trj_y, trj_d, trj_t = stimulus.read_trajectory(target_rate_map_path, target_rate_map_arena, target_rate_map_trajectory)

    time_range = (0., min(np.max(trj_t), tstop))
    
    time_bins = np.arange(time_range[0], time_range[1], time_step)
    target_rate_vector = np.interp(time_bins, trj_t, trj_rate_map)

    
    param_bounds, param_names, param_initial_dict, param_range_tuples = \
      optimize_params(env, population, param_type, param_config_name)
    
    def from_param_dict(params_dict):
        result = []
        for param_pattern, (update_operator, population, source, sec_type, syn_name, param_name, param_range) in zip(param_names, param_range_tuples):
            result.append((update_operator, population, source, sec_type, syn_name, param_name, params_dict[param_pattern]))
        return result

    def gid_firing_rate_vector(spkdict, gid):
        if gid in spkdict[population]:
            spkdict1 = {gid: spkdict[population][gid]}
        else:
            spkdict1 = {gid: np.asarray([], dtype=np.float32)}
        rate_dict = spikedata.spike_rates(spkdict1)
        spike_density_dict = spikedata.spike_density_estimate (population, spkdict1, time_bins)
        if gid in spkdict[population]:
            rate = spike_density_dict[gid]['rate']
            logger.info('firing rate objective: spike times of gid %i: %s' % (gid, str(spkdict[population][gid])))
            logger.info('firing rate objective: firing rate of gid %i: %s' % (gid, str(rate)))
            logger.info('firing rate objective: min/max rates of gid %i are %.2f / %.2f Hz' % (gid, np.min(rate), np.max(rate)))
        return spike_density_dict[gid]['rate']
    logger.info("firing rate objective: target time bins: %s" % str(time_bins))
    logger.info("firing rate objective: target vector: %s" % str(target_rate_vector))
    logger.info("firing rate objective: target rate vector min/max is %.2f Hz (%.2f ms) / %.2f Hz (%.2f ms)" % (np.min(target_rate_vector), time_bins[np.argmin(target_rate_vector)], np.max(target_rate_vector), time_bins[np.argmax(target_rate_vector)]))
    f = lambda **v: (-np.square(np.subtract(gid_firing_rate_vector(run_with(env, {population: {gid: from_param_dict(v)}}), gid), target_rate_vector)).mean())
    
    return f


def optimize_run(env, pop_name, param_config_name, init_objfun,
                 opt_iter=10, solver_epsilon=1e-2, param_type='synaptic', init_params={}, 
                 results_file=None, verbose=False):
    import distgfs

    param_bounds, param_names, param_initial_dict, param_range_tuples = \
      optimize_params(env, pop_name, param_type, param_config_name)
    
    hyperprm_space = { param_pattern: [param_range[0], param_range[1]]
                       for param_pattern, (update_operator, pop_name, source, sec_type, syn_name, _, param_range) in
                           zip(param_names, param_range_tuples) }

    # Create an optimizer parameter set
    if results_file is None:
        if env.results_path is not None:
            file_path = '%s/distgfs.network_clamp.%s.h5' % (env.results_path, str(env.results_file_id))
        else:
            file_path = 'distgfs.network_clamp.%s.h5' % (str(env.results_file_id))
    else:
        file_path = '%s/%s' % (env.results_path, results_file)
    distgfs_params = {'opt_id': 'network_clamp.rate_dist',
                      'obj_fun_init_name': init_objfun, 
                      'obj_fun_init_module': 'dentate.network_clamp',
                      'obj_fun_init_args': init_params,
                      'reduce_fun_name': 'distgfs_reduce_fun',
                      'reduce_fun_module': 'dentate.network_clamp',
                      'problem_parameters': {},
                      'space': hyperprm_space,
                      'file_path': file_path,
                      'save': True,
                      'n_iter': opt_iter,
                      'solver_epsilon': solver_epsilon }

    opt_params, outputs = distgfs.run(distgfs_params, spawn_workers=True, verbose=verbose)
    
    logger.info('Optimized parameters: %s' % pprint.pformat(opt_params))
    logger.info('Optimized objective function: %s' % pprint.pformat(outputs))

    return opt_params, outputs

    
def netclamp_dist_ctrl(controller, init_params, cell_index_set):
    """Controller for distributed network clamp runs."""
    gfsopt = gfsinit(gfsopt_params)
    logger.info("Optimizing for %d iterations..." % gfsopt.n_iter)
    iter_count = 0
    task_ids = []
    while iter_count < gfsopt.n_iter:
        controller.recv()
        
        if (iter_count > 0) and gfsopt.save and (iter_count % gfsopt.save_iter == 0):
            gfsopt.save_evals()

        if len(task_ids) > 0:
            task_id, res = controller.get_next_result()
            
            if gfsopt.reduce_fun is None:
                rres = res
            else:
                rres = gfsopt.reduce_fun(res)
            eval_req = gfsopt.evals[task_id]
            vals = list(eval_req.x)
            eval_req.set(rres)
            task_ids.remove(task_id)
            iter_count += 1
            logger.info("optimization iteration %d: parameter coordinates %s: %s" % (iter_count, str(vals), str(rres)))
            
        while (len(controller.ready_workers) > 0) and (len(gfsopt.evals) < gfsopt.n_iter):
            eval_req = gfsopt.optimizer.get_next_x()
            vals = list(eval_req.x)
            task_id = controller.submit_call("eval_fun", module_name="distgfs",
                                             args=(gfsopt.opt_id, iter_count, vals,))
            task_ids.append(task_id)
            gfsopt.evals[task_id] = eval_req
                
    if gfsopt.save:
        gfsopt.save_evals()
    controller.info()

def netclamp_dist_work(worker, gfsopt_params):
    """Initialize workers for distributed network clamp runs."""
    gfsinit(gfsopt_params)
    

def write_output(env):
    rank = env.comm.rank
    if rank == 0:
        logger.info("*** Writing spike data")
    io_utils.spikeout(env, env.results_file_path)
    if rank == 0:
        logger.info("*** Writing intracellular data")
    io_utils.recsout(env, env.results_file_path)


@click.group()
def cli():
    pass


@click.command()
@click.option("--config-file", '-c', required=True, type=str, help='model configuration file name')
@click.option("--population", '-p', required=True, type=str, default='GC', help='target population')
@click.option("--gid", '-g', required=True, type=int, default=0, help='target cell gid')
@click.option("--template-paths", type=str, required=True,
              help='colon-separated list of paths to directories containing hoc cell templates')
@click.option("--dataset-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='path to directory containing required neuroh5 data files')
@click.option("--config-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              default='config',
              help='path to directory containing network and cell mechanism config files')
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), \
              help='path to directory where output files will be written')
@click.option("--spike-events-path", '-s', required=True, type=click.Path(),
              help='path to neuroh5 file containing spike times')
@click.option("--spike-events-namespace", type=str, default='Spike Events',
              help='namespace containing spike times')
@click.option("--spike-events-t", required=False, type=str, default='t',
              help='name of variable containing spike times')
@click.option('--plot-cell', is_flag=True, help='plot the distribution of weight and g_unit synaptic parameters')
@click.option('--write-cell', is_flag=True, help='write out selected cell tree morphology and connections')
@click.option('--profile-memory', is_flag=True, help='calculate and print heap usage after the simulation is complete')
@click.option('--recording-profile', type=str, default='Network clamp default', help='recording profile to use')

def show(config_file, population, gid, template_paths, dataset_prefix, config_prefix, results_path,
         spike_events_path, spike_events_namespace, spike_events_t, plot_cell, write_cell, profile_memory, recording_profile):
    """
    Show configuration for the specified cell.
    """

    np.seterr(all='raise')

    verbose = True
    init_params = dict(locals())
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        comm0 = env.comm.Split(2 if rank == 0 else 1, 0)
    
        env = Env(**init_params, comm=comm0)
        configure_hoc_env(env)

        init(env, population, gid, spike_events_path, \
            spike_events_namespace=spike_events_namespace, \
            t_var=spike_events_t, plot_cell=plot_cell, write_cell=write_cell)
            
        if env.profile_memory:
            profile_memory(logger)
            
    comm.barrier()

@click.command()
@click.option("--config-file", '-c', required=True, type=str, help='model configuration file name')
@click.option("--population", '-p', required=True, type=str, default='GC', help='target population')
@click.option("--gid", '-g', required=False, type=int, default=0, help='target cell gid')
@click.option("--generate-inputs", '-e', required=False, type=str, multiple=True,
              help='generate spike trains for the given presynaptic population')
@click.option("--generate-weights", '-w', required=False, type=str, multiple=True,
              help='generate weights for the given presynaptic population')
@click.option("--tstop", '-t', type=float, default=150.0, help='simulation end time')
@click.option("--t-max", type=float)
@click.option("--t-min", type=float)
@click.option("--template-paths", type=str, required=True,
              help='colon-separated list of paths to directories containing hoc cell templates')
@click.option("--dataset-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='path to directory containing required neuroh5 data files')
@click.option("--config-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              default='config',
              help='path to directory containing network and cell mechanism config files')
@click.option("--spike-events-path", '-s', required=True, type=click.Path(),
              help='path to neuroh5 file containing spike times')
@click.option("--spike-events-namespace", type=str, default='Spike Events',
              help='namespace containing spike times')
@click.option("--spike-events-t", required=False, type=str, default='t',
              help='name of variable containing spike times')
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), \
              help='path to directory where output files will be written')
@click.option("--results-file-id", type=str, required=False, default=None, \
              help='identifier that is used to name neuroh5 files that contain output spike and intracellular trace data')
@click.option("--results-namespace-id", type=str, required=False, default=None, \
              help='identifier that is used to name neuroh5 namespaces that contain output spike and intracellular trace data')
@click.option('--plot-cell', is_flag=True, help='plot the distribution of weight and g_unit synaptic parameters')
@click.option('--write-cell', is_flag=True, help='write out selected cell tree morphology and connections')
@click.option('--profile-memory', is_flag=True, help='calculate and print heap usage after the simulation is complete')
@click.option('--recording-profile', type=str, default='Network clamp default', help='recording profile to use')

def go(config_file, population, gid, generate_inputs, generate_weights, tstop, t_max, t_min,
       template_paths, dataset_prefix, config_prefix, spike_events_path, spike_events_namespace, spike_events_t,
       results_path, results_file_id, results_namespace_id, plot_cell, write_cell,
       profile_memory, recording_profile):

    """
    Runs network clamp simulation for the specified gid, or for all gids found in the input data file.
    """
    
    if results_file_id is None:
        results_file_id = uuid.uuid4()
    init_params = dict(locals())
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    np.seterr(all='raise')
    verbose = True

    if results_file_id is None:
        if rank == 0:
            results_file_id = uuid.uuid4()
        
        results_file_id = comm.bcast(results_file_id, root=0)
        init_params['results_file_id'] = results_file_id
    
    cell_index_set = set([])
    if gid is None:
        cell_index_data = None
        if rank == 0:
            comm0 = env.comm.Split(2 if rank == 0 else 1, 0)
            env = Env(**init_params, comm=comm0)
            attr_info_dict = read_cell_attribute_info(env.data_file_path, populations=[population],
                                                      read_cell_index=True, comm=comm0)
            cell_index = None
            attr_name, attr_cell_index = next(iter(attr_info_dict[pop_name]['Trees']))
            cell_index_set = set(attr_cell_index)
        cell_index_set = comm.bcast(cell_index_set, root=0)
    else:
        cell_index_set.add(gid)

    comm.barrier()
        
    if size > 1:
        import distwq
        if distwq.is_controller:
            distwq.run(fun_name="netclamp_dist_ctrl", module_name="dentate.network_clamp",
                       verbose=True, args=(init_params, cell_index_set),
                       spawn_workers=True, nprocs_per_worker=1)

        else:
            distwq.run(fun_name="netclamp_dist_work", module_name="dentate.network_clamp",
                       verbose=True, args=(init_params, cell_index_set),
                       spawn_workers=True, nprocs_per_worker=1)
    else:
        env = Env(**init_params, comm=comm)
        for gid in cell_index:
            init(env, population, gid, spike_events_path, generate_inputs_pops=set(generate_inputs),
                generate_weights_pops=set(generate_weights), spike_events_namespace=spike_events_namespace,
                t_var=spike_events_t, t_min=t_min, t_max=t_max,
                plot_cell=plot_cell, write_cell=write_cell)
            run(env)
            write_output(env)
        if env.profile_memory:
            profile_memory(logger)


@click.command()
@click.option("--config-file", '-c', required=True, type=str, help='model configuration file name')
@click.option("--population", '-p', required=True, type=str, default='GC', help='target population')
@click.option("--gid", '-g', required=True, type=int, default=0, help='target cell gid')
@click.option("--generate-inputs", '-e', required=False, type=str, multiple=True,
              help='generate spike trains for the given presynaptic population')
@click.option("--generate-weights", '-w', required=False, type=str, multiple=True,
              help='generate weights for the given presynaptic population')
@click.option("--tstop", '-t', type=float, default=150.0, help='simulation end time')
@click.option("--t-max", type=float)
@click.option("--t-min", type=float)
@click.option("--opt-iter", type=int, default=10, help='number of optimization iterations')
@click.option("--template-paths", type=str, required=True,
              help='colon-separated list of paths to directories containing hoc cell templates')
@click.option("--dataset-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              help='path to directory containing required neuroh5 data files')
@click.option("--config-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              default='config',
              help='path to directory containing network and cell mechanism config files')
@click.option("--param-config-name", type=str, 
              help='parameter configuration name to use for optimization (defined in config file)')
@click.option("--param-type", type=str, default='synaptic',
              help='parameter type to use for optimization (synaptic)')
@click.option('--recording-profile', type=str, help='recording profile to use')
@click.option("--results-file", required=False, type=str, help='optimization results file')
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), \
              help='path to directory where output files will be written')
@click.option("--spike-events-path", '-s', required=True, type=click.Path(),
              help='path to neuroh5 file containing spike times')
@click.option("--spike-events-namespace", type=str, required=False, default='Spike Events',
              help='namespace containing spike times')
@click.option("--spike-events-t", required=False, type=str, default='t',
              help='name of variable containing spike times')
@click.option("--target-rate-map-path", required=False, type=click.Path(),
              help='path to neuroh5 file containing target rate maps used for rate optimization')
@click.option("--target-rate-map-namespace", type=str, required=False, default='Input Spikes',
              help='namespace containing target rate maps used for rate optimization')
@click.option("--target-rate-map-arena", type=str, required=False, 
              help='name of arena used for rate optimization')
@click.option("--target-rate-map-trajectory", type=str, required=False, 
              help='name of trajectory used for rate optimization')
@click.option("--target-state-variable", type=str, required=False, 
              help='name of state variable used for state optimization')
@click.argument('target')# help='rate, rate_dist, state'


def optimize(config_file, population, gid, generate_inputs, generate_weights, t_max, t_min, tstop, opt_iter, 
             template_paths, dataset_prefix, config_prefix, spike_events_path, spike_events_namespace, spike_events_t,
             param_config_name, param_type, recording_profile, results_file, results_path, target_rate_map_path, target_rate_map_namespace, target_rate_map_arena, target_rate_map_trajectory, target_state_variable, target):
    """
    Optimize the firing rate of the specified cell in a network clamp configuration.
    """
    init_params = dict(locals())

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    results_file_id = None
    if rank == 0:
        results_file_id = uuid.uuid4()
        
    results_file_id = comm.bcast(results_file_id, root=0)
    
    np.seterr(all='raise')
    verbose = False
    if rank == 0 or rank == 1:
        verbose = True
    cache_queries = True
    params = dict(locals())
    env = Env(**params)

    if size == 1:
        configure_hoc_env(env)
        init(env, population, gid, spike_events_path, 
            generate_inputs_pops=set(generate_inputs), 
            generate_weights_pops=set(generate_weights), 
            spike_events_namespace=spike_events_namespace, 
            t_var=spike_events_t, t_min=t_min, t_max=t_max)
        
    if (population in env.netclamp_config.optimize_parameters[param_type]):
        opt_params = env.netclamp_config.optimize_parameters[param_type][population]
    else:
        raise RuntimeError(
            "network_clamp.optimize: population %s does not have optimization configuration" % population)

    if target == 'rate':
        opt_target = opt_params['Targets']['firing rate']
        init_params['target_rate'] = opt_target
        init_objfun_name = 'init_rate_dist_objfun'
    elif target == 'ratedist' or target == 'rate_dist':
        init_objfun_name = 'init_rate_dist_objfun'
    elif target == 'state':
        opt_target = opt_params['Targets']['state'][state_variable]
        init_params['target_value'] = opt_target
        init_params['state_variable'] = target_state_variable
        init_objfun_name = 'init_state_objfun'
    else:
        raise RuntimeError('network_clamp.optimize: unknown optimization target %s' % target) 
        
    optimize_run(env, population, param_config_name, init_objfun_name,
                 opt_iter=opt_iter, param_type=param_type,
                 init_params=init_params, results_file=results_file,
                 verbose=verbose)


cli.add_command(show)
cli.add_command(go)
cli.add_command(optimize)

if __name__ == '__main__':

    cli(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):],
        standalone_mode=False)
