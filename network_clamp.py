"""
Routines for Network Clamp simulation.
"""
import click
from collections import defaultdict
from dentate.utils import *
from dentate.neuron_utils import *
from dentate.env import Env
from dentate.cells import *
from dentate import spikedata, io_utils, synapses


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
                seed_offset = int(env.modelConfig['Random Seeds']['Sparse Weights'])
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
                seed_offset = int(env.modelConfig['Random Seeds']['GC Log-Normal Weights 1'])
                seed = int(seed_offset + 1)
                weights_dict[presyn_id] = \
                    synapses.generate_log_normal_weights(weights_name, mu, sigma, seed, source_syn_dict, clip=clip)
            elif weight_rule['class'] == 'Normal':
                weights_name = weight_rule['name']
                rule_params = weight_rule['params']
                mu = rule_params['mu']
                sigma = rule_params['sigma']
                seed_offset = int(env.modelConfig['Random Seeds']['GC Normal Weights'])
                seed = int(seed_offset + 1)
                weights_dict[presyn_id] = \
                    synapses.generate_normal_weights(weights_name, mu, sigma, seed, source_syn_dict)
            else:
                raise RuntimeError('network_clamp.generate_weights: unknown weight generator rule class %s' % \
                                   weight_rule['class'])

    return weights_dict


def load_cell(env, pop_name, gid, mech_file_path=None, correct_for_spines=False, load_edges=True, tree_dict=None,
              load_synapses=True, synapses_dict=None):
    """
    Instantiates the mechanisms of a single cell.

    :param env: env.Env
    :param pop_name: str
    :param gid: int
    :param mech_file_path: str; path to cell mechanism config file
    :param correct_for_spines: bool

    Environment can be instantiated as:
    env = Env(config_file, template_paths, dataset_prefix, config_prefix)
    :param template_paths: str; colon-separated list of paths to directories containing hoc cell templates
    :param dataset_prefix: str; path to directory containing required neuroh5 data files
    :param config_prefix: str; path to directory containing network and cell mechanism config files
    
    """
    configure_hoc_env(env)

    cell = get_biophys_cell(env, pop_name, gid, load_edges=load_edges, \
                            tree_dict=tree_dict, load_synapses=load_synapses,
                            synapses_dict=synapses_dict)

    if mech_file_path is not None:
        init_biophysics(cell, reset_cable=True, from_file=True, mech_file_path=mech_file_path,
                        correct_cm=correct_for_spines, correct_g_pas=correct_for_spines, env=env)
    synapses.init_syn_mech_attrs(cell, env)

    env.biophys_cells[pop_name][gid] = cell

    return cell


def register_cell(env, population, gid, cell):
    """
    Registers a cell in a ParallelContext network environment.

    :param env: an instance of env.Env
    :param population: population name
    :param gid: gid
    :param cell: cell instance
    """
    rank = env.comm.rank
    env.gidset.add(gid)
    env.cells.append(cell)
    env.pc.set_gid2node(gid, rank)
    # Tell the ParallelContext that this cell is a spike source
    # for all other hosts. NetCon is temporary.
    hoc_cell = getattr(cell, "hoc_cell", None)
    if hoc_cell is None:
        nc = cell.connect2target(h.nil)
    else:
        nc = cell.hoc_cell.connect2target(h.nil)
    env.pc.cell(gid, nc, 1)
    # Record spikes of this cell
    env.pc.spike_record(gid, env.t_vec, env.id_vec)


def init_cell(env, pop_name, gid, load_edges=True):
    """
    Instantiates a cell and all its synapses.

    :param env: an instance of env.Env
    :param pop_name: population name
    :param gid: gid
    """

    ## Determine if a mechanism configuration file exists for this cell type
    if 'mech_file_path' in env.celltypes[pop_name]:
        mech_file_path = env.celltypes[pop_name]['mech_file_path']
    else:
        mech_file_path = None

    ## Determine if correct_for_spines flag has been specified for this cell type
    synapse_config = env.celltypes[pop_name]['synapses']
    if 'correct_for_spines' in synapse_config:
        correct_for_spines_flag = synapse_config['correct_for_spines']
    else:
        correct_for_spines_flag = False

    ## Determine presynaptic populations that connect to this cell type
    presyn_names = env.projection_dict[pop_name]

    ## Load cell gid and its synaptic attributes and connection data
    cell = load_cell(env, pop_name, gid, mech_file_path=mech_file_path, \
                     correct_for_spines=correct_for_spines_flag, \
                     load_edges=load_edges)
    register_cell(env, pop_name, gid, cell)

    rec = make_rec(0, pop_name, gid, cell, \
                   sec=cell.soma[0].sec, loc=0.5, param='v', \
                   dt=h.dt, description='Soma')
    env.recs_dict[pop_name]['Soma'].append(rec)
    if len(cell.hillock) > 0:
        rec = make_rec(1, pop_name, gid, cell, \
                       sec=cell.hillock[0].sec, loc=0.5, param='v', \
                       dt=h.dt, description='Axon hillock')
        env.recs_dict[pop_name]['Axon hillock'].append(rec)
    if len(cell.apical) > 0:
        n = len(cell.apical)
        for i in range(n):
            rec = make_rec(i + 2, pop_name, gid, cell, \
                           sec=cell.apical[i].sec, loc=0.5, param='v', \
                           dt=h.dt, description='Apical dendrite')
            env.recs_dict[pop_name]['Apical dendrite'].append(rec)

    report_topology(cell, env)

    return cell


def init(env, pop_name, gid, spike_events_path, generate_inputs_pops=set([]), generate_weights_pops=set([]),
         spike_events_namespace='Spike Events', t_var='t', t_min=None, t_max=None):
    """
    Instantiates a cell and all its synapses and connections and loads
    or generates spike times for all synaptic connections.

    :param env: an instance of env.Env
    :param pop_name: population name
    :param gid: gid
    :param spike_events_path:

    """
    io_utils.mkout(env, env.results_file_path)

    ## If specified, presynaptic spikes that only fall within this time range
    ## will be loaded or generated
    if t_max is None:
        t_range = None
    else:
        if t_min is None:
            t_range = [0.0, t_max]
        else:
            r_range = [t_min, t_max]

    ## Attribute namespace that contains recorded spike events
    if env.results_id is None:
        namespace_id = spike_events_namespace
    else:
        namespace_id = "%s %s" % (spike_events_namespace, str(env.results_id))

    ## Determine presynaptic populations that connect to this cell type
    presyn_names = env.projection_dict[pop_name]

    ## Load cell gid and its synaptic attributes and connection data
    cell = init_cell(env, pop_name, gid)

    ## Load spike times of presynaptic cells
    spkdata = spikedata.read_spike_events(spike_events_path, \
                                          presyn_names, \
                                          spike_events_namespace, \
                                          time_variable=t_var, \
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
    source_syn_dict = defaultdict(lambda: defaultdict(list))
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

    cell = env.pc.gid2cell(gid)
    for sec in list(cell.all):
        h.psection(sec=sec)

    env.pc.set_maxstep(10)
    h.stdinit()


def run(env):
    """
    Runs network clamp simulation. Assumes that procedure `init` has been
    called with the network configuration provided by the `env`
    argument.

    :param env:
    """

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    env.t_vec.resize(0)
    env.id_vec.resize(0)

    h.t = 0.0
    h.tstop = env.tstop

    h.finitialize(env.v_init)

    if rank == 0:
        logger.info("*** Running simulation with dt = %f and tstop = %f" % (h.dt, h.tstop))

    env.pc.barrier()
    env.pc.psolve(h.tstop)

    if rank == 0:
        logger.info("*** Simulation completed")
    env.pc.barrier()

    comptime = env.pc.step_time()
    cwtime = comptime + env.pc.step_wait()
    maxcw = env.pc.allreduce(cwtime, 2)
    avgcomp = old_div(env.pc.allreduce(comptime, 1), nhosts)
    maxcomp = env.pc.allreduce(comptime, 2)

    if rank == 0:
        logger.info("Host %i  ran simulation in %g seconds" % (rank, comptime))

    env.pc.runworker()
    env.pc.done()


def run_with(env, param_dict):
    """
    Runs network clamp simulation with the specified parameters for
    the given gid(s).  Assumes that procedure `init` has been called with
    the network configuration provided by the `env` argument.

    :param env:
    :param param_dict: dictionary { gid: params }

    """

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    for pop_name, gid_param_dict in viewitems(param_dict):
        biophys_cell_dict = env.biophys_cells[pop_name]
        for gid, params_tuples in viewitems(gid_param_dict):
            biophys_cell = biophys_cell_dict[gid]
            for source, sec_type, syn_name, param_name, param_value in params_tuples:
                synapses.modify_syn_param(biophys_cell, env, sec_type, syn_name,
                                          param_name=param_name, value=param_value,
                                          filters={'sources': [source]},
                                          origin='soma', update_targets=True)
            cell = env.pc.gid2cell(gid)

    env.t_vec.resize(0)
    env.id_vec.resize(0)

    # h.cvode_active(1)

    h.t = 0.0
    h.tstop = env.tstop

    h.finitialize(env.v_init)

    if rank == 0:
        logger.info("*** Running simulation with dt = %f and tstop = %f" % (h.dt, h.tstop))
        logger.info("*** Parameters: %s" % str(param_dict))

    env.pc.barrier()
    env.pc.psolve(h.tstop)

    if rank == 0:
        logger.info("*** Simulation completed")
    env.pc.barrier()

    comptime = env.pc.step_time()
    cwtime = comptime + env.pc.step_wait()
    maxcw = env.pc.allreduce(cwtime, 2)
    avgcomp = old_div(env.pc.allreduce(comptime, 1), nhosts)
    maxcomp = env.pc.allreduce(comptime, 2)

    if rank == 0:
        logger.info("Host %i  ran simulation in %g seconds" % (rank, comptime))

    env.pc.runworker()
    env.pc.done()

    return spikedata.get_env_spike_dict(env)


def make_firing_rate_target(env, pop_name, gid, target_rate, from_param_vector):
    def gid_firing_rate(spkdict, gid):
        if gid in spkdict[pop_name]:
            spkdict1 = {gid: spkdict[pop_name][gid]['t']}
        else:
            spkdict1 = {gid: np.asarray([], dtype=np.float32)}
        rate_dict = spikedata.spike_rates(spkdict1)
        if gid in spkdict[pop_name]:
            logger.info('firing rate objective: spikes times of gid %i: %s' % (gid, str(spkdict[pop_name][gid]['t'])))
        logger.info('firing rate objective: rate of gid %i is %.2f' % (gid, rate_dict[gid]))
        return rate_dict[gid]

    f = lambda *v: (abs(gid_firing_rate(run_with(env, {pop_name: {gid: from_param_vector(v)}}), gid) - target_rate))

    return f


def optimize_rate(env, pop_name, gid, opt_iter=10):
    import dlib

    if (pop_name in env.netclamp_config.optimize_parameters):
        opt_params = env.netclamp_config.optimize_parameters[pop_name]
        param_ranges = opt_params['Parameter ranges']
        opt_target = opt_params['Target firing rate']
    else:
        raise RuntimeError(
            "network_clamp.optimize_rate: population %s does not have optimization configuration" % pop_name)

    param_range_tuples = []
    for source, source_dict in sorted(viewitems(param_ranges), key=lambda k_v3: k_v3[0]):
        for sec_type, sec_type_dict in sorted(viewitems(source_dict), key=lambda k_v2: k_v2[0]):
            for syn_name, syn_mech_dict in sorted(viewitems(sec_type_dict), key=lambda k_v1: k_v1[0]):
                for param_name, param_range in sorted(viewitems(syn_mech_dict), key=lambda k_v: k_v[0]):
                    param_range_tuples.append((source, sec_type, syn_name, param_name, param_range))

    min_values = [(source, sec_type, syn_name, param_name, param_range[0]) for
                  source, sec_type, syn_name, param_name, param_range in param_range_tuples]
    max_values = [(source, sec_type, syn_name, param_name, param_range[1]) for
                  source, sec_type, syn_name, param_name, param_range in param_range_tuples]

    def from_param_vector(params):
        result = []
        assert (len(params) == len(param_range_tuples))
        for i, (source, sec_type, syn_name, param_name, param_range) in enumerate(param_range_tuples):
            result.append((source, sec_type, syn_name, param_name, params[i]))
        return result

    def to_param_vector(params):
        result = []
        for (source, sec_type, syn_name, param_name, param_value) in params:
            result.append(param_value)
        return result

    f_firing_rate = make_firing_rate_target(env, pop_name, gid, opt_target, from_param_vector)
    opt_params, outputs = dlib.find_min_global(f_firing_rate, to_param_vector(min_values), to_param_vector(max_values),
                                               opt_iter)

    logger.info('Optimized parameters: %s' % str(from_param_vector(opt_params)))
    logger.info('Optimized objective function: %s' % str(outputs))

    return opt_params, outputs


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
@click.option("--spike-events-path", '-s', required=True, type=click.Path(),
              help='path to neuroh5 file containing spike times')
@click.option("--spike-events-namespace", type=str, default='Spike Events',
              help='namespace containing spike times')
@click.option('--profile-memory', is_flag=True, help='calculate and print heap usage after the simulation is complete')
def show(config_file, population, gid, tstop, template_paths, dataset_prefix, config_prefix, spike_events_path,
         spike_events_namespace, profile_memory):
    """
    Show configuration for the specified cell.
    """

    comm = MPI.COMM_WORLD
    np.seterr(all='raise')

    verbose = True
    params = dict(locals())
    env = Env(**params)
    configure_hoc_env(env)

    init_cell(env, population, gid, load_edges=False)

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
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), \
              help='path to directory where output files will be written')
@click.option("--results-id", type=str, required=False, default=None, \
              help='identifier that is used to name neuroh5 namespaces that contain output spike and intracellular trace data')
@click.option('--profile-memory', is_flag=True, help='calculate and print heap usage after the simulation is complete')
def go(config_file, population, gid, generate_inputs, generate_weights, tstop, template_paths, dataset_prefix,
       config_prefix, spike_events_path, spike_events_namespace, results_path, results_id, profile_memory):
    """
    Runs network clamp simulation for the specified cell.
    """

    comm = MPI.COMM_WORLD
    np.seterr(all='raise')
    verbose = True
    params = dict(locals())
    env = Env(**params)
    configure_hoc_env(env)

    init(env, population, gid, spike_events_path, \
         generate_inputs_pops=set(generate_inputs), \
         generate_weights_pops=set(generate_weights), \
         spike_events_namespace=spike_events_namespace, \
         t_var='t', t_min=None, t_max=None)

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
@click.option("--opt-iter", type=int, default=10, help='number of optimization iterations')
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
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True), \
              help='path to directory where output files will be written')
@click.argument('target')
def optimize(config_file, population, gid, generate_inputs, generate_weights, tstop, opt_iter, template_paths,
             dataset_prefix, config_prefix, spike_events_path, spike_events_namespace, results_path, target):
    """
    Optimize the firing rate of the specified cell in a network clamp configuration.
    """

    comm = MPI.COMM_WORLD
    np.seterr(all='raise')
    verbose = True
    cache_queries = True
    params = dict(locals())
    env = Env(**params)
    configure_hoc_env(env)

    init(env, population, gid, spike_events_path, \
         generate_inputs_pops=set(generate_inputs), \
         generate_weights_pops=set(generate_weights), \
         spike_events_namespace=spike_events_namespace, \
         t_var='t', t_min=None, t_max=None)

    if target == 'rate':
        optimize_rate(env, population, gid, opt_iter=opt_iter)
    else:
        raise RuntimeError('network_clamp.optimize: unknown optimization target %s' % \
                           target)


cli.add_command(show)
cli.add_command(go)
cli.add_command(optimize)

if __name__ == '__main__':
    cli(args=sys.argv[(list_find(lambda s: s.find(os.path.basename(__file__)) != -1, sys.argv) + 1):],
        standalone_mode=False)
