
import click, gc, math
from neuroh5.io import NeuroH5CellAttrGen, append_cell_attributes, read_population_ranges
from dentate.env import Env
from dentate.stimulus import InputSelectivityConfig, choose_input_selectivity_type, get_2D_arena_spatial_mesh, \
    get_2D_arena_bounds, get_2D_arena_extents, generate_input_selectivity_features
from dentate.utils import *
from mpi4py import MPI
import h5py

logger = get_script_logger(os.path.basename(__file__))

context = Struct(**dict(locals()))

sys_excepthook = sys.excepthook
def mpi_excepthook(type, value, traceback):
    sys_excepthook(type, value, traceback)
    sys.stdout.flush()
    sys.stderr.flush()
    if MPI.COMM_WORLD.size > 1:
        MPI.COMM_WORLD.Abort(1)
sys.excepthook = mpi_excepthook

def roundf( n, f ):
 
    # Smaller multiple
    a = (n // f) * f
     
    # Larger multiple
    b = a + f
     
    # Return of closest of two
    return (b if n - a > b - n else a)

def debug_callback(context):
    from dentate.plot import plot_2D_rate_map, close_figure
    fig_title = '%s %s cell %i' % (context.population, context.this_selectivity_type_name, context.gid)
    fig_options = copy.copy(context.fig_options)
    if context.save_fig is not None:
        fig_options.saveFig = '%s %s' % (context.save_fig, fig_title)
    fig = plot_2D_rate_map(x=context.arena_x_mesh, y=context.arena_y_mesh, rate_map=context.rate_map,
                           peak_rate = context.env.stimulus_config['Peak Rate'][context.population][context.this_selectivity_type],
                           title='%s\nNormalized cell position: %.3f' % (fig_title, context.norm_u_arc_distance),
                           **fig_options())
    close_figure(fig)
    
def merge_dict(d1, d2, datatype):
    dd = {}
    for d in (d1, d2):
        dd.update(d1)
    return dd

def merge_count_dict(d1, d2, datatype):
    dd = defaultdict(lambda: 0)
    for d in (d1, d2): 
        for key, value in d.items():
            dd[key] += value
    return dict(dd.items())

def concatenate_ndarray_dict(a, b, datatype):
    merged_ndarray_dict = {}
    for k in a:
        if k not in b:
            merged_rate_map_dict[k] = a[k]
        else:
            if len(b[k]) == 0:
                merged_ndarray_dict[k] = a[k]
            elif len(a[k]) == 0:
                merged_ndarray_dict[k] = b[k]
            else:
                merged_ndarray_dict[k] = np.concatenate((a[k], b[k]), axis=None)
    for k in b:
        if k not in a:
            merged_rate_map_dict[k] = b[k]
    return merged_ndarray_dict

def merge_rate_map_dict(m1, m2, datatype):
    merged_rate_map_dict = {}
    for population in m1:
        if population not in merged_rate_map_dict:
            merged_rate_map_dict[population] = {}
        for selectivity_type_name in m1[population]:
            srm = m1[population][selectivity_type_name]
            if selectivity_type_name in merged_rate_map_dict[population]:
                merged_rate_map_dict[population][selectivity_type_name] += srm
            else:
                merged_rate_map_dict[population][selectivity_type_name] = srm.copy()
    for population in m2:
        if population not in merged_rate_map_dict:
            merged_rate_map_dict[population] = {}
        for selectivity_type_name in m2[population]:
            srm = m2[population][selectivity_type_name]
            if selectivity_type_name in merged_rate_map_dict[population]:
                merged_rate_map_dict[population][selectivity_type_name] += srm
            else:
                merged_rate_map_dict[population][selectivity_type_name] = srm.copy()
    return merged_rate_map_dict

mpi_op_merge_dict = MPI.Op.Create(merge_dict, commute=True)
mpi_op_merge_count_dict = MPI.Op.Create(merge_count_dict, commute=True)
mpi_op_merge_rate_map_dict = MPI.Op.Create(merge_rate_map_dict, commute=True)
mpi_op_concatenate_ndarray_dict = MPI.Op.Create(concatenate_ndarray_dict, commute=True)
    

@click.command()
@click.option("--config", required=True, type=str)
@click.option("--config-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True),
              default='config')
@click.option("--coords-path", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--distances-namespace", '-n', type=str, default='Arc Distances')
@click.option("--output-path", type=click.Path(file_okay=True, dir_okay=False), default=None)
@click.option("--arena-id", type=str, default='A')
@click.option("--populations", '-p', type=str, multiple=True)
@click.option("--use-noise-gen", is_flag=True, default=False)
@click.option("--io-size", type=int, default=-1)
@click.option("--chunk-size", type=int, default=1000)
@click.option("--value-chunk-size", type=int, default=1000)
@click.option("--cache-size", type=int, default=50)
@click.option("--write-size", type=int, default=10000)
@click.option("--verbose", '-v', is_flag=True)
@click.option("--gather", is_flag=True)
@click.option("--interactive", is_flag=True)
@click.option("--debug", is_flag=True)
@click.option("--plot", is_flag=True)
@click.option("--show-fig", is_flag=True)
@click.option("--save-fig", required=False, type=str, default=None)
@click.option("--save-fig-dir", type=click.Path(exists=True, file_okay=False, dir_okay=True), default=None)
@click.option("--font-size", type=float, default=14)
@click.option("--fig-format", required=False, type=str, default='svg')
@click.option("--dry-run", is_flag=True)
def main(config, config_prefix, coords_path, distances_namespace, output_path, arena_id, populations, use_noise_gen, io_size,
         chunk_size, value_chunk_size, cache_size, write_size, verbose, gather, interactive, debug, plot, show_fig,
         save_fig, save_fig_dir, font_size, fig_format, dry_run):
    """

    :param config: str (.yaml file name)
    :param config_prefix: str (path to dir)
    :param coords_path: str (path to file)
    :param distances_namespace: str
    :param output_path: str
    :param arena_id: str
    :param populations: tuple of str
    :param io_size: int
    :param chunk_size: int
    :param value_chunk_size: int
    :param cache_size: int
    :param write_size: int
    :param verbose: bool
    :param gather: bool; whether to gather population attributes to rank 0 for interactive analysis or plotting
    :param interactive: bool
    :param debug: bool
    :param plot: bool
    :param show_fig: bool
    :param save_fig: str (base file name)
    :param save_fig_dir:  str (path to dir)
    :param font_size: float
    :param fig_format: str
    :param dry_run: bool
    """
    comm = MPI.COMM_WORLD
    rank = comm.rank

    config_logging(verbose)

    env = Env(comm=comm, config_file=config, config_prefix=config_prefix, template_paths=None)
    if io_size == -1:
        io_size = comm.size
    if rank == 0:
        logger.info(f'{comm.size} ranks have been allocated')

    if save_fig is not None:
        plot = True

    if plot:
        import matplotlib.pyplot as plt
        from dentate.plot import plot_2D_rate_map, default_fig_options, save_figure, clean_axes, close_figure

        fig_options = copy.copy(default_fig_options)
        fig_options.saveFigDir = save_fig_dir
        fig_options.fontSize = font_size
        fig_options.figFormat = fig_format
        fig_options.showFig = show_fig

    if save_fig is not None:
        save_fig = '%s %s' % (save_fig, arena_id)
        fig_options.saveFig = save_fig

    if not dry_run and rank == 0:
        if output_path is None:
            raise RuntimeError('generate_input_selectivity_features: missing output_path')
        if not os.path.isfile(output_path):
            input_file = h5py.File(coords_path, 'r')
            output_file = h5py.File(output_path, 'w')
            input_file.copy('/H5Types', output_file)
            input_file.close()
            output_file.close()
    comm.barrier()
    population_ranges = read_population_ranges(coords_path, comm)[0]

    if len(populations) == 0:
        populations = sorted(population_ranges.keys())

    reference_u_arc_distance_bounds_dict = {}
    if rank == 0:
        for population in sorted(populations):
            if population not in population_ranges:
                raise RuntimeError('generate_input_selectivity_features: specified population: %s not found in '
                                   'provided coords_path: %s' % (population, coords_path))
            if population not in env.stimulus_config['Selectivity Type Probabilities']:
                raise RuntimeError('generate_input_selectivity_features: selectivity type not specified for '
                                   'population: %s' % population)
            with h5py.File(coords_path, 'r') as coords_f:
                pop_size = population_ranges[population][1]
                unique_gid_count = len(set(
                    coords_f['Populations'][population][distances_namespace]['U Distance']['Cell Index'][:]))
                if pop_size != unique_gid_count:
                    raise RuntimeError('generate_input_selectivity_features: only %i/%i unique cell indexes found '
                                       'for specified population: %s in provided coords_path: %s' %
                                       (unique_gid_count, pop_size, population, coords_path))
                try:
                    reference_u_arc_distance_bounds_dict[population] = \
                      coords_f['Populations'][population][distances_namespace].attrs['Reference U Min'], \
                      coords_f['Populations'][population][distances_namespace].attrs['Reference U Max']
                except Exception:
                    raise RuntimeError('generate_input_selectivity_features: problem locating attributes '
                                       'containing reference bounds in namespace: %s for population: %s from '
                                       'coords_path: %s' % (distances_namespace, population, coords_path))
    comm.barrier()
    reference_u_arc_distance_bounds_dict = comm.bcast(reference_u_arc_distance_bounds_dict, root=0)

    selectivity_type_names = dict([ (val, key) for (key, val) in viewitems(env.selectivity_types) ])
    selectivity_type_namespaces = dict()
    for this_selectivity_type in selectivity_type_names:
        this_selectivity_type_name = selectivity_type_names[this_selectivity_type]
        chars = list(this_selectivity_type_name)
        chars[0] = chars[0].upper()
        selectivity_type_namespaces[this_selectivity_type_name] = ''.join(chars) + ' Selectivity %s' % arena_id

    if arena_id not in env.stimulus_config['Arena']:
        raise RuntimeError(f'Arena with ID: {arena_id} not specified by configuration at file path: {config_prefix}/{config}')
    arena = env.stimulus_config['Arena'][arena_id]
    arena_x_mesh, arena_y_mesh = None, None
    if rank == 0:
        arena_x_mesh, arena_y_mesh = \
             get_2D_arena_spatial_mesh(arena=arena, spatial_resolution=env.stimulus_config['Spatial Resolution'])
    arena_x_mesh = comm.bcast(arena_x_mesh, root=0)
    arena_y_mesh = comm.bcast(arena_y_mesh, root=0)

    local_random = np.random.RandomState()
    selectivity_seed_offset = int(env.model_config['Random Seeds']['Input Selectivity'])
    local_random.seed(selectivity_seed_offset - 1)

    selectivity_config = InputSelectivityConfig(env.stimulus_config, local_random)
    if plot and rank == 0:
        selectivity_config.plot_module_probabilities(**fig_options())

    if (debug or interactive) and rank == 0:
        context.update(dict(locals()))

    pop_norm_distances = {}
    rate_map_sum = {}
    x0_dict = {}
    y0_dict = {}
    write_every = max(1, int(math.floor(write_size / comm.size)))
    for population in sorted(populations):
        if rank == 0:
            logger.info(f'Generating input selectivity features for population {population}...')

        reference_u_arc_distance_bounds = reference_u_arc_distance_bounds_dict[population]

        modular = True
        if population in env.stimulus_config['Non-modular Place Selectivity Populations']:
            modular = False

        noise_gen_dict = None
        if use_noise_gen:
            noise_gen_dict = {}
            if modular:
                for module_id in range(env.stimulus_config['Number Modules']):
                    extent_x, extent_y = get_2D_arena_extents(arena)
                    margin = round(selectivity_config.place_module_field_widths[module_id] / 2.)
                    arena_x_bounds, arena_y_bounds = get_2D_arena_bounds(arena, margin=margin)
                    noise_gen = MPINoiseGenerator(comm=comm, bounds=(arena_x_bounds, arena_y_bounds),
                                                  tile_rank=comm.rank, bin_size=0.5, mask_fraction=0.99,
                                                  seed=selectivity_seed_offset)
                    noise_gen_dict[module_id] = noise_gen
            else:
                margin = round(np.mean(selectivity_config.place_module_field_widths) / 2.)
                arena_x_bounds, arena_y_bounds = get_2D_arena_bounds(arena, margin=margin)
                noise_gen_dict[-1] = MPINoiseGenerator(comm=comm, bounds=(arena_x_bounds, arena_y_bounds),
                                                       tile_rank=comm.rank, bin_size=0.5, mask_fraction=0.99,
                                                       seed=selectivity_seed_offset)

        
        this_pop_norm_distances = {}
        this_rate_map_sum = defaultdict(lambda: np.zeros_like(arena_x_mesh))
        this_x0_list = []
        this_y0_list = []
        start_time = time.time()
        gid_count = defaultdict(lambda: 0)
        distances_attr_gen = NeuroH5CellAttrGen(coords_path, population, namespace=distances_namespace,
                                                comm=comm, io_size=io_size, cache_size=cache_size)

        selectivity_attr_dict = dict((key, dict()) for key in env.selectivity_types)
        for iter_count, (gid, distances_attr_dict) in enumerate(distances_attr_gen):
            req = comm.Ibarrier()
            if gid is None:
                if noise_gen_dict is not None:
                    all_module_ids = [-1]
                    if modular:
                        all_module_ids = comm.allreduce(set([]), op=mpi_op_set_union)
                    for module_id in all_module_ids:
                        this_noise_gen = noise_gen_dict[module_id]
                        global_num_fields = this_noise_gen.sync(0)
                        for i in range(global_num_fields):
                            this_noise_gen.add(np.empty( shape=(0, 0), dtype=np.float32 ), None)
            else:
                if rank == 0:
                    logger.info(f'Rank {rank} generating selectivity features for gid {gid}...')
                u_arc_distance = distances_attr_dict['U Distance'][0]
                v_arc_distance = distances_attr_dict['V Distance'][0]
                norm_u_arc_distance = ((u_arc_distance - reference_u_arc_distance_bounds[0]) /
                                       (reference_u_arc_distance_bounds[1] - reference_u_arc_distance_bounds[0]))

                this_pop_norm_distances[gid] = norm_u_arc_distance

                this_selectivity_type_name, this_selectivity_attr_dict = \
                 generate_input_selectivity_features(env, population, arena,
                                                     arena_x_mesh, arena_y_mesh,
                                                     gid, (norm_u_arc_distance, v_arc_distance),
                                                     selectivity_config, selectivity_type_names,
                                                     selectivity_type_namespaces,
                                                     noise_gen_dict=noise_gen_dict,
                                                     rate_map_sum=this_rate_map_sum,
                                                     debug= (debug_callback, context) if debug else False)
                if 'X Offset' in this_selectivity_attr_dict:
                    this_x0_list.append(this_selectivity_attr_dict['X Offset'])
                    this_y0_list.append(this_selectivity_attr_dict['Y Offset'])
                selectivity_attr_dict[this_selectivity_type_name][gid] = this_selectivity_attr_dict
                gid_count[this_selectivity_type_name] += 1
            if noise_gen_dict is not None:
                for m in noise_gen_dict:
                    noise_gen_dict[m].tile_rank = (noise_gen_dict[m].tile_rank + 1) % comm.size
            req.wait()

            if (iter_count > 0 and iter_count % write_every == 0) or (debug and iter_count == 10):
                total_gid_count = 0
                gid_count_dict = dict(gid_count.items())
                req = comm.Ibarrier() 
                selectivity_gid_count = comm.reduce(gid_count_dict, root=0, op=mpi_op_merge_count_dict)
                req.wait()
                if rank == 0:
                    for selectivity_type_name in selectivity_gid_count:
                        total_gid_count += selectivity_gid_count[selectivity_type_name]
                    for selectivity_type_name in selectivity_gid_count:
                        logger.info('generated selectivity features for %i/%i %s %s cells in %.2f s' %
                                    (selectivity_gid_count[selectivity_type_name], total_gid_count, population, selectivity_type_name,
                                    (time.time() - start_time)))

                if not dry_run:
                    for selectivity_type_name in sorted(selectivity_attr_dict.keys()):
                       req = comm.Ibarrier() 
                       if rank == 0:
                            logger.info(f'writing selectivity features for {population} [{selectivity_type_name}]...')
                       selectivity_type_namespace = selectivity_type_namespaces[selectivity_type_name]
                       append_cell_attributes(output_path, population, selectivity_attr_dict[selectivity_type_name],
                                               namespace=selectivity_type_namespace, comm=comm, io_size=io_size,
                                               chunk_size=chunk_size, value_chunk_size=value_chunk_size)
                       req.wait()
                    del selectivity_attr_dict
                    selectivity_attr_dict = dict((key, dict()) for key in env.selectivity_types)
                    gc.collect()


            if debug and iter_count >= 10:
                break


            
        pop_norm_distances[population] = this_pop_norm_distances
        rate_map_sum[population] = dict(this_rate_map_sum)
        if len(this_x0_list) > 0:
            x0_dict[population] = np.concatenate(this_x0_list, axis=None)
            y0_dict[population] = np.concatenate(this_y0_list, axis=None)
            
        
        total_gid_count = 0
        gid_count_dict = dict(gid_count.items())
        req = comm.Ibarrier()
        selectivity_gid_count = comm.reduce(gid_count_dict, root=0, op=mpi_op_merge_count_dict)
        req.wait()
        
        if rank == 0:
            for selectivity_type_name in selectivity_gid_count:
                total_gid_count += selectivity_gid_count[selectivity_type_name]
            for selectivity_type_name in selectivity_gid_count:
                logger.info('generated selectivity features for %i/%i %s %s cells in %.2f s' %
                            (selectivity_gid_count[selectivity_type_name], total_gid_count, population, selectivity_type_name,
                             (time.time() - start_time)))

        if not dry_run:
            for selectivity_type_name in sorted(selectivity_attr_dict.keys()):
                req = comm.Ibarrier()
                if rank == 0:
                    logger.info(f'writing selectivity features for {population} [{selectivity_type_name}]...')
                selectivity_type_namespace = selectivity_type_namespaces[selectivity_type_name]
                append_cell_attributes(output_path, population, selectivity_attr_dict[selectivity_type_name],
                                       namespace=selectivity_type_namespace, comm=comm, io_size=io_size,
                                       chunk_size=chunk_size, value_chunk_size=value_chunk_size)
                req.wait()
            del selectivity_attr_dict
            gc.collect()
        req = comm.Ibarrier()
        req.wait()
                
    if gather:
        merged_pop_norm_distances = {}
        for population in sorted(populations):
            merged_pop_norm_distances[population] = \
              comm.reduce(pop_norm_distances[population], root=0,
                          op=mpi_op_merge_dict)
        merged_rate_map_sum = comm.reduce(rate_map_sum, root=0, op=mpi_op_merge_rate_map_dict)
        merged_x0 = comm.reduce(x0_dict, root=0, op=mpi_op_concatenate_ndarray_dict)
        merged_y0 = comm.reduce(y0_dict, root=0, op=mpi_op_concatenate_ndarray_dict)
        if rank == 0:
            if plot:
                for population in merged_pop_norm_distances:
                    norm_distance_values = np.asarray(list(merged_pop_norm_distances[population].values()))
                    hist, edges = np.histogram(norm_distance_values, bins=100)
                    fig, axes = plt.subplots(1)
                    axes.plot(edges[1:], hist)
                    axes.set_title(f'Population: {population}')
                    axes.set_xlabel('Normalized cell position')
                    axes.set_ylabel('Cell count')
                    clean_axes(axes)
                    if save_fig is not None:
                        save_figure(f'{save_fig} {population} normalized distances histogram',
                                    fig=fig, **fig_options())
                    if fig_options.showFig:
                        fig.show()
                    close_figure(fig)
                for population in merged_rate_map_sum:
                    for selectivity_type_name in merged_rate_map_sum[population]:
                        fig_title = f'{population} {this_selectivity_type_name} summed rate maps'
                        if save_fig is not None:
                            fig_options.saveFig = f'{save_fig} {fig_title}'
                        plot_2D_rate_map(x=arena_x_mesh, y=arena_y_mesh,
                                         rate_map=merged_rate_map_sum[population][selectivity_type_name],
                                         title=f'Summed rate maps\n{population} {selectivity_type_name} cells',
                                               **fig_options())
                for population in merged_x0:
                    fig_title = f'{population} field offsets'
                    if save_fig is not None:
                        fig_options.saveFig = f'{save_fig} {fig_title}'
                    x0 = merged_x0[population]
                    y0 = merged_y0[population]
                    fig, axes = plt.subplots(1)
                    axes.scatter(x0, y0)
                    if save_fig is not None:
                        save_figure(f'{save_fig} {fig_title}', fig=fig, **fig_options())
                    if fig_options.showFig:
                        fig.show()
                    close_figure(fig)

    if interactive and rank == 0:
        context.update(locals())


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda x: os.path.basename(x) == os.path.basename(__file__), sys.argv) + 1):],
         standalone_mode=False)
