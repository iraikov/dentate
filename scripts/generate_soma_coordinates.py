##
## Generate soma coordinates within layer-specific volume.
##

import sys, itertools, os.path, math, random, click, logging
from mpi4py import MPI
import numpy as np
from neuroh5.io import read_population_ranges, append_cell_attributes
import h5py
from dentate.utils import *
from dentate.env import Env
from dentate.geometry import make_volume, DG_volume, make_uvl_distance
from dentate.alphavol import alpha_shape
import dlib, rbf
from rbf.pde.nodes import min_energy_nodes
from rbf.pde.geometry import contains

script_name = os.path.basename(__file__)
logger = get_script_logger(script_name)

def mpi_excepthook(type, value, traceback):
    """

    :param type:
    :param value:
    :param traceback:
    :return:
    """
    sys_excepthook(type, value, traceback)
    if MPI.COMM_WORLD.size > 1:
        MPI.COMM_WORLD.Abort(1)


sys_excepthook = sys.excepthook
sys.excepthook = mpi_excepthook

def random_subset( iterator, K ):
    result = []
    N = 0

    for item in iterator:
        N += 1
        if len( result ) < K:
            result.append( item )
        else:
            s = int(random.random() * N)
            if s < K:
                result[ s ] = item

    return result

def uvl_in_bounds(uvl_coords, pop_min_extent, pop_max_extent):
    result = (uvl_coords[0] <= pop_max_extent[0]) and \
      (uvl_coords[0] > pop_min_extent[0]) and \
      (uvl_coords[1] < pop_max_extent[1]) and \
      (uvl_coords[1] > pop_min_extent[1]) and \
      (uvl_coords[2] < pop_max_extent[2]) and \
      (uvl_coords[2] > pop_min_extent[2])
    return result


@click.command()
@click.option("--config", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--types-path", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--template-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option("--output-path", required=True, type=click.Path(exists=False, file_okay=True, dir_okay=False))
@click.option("--output-namespace", type=str, default='Generated Coordinates')
@click.option("--populations", '-i', type=str, multiple=True)
@click.option("--resolution", type=(int,int,int), default=(30,30,10))
@click.option("--alpha-radius", type=float, default=120.)
@click.option("--nodeiter", type=int, default=10)
@click.option("--optiter", type=int, default=200)
@click.option("--io-size", type=int, default=-1)
@click.option("--chunk-size", type=int, default=1000)
@click.option("--value-chunk-size", type=int, default=1000)
@click.option("--verbose", '-v', type=bool, default=False, is_flag=True)
def main(config, types_path, template_path, output_path, output_namespace, populations, resolution, alpha_radius, nodeiter, optiter, io_size, chunk_size, value_chunk_size, verbose):

    config_logging(verbose)
    logger = get_script_logger(script_name)

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    
    if io_size == -1:
        io_size = comm.size
    if rank == 0:
        logger.info('%i ranks have been allocated' % comm.size)
    sys.stdout.flush()

    if rank==0:
        if not os.path.isfile(output_path):
            input_file  = h5py.File(types_path,'r')
            output_file = h5py.File(output_path,'w')
            input_file.copy('/H5Types',output_file)
            input_file.close()
            output_file.close()
    comm.barrier()

    env = Env(comm=comm, config_file=config)

    layer_min_extents = env.geometry['Parametric Surface']['Minimum Extent']
    layer_max_extents = env.geometry['Parametric Surface']['Maximum Extent']
    rotate = env.geometry['Parametric Surface']['Rotation']

    random_seed = int(env.modelConfig['Random Seeds']['Soma Locations'])
    
    population_ranges = read_population_ranges(output_path, comm)[0]

    for population in populations:

        if verbose and (rank == 0):
            logger.info( 'population: %s' % population )

        (population_start, population_count) = population_ranges[population]

        pop_min_extent = env.geometry['Cell Layers']['Minimum Extent'][population]
        pop_max_extent = env.geometry['Cell Layers']['Maximum Extent'][population]

        if verbose and (rank == 0):
            logger.info('min extent: %f %f %f' % (pop_min_extent[0],pop_min_extent[1],pop_min_extent[2]))
            logger.info('max extent: %f %f %f' % (pop_max_extent[0],pop_max_extent[1],pop_max_extent[2]))

        xyz_coords = None
        xyz_coords_interp = None
        uvl_coords_interp = None
        if rank == 0:
            if verbose:
                logger.info("Constructing volume...")
                
            vol = make_volume((pop_min_extent[0], pop_max_extent[0]), \
                              (pop_min_extent[1], pop_max_extent[1]), \
                              (pop_min_extent[2], pop_max_extent[2]), \
                              rotate=rotate, resolution=resolution)
            
            if verbose:
                logger.info("Constructing volume triangulation...")
            tri = vol.create_triangulation()

            if verbose:
                logger.info("Constructing alpha shape...")
            alpha = alpha_shape([], alpha_radius, tri=tri)
    
            vert = alpha.points
            smp  = np.asarray(alpha.bounds, dtype=np.int64)

            N = int(population_count*2) # total number of nodes
            node_count = 0

            if verbose:
                logger.info("Generating %i nodes..." % N)

            if verbose:
                rbf_logger = logging.Logger.manager.loggerDict['rbf.pde.nodes']
                rbf_logger.setLevel(logging.DEBUG)

            while node_count < population_count:
                # create N quasi-uniformly distributed nodes
                out = min_energy_nodes(N,(vert,smp),iterations=nodeiter)
                nodes = out[0]
        
                # remove nodes outside of the domain
                in_nodes = nodes[contains(nodes,vert,smp)]
                
                node_count = len(in_nodes)
                N = int(1.5*N)
            
                if verbose:
                    logger.info("%i interior nodes out of %i nodes generated" % (node_count, len(nodes)))

            if verbose:
                logger.info("Inverse interpolation of %i nodes..." % node_count)

            xyz_coords = in_nodes.reshape(-1,3)
            uvl_coords_interp = vol.inverse(xyz_coords)
            xyz_coords_interp = vol(uvl_coords_interp[:,0],uvl_coords_interp[:,1],uvl_coords_interp[:,2],mesh=False).reshape(3,-1).T

            if verbose:
                logger.info("Broadcasting generated nodes...")

            
        xyz_coords = comm.bcast(xyz_coords, root=0)
        xyz_coords_interp = comm.bcast(xyz_coords_interp, root=0)
        uvl_coords_interp = comm.bcast(uvl_coords_interp, root=0)

        coords = []
        coords_dict = {}
        xyz_error = np.asarray([0.0, 0.0, 0.0])

        if verbose:
            if rank == 0:
                logger.info("Computing UVL coordinates...")

        for i in range(0,xyz_coords.shape[0]):

            coord_ind = i
            if i % size == rank:

                xyz_error_interp  = np.abs(np.subtract(xyz_coords[coord_ind,:], xyz_coords_interp[coord_ind,:]))

                f_uvl_distance = make_uvl_distance(xyz_coords[coord_ind,:],rotate=rotate)
                uvl_coords_opt,dist = dlib.find_min_global(f_uvl_distance, pop_min_extent, pop_max_extent, optiter)
                xyz_coords_opt = DG_volume(uvl_coords_opt[0], uvl_coords_opt[1], uvl_coords_opt[2], rotate=rotate)[0]
                xyz_error_opt  = np.abs(np.subtract(xyz_coords[coord_ind,:], xyz_coords_opt))

                
                if uvl_in_bounds(uvl_coords_opt, pop_min_extent, pop_max_extent) and \
                   np.all (np.less (xyz_error_opt, xyz_error_interp)):
                    uvl_coords  = uvl_coords_opt
                    xyz_coords1 = xyz_coords_opt
                elif uvl_in_bounds(uvl_coords_interp[coord_ind,:], pop_min_extent, pop_max_extent):
                    uvl_coords  = uvl_coords_interp[coord_ind,:].ravel()
                    xyz_coords1 = xyz_coords_interp[coord_ind,:].ravel()
                else:
                    uvl_coords = None
                    xyz_coords1 = None

                if uvl_coords is not None:

                    xyz_error   = np.add(xyz_error, np.abs(np.subtract(xyz_coords[coord_ind,:], xyz_coords1)))

                    if verbose:
                        logger.info('Rank %i: cell %i: %f %f %f' % (rank, i, uvl_coords[0], uvl_coords[1], uvl_coords[2]))

                    coords.append((xyz_coords1[0],xyz_coords1[1],xyz_coords1[2],\
                                       uvl_coords[0],uvl_coords[1],uvl_coords[2]))
                                       
        
        total_xyz_error = np.zeros((3,))
        comm.Allreduce(xyz_error, total_xyz_error, op=MPI.SUM)

        coords_count = 0
        coords_count = np.sum(np.asarray(comm.allgather(len(coords))))

        if verbose:
            if rank == 0:
                logger.info('Total %i coordinates generated' % coords_count)

        mean_xyz_error = np.asarray([old_div(total_xyz_error[0], coords_count), \
                                     old_div(total_xyz_error[1], coords_count), \
                                     old_div(total_xyz_error[2], coords_count)])

        
        if verbose:
            if rank == 0:
                logger.info("mean XYZ error: %f %f %f " % (mean_xyz_error[0], mean_xyz_error[1], mean_xyz_error[2]))

        if rank == 0:
            color = 1
        else:
            color = 0

        ## comm0 includes only rank 0
        comm0 = comm.Split(color, 0)

        coords_lst = comm.gather(coords, root=0)
        if rank == 0:
            all_coords = []
            for sublist in coords_lst:
                for item in sublist:
                    all_coords.append(item)


            if coords_count < population_count:
                logger.warning("Generating additional %i coordinates " % (population_count - len(all_coords)))

                safety = 0.01
                sampled_coords = all_coords
                for i in range(population_count - len(all_coords)):
                    coord_u = np.random.uniform(pop_min_extent[0] + safety, pop_max_extent[0] - safety)
                    coord_v = np.random.uniform(pop_min_extent[1] + safety, pop_max_extent[1] - safety)
                    coord_l = np.random.uniform(pop_min_extent[2] + safety, pop_max_extent[2] - safety)
                    xyz_coords = DG_volume(coord_u, coord_v, coord_l, rotate=rotate).ravel()
                    sampled_coords.append((xyz_coords[0],xyz_coords[1],xyz_coords[2],\
                                           coord_u, coord_v, coord_l))
            else:
                sampled_coords = random_subset(all_coords, int(population_count))

            
            sampled_coords.sort(key=lambda coord: coord[3]) ## sort on U coordinate
            coords_dict = { population_start+i :  { 'X Coordinate': np.asarray([x_coord],dtype=np.float32),
                                    'Y Coordinate': np.asarray([y_coord],dtype=np.float32),
                                    'Z Coordinate': np.asarray([z_coord],dtype=np.float32),
                                    'U Coordinate': np.asarray([u_coord],dtype=np.float32),
                                    'V Coordinate': np.asarray([v_coord],dtype=np.float32),
                                    'L Coordinate': np.asarray([l_coord],dtype=np.float32) }
                            for (i,(x_coord,y_coord,z_coord,u_coord,v_coord,l_coord)) in enumerate(sampled_coords) }

            append_cell_attributes(output_path, population, coords_dict,
                                    namespace=output_namespace,
                                    io_size=io_size, chunk_size=chunk_size,
                                    value_chunk_size=value_chunk_size,comm=comm0)

        comm.Barrier()
        

if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda x: os.path.basename(x) == os.path.basename(__file__), sys.argv)+1):])
