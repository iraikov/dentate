from specify_cells import *
from mpi4py import MPI
from neurotrees.io import NeurotreeGen
from neurotrees.io import append_cell_attributes
# import mkl

# mkl.set_num_threads(1)

comm = MPI.COMM_WORLD
rank = comm.rank

if rank == 0:
    print '%i ranks have been allocated' % comm.size
sys.stdout.flush()

neurotrees_dir = '../morphologies/'
# forest_file = '122016_DGC_forest_test_copy.h5'
# neurotrees_dir = os.environ['PI_SCRATCH']+'/DGC_forest/hdf5/'
# neurotrees_dir = os.environ['PI_HOME']+'/'
# forest_file = 'DGC_forest_full.h5'
# forest_file = 'DGC_forest_syns_012717.h5'
forest_file = 'DGC_forest_syn_locs_test_041217.h5'

# forest_file = 'DGC_forest_test.h5'

population = 'GC'
g = NeurotreeGen(MPI._addressof(comm), neurotrees_dir+forest_file, population, io_size=comm.size)

sys.stdout.flush()

count = 0

start_time = time.time()
for gid, morph_dict in g:
    local_time = time.time()
    # mismatched_section_dict = {}
    synapse_dict = {}
    cell = DG_GC(neurotree_dict=morph_dict, gid=gid, full_spines=False)
    # this_mismatched_sections = cell.get_mismatched_neurotree_sections()
    # if this_mismatched_sections is not None:
    #    mismatched_section_dict[gid] = this_mismatched_sections
    synapse_dict[gid] = cell.export_neurotree_synapse_attributes()
    del cell
    append_cell_attributes(MPI._addressof(comm), neurotrees_dir+forest_file, population, synapse_dict,
                           namespace='Synapse_Attributes', value_chunk_size=48000)
    print 'Rank %i took %i s to compute syn_locs for %s gid: %i' % (rank, time.time() - local_time, population, gid)
    count += 1
    sys.stdout.flush()
    del synapse_dict
    gc.collect()

# len_mismatched_section_dict_fragments = comm.gather(len(mismatched_section_dict), root=0)
# len_GID_fragments = comm.gather(len(GID), root=0)
global_count = comm.gather(count, root=0)
if rank == 0:
    print '%i ranks took %i s to compute synapse locations for %i morphologies' % (comm.size,
                                                                                   time.time() - start_time,
                                                                                   np.sum(global_count))
    # print '%i morphologies have mismatched section indexes' % np.sum(len_mismatched_section_dict_fragments)
