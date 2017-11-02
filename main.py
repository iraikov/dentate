##
##  Dentate Gyrus model initialization script
##

import sys, os
import os.path
import click
import itertools
from collections import defaultdict
from datetime import datetime
import numpy as np
from mpi4py import MPI # Must come before importing NEURON
from neuron import h
from neuroh5.io import read_projection_names, scatter_read_graph, bcast_graph, scatter_read_trees, scatter_read_cell_attributes, write_cell_attributes
import h5py
from env import Env
import lpt, utils, synapses, cells

## Estimate cell complexity. Code by Michael Hines from the discussion thread
## https://www.neuron.yale.edu/phpBB/viewtopic.php?f=31&t=3628
def cx(env):
  rank   = int(env.pc.id())
  lb = h.LoadBalance()
  if os.path.isfile("mcomplex.dat"):
    lb.read_mcomplex() 
  cxvec = h.Vector(len(env.gidlist))
  for i, gid in enumerate(env.gidlist):
    cxvec.x[i] = lb.cell_complexity(env.pc.gid2cell(gid))
  env.cxvec = cxvec
  return cxvec

# for given cxvec on each rank what is the fractional load balance.
def ld_bal(env):
  rank   = int(env.pc.id())
  nhosts = int(env.pc.nhost())
  cxvec  = env.cxvec
  sum_cx = sum(cxvec)
  max_sum_cx = env.pc.allreduce(sum_cx, 2)
  sum_cx = env.pc.allreduce(sum_cx, 1)
  if rank == 0:
    print ("*** expected load balance %.2f" % (sum_cx / nhosts / max_sum_cx))

# Each rank has gidvec, cxvec: gather everything to rank 0, do lpt
# algorithm and write to a balance file.
def lpt_bal(env):
  rank   = int(env.pc.id())
  nhosts = int(env.pc.nhost())

  cxvec  = env.cxvec
  gidvec = env.gidlist
  #gather gidvec, cxvec to rank 0
  src    = [None]*nhosts
  src[0] = zip(cxvec.to_python(), gidvec)
  dest   = env.pc.py_alltoall(src)
  del src

  if rank == 0:
    lb = h.LoadBalance()
    allpairs = sum(dest,[])
    del dest
    parts = lpt.lpt(allpairs, nhosts)
    lpt.statistics(parts)
    part_rank = 0
    with open('parts.%d' % nhosts, 'w') as fp:
      for part in parts:
        for x in part[1]:
          fp.write('%d %d\n' % (x[1],part_rank))
        part_rank = part_rank+1

        
def mkspikeout (env, spikeout_filename):
    datasetPath     = os.path.join(env.datasetPrefix,env.datasetName)
    forestFilePath  = os.path.join(datasetPath,env.modelConfig['Cell Data'])
    forestFile      = h5py.File(forestFilePath)
    spikeoutFile    = h5py.File(spikeout_filename)
    forestFile.copy('/H5Types',spikeoutFile)
    forestFile.close()
    spikeoutFile.close()

    
def spikeout (env, output_path, t_vec, id_vec):
    binlst  = []
    typelst = env.celltypes.keys()
    for k in typelst:
        binlst.append(env.celltypes[k]['start'])
        
    binvect  = np.array(binlst)
    sort_idx = np.argsort(binvect,axis=0)
    bins     = binvect[sort_idx][1:]
    types    = [ typelst[i] for i in sort_idx ]
    inds     = np.digitize(id_vec, bins)

    if not str(env.resultsId):
        namespace_id = "Spike Events" 
    else:
        namespace_id = "Spike Events %s" % str(env.resultsId)
    
    for i in range(0,len(types)):
        if i > 0:
            start = bins[i-1]
        else:
            start = 0
        spkdict  = {}
        sinds    = np.where(inds == i)
        if len(sinds) > 0:
            ids      = id_vec[sinds]
            ts       = t_vec[sinds]
            for j in range(0,len(ids)):
                id = ids[j] - start
                t  = ts[j]
                if spkdict.has_key(id):
                    spkdict[id]['t'].append(t)
                else:
                    spkdict[id]= {'t': [t]}
            for j in spkdict.keys():
                spkdict[j]['t'] = np.array(spkdict[j]['t'])
        pop_name = types[i]
        write_cell_attributes(env.comm, output_path, pop_name, spkdict, namespace=namespace_id)
        

def connectcells(env):
    datasetPath = os.path.join(env.datasetPrefix,env.datasetName)
    connectivityFilePath = os.path.join(datasetPath,env.modelConfig['Connection Data'])
    forestFilePath = os.path.join(datasetPath,env.modelConfig['Cell Data'])

    if env.verbose:
      if env.pc.id() == 0:
        print '*** Connectivity file path is %s' % connectivityFilePath

    prj_dict = defaultdict(list)
    for (src,dst) in read_projection_names(env.comm,connectivityFilePath):
      prj_dict[dst].append(src)

    if env.verbose:
      if env.pc.id() == 0:
        print '*** Reading projections: ', prj_dict.items()
      
    for (postsyn_name, presyn_names) in prj_dict.iteritems():

      synapse_config = env.celltypes[postsyn_name]['synapses']
      if synapse_config.has_key('spines'):
        spines = synapse_config['spines']
      else:
        spines = False

      if synapse_config.has_key('unique'):
        unique = synapse_config['unique']
      else:
        unique = False
      
      if env.verbose:
          if env.pc.id() == 0:
            print '*** Reading synapse attributes of population %s' % (postsyn_name)

      if env.nodeRanks is None:
        cell_attributes_dict = scatter_read_cell_attributes(env.comm, forestFilePath, postsyn_name, 
                                                            namespaces=['Synapse Attributes'],
                                                            io_size=env.IOsize)
      else:
        cell_attributes_dict = scatter_read_cell_attributes(env.comm, forestFilePath, postsyn_name, 
                                                            namespaces=['Synapse Attributes'],
                                                            node_rank_map=env.nodeRanks,
                                                            io_size=env.IOsize)
      cell_synapses_dict = { k : v for (k,v) in cell_attributes_dict['Synapse Attributes'] }
      del cell_attributes_dict

      for presyn_name in presyn_names:

        edge_count = 0

        if env.verbose:
          if env.pc.id() == 0:
            print '*** Connecting %s -> %s' % (presyn_name, postsyn_name)

        if env.nodeRanks is None:
          (graph, a) = scatter_read_graph(env.comm,connectivityFilePath,io_size=env.IOsize,
                                          projections=[(presyn_name, postsyn_name)],
                                          namespaces=['Synapses','Connections'])
        else:
          (graph, a) = scatter_read_graph(env.comm,connectivityFilePath,io_size=env.IOsize,
                                          node_rank_map=env.nodeRanks,
                                          projections=[(presyn_name, postsyn_name)],
                                          namespaces=['Synapses','Connections'])
        
        edge_iter = graph[postsyn_name][presyn_name]


        connection_dict = env.connection_generator[postsyn_name][presyn_name].connection_properties
        kinetics_dict = env.connection_generator[postsyn_name][presyn_name].synapse_kinetics

        syn_id_attr_index = a[postsyn_name][presyn_name]['Synapses']['syn_id']
        distance_attr_index = a[postsyn_name][presyn_name]['Connections']['distance']

        for (postsyn_gid, edges) in edge_iter:

          postsyn_cell   = env.pc.gid2cell(postsyn_gid)
          cell_syn_dict  = cell_synapses_dict[postsyn_gid]

          presyn_gids    = edges[0]
          edge_syn_ids   = edges[1]['Synapses'][syn_id_attr_index]
          edge_dists     = edges[1]['Connections'][distance_attr_index]

          cell_syn_types    = cell_syn_dict['syn_types']
          cell_swc_types    = cell_syn_dict['swc_types']
          cell_syn_locs     = cell_syn_dict['syn_locs']
          cell_syn_sections = cell_syn_dict['syn_secs']

          edge_syn_ps_dict  = synapses.mksyns(postsyn_gid,
                                              postsyn_cell,
                                              edge_syn_ids,
                                              cell_syn_types,
                                              cell_swc_types,
                                              cell_syn_locs,
                                              cell_syn_sections,
                                              kinetics_dict, env,
                                              add_synapse=synapses.add_unique_synapse if unique else synapses.add_shared_synapse,
                                              spines=spines)

          if env.verbose:
            if edge_count == 0:
              for sec in list(postsyn_cell.all):
                  h.psection(sec=sec)
          
          for (presyn_gid, edge_syn_id, distance) in itertools.izip(presyn_gids, edge_syn_ids, edge_dists):
            syn_ps_dict = edge_syn_ps_dict[edge_syn_id]
            for (syn_mech, syn_ps) in syn_ps_dict.iteritems():
              connection_syn_mech_config = connection_dict[syn_mech]
              weight = connection_syn_mech_config['weight']
              delay  = distance / connection_syn_mech_config['velocity']
              if type(weight) is float:
                h.nc_appendsyn (env.pc, h.nclist, presyn_gid, postsyn_gid, syn_ps, weight, delay)
              else:
                h.nc_appendsyn_wgtvector (env.pc, h.nclist, presyn_gid, postsyn_gid, syn_ps, weight, delay)

            edge_count += len(presyn_gids)

def connectgjs(env):
    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    datasetPath = os.path.join(env.datasetPrefix,env.datasetName)
    
    gapjunctions = env.gapjunctions
    if env.gapjunctionsFile is None:
        gapjunctionsFilePath = None
    else:
        gapjunctionsFilePath = os.path.join(datasetPath,env.gapjunctionsFile)

    if gapjunctions is not None:

        h('objref gjlist')
        h.gjlist = h.List()
        if env.verbose:
            if env.pc.id() == 0:
                print gapjunctions
        datasetPath = os.path.join(env.datasetPrefix,env.datasetName)
        (graph, a) = bcast_graph(env.comm,gapjunctionsFilePath,attributes=True)

        ggid = 2e6
        for name in gapjunctions.keys():
            if env.verbose:
                if env.pc.id() == 0:
                    print "*** Creating gap junctions %s" % name
            prj = graph[name]
            attrmap = a[name]
            weight_attr_idx = attrmap['Weight']+1
            dstbranch_attr_idx = attrmap['Destination Branch']+1
            dstsec_attr_idx = attrmap['Destination Section']+1
            srcbranch_attr_idx = attrmap['Source Branch']+1
            srcsec_attr_idx = attrmap['Source Section']+1
            for destination in sorted(prj.keys()):
                edges = prj[destination]
                sources      = edges[0]
                weights      = edges[weight_attr_idx]
                dstbranches  = edges[dstbranch_attr_idx]
                dstsecs      = edges[dstsec_attr_idx]
                srcbranches  = edges[srcbranch_attr_idx]
                srcsecs      = edges[srcsec_attr_idx]
                for i in range(0,len(sources)):
                    source    = sources[i]
                    srcbranch = srcbranches[i]
                    srcsec    = srcsecs[i]
                    dstbranch = dstbranches[i]
                    dstsec    = dstsecs[i]
                    weight    = weights[i]
                    if env.pc.gid_exists(source):
                        h.mkgap(env.pc, h.gjlist, source, srcbranch, srcsec, ggid, ggid+1, weight)
                    if env.pc.gid_exists(destination):
                        h.mkgap(env.pc, h.gjlist, destination, dstbranch, dstsec, ggid+1, ggid, weight)
                    ggid = ggid+2

            del graph[name]


def mkcells(env):

    h('objref templatePaths, templatePathValue')

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    datasetPath  = os.path.join(env.datasetPrefix, env.datasetName)

    h.templatePaths = h.List()
    for path in env.templatePaths:
        h.templatePathValue = h.Value(1,path)
        h.templatePaths.append(h.templatePathValue)
    popNames = env.celltypes.keys()
    popNames.sort()
    for popName in popNames:
        templateName = env.celltypes[popName]['template']
        h.find_template(env.pc, h.templatePaths, templateName)

    dataFilePath = os.path.join(datasetPath,env.modelConfig['Cell Data'])

    if rank == 0:
        print 'cell attributes: ', env.cellAttributeInfo

    for popName in popNames:

        if env.verbose:
            if env.pc.id() == 0:
                print "*** Creating population %s" % popName

        templateName = env.celltypes[popName]['template']
        templateClass = eval('h.%s' % templateName)

        if env.celltypes[popName].has_key('synapses'):
            synapses = env.celltypes[popName]['synapses']
        else:
            synapses = {}

        if env.cellAttributeInfo.has_key(popName) and env.cellAttributeInfo[popName].has_key('Trees'):
            if env.verbose:
                if env.pc.id() == 0:
                    print "*** Reading trees for population %s" % popName

            if env.nodeRanks is None:
                (trees, forestSize) = scatter_read_trees(env.comm, dataFilePath, popName, io_size=env.IOsize)
            else:
                (trees, forestSize) = scatter_read_trees(env.comm, dataFilePath, popName, io_size=env.IOsize,
                                                            node_rank_map=env.nodeRanks)
            if env.verbose:
                if env.pc.id() == 0:
                    print "*** Done reading trees for population %s" % popName

            h.numCells = 0
            i=0
            for (gid, tree) in trees:
                if env.verbose:
                    if env.pc.id() == 0:
                        print "*** Creating gid %i" % gid
            
                verboseflag = 0
                model_cell = cells.make_neurotree_cell(templateClass, neurotree_dict=tree, gid=gid, local_id=i, dataset_path=datasetPath)
                if env.verbose:
                    if (rank == 0) and (i == 0):
                        for sec in list(model_cell.all):
                            h.psection(sec=sec)
                env.gidlist.append(gid)
                env.cells.append(model_cell)
                env.pc.set_gid2node(gid, int(env.pc.id()))
                ## Tell the ParallelContext that this cell is a spike source
                ## for all other hosts. NetCon is temporary.
                nc = model_cell.connect2target(h.nil)
                env.pc.cell(gid, nc, 1)
                ## Record spikes of this cell
                env.pc.spike_record(gid, env.t_vec, env.id_vec)
                i = i+1
                h.numCells = h.numCells+1
            if env.verbose:
                if env.pc.id() == 0:
                    print "*** Created %i cells" % i
                    
        elif env.cellAttributeInfo.has_key(popName) and env.cellAttributeInfo[popName].has_key('Coordinates'):
            if env.verbose:
                if env.pc.id() == 0:
                    print "*** Reading coordinates for population %s" % popName
            
            if env.nodeRanks is None:
                cell_attributes_dict = scatter_read_cell_attributes(env.comm, dataFilePath, popName, 
                                                                    namespaces=['Coordinates'],
                                                                    io_size=env.IOsize)
            else:
                cell_attributes_dict = scatter_read_cell_attributes(env.comm, dataFilePath, popName, 
                                                                    namespaces=['Coordinates'],
                                                                    node_rank_map=env.nodeRanks,
                                                                    io_size=env.IOsize)
            if env.verbose:
                if env.pc.id() == 0:
                    print "*** Done reading coordinates for population %s" % popName
                    
            coords = cell_attributes_dict['Coordinates']

            h.numCells = 0
            i=0
            for (gid, _) in coords:
                if env.verbose:
                    if env.pc.id() == 0:
                        print "*** Creating gid %i" % gid
            
                verboseflag = 0
                model_cell = cells.make_cell(templateClass, gid=gid, local_id=i, dataset_path=datasetPath)
                env.gidlist.append(gid)
                env.cells.append(model_cell)
                env.pc.set_gid2node(gid, int(env.pc.id()))
                ## Tell the ParallelContext that this cell is a spike source
                ## for all other hosts. NetCon is temporary.
                nc = model_cell.connect2target(h.nil)
                env.pc.cell(gid, nc, 1)
                ## Record spikes of this cell
                env.pc.spike_record(gid, env.t_vec, env.id_vec)
                i = i+1
                h.numCells = h.numCells+1
            
             
def mkstim(env):

    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    datasetPath  = os.path.join(env.datasetPrefix, env.datasetName)
    
    inputFilePath = os.path.join(datasetPath,env.modelConfig['Cell Data'])

    popNames = env.celltypes.keys()
    popNames.sort()
    for popName in popNames:
        if env.celltypes[popName].has_key('vectorStimulus'):
            vecstim_namespace = env.celltypes[popName]['vectorStimulus']

            if env.nodeRanks is None:
              cell_attributes_dict = scatter_read_cell_attributes(env.comm, inputFilePath, popName, 
                                                                  namespaces=[vecstim_namespace],
                                                                  io_size=env.IOsize)
            else:
              cell_attributes_dict = scatter_read_cell_attributes(env.comm, inputFilePath, popName, 
                                                                  namespaces=[vecstim_namespace],
                                                                  node_rank_map=env.nodeRanks,
                                                                  io_size=env.IOsize)
            cell_vecstim = cell_attributes_dict[vecstim_namespace]
            for (gid, vecstim_dict) in cell_vecstim:
              if env.verbose:
                    if env.pc.id() == 0:
                        if len(vecstim_dict['spiketrain']) > 0:
                            print "*** Spike train for gid %i is of length %i (first spike at %g ms)" % (gid, len(vecstim_dict['spiketrain']),vecstim_dict['spiketrain'][0])
                        else:
                            print "*** Spike train for gid %i is of length %i" % (gid, len(vecstim_dict['spiketrain']))
                        
              cell = env.pc.gid2cell(gid)
              cell.play(h.Vector(vecstim_dict['spiketrain']))

            

            
def init(env):

    h.load_file("nrngui.hoc")
    h.load_file("loadbal.hoc")
    h('objref fi_status, fi_checksimtime, pc, nclist, nc, nil')
    h('strdef datasetPath')
    h('numCells = 0')
    h('totalNumCells = 0')
    h('max_walltime_hrs = 0')
    h('mkcellstime = 0')
    h('mkstimtime = 0')
    h('connectcellstime = 0')
    h('connectgjstime = 0')
    h('results_write_time = 0')
    h.nclist = h.List()
    datasetPath  = os.path.join(env.datasetPrefix, env.datasetName)
    h.datasetPath = datasetPath
    ##  new ParallelContext object
    h.pc   = h.ParallelContext()
    env.pc = h.pc
    ## polymorphic value template
    h.load_file("./templates/Value.hoc")
    ## randomstream template
    h.load_file("./templates/ranstream.hoc")
    ## stimulus cell template
    h.load_file("./templates/StimCell.hoc")
    h.xopen("./lib.hoc")
    h.dt = env.dt
    h.tstop = env.tstop
    if env.optldbal or env.optlptbal:
        lb = h.LoadBalance()
        if not os.path.isfile("mcomplex.dat"):
            lb.ExperimentalMechComplex()
    h.startsw()
    mkcells(env)
    env.mkcellstime = h.stopsw()
    env.pc.barrier()
    if (env.pc.id() == 0):
        print "*** Cells created in %g seconds" % env.mkcellstime
    print "*** Rank %i created %i cells" % (env.pc.id(), len(env.cells))
    h.startsw()
    mkstim(env)
    env.mkstimtime = h.stopsw()
    if (env.pc.id() == 0):
        print "*** Stimuli created in %g seconds" % env.mkstimtime
    env.pc.barrier()
    h.startsw()
    connectcells(env)
    env.connectcellstime = h.stopsw()
    env.pc.barrier()
    if (env.pc.id() == 0):
        print "*** Connections created in %g seconds" % env.connectcellstime
    print "*** Rank %i created %i connections" % (env.pc.id(), int(h.nclist.count()))
    h.startsw()
    #connectgjs(env)
    env.connectgjstime = h.stopsw()
    if (env.pc.id() == 0):
        print "*** Gap junctions created in %g seconds" % env.connectgjstime
    env.pc.setup_transfer()
    env.pc.set_maxstep(10.0)
    h.max_walltime_hrs   = env.max_walltime_hrs
    h.mkcellstime        = env.mkcellstime
    h.mkstimtime         = env.mkstimtime
    h.connectcellstime   = env.connectcellstime
    h.connectgjstime     = env.connectgjstime
    h.results_write_time = env.results_write_time
    h.fi_checksimtime    = h.FInitializeHandler("checksimtime(pc)")
    if (env.pc.id() == 0):
        print "dt = %g" % h.dt
        print "tstop = %g" % h.tstop
        h.fi_status = h.FInitializeHandler("simstatus()")
    h.v_init = env.v_init
    h.stdinit()
    h.finitialize(env.v_init)
    env.pc.barrier()
    if env.optldbal or env.optlptbal:
        cx(env)
        ld_bal(env)
        if env.optlptbal:
            lpt_bal(env)

# Run the simulation
def run (env):
    rank = int(env.pc.id())
    nhosts = int(env.pc.nhost())

    env.pc.barrier()
    env.pc.psolve(h.tstop)

    if (rank == 0):
        print "*** Simulation completed"

    spikeoutPath = "%s/%s_spikeout.h5" % (env.resultsPath, env.modelName)
    if (rank == 0):
        try:
            mkspikeout (env, spikeoutPath)
        except:
            pass
    spikeout(env, spikeoutPath, np.array(env.t_vec, dtype=np.float32), np.array(env.id_vec, dtype=np.uint32))

    # TODO:
    #if (env.vrecordFraction > 0):
    #    h.vrecordout("%s/%s_vrecord_%d.dat" % (env.resultsPath, env.modelName, env.pc.id(), env.indicesVrecord))

    comptime = env.pc.step_time()
    cwtime   = comptime + env.pc.step_wait()
    maxcw    = env.pc.allreduce(cwtime, 2)
    avgcomp  = env.pc.allreduce(comptime, 1)/nhosts
    maxcomp  = env.pc.allreduce(comptime, 2)

    if (env.pc.id() == 0):
        print "Execution time summary for host 0:"
        print "  created cells in %g seconds" % env.mkcellstime
        print "  connected cells in %g seconds" % env.connectcellstime
        print "  created gap junctions in %g seconds" % env.connectgjstime
        print "  ran simulation in %g seconds" % comptime
        if (maxcw > 0):
            print "  load balance = %g" % (avgcomp/maxcw)

    env.pc.runworker()
    env.pc.done()
    h.quit()

    
@click.command()
@click.option("--config-file", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--template-paths", type=str)
@click.option("--dataset-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option("--results-path", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option("--results-id", type=str, required=False, default='')
@click.option("--node-rank-file", required=False, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--io-size", type=int, default=1)
@click.option("--coredat", is_flag=True)
@click.option("--vrecord-fraction", type=float, default=0.0)
@click.option("--tstop", type=int, default=1)
@click.option("--v-init", type=float, default=-75.0)
@click.option("--max-walltime-hours", type=float, default=1.0)
@click.option("--results-write-time", type=float, default=30.0)
@click.option("--dt", type=float, default=0.025)
@click.option("--ldbal", is_flag=True)
@click.option("--lptbal", is_flag=True)
@click.option('--verbose', '-v', is_flag=True)
def main(config_file, template_paths, dataset_prefix, results_path, results_id, node_rank_file, io_size, coredat, vrecord_fraction, tstop, v_init, max_walltime_hours, results_write_time, dt, ldbal, lptbal, verbose):
    np.seterr(all='raise')
    env = Env(MPI.COMM_WORLD, config_file, 
              template_paths, dataset_prefix, results_path, results_id,
              node_rank_file, io_size,
              vrecord_fraction, coredat, tstop, v_init,
              max_walltime_hours, results_write_time,
              dt, ldbal, lptbal, verbose)
    init(env)
    run(env)

if __name__ == '__main__':
    main(args=sys.argv[(sys.argv.index("main.py")+1):])

