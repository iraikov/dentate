##
## This file computes an approximation to the LFP generated by the
## pyramidal cells in the network, based on the formula by Schomburg et
## al., J Neurosci 2012.
## 
## The approximate LFP is calculated as the sum of current contributions
## of all compartments, scaled by the distances to the recording
## electrode and extracellular medium resistivity.  The time resolution
## of the LFP calculation may be lower than that of the simulation by
## setting lfp_dt.

from neuron import h
import itertools
import math

## 
## Computes xyz coords of nodes in a model cell  whose topology & geometry are defined by pt3d data.
## Code by Ted Carnevale.
## 

def interpxyz(nn, nsegs, xx, yy, zz, ll, xint, yint, zint):
    
    ## To use Vector class's .interpolate() 
    ## must first scale the independent variable
    ## i.e. normalize length along centroid
    ll.div(ll.x[nn-1])
    
    ## initialize the destination "independent" vector
    range = h.Vector(nsegs+2)
    range.indgen(1/nsegs)
    range.sub(1/(2*nsegs))
    range.x[0]=0
    range.x[nsegs+1]=1
    
    ## length contains the normalized distances of the pt3d points 
    ## along the centroid of the section.  These are spaced at 
    ## irregular intervals.
    ## range contains the normalized distances of the nodes along the 
    ## centroid of the section.  These are spaced at regular intervals.
    ## Ready to interpolate.
    
    xint.interpolate(range, ll, xx)
    yint.interpolate(range, ll, yy)
    zint.interpolate(range, ll, zz)


class LFP:

    def __init__(self, pc, celltypes, pos, rho = 333.0, fdst = 0.1, MaxEDist=100., dt_lfp=0.5, seed=1):
        self.pc         = pc
        self.dt_lfp     = dt_lfp
        self.seed       = seed
        self.epoint     = pos
        self.MaxEDist   = MaxEDist
        self.rho        = rho ## extracellular resistivity, [ohm cm]
        self.fdst       = fdst ## percent of distant cells to include in the computation
        self.lfplist    = h.List()
        self.lfp_ids    = {}
        self.lfp_types  = {}
        self.lfpkmatrix = {}
        self.celltypes  = celltypes
        if (int(pc.id()) == 0):
            self.fih_lfp = h.FInitializeHandler(1, self.sample_lfp)
            
    def setup_lfpkmatrix():

        ex, ey, ez = self.epoint
        for pop_name in self.celltypes.keys():
            
            lfp_ids = self.lfp_ids[pop_name]
            lfpkmatrix = self.lfpkmatrix[pop_name]

            for i in xrange(0, lfp_ids.size()):
                ## Iterates over all cells chosen for the LFP computation
	    
                gid  = lfp_ids.x[i]
                cell = pc.gid2cell(gid)
	    
                ## Iterates over each compartment of the cell
                for sec in list(cell.all):
                    if h.ismembrane('extracellular',sec=sec):

                        nn = sec.n3d()
		    
                        xx = h.Vector(nn)
                        yy = h.Vector(nn)
                        zz = h.Vector(nn)
                        ll = h.Vector(nn)
		    
                        for ii in xrange(0,nn):
                            xx.x[ii] = sec.x3d(ii)
                            yy.x[ii] = sec.y3d(ii)
                            zz.x[ii] = sec.z3d(ii)
                            ll.x[ii] = sec.arc3d(ii)
		    
                        xint = h.Vector(sec.nseg+2)
                        yint = h.Vector(sec.nseg+2)
                        zint = h.Vector(sec.nseg+2)
                    
                        interpxyz(nn,nseg,xx,yy,zz,ll,xint,yint,zint)
		    
                        j = 0
                        for seg in sec:
                        
                            sx = xint.x[j]
                            sy = yint.x[j]
                            sz = zint.x[j]
			
                            ## l = L/nseg is compartment length 
                            ## r is the perpendicular distance from the electrode to a line through the compartment
                            ## h is longitudinal distance along this line from the electrode to one end of the compartment
                            ## s = l + h is longitudinal distance to the other end of the compartment
                            l = sec.L/sec.nseg
                            r = math.sqrt((ex-sx)*(ex-sx) + (ey-sy)*(ey-sy) + (ez-sz)*(ez-sz))
                            h = l/2
                            s = l + h
                            k = 0.0001 * sec.area(seg) * (self.rho / (4.0 * math.pi * l)) * math.abs(math.log((math.sqrt(h*h + r*r) - h) / (math.sqrt(s*s + r*r) - s)))
                            if math.isnan(k):
                                k = 0.
                            ## Distal cell
                            if (lfp_types.x[i] == 2):
                                k = (1.0/fdst)*k
                            ##printf ("host %d: npole_lfp: gid = %d i = %d j = %d r = %g h = %g k = %g\n", pc.id, gid, i, j, r, h, k)
                            lfpkmatrix.x[i][j] = k
                            j = j + 1

                            
    def setup_lfp(self):
    
        ex, ey, ez = self.epoint
    
        ##printf ("host %d: entering setup_npole_lfp" % int(self.pc.id()))
    
        ## Vector for storing longitudinal and perpendicular distances
        ## between recording electrode and compartments
        dists = h.Vector(2)
    
        ## Determine which cells will be used for the LFP computation and the sizes of their compartments
        for (ipop, pop_name) in enumerate(self.celltypes.keys()):

            ranlfp     = h.Random(self.seed + i)
            ranlfp.uniform(0, 1)

            lfp_ids    = h.Vector()
            lfp_types  = h.Vector()
            m = 0
            n = 0

            pop_start = self.celltypes[pop_name]['start']
            pop_num = self.celltypes[pop_name]['num']
            
            for gid in xrange(pop_start, pop_start+num):

                ransample = ranlfp.repick()

                if not env.pc.gid_exists(gid):
                    continue

                cell = pc.gid2cell(gid)
                ## Relative to the recording electrode position
                if (math.sqrt((cell.x-ex)*(cell.x-ex) + (cell.y-ey)*(cell.y-ey) + (cell.z-ez)*(cell.z-ez)) < self.MaxEDist):
                    lfptype = 1 ## proximal cell; compute extracellular potential
                else:
                    if (ransample < fdst):
                        lfptype = 2 ## distal cell -- compute extracellular potential only for fdst fraction of total
                    else:
                        lfptype = 0 ## do not compute extracellular potential

                if (lfptype > 0):
                    lfp_ids.append(gid)
                    lfp_types.append(lfptype)
                    m = m+1
                    if (n == 0):
                        for sec in list(cell.all):
                            sec.insert('extracellular')
                            n = n + sec.nseg

            self.lfp_ids[pop_name] = lfp_ids
            self.lfp_types[pop_name] = lfp_types
            if (m > 0):
                self.lfpkmatrix[pop_name] = h.Matrix(m,n)

            self.setup_lfpkmatrix()

            
    def pos_lfp(self):
        ## Calculate the average LFP of select pyramidal cells in the network,
        ##  only including pyramidal cells whose somata are within MaxEDist
        ##  microns of the (x,y,z) recording electrode location
	
        sumcell = 0
        vlfp = 0 

        for pop_name in self.celltypes.keys():
            lfp_ids = self.lfp_ids[pop_name]
            lfpkmatrix = self.lfpkmatrix[pop_name]
            ## Iterate over all cell types
            for i in xrange(0, lfp_ids.size()):
                ## Iterate over all cells chosen for the LFP computation
                gid  = lfp_ids.x[i]
                cell = pc.gid2cell(gid)

                for sec in list(cell.all):
                    if h.ismembrane('extracellular',sec=sec):
                        j = 0
                        for seg in sec:
                            vlfp = vlfp + (seg._ref_i_membrane * lfpkmatrix.x[i][j])
                            j = j + 1

        meanlfp = self.pc.allreduce(h.vlfp, 1)
        return meanlfp


    def sample_lfp(self):
    
        ## recording electrode position (um)
        ex, ey, ez = self.epoint
    
        ## At t=0, calculate distances from recording electrode to all
        ## compartments of all pyramidal cells, calculate scaling
        ## coefficients for the LFP calculation, and save them in
        ## lfpkmatrix.
        
        if (h.t == 0.):
            self.setup_lfp()
    
        ## Compute LFP across a subset of cells within a certain distance
        ## from the recording electrode:

        meanlfp = self.pos_lfp()
    
        if (int(self.pc.id()) == 0):
            vec = h.Vector()
            ## For this time step, create a vector with entries of time and average LFP
            vec.append(h.t, meanlfp)			
            ## Append the vector for this time step to the list
            self.lfplist.append(vec.c)			
    
        ## Add another event to the event queue, to 
        ## execute sample_lfp again, lfp_dt ms from now
        h.cvode.event(h.t + self.lfp_dt, self.sample_lfp)

