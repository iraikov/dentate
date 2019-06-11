"""An approximation to the LFP generated by the cells in the network, based on the method developed by Schomburg et al., J Neurosci 2012.
 
The approximate LFP is calculated as the sum of current contributions
of all compartments, scaled by the distances to the recording
electrode and extracellular medium resistivity.  The time resolution
of the LFP calculation may be lower than that of the simulation by
setting dt_lfp.
"""
from __future__ import division

from builtins import range
from builtins import object
from past.utils import old_div
from neuron import h
import itertools
import math

h('objref strfun')
h.strfun = h.StringFunctions()


def interpxyz(nn, nsegs, xx, yy, zz, ll, xint, yint, zint):
    """Computes xyz coords of nodes in a model cell  whose topology & geometry are defined by pt3d data.
    Code by Ted Carnevale.
    """

    ## To use Vector class's .interpolate() 
    ## must first scale the independent variable
    ## i.e. normalize length along centroid
    ll.div(ll.x[nn - 1])

    ## initialize the destination "independent" vector
    rangev = h.Vector(nsegs + 2)
    rangev.indgen(old_div(1, nsegs))
    rangev.sub(old_div(1, (2 * nsegs)))
    rangev.x[0] = 0
    rangev.x[nsegs + 1] = 1

    ## length contains the normalized distances of the pt3d points 
    ## along the centroid of the section.  These are spaced at 
    ## irregular intervals.
    ## range contains the normalized distances of the nodes along the 
    ## centroid of the section.  These are spaced at regular intervals.
    ## Ready to interpolate.
    xint.interpolate(rangev, ll, xx)
    yint.interpolate(rangev, ll, yy)
    zint.interpolate(rangev, ll, zz)


class LFP(object):

    def __init__(self, label, pc, celltypes, pos, rho=333.0, fdst=0.1, maxEDist=100., dt_lfp=0.5, seed=1):
        self.label = label
        self.pc = pc
        self.dt_lfp = dt_lfp
        self.seed = seed
        self.epoint = pos
        self.maxEDist = maxEDist
        self.rho = rho  ## extracellular resistivity, [ohm cm]
        self.fdst = fdst  ## percent of distant cells to include in the computation
        self.meanlfp = []
        self.t = []
        self.lfp_ids = {}
        self.lfp_types = {}
        self.lfp_coeffs = {}
        self.celltypes = celltypes
        self.fih_lfp = h.FInitializeHandler(1, self.sample_lfp)

    def setup_lfp_coeffs(self):

        ex, ey, ez = self.epoint
        for pop_name in self.celltypes:

            lfp_ids = self.lfp_ids[pop_name]
            lfp_types = self.lfp_types[pop_name]
            lfp_coeffs = self.lfp_coeffs[pop_name]

            for i in range(0, int(lfp_ids.size())):
                ## Iterates over all cells chosen for the LFP computation

                gid = lfp_ids.x[i]
                cell = self.pc.gid2cell(gid)

                ## Iterates over each compartment of the cell
                for sec in list(cell.all):
                    if h.ismembrane('extracellular', sec=sec):

                        nn = sec.n3d()

                        xx = h.Vector(nn)
                        yy = h.Vector(nn)
                        zz = h.Vector(nn)
                        ll = h.Vector(nn)

                        for ii in range(0, nn):
                            xx.x[ii] = sec.x3d(ii)
                            yy.x[ii] = sec.y3d(ii)
                            zz.x[ii] = sec.z3d(ii)
                            ll.x[ii] = sec.arc3d(ii)

                        xint = h.Vector(sec.nseg + 2)
                        yint = h.Vector(sec.nseg + 2)
                        zint = h.Vector(sec.nseg + 2)

                        interpxyz(nn, sec.nseg, xx, yy, zz, ll, xint, yint, zint)

                        j = 0
                        sx0 = xint.x[0]
                        sy0 = yint.x[0]
                        sz0 = zint.x[0]
                        for seg in sec:

                            sx = xint.x[j]
                            sy = yint.x[j]
                            sz = zint.x[j]

                            ## l = L/nseg is compartment length 
                            ## rd is the perpendicular distance from the electrode to a line through the compartment
                            ## ld is longitudinal distance along this line from the electrode to one end of the compartment
                            ## sd = l - ld is longitudinal distance to the other end of the compartment
                            l = old_div(sec.L, sec.nseg)
                            rd = math.sqrt((ex - sx) * (ex - sx) + (ey - sy) * (ey - sy) + (ez - sz) * (ez - sz))
                            ld = math.sqrt((sx - sx0) * (sx - sx0) + (sy - sy0) * (sy - sy0) + (sz - sz0) * (sz - sz0))
                            sd = l - ld
                            k = 0.0001 * h.area(seg.x) * (old_div(self.rho, (4.0 * math.pi * l))) * abs(math.log(
                                old_div((math.sqrt(ld * ld + rd * rd) - ld), (math.sqrt(sd * sd + rd * rd) - sd))))
                            if math.isnan(k):
                                k = 0.
                            ## Distal cell
                            if (lfp_types.x[i] == 2):
                                k = (1.0 / self.fdst) * k
                            ##printf ("host %d: npole_lfp: gid = %d i = %d j = %d r = %g h = %g k = %g\n", pc.id, gid, i, j, r, h, k)
                            lfp_coeffs.o(i).x[j] = k
                            j = j + 1

    def setup_lfp(self):

        ex, ey, ez = self.epoint

        ##printf ("host %d: entering setup_npole_lfp" % int(self.pc.id()))

        ## Determine which cells will be used for the LFP computation and the sizes of their compartments
        for (ipop, pop_name) in enumerate(self.celltypes):

            ranlfp = h.Random(self.seed + ipop)
            ranlfp.uniform(0, 1)

            lfp_ids = h.Vector()
            lfp_types = h.Vector()
            lfp_coeffs = h.List()

            pop_start = self.celltypes[pop_name]['start']
            pop_num = self.celltypes[pop_name]['num']

            for gid in range(pop_start, pop_start + pop_num):

                ransample = ranlfp.repick()

                if not self.pc.gid_exists(gid):
                    continue

                cell = self.pc.gid2cell(gid)

                if h.strfun.is_artificial(cell) > 0:
                    continue

                somasec = list(cell.soma)
                x = somasec[0].x3d(0)
                y = somasec[0].y3d(0)
                z = somasec[0].z3d(0)

                ## Relative to the recording electrode position
                if (math.sqrt((x - ex) ** 2 + (y - ey) ** 2 + (z - ez) ** 2) < self.maxEDist):
                    lfptype = 1  ## proximal cell; compute extracellular potential
                else:
                    if (ransample < self.fdst):
                        lfptype = 2  ## distal cell -- compute extracellular potential only for fdst fraction of total
                    else:
                        lfptype = 0  ## do not compute extracellular potential

                if (lfptype > 0):
                    lfp_ids.append(gid)
                    lfp_types.append(lfptype)
                    n = 0
                    for sec in list(cell.all):
                        sec.insert('extracellular')
                        n = n + sec.nseg
                    vec = h.Vector()
                    vec.resize(n)
                    lfp_coeffs.append(vec)

            self.lfp_ids[pop_name] = lfp_ids
            self.lfp_types[pop_name] = lfp_types
            self.lfp_coeffs[pop_name] = lfp_coeffs

        self.setup_lfp_coeffs()

    def pos_lfp(self):
        ## Calculate the average LFP of select cells in the network,
        ##  only including cells whose somata are within maxEDist
        ##  microns of the (x,y,z) recording electrode location

        vlfp = 0.

        for pop_name in self.celltypes:
            lfp_ids = self.lfp_ids[pop_name]
            lfp_coeffs = self.lfp_coeffs[pop_name]
            ## Iterate over all cell types
            for i in range(0, int(lfp_ids.size())):
                ## Iterate over the cells chosen for the LFP computation
                gid = lfp_ids.x[i]
                cell = self.pc.gid2cell(gid)

                for sec in list(cell.all):
                    if h.ismembrane('extracellular', sec=sec):
                        j = 0
                        for seg in sec:
                            vlfp = vlfp + (seg._ref_i_membrane[0] * lfp_coeffs.o(i).x[j])
                            j = j + 1

        meanlfp = self.pc.allreduce(vlfp, 1)
        return meanlfp

    def sample_lfp(self):

        ## recording electrode position (um)
        ex, ey, ez = self.epoint

        ## At t=0, calculate distances from recording electrode to all
        ## compartments of all cells, calculate scaling coefficients
        ## for the LFP calculation, and save them in lfp_coeffs.

        if (h.t == 0.):
            self.setup_lfp()

        ## Compute LFP across the subset of cells:
        meanlfp = self.pos_lfp()

        if (int(self.pc.id()) == 0):
            ## For this time step, append to lists with entries of time and average LFP
            self.meanlfp.append(meanlfp)
            self.t.append(h.t)

        ## Add another event to the event queue, to 
        ## execute sample_lfp again, dt_lfp ms from now
        h.cvode.event(h.t + self.dt_lfp, self.sample_lfp)
