import sys, os, time, copy 
from mpi4py import MPI
import numpy as np
import itertools as it
import scipy as sp
import scipy.stats as spst
from sklearn.gaussian_process import GaussianProcessRegressor as GPR
import h5py, yaml
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.backends.backend_pdf import PdfPages
import pprint

class NetClampParam:
    def __init__(self, fils, fil_dir, prefix=None):
        self.N_fils = len(fils)
        self.fil_arr = np.empty(shape=(self.N_fils,2), dtype='O')
        for idx, fil in enumerate(fils):
            self.fil_arr[idx,0] = h5py.File(os.path.join(fil_dir, fil), 'r')
            head_group = [i for i in self.fil_arr[0,0]][0]
            self.fil_arr[idx,1] = np.sort(self.fil_arr[idx,0][head_group]['problem_ids']) 
        
        self.head_group = head_group 
        self.ref_point = self.fil_arr[0,0][self.head_group]['{:d}'.format(self.fil_arr[0,1][0])]
        ref_features = np.array(self.ref_point['features'])
        self.N_cells = self.fil_arr[idx,1].shape[0]
        self.N_params, self.N_objectives, self.N_trials = ref_features['trial_objs'].shape

        self.get_params()

        self.target = np.array(self.fil_arr[0,0][self.head_group]['metadata'])[0][0]

        ts = time.strftime("%Y%m%d_%H%M%S")
        prefix = '' if prefix is None else '{!s}_'.format(prefix)
        self.suffix = 'Combined' if self.N_fils > 1 else '{:08d}'.format(int(self.fil_arr[idx,0].filename.split('_')[-1][:-3]))
        plot_save_ext = 'pdf'
        
        self.filnam = '{!s}{!s}_{!s}_{!s}'.format(prefix, self.population, ts, self.suffix)
        self.plot_filnam = '{!s}.{!s}'.format(self.filnam, plot_save_ext) 

        self.get_best_arrs([self.target])

        self.pop_props = {
                        'AAC':  ('AxoAxonic', -65, -42),
                        'BC':   ('PV+ Basket', -62, -43.1),
                        'HC':   ('HIPP', -67, -28.1),
                        'HCC':  ('HICAP', -65, -36.6),
                        'IS':   ('IS', -69.7, -44.3),
                        'MOPP': ('MOPP', -75.6, -30),
                        'NGFC': ('NGFC', -75.6, -30),
                        'MC':   ('Mossy'),
                        }

        self.plot_best()

#        self.get_param_obj_tol()
    #    self.get_param_criteria_arr()


       # self.generate_yaml()

    def get_params(self):
        raw_param_dtype = np.array(self.ref_point['parameters']).dtype
        self.param_names = np.array(raw_param_dtype.names)

        self.param_dtype = np.dtype([(nam, raw_param_dtype.fields[nam][0]) for nam in self.param_names])
        self.param_dim = self.param_names.shape[0] 

        param_dtype = np.dtype([('population', 'U8'), ('presyn', list), ('loc', list), ('syn', 'U8'), ('prop', 'U12'), ('lo_bound', 'f8'), ('up_bound', 'f8'), ('presyn_lab', 'U64')])
        self.param_props = np.empty(shape=self.param_dim, dtype=param_dtype)

        self.param_tup = []
        for pidx, parm in enumerate(self.param_names):
            tmp = parm.split('.')
            popu = tmp[0]
            presyn = tmp[1][1:-1].replace("'", "").split(', ') if tmp[1].startswith('(') else [tmp[1]] 
            loc = tmp[2][1:-1].replace("'", "").split(', ') if tmp[2].startswith('(') else [tmp[2]] 
            syn = tmp[3]
            prop = tmp[4]
            lo_bnd = self.fil_arr[0,0][self.head_group]['parameter_spec']['lower'][pidx]
            up_bnd = self.fil_arr[0,0][self.head_group]['parameter_spec']['upper'][pidx]
            presyn_lab = ', '.join(psyn for psyn in presyn) 

            self.param_props[pidx] = popu, presyn, loc, syn, prop, lo_bnd, up_bnd, presyn_lab

        self.param_arr = np.empty(shape=(self.N_fils, self.N_cells, self.N_params), dtype=self.param_dtype) 

        self.population = self.param_props[0]['population'] 

    def get_best_arrs(self, target):
        self.obj_val_mean_arr = np.empty(shape=(self.N_fils, self.N_cells, self.N_params, self.N_objectives))

        self.best_idx = np.empty(shape=(self.N_fils, self.N_cells, self.N_objectives, self.N_trials), dtype=np.uint32)
        self.best_arr = np.empty(shape=(self.N_fils, self.N_cells, self.N_objectives, self.N_trials))
        self.best_Vmean = np.empty(shape=(self.N_fils, self.N_cells, self.N_objectives, self.N_trials))
        self.best_prm = np.empty(shape=(self.N_fils, self.N_cells, self.N_objectives, self.N_trials), dtype=self.param_dtype)

        self.bestmean_idx = np.empty(shape=(self.N_fils, self.N_cells, self.N_objectives), dtype=np.uint32)
        self.bestmean_arr = np.empty(shape=(self.N_fils, self.N_cells, self.N_objectives))
        self.bestmean_prm = np.empty(shape=(self.N_fils, self.N_cells, self.N_objectives), dtype=self.param_dtype)
        self.bestmean_V_mean = np.empty(shape=(self.N_fils, self.N_cells))

        for fidx, filobj in enumerate(self.fil_arr[:,0]):
            for cidx, cell in enumerate(self.fil_arr[fidx, 1]):
                ref_cell = filobj[self.head_group]['{:d}'.format(cell)] 
                for prm in self.param_names:
                    self.param_arr[fidx, cidx, :][prm] = ref_cell['parameters'][prm] 
                ref_features = ref_cell['features'] 
                self.obj_val_mean_arr[fidx, cidx] = np.array(ref_features['mean_rate']) 
                for i in range(self.N_objectives):
                    np.argmin(np.abs(ref_features['trial_objs'][:, i ,:]-target[i]), axis=0, out=self.best_idx[fidx,cidx, i])
                    self.bestmean_idx[fidx, cidx, i] = np.argmin(np.abs(self.obj_val_mean_arr[fidx, cidx, :, i]-target[i]))

                for pidx, val in np.ndenumerate(self.best_idx[fidx,cidx]):
                    self.best_arr[fidx, cidx, pidx[0], pidx[1]] = ref_features['trial_objs'][val, pidx[0], pidx[1]]    
                    self.best_prm[fidx, cidx, pidx[0], pidx[1]] = self.param_arr[fidx, cidx, val]
                    self.best_Vmean[fidx, cidx, pidx[0], pidx[1]] = ref_features['mean_v'][val, pidx[1]] 

                for pidx, val in np.ndenumerate(self.bestmean_idx[fidx,cidx]):
                    self.bestmean_arr[fidx, cidx, pidx] = self.obj_val_mean_arr[fidx, cidx, val, pidx]
                    self.bestmean_prm[fidx, cidx, pidx] = self.param_arr[fidx, cidx, val] 



    def generate_yaml(self, CriteriaList=None, gid_val_arr=None):
        main_dict = {self.population.item(): {gid.item(): [[str(parm['population']), 
                           [str(par).strip(',') for par in  parm['presyn']] if len(parm['presyn'])>1 else str(parm['presyn'][0]).strip(','), 
                           [str(loc).strip(',') for loc in parm['loc']] if len(parm['loc'])>1 else str(parm['loc'][0]).strip(','), 
                           str(parm['syn']).strip(','), str(parm['prop']).strip(','), None] 
                     for parm in self.param_props]
                for gid in self.fil_arr[0,1]}}

        if CriteriaList is None:
            CriteriaList, gid_val_arr = self.get_param_criteria_arr(CriteriaList=CriteriaList)

        N_crit = len(CriteriaList)

        crit_dict_arr = np.empty(shape=N_crit, dtype='O')
        crit_yaml_arr = np.empty(shape=N_crit, dtype='O')

        for critidx, crit in enumerate(CriteriaList):
            crit_dict_arr[critidx] = copy.deepcopy(main_dict)
            for gididx, gid in enumerate(self.fil_arr[0,1]):
                for val, prmidx  in zip(gid_val_arr[critidx, 0, :, gididx], range(self.param_dim)): 
                    crit_dict_arr[critidx][self.population][gid][prmidx][-1] = val.item()

            crit_yaml_arr[critidx] = '{!s}_{!s}.yaml'.format(self.filnam, crit)

        for datadict, yamlfil in zip(crit_dict_arr, crit_yaml_arr):
            yaml.dump(datadict, open(yamlfil, 'w'), default_flow_style=True, Dumper=yaml.CDumper)


    def get_param_criteria_arr(self, CriteriaList=None):
        if not hasattr(self, 'best_prm'):
            self.get_best_arrs([self.target]) 

    
        criteria_dict = {
                            'UniformBestMean': "get_best()",
                            'SpecificBestMean': "get_best(uniform=False)",
                            'UniformTrialMean': "get_trial()",
                            'SpecificTrialMean': "get_trial(uniform=False)",
                            'UniformBestMedian': "get_best(fn='median')",
                            'SpecificBestMedian': "get_best(fn='median', uniform=False)",
                            'UniformTrialMedian': "get_trial(fn='median')",
                            'SpecificTrialMedian': "get_trial(fn='median', uniform=False)",
                            'UniformBestMode': "get_best_mode()",
                            'SpecificBestMode': "get_best_mode(uniform=False)",
                            'UniformTrialMode': "get_trial_mode()",
                            'SpecificTrialMode': "get_trial_mode(uniform=False)",
            #                'UniformBestSpeciifc': "get_specific()",
                    }

        def get_best(fn='mean', uniform=True):
            func = getattr(np, fn)
            axis = (0,) if not uniform else None           
            tmp_arr = np.empty(shape=(self.N_objectives, self.param_dim, self.N_cells))
            for prmidx, parm in enumerate(self.raw_param_dtype.names):
               # for objidx in range(self.N_objectives):
               #     vparm = func(self.bestmean_prm[parm][:,:,objidx,:], axis=axis)
               #     tmp_arr[objidx, prmidx,:] = vparm
                vparm = func(self.bestmean_prm[parm], axis=axis)
                tmp_arr[0, prmidx,:] = vparm
            return tmp_arr 

        def get_trial(fn='mean', uniform=True):
            func = getattr(np, fn)
            axis = (0,-1) if not uniform else None           
            tmp_arr = np.empty(shape=(self.N_objectives, self.param_dim, self.N_cells))
            for prmidx, parm in enumerate(self.raw_param_dtype.names):
                for objidx in range(self.N_objectives):
                    vparm = func(self.best_prm[parm][:,:,objidx,:], axis=axis)
                    tmp_arr[objidx, prmidx,:] = vparm
            return tmp_arr 

        def get_best_mode(uniform=True):
            axis = 0 if not uniform else None           
            tmp_arr = np.empty(shape=(self.N_objectives, self.param_dim, self.N_cells))
            for prmidx, parm in enumerate(self.raw_param_dtype.names):
               # for objidx in range(self.N_objectives):
               #     vparm = func(self.bestmean_prm[parm][:,:,objidx,:], axis=axis)
               #     tmp_arr[objidx, prmidx,:] = vparm
                vparm, _ = spst.mode(self.bestmean_prm[parm], axis=axis)
                tmp_arr[0, prmidx,:] = vparm
            return tmp_arr 

        def get_trial_mode(uniform=True):
            tmp_arr = np.empty(shape=(self.N_objectives, self.param_dim, self.N_cells))
            if uniform:
                for prmidx, parm in enumerate(self.raw_param_dtype.names):
                    for objidx in range(self.N_objectives):
                        vparm, _ = spst.mode(self.best_prm[parm][:,:,objidx,:], axis=None)
                        tmp_arr[objidx, prmidx,:] = vparm
            else:
                for prmidx, parm in enumerate(self.raw_param_dtype.names):
                    for objidx in range(self.N_objectives):
                        for cidx in range(self.N_cells):
                            vparm, _ = spst.mode(np.ravel(self.best_prm[parm][:,cidx,objidx,:]))
                            tmp_arr[objidx, prmidx, cidx] = vparm
            return tmp_arr 

        def generate_yaml(Uniform=True):
            Critlist = []
            gid_val_arr = np.empty(shape=(self.N_cells, self.N_objectives, self.param_dim, self.N_cells))
            for cidx, cell in enumerate(self.fil_arr[fidx, 1]):
                ref_cell = filobj[self.head_group]['{:d}'.format(cell)] 
                Critlist.append('Uniform_BestGid_{!s}'.format(cell))
                bestoversam_idx = np.argmin(np.abs(self.bestmean_arr[:,cidx]-target[0]))  
                bestoversam_prm = self.bestmean_prm[bestoversam_idx, cidx] 

                parms = np.array([bestoversam_prm[prm] for prm in self.raw_param_dtype.names])
                for i in range(self.N_cells):
                    gid_val_arr[cidx,0,:,i] = parms 
                print(self.population, cell, bestoversam_idx, parms)
            self.generate_yaml(Critlist, gid_val_arr)

    #    generate_yaml()

        CriteriaList = list(criteria_dict.keys()) if CriteriaList is None else CriteriaList
        N_crit = len(CriteriaList)

        gid_val_arr = np.empty(shape=(N_crit, self.N_objectives, self.param_dim, self.N_cells))

        for critidx, crit in enumerate(CriteriaList):
            gid_val_arr[critidx] = eval(criteria_dict[crit]) 

        return (CriteriaList, gid_val_arr) 

    def get_obj_axes(self):
    
        fig_h = self.N_objectives * 4 + 5
        hratio = [self.N_objectives, 1]
        fig = plt.figure(figsize=(20,fig_h), constrained_layout=True)

        spec = gs.GridSpec(nrows=2, ncols=1, figure=fig, height_ratios=hratio)

        optobjspec = spec[0,0].subgridspec(self.N_objectives, self.N_cells)
        objspec = spec[1,0].subgridspec(1, 2)

        opt_axes = np.empty(shape=optobjspec.get_geometry(), dtype='O')
        obj_axes = np.empty(shape=self.N_objectives, dtype='O')
        fea_axes = np.empty(shape=1, dtype='O')

        for idx, ax in np.ndenumerate(opt_axes): 
            opt_axes[idx] = fig.add_subplot(optobjspec[idx])

        for idx, ax in np.ndenumerate(obj_axes): 
            obj_axes[idx] = fig.add_subplot(objspec[0,idx[0]])

        fea_axes[0] = fig.add_subplot(objspec[0, 1])

        return fig, opt_axes, obj_axes, fea_axes

    def get_param_axes(self):
        self.N_comb = (self.param_dim**2-self.param_dim)//2
        corprmgrd, corprmcrds = self.get_grid_shape(self.N_comb)
        hratio = [1,1, corprmgrd[0]]

        fig = plt.figure(figsize=(20, 6+corprmgrd[0]*3), constrained_layout=True)
        spec = gs.GridSpec(nrows=3, ncols=1, figure=fig, height_ratios=hratio)

        prmspec = spec[0,0].subgridspec(1, self.param_dim)
        optprmspec = spec[1,0].subgridspec(1, self.param_dim)
      #  corprmspec = spec[2,0].subgridspec(1, (self.param_dim**2-self.param_dim)//2)

        corprmspec = spec[2,0].subgridspec(*corprmgrd)

        prm_axes = np.empty(prmspec.get_geometry(), dtype='O')
        opr_axes = np.empty(optprmspec.get_geometry(), dtype='O')
        cpr_axes = np.empty(corprmspec.get_geometry(), dtype='O')

        for idx, ax in np.ndenumerate(prm_axes): 
            prm_axes[idx] = fig.add_subplot(prmspec[idx])

        for idx, ax in np.ndenumerate(opr_axes): 
            opr_axes[idx] = fig.add_subplot(optprmspec[idx])

     #   for idx, ax in np.ndenumerate(cpr_axes): 
        for idx in corprmcrds:
            cpr_axes[idx] = fig.add_subplot(corprmspec[idx])


        return fig, prm_axes, opr_axes, cpr_axes, corprmcrds 


    def plot_best(self):
    #    param_names = np.array(self.param_dtype.names)
    #    param_dim = self.param_dim 
        N_iter_arr = np.arange(self.N_params)

# TODO send to summary page
 #       prm_str = ''
 #       for idx, prm in enumerate(param_names):
 #           prm_str += '{!s}: {:f}, {:f} '.format(self.param_props['syn'][idx], np.mean(self.best_prm[prm]), np.mean(self.bestmean_prm[prm]))

        self.prmcomb = np.array([i for i in it.combinations(range(self.param_dim), r=2)], dtype=np.uint32)
        N_prmcomb = self.prmcomb.shape[0]

        xcat = self.fil_arr[0,1]
        boxwidth = (xcat[-1]-xcat[0])/10
        plotshap = ['o', 's', 'D', 'X', '^']
        N_fil_cell_max = max(self.N_fils, self.N_cells)
        colors = plt.get_cmap('tab10').colors if N_fil_cell_max <=10 else mpl.cm.get_cmap('tab20').colors
        self.colors = colors

        obj_fig, opt_axes, obj_axes, fea_axes = self.get_obj_axes() 
        prm_fig, prm_axes , opr_axes, cpr_axes, prmcoords  = self.get_param_axes() 

        titnote='Marker shape: sampling set; when axis is specific for gid, colors represent trials, else gids; black represents mean_rate (objective)'
        obj_fig.suptitle('{!s}\n{!s}\n{!s}'.format(self.pop_props[self.population][0], self.suffix, titnote))
#        prm_fig.suptitle('{!s}\n{!s}\n{!s}\n{!s}'.format(self.pop_props[self.population][0], self.suffix, titnote, prm_str))

        for i in range(self.N_objectives): 
            ytop = self.target[i]*5 if self.target[i]<10 else self.target[i]*3
            for j in range(self.N_fils):
                obj_axes[i].scatter(xcat, self.bestmean_arr[j, :], marker=plotshap[j], color='k')
                obj_axes[i].axhline(self.target[0], color='k', ls='--')
                for tidx in range(self.N_trials):
                    obj_axes[i].scatter(xcat, self.best_arr[j,:, i,tidx], marker=plotshap[j], color=colors[tidx], facecolors='none', lw=0.5)

                    fea_axes[0].scatter(xcat, self.best_Vmean[j,:, i,tidx], marker=plotshap[j], color=colors[tidx], facecolors='none', lw=0.5)

            obj_axes[i].set_ylabel(r'Firing Rate [spikes/s]')

        fea_axes[0].set_ylabel(r'Mean voltage [mV]')
        fea_axes[0].axhline(self.pop_props[self.population][1], color='k', lw=0.5)
        fea_axes[0].axhline(self.pop_props[self.population][2], color='k', lw=0.5)

        for i, prm in enumerate(self.param_names): 
            prm_axes[0,i].axhline(self.param_props['lo_bound'][i], color='k', ls='--')
            prm_axes[0,i].axhline(self.param_props['up_bound'][i], color='k', ls='--')
            prm_axes[0,i].set_ylabel(r'${!s}$ weight'.format(self.param_props['syn'][i]))
            tit = self.param_props['presyn_lab'][i] 
            prm_axes[0,i].set_title(tit)

            opr_axes[0,i].axhline(self.param_props['lo_bound'][i], color='k', ls='--')
            opr_axes[0,i].axhline(self.param_props['up_bound'][i], color='k', ls='--')
            opr_axes[0,i].set_ylabel(r'${!s}$ weight'.format(self.param_props['syn'][i]))
            opr_axes[0,i].set_xlabel(r'Iteration Index')

            for cidx, cell in enumerate(xcat): 
                prm_axes[0,i].boxplot(np.ravel(self.best_prm[prm][:,cidx, 0,:]), vert=True, positions=[cell], showmeans=True, widths=boxwidth)

            for j in range(self.N_fils):
                prm_axes[0,i].scatter(xcat, self.bestmean_prm[prm][j, :], marker=plotshap[j], color='k')
                for tidx in range(self.N_trials):
                    prm_axes[0,i].scatter(xcat, self.best_prm[prm][j,:, 0,tidx], marker=plotshap[j], color=colors[tidx], facecolors='none', lw=0.5)
                
                for cidx, cell in enumerate(xcat): 
                    opr_axes[0,i].scatter(N_iter_arr, self.param_arr[j, cidx][prm],  marker=plotshap[j], color=colors[cidx], facecolors='none', lw=0.01, s=1)
                    opr_axes[0,i].scatter(self.bestmean_idx[j,cidx], self.bestmean_prm[prm][j,cidx],  marker=plotshap[j], color=colors[cidx], edgecolor='k', lw=0.5, s=15, zorder=100)

                    for tidx in range(self.N_trials):
                        opr_axes[0,i].scatter(self.best_idx[j,cidx,0,tidx], self.best_prm[prm][j,cidx, 0, tidx],  marker=plotshap[j], color=colors[cidx], s=4)


        opt_val_arr = np.empty(shape=(self.N_fils, self.N_cells, self.N_params, self.N_objectives, self.N_trials))
        for fidx, filobj in enumerate(self.fil_arr[:,0]):
            for cidx, cell in enumerate(self.fil_arr[fidx, 1]):
                ref_cell = filobj[self.head_group]['{:d}'.format(cell)] 
                opt_val_arr[fidx,cidx] = ref_cell['features']['trial_objs']
                ref_ax = opt_axes[0,cidx]
                ref_ax.scatter(N_iter_arr, self.obj_val_mean_arr[fidx, cidx,:,0], marker=plotshap[fidx], color='k', facecolors='none', lw=0.01, s=1, zorder=200)
                ref_ax.axhline(self.target[0], color='k', ls='--', zorder=0.25)
                for tidx in range(self.N_trials):
                    ref_ax.scatter(N_iter_arr, opt_val_arr[fidx,cidx,:,0,tidx], marker=plotshap[fidx], color=colors[tidx], facecolors='none', lw=0.01, s=1, zorder=1)

                for combidx, coord in zip(range(N_prmcomb), prmcoords):
                    comb  = self.prmcomb[combidx, : ]
                    i, j = comb

                    cpr_axes[coord].axvline(self.param_props['lo_bound'][i], color='k', ls='--')
                    cpr_axes[coord].axvline(self.param_props['up_bound'][i], color='k', ls='--')
                    cpr_axes[coord].axhline(self.param_props['lo_bound'][j], color='k', ls='--')
                    cpr_axes[coord].axhline(self.param_props['up_bound'][j], color='k', ls='--')

                    cpr_axes[coord].set_xlabel(r'${!s}$ {!s}'.format(self.param_props['syn'][i], self.param_props['presyn_lab'][i]))
                    cpr_axes[coord].set_ylabel(r'${!s}$ {!s}'.format(self.param_props['syn'][j], self.param_props['presyn_lab'][j]))
                    cpr_axes[coord].scatter(self.param_arr[fidx, cidx][self.param_names[i]],  self.param_arr[fidx, cidx][self.param_names[j]], marker=plotshap[fidx], color=colors[cidx], facecolors='none', lw=0.01, s=1)

                    cpr_axes[coord].scatter(np.ravel(self.best_prm[self.param_names[i]][fidx,cidx, 0,:]), np.ravel(self.best_prm[self.param_names[j]][fidx,cidx, 0,:]), marker=plotshap[fidx], color=colors[cidx], s=4)
                    cpr_axes[coord].scatter(np.ravel(self.bestmean_prm[self.param_names[i]][fidx,cidx]), np.ravel(self.bestmean_prm[self.param_names[j]][fidx,cidx]), marker=plotshap[fidx], color=colors[cidx], edgecolor='k', lw=0.5, s=15, zorder=100)



        for cidx, cell in enumerate(self.fil_arr[0, 1]):
            ref_ax = opt_axes[0,cidx]
            perc = np.quantile(self.obj_val_mean_arr[:,cidx,:,0], [0.25,0.5, 0.75], axis=(0,), keepdims=False)
            ref_ax.plot(N_iter_arr, perc[1,:], color='saddlebrown', lw=0.4, zorder=0.75)
            ref_ax.fill_between(N_iter_arr, y1=perc[2,:], y2=perc[0,:], color='peachpuff', alpha=0.5, zorder=0.5)
            ref_ax.set_title('{:d}'.format(int(cell)), color=colors[cidx])
            ref_ax.set_ylim(top=ytop, bottom=0)

        opt_axes[0,0].set_ylabel(r'Firing Rate [spikes/s]')
        opt_axes[0,0].set_xlabel(r'Iteration Index')

        summ_fig = self.get_summ_fig()
    #    tol_fig = self.get_param_obj_tol()
        #fig.savefig(self.plot_filnam, transparent=True)
     #   self.plot_filnam = 'Debug.pdf'
        with PdfPages(self.plot_filnam) as pdf:
            pdf.savefig(obj_fig, transparent=True)
            pdf.savefig(prm_fig, transparent=True)
            pdf.savefig(summ_fig, transparent=True)
        #    pdf.savefig(tol_fig, transparent=True)

    def get_param_obj_tol(self, tol=0.01):
        prm_idx = np.empty(shape=(self.N_cells, self.N_objectives), dtype=np.ndarray) 
        prm_arr = np.empty(shape=(self.N_cells, self.N_objectives), dtype=np.ndarray) 
        fea_arr = np.empty(shape=(self.N_cells, self.N_objectives), dtype=np.ndarray) 

        colors = plt.get_cmap('tab10').colors if self.N_cells <=10 else mpl.cm.get_cmap('tab20').colors

        for i in range(self.N_objectives):
            for cidx in range(self.N_cells):
                tgt = np.array([1.-tol, 1.+tol]) * self.target[i] 
                objs = self.obj_val_mean_arr[:,cidx,:,i]
                prm_idx[cidx, i] = np.where((objs>= tgt[0]) & (objs<=tgt[1]))
                N_prm = prm_idx[cidx, i][0].shape[0]
                prm_arr[cidx, i] = np.empty(shape=N_prm, dtype=self.param_dtype) 
                fea_arr[cidx, i] = np.empty(shape=N_prm)

                for idx in range(N_prm):
                    filidx = prm_idx[cidx, i][0][idx]
                    prmidx = prm_idx[cidx, i][1][idx]

                    prm_arr[cidx,i][idx] = self.param_arr[filidx, cidx, prmidx]
                    fea_arr[cidx,i][idx] = self.obj_val_mean_arr[filidx,cidx,idx,i] 

        results = self.get_param_reg(prm_arr, fea_arr)

        corprmgrd, prmcoords = self.get_grid_shape(self.N_comb)
        N_prmcomb = len(prmcoords)

        self.N_obj_groups =  self.N_objectives

        fig = plt.figure(figsize=(20, self.N_obj_groups * corprmgrd[0]*3), constrained_layout=True)
        spec = gs.GridSpec(nrows=self.N_obj_groups, ncols=1, figure=fig)

        corprmspec = np.empty(shape=(self.N_obj_groups,2), dtype='O') 
        for gidx in range(self.N_obj_groups):
            corprmspec[gidx,0] = spec[gidx,0].subgridspec(*corprmgrd)
            corprmspec[gidx,1] = np.empty(corprmspec[gidx,0].get_geometry(), dtype='O') 
            cpr_axes = corprmspec[gidx,1]
            for idx in prmcoords:
                cpr_axes[idx] = fig.add_subplot(corprmspec[gidx,0][idx])

            for cidx in range(self.N_cells):
                for combidx, coord in zip(range(N_prmcomb), prmcoords):
                    comb  = self.prmcomb[combidx, : ]
                    i, j = comb

                    cpr_axes[coord].axvline(self.param_props['lo_bound'][i], color='k', ls='--')
                    cpr_axes[coord].axvline(self.param_props['up_bound'][i], color='k', ls='--')
                    cpr_axes[coord].axhline(self.param_props['lo_bound'][j], color='k', ls='--')
                    cpr_axes[coord].axhline(self.param_props['up_bound'][j], color='k', ls='--')

                    cpr_axes[coord].set_xlabel(r'${!s}$ {!s}'.format(self.param_props['syn'][i], self.param_props['presyn_lab'][i]))
                    cpr_axes[coord].set_ylabel(r'${!s}$ {!s}'.format(self.param_props['syn'][j], self.param_props['presyn_lab'][j]))
                    cpr_axes[coord].scatter(prm_arr[cidx, gidx][self.param_names[i]], prm_arr[cidx, gidx][self.param_names[j]], color=colors[cidx], facecolors=colors[cidx], lw=0.1, s=2)

        fig.suptitle('Tolerance level: {:.3f}'.format(tol))

        return fig

    def get_param_reg(self, prm_arr, fea_arr, multi_obj=False):

        if not multi_obj:
            GPR_objs = np.empty(shape=(self.N_cells, self.N_objectives), dtype='O') 
            for cidx in range(self.N_cells):
                N_X = prm_arr[cidx, 0].shape[0]
                X = np.empty(shape=(N_X, self.param_dim))
                for i in range(self.N_objectives):
                    Y = fea_arr[cidx, i]
                    for prmidx, prm in enumerate(self.param_names):
                        X[:,prmidx] = prm_arr[cidx,i][prm] 
                    GPR_objs[cidx, i] = GPR().fit(X, Y)
                  #  GPR_params = GPR_objs[cidx, i].get_params()
                    GPR_params = GPR_objs[cidx, i].log_marginal_likelihood_value_
                    y_predicted = GPR_objs[cidx, i].predict(X, return_cov=True)
                   # scor = GPR_objs[cidx, i].score(X, Y)
#                    print(N_X, GPR_params.shape)
                    print(y_predicted[1].shape)
                        

    def get_summ_fig(self):
        fig = plt.figure(figsize=(20, 10), constrained_layout=True)
        spec = gs.GridSpec(nrows=1, ncols=2, figure=fig)

        ax_summ = fig.add_subplot(spec[0,0])
        ax_summ.axis('off')

        nz_masked = np.ma.masked_equal(self.obj_val_mean_arr[:,:,:,0], 0.0, copy=False)
        minfir_idx = np.unravel_index(np.argmin(np.abs(nz_masked)), nz_masked.shape)

        param = self.param_arr[minfir_idx] 

        table = [[i['presyn'], i['syn'], param[idx]] for idx, i in enumerate(self.param_props)]
        collab = ['PreSynaptic', 'Synapse', 'Min Fir.']

        best_table = mpl.table.table(ax=ax_summ, cellText=table, colLabels=collab, loc=0) 
    
        return fig 

    def get_grid_shape(self, num, base_ncols=6):
        n_rows = int(np.ceil(num/base_ncols))
        n_cols = int(np.ceil(num/n_rows))
        if num == n_rows*n_cols: 
            if n_rows == 1:
                coords = [(0,j) for j in range(n_cols)]
            else:
                coords = [(i,j) for i in range(n_rows) for j in range(n_cols)]
        else:
            coords = [(i,j) for i in range(n_rows-1) for j in range(n_cols)]
            rem = num - (n_rows-1) * n_cols
            for j in range(rem):
                coords.append((n_rows-1, j))

        return (n_rows, n_cols), coords

def distribute_chores(fil_list, fil_dir, Combined=True, prefix=None):
    N_fils = len(fil_list)

    if Combined:
        chores = fil_list
        N_chores = N_fils
    else:
        chores = []
        for fil_set in fil_list:
            for fil in fil_set:
                chores.append([fil])
        N_chores = len(chores)
#
    mpicom = MPI.COMM_WORLD

    N_hosts = mpicom.Get_size()
    Rank = mpicom.Get_rank()
    
    my_chores = np.arange(Rank, N_chores, N_hosts)

    for idx in my_chores:
        NetClampParam(chores[idx], fil_dir=fil_dir, prefix=prefix)

def get_list(fil_dir, tag, optimizer='distgfs', ext='h5', groupSubpop=True):
    fil_lst_raw = []
    fil_lst = []
    subpop_lst = ['GC_', 'MC_', 'AAC_', 'BC_', 'HC_', 'HCC_', 'IS_', 'MOPP_', 'NGFC_']
    subpop_dict = {k:[] for k in subpop_lst}

    if type(tag) is not list:
        print('Invalid tag type, must be a list of strings.')
    else:
        fil_lst_raw = [f for f in os.listdir(fil_dir) if f.startswith(optimizer) and f.endswith(ext) and any(tg in f for tg in tag)]

    fil_lst_raw = list(set(fil_lst_raw))

    if groupSubpop:
        for f in fil_lst_raw:
            subpop = f.split('.')[2].split('_')[0]
            subpop_dict[subpop+'_'].append(f)
        for k, v in subpop_dict.items():
            if len(v):
                fil_lst.append(v)
    else:
        fil_lst = fil_lst_raw

    return fil_lst

    

    
if __name__ == '__main__':

    remote = False if sys.platform == 'darwin' else True

    fil_dir='/scratch1/04119/pmoolcha/HDM/dentate/results/netclamp' if remote else '/Volumes/Work/SolteszLab/HDM/dentate/results/netclamp'
        

    interneuron_opt  = [ 
[
            'distgfs.network_clamp.AAC_20210422_193222_04893658.h5',
            'distgfs.network_clamp.AAC_20210422_193222_27137089.h5',
            'distgfs.network_clamp.AAC_20210422_193223_36010476.h5',
            'distgfs.network_clamp.AAC_20210422_193223_49937004.h5',
            'distgfs.network_clamp.AAC_20210422_193223_53499406.h5',
],
[
            'distgfs.network_clamp.BC_20210422_193223_93454042.h5',
            'distgfs.network_clamp.BC_20210422_193227_52357252.h5',
            'distgfs.network_clamp.BC_20210422_193228_01503844.h5',
            'distgfs.network_clamp.BC_20210422_193228_28135771.h5',
],
[
            'distgfs.network_clamp.HC_20210422_193233_15879716.h5',
            'distgfs.network_clamp.HC_20210422_193236_28682721.h5',
            'distgfs.network_clamp.HC_20210422_193236_45419272.h5',
            'distgfs.network_clamp.HC_20210422_193236_53736785.h5',
            'distgfs.network_clamp.HC_20210422_193236_63599789.h5',
],
[
            'distgfs.network_clamp.HCC_20210422_193228_33236209.h5',
            'distgfs.network_clamp.HCC_20210422_193232_12260638.h5',
            'distgfs.network_clamp.HCC_20210422_193232_17609813.h5',
            'distgfs.network_clamp.HCC_20210422_193233_71407528.h5',
            'distgfs.network_clamp.HCC_20210422_193233_92055940.h5',
],
[
            'distgfs.network_clamp.IS_20210422_213208_75940072.h5',
],
[
            'distgfs.network_clamp.MOPP_20210422_213404_31571230.h5',
            'distgfs.network_clamp.MOPP_20210422_213404_45373570.h5',
            'distgfs.network_clamp.MOPP_20210422_213512_68839073.h5',
            'distgfs.network_clamp.MOPP_20210422_213512_85763600.h5',
            'distgfs.network_clamp.MOPP_20210422_214312_29079471.h5',
],
[
            'distgfs.network_clamp.NGFC_20210422_235515_12740157.h5',
            'distgfs.network_clamp.NGFC_20210422_235518_97895890.h5',
            'distgfs.network_clamp.NGFC_20210422_235800_93872787.h5',
            'distgfs.network_clamp.NGFC_20210422_235801_96772370.h5',
            'distgfs.network_clamp.NGFC_20210422_235802_95844113.h5',
],
]
    newint_opt = [
[
    'distgfs.network_clamp.AAC_20210423_205858_36010476.h5',
    'distgfs.network_clamp.AAC_20210423_205858_53499406.h5',
    'distgfs.network_clamp.AAC_20210423_205858_04893658.h5',
    'distgfs.network_clamp.AAC_20210423_205858_49937004.h5',
    'distgfs.network_clamp.AAC_20210423_205858_27137089.h5',
],
[
    'distgfs.network_clamp.BC_20210423_210051_01503844.h5',
    'distgfs.network_clamp.BC_20210423_205910_93454042.h5',
    'distgfs.network_clamp.BC_20210423_205955_74865768.h5',
    'distgfs.network_clamp.BC_20210423_210051_52357252.h5',
    'distgfs.network_clamp.BC_20210423_210803_28135771.h5',
],
[
    'distgfs.network_clamp.HC_20210423_211716_15879716.h5',
    'distgfs.network_clamp.HC_20210423_211528_45419272.h5',
    'distgfs.network_clamp.HC_20210423_211839_28682721.h5',
    'distgfs.network_clamp.HC_20210423_211900_53736785.h5',
    'distgfs.network_clamp.HC_20210423_212004_63599789.h5',
],
[
    'distgfs.network_clamp.HCC_20210423_211101_92055940.h5',
    'distgfs.network_clamp.HCC_20210423_211035_33236209.h5',
    'distgfs.network_clamp.HCC_20210423_211128_17609813.h5',
    'distgfs.network_clamp.HCC_20210423_211127_71407528.h5',
    'distgfs.network_clamp.HCC_20210423_211330_12260638.h5',
],
[
    'distgfs.network_clamp.MOPP_20210423_225834_45373570.h5',
    'distgfs.network_clamp.MOPP_20210423_225855_85763600.h5',
    'distgfs.network_clamp.MOPP_20210423_225756_68839073.h5',
    'distgfs.network_clamp.MOPP_20210423_225756_31571230.h5',
    'distgfs.network_clamp.MOPP_20210423_230054_29079471.h5',
],
[
    'distgfs.network_clamp.NGFC_20210423_230155_95844113.h5',
    'distgfs.network_clamp.NGFC_20210423_230350_93872787.h5',
    'distgfs.network_clamp.NGFC_20210423_231033_96772370.h5',
    'distgfs.network_clamp.NGFC_20210423_230002_12740157.h5',
    'distgfs.network_clamp.NGFC_20210423_230142_97895890.h5',
],
]

#    testfils= [
#[
#        'distgfs.network_clamp.AAC_20210402_163643_27137089.h5',
##        'distgfs.network_clamp.AAC_20210402_163643_36010476.h5',
#],
#]
    
    pairwise = [
#[
#        'distgfs.network_clamp.AAC_20210404_122356_49937004.h5',
#        'distgfs.network_clamp.AAC_20210404_122356_53499406.h5',
#        'distgfs.network_clamp.AAC_20210404_122356_36010476.h5',
#        'distgfs.network_clamp.AAC_20210404_122356_27137089.h5',
#        'distgfs.network_clamp.AAC_20210404_122356_04893658.h5',
#],
#[
#        'distgfs.network_clamp.BC_20210404_122356_74865768.h5',
#        'distgfs.network_clamp.BC_20210404_122356_52357252.h5',
#        'distgfs.network_clamp.BC_20210404_122356_28135771.h5',
#        'distgfs.network_clamp.BC_20210404_122356_93454042.h5',
#        'distgfs.network_clamp.BC_20210405_130348_01503844.h5',
#],
#[
#        'distgfs.network_clamp.HCC_20210404_124301_33236209.h5',
#        'distgfs.network_clamp.HCC_20210404_124358_71407528.h5',
#        'distgfs.network_clamp.HCC_20210404_124301_92055940.h5',
#        'distgfs.network_clamp.HCC_20210404_124358_17609813.h5',
#],
[
        'distgfs.network_clamp.HC_20210404_203555_15879716.h5',
        'distgfs.network_clamp.HC_20210404_204057_45419272.h5',
        'distgfs.network_clamp.HC_20210404_205135_28682721.h5',
        'distgfs.network_clamp.HC_20210404_204313_53736785.h5',
],
#[
#        'distgfs.network_clamp.MOPP_20210404_221007_31571230.h5',
#        'distgfs.network_clamp.MOPP_20210404_222410_85763600.h5',
#        'distgfs.network_clamp.MOPP_20210404_222410_68839073.h5',
#        'distgfs.network_clamp.MOPP_20210404_222410_45373570.h5',
#        'distgfs.network_clamp.MOPP_20210405_051112_29079471.h5',
#],
#[
#        'distgfs.network_clamp.NGFC_20210405_051412_97895890.h5',
#        'distgfs.network_clamp.NGFC_20210404_224448_12740157.h5',
#        'distgfs.network_clamp.NGFC_20210405_052414_93872787.h5',
#        'distgfs.network_clamp.NGFC_20210405_052833_95844113.h5',
#        'distgfs.network_clamp.NGFC_20210405_053112_96772370.h5',
#],
#[
#        'distgfs.network_clamp.IS_20210404_144931_49627038.h5',
#],
]

 #   ModulatedMicroResRange = get_list(fil_dir, ['20210525'])
    ModuMicroBigRange = get_list(fil_dir, ['20210622_163430', '20210622_163845', '20210622_164202', '20210622_164310'])

#    distribute_chores(interneuron_opt, fil_dir, Combined=True, prefix='GlobalExc')
#    distribute_chores(newint_opt, fil_dir, Combined=True, prefix='CombinedCompartments')
#    distribute_chores(pairwise, fil_dir, Combined=True, prefix='Pairwise')
#    distribute_chores(testfils, fil_dir, Combined=True, prefix='Test')

#    distribute_chores(ModulatedMicroResRange, fil_dir, Combined=True, prefix='ModMctResRng')
    distribute_chores(ModuMicroBigRange, fil_dir, Combined=True, prefix='ModMctBigRng')
