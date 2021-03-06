import os, sys, click
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from neuroh5.io import NeuroH5CellAttrGen, read_cell_attributes
from mpi4py import MPI
from pprint import pprint

from dentate.utils import *
from dentate.stimulus import gid2module_dictionary, module2gid_dictionary, fraction_active

script_name = os.path.basename(__file__)


def add_colorbar(img, ax):
    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(img, cax=cax)


def turn_off_xy(axes):
    turn_off_x(axes)
    turn_off_y(axes)


def turn_off_x(axes):
    axes.xaxis.set_ticks_position('none')
    axes.get_xaxis().set_visible(False)


def turn_off_y(axes):
    axes.yaxis.set_ticks_position('none')
    axes.get_yaxis().set_visible(False)


def plot_rate_maps_single_module(cells, plot=False, save=True, **kwargs):

    rate_maps = [] 
    for gid in cells:
        cell = cells[gid]
        rate_maps.append(cell['Rate Map'].reshape(cell['Nx'][0], cell['Ny'][0]))
    rate_maps = np.asarray(rate_maps, dtype='float32')
    
    summed_map = np.sum(rate_maps, axis=0)
    mean_map   = np.mean(rate_maps,axis=0)
    var_map    = np.var(rate_maps, axis=0)
    images = [summed_map, mean_map, var_map]

    ctype = kwargs.get('ctype', 'grid')
    module = kwargs.get('module', 0)

    fig, axes = plt.subplots(3,2)
    for j in range(2):
        for i in range(3):
            img = axes[i,j].imshow(images[i], cmap='inferno')
            plt.colorbar(img, ax=axes[i,j])
            if j == 1:
                img.set_clim(0, np.max(images[i]))
            if i == 0:
                axes[i,j].set_title('Cells: %s. Module: %d. Summed RM' % (ctype, module))
            elif i == 1:
                axes[i,j].set_title('Cells: %s. Module: %d. Mean RM' % (ctype, module))
            else:
                axes[i,j].set_title('Cells: %s. Module: %d. Var RM' % (ctype, module))
    if save:
        plt.savefig('%s-module-%d-ratemap.svg' % (ctype, module), {'format': 'svg'})
    if plot:
        plt.show()


def plot_rate_maps_multiple_modules(module_dictionary, modules, plot=False, save=True, **kwargs):
    assert(len(modules) == 10)
    fig_sum, axes_sum  = plt.subplots(2,5, figsize=[16., 6.])
    fig_mean, axes_mean = plt.subplots(2,5, figsize=[16., 6.])
    fig_var, axes_var  = plt.subplots(2,5, figsize=[16., 6.])

    ctype = kwargs.get('ctype', 'place')
    ax1_count = 0
    ax2_count = 0
    for (i,mod) in enumerate(module_dictionary):
        if i % 2 == 0 and i > 0:
            ax2_count += 1
            ax1_count = 0
        cells = module_dictionary[mod]
        rate_maps = []
        for gid in cells:
            cell = cells[gid]
            rate_maps.append(cell.rate_map.reshape(cell.nx, cell.ny))
        rate_maps = np.asarray(rate_maps, dtype='float32')
        
        summed_map = np.sum(rate_maps, axis=0)
        mean_map   = np.mean(rate_maps, axis=0)
        var_map    = np.var(rate_maps, axis=0)

        img_sum    = axes_sum[ax1_count,ax2_count].imshow(summed_map, cmap='Greys_r', vmin=0., vmax=np.max(summed_map))
        axes_sum[ax1_count, ax2_count].set_title('Cells: %s. Module: %i. Summed RM' % (ctype, mod))
        add_colorbar(img_sum, axes_sum[ax1_count, ax2_count])

        img_mean   = axes_mean[ax1_count, ax2_count].imshow(mean_map, cmap='Greys_r', vmin=0., vmax=np.max(mean_map) )
        axes_mean[ax1_count, ax2_count].set_title('Cells: %s. Module: %i. Mean RM' % (ctype, mod))
        add_colorbar(img_mean, axes_mean[ax1_count, ax2_count])
 
        img_var    = axes_var[ax1_count, ax2_count].imshow(var_map, cmap='Greys_r', vmin=0., vmax=np.max(var_map))
        axes_var[ax1_count, ax2_count].set_title('Cells: %s. Module: %i. Var RM' % (ctype, mod))
        add_colorbar(img_var, axes_var[ax1_count, ax2_count])

        turn_off_xy(axes_sum[ax1_count, ax2_count])
        turn_off_xy(axes_mean[ax1_count, ax2_count])
        turn_off_xy(axes_var[ax1_count, ax2_count])

        ax1_count += 1

    if save:
        fig_sum.savefig('%s-summed-ratemap.svg' % (ctype), format='svg')
        fig_mean.savefig('%s-mean-ratemap.svg' % (ctype), format='svg')
        fig_var.savefig('%s-var-ratemap.svg' % (ctype), format='svg')
    if plot:
        plt.show()
        

def plot_xy_offsets_single_module(cells, plot=False, save=True, **kwargs):

    pop_xy_offsets = []
    for gid in cells:
        cell = cells[gid]
        if 'X Offset Scaled' not in cell or 'Y Offset Scaled' not in cell:
            offsets = list(zip(cell['X Offset'], cell['Y Offset']))
        else:

            offsets = list(zip(cell['X Offset Scaled'], cell['Y Offset Scaled']))
        for (x_offset, y_offset) in offsets:
            pop_xy_offsets.append((x_offset, y_offset))
    pop_xy_offsets = np.asarray(pop_xy_offsets, dtype='float32')

    ctype  = kwargs.get('ctype', 'grid')
    module = kwargs.get('module', 0)

    plt.figure()
    plt.scatter(pop_xy_offsets[:,0], pop_xy_offsets[:,1])
    plt.title('Cells: %s. Module: %d. xy-offsets' % (ctype, module))
    
    if save:
        plt.savefig('%s-module-%d-xyoffsets.png' % (ctype, module))
    if plot:
        plt.show()


def plot_xy_offsets_multiple_modules(modules_dictionary, modules, plot=False, save=False, **kwargs):
    assert(len(modules) == 10)

    fig, axes = plt.subplots(2,5, figsize=[20., 12.])
    ctype     = kwargs.get('ctype', 'place')
    
    ax1_count = 0
    ax2_count = 0
    for (i, module) in enumerate(modules_dictionary):
        if i % 2 == 0 and i > 0:
            ax2_count += 1
            ax1_count = 0
        pop_xy_offsets = []
        cells = modules_dictionary[module]
        for gid in cells:
            cell = cells[gid]
            offsets = list(zip(cell.x_offset, cell.y_offset))
            for (x_offset, y_offset) in offsets:
                pop_xy_offsets.append((x_offset, y_offset))
        pop_xy_offsets = np.asarray(pop_xy_offsets, dtype='float32')
        
        axes[ax1_count, ax2_count].scatter(pop_xy_offsets[:,0], pop_xy_offsets[:,1])
        axes[ax1_count, ax2_count].set_title('Cells: %s. Module: %i' % (ctype, module))
        ax1_count += 1
    if save:
        fig.savefig('%s-xy-offsets.svg' % (ctype), format='svg')
    if plot:
        plt.show()


def plot_fraction_active_single_module(cells, nxny, plot=False, save=True, **kwargs):
    factive = fraction_active(cells, 2.)
    fraction_active_img = np.zeros(nxny)
    for (i,j) in factive:
        fraction_active_img[i,j] = factive[(i,j)]

    ctype  = kwargs.get('ctype', 'grid')
    module = kwargs.get('module',0)
    target = kwargs.get('target', 0.05)
   
    fig, axes = plt.subplots(1,2)

    img = axes[0].imshow(fraction_active_img, cmap='inferno')
    plt.colorbar(img, ax=axes[0])
    #img.set_clim(0, 1.0)
    axes[0].set_title('Fraction Active')

    img = axes[1].imshow(fraction_active_img - target, cmap='inferno')
    plt.colorbar(img, ax=axes[1])
    #img.set_clim(0, 1.0)
    axes[1].set_title('Fraction active distance from target')
 
    if save:
        plt.savefig('%s-module-%d-fractionactive.png' % (ctype, module))
    if plot:
        plt.show()


def plot_fraction_active_multiple_modules(modules_dictionary, modules, plot=False, save=True, **kwargs):
    assert(len(modules) == 10)
    ctype     = kwargs.get('ctype', 'place')
    fig, axes = plt.subplots(2,5, figsize=[16., 6.])

    ax_count1 = 0
    ax_count2 = 0
    for (i, module) in enumerate(modules_dictionary):
        if i % 2 == 0 and i > 0:
            ax_count2 += 1
            ax_count1 = 0

        cells   = modules_dictionary[module]
        factive = fraction_active(cells, 2.)
        nx = np.max([x[0] for x in list(factive.keys())]) + 1
        ny = np.max([x[1] for x in list(factive.keys())]) + 1
        fraction_active_img = np.zeros((nx, ny))
        for (i,j) in factive:
            fraction_active_img[i,j] = factive[(i,j)]
        img = axes[ax_count1, ax_count2].imshow(fraction_active_img, cmap='inferno')
        add_colorbar(img, axes[ax_count1, ax_count2])
        axes[ax_count1, ax_count2].set_title('Cells: %s. Module: %i. FA' % (ctype, module))
        turn_off_xy(axes[ax_count1, ax_count2])
        ax_count1 += 1

    if save:
        fig.savefig('%s-fraction-active.svg' % (ctype), format='svg')
    if plot:
        plt.show()


def plot_rate_histogram_single_module(cells, plot=False, save=True, **kwargs):
    ctype  = kwargs.get('ctype', 'grid')
    module = kwargs.get('module', 0)
    peak_firing_rate = kwargs.get('peak firing rate', 20.0)
    bins = kwargs.get('nbins', 40)

    rate_maps = []
    for gid in cells:
        cell = cells[gid]
        nx, ny = cell['Nx'][0], cell['Ny'][0]
        rate_maps.append(cell['Rate Map'].reshape(nx, ny))
    rate_maps = np.asarray(rate_maps, dtype='float32')

    N, nx, ny = rate_maps.shape
    weights = np.ones(N) / float(N)
    hists, edges_list = [], []
    for i in range(nx):
        for j in range(ny):
            hist, edges = np.histogram(rate_maps[:,i,j], bins=bins, range=(0.0, peak_firing_rate), weights=weights)
            hists.append(hist)
            edges_list.append(edges)
    hists = np.asarray(hists)
    hist_mean = np.mean(hists, axis=0)
    hist_std = np.sqrt(np.var(hists, axis=0))

    fig, axes = plt.subplots(2, 1)
    axes[0].bar(edges_list[0][1:], hist_mean, alpha=0.5, log=False, yerr=hist_std)
    axes[0].set_title('Firing rate histogram')
    axes[0].set_ylabel('Probability')

    axes[1].bar(edges_list[0][1:], hist_mean, alpha=0.5, log=True, yerr=hist_std)
    axes[1].set_ylabel('Log Probability')
    axes[1].set_xlabel('Firing rate (Hz')
 
    if save:
        plt.savefig('%s-module-%d-rate-histogram.png' % (ctype, module))
    if plot:
        plt.show()


def plot_rate_histogram_multiple_modules(module_dictionary, modules, plot=False, save=True, **kwargs):
    ctype   = kwargs.get('ctype', 'place')
    bins    = kwargs.get('bins', 40)
    peak_firing_rate = kwargs.get('peak firing rate', 20.0)

    fig, axes = plt.subplots(2,5, figsize=[20., 12.])
    ax1_count = 0
    ax2_count = 0

    for (i, module) in enumerate(module_dictionary):
        if i % 2 == 0 and i > 0:
            ax1_count = 0
            ax2_count += 1

        cells = module_dictionary[module]
        rate_maps = []
        for gid in cells:
            cell   = cells[gid]
            nx, ny = cell.nx, cell.ny
            rate_maps.append(cell.rate_map.reshape(nx, ny))
        rate_maps = np.asarray(rate_maps, dtype='float32')
        
        N, nx, ny = rate_maps.shape
        weights = np.ones(N) / float(N)
        hists, edges_list = [], []
        for i in range(nx):
            for j in range(ny):
                hist, edges = np.histogram(rate_maps[:,i,j], bins=bins, range=(0.0, peak_firing_rate), weights=weights)
                hists.append(hist)
                edges_list.append(edges)
        hists = np.asarray(hists)
        hist_mean = np.mean(hists, axis=0)
        hist_std  = np.std(hists, axis=0)
        
        axes[ax1_count, ax2_count].bar(edges_list[0][1:], hist_mean, alpha=0.5, log=True, yerr=hist_std)
        axes[ax1_count, ax2_count].set_ylabel('log Probability')
        axes[ax1_count, ax2_count].set_xlabel('firing rate (hz)')
        axes[ax1_count, ax2_count].set_title('Cells: %s. Module: %i' % (ctype, module))
        ax1_count += 1
    if save: 
        fig.savefig('%s-rate-histogram.svg' % (ctype), format='svg')
    if plot:
        plt.show()


def plot_lambda_activity_histograms(module_dictionary, modules, save=False, plot=False, threshold=2.0, **kwargs):
    
    fig, axes   = plt.subplots(2)
    ctype     = kwargs.get('ctype', 'grid')
    nx, ny    = 20, 20
    nxx, nyy  = np.meshgrid(np.arange(nx), np.arange(ny))
    coords    = list(zip(nxx.reshape(-1,), nyy.reshape(-1,)))
    active_lambda_maps = {(i,j): [] for (i,j) in coords}
    for (i, module) in enumerate(module_dictionary):
        cells = module_dictionary[module]
        for gid in cells:
            cell = cells[gid]
            rate_map = cell['Rate Map'].reshape(cell['Nx'][0], cell['Ny'][0])
            module = cell['Module'][0]
            for x in range(nx):
                for y in range(ny):
                    response = rate_map[(x,y)]
                    if response >= threshold:
                        active_lambda_maps[(x,y)].append(module)

    hist_lst, edges_lst = [], []
    for (i,j) in active_lambda_maps:
        pos_activity = active_lambda_maps[(i, j)]
        hist, edges = np.histogram(pos_activity,bins=10)
        hist_lst.append(hist)
        edges_lst.append(edges)
    hists = np.asarray(hist_lst)

    hist_mean = np.mean(hists,axis=0)
    hist_std  = np.std(hists,axis=0)
    hist_mean_normalized = old_div(hist_mean, [len(module_dictionary[i]) for i in module_dictionary])
    hist_std_normalized = old_div(hist_std, [len(module_dictionary[i]) for i in module_dictionary])

    axes[0].bar(edges_lst[0][1:], hist_mean, alpha=0.5, yerr=hist_std)
    axes[0].set_title('Activity per module')

    axes[1].bar(edges_lst[0][1:], hist_mean_normalized, alpha=0.5, yerr=hist_std_normalized)
    axes[1].set_title('Normalized acitivity per module')

    if save:
        fig.savefig('%s-module-activity.svg' % (ctype), format='svg')
    if plot:
        plt.show()


def read_population_storage(file_path):
    try:
        storage = PopulationStorage(file_path=file_path)
    except:
        raise Exception('Error occured when loading PopulationStorage object')
    return storage


def plot_population_storage(storage):
    try:
        storage.plot()
    except:
        raise Exception('Error occured when loading PopulationStorage object')


def plot_population_input(storage, bounds=None):
    individuals = np.asarray([generation for generation in storage.history]).flatten()
    inputs      = np.asarray([individual.x for individual in individuals])
    frac_active = np.asarray([individual.features[0] for individual in individuals])

    if bounds is not None:
        fig, (ax1, ax2) = plt.subplots(2,1)
        valid_inputs      = inputs[(frac_active >= bounds[0]) & (frac_active <= bounds[1])]
        valid_frac_active = frac_active[(frac_active >= bounds[0]) & (frac_active <= bounds[1])]
    else:
        fig, ax1 = plt.subplots(1,1)

    full_scatter = ax1.scatter(inputs[:,0], inputs[:,1], c=frac_active, cmap='viridis')
    plt.colorbar(full_scatter, ax=ax1)
    ax1.set_xlabel('p_inactive')
    ax1.set_ylabel('p_r')

    if bounds is not None:
        bounded_scatter = ax2.scatter(valid_inputs[:,0], valid_inputs[:,1], c=valid_frac_active, cmap='viridis')
        plt.colorbar(bounded_scatter, ax=ax2) 
        ax2.set_xlabel('p_inactive')
        ax2.set_ylabel('p_r')


def plot_group(module_dictionary, modules, plot=False, save=False, **kwargs):
    plot_rate_maps_multiple_modules(module_dictionary, modules, plot=False, save=save, **kwargs)
    plot_fraction_active_multiple_modules(module_dictionary, modules, plot=False, save=save, **kwargs)
    plot_xy_offsets_multiple_modules(module_dictionary, modules, plot=False, save=save, **kwargs)
    plot_rate_histogram_multiple_modules(module_dictionary, modules, plot=plot, save=save, **kwargs)


@click.command()
@click.option("--features-path", required=True, type=click.Path(exists=True, file_okay=True, dir_okay=False))
@click.option("--cell-type", required=True, type=str)
@click.option("--arena-id", type=str, default='A')
@click.option("--show-fig", is_flag=True)
@click.option("--save-fig", is_flag=True)
def main(features_path, cell_type, arena_id, show_fig, save_fig):

    font = {'family': 'normal', 'weight': 'bold', 'size': 6}
    matplotlib.rc('font', **font)

    comm = MPI.COMM_WORLD
    modules = np.arange(10) + 1

    if cell_type == 'grid':
        mpp_grid = read_cell_attributes(features_path, 'MPP', 'Grid Selectivity %s' % arena_id)
        cells_modules_dictionary = gid2module_dictionary([mpp_grid], modules)
    elif cell_type == 'place':
        lpp_place = read_cell_attributes(features_path, 'LPP', 'Place Selectivity %s' % arena_id)
        mpp_place = read_cell_attributes(features_path, 'MPP', 'Place Selectivity %s' % arena_id )
        cells_modules_dictionary = gid2module_dictionary([mpp_place, lpp_place], modules)
    elif cell_type == 'both':
        lpp_place = read_cell_attributes(features_path, 'LPP', 'Place Selectivity %s' % arena_id)
        mpp_place = read_cell_attributes(features_path, 'MPP', 'Place Selectivity %s' % arena_id)
        mpp_grid  = read_cell_attributes(features_path, 'MPP', 'Grid Selectivity %s' % arena_id)
        cells_modules_dictionary = gid2module_dictionary([mpp_grid, mpp_place, lpp_place], modules)

    kwargs = {'ctype': cell_type}
    plot_group(cells_modules_dictionary, modules, plot=show_fig, save=save_fig, **kwargs)


if __name__ == '__main__':
    main(args=sys.argv[(list_find(lambda s: s.find(script_name) != -1, sys.argv) + 1):], standalone_mode=False)
