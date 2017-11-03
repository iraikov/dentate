
import sys, time, gc
import numpy as np
from neuroh5.io import read_cell_attributes, read_population_ranges

#  custom data type for type of feature selectivity
selectivity_grid = 0
selectivity_place_field = 1

def generate_trajectory(arena_dimension = 100., velocity = 30., spatial_resolution = 1.):  # cm

    # arena_dimension - minimum distance from origin to boundary (cm)

    x = np.arange(-arena_dimension, arena_dimension, spatial_resolution)
    y = np.arange(-arena_dimension, arena_dimension, spatial_resolution)
    distance = np.insert(np.cumsum(np.sqrt(np.sum([np.diff(x) ** 2., np.diff(y) ** 2.], axis=0))), 0, 0.)
    interp_distance = np.arange(distance[0], distance[-1], spatial_resolution)
    t = interp_distance / velocity * 1000.  # ms
    t_stop = t[-1]
    interp_x = np.interp(interp_distance, distance, x)
    interp_y = np.interp(interp_distance, distance, y)
    d = interp_distance

    return (interp_x, interp_y)



def generate_spatial_ratemap (selectivity_type, features_dict, interp_x, interp_y,
                              grid_peak_rate = 40., place_peak_rate = 40.): # Hz

    response = np.zeros_like(d, dtype='float32')

    a = 0.3
    b = -1.5
    u = lambda ori: (np.cos(ori), np.sin(ori))
    ori_array = 2. * np.pi * np.array([-30., 30., 90.]) / 360.  # rads
    g = lambda x: np.exp(a * (x - b)) - 1.
    scale_factor = g(3.)
    grid_rate = lambda grid_spacing, ori_offset, x_offset, y_offset: \
      lambda x, y: grid_peak_rate / scale_factor * \
      g(np.sum([np.cos(4. * np.pi / np.sqrt(3.) /
                           grid_spacing * np.dot(u(theta - ori_offset), (x - x_offset, y - y_offset)))
                    for theta in ori_array]))

    place_rate = lambda field_width, x_offset, y_offset: \
      lambda x, y: place_peak_rate * np.exp(-((x - x_offset) / (field_width / 3. / np.sqrt(2.))) ** 2.) * \
      np.exp(-((y - y_offset) / (field_width / 3. / np.sqrt(2.))) ** 2.)

    if selectivity_type == selectivity_grid:
        ori_offset = features_dict['Grid Orientation'][0]
        grid_spacing = features_dict['Grid Spacing'][0]
        x_offset = features_dict['X Offset'][0]
        y_offset = features_dict['Y Offset'][0]
        rate = np.vectorize(grid_rate(grid_spacing, ori_offset, x_offset, y_offset))
    elif selectivity_type == selectivity_place_field:
        field_width = features_dict['Field Width'][0]
        x_offset = features_dict['X Offset'][0]
        y_offset = features_dict['Y Offset'][0]
        rate = np.vectorize(place_rate(field_width, x_offset, y_offset))

    response = rate(interp_x, interp_y).astype('float32', copy=False)

    return response



def read_stimulus (comm, stimulus_path, stimulus_namespace, population, verbose=False):
        ratemap_lst = []
        attr_gen = read_cell_attributes(comm, stimulus_path, population, namespace=stimulus_namespace)
        for gid, stimulus_dict in attr_gen:
            rate = stimulus_dict['rate']
            spiketrain = stimulus_dict['spiketrain']
            modulation = stimulus_dict['modulation']
            peak_index = stimulus_dict['peak_index']
            ratemap_lst.append((gid, rate, spiketrain, peak_index))

        ## sort by peak_index
        ratemap_lst.sort(key=lambda item: item[3])

        return ratemap_lst
        


            