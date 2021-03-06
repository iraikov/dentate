
import gc
import os
import sys

import click
import dentate
from dentate import plot
from dentate import utils
from mpi4py import MPI

script_name = os.path.basename(__file__)

@click.command()
@click.option("--spike-events-path", '-p', required=True, type=click.Path())
@click.option("--spike-events-namespace", '-n', type=str, default='Spike Events')
@click.option("--populations", '-i', type=str, multiple=True)
@click.option("--include-artificial/--exclude-artificial", type=bool, default=True, is_flag=True)
@click.option("--bin-size", type=float, default=50.0)
@click.option("--t-variable", type=str, default='t')
@click.option("--t-max", type=float)
@click.option("--t-min", type=float)
@click.option("--quantity", type=str, default='rate')
@click.option("--graph-type", type=str, default='point')
@click.option("--font-size", type=float, default=14)
@click.option("--overlay", type=bool, default=False, is_flag=True)
@click.option("--unit", type=str, default='cell')
@click.option("--verbose", "-v", type=bool, default=False, is_flag=True)
def main(spike_events_path, spike_events_namespace, populations, include_artificial, bin_size, t_variable, t_max, t_min, quantity, graph_type, font_size, overlay, unit, verbose):

    utils.config_logging(verbose)
    
    if t_max is None:
        time_range = None
    else:
        if t_min is None:
            time_range = [0.0, t_max]
        else:
            time_range = [t_min, t_max]

    if not populations:
        populations = ['eachPop']

    if unit == 'cell':
        plot.plot_spike_distribution_per_cell (spike_events_path, spike_events_namespace, include=populations, include_artificial=include_artificial,
                                               time_variable=t_variable, time_range=time_range, quantity=quantity, fontSize=font_size,
                                               graph_type = graph_type, overlay=overlay, saveFig=True)
    elif unit == 'time':
        plot.plot_spike_distribution_per_time (spike_events_path, spike_events_namespace, include=populations, include_artificial=include_artificial,
                                               time_variable=t_variable, time_range=time_range, time_bin_size=bin_size, quantity=quantity,
                                               fontSize=font_size, overlay=overlay, saveFig=True)

if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == script_name, sys.argv)+1):])
