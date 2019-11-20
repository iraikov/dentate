
import os, sys
import click
import dentate
from dentate import plot
from dentate import utils

script_name = os.path.basename(__file__)

@click.command()
@click.option("--config-path", '-c', required=True, type=click.Path())
@click.option("--input-path", '-p', required=True, type=click.Path())
@click.option("--t-max", type=float)
@click.option("--t-min", type=float)
@click.option("--window-size", type=int, default=4096)
@click.option("--overlap", type=float, default=0.9)
@click.option("--frequency-range", type=(float, float), default=(0., 500.))
@click.option("--font-size", type=float, default=14)
@click.option("--verbose", "-v", type=bool, default=False, is_flag=True)
def main(config_path, input_path, t_max, t_min, window_size, overlap, frequency_range, font_size, verbose):

    utils.config_logging(verbose)

    if t_max is None:
        time_range = None
    else:
        if t_min is None:
            time_range = [0.0, t_max]
        else:
            time_range = [t_min, t_max]

    plot.plot_lfp_spectrogram (config_path, input_path, time_range=time_range, \
                   window_size=window_size, overlap=overlap, frequency_range=frequency_range, \
                   fontSize=font_size, saveFig=True)
    

if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == script_name, sys.argv)+1):])
