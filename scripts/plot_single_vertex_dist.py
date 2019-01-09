
import sys, os, gc, math
from mpi4py import MPI
import click
import dentate
from dentate import utils, plot
from dentate.env import Env

script_name = os.path.basename(__file__)

@click.command()
@click.option("--config", required=True, type=str)
@click.option("--config-prefix", required=True, type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option("--connectivity-path", '-p', required=True, type=click.Path())
@click.option("--coords-path", '-c', required=True, type=click.Path())
@click.option("--distances-namespace", '-t', type=str, default='Arc Distances')
@click.option("--destination-gid", '-g', type=int)
@click.option("--destination", '-d', type=str)
@click.option("--source", '-s', type=str)
@click.option("--extent-type", type=str, default='local')
@click.option("--bin-size", type=float, default=20.0)
@click.option("--font-size", type=float, default=14)
@click.option("--verbose", "-v", type=bool, default=False, is_flag=True)
def main(config, config_prefix, connectivity_path, coords_path, distances_namespace, destination_gid, destination, source, extent_type, bin_size, font_size, verbose):

    utils.config_logging(verbose)
    logger = utils.get_script_logger(os.path.basename(script_name))

    env = Env(config_file=config, config_prefix=config_prefix)

    plot.plot_single_vertex_dist (env, connectivity_path, coords_path, distances_namespace, \
                                  destination_gid, destination, source, \
                                  extent_type=extent_type, bin_size=bin_size, fontSize=font_size, \
                                  saveFig=True)
    

if __name__ == '__main__':
    main(args=sys.argv[(utils.list_find(lambda x: os.path.basename(x) == script_name, sys.argv)+1):])
