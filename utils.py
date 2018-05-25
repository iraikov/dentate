from collections import defaultdict, Iterable, namedtuple
import sys, os.path, string, time, gc, math, datetime, numbers, itertools
import copy, pprint, logging
import numpy as np
import yaml
import pprint
import logging


class IncludeLoader(yaml.Loader):

    def __init__(self, stream):

        self._root = os.path.split(stream.name)[0]

        yaml.Loader.__init__(self, stream)

    def include(self, node):

        filename = os.path.join(self._root, self.construct_scalar(node))

        with open(filename, 'r') as f:
            return yaml.load(f, IncludeLoader)

IncludeLoader.add_constructor('!include', IncludeLoader.include)


def write_to_yaml(file_path, data, convert_scalars=False):
    """

    :param file_path: str (should end in '.yaml')
    :param data: dict
    :param convert_scalars: bool
    :return:
    """
    import yaml
    with open(file_path, 'w') as outfile:
        if convert_scalars:
            data = nested_convert_scalars(data)
        yaml.dump(data, outfile, default_flow_style=False)


def read_from_yaml(file_path):
    """

    :param file_path: str (should end in '.yaml')
    :return:
    """
    import yaml
    if os.path.isfile(file_path):
        with open(file_path, 'r') as stream:
            data = yaml.load(stream)
        return data
    else:
        raise Exception('File: {} does not exist.'.format(file_path))


def nested_convert_scalars(data):
    """
    Crawls a nested dictionary, and converts any scalar objects from numpy types to python types.
    :param data: dict
    :return: dict
    """
    if isinstance(data, dict):
        for key in data:
            data[key] = nested_convert_scalars(data[key])
    elif isinstance(data, Iterable) and not isinstance(data, (str, tuple)):
        for i in xrange(len(data)):
            data[i] = nested_convert_scalars(data[i])
    elif hasattr(data, 'item'):
        try:
            data = np.asscalar(data)
        except TypeError:
            pass
    return data


def list_find(f, lst):
    """

    :param f:
    :param lst:
    :return:
    """
    i=0
    for x in lst:
        if f(x):
            return i
        else:
            i=i+1
    return None


def list_argsort(f, seq):
    """
    http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3383106#3383106
    lambda version by Tony Veijalainen
    :param f:
    :param seq:
    :return:
    """
    return [i for i,x in sorted(enumerate(seq), key = lambda x: f(x[1]))]

def list_index (element, lst):
    try:
        index_element = lst.index(element)
        return index_element
    except ValueError:
        return None


def random_choice_w_replacement(ranstream,n,p):
    return ranstream.multinomial(n,p.ravel())


def make_random_clusters(centers, n_samples_per_center, n_features=2, cluster_std=1.0,
                         center_ids=None, center_box=(-10.0, 10.0), random_seed=None):
    """Generate isotropic Gaussian blobs for clustering.

    Parameters
    ----------
    centers : int or array of shape [n_centers, n_features]
        The number of centers to generate, or the fixed center locations.
    n_samples_per_center : int array
        Number of points for each cluster.
    n_features : int, optional (default=2)
        The number of features for each sample.
    cluster_std : float or sequence of floats, optional (default=1.0)
        The standard deviation of the clusters.
    center_ids : array of integer center ids, if None then centers will be numbered 0 .. n_centers-1
    center_box : pair of floats (min, max), optional (default=(-10.0, 10.0))
        The bounding box for each cluster center when centers are
        generated at random.
    random_seed : int or None, optional (default=None)
        If int, random_seed is the seed used by the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.
    Returns
    -------
    X : array of shape [n_samples, n_features]
        The generated samples.
    y : array of shape [n_samples]
        The integer labels for cluster membership of each sample.
    Examples
    --------
    >>> X, y = make_random_clusters (centers=6, n_samples_per_center=np.array([1,3,10,15,7,9]), n_features=1, \
                                     center_ids=np.array([10,13,21,25,27,29]).reshape(-1,1), cluster_std=1.0, \
                                     center_box=(-10.0, 10.0))
    >>> print(X.shape)
    (45, 1)
    >>> y
    array([10, 13, 13, 13, ..., 29, 29, 29])
    """
    rng = np.random.RandomState(random_seed)

    if isinstance(centers, numbers.Integral):
        centers = np.sort(rng.uniform(center_box[0], center_box[1], \
                                      size=(centers, n_features)), axis=0)
    else:
        assert(isinstance(centers, np.ndarray))
        n_features = centers.shape[1]

    if center_ids is None:
        center_ids = np.arange(0, centers.shape[0])
        
    if isinstance(cluster_std, numbers.Real):
        cluster_std = np.ones(len(centers)) * cluster_std

    X = []
    y = []

    n_centers = centers.shape[0]

    for i, (cid, n, std) in enumerate(itertools.izip(center_ids, n_samples_per_center, cluster_std)):
        if n > 0:
            X.append(centers[i] + rng.normal(scale=std, size=(n, n_features)))
            y += [cid] * n

    X = np.concatenate(X)
    y = np.array(y)

    return X, y


def random_clustered_shuffle(centers, n_samples_per_center, center_ids=None, cluster_std=1.0, center_box=(-1.0, 1.0), random_seed=None):
    """Generates a Gaussian random clustering given a number of cluster
    centers, samples per each center, optional integer center ids, and
    cluster standard deviation.

    Parameters
    ----------
    centers : int or array of shape [n_centers]
        The number of centers to generate, or the fixed center locations.
    n_samples_per_center : int array
        Number of points for each cluster.
    cluster_std : float or sequence of floats, optional (default=1.0)
        The standard deviation of the clusters.
    center_ids : array of integer center ids, if None then centers will be numbered 0 .. n_centers-1
    random_seed : int or None, optional (default=None)
        If int, random_seed is the seed used by the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    >>> x = random_clustered_shuffle(centers=6,center_ids=np.array([10,13,21,25,27,29]).reshape(-1,1), \
                                     n_samples_per_center=np.array([1,3,10,15,7,9]))
    >>> array([10, 13, 13, 25, 13, 29, 21, 25, 27, 21, 27, 29, 25, 25, 25, 21, 29,
               27, 25, 21, 29, 25, 25, 25, 25, 29, 21, 25, 21, 29, 29, 29, 21, 25,
               29, 21, 27, 27, 21, 27, 25, 21, 25, 27, 25])
    """

    if isinstance(centers, numbers.Integral):
        n_centers = centers
    else:
        assert(isinstance(centers, np.ndarray))
        n_centers = len(centers)
    
    X, y = make_random_clusters (centers, n_samples_per_center, n_features=1, \
                                 center_ids=center_ids, cluster_std=cluster_std, center_box=center_box, \
                                 random_seed=random_seed)
    s = np.argsort(X,axis=0).ravel()
    return y[s].ravel()

