import numpy as np
from scipy.spatial.distance import pdist, squareform


def nearest_neighbor_distance(X, Y, Z):
    """
    Determines the nearest neighbor distance (center of mass distance)
    from an array of centers of mass at positions X, Y, and Z.

    :param X: X-coordinates of the centers of mass
    :param Y: Y-coordinates of the centers of mass
    :param Z: Z-coordinates of the centers of mass
    :return: Distance to the nearest neighboring object (pore).
    :rtype: numpy.ndarray
    """
    xyz = np.array((X, Y, Z)).T
    distances = squareform(pdist(xyz))
    distances[np.diag_indices_from(distances)] = distances.max()
    return distances.min(axis=1)


def sphere_equivalent_diameter(volume):
    """
    Returns the sphere-equivalent diameter given a vector of volumes.

    :param volume: Volumes of the objects.
    :return: Sphere-equivalent diameters
    """
    volume = np.asarray(volume)
    return (6*volume/np.pi)**(1/3)


def median_pore_diameter(volume):
    """
    Calculates the median pore diameter from the pore volumes.

    :param volume: Measured pore volumes.
    :return: Median pore diameters.
    """
    return np.median(sphere_equivalent_diameter(volume))


def median_pore_spacing(X, Y, Z):
    """
    Calculates the median pore spacing between pores at the specified
    centers of mass.

    :param X: X-coordinates of the centers of mass.
    :param Y: Y-coordinates of the centers of mass.
    :param Z: Z-coordinates of the centers of mass
    :return: The median pore spacing.
    """
    return np.median(nearest_neighbor_distance(X, Y, Z))


def mean_pore_spacing(X, Y, Z):
    """
    Calculates the mean pore spacing for pores at the specified centers
    of mass.

    :param X: X-coordinates of the centers of mass.
    :param Y: Y-coordinates of the centers of mass
    :param Z: Z-coordinates of the centers of mass
    :return: The mean pore spacing.
    """
    return np.mean(nearest_neighbor_distance(X, Y, Z))

def max_pore_diameter(volume):
    """
    Calculates the max pore diameter for pores at the specified centers of mass.

    :param X: X-coordinates of the centers of mass
    :param Y: Y-coordinates of the centers of mass
    :param Z: Z-coordinates of the centers of mass
    :return: The max pore diameter
    """

    return np.max(sphere_equivalent_diameter(volume))