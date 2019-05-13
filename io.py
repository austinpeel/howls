import os
import numpy as np
from astropy.io import fits
from common import home
from .constants import __NPIX__


def datapath(subdir=None):
    """Retrieve the system path to the HOWLS convergence maps.

    It is assumed that the kappa maps are fits files located in
    $HOME/Data/HOWLS/...

    Parameters
    ----------
    subdir : str, optional
        Subdirectory(ies) to add after $HOME/Data/HOWLS/

    """

    path = os.path.join(home, 'Data/HOWLS/')
    if subdir is not None:
        path = os.path.join(path, str(subdir))
    return path


def kappapath(map_id, noisy=False, warnings=True):
    """Retrieve path to HOWLS convergence map.

    Parameters
    ----------
    map_id : int
        Map realization number.
    noisy : bool
        Get the noisy version of the map if True.

    """
    option = ['kappa', 'kappa_noise'][noisy]
    filename = option + '_LCDM_zs2_step1_{:03d}.fits'.format(map_id)
    filepath = os.path.join(datapath(), option, filename)

    if warnings:
        if not os.path.exists(filepath):
            print("Warning: could not find {}".format(filepath))
    return filepath


def fetch_kappa(map_id, noisy=False, subsize=None):
    """Retrieve a HOWLS convergence map.

    Parameters
    ----------
    map_id : int
        Map realization number.
    noisy : bool
        Get the noisy version of the map if True.
    subsize : int, optional
        If provided, return a random (subsize, subsize) section of the map.

    """
    kappa = None
    try:
        kappa = fits.getdata(kappapath(map_id, noisy))
        if subsize is not None:
            subsize = int(subsize)
            x1 = np.random.randint(0, __NPIX__ - subsize)
            x2 = x1 + subsize
            y1 = np.random.randint(0, __NPIX__ - subsize)
            y2 = y1 + subsize
            kappa = kappa[x1:x2, y1:y2]
    except:
        pass
    return kappa


def fetch_saved_stats(filename, verbose=False):
    stats = None
    filepath = os.path.join(datapath("output"), filename)
    if not filepath.endswith(".npy"):
        filepath = filepath + ".npy"
    try:
        stats = np.load(filepath)
        if verbose:
            print("Loaded " + filepath)
    except:
        print("Error: could not open " + filepath)
    return stats
