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


def kappapath(step, Om, map_id, noisy, warnings=True):
    """Retrieve path to HOWLS convergence map.

    Parameters
    ----------
    step : int
        Step number of the HOWLS project. Options so far are only {1, 2}.
    Om : float or None, optional
        Omega_m density parameter of the simulation (used in Step 2).
        If None and step is 2, use the fiducial model.
    map_id : int
        Map realization number.
    noisy : bool
        Get the noisy version of the map if True.

    """
    # Build data path to kappa maps
    if step == 1:
        option = ['kappa', 'kappa_noise'][noisy]
        filename = option + '_LCDM_zs2_step1_{:03d}.fits'.format(map_id)
        filepath = os.path.join(datapath('step1'), option, filename)
    elif step == 2:
        filename = 'kappa_noise_LCDM_'
        if Om == None:
            dir = 'LCDM_fiducial'
            opt = ''
        else:
            dir = 'LCDM_Om_{}'.format(Om)
            opt = 'Om_{}_'.format(Om)
        filename += opt + 'zs2_step2_{:03d}.fits'.format(map_id)
        filepath = os.path.join(datapath('step2'), dir, filename)
    else:
        if warnings:
            print("Warning: invalid step number. Using 2.")
            step = 2
        return kappapath(step=2, Om=Om, map_id=map_id, noisy=noisy)

    if warnings:
        if not os.path.exists(filepath):
            print("Warning: could not find {}".format(filepath))
    return filepath


def fetch_kappa(step, Om=None, map_id=0, noisy=False, subsize=None):
    """Retrieve a HOWLS convergence map.

    Parameters
    ----------
    step : int
        Step number of the HOWLS project. Options so far are only {1, 2}.
    Om : float or None, optional
        Omega_m density parameter of the simulation (used in Step 2).
        If None and step is 2, use the fiducial model.
    map_id : int
        Map realization number.
    noisy : bool
        Get the noisy version of the map if True (used in Step 1).
    subsize : int, optional
        If provided, return a random (subsize, subsize) section of the map.

    """
    kappa = None
    try:
        kappa = fits.getdata(kappapath(step, Om, map_id, noisy))
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
