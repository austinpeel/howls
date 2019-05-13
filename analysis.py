from __future__ import print_function
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from .io import fetch_kappa, datapath
from .constants import __MAP_SIZE__, __PIX_AREA__
from .conversions import j2arcmin, j2pix
from wltools.transforms import mr_transform
from wltools.peaks import find_peaks2d
from wltools.plotting import plot2d
from wltools.utils import wavelet_norm, print_time
from scipy.stats import skew, moment
from scipy.stats import kurtosis as kurt
from lenstools import ConvergenceMap


def compute_peaks():
    """Compute peak counts as a function of scale and SNR.

    Definitions
    -----------
    Scales : wavelet (starlet) j = 1, ..., 7
             corresponding to aperture filter theta = 0.29, ..., 18.8 arcmin

    SNR    : 0.0 <= nu <= 5.0 with dnu = 0.1
             This gives 50 bins in total.

    """
    start_time = time.time()
    # Parameters
    nmaps = 256
    jmin = 3
    jmax = 7
    nbins = 100
    bin_edges, dnu = np.linspace(0, 10, nbins + 1, retstep=True)
    # Add final high SNR bin
    bin_edges = np.array(list(bin_edges) + [1e5])
    nbins = nbins + 1
    peaks_clean = np.zeros((nmaps, jmax - jmin + 1, nbins))
    peaks_noisy = np.zeros((nmaps, jmax - jmin + 1, nbins))

    # True kappa noise
    sigma_kappa = 0.3 / np.sqrt(30 * __PIX_AREA__.value)

    print("Setup")
    print("-----")
    print("Maps")
    print("  {}".format(nmaps))
    print("Scales")
    for ii in range(jmin, jmax + 1):
        print("  {} : {:.2f}".format(ii, j2arcmin(ii)))
    print("SNR")
    print("  {} bins".format(nbins))
    print("  {:.1f} to {:.1f} with delta = {:.3f}".format(bin_edges[0],
          bin_edges[-2], dnu))
    
    print("\nDoing clean maps")
    print("----------------")
    for ii in range(nmaps):
        print("\r{}".format(ii + 1), end='')
        kappa = fetch_kappa(map_id=ii, noisy=False)
        wt = mr_transform(kappa, nscales=jmax)
        for jj, w in zip(range(jmin - 1, jmax), wt[(jmin - 1):-1]):
            b = int(j2pix(jj + 1)) # border pixels to exclude
            wcrop = w[b:-b, b:-b]
            x, y, h = find_peaks2d(wcrop / np.std(w))
            hist, bins = np.histogram(h, bin_edges)
            peaks_clean[ii, jj - (jmin - 1)] = hist
        sys.stdout.flush()

    print("\n\nDoing noisy maps")
    print("----------------")
    for ii in range(nmaps):
        print("\r{}".format(ii + 1), end='')
        kappa = fetch_kappa(map_id=ii, noisy=True)
        wt = mr_transform(kappa, nscales=jmax)
        sigma1 = np.std(wt[0])
        for jj, w in zip(range(jmin - 1, jmax), wt[(jmin - 1):-1]):
            b = int(j2pix(jj + 1)) # border pixels to exclude
            wcrop = w[b:-b, b:-b]
            # Estimate the true noise level of this scale
            # sigma = sigma1 * wavelet_norm(jj + 1) / wavelet_norm(1)
            sigma = sigma_kappa * wavelet_norm(jj + 1) # True sigma
            x, y, h = find_peaks2d(wcrop / sigma)
            hist, bins = np.histogram(h, bin_edges)
            peaks_noisy[ii, jj - (jmin - 1)] = hist
        sys.stdout.flush()

    # Save the result
    outpath = datapath("output")
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    np.save(os.path.join(outpath, "peaks_clean"), peaks_clean)
    np.save(os.path.join(outpath, "peaks_noisy"), peaks_noisy)

    print("\n")
    print_time(time.time() - start_time)


def compute_moments():
    """Compute moments (2, 3, 4) as a function of scale.

    Definitions
    -----------
    Scales : wavelet (starlet) j = 1, ..., 7
             corresponding to aperture filter theta = 0.29, ..., 18.8 arcmin

    Moments : variance - Second central moment
              skewness - Third standardized central moment
              kurtosis - Fourth standardized central moment (Fisher's def.)
              mu3      - Third central moment (not standardized)
              mu4      - Fourth central moment (not standardized)

    Notes
    -----
    np.var(ddof=1)       agrees with Mathematica Variance[] function.
    scipy.stats.skew     agrees with Mathematica Skewness[] function.
    scipy.stats.kurtosis agrees with Mathematica Kurtosis[] function.
    scipy.stats.moment   agrees with Mathematica CentralMoment[] function.

    """
    start_time = time.time()
    # Parameters
    nmaps = 256
    jmin = 3
    jmax = 7
    var_clean = np.zeros((nmaps, jmax - jmin + 1))
    skew_clean = np.zeros_like(var_clean)
    kurt_clean = np.zeros_like(var_clean)
    mu3_clean = np.zeros_like(var_clean)
    mu4_clean = np.zeros_like(var_clean)
    var_noisy = np.zeros_like(var_clean)
    skew_noisy = np.zeros_like(var_clean)
    kurt_noisy = np.zeros_like(var_clean)
    mu3_noisy = np.zeros_like(var_clean)
    mu4_noisy = np.zeros_like(var_clean)

    print("Setup")
    print("-----")
    print("Maps")
    print("  {}".format(nmaps))
    print("Scales")
    for ii in range(jmin, jmax + 1):
        print("  {} : {:.2f}".format(ii, j2arcmin(ii)))

    print("\nDoing clean maps")
    print("----------------")
    for ii in range(nmaps):
        print("\r{}".format(ii + 1), end='')
        kappa = fetch_kappa(map_id=ii, noisy=False)
        wt = mr_transform(kappa, nscales=jmax)
        for jj, w in zip(range(jmin - 1, jmax), wt[(jmin - 1):-1]):
            b = int(j2pix(jj + 1)) # border pixels to exclude
            wcrop = w[b:-b, b:-b]
            var_clean[ii, jj - (jmin - 1)] = np.var(wcrop, axis=None, ddof=1)
            skew_clean[ii, jj - (jmin - 1)] = skew(wcrop, axis=None)
            kurt_clean[ii, jj - (jmin - 1)] = kurt(wcrop, axis=None)
            mu3_clean[ii, jj - (jmin - 1)] = moment(wcrop, 3, axis=None)
            mu4_clean[ii, jj - (jmin - 1)] = moment(wcrop, 4, axis=None)
        sys.stdout.flush()

    print("\n\nDoing noisy maps")
    print("----------------")
    for ii in range(nmaps):
        print("\r{}".format(ii + 1), end='')
        kappa = fetch_kappa(map_id=ii, noisy=True)
        wt = mr_transform(kappa, nscales=jmax)
        for jj, w in zip(range(jmin - 1, jmax), wt[(jmin - 1):-1]):
            b = int(j2pix(jj + 1)) # border pixels to exclude
            wcrop = w[b:-b, b:-b]
            var_noisy[ii, jj - (jmin - 1)] = np.var(wcrop, axis=None, ddof=1)
            skew_noisy[ii, jj - (jmin - 1)] = skew(wcrop, axis=None)
            kurt_noisy[ii, jj - (jmin - 1)] = kurt(wcrop, axis=None)
            mu3_noisy[ii, jj - (jmin - 1)] = moment(wcrop, 3, axis=None)
            mu4_noisy[ii, jj - (jmin - 1)] = moment(wcrop, 4, axis=None)
        sys.stdout.flush()

    # Save the result
    outpath = datapath("output")
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    np.save(os.path.join(outpath, "var_clean"), var_clean)
    np.save(os.path.join(outpath, "skew_clean"), skew_clean)
    np.save(os.path.join(outpath, "kurt_clean"), kurt_clean)
    np.save(os.path.join(outpath, "mu3_clean"), mu3_clean)
    np.save(os.path.join(outpath, "mu4_clean"), mu4_clean)
    np.save(os.path.join(outpath, "var_noisy"), var_noisy)
    np.save(os.path.join(outpath, "skew_noisy"), skew_noisy)
    np.save(os.path.join(outpath, "kurt_noisy"), kurt_noisy)
    np.save(os.path.join(outpath, "mu3_noisy"), mu3_noisy)
    np.save(os.path.join(outpath, "mu4_noisy"), mu4_noisy)

    print("\n")
    print_time(time.time() - start_time)


def power_spectrum(noisy=False):
    pass


def compare_peaks(noisy=False):
    """Compare my peak counts results to lenstools."""
    # Load a random LOS (out of 256) for zs = 2.
    lcdm = fetch_kappa('LCDM', map_id=np.random.randint(0, 256))
    cmap = ConvergenceMap(data=lcdm[:256, :256], angle=(MAP_SIZE / 8))

    # Compute wavelet transform with 5 scales
    wt = starlet2d(lcdm, nscales=5)
    cmap4 = ConvergenceMap(data=wt[3, :256, :256], angle=(MAP_SIZE / 8))

    # Find peaks with lenstools
    bin_edges = np.linspace(cmap4.data.min(), cmap4.data.max() * 1.01, 33)
    height, position = cmap4.locatePeaks(thresholds=bin_edges)
    xpeak, ypeak = position.value.T

    # Find peaks with wltools
    x, y, h = find_peaks2d(cmap4.data, include_border=True)

    # wlt.plotting.plot2d(lcdm, title="Austin's kappa map")
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 6))
    cmap4.visualize(cmap='gist_stern', colorbar=True, fig=fig, ax=ax0)
    ax0.set_title("lenstools : {} peaks".format(len(height)))
    ax0.scatter(xpeak, ypeak, s=25, c='c')
    plot2d(cmap4.data, fig=fig, ax=ax1)
    ax1.set_title("wltools : {} peaks".format(len(h)))
    ax1.scatter(y, x, s=25, c='c')
    plt.show()
