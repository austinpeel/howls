import numpy as np
import matplotlib.pyplot as plt
from .io import fetch_saved_stats
from scipy.stats import gaussian_kde


def plot_peaks(jval=3, noisy=False, ax=None):
    """Plot the aperture mass peaks function at a given scale.

    Parameters
    ----------
    jval : int, optional
        Wavelet scale at which aperture mass statistics have been computed.
    noisy : bool
        Whether to plot statistics of the noisy or clean maps.

    """
    # Load values
    filename = 'peaks_' + ['clean', 'noisy'][noisy]
    all_vals = fetch_saved_stats(filename, verbose=True)
    vals = all_vals[:, jval - 3, :]
    mean = vals.mean(axis=0)
    std = vals.std(axis=0, ddof=1) / np.sqrt(256)
    bin_edges, dnu = np.linspace(0, 10.1, 102, retstep=True)
    bins = (bin_edges[1:] + bin_edges[:-1]) / 2

    # Plot
    if ax is None:
        fig, ax = plt.subplots(1, 1, facecolor='w')
    ax.errorbar(bins, mean, yerr=std, label=['clean', 'noisy'][noisy])
    ax.set_xlabel("S/N")
    ax.set_ylabel("peak count")
    ax.set_title("scale $j={}$".format(jval))
    ax.legend()


def plot_histogram(stat, jval=3, noisy=False, nbins=16, ax=None):
    """Plot the histogram of aperture mass statistics at a given scale.

    Parameters
    ----------
    stat : str
        Statistic to plot, one of {'var', 'skew', 'kurt', 'mu3', 'mu4'}
    jval : int, optional
        Wavelet scale at which aperture mass statistics have been computed.
    noisy : bool
        Whether to plot statistics of the noisy or clean maps.
    nbins : int, optional
        Number of histogram bins. Default is 12.

    """
    # Load values
    filename = stat + '_' + ['clean', 'noisy'][noisy]
    all_vals = fetch_saved_stats(filename, verbose=True)
    vals = all_vals.T[jval - 3]

    # Approximate the pdf by kernel density estimation
    kde = gaussian_kde(vals, bw_method='silverman')

    # Plot
    if ax is None:
        fig, ax = plt.subplots(1, 1, facecolor='w')
    ax.hist(vals, nbins, density=True, fc='k', alpha=0.4)
    x = np.linspace(0.9 * vals.min(), 1.1 * vals.max(), 10**4)
    ax.plot(x, kde.pdf(x), c='k', linewidth=1.2, label=filename)
    ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax.set_xlabel(stat)
    ax.set_ylabel("PDF")
    ax.set_title("wavelet scale $j={}$".format(jval))
    ax.legend()
