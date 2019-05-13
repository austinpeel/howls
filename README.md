# howls
Convergence map analysis for HOWLS (Higher-Order Weak-Lensing Statistics) project

## Step 2 (in progress)

## Step 1 (16 Nov. 2018)

Peaks have been computed as a function of wavelet scale and signal-to-noise (SNR).
A peak is defined as a pixel in a filtered convergence map (equivalently aperture mass)
whose value is larger than its eight neighbors.

Scales
------
I considered 5 scales (j = [3, 4, 5, 6, 7] in the wavelet formalism), which correspond
to aperture mass filters of angular size theta = [1.17, 2.34, 4.69, 9.38, 18.8] arcmin.

Signal-to-noise
---------------
SNR bins are defined as follows: 0.0 to 10.0 with a spacing of 0.1. There is also a final
highest bin collecting peaks with SNR >= 10.0. This gives a total of 101 bins.

The bin edges are treated so that the lower limit is inclusive, while the right limit is
exclusive. For example, a peak with SNR = 3.97 falls in the [3.9, 4.0) bin, while a peak
with SNR = 4.0 falls in the [4.0, 4.1) bin.

In the clean data case, the noise at a given scale is computed as the standard deviation
of the wavelet map at that scale. In the noisy data case, the noise is computed as the
standard deviation, at the given scale, of a gaussian random field with dispersion
corresponding to the true kappa pixel noise added.

Notes
-----
Pixels in filtered maps closer to the border than the filter size theta were excluded
from the peak finding process.

The wide SNR range was chosen to account for the large number of high SNR peaks detected
at smaller scales (mostly in the clean case). The fine spacing of 0.1 was chosen with the
idea that you can always coarse-grain the datavector to get fewer bins without losing
information, but the opposite isn't true.

Data access
-----------
Peak counts for both the clean and noisy cases have been stored as 3D numpy arrays with
shape (256, 5, 101) corresponding to (map number, scale, SNR bin).
