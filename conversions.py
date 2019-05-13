from .constants import __PIX_SCL__
from wltools.conversions import to_arcmin


def arcmin2pix(theta):
    """Compute pixels based on 25 deg^2 HOWLS kappa maps of 2048^2 pixels"""
    return (to_arcmin(theta) / __PIX_SCL__).value


def pix2arcmin(n):
    """Compute angular size in arcmin corresponding to n pixels."""
    return n * __PIX_SCL__


def j2arcmin(j):
    """Compute aperture in arcmin corresponding to wavelet scale j=1,2,..."""
    if j < 1:
        print("Invalid j: must be 1, 2, ...")
        return 0
    return 2**int(j) * __PIX_SCL__


def j2pix(j):
    """Compute pixel size corresponding to wavelet scale j."""
    if j < 1:
        print("Invalid j: must be 1, 2, ...")
        return 0
    return float(2**int(j))
