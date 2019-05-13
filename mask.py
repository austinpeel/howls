import numpy as np
from dustgrain.constants import NPIX
from wltools.utils import radial_distance_matrix


def make_mask(nholes=10, size=50):
    mask = np.ones((NPIX, NPIX)).astype(bool)
    for ii in range(nholes):
        x = np.random.randint(0, NPIX)
        y = np.random.randint(0, NPIX)
        holy = radial_distance_matrix(NPIX, center=(x, y)) > size
        mask = mask & holy
    return mask.astype(int8)
