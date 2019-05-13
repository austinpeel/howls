# Constant descriptors of DUSTGRAIN-pathfinder maps
import astropy.units as u


__MAP_SIZE__ = 5 * u.deg
__MAP_AREA__ = __MAP_SIZE__**2
__NPIX__ = 2048
__PIX_SCL__ = (__MAP_SIZE__ / __NPIX__).to(u.arcmin)
__PIX_AREA__ = __PIX_SCL__**2
__ZS__ = 2.0
