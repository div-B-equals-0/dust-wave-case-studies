from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt

def set_wcs_limits(ax, corner_coords, wcs):
    xp, yp = skycoord_to_pixel(corner_coords, wcs)
    ax.set_xlim(*xp)
    ax.set_ylim(*yp)

if __name__ == "__main__":

   datadir = "data/Sigma-Ori/"
   hdu = fits.open(datadir + "sigma-ori-Tcol-24-70-multibin.fits")[0]

   c = SkyCoord([
       # Bottom left corner
       [84.8625309*u.deg, -2.8242803*u.deg],
       # Top right corner
       [84.5378848*u.deg, -2.2873997*u.deg],
   ])
   fig, ax = plt.subplots()
   set_wcs_limits(ax, c, WCS(hdu.header))
   print(ax.get_xlim(), ax.get_ylim())
