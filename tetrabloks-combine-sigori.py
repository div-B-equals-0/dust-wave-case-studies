import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import regions as rg

datadir = "data/Sigma-Ori/"
nlist = [1, 2, 4, 8, 16, 32]

win = rg.read_ds9("data/Sigma-Ori/sigma-ori-window-70.reg")[0]

for wavA, wavB, minA, minB in [
        ["24", "70", 3.0, 200.0],
        ["160", "350", 30.0, 5.0]
]:
    outimA, outimB = None, None
    for n in reversed(nlist):
        fnA = datadir + f"sigma-ori-F-nu-MJy-sr-{wavA}-bin{n:03d}.fits"
        fnB = datadir + f"sigma-ori-F-nu-MJy-sr-{wavB}-bin{n:03d}.fits"

        hduA = fits.open(fnA)["SCALED"]
        hduB = fits.open(fnB)["SCALED"]
        if outimA is None:
            hdr = hduA.header
            outimA = np.empty_like(hduA.data)
            outimB = np.empty_like(hduA.data)

        # Minimum flux requirement
        mask = (n**2 * hduA.data > minA) & (n**2 * hduB.data > minB)
        outimA[mask] = hduA.data[mask]
        outimB[mask] = hduB.data[mask]


    # Cut off the image at the border of the original diamond window
    winmask = win.to_pixel(WCS(hduA.header)).to_mask()
    mask = ~winmask.to_image(hduA.data.shape).astype(bool)
    outimA[mask] = np.nan
    outimB[mask] = np.nan

    out_fnA =  datadir + f"sigma-ori-F-nu-MJy-sr-{wavA}-multibin.fits"
    out_fnB =  datadir + f"sigma-ori-F-nu-MJy-sr-{wavB}-multibin.fits"
    fits.PrimaryHDU(header=hdr, data=outimA).writeto(out_fnA, overwrite=True)
    fits.PrimaryHDU(header=hdr, data=outimB).writeto(out_fnB, overwrite=True)
