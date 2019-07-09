import sys
import os
import numpy as np
from astropy.io import fits
sys.path.append("../multibin-maps")
from rebin_utils import downsample, oversample, pad_array

nlist = [1, 2, 4, 8, 16, 32]
mingoods = [2]*len(nlist)

try:
    infile = sys.argv[1]
except IndexError:
    sys.exit(f"Usage: {sys.argv[0]} FITSFILE")

hdu = fits.open(infile)[0]
if hdu.data is None:
    hdu = fits.open(infile)[1]
hdr = hdu.header

# Maximum binning
nmax = nlist[-1]
# Pad arrays to nearest multiple of nmax
im = pad_array(hdu.data, nmax)

# First version does not use a separate weight array, and uses NaN as mask
w = np.ones_like(im)
# Maybe do a star mask later
m = np.isfinite(im)

for n, mingood in zip(nlist, mingoods):
    im[~m] = 0.0
    outfile = infile.replace(".fits", f"-bin{n:03d}.fits")
    print("Saving", outfile)
    # Save both the scaled image and the weights after scaling back up to full res
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(data=oversample(im, n), header=hdr, name="scaled"),
        fits.ImageHDU(data=oversample(w, n), header=hdr, name="weight"),
    ]).writeto(outfile, overwrite=True)
    # Now do the rebinning by a factor of two
    [im,], m, w = downsample([im,], m, weights=w, mingood=mingood)
