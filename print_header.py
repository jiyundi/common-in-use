from astropy.io import fits

hdulist = fits.open('spec1d.m46.034.slit.fits')
hdr01   = hdulist[1].header
print(hdr01)
