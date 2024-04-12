from astropy.io import fits

outfile = 'header.txt'

hdulist = fits.open('spec1d.m46.034.slit.fits')
hdr01   = hdulist[1].header
hdr01.tofile(outfile, sep='\n', endcard=False,
             padding=False, overwrite=True)
