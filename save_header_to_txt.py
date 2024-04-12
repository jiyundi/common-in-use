from astropy.io import fits

outfile = open('header_.txt', 'w')

hdulist = fits.open('spec1d.m46.034.slit.fits')
hdr01   = hdulist[1].header
outfile.write(hdr01)
outfile.close()

infile = open('header_.txt',  'r')
oufile = open('header.txt', 'w')
for line in infile:
    for j in range(len(line)//80 - 1):
        oufile.write(line[j*80:(j+1)*80]+'\n')
infile.close()
oufile.close()
