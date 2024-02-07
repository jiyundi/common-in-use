from astropy.io import fits
hdulist=fits.open('spec1d.648.001.90055.fits')
print('Your .fits file has been opened...')

#### ============ TUTORIAL ============ ####
##1. Each .fits contains many headers (HDUs)



##2. Every .fits will firstly give primary header 
#    ("Table of Contents") and then data extensions
#    Note: each extension embedded bintable and keywords.

# FITS FILE
#   '=====> hdul[0]: "Table of Contents"
#       '------->>>>> primary HDU: 
#                     "No./Name/Ver/Type///////Cards/Dimensions//Format "
#                     "0   IM1  1   PrimaryHDU 152   ()                 "
#                     "1   FLUX 1   ImageHDU   987   (5632, 190) float64" 
#   '=====> hdul[1]: Wavelength array vs. slit ID with keywords
#       '------->>>>> wavelength array vs. slit ID
#       '------->>>>> keywords: 
#                     "NAXIS = 2 / Number of data axes"

# Primary header (Table of Contents)
hdr00info = hdulist.info()
print(hdr00info)



##3. For those unfamiliar with FITS headers, 
#    they consist of a list of 80 byte “cards”,
#    where a card contains a keyword, a value, 
#    and a comment. I.e.
#      "Cards" = "how many keywords there", 
#    e.g.
#      MASKNAME= 'Egami2019a_1_476'   / Mask name from the database

# Primary header keywords (Table of Contents)
# hdr00kwds = hdulist[0].header
# print(hdr00kwds)

# 1st-item header keywords (data's)
hdr01kwds = hdulist[1].header
print(hdr01kwds)

# 1st-item header data (actual data, in (1D arrays x slits).)
hdr01datacompare = hdulist[1].data

# w = hdulist[1].data["LAMBDA"][0]
# f = hdulist[1].data["FLUX"][0]
# e = hdulist[1].data["IVAR"][0]