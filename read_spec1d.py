import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits

def binning_1d(arr, binfac):
    """ Bins up two or three column spectral data by a specified factor. """
    binfac = int(binfac)
    nbins  = len(arr) // binfac
    binned = np.zeros((nbins, 1))
    for i in range(len(binned)):
        spec_slice = arr[i*binfac:(i+1)*binfac]
        binned[i] = np.mean(spec_slice)
    return binned[:,0]

def read_spec1d(spec1dfilepath, binfactor):
    # **observed** spectrum from Keck/DEIMOS
    # Reduced by PypeIt and pypeittospecpro.py
    hdulist = fits.open(spec1dfilepath)
    print('Opened...         ' + spec1dfilepath)
    hdr01data = hdulist[1].data
    w = hdr01data["LAMBDA"][0]
    f = hdr01data["FLUX"][0]
    e = hdr01data["IVAR"][0]
    arr_obs_wave = binning_1d(w, binfactor)
    arr_obs_flux = binning_1d(f, binfactor)
    arr_obs_erro = binning_1d(e, binfactor)
    
    return arr_obs_wave, arr_obs_flux, arr_obs_erro

binfactor      = 10
spec1dfilepath = 'spec1d.m46.034.slit.fits'
arr_obs_wave, arr_obs_flux, arr_obs_erro = read_spec1d(spec1dfilepath, binfactor)
arr_obs = np.concatenate(([arr_obs_wave],
                          [arr_obs_flux],
                          [arr_obs_erro]), axis=0)
