from astropy.cosmology import FlatLambdaCDM

def scale_calc(z, H0, Omega_M, Omega_Lam):
    cosmo = FlatLambdaCDM(H0, Om0=Omega_M) # a selected cosmological model
    d_A = cosmo.kpc_proper_per_arcmin(z)
    d_A = d_A.value/60 # kpc/arcmin --> kpc/arcsec
    return d_A # kpc/arcsec
