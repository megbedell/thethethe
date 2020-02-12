import numpy as np
from astropy.io import fits
import dace
from dace.spectroscopy import Spectroscopy
import tarfile
import os

speed_of_light = 3.e8 # m/s

def download_spectra(starname, n_files=2):
    # query:
    observedTargets = Spectroscopy.query_database(limit=100,
                                              filters={'public': {'is':True},
                                                       'obj_id_catname': {'contains': starname},
                                                       'ins_name': {'contains': 'HARPS'},
                                                       'ins_mode': {'contains': 'HARPS'}})
    to_download = []
    j = 0
    for i,n in enumerate(observedTargets['obj_id_catname']):
        if n == starname:
            to_download.append(observedTargets['file_rootpath'][i].strip('.fits'))
            j += 1
        if j >= n_files:
            break
    if len(to_download) < n_files:
        print("Insufficient data returned for {0}".format(starname))
        return None
    # download:
    Spectroscopy.download_files('s1d', to_download, '{0}.tar.gz'.format(starname))
    # unpack:
    tar = tarfile.open('{0}.tar.gz'.format(starname), "r:gz")
    files = []
    for tarinfo in tar:
        files.append('harps-data/'+tarinfo.name)
    tar.extractall(path='./harps-data/')
    tar.close()
    os.remove('{0}.tar.gz'.format(starname))
    return files

def read_spectrum(specfile):
    hdul = fits.open(specfile)
    flux = hdul[0].data
    wave = np.arange(len(flux)) * hdul[0].header['CDELT1'] + hdul[0].header['CRVAL1']
    return wave, flux

def calc_rv_err(wave, flux, dw=0.01, perpix=False):
    flux = np.abs(flux) # TOTAL HACK relying on negative fluxes being v v small anyway
    err_flux = np.sqrt(flux) # Poisson noise
    df_dw = (np.roll(flux, -1) - np.roll(flux, 1)) / (2. * dw) # local slope - will fail for first/last pixels
    df_dv = df_dw * wave / speed_of_light
    err_rv_perpix = err_flux / df_dv
    mask = (wave > 5300.) & (wave < 5340.)
    err_rv_perpix = err_rv_perpix[~mask] # trim between chips
    err_rv_perpix = err_rv_perpix[1:-2] # trim ends
    if perpix:
        flux_perpix = flux[~mask]
        flux_perpix = flux_perpix[1:-2]
        return err_rv_perpix, flux_perpix
    err_rv = 1./np.sqrt(np.sum(1./err_rv_perpix**2))
    return err_rv

