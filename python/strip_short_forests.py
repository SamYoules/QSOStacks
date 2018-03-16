# Date: 8 Jan 2018
# Author: Sam Youles
# Package: picca
# Deselect all spectra with short forests
# Spectra obtained from https://dr14.sdss.org/optical/spectrum/search

from astropy.io import fits
from picca import prep_del, io
from picca.data import forest, delta
import fitsio
import time
import numpy as N
import pylab as P
import scipy as sp
from scipy.interpolate import interp1d
import sys

def read_from_spec(in_dir,thid,ra,dec,zqso,plate,mjd,fid,order,mode):
    pix_data = []
    for t,r,d,z,p,m,f in zip(thid,ra,dec,zqso,plate,mjd,fid):
        try:
            fin = in_dir + '/spec-%d-%d-%04d.fits'%(p, m, f)
            h = fitsio.FITS(fin)
        except IOError:
            print("error reading {}\n".format(fin))
            continue

        print("{} read\n".format(fin))
        ll = h[1]["loglam"][:]
        fl = h[1]["flux"][:]
        iv = h[1]["ivar"][:]*(h[1]["and_mask"][:]==0)
        d = forest(ll,fl,iv, t, r, d, z, p, m, f,order)
        pix_data.append(d)
        h.close()
    return pix_data

def read_data(in_dir, thid, ra, dec, z, plate, mjd, fiberid, order, mode):

    # Read forest pixel arrays from each spec file into pix_data list
    pix_data = read_from_spec(in_dir,thid, ra, dec, z, plate, mjd, fiberid, order, mode=mode)

    if not pix_data is None:
        print("{} read from pix\n".format(len(pix_data)))
    if not pix_data is None and len(pix_data)>0:
        data = pix_data

    # strip out short forests and non-numeric values from data
    for d in data:
        if not hasattr(d,'ll') or len(d.ll) < npix_min:
            print "{} forest too short\n".format(d)
            continue

        if sp.isnan((d.fl*d.iv).sum()):
            print "{} nan found\n".format(d)
            continue

        k = (str(d.thid) + ' ' + str(d.plate) + ' ' + str(d.mjd) + ' ' + str(d.fid) + ' ' + str(d.ra) + ' ' + str(d.dec) + ' ' + str(d.zqso))
        print >>Fout, k

# Initialise
forest.lmin = sp.log10(3600.)
forest.lmax = sp.log10(5500.)
forest.lmin_rest = sp.log10(1040.)
forest.lmax_rest = sp.log10(1200.)
forest.rebin = 3
forest.dll = 3*1e-4
forest.dla_mask = 0.8

a = fits.open('iter.fits.gz')

iter_data = a[2].data
forest.var_lss = interp1d(iter_data.loglam, iter_data.var_lss, fill_value = 'extrapolate', kind = 'nearest')
forest.eta = interp1d(iter_data.loglam, iter_data.eta, fill_value = 'extrapolate', kind = 'nearest')
forest.fudge = interp1d(iter_data.loglam, iter_data.fudge, fill_value = 'extrapolate', kind = 'nearest')
iter_data = a[3].data
forest.mean_cont = interp1d(iter_data.loglam_rest, iter_data.mean_cont, fill_value = 'extrapolate', kind = 'nearest')

npix_min = 50

DR = sys.argv[1]

# Create arrays for THID etc from input file
THID, PLATE, MJD, FIBERID, RA, DEC, Z = N.loadtxt('{}qsoinfo.txt'.format(DR), unpack = 1)

THID = THID.astype(int)
PLATE = PLATE.astype(int)
MJD = MJD.astype(int)
FIBERID = FIBERID.astype(int)
order = 1
mode = 'spec'

Fout = open("shortqso{}.txt".format(DR[2:]), "w")

# Read data from relevant spec files
read_data(DR, THID, RA, DEC, Z, PLATE, MJD, FIBERID, order, mode)
  
Fout.close()

