# Date: 31 Jan 2018
# Author: Sam Youles
# Package: picca
# Stack_deltas_clip.py
# Select all spectra from shortqso12.txt [14] (contains DR12 [14] quasars from plates
# 7338, 7339 and 7340) & create fits files of deltas from a stack of the
#   forests for each quasar.
#   Use sigma-clipping to exclude outliers.
#   i.e. fit the individual continua, stack weighted flux, fit the continuum for
#   the stack, and then make the deltas.
#   Also create a stack of the deltas for each of the forests.
# User input example: run stack_deltas.py DR14
# Spectra obtained from https://dr12.sdss.org/optical/spectrum/search [14]
#        and https://dr12.sdss.org/optical/spectrum/search
#        and https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/test/bautista/v5_10_7/spectra/lite/

from astropy.io import fits
from picca import prep_del, io
from picca.data import forest, delta
from func import read_from_spec, read_data, stack_flux, init
import fitsio
import numpy as N
import scipy as sp
from scipy.interpolate import interp1d
import sys

print sys.argv[1]

# Initialise
init()
order = 1
mode = 'spec'
npix_min = 50

# Create arrays for THID etc from input files
DR = sys.argv[1]
fname = 'shortqso{}.txt'.format(DR[2:])
THID, PLATE, MJD, FIBERID, RA, DEC, Z = N.loadtxt(fname, unpack = 1)
THID = THID.astype(int)
PLATE = PLATE.astype(int)
MJD = MJD.astype(int)
FIBERID = FIBERID.astype(int)

# Create a shorter array containing only unique Thing_ids
UTHID = N.unique(THID)

# Create shorter arrays for other variables, corresponding to unique THID
for T in UTHID:
    w = (THID == T)
    thid = THID[w]
    plate = PLATE[w]
    mjd = MJD[w]
    fiberid = FIBERID[w]
    ra = RA[w]
    dec = DEC[w]
    z = Z[w]

    # Read data from relevant spec files
    data = read_data(DR, thid, ra, dec, z, plate, mjd, fiberid, order, mode)

    # Find the mean flux value of the forest for each epoch
    epochflux = []
    g = data[0]
    for f in data:
        epochflux.append(N.sum(f.fl * f.iv)/N.sum(f.iv))

    # Find median and 3-sigma value for epochflux
    flux = N.asarray(epochflux)
    median = N.median(flux)
    sixteen = N.percentile(flux, 16)
    eightyfour = N.percentile(flux, 84)
    sigma = 0.5 * (eightyfour - sixteen)
    sigma3minus = median - sigma * 3.
    sigma3plus = median + sigma * 3.

    # Discard any epochs that are 3 sigma or more away from the mean
    i = 0
    l = []
    for f in data:
        if (flux[i] < sigma3plus and flux[i] > sigma3minus):
            l.append(f)
        i += 1
    data = N.asarray(l)

    # Fit continuum for all spectra
    for d in data:
        d.cont_fit()

    # Define g in order to retrieve thid, etc for plot
    try:
        g = data[0]
    except:
        continue

    # Create a stack of all the deltas
    ll2, de2, wst2  = stack_flux(data, 1)

    # Create deltas of the stack (of weighted flux)
    # Stack the weighted flux
    ll1, st1, wst1  = stack_flux(data, 0)
    # Fit the continuum for the stack
    d = forest(ll1, st1, wst1, g.thid, g.ra, g.dec, g.zqso, g.plate, g.mjd, g.fid, g.order)
    d.cont_fit()
    # Create deltas in same bin widths as before
    nstack = int((forest.lmax - forest.lmin) / forest.dll) + 1
    de1 = sp.zeros(nstack)
    bins=((d.ll - d.lmin) / d.dll + 0.5).astype(int)
    c = sp.bincount(bins, weights = d.fl / d.co - 1)
    de1[:len(c)] += c


    # Find start and end of forest
    lambdamin = []
    lambdamax = []

    for f in data:
        lambdamin.append(min(f.ll))
        lambdamax.append(max(f.ll))

    llmin = min(lambdamin)
    llmax = max(lambdamax)

    j = 0
    k = 0

    while (ll1[j] < llmin):
        j += 1

    k = len(ll1) - 1
    while (ll1[k] > llmax):
        k -= 1
 
    # Write to fits files
    out = fitsio.FITS("{}stacks_clip/delta-{}".format(DR, d.thid)+".fits.gz",'rw',clobber=True)
    hd={}
    hd["RA"]=d.ra
    hd["DEC"]=d.dec
    hd["Z"]=d.zqso
    hd["PMF"]="{}-{}-{}".format(d.plate,d.mjd,d.fid)
    hd["THING_ID"]=d.thid
    hd["PLATE"]=d.plate
    hd["MJD"]=d.mjd
    hd["FIBERID"]=d.fid
    hd["ORDER"]=d.order

    cols=[ll1[j:k], de1[j:k], wst1[j:k], ll2[j:k], de2[j:k], wst2[j:k]]
    names=['LOGLAM1','DELTA1','WEIGHT1','LOGLAM2','DELTA2','WEIGHT2']
    out.write(cols,names=names,header=hd)
    out.close()

