# Date: 7 Mar 2018
# Author: Sam Youles
# Package: picca
# delta_errors.py
# Select all spectra for each quasar from shortqso12.txt (0r 14 or 15) & plot
#   errors for deltas-of-stack and stacks-of-deltas.
#   i.e. fit the individual continua, stack weighted flux, fit the continuum for
#   the stack, and then make the deltas.
#   Also create a stack of the deltas for each of the forests.
#   Plot errors against wavelength for both stacks.
#   Used sigma-clipped data sample.
# User input example: run delta_errors_clip.py DR14
# Spectra obtained from https://dr12.sdss.org/optical/spectrum/search
#        and https://dr12.sdss.org/optical/spectrum/search
#        and https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/test/bautista/v5_10_7/spectra/lite/

from astropy.io import fits
from picca import prep_del, io
from picca.data import forest, delta
from func import read_from_spec, read_data, stack_flux, init
import fitsio
import numpy as N
import pylab as P
import scipy as sp
from scipy.interpolate import interp1d
import sys


# Initialise
init()
order = 1
mode = 'spec'

# Set up the figure
P.rcParams.update({'font.size':16})
fig = P.figure(figsize=(10,8), edgecolor='none', facecolor='white', dpi=80)

# Create arrays for THID etc from input files
DR = sys.argv[1]
fname = 'short{}.txt'.format(DR[2:])
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

    # Find start and end of forest
    lambdamin = []
    lambdamax = []
    j = 0
    while j < len(data):
        f = data[j]
        lambdamin.append(min(f.ll))
        lambdamax.append(max(f.ll))
        j += 1
    llmin = min(lambdamin)
    llmax = max(lambdamax)

    # Create 2D arrays for the weights on each pixel in each resample
    wei1 = []
    wei2 = []

    # Create jackknife resamples from data
    for i in range(len(data)):
        l = []
        for j in range(len(data)):
            if j != i:
                l.append(data[j])
        resample_data = N.asarray(l)

        # Fit continuum for all spectra
        for d in resample_data:
            d.cont_fit()

        # Create a stack of all the deltas
        ll2, de2, wst2  = stack_flux(resample_data, 1)

        ## Create deltas of the stack (of weighted flux)
        # Stack the weighted flux
        ll1, st1, wst1  = stack_flux(resample_data, 0)

        # Fit the continuum for the stack
        g = data[0]
        d = forest(ll1, st1, wst1, g.thid, g.ra, g.dec, g.zqso, g.plate, g.mjd, g.fid, g.order)
        try:
            d.cont_fit()
        except:
            print 'Error fitting continuum: ' + str(g.thid)
            break
        # Create deltas and weights in same bin widths as before
        nstack = int((forest.lmax - forest.lmin) / forest.dll) + 1
        de1 = sp.zeros(nstack)
        bins=((d.ll - d.lmin) / d.dll + 0.5).astype(int)
        c = sp.bincount(bins, weights = d.fl / d.co - 1)
        de1[:len(c)] += c
        eta = forest.eta(d.ll)
        var_lss = forest.var_lss(d.ll)
        iv = d.iv / eta
        we = iv * d.co**2 / (iv * d.co**2 * var_lss + 1)
        c = sp.bincount(bins, weights = we)
        wst1[:len(c)] += c

        # Get rid of leading and trailing zeros
        w = (wst1 != 0.)
        we1 = wst1[w]
        wei1.append(we1)
        w = (wst2 != 0.)
        we2 = wst2[w]
        wei2.append(we2)

    # Get rid of leading and trailing zeros
    w = (ll1 >= llmin)
    l1 = ll1[w]
    w = (l1 <= llmax)
    loglam1 = l1[w]

    # Get the mean of the weights for each pixel in the forest (each column in the 2D arrays)
    try:
        we1 = N.mean(wei1, axis=0)
        we2 = N.mean(wei2, axis=0)
    except:
        print T, ' arrays of different size'
        continue

    P.plot(10**loglam1, 1./N.sqrt(we1), lw = 2, color = 'royalblue', label=('Deltas of Stack'))
    P.plot(10**loglam1, 1./N.sqrt(we2), lw = 2, color = 'darkred', label=('Stack of Deltas'))
    P.xlabel('wavelength [Angstrom]')
    P.ylabel('Errors')
    P.legend()
    P.title(r'Errors for {} QSO: {}'.format(DR, T), color='black')
    P.savefig('{}JKclip/errors_{}.png'.format(DR, T))
    P.clf()
P.close()


