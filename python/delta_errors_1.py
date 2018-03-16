# Date: 31 Jan 2018
# Author: Sam Youles
# Package: picca
# delta_errors_1.py
# Select all spectra for the same quasar from shortqso12.txt (or 14 or 15) &
#   plot errors for deltas-of-stack and stacks-of-deltas.
#   i.e. fit the individual continua, stack weighted flux, fit the continuum for
#   the stack, and then make the deltas.
#   Also create a stack of the deltas for each of the forests.
#   Plot errors against wavelength for both stacks.
# User input example: run delta_errors_1.py DR14 496389662
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

# Create shorter arrays for other variables, corresponding to specific THID
T = int(sys.argv[2])
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

# Set up the figure
P.rcParams.update({'font.size':16})
fig = P.figure(figsize=(10,8), edgecolor='none', facecolor='white', dpi=80)

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
    d = forest(ll1, st1, wst1, d.thid, d.ra, d.dec, d.zqso, d.plate, d.mjd, d.fid, d.order)
    d.cont_fit()

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
we1 = N.mean(wei1, axis=0)
we2 = N.mean(wei2, axis=0)

P.plot(10**loglam1, 1./N.sqrt(we1), lw = 2, color = 'royalblue', label=('Deltas of Stack'))
P.plot(10**loglam1, 1./N.sqrt(we2), lw = 2, color = 'darkred', label=('Stack of Deltas'))
P.xlabel('wavelength [Angstrom]')
P.ylabel('Errors')
P.legend()
P.title(r'Errors for {} QSO: {}'.format(DR, T), color='black')
P.show()
P.savefig('{}JK/errors_{}.png'.format(DR, T))
P.close()



