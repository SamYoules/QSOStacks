# Date: 9 Feb 2018
# Author: Sam Youles
# Package: picca
# stack_flux.py
# Select all spectra for the same quasar from shortqso12.txt (shortqso14.txt)
#   & create a static plot of changes in the forest over time. For comparison,
#   the stacked, weighted flux and continuum is also plotted. Creates stacks for whole data set.
# User input example: run static_stack_1.py DR14
# Spectra obtained from https://dr14.sdss.org/optical/spectrum/search
#        and https://dr12.sdss.org/optical/spectrum/search
#        and https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/test/bautista/v5_10_7/spectra/lite/

from astropy.io import fits
from picca import prep_del, io
from picca.data import forest, delta
from func import read_from_spec, read_data, stack_flux, init
import fitsio
import time
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

# Set up the figure
P.rcParams.update({'font.size':16})
fig = P.figure(figsize=(10,8), edgecolor='none', facecolor='white', dpi=80)

# Create arrays for THID etc from input files
DR = sys.argv[1]
if DR == 'DR12':
   fname = 'shortqso12.txt'
elif DR ==  'DR14':
   fname = 'shortqso14.txt'
elif DR ==  'DR15':
   fname = 'shortqso15.txt'
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

    # Fit the continua
    for d in data:
        d.cont_fit()

    # Define g in order to retrieve thid, etc for plot
    g = data[0]

    # Stack the weighted flux
    ll1, st1, wst1  = stack_flux(data, 0)
    # Fit the continuum for the stack
    d = forest(ll1, st1, wst1, g.thid, g.ra, g.dec, g.zqso, g.plate, g.mjd, g.fid, g.order)
    d.cont_fit()
    # Create continuum in same bin widths as before
    nstack = int((forest.lmax - forest.lmin) / forest.dll) + 1
    co1 = sp.zeros(nstack)
    bins=((d.ll - d.lmin) / d.dll + 0.5).astype(int)
    c = sp.bincount(bins, weights = d.co)
    co1[:len(c)] += c

    # Set limits for axes
    llmin = []
    llmax = []
    flmin = []
    flmax = []
    j = 0
    while j < len(data):
        f = data[j]
        llmin.append(min(f.ll))
        llmax.append(max(f.ll))
        flmin.append(min(f.fl))
        flmax.append(max(f.fl))
        j += 1
    xmin = 10**min(llmin)
    xmax = 10**max(llmax)
    ymin = min(flmin)
    ymax = max(flmax)
    ax = P.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
    P.xlabel(r'$\lambda \ [\AAngstr\"{o}m]$')
    P.ylabel(r'$f_\lambda \ [10^{-19} \ W m^{-2} \ nm^{-1}]$')
    for f in data:
        P.plot(10**f.ll, f.fl, lw=1, color='silver')
    P.plot(10**ll1, st1, lw=2, color='crimson')
    P.plot(10**ll1, co1, lw=2, color='crimson')

    # Add title:
    P.title(r'Stacked Flux for {} QSO: {}'.format(DR, g.thid), color='black')
    P.legend()
    P.savefig('{}static/stack_{}.png'.format(DR, g.thid))
    P.clf()

P.close()    
    

