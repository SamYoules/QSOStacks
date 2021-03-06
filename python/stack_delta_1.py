# Date: 31 Jan 2018
# Author: Sam Youles
# Package: picca
# stack_delta_1.py
# Select all spectra for the same quasar from shortqso12.txt (shortqso14, shortqso15) & create fits files of deltas from a stack of the forests.
#   i.e. fit the individual continua, stack weighted flux, fit the continuum for
#   the stack, and then make the deltas.
#   Also create a stack of the deltas for each of the forests.
# User input example: run stack_delta_1.py DR14 496389662
# Spectra obtained from https://dr12.sdss.org/optical/spectrum/search
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
T = sys.argv[2]
w = (THID == T)
thid = THID[w]
plate = PLATE[w]
mjd = MJD[w]
fiberid = FIBERID[w]
ra = RA[w]
dec = DEC[w]
z = Z[w]

order = 1
mode = 'spec'

# Read data from relevant spec files
data = read_data(DR, thid, ra, dec, z, plate, mjd, fiberid, order, mode)

# Fit continuum for all spectra
for d in data:
    d.cont_fit()
# Define g in order to retrieve thid, etc for plot
g = data[0]

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
j = 0
while j < len(data):
    f = data[j]
    lambdamin.append(min(f.ll))
    lambdamax.append(max(f.ll))
    j += 1
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
out = fitsio.FITS("{}stacks/delta-{}".format(DR, d.thid)+".fits.gz",'rw',clobber=True)
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



