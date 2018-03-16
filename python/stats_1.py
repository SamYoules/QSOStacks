# Date: 9 Feb 2018
# Author: Sam Youles
# Package: picca
# stats_1.py
# Select all spectra for the same quasar from shortqso12.txt (shortqso14.txt)
#   & pull out values for flux at a specific wavelength for all spectra and
#   create histogram of flux, with mean, median, and sigma values.
# User input example: run stats_1.py DR14 496389662
# Use the file matching_thids.txt to find THING_IDs for quasars in DR12 & DR14/DR15.
# If you need to see plate-mjd-fid, use commonqsos.txt instead.
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
#T = input("Please enter THING_ID for {}: ".format(DR))
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


# Fit the continua
newdata = []
for d in data:
    d.cont_fit()
    if d.bad_cont is not None:
        print 'bad cont: ' + str(T)
        continue
    else:
        newdata.append(d)
data = N.asarray(newdata)

# Define g in order to retrieve thid, etc for stats, continuum fitting and plot
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

# Get mean flux for each epoch
epochflux = []
for f in data:
    epochflux.append(N.sum(f.fl * f.iv)/N.sum(f.iv))

# Find mean, standard deviations, median
flux = N.asarray(epochflux)
meanflux = N.mean(epochflux)
#n = flux.size
#sq_mean = (N.sum(flux**2)) / n
#mean_sq = (N.sum(flux) / n)**2
#sigma = N.sqrt(abs(sq_mean - mean_sq))
median = N.median(flux)
sixteen = N.percentile(flux, 16)
eightyfour = N.percentile(flux, 84)

sigma = 0.5 * (eightyfour - sixteen)
sigma_1_minus = median - sigma
sigma_1_plus = median + sigma
sigma_2_minus = median - sigma * 2.
sigma_2_plus = median + sigma * 2.
sigma_3_minus = median - sigma * 3.
sigma_3_plus = median + sigma * 3.
sigma_4_minus = median - sigma * 4.
sigma_4_plus = median + sigma * 4.
sigma_5_minus = median - sigma * 5.
sigma_5_plus = median + sigma * 5.

# Set bin widths and axes
xmin = min(flux)
xmax = max(flux)
binwdth = (xmax - xmin)/150
P.axes(xlim=(xmin - binwdth,xmax + binwdth))

mybins = N.arange(xmin, xmax, binwdth)
P.hist(flux, bins = mybins, rwidth = 0.9, align = 'mid')

# Shade up to 5 sigma from median 
P.axvspan(sigma_1_minus, median, alpha=0.2, color='blue')
P.axvspan(median, sigma_1_plus, alpha=0.2, color='blue')

P.axvspan(sigma_2_minus, sigma_1_minus, alpha=0.14, color='blue')
P.axvspan(sigma_1_plus, sigma_2_plus, alpha=0.14, color='blue')

P.axvspan(sigma_3_minus, sigma_2_minus, alpha=0.12, color='blue')
P.axvspan(sigma_2_minus, sigma_3_plus, alpha=0.12, color='blue')

P.axvspan(sigma_4_minus, sigma_3_minus, alpha=0.1, color='blue')
P.axvspan(sigma_3_minus, sigma_4_plus, alpha=0.1, color='blue')

P.axvspan(sigma_5_minus, sigma_4_minus, alpha=0.08, color='blue')
P.axvspan(sigma_4_minus, sigma_5_plus, alpha=0.08, color='blue')

P.axvline(x=meanflux, color = 'red', label = 'mean', lw = 2)

P.axvline(x=median, color = 'red', label = 'median', lw = 2, ls ='--')
P.axvline(x=sixteen, color = 'yellow', label = '16', lw = 2, ls ='--')
P.axvline(x=eightyfour, color = 'green', label = '84', lw = 2, ls ='--')

P.legend()
P.xlabel(r'$Average flux per forest \ [10^{-19} \ W m^{-2} \ nm^{-1}]$')
P.title(r'{} QSO: {}'.format(DR, T), color='black')
P.show()
P.savefig('{}stats/stats_{}.png'.format(DR, T))


