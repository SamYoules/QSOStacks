# Date: 30 Jan 2018
# Author: Sam Youles
# Package: picca
# Name: plot_delta_1.py
# Create plots from stack files, comparing stacks of deltas with deltas of stacks 
# User input example: run plot_delta_1.py DR14 496389662
# Stack files (DR14stacks/delta-nnnnnnnn.fits.gz) contain:
# 	header: RA, DEC,Z,PMF, THING_ID, PLATE, MJD, FIBERID, ORDER
# 	DATA: LOGLAM1, DELTA1, WEIGHT1, LOGLAM2, DELTA2, WEIGHT2

from astropy.io import fits
import glob
import numpy as N
import pylab as P

# Set up the figure
P.rcParams.update({'font.size':16})
fig = P.figure(figsize=(10,8), edgecolor='none', facecolor='white', dpi=80)

# Open each fits file and make plot
DR = sys.argv[1]
T = sys.argv[2]
a = fits.open('{}stacks/delta-{}.fits.gz'.format(DR, T))
d = a[1].data

P.plot(10**d.LOGLAM1, d.DELTA1, lw = 2, color = 'royalblue', label=('Deltas of Stack'))
P.plot(10**d.LOGLAM2, d.DELTA2, lw = 2, color = 'darkred', label=('Stack of Deltas'))
P.xlabel('wavelength [Angstrom]')
P.ylabel('Deltas')
P.legend()
P.title(r'Ly-$\alpha$ {} Deltas for QSO: {}'.format(DR, T), color='black')
P.savefig('{}deltas/deltas_{}.png'.format(DR, T))
P.close()

