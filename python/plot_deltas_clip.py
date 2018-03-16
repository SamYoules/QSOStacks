# Date: 30 Jan 2018
# Author: Sam Youles
# Package: picca
# Name: plot_deltas_clip.py
# Create plots from sigma-clipped stack files, comparing stacks of deltas with deltas of stacks 
# User input example: run plot_deltas.py DR14
# Stack files (DR12clippedstacks/delta-nnnnnnnn.fits.gz) contain:
# 	header: RA, DEC,Z,PMF, THING_ID, PLATE, MJD, FIBERID, ORDER
# 	DATA: LOGLAM1, DELTA1, WEIGHT1, LOGLAM2, DELTA2, WEIGHT2
	
from astropy.io import fits
import glob
import numpy as N
import pylab as P
import sys

# Set up the figure
P.rcParams.update({'font.size':16})
fig = P.figure(figsize=(10,8), edgecolor='none', facecolor='white', dpi=80)

# Open each fits file and make plot
DR = sys.argv[1]
allfiles = glob.glob('{}stacks_clip/delta*.fits.gz'.format(DR))

for f in allfiles:
    try:
        a = fits.open(f)
        d = a[1].data
        h = a[1].header
        thid = h[23]
    except:
        print ("Couldn't open ", f)
        continue

    P.plot(10**d.LOGLAM1, d.DELTA1, lw = 2, color = 'royalblue', label=('Deltas of Stack'))
    P.plot(10**d.LOGLAM2, d.DELTA2, lw = 2, color = 'darkred', label=('Stack of Deltas'))
    P.xlabel('wavelength [Angstrom]')
    P.ylabel('Deltas')
    P.legend()
    P.title(r'Ly-$\alpha$ {} Deltas for QSO: {}'.format(DR, thid), color='black')
    P.savefig('{}deltas_clip/deltas_{}.png'.format(DR, thid))
    P.clf()

P.close()

