# Date: 11 Jan 2018
# Author: Sam Youles
# Package: picca
# make_qsoinfo.py
# Create qsoinfo.txt file containing: plate, mjd, fiberid, thingid, RA, Dec, Z
# Spectra obtained from https://dr14.sdss.org/optical/spectrum/search
#        and https://dr12.sdss.org/optical/spectrum/search
#        and https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/test/bautista/v5_10_7/spectra/lite/
# Z found from DR14Q_v3_0.fits
# Purpose: Allows multiple spectra for the same object (THING_ID) to be stacked

import glob
from astropy.io import fits
import numpy as N
import sys

DR = sys.argv[1]

allfiles = glob.glob('{}/spec*.fits'.format(DR))
fout = open('{}qsoinfo.txt'.format(DR), 'w')

m = 0
drq = fits.open('../picca/data/DR14Q_v3_0.fits')[1].data
for f in allfiles:
    try:
        tab = fits.open(f)[2].data
    except:
        print ("Couldn't open ", f)
        continue
    try:
        j = N.where(drq.THING_ID == tab.THING_ID[0])
        zqso = drq.Z[j][0]
    except:
        zqso = tab.Z[0]
    print >> fout, tab.THING_ID[0], tab.PLATE[0], tab.MJD[0], tab.FIBERID[0], tab.RA[0], tab.DEC[0], zqso
    # Counter to see progress
    m +=1
    if (m % 100 == 0):
        print (m)

fout.close()
