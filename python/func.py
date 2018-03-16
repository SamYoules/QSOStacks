from astropy.io import fits
from picca import prep_del, io
from picca.data import forest, delta
import fitsio
import numpy as N
import scipy as sp
from scipy.interpolate import interp1d
import sys

def init():
    # Initialise
    forest.lmin = sp.log10(3600.)
    forest.lmax = sp.log10(5500.)
    forest.lmin_rest = sp.log10(1040.)
    forest.lmax_rest = sp.log10(1200.)
    forest.rebin = 3
    forest.dll = 3*1e-4
    forest.dla_mask = 0.8

    a = fits.open('iter.fits.gz')
    fo = open('rejects.txt', 'a')

    iter_data = a[2].data
    forest.var_lss = interp1d(iter_data.loglam, iter_data.var_lss, fill_value = 'extrapolate', kind = 'nearest')
    forest.eta = interp1d(iter_data.loglam, iter_data.eta, fill_value = 'extrapolate', kind = 'nearest')
    forest.fudge = interp1d(iter_data.loglam, iter_data.fudge, fill_value = 'extrapolate', kind = 'nearest')
    iter_data = a[3].data
    forest.mean_cont = interp1d(iter_data.loglam_rest, iter_data.mean_cont, fill_value = 'extrapolate', kind = 'nearest')


def read_from_spec(in_dir,thid,ra,dec,zqso,plate,mjd,fid,order,mode):
    # Open spec fits file and write data to pix_data
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
    npix_min = 50
    pix_data = read_from_spec(in_dir,thid, ra, dec, z, plate, mjd, fiberid, order, mode=mode)

    if not pix_data is None:
        print("{} read from pix\n".format(len(pix_data)))
    if not pix_data is None and len(pix_data)>0:
        data = pix_data
    else:
        sys.exit()

    # strip out short forests and non-numeric values from data
    l = []
    for d in data:
        if not hasattr(d,'ll') or len(d.ll) < npix_min:
            print "{} forest too short\n".format(d)
            sys.exit()

        if sp.isnan((d.fl*d.iv).sum()):
            print "{} nan found\n".format(d)
            sys.exit()

        l.append(d)

    data = N.asarray(l)
    return data

def stack_flux(data, delta):
    '''Make a weighted sum of flux/delta values in wavelength bins.'''

    nstack = int((forest.lmax - forest.lmin) / forest.dll) + 1
    ll = forest.lmin + sp.arange(nstack) * forest.dll
    st = sp.zeros(nstack)
    wst = sp.zeros(nstack)
    data_bad_cont = []

    # Stack flux & weights, or deltas & weights
    for d in data:
        if d.bad_cont is not None:
            data_bad_cont.append(d)
            continue

        bins=((d.ll - d.lmin) / d.dll + 0.5).astype(int)
        eta = forest.eta(d.ll)
        var_lss = forest.var_lss(d.ll)
        fudge = forest.fudge(d.ll)

        if (delta == 0):
            # convert ivar into normalized ivar (going from flux units to F units)
            ivar_F = d.iv * d.co**2

            # correct this variance, adding the var_lss and eta factors
            var_F = 1./ivar_F
            var_F_tot = var_F*eta + var_lss + fudge/var_F

            # convert back to flux units
            var_flux_tot = var_F_tot * d.co**2 
            we = 1./var_flux_tot
            c = sp.bincount(bins, weights = d.fl * we)
        else:
            iv = d.iv / eta
            we = iv * d.co**2 / (iv * d.co**2 * var_lss + 1)
            c = sp.bincount(bins, weights = (d.fl/d.co - 1) * we)

        st[:len(c)] += c
        c = sp.bincount(bins, weights = we)
        wst[:len(c)] += c

    w = wst>0
    st[w] /= wst[w]
    for d in data_bad_cont:
        print ("rejected {} due to {}\n".format(d.thid,d.bad_cont))

    return ll, st, wst

