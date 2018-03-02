#! /usr/bin/env python
# The program to re-write the jla supernova data in sncosmo format
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from matplotlib import pyplot as plt
import re
import os
import sncosmo
import builtins_jla
from astropy.table import Table
#from ._extinction import ccm89

# wavelength limits for salt2 model
wl_min_sal = 3000
wl_max_sal = 7000

# wavelength limits for sugar model
wl_min_sug = 3341.41521
wl_max_sug = 8576.61898

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def bandpass_interpolators(name, radius, fig=False):

    # we'll figure out min and max wave as we go.
    minwave = float('inf')
    maxwave = 0.

    bi = sncosmo.bandpasses._BANDPASS_INTERPOLATORS.retrieve(name)
    r = radius
    band = bi.at(r)

    # update min,max wave
    minwave = min(minwave, band.minwave())
    maxwave = max(maxwave, band.maxwave())
    wave = np.linspace(band.minwave(), band.maxwave(), 1000)
    trans = band(wave)
    if fig == True:
        plt.plot(wave,trans,'r')
        plt.show()
    return wave, trans

# Extinction law (Cardelli, 1989 ApJ.345.245C) :
# delta_mag = A(lambda) = [a(x)Rv+b(x)]E(B-V)

def wave_eff_A(wave, trans, Rv=3.1):
    sn_name = 'lc-SDSS6108.list'
    EBV = read_lcSDSS(sn_name)['@MWEBV']
    z = read_lcSDSS(sn_name)['@Z_HELIO']
    rest_wave = wave*(1+np.array(z))
    spl = Spline1d(rest_wave, trans, k=1, ext = 1)
    dt = 10000
    xs = np.linspace(min(rest_wave), max(rest_wave), dt)
    dxs = ((max(rest_wave)-min(rest_wave))/dt)*np.ones(len(xs))
    w_eff =  (np.sum((10**(-ccm89(xs,Rv*EBV,Rv)/2.5))*spl(xs)*xs*dxs)/np.sum((10**(-ccm89(xs,Rv*EBV,Rv)/2.5))*spl(xs)*dxs))/(1+z)
    return w_eff

def A_l(xx, EBV, Rv=3.1):
    A_l = []
    xx = np.array(xx)
    x = 1./(xx/10**4.) # xx in AA
    for i in x:
        if i < 1.1:
            y = i**1.61
            a = 0.574*y
            b = -0.527*y
        elif 1.1 <= i <= 3.3:
            y = i - 1.82
            a = 1 + y*(0.17699 + y*(-0.50447 + y*(-0.02427+ y*(0.72085 + y*(0.01979 + y*(-0.77530 + 0.32999*y))))))
            b =  y*(1.41338 + y*(2.28305 + y*(1.07233 +y*(-5.38434 + y*(-0.62251 + y*(5.30260 - 2.09002*y))))))
        elif 3.3 < i < 5.9:
            a = 1.752 - 0.316*i - 0.104/((i-4.67)*(i-4.67) + 0.341)
            b = -3.09 + 1.825*i + 1.206/((i-4.62)*(i-4.62) + 0.263)
        elif 5.9 <= i <= 8.:
            Fa = -(i-5.9)*(i-5.9) * (0.04473 + 0.009779*(i-5.9))
            a = 1.752 - 0.316*i - 0.104/((i-4.67)*(i-4.67) + 0.341) + Fa
            Fb = (i-5.9)*(i-5.9) * (0.213 + 0.1207*(i-5.9))
            b = -3.09 + 1.825*i + 1.206/((i-4.62)*(i-4.62) + 0.263) + Fb
        A = (a*Rv + b)*EBV
        A_l.append(A)
    A_l = np.array(A_l)
    return A_l

def wl_cut_salt2(fname, EBV, z):
    filt = sncosmo.get_bandpass(fname)
    wlen = filt.wave
    tran = filt.trans
    dt = 10000
    spl = Spline1d(wlen, tran, k=1, ext = 1)

    xs = np.linspace(min(wlen), max(wlen), dt)
    dxs = ((max(wlen)-min(wlen))/(dt-1))
    wlen_eff = np.sum((10**(-A_l(xs,EBV)/2.5))*spl(xs)*xs*dxs)/np.sum((10**(-A_l(xs,EBV)/2.5))*spl(xs)*dxs)
    if wl_min_sal >= wlen_eff/(1+z) or wlen_eff/(1+z) >= wl_max_sal:
        return ('False', wlen_eff/(1+z))
    else:
        return('True', wlen_eff/(1+z))

def wl_cut_sugar(fname, z):
        filt = sncosmo.get_bandpass(fname)
        wlen = filt.wave
        tran = filt.trans
        dt = 10000
        wlen_shift = wlen/(1+z)
        spl = Spline1d(wlen_shift, tran, k=1, ext = 1)
        xs = np.linspace(min(wlen_shift), max(wlen_shift), dt)
        dxs = ((max(wlen_shift)-min(wlen_shift))/(dt-1))
        area_full = np.sum(spl(xs)*dxs) # full area under the filter

        xs_cut = np.linspace(wl_min_sug, wl_max_sug, dt)
        dxs_cut = ((wl_max_sug-wl_min_sug)/(dt-1))
        area_cut = np.sum(spl(xs_cut)*dxs_cut)
        r = 1.-area_cut/area_full # area under the filter outside the model

        if r < 0.1:
            return ('True',r)
        else:
            return ('False',r)

def read_lc_jla(sn_name, model = None):
    infile = open('jla_data/jla_light_curves/'+ sn_name, 'r') # download the sn data
    photometry = []
    time = []
    band = []
    flux = []
    fluxerr = []
    zp = []
    zpsys = []
    head = {}
    for line in infile:
        if line[0] == '@':
            d = line.split()
            if is_number(d[1]):
                head[d[0]] = float(d[1])
            else:
                head[d[0]] = d[1]
            continue
        elif line[0] == '#': # miss the lines started from #
            continue
        elif len(line) == 1: # miss empty lines
            continue
        photometry = line.split()
        time.append(float(photometry[0]))
        flux.append(float(photometry[1]))
        fluxerr.append(float(photometry[2]))
        zp.append(float(photometry[3]))
        band.append('jla_' + photometry[4])
        zpsys.append('jla_' + photometry[5])
    infile.close()

    data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})

    cov = 'covmat_' + sn_name.rsplit('.')[0] + '.dat'
    if cov in os.listdir('jla_data/jla_light_curves/'):
        # with open('jla_data/jla_light_curves/' + cov, 'r') as table:
        #     size = table.readline()
        #     cov_file = np.genfromtxt(table)

        infile = open('jla_data/jla_light_curves/' + cov, 'r')
        size = infile.readline()
        table = [line.encode('utf-8') for line in infile]
        cov_file = np.genfromtxt(table)
        
        data['fluxcov'] = cov_file

    dic = {}
    for x in data:
        if x[1] in dic.keys():
            dic[x[1]] += 1
        else:
            dic[x[1]] = 1

    if '@X_FOCAL_PLANE' in head.keys():
        radius = np.sqrt(head['@X_FOCAL_PLANE']**2. + head['@Y_FOCAL_PLANE']**2.)

    for fname in dic.keys():
        if fname.startswith('jla_MEGACAMPSF::'):
            name = fname[4:]
            filt = bandpass_interpolators(name,radius)
            wlen = filt[0]
            tran = filt[1]
            band = sncosmo.Bandpass(wlen, tran, name=fname)
            sncosmo.registry.register(band, force=True)

    bands_ab = {'jla_SDSS::u': ('jla_AB_B12_0',  0.06791),
            'jla_SDSS::g': ('jla_AB_B12_0', -0.02028),
            'jla_SDSS::r': ('jla_AB_B12_0', -0.00493),
            'jla_SDSS::i': ('jla_AB_B12_0', -0.01780),
            'jla_SDSS::z': ('jla_AB_B12_0', -0.01015),
            'jla_MEGACAMPSF::u': ('jla_AB_B12_0',  0),
            'jla_MEGACAMPSF::g': ('jla_AB_B12_0', 0),
            'jla_MEGACAMPSF::r': ('jla_AB_B12_0', 0),
            'jla_MEGACAMPSF::i': ('jla_AB_B12_0', 0),
            'jla_MEGACAMPSF::z': ('jla_AB_B12_0', 0)}
    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_ab),'jla_AB_B12', force=True)


    f_in = {}
    f_out = {}
    for i in dic.keys():
        if model == 'salt2':
            res = wl_cut_salt2(i, head['@MWEBV'], head['@Z_HELIO'])
            if res[0] == 'True':
                f_in[i] = res[1]
            else:
                f_out[i] = res[1]
                #print('We excluded passband %s (%d points) because restframewavelength = %7.3f does not belong to the interval [%d,%d]' % (i,dic[i],res[1],wl_min_sal,wl_max_sal))
        elif model == 'sugar':
            res = wl_cut_sugar(i, head['@Z_HELIO'])
            if res[0] == 'True':
                f_in[i] = res[1]
            else:
                f_out[i] = res[1]
                print('We excluded passband %s (%d points) because it does not belong to the interval [%d,%d]' % (i,dic[i],wl_min_sug,wl_max_sug))
        else:
            print('ERROR: model name has to be salt2 or sugar')

    mask = []
    for row in data:
        if row[1] in f_in.keys():
            mask.append(True)
        else:
            mask.append(False)
    mask = np.array(mask)

    data_cut = sncosmo.select_data(data, mask)

    return head, data_cut
