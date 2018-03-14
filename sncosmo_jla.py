#! /usr/bin/env python
import sncosmo
import builtins_jla
import os
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from read_jla import read_lc_jla
import numpy as np
import copy

t_min = -15
t_max = 45

t_min_sug = -12
t_max_sug = 42

wl_min_sal = 3000.
wl_max_sal = 7000.

def fit_salt2(results=False):
    dic = {}
    jla_file = np.loadtxt('./jla_data/jla_lcparams.txt',dtype='str')
    for line in jla_file:
        dic['lc-' + line[0] + '.list'] = float(line[1]), float(line[20])

    outfile = open('res_salt22.txt', 'w')
    outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor')
    # list_jla = ['lc-03D1co.list','lc-04D4jy.list','lc-05D3ax.list','lc-03D1dt.list','lc-04D3dd.list','lc-06D2cd.list','lc-03D3cd.list','lc-05D3ha.list']
    list_jla = os.listdir('jla_data/jla_light_curves/')
    fitfail = []
    for filename in list_jla:
        if filename.startswith('lc-'):
            sn_name = filename
            print(sn_name)
            head, data = read_lc_jla(sn_name, model = 'salt2')
            source = sncosmo.get_source('salt2', version='2.4')
            source.EBV_snfit = head['@MWEBV']
            source.z_snfit = head['@Z_HELIO']
            source.Rv_snfit = 3.1

            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])
            model.set(mwebv=head['@MWEBV'])
            model.set(z=head['@Z_HELIO'])

            try:
                zperr_U = head['@ZPERR_STANDARD::U']
                if zperr_U == 0.1:
                    apply_ZPERR = True
                else:
                    apply_ZPERR = False
            except:
                apply_ZPERR = False

            try:
                # First iteration (x1 is fixed)
                model.set(x1=0.01)
                res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'c'], modelcov=False, guess_t0=True, apply_ZPERR=False)
                chi2 = res.chisq
                chi2p = chi2+1
                m = 0
                while chi2 < chi2p:
                    if m > 0:
                        resp = res
                        fitted_modelp = fitted_model

                    t_peak = fitted_model.parameters[1]

                    t1 = t_peak + t_min*(1 + model.get('z'))
                    t2 = t_peak + t_max*(1 + model.get('z'))

                    A=[]
                    data_new = copy.deepcopy(data)
                    for i in range(len(data_new)):
                        if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                            A.append(i)
                    A=np.array(A)
                    for i in range(len(A)):
                        #print('We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2))
                        data_new.remove_row(A[i])
                        A-=1

                    res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, apply_ZPERR=apply_ZPERR)
                    chi2p = chi2
                    chi2 = res.chisq
                    m += 1

                # Final results
                res=resp
                fitted_model=fitted_modelp

                if results == True:
                    print(res, '\nNumber of iterations: ',m)
                    sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
                    #plt.savefig(sn_name[:9]+'.pdf')
                    plt.show()

                # Calculation of m_b
                source2 = sncosmo.get_source('salt2', version='2.4')
                model2 = sncosmo.Model(source=source2)
                model2.set(t0=res.parameters[1], x0=res.parameters[2], x1=res.parameters[3], c=res.parameters[4])
                mb_without_corr = model2.bandmag('jla_STANDARD::B', 'jla_VEGA2_mb', [res.parameters[1]])
                zcmb = dic[sn_name][0]
                dbmag_pecvel = -5*np.log10((1+head['@Z_HELIO'])/(1+zcmb))
                biascor = dic[sn_name][1]

                mb = mb_without_corr + biascor + dbmag_pecvel

                # m_b uncertainty
                dmbfit, cov_mb_c, cov_mb_x1, cov_mb_x0 = mb_uncertainty(res)
#                print(dmbfit, cov_mb_c, cov_mb_x1, cov_mb_x0)


                outfile.write('\n %s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5], res.parameters[0], mb, dbmag_pecvel, res.parameters[3], res.errors['x1'], res.parameters[4], res.errors['c'], res.parameters[1], res.errors['t0']))

            except:
                fitfail.append(sn_name)
                print('Error: fit fail for: ',sn_name)

            print(fitfail)

    outfile.close()

def fit_sugar(results=False):
    dic = {}
    # outfile = open('res_sugar.txt', 'w')
    # outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor')
#     list_jla = ['lc-03D1co.list','lc-04D4jy.list','lc-05D3ax.list','lc-03D1dt.list','lc-04D3dd.list','lc-06D2cd.list','lc-03D3cd.list','lc-05D3ha.list']
    list_jla = ['lc-sn2006n.list']
    jla_file = np.loadtxt('./jla_data/jla_lcparams.txt',dtype='str')
    for line in jla_file:
        dic['lc-' + line[0] + '.list'] = float(line[12])
#    list_jla = os.listdir('jla_data/jla_light_curves/')
    fitfail = []
    for filename in list_jla:
#        if filename.startswith('lc-SDSS'):
            sn_name = filename
            print(sn_name)
            head, data = read_lc_jla(sn_name, model = 'sugar')
            source = sncosmo.get_source('sugar')
            source.EBV_snfit = head['@MWEBV']
            source.z_snfit = head['@Z_HELIO']
            source.Rv_snfit = 3.1
            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])
            model.set(mwebv=head['@MWEBV'])
            model.set(z=head['@Z_HELIO'])

            try:
                zperr_U = head['@ZPERR_STANDARD::U']
                if zperr_U == 0.1:
                    apply_ZPERR = True
                else:
                    apply_ZPERR = False
            except:
                apply_ZPERR = False

            try:            
                # First iteration (x1 is fixed)
                model.set(q1=1)
                res, fitted_model = sncosmo.fit_lc(data, model, ['t0','q1','q2', 'q3', 'A', 'Mgr'],
                                    bounds={'t0':(dic[sn_name]-5,dic[sn_name]+5)}, modelcov=False, guess_t0=True, apply_ZPERR=False)                
                # res, fitted_model = sncosmo.fit_lc(data, model, ['t0','q1','q2', 'q3', 'A', 'Mgr'],
                #                     bounds={'t0':(head['@DayMax']-5,head['@DayMax']+5)}, modelcov=False, guess_t0=True,phase_range=(-12,42), apply_ZPERR=False)
    
                    
                chi2 = res.chisq
                chi2p = chi2+1
                m = 0
                while chi2 < chi2p:
                    if m > 0:
                        resp = res
                        fitted_modelp = fitted_model
    
                    t_peak = fitted_model.parameters[1]
    
                    t1 = t_peak + t_min*(1 + model.get('z'))
                    t2 = t_peak + t_max*(1 + model.get('z'))
    
                    A=[]
                    data_new = copy.deepcopy(data)
                    for i in range(len(data_new)):
                        if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                            A.append(i)
                    A=np.array(A)
                    for i in range(len(A)):
                        # print('We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2))
                        data_new.remove_row(A[i])
                        A-=1
    
                    res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0', 'q1', 'q2', 'q3', 'A', 'Mgr'], modelcov=True, apply_ZPERR=False)
                    chi2p = chi2
                    chi2 = res.chisq
                    m += 1

                # Final results
                res=resp
                fitted_model=fitted_modelp

                if results == True:
                    print(res, '\nNumber of iterations: ',m)
                    sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
                    #plt.savefig(sn_name[:9]+'.pdf')
                    plt.show()

                # outfile.write('\n %s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5], res.parameters[0], mb, dbmag_pecvel, res.parameters[3], res.errors['x1'], res.parameters[4], res.errors['c'], res.parameters[1], res.errors['t0']))

            except:
                fitfail.append(sn_name)
                print('Error: fit fail for: ',sn_name)

            print(fitfail)

    # outfile.close()

def mB_calc():
    #interpolation of TB and Trest
    filt=np.genfromtxt('jla_data/Instruments/SNLS3-Landolt-model/sb-shifted.dat')
    wlen=filt[:,0]
    tran=filt[:,1]
    splB = Spline1d(wlen, tran, k=1, ext=1)

    #interpolation of ref spectrum (BD17)
    data = np.genfromtxt('jla_data/MagSys/bd_17d4708_stisnic_003.ascii')
    wlen_ref = data[:,0]
    flux_ref = data[:,1]
    splref = Spline1d(wlen_ref, flux_ref, k=1, ext=1)

    plt.plot(wlen_ref,flux_ref,'ro',np.arange(1000.,10000.,1.),splref(np.arange(1000.,10000.,1.)))
    plt.show()
#return mb

def chi2(name=None):
    sn_name='lc-SDSS762.list'
    head, data = read_lc_jla(sn_name, model = 'salt2')
    source = sncosmo.get_source('salt2', version='2.4')
    source.EBV_snfit = head['@MWEBV']
    source.z_snfit = head['@Z_HELIO']
    source.Rv_snfit = 3.1

    dust = sncosmo.CCM89Dust()

    model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
    model.set(mwebv=head['@MWEBV'], z=head['@Z_HELIO'], t0=53625.149, x0=0.0001,x1=1.1259,c=-0.0388) #

    t_peak = model.parameters[1]

    t1 = t_peak + t_min*(1 + model.get('z'))
    t2 = t_peak + t_max*(1 + model.get('z'))

    A=[]
    data_new = copy.deepcopy(data)
    for i in range(len(data_new)):
        if data_new[i][0] <= t1 or data_new[i][0] >= t2:
            A.append(i)
    A=np.array(A)
    for i in range(len(A)):
        #print('We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2))
        data_new.remove_row(A[i])
        A-=1
    data_new = sncosmo.fitting.photometric_data(data_new)
    chi = sncosmo.fitting.generate_chisq(data_new, model, modelcov=False, apply_ZPERR=False) # to be careful with z what 'zperr_u = True'
    #print(model.bandflux('jla_STANDARD::U', [51072.12,51077.14,51078.14,51083.23,51136.1], zp=14.205682, zpsys='jla_VEGA2'))
    #for i in range(13):
    #   print(model.bandflux(data[i][1], float(data[i][0]), zp=float(data[i][4]), zpsys='AB_jla'))
    sncosmo.plot_lc(data_new, model=model)
    plt.savefig('lc-SDSS14318_snfit.pdf')
    plt.show()

def results(filename):
    salt2_param = open(filename, 'r')
    lines = salt2_param.readlines()
    salt2_param.close()
    data = {}
    first_line = lines[0]
    properties = first_line.split()

    for line in lines[1:]:
        words = line.split()
        i = words[0]
        values = words[1:]
        data[i] = {}
        for p, v in zip(properties[1:], values):
            data[i][p] = float(v)
    return data

def comparison_plot(par='x1', er_par='dx1'):
    # Plot Setup
    rcParams['font.size'] = 16.
    font = {'family': 'normal', 'size': 16}
    rc('axes', linewidth=1.5)
    rc("text", usetex=True)
    rc('font', family='serif')
    rc('font', serif='Times')
    rc('legend', fontsize=25)
    rc('xtick.major', size=5, width=1.5)
    rc('ytick.major', size=5, width=1.5)
    rc('xtick.minor', size=3, width=1)
    rc('ytick.minor', size=3, width=1)
    fig = plt.figure(figsize=(8.,8.))
    x = 0
    for key in results('res_salt2.txt'):
        dif = results('res_salt2.txt')[key][par]-results('jla_data/jla_lcparams.txt')[key][par]
        dif_er = results('res_salt2.txt')[key][er_par] - results('jla_data/jla_lcparams.txt')[key][er_par]
        #if abs(dif) > np.sqrt(results('res_salt2.txt')[key][er_par]**2. + results('jla_data/jla_lcparams.txt')[key][er_par]**2.):
        #if abs(dif) > results('jla_data/jla_lcparams.txt')[key][er_par]:
        #if -0.01 < dif_er < 0.01:
        #    print(key)
        plt.errorbar(results('res_salt2.txt')[key][par],results('jla_data/jla_lcparams.txt')[key][par],xerr=results('res_salt2.txt')[key][er_par],yerr=results('jla_data/jla_lcparams.txt')[key][er_par], color='red', fmt='o', mfc='red', zorder=1)
        #plt.errorbar(results('res_salt2.txt')[key][par],results('jla_data/jla_lcparams.txt')[key][par],xerr=None,yerr=results('jla_data/jla_lcparams.txt')[key][er_par], color='red', fmt='o', mfc='red', zorder=1)

        x = x + 1

    ax = plt.subplot(111)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xlim([-4,4])
    ax.set_ylim([-4,4])
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls = '--', color = 'black', alpha = 0.8)
    plt.ylabel('snfit',fontsize=25)
    plt.xlabel('sncosmo',fontsize=25)
    plt.figtext(0.2, 0.8, par, fontsize=25)
    pdffile = 'x1_plot.pdf'
    plt.savefig(pdffile, bbox_inches='tight')
    #plt.savefig('plot.png')
    plt.show()

def comparison_hist(par='x1', er_par='dx1'):
    # Plot Setup
    rcParams['font.size'] = 16.
    font = {'family': 'normal', 'size': 16}
    rc('axes', linewidth=1.5)
    rc("text", usetex=True)
    rc('font', family='serif')
    rc('font', serif='Times')
    rc('legend', fontsize=25)
    rc('xtick.major', size=5, width=1.5)
    rc('ytick.major', size=5, width=1.5)
    rc('xtick.minor', size=3, width=1)
    rc('ytick.minor', size=3, width=1)
    #fig = plt.figure(figsize=(8.,8.))

    gs = gridspec.GridSpec(2, 1) #subplots ratio
    f, (ax1, ax2) = plt.subplots(2, sharex=True)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1, 0])
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    dif_par = []
    dif_er_par = []
    for key in results('res_salt2.txt'):
        dif_par.append(results('res_salt2.txt')[key][par]-results('jla_data/jla_lcparams.txt')[key][par])
        dif_er_par.append(results('res_salt2.txt')[key][er_par]-results('jla_data/jla_lcparams.txt')[key][er_par])
        #if abs(dif_par) > np.sqrt(results('res_salt2.txt')[key][er_par]**2. + results('jla_data/jla_lcparams.txt')[key][er_par]**2.):
           #print(key)
    ax1.hist(dif_par,25, color = 'red')
    ax2.hist(dif_er_par,25, color = 'red')

    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax1.set_ylabel('N',fontsize=18)
    ax2.set_ylabel('N',fontsize=18)
    ax1.set_xlabel(par,fontsize=18)
    ax2.set_xlabel(er_par,fontsize=18)
    ax1.set_yscale('log',nonposy='clip')
    ax2.set_yscale('log',nonposy='clip')
    ax1.set_ylim(bottom=0.1)
    ax2.set_ylim(bottom=0.1)
    #pdffile = 't_hist.pdf'
    #plt.savefig(pdffile, bbox_inches='tight')
    plt.show()

#########################
# Calculation of mb without sncosmo

def mB_determination(res):

    scale_factor = 10**-12

    #interpolation of TB and Trest
    filt2 = np.genfromtxt('jla_data/Instruments/SNLS3-Landolt-model/sb-shifted.dat')
    wlen = filt2[:,0]
    tran = filt2[:,1]
    splB = Spline1d(wlen, tran, k=1,ext = 1)



    #interpolation of ref spectrum
    data = np.genfromtxt('jla_data/MagSys/bd_17d4708_stisnic_002.ascii')
    dispersion = data[:,0]
    flux_density = data[:,1]
    splref = Spline1d(dispersion, flux_density, k=1,ext = 1)


    #interpolation of the spectrum model
    template_0 = np.genfromtxt('jla_data/salt2-4/salt2_template_0.dat')
    template_1 = np.genfromtxt('jla_data/salt2-4/salt2_template_1.dat')
#    salt2source=sncosmo.SALT2Source('/users/divers/lsst/mondon/hubblefit/sncosmo_jla/salt2-4')

    wlM0 = []
    M0 = []
    for i in range(len(template_0[:,0])):
        if template_0[:,0][i] == 0.0:
             wlM0.append(template_0[:,1][i])
             M0.append(template_0[:,2][i])
    splM0 = Spline1d(wlM0, M0, k=1,ext = 1)

    wlM1 = []
    M1 = []
    for i in range(len(template_1[:,0])):
        if template_1[:,0][i] == 0.0:
            wlM1.append(template_1[:,1][i])
            M1.append(template_1[:,2][i])
    splM1 = Spline1d(wlM1, M1, k=1,ext = 1)

    #computation of the integral
    dt = 100000
    xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
    dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))

#    I1=np.sum((splM0(xs)*10**-12+res.parameters[3]*splM1(xs)*10**-12)*(10**(-0.4*salt2source.colorlaw(xs)*res.parameters[4]))*xs*splB(xs)*dxs)
    I1 = np.sum((splM0(xs)*scale_factor + res.parameters[3]*splM1(xs)*scale_factor)*(10**(-0.4*color_law_salt2(xs)*res.parameters[4]))*xs*splB(xs)*dxs)
    I2 = np.sum(splref(xs)*xs*splB(xs)*dxs)
#    print(I1, I2)



    #computation of mb
    mref = 9.907
    mb = -2.5*np.log10(res.parameters[2]*(I1/I2))+mref

    return mb

def color_law_salt2(wl):
        B_wl = 4302.57
        V_wl = 5428.55
        l = (wl-B_wl)/(V_wl-B_wl)
        l_lo = (2800.-B_wl)/(V_wl-B_wl)
        l_hi = (7000.-B_wl)/(V_wl-B_wl)
        a = -0.504294
        b = 0.787691
        c = -0.461715
        d = 0.0815619
        cst = 1-(a+b+c+d)
        cl = []
        for i in range (len(l)):
            if l[i] > l_hi:
                cl.append(-(cst*l_hi+l_hi**2*a+l_hi**3*b+l_hi**4*c+l_hi**5*d+(cst+2*l_hi*a+3*l_hi**2*b+4*l_hi**3*c+5*l_hi**4*d)*(l[i]-l_hi)))
            if l[i] < l_lo:
                cl.append(-(cst*l_lo+l_lo**2*a+l_lo**3*b+l_lo**4*c+l_lo**5*d+(cst+2*l_lo*a+3*l_lo**2*b+4*l_lo**3*c+5*l_lo**4*d)*(l[i]-l_lo)))
            if l[i]>= l_lo and l[i]<= l_hi:
                cl.append(-(cst*l[i]+l[i]**2*a+l[i]**3*b+l[i]**4*c+l[i]**5*d))
        return np.array(cl)

def mb_uncertainty(res):

    h = 10**-9

    #build mb derivative for all parameters
    resx0 = copy.deepcopy(res)
    resx0.parameters[2] = resx0.parameters[2] + h
    dmb_dx0 = (mB_determination(resx0)- mB_determination(res))/h

    resx1 = copy.deepcopy(res)
    resx1.parameters[3] = resx0.parameters[3] + h
    dmb_dx1 = (mB_determination(resx1)- mB_determination(res))/h

    resc = copy.deepcopy(res)
    resc.parameters[4] = resx0.parameters[4] + h
    dmb_dc = (mB_determination(resc)- mB_determination(res))/h

    #build vectors for all mb derivative
    vect = np.array([dmb_dx0, dmb_dx1, dmb_dc])

    #build the covariance matrix for salt2 parameters
    mat = np.delete(res.covariance, (0), axis=0)
    mat = mat.T
    mat = np.delete(mat, (0), axis=0)
    mat = mat.T
    print(len(mat))

    dmb = np.sqrt(np.dot(np.dot(vect.T, mat), vect))
    cov_mb_c = mat[2,2] * dmb_dc + mat[2,1] * dmb_dx1 + mat[2,0] * dmb_dx0
    cov_mb_x1 = mat[1,1] * dmb_dx1 + mat[1,0] * dmb_dx0 + mat[1,2] * dmb_dc
    cov_mb_x0 = mat[0,0] * dmb_dx0 +  mat[0,1] * dmb_dx1 + mat[0,2] * dmb_dc

    return dmb, cov_mb_c, cov_mb_x1, cov_mb_x0
