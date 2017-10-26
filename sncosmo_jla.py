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
#from decimal import Decimal







t_min = -15
t_max = 45

t_min_sug = -12
t_max_sug = 42

# wavelength limits for salt2 model
wl_min_sal = 3000
wl_max_sal = 7000
    
def fit_salt2(results=False):
    jla_file = np.loadtxt('./jla_data/jla_lcparams.txt',dtype='str')
    jla_zcmb = np.array(jla_file[:,1],float)
    jla_sn_name = np.array(jla_file[:,0])
    jla_biascor = np.array(jla_file[:,20],float)
#    outfile = open('results/res_salt2.txt', 'w')
#    outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor x0 dx0')
#    list_SDSS = ['lc-03D1au.list']
    list_SDSS = ['lc-sn1998dx.list']    
#    list_SDSS = ['lc-SDSS13425.list']
#    list_SDSS = ['lc-SDSS14318.list']
#    list_SDSS = [raw_input('Name of SN: ')]    
    
    
    fitfail=[]   
    zperr=[]     
#    for filename in os.listdir('./jla_data/jla_light_curves/'):
    for filename in list_SDSS:
#        if 'lc-SDSS' == filename[:7] or 'lc-sn'==filename[:5]:
#        if 'lc-SDSS' == filename[:7]:   
        if filename.startswith('lc-'):
            sn_name = filename
            
            print sn_name
            head, data = read_lc_jla(sn_name, model = 'salt2')            
            source = sncosmo.get_source('salt2', version='2.4')
            source.EBV_snfit = head['@MWEBV']    
            source.z_snfit = head['@Z_HELIO']
            source.Rv_snfit = 3.1
    
            dust = sncosmo.CCM89Dust()
            #model = sncosmo.Model(source=source)
            model = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])
            model.set(mwebv=head['@MWEBV'])
            model.set(z=head['@Z_HELIO'])
           
#            
            try:
            #note for me: remove shift if you put errors 
                zperr_U = head['@ZPERR_STANDARD::U']
                if zperr_U == 0.1:
                    apply_ZPERR = True
                    zperr.append(sn_name)
                else: 
                    apply_ZPERR = False
            except:
                apply_ZPERR = False
           

            try:     
                #initial iteration x1 fix
                model.set(x1=0.0) #with 0.01 the  fit work for all SDSS and nearby   
                    
                print 'initialisation'
                res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'c'], modelcov = False, guess_t0=True,apply_ZPERR=False)
                print 'first iteration'
                chi2 = res.chisq
                chi2p = chi2 +1
                m=0
                while chi2 < chi2p:
                        if m > 0:
                            resp = res
                            fitted_modelp = fitted_model
    
                        t_peak = fitted_model.parameters[1]
                        #print t_peak,fitted_model.parameters[4]
            
                        t1 = t_peak + t_min*(1 + model.get('z'))
                        t2 = t_peak + t_max*(1 + model.get('z'))
                                    
                        A=[]
                        data_new = copy.deepcopy(data)
                        for i in range(len(data_new)):                    
                            if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                                A.append(i)
                        A=np.array(A)
                        for i in range(len(A)):
                            #print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2)
                            data_new.remove_row(A[i])
                            A-=1    
                        #model.set(t0=t_peak)
                        #er_t_peak = res.errors['t0']
                        
                        #res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'x1', 'c'], modelcov = True, guess_t0=True)
                                          
#                        res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, bounds={'t0':(t_peak-5, t_peak+5),'x0':(x0_guess-5, x0_guess+5),'x1':(x1_guess-5, x1_guess+5),'c':(c_guess-5, c_guess+5)},apply_ZPERR=apply_ZPERR)       
                        res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, apply_ZPERR=apply_ZPERR) 
                        chi2p = chi2
    
                        chi2 = res.chisq
                        print chi2p, chi2
                        m += 1
                        
                #final results
                res = resp
                fitted_model = fitted_modelp
                print (res.chisq)
    
    
                    
                    #print fitted_model.get('t0')
                    #print fitted_model.maxtime(), fitted_model.mintime()
                    #print res
    
                plt.show()
                mB_without_corr = mB_determination(res)
                
                for i in range (len(jla_sn_name)):
                    
                    if sn_name == 'lc-' + jla_sn_name[i] + '.list':
                       biascor = jla_biascor[i]
                       zcmb = jla_zcmb[i]
                        
                mb = mB_without_corr + biascor - 5*np.log10((1 + head['@Z_HELIO'])/(1 + zcmb))
#                dmbfit, cov_mb_c, cov_mb_x1, cov_mb_x0 = mb_uncertainty(res)
#                cor_dmb = ()
                print mb
                
                
                if results == True:    
                    print res, '\nNumber of iterations: ', m
                    sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
                    #plt.savefig('lc.eps')
                    plt.show()
    
#                outfile.write('\n')        
#                outfile.write('%s 999 %f 999 %f 999 %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999 %e %e' %(sn_name[3:-5], res.parameters[0], mb, res.parameters[3], res.errors['x1'], res.parameters[4], res.errors['c'], res.parameters[1], res.errors['t0'], res.parameters[2], res.errors['x0']))          

#
            except:
                fitfail.append(sn_name)
                print 'Error: fit fail for: ',sn_name          
#            
#            
            print fitfail
            print zperr        
#    outfile.close()
    return res
    
    
def fit_sugar():

    #outfile = open('res_sugar.txt', 'w')
    for filename in os.listdir('/home/maria/Dropbox/Science/Supernovae/JLA_SALT2/JLA_fit/jla_data'):
    #list_SDSS = ['lc-SDSS18468.list']
    #list_SDSS = ['b_lc-SDSS14318.list','b_lc-SDSS17274.list','b_lc-SDSS20470.list','b_lc-SDSS21042.list','b_lc-SDSS16402.list','b_lc-SDSS11300.list','b_lc-SDSS16206.list','b_lc-SDSS20575.list','b_lc-SDSS16281.list','b_lc-SDSS18617.list','b_lc-SDSS17220.list']    
    #for filename in list_SDSS:
        if 'lc-SDSS' == filename[:7]:
            sn_name = filename
            print sn_name
            source = sncosmo.get_source('sugar')    
            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])
            
            head = read_lc_jla(sn_name)
            model.set(mwebv=head['@MWEBV'])
            model.set(z=head['@Z_HELIO'])
    
            data = sncosmo.read_lc('/home/maria/Dropbox/Science/Supernovae/sncosmo/test_jla/' + sn_name)
    
            print head['@DayMax']
            res, fitted_model = sncosmo.fit_lc(data, model, ['t0','q1','q2', 'q3', 'A', 'Mgr'], bounds={'t0':(head['@DayMax']-3,head['@DayMax']+3)},guess_t0=True,phase_range=(-12,42))

            print res
            sncosmo.plot_lc(data, model=fitted_model, errors=res.errors)
            plt.show()

            #outfile.write('%s \n' %(sn_name)) 
            #outfile.write('%s \n' %(res)) 
        
    #outfile.close()    

def chi2(sn_name):

    head, data = read_lc_jla(sn_name, model='salt2')

    source = sncosmo.get_source('salt2', version='2.4')
    source.EBV_snfit = head['@MWEBV']    
    source.z_snfit = head['@Z_HELIO']
    source.Rv_snfit = 3.1
    dust = sncosmo.CCM89Dust()
    
    #load snfit parameters 
    snfit = np.loadtxt('./x0_mb_JLA.txt', dtype='str')
    name_snfit = np.array(snfit[:,0])
    mb_snfit = np.array(snfit[:,3],float)
    dmb_snfit = np.array(snfit[:,4],float)
    x0_snfit = np.array(snfit[:,5],float)
    x1_snfit = np.array(snfit[:,7],float)
    c_snfit = np.array(snfit[:,9],float)
    t0_snfit = np.array(snfit[:,11],float)
    z_snfit = np.array(snfit[:,13],float)
    for i in range(len(name_snfit)):
        sn_name_snfit = 'lc-'+name_snfit[i]+'.list'
        if sn_name_snfit == sn_name:
            x0_snfit = x0_snfit[i]
            x1_snfit = x1_snfit[i]
            c_snfit = c_snfit[i]
            t0_snfit = t0_snfit[i]
            mb_snfit = mb_snfit[i]
            dmb_snfit = dmb_snfit[i]
            z_snfit = z_snfit[i]
            
    model = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])
    model.set(mwebv=head['@MWEBV'], z=head['@Z_HELIO'], t0=t0_snfit, x0=x0_snfit ,x1=x1_snfit,c=c_snfit) #

    print model.parameters
    print sncosmo.chisq(data, model, modelcov=True)
    return model, mb_snfit, dmb_snfit,z_snfit
#    print model.bandflux('jla_SDSS::g', [54277.5, 54285.0, 54292.5, 54300., 54307.5, 54315., 54322.5], zp=25, zpsys='jla_AB_B12')
    #for i in range(13):
    #    print model.bandflux(data[i][1], float(data[i][0]), zp=float(data[i][4]), zpsys='AB_jla')
#    sncosmo.plot_lc(data, model=model)
#    plt.show()
    
def chi2_maria(name=None):
#    sn_name='lc-sn1998dx.list'
    sn_name='lc-SDSS14318.list'
    head, data = read_lc_jla(sn_name, model = 'salt2')
    source = sncosmo.get_source('salt2', version='2.4')
    source.EBV_snfit = head['@MWEBV']
    source.z_snfit = head['@Z_HELIO']
    source.Rv_snfit = 3.1
    
    dust = sncosmo.CCM89Dust()

    model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
#    model.set(mwebv=head['@MWEBV'], z=head['@Z_HELIO'], t0=51071.92 ,x1=-1.5973806882,c=-0.118311996111, x0=0.00173409432723) #sn1998dx
    model.set(mwebv=head['@MWEBV'], z=head['@Z_HELIO'], t0=54072.073569 ,x1=0.864354,c=0.028353, x0=0.00138909129413) #SDSSS14318
#    res, fitted_model = sncosmo.fit_lc(data, model, [], modelcov = True, bounds={'x0':(0.001,0.002)} ,apply_ZPERR=True)

    t_peak = model.parameters[1]
    print model.parameters
    t1 = t_peak + t_min*(1 + model.get('z'))
    t2 = t_peak + t_max*(1 + model.get('z'))
                
    A=[]
    data_new = copy.deepcopy(data)
    for i in range(len(data_new)):                    
        if data_new[i][0] <= t1 or data_new[i][0] >= t2:
            A.append(i)
    A=np.array(A)
    for i in range(len(A)):
        #print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2)
        data_new.remove_row(A[i])
        A-=1 
#    print sncosmo.chisq(data,data,model,modelcov=True)

    #for i in range(13):
    # print model.bandflux(data[i][1], float(data[i][0]), zp=float(data[i][4]), zpsys='AB_jla')
    sncosmo.plot_lc(data_new, model=model)
    plt.savefig('lc-SDSS14318_snfit.pdf')
    plt.show()

def mB_determination(res):
    
    scale_factor = 10**-12
    
    #interpolation of TB and Trest
    filt2 = np.genfromtxt('./jla_data/Instruments/SNLS3-Landolt-model/sb-shifted.dat')
    wlen = filt2[:,0]
    tran = filt2[:,1]
    splB = Spline1d(wlen, tran, k=1,ext = 1)



    #interpolation of ref spectrum
    data = np.genfromtxt('./jla_data/MagSys/bd_17d4708_stisnic_002.ascii')
    dispersion = data[:,0]
    flux_density = data[:,1]
    splref = Spline1d(dispersion, flux_density, k=1,ext = 1)

  
    #interpolation of the spectrum model
    template_0 = np.genfromtxt('./salt2-4/salt2_template_0.dat')    
    template_1 = np.genfromtxt('./salt2-4/salt2_template_1.dat')
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
#    print I1, I2   
    
    
    
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
#        return -(cst*l+l**2*a+l**3*b+l**4*c+l**5*d)
def Plot_diff_mb():
    
    mb_sncosmo = np.array(fit_salt2())
    mb_snfittxt = np.loadtxt('/users/divers/lsst/mondon/hubblefit/sncosmo_jla/x0_mb_JLA.txt',dtype='str')
    mb_snfit = np.array(mb_snfittxt[:,3],float)
    dmb_snfit = np.array(mb_snfittxt[:,4],float)
    plt.errorbar(mb_sncosmo, mb_snfit, xerr=np.nan, yerr=dmb_snfit, marker='o', color='red')
    plt.xlabel('sncosmo')
    plt.ylabel('snfit')
    plt.show()    
    
def mb_using_snfit_parameters():
    diff = []
    z = []
    for filename in os.listdir('./jla_data/jla_light_curves/'):
#    for filename in list_SDSS:
        if 'lc-SDSS' == filename[:7]:
            sn_name = filename
            res,mb_snfit1,dmb_snfit1,z1 = chi2(sn_name)
            mb_snfit_code=mB_determination(res)
            diff.append(mb_snfit_code-mb_snfit1)
            z.append(z1)
    plt.scatter(z,diff,marker = '.',color='red',label='mb')


def comparison_plot(par='x1', er_par='dx1'):
    # Plot Setup
#    res_snfit = 'jla_data/jla_lcparams.txt'
    res_snfit = 'x0_mb_JLA.txt'
    res_sncosmo = 'results/res_salt2.txt'
#    res_sncosmo = 'results/res_salt2.0.txt'
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
    for key in results(res_sncosmo):
        dif = results(res_sncosmo)[key][par]-results(res_snfit)[key][par]
        if abs(dif) > np.sqrt(results(res_sncosmo)[key][er_par]**2. + results(res_snfit)[key][er_par]**2.):
            print key, results(res_sncosmo)[key][par], results(res_snfit)[key][par]
        plt.errorbar(results(res_sncosmo)[key][par],results(res_snfit)[key][par],xerr=results(res_snfit)[key][er_par],yerr=results(res_snfit)[key][er_par],marker = 'o',color='red')

        x = x + 1
    #plt.plot(range(x),np.zeros(x),'k')

    ax = plt.subplot(111)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
#    ax.set_xlim([-0.4,0.4])
#    ax.set_ylim([-0.4,0.4])
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls = '--', color = 'black', alpha = 0.8) 
    plt.ylabel('snfit',fontsize=25)
    plt.xlabel('sncosmo',fontsize=25)
    plt.figtext(0.2, 0.8, par, fontsize=25)

#    pdffile = 'x0_plot.pdf'
#    plt.savefig(pdffile, bbox_inches='tight')
#    plt.savefig(par+'_plot_f.png')
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
    
def comparison_plot_difference(par='x1', er_par='dx1'):
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
        dif = results('jla_data/jla_lcparams.txt')[key]['mb']+results('jla_data/jla_lcparams.txt')[key]['biascor']-results('res_salt2.txt')[key]['dmb']-results('res_salt2.txt')[key]['mb']
        dif2 = results('jla_data/jla_lcparams2.txt')[key]['mB']-results('res_salt2.txt')[key]['mb']
        #dif = results('jla_data/jla_lcparams.txt')[key]['mb']+results('jla_data/jla_lcparams.txt')[key]['biascor']-results('res_salt2.txt')[key]['dmb']-results('jla_data/jla_lcparams2.txt')[key]['mB']
        #dif = results('res_salt2.txt')[key][par]-results('jla_data/jla_lcparams.txt')[key][par]-results('jla_data/jla_lcparams.txt')[key]['biascor']
        #if abs(dif) > np.sqrt(results('res_salt2.txt')[key][er_par]**2. + results('jla_data/jla_lcparams.txt')[key][er_par]**2.):
        #if abs(dif) > 1e-5:
            #print key, results('jla_data/jla_lcparams.txt')[key]['mb'],results('jla_data/jla_lcparams.txt')[key]['biascor'],results('res_salt2.txt')[key]['dmb'],results('jla_data/jla_lcparams2.txt')[key]['mB']
        plt.plot(x,dif,marker = 'o',color='red')
        plt.plot(x,dif2,marker = '+',color='blue')

        x = x + 1

#   source = sncosmo.get_source('salt2', version='2.4')
#   model = sncosmo.Model(source=source)
#   for key in results('res_salt2.txt'):
#       t0 = results('jla_data/jla_lcparams2.txt')[key]['tmax']
#       x0 = results('jla_data/jla_lcparams2.txt')[key]['x0']
#       x1 = results('jla_data/jla_lcparams2.txt')[key]['x1']
#       c = results('jla_data/jla_lcparams2.txt')[key]['c']
#       model.set(t0=t0,x0=x0,x1=x1,c=c)
#       mb_sncosmo = model.bandmag('standard::b','jla1',[t0])
#       plt.plot(x,mb_sncosmo-results('jla_data/jla_lcparams2.txt')[key]['mBc'],marker = 'o',color='red')
#       x = x + 1

    ax = plt.subplot(111)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.ylabel('snfit-sncosmo',fontsize=25)
    plt.xlabel('n',fontsize=25)

    #plt.savefig('plot.png')
    plt.show()

def comparison_hist(par='tmax', er_par='dtmax'):
    # Plot Setup
#    res_snfit = 'jla_data/jla_lcparams.txt'
    res_snfit = 'x0_mb_JLA.txt'
#    res_sncosmo = 'results/res_salt2_without_zperr.txt'
    res_sncosmo = 'results/res_salt2.txt'
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
    for key in results(res_sncosmo):
        dif_par.append(results(res_sncosmo)[key][par]-results(res_snfit)[key][par])
        dif_er_par.append(results(res_sncosmo)[key][er_par]-results(res_snfit)[key][er_par])
        #if abs(dif_par) > np.sqrt(results('res_salt2.txt')[key][er_par]**2. + results('jla_data/jla_lcparams.txt')[key][er_par]**2.):
        #   print key
    ax1.hist(dif_par,25, color = 'red')
    ax2.hist(dif_er_par,25, color = 'red')

    ax1.xaxis.set_minor_locator(AutoMinorLocator(2))
#    ax2.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax2.xaxis.set_major_locator(MultipleLocator(0.00004))
    ax1.set_ylabel('N',fontsize=18)
    ax2.set_ylabel('N',fontsize=18)
    ax1.set_xlabel(par,fontsize=18)
    ax2.set_xlabel(er_par,fontsize=18)
    ax1.set_yscale('log',nonposy='clip')
    ax2.set_yscale('log',nonposy='clip')
    ax1.set_ylim(bottom=0.1)
    ax2.set_ylim(bottom=0.1)
#    pdffile = 'x0_hist.pdf'
#    plt.savefig(pdffile, bbox_inches='tight')
#    plt.savefig('x0_hist.pdf')
    plt.show()

  


    
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
#    mat = np.delete(res.covariance, (0), axis=0)
#    mat = mat.T
#    mat = np.delete(mat, (0), axis=0)
#    mat = mat.T
    mat = res.covariance
    print len(mat)
#    print mat
    
    dmb = np.sqrt(np.dot(np.dot(vect.T, mat), vect))
    cov_mb_c = res.covariance[2,2] * dmb_dc + res.covariance[2,1] * dmb_dx1 + res.covariance[2,0] * dmb_dx0
    cov_mb_x1 = res.covariance[1,1] * dmb_dx1 + res.covariance[1,0] * dmb_dx0 + res.covariance[1,2] * dmb_dc
    cov_mb_x0 = res.covariance[0,0] * dmb_dx0 +  res.covariance[0,1] * dmb_dx1 + res.covariance[0,2] * dmb_dc
    
    

    return dmb, cov_mb_c, cov_mb_x1, cov_mb_x0
    
    
class fake_res():    
    #SDSS12855
    def __init__(self, dict=None):             
        self.parameters = [999, 999, 8.75348389959e-05, -1.91819718101, 0.040562656067]
        self.covariance = np.array([[3.04977916498e-06**2, -2.26064049098e-07, -7.96590000957e-08],
               [-2.26064049098e-07, 0.297270343743**2, 0.000629305100058],
               [-7.96590000957e-08, 0.000629305100058, 0.0348134181789**2]])      


#        

    
    
    
    
    
    