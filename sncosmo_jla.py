#! /usr/bin/env python
import sncosmo
import os
from matplotlib import pyplot as plt
from read_jla import read_lc_jla
import numpy as np
import copy
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
import astropy.units as u

clight=299792.458
h=6.626070040*10**-34
t_min = -15
t_max = 45

t_min_sug = -12
t_max_sug = 42

# wavelength limits for salt2 model
wl_min_sal = 3000
wl_max_sal = 7000
    
def fit_salt2():
#    outfile = open('results/res_salt2.txt', 'w')
#    outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor mb_sncosmo x0 dx0')
    list_SDSS = ['lc-SDSS3377.list']
#    list_SDSS = [raw_input('Name of SN: ')]
    fitfail=[]        
#    for filename in os.listdir('./jla_data/jla_light_curves/'):
    for filename in list_SDSS:
        if 'lc-SDSS' == filename[:7]:
            sn_name = filename
            print sn_name
            head, data = read_lc_jla(sn_name, model = 'salt2')

            source = sncosmo.get_source('salt2', version='2.4')
            source.EBV_snfit = head['@MWEBV']    
            source.z_snfit = head['@Z_HELIO']
            source.Rv_snfit = 3.1
    
            dust = sncosmo.CCM89Dust()
            #model = sncosmo.Model(source=source)
            model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
            model.set(mwebv=head['@MWEBV'])
            model.set(z=head['@Z_HELIO'])
            data_original=copy.copy(data)
    
#            res, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], modelcov = True,bounds={'x1':(-10, 10),'c':(-5, 5)})
            try:     
                #initial iteration x1 fix
                model.set(x1=0.0)
                res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'c'], modelcov = False, guess_t0=True)
                chi2=res.chisq
                chi2p=chi2+1
                m=0
                while chi2 < chi2p:
                    t_peak = fitted_model.parameters[1]
                    x0_guess = fitted_model.parameters[2]
                    x1_guess = fitted_model.parameters[3]
                    c_guess = fitted_model.parameters[4]
#                    print t_peak,fitted_model.parameters[4]
        
                    t1 = t_peak + t_min*(1 + model.get('z'))
                    t2 = t_peak + t_max*(1 + model.get('z'))
                                

                    A=[]
                    data=copy.copy(data_original)
                    for i in range(len(data)):                    
                        if data[i][0] <= t1 or data[i][0] >= t2:
                            A.append(i)

                    A=np.array(A)    
                    for i in range(len(A)):
#                        print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data[A[i]][0],t1,t2)
                        data.remove_row(A[i])
                        A-=1

                    #model.set(t0=t_peak)
                    #er_t_peak = res.errors['t0']
        
#                    res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'x1', 'c'], modelcov = True, guess_t0=True)
                    res, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, bounds={'t0':(t_peak-5, t_peak+5),'x0':(x0_guess-5, x0_guess+5),'x1':(x1_guess-5, x1_guess+5),'c':(c_guess-5, c_guess+5)})       
                    chi2p=chi2
                    resp=res
                    fitted_modelp=fitted_model
                    chi2=res.chisq
                    print chi2p,chi2
                    m+=1

                #final results
                res=resp
                fitted_model=fitted_modelp



                
                #print fitted_model.get('t0')
                #print fitted_model.maxtime(), fitted_model.mintime()
#                print res
                print 'number of iteration ',m
#                sncosmo.plot_lc(data, model=fitted_model, errors=res.errors)
#                plt.savefig('SDSS9032_sugar.eps')
                plt.show()
#                mb=mB_determination(res)
#                print mb
                mb_sncosmo=fitted_model.bandmag('jla_STANDARD::B','jla_VEGA2', [res.parameters[1]])[0]
                print fitted_model.bandmag('STANDARD::B','VEGA2', [res.parameters[1]])
            except:
                fitfail.append(sn_name)
                print 'Error: fit fail for ',sn_name
            
            
            print(res.parameters)
            print fitfail

#    return res
#            outfile.write('%s %s \n' %(sn_name,fitted_model.bandmag('standard::b','jla1',fitted_model.parameters[1]))) 
#            outfile.write('%s \n' %(res)) 
        

#            outfile.write('\n')
        
#            outfile.write('%s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5],res.parameters[0],res.parameters[2],res.errors['x0'],res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c'],res.parameters[1],res.errors['t0'])) 
#            outfile.write('%s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5],res.parameters[0],mb,0,res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c'],res.parameters[1],res.errors['t0'])) 
            #outfile.write('%s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5],res.parameters[0],res.parameters[2],res.errors['x0'],res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c'],t_peak,er_t_peak))  
#            outfile.write('%s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999 %f %f %f' %(sn_name[3:-5],res.parameters[0],mb,0,res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c'],res.parameters[1],res.errors['t0'],mb_sncosmo,res.parameters[2],res.errors['x0'])) 
            
#    outfile.close()

def fit1():
    #par_c = {}

    outfile = open('res.txt', 'w')
    outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor')
    list_SDSS = ['lc-SDSS21062.list']
    #for filename in os.listdir('/home/maria/Dropbox/Science/Supernovae/JLA_SALT2/JLA_fit/jla_data'):

    for filename in list_SDSS:
        if 'lc-SDSS' == filename[:7] or 'lc-sn' == filename[:5]:
            sn_name = filename
            print sn_name
            source = sncosmo.get_source('salt2', version='2.4')
            source.EBV_abs = read_lcSDSS(sn_name)['@MWEBV']    
            source.z_abs = read_lcSDSS(sn_name)['@Z_HELIO']
            source.Rv_abs = 3.1
    
            dust = sncosmo.CCM89Dust()
            #model = sncosmo.Model(source=source)
            model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
            model.set(mwebv=read_lcSDSS(sn_name)['@MWEBV'])
            model.set(z=read_lcSDSS(sn_name)['@Z_HELIO'])
    
            data = sncosmo.read_lc('/home/maria/Dropbox/Science/Supernovae/sncosmo/test_jla/' + sn_name)
    
            # first iteration        
            model.set(x0=0.00000000001)
            res, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x1', 'c'], modelcov = False)

            t_peak = fitted_model.parameters[1]
            #x0_guess = fitted_model.parameters[2]
            #x1_guess = fitted_model.parameters[3]
            #c_guess = fitted_model.parameters[4]
            #print t_peak, x0_guess, x1_guess, c_guess
            t1 = t_peak + t_min*(1 + model.get('z'))
            t2 = t_peak + t_max*(1 + model.get('z'))

            data2 = copy.deepcopy(data)

            A=[]
            for i in range(len(data)):                    
                if data[i][0] <= t1 or data[i][0] >= t2:
                    A.append(i)
            A=np.array(A)    
            for i in range(len(A)):
                print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data2[A[i]][0],t1,t2)
                data2.remove_row(A[i])
                A-=1
            # second iteration
    
            res2, fitted_model2 = sncosmo.fit_lc(data2, model, ['t0', 'x0', 'x1', 'c'], modelcov = False)
    
            t_peak2 = fitted_model2.parameters[1]
    
            t12 = t_peak2 + t_min*(1 + model.get('z'))
            t22 = t_peak2 + t_max*(1 + model.get('z'))
    
            data3 = copy.deepcopy(data)
    
            A=[]
            for i in range(len(data)):                    
                if data[i][0] <= t12 or data[i][0] >= t22:
                    A.append(i)
            A=np.array(A)    
            for i in range(len(A)):
                print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data3[A[i]][0],t1,t2)
                data3.remove_row(A[i])
                A-=1
    
            # 3d iteration
            res3, fitted_model3 = sncosmo.fit_lc(data3, model, ['t0', 'x0', 'x1', 'c'], modelcov = True)
    
    
            #res, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, bounds={'t0':(t_peak-5, t_peak+5),'x0':(x0_guess-5, x0_guess+5),'x1':(x1_guess-5, x1_guess+5),'c':(c_guess-5, c_guess+5)})
    
            
            sncosmo.plot_lc(data3, model=fitted_model3, errors=res.errors)
            plt.show()    
            print res3
    
            #dic_param = {}
        #dic_param[sn_name]=
        
        #c = fitted_model.parameters[4]
        #c = res
        #par_c[filename] = c

            outfile.write('\n')

            #outfile.write('%s 999 %f 999 %f %f %f %f %f %f 999 999 999 999 999 999 999 999 999 999 999' %(sn_name[3:-5],res.parameters[0],res.parameters[2],res.errors['x0'],res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c']))  
    outfile.close()


    #color = np.genfromtxt('/home/maria/Dropbox/Science/Supernovae/sncosmo/test_color.dat')
    #color2 = np.genfromtxt('/home/maria/Dropbox/Science/Supernovae/sncosmo/test_color2.dat')
    #c_salt = color[:,1]
    #c_snc = color[:,2]
    #plt.plot((c_snc-c_salt)/c_salt*100,c_salt,'+', (color2[:,2]-color2[:,1])/color2[:,1]*100,color2[:,1],'o')
    #plt.show()


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
            model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
            
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
    snfit=np.loadtxt('./x0_mb_JLA.txt',dtype='str')
    name_snfit=np.array(snfit[:,0])
    mb_snfit=np.array(snfit[:,3],float)
    dmb_snfit=np.array(snfit[:,4],float)
    x0_snfit=np.array(snfit[:,5],float)
    x1_snfit=np.array(snfit[:,7],float)
    c_snfit=np.array(snfit[:,9],float)
    t0_snfit=np.array(snfit[:,11],float)
    z_snfit=np.array(snfit[:,13],float)
    for i in range(len(name_snfit)):
        sn_name_snfit='lc-'+name_snfit[i]+'.list'
        if sn_name_snfit==sn_name:
            x0_snfit=x0_snfit[i]
            x1_snfit=x1_snfit[i]
            c_snfit=c_snfit[i]
            t0_snfit=t0_snfit[i]
            mb_snfit=mb_snfit[i]
            dmb_snfit=dmb_snfit[i]
            z_snfit=z_snfit[i]
            
    model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
    model.set(mwebv=head['@MWEBV'], z=head['@Z_HELIO'], t0=t0_snfit, x0=x0_snfit ,x1=x1_snfit,c=c_snfit) #

    print model.parameters
    print sncosmo.chisq(data,model,modelcov=True)
    return model, mb_snfit,dmb_snfit,z_snfit
#    print model.bandflux('jla_SDSS::g', [54277.5, 54285.0, 54292.5, 54300., 54307.5, 54315., 54322.5], zp=25, zpsys='jla_AB_B12')
    #for i in range(13):
    #    print model.bandflux(data[i][1], float(data[i][0]), zp=float(data[i][4]), zpsys='AB_jla')
#    sncosmo.plot_lc(data, model=model)
#    plt.show()

def mB_determination(res):
    
    
    #interpolation of TB and Trest
    filt2=np.genfromtxt('./jla_data/Instruments/SNLS3-Landolt-model/sb-shifted.dat')
    wlen=filt2[:,0]
    tran=filt2[:,1]
    splB = Spline1d(wlen, tran, k=1,ext = 1)



    #interpolation of ref spectrum
    data = np.genfromtxt('./jla_data/MagSys/bd_17d4708_stisnic_003.ascii')
    dispersion = data[:,0]
    flux_density = data[:,1]
    splref = Spline1d(dispersion, flux_density, k=1,ext = 1)

  
    #interpolation of the spectrum model
    template_0=np.genfromtxt('./salt2-4/salt2_template_0.dat')    
    template_1=np.genfromtxt('./salt2-4/salt2_template_1.dat')
#    salt2source=sncosmo.SALT2Source('/users/divers/lsst/mondon/hubblefit/sncosmo_jla/salt2-4')
    
    wlM0=[]
    M0=[]
    for i in range(len(template_0[:,0])):
        if template_0[:,0][i]==0.0:
           wlM0.append(template_0[:,1][i]) 
           M0.append(template_0[:,2][i])
    splM0 = Spline1d(wlM0, M0, k=1,ext = 1)

    wlM1=[]
    M1=[]
    for i in range(len(template_1[:,0])):
        if template_1[:,0][i]==0.0:
           wlM1.append(template_1[:,1][i]) 
           M1.append(template_1[:,2][i])
    splM1 = Spline1d(wlM1, M1, k=1,ext = 1)

    #computation of the integral
    dt = 100000
    xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
    dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))
#    cl=[]
#    for i in xs:
#        cl.append(color_law_salt2(xs[i]))
#    I1=np.sum((splM0(xs)*10**-12+res.parameters[3]*splM1(xs)*10**-12)*(10**(-0.4*cl*res.parameters[4]))*xs*splB(xs)*dxs)

    I1=np.sum((splM0(xs)*10**-12+res.parameters[3]*splM1(xs)*10**-12)*(10**(-0.4*color_law_salt2(xs)*res.parameters[4]))*xs*splB(xs)*dxs)    
    I2=np.sum(splref(xs)*xs*splB(xs)*dxs)
    print I1, I2   
    
    
    
    #computation of mb
    mref=9.907
    mb=-2.5*np.log10(res.parameters[2]*(I1/I2))+mref

    return mb

def color_law_salt2(wl):
        B_wl=4302.57
        V_wl=5428.55
        l=(wl-B_wl)/(V_wl-B_wl)
        l_lo=(2800.-B_wl)/(V_wl-B_wl)
        l_hi=(7000.-B_wl)/(V_wl-B_wl)
        a=-0.504294
        b=0.787691
        c=-0.461715
        d=0.0815619
        cst= 1-(a+b+c+d)
        cl=[]
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
    
    mb_sncosmo=np.array(fit_salt2())
    mb_snfittxt=np.loadtxt('/users/divers/lsst/mondon/hubblefit/sncosmo_jla/x0_mb_JLA.txt',dtype='str')
    mb_snfit= np.array(mb_snfittxt[:,3],float)
    dmb_snfit=np.array(mb_snfittxt[:,4],float)
    plt.errorbar(mb_sncosmo,mb_snfit,xerr=np.nan,yerr=dmb_snfit,marker = 'o',color='red')
    plt.xlabel('sncosmo')
    plt.ylabel('snfit')
    plt.show()    
    
def mb_using_snfit_parameters():
    diff=[]
    z=[]
    for filename in os.listdir('./jla_data/jla_light_curves/'):
#    for filename in list_SDSS:
        if 'lc-SDSS' == filename[:7]:
            sn_name = filename
            res,mb_snfit1,dmb_snfit1,z1 =chi2(sn_name)
            mb_snfit_code=mB_determination(res)
            diff.append(mb_snfit_code-mb_snfit1)
            z.append(z1)
    plt.scatter(z,diff,marker = '.',color='red',label='mb')

            
            
          