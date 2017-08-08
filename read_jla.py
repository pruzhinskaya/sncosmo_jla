#! /usr/bin/env python
# The program to re-write the jla supernova data in sncosmo format
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import re
import sncosmo
import builtins_jla
from astropy.table import Table
#from ._extinction import ccm89


#path to jla data
p = './jla_data/'


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

def wl_cut_salt2(fname,EBV,z):
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

def wl_cut_sugar(fname,z):		
		filt = sncosmo.get_bandpass(fname)
		wlen = filt.wave
        	tran = filt.trans
		dt = 10000
		wlen_shift = wlen/(1+z)
		spl = Spline1d(wlen_shift, tran, k=1, ext = 1)
		#print
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

def read_lc_jla(sn_name):
<<<<<<< HEAD
	infile = open(p+'jla_light_curves/'+ sn_name, 'r') # download the sn data
=======
	infile = open('jla_data/jla_light_curves'+ sn_name, 'r') # download the sn data
>>>>>>> 812c849a1ded399866613d2a2cbd733e24aa1472
	numbers = []
	d = []
	words = []
	words_new = []
	head = {}
	for line in infile:
		if line[0] == '@': 
			d = line.split()
			if is_number(d[1]):
				head[d[0]] = float(d[1])
			else:
				head[d[0]] = d[1]
			continue
		elif line[0] == '#':  # miss the lines started from #
			continue	
		elif len(line) == 1:  # miss empty lines
			continue
		words = line.split()	
		words_new.append(float(words[0]))
		words_new.append(float(words[1]))
		words_new.append(float(words[2]))
		words_new.append(float(words[3]))
		words_new.append('jla_' + words[4])
		words_new.append('jla_' + words[5])
		numbers.append(words_new)
		words_new = []
	
#	numbers.append(words)
	infile.close()

	dic = {}
	for x in numbers:
		if x[4] in dic.keys():
			dic[x[4]]+=1
		else:
			dic[x[4]] = 1
      
	f_in = {}	
	f_out = {}
	for i in dic.keys():
		res = wl_cut_sugar(i,head['@Z_HELIO']) # sugar
		#res = wl_cut_salt2(i,head['@MWEBV'],head['@Z_HELIO']) # salt2
		if res[0] == 'True':
			f_in[i] = res[1]
		else:
			f_out[i] = res[1]
<<<<<<< HEAD
#			print 'We excluded passband %s (%d points) because it does not belong to the interval [%d,%d]' % ((i, dic[i],wl_min_sug,wl_max_sug) # sugar
			print 'We excluded passband %s (%d points) because restframewavelength = %7.3f does not belong to the interval [%d,%d]' % (i, dic[i],res[1],wl_min_sal,wl_max_sal) # salt
                 
=======
			print 'We excluded passband %s (%d points) because it does not belong to the interval [%d,%d]' % (i,len(dic[i]),wl_min_sug,wl_max_sug) # sugar
			#print 'We excluded passband %s (%d points) because restframewavelength = %7.3f does not belong to the interval [%d,%d]' % (i,len(dic[i]),res[1],wl_min_sal,wl_max_sal) # salt2
>>>>>>> 812c849a1ded399866613d2a2cbd733e24aa1472
	time = []
	band = []
	flux = []
	fluxerr = []
	zp = []
	zpsys = []
	for x in numbers:
		if x[4] in f_in.keys():
			time.append(x[0])
			band.append(x[4])
			flux.append(x[1])
			fluxerr.append(x[2])
			zp.append(x[3])
			zpsys.append(x[5])
		else:
			continue
	data=Table([time,band,flux,fluxerr,zp,zpsys], names=('time','band','flux','fluxerr','zp','zpsys'),meta={'name':'data'})
        return head, data	
		
def results(filename):
	salt2_param = open('/home/maria/Dropbox/Science/Supernovae/'+filename, 'r')
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

def comparison_plot(par='tmax',er_par = 'dtmax'):
	x = 0
	for key in results('sncosmo/res_com'):
		#dif = results('sncosmo/res_SDSS.txt')[key][par]-results('JLA_SALT2/jla_lcparams.txt')[key][par]
		#if abs(dif) > 2.3:
		#if results('sncosmo/res_com')[key][er_par] > 1:
		#	print key
		#plt.errorbar(results('sncosmo/res_SDSS.txt')[key][par],results('JLA_SALT2/jla_lcparams.txt')[key][par],xerr=results('sncosmo/res_SDSS.txt')[key][er_par],yerr=results('JLA_SALT2/jla_lcparams.txt')[key][er_par],marker = 'o',color='red')
		plt.errorbar(x,results('sncosmo/res_com')[key][par]-results('JLA_SALT2/jla_lcparams.txt')[key][par],results('JLA_SALT2/jla_lcparams.txt')[key][er_par],marker = 'o',color='red')
		#plt.plot(results('sncosmo/ALLSDSSres.txt')[key]['color'],results('JLA_SALT2/jla_lcparams.txt')[key]['color'],'ro')
		#plt.plot(x,results('sncosmo/ALLSDSSres.txt')[key][par]-results('JLA_SALT2/jla_lcparams.txt')[key][par],'ro')
		x = x + 1
	#plt.plot(range(x),np.zeros(x),'k')
	#plt.plot([53540,54650],[53540,54650],'k')
	#plt.xlim([53540,54650])
	#plt.ylim([53540,54650])
	#plt.xlabel('sncosmo')
	#plt.ylabel('snfit')
	plt.ylabel('$t_{snfit}-t_{sncosmo}$',fontsize=25)
	plt.xlabel('N')
	plt.savefig('all_tmax_n.png')
	plt.show()

def comparison_hist(par='tmax',er_par = 'dtmax'):

	gs = gridspec.GridSpec(2, 1) #subplots ratio
	f, (ax1, ax2) = plt.subplots(2, sharex=True)

	ax1 = plt.subplot(gs[0, 0])
	ax2 = plt.subplot(gs[1, 0])

	dif_par = []
	dif_er_par = []	
	for key in results('sncosmo/NBdata'):
		dif_par.append(results('sncosmo/NBdata')[key][par]-results('JLA_SALT2/jla_lcparams.txt')[key][par])
		dif_er_par.append(results('sncosmo/NBdata')[key][er_par]-results('JLA_SALT2/jla_lcparams.txt')[key][er_par])
	ax1.hist(dif_par,25,label=par)
	ax2.hist(dif_er_par,25,label=er_par)
	ax1.legend()
	ax2.legend()
	ax1.set_ylabel('N')
	ax2.set_ylabel('N')
	plt.savefig('nb_tmax.png')	
	plt.show()
