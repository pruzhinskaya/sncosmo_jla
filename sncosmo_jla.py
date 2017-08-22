#! /usr/bin/env python
import sncosmo
import builtins_jla
import os
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from read_jla import read_lc_jla
import numpy as np
import copy

t_min = -15
t_max = 45

t_min_sug = -12
t_max_sug = 42

def fit_salt2():
	outfile = open('res_salt2.txt', 'w')
	outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor')
	list_jla = ['lc-SDSS17240.list','lc-sn2006az.list']
	#list_jla = os.listdir('jla_data/jla_light_curves/')
	fitfail = []
	for filename in list_jla:
		#if 'lc-SDSS' == filename[:7]:
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

			try:     
                		#first iteration (x1 is fixed)
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

					t1 = t_peak + t_min*(1 + model.get('z'))
					t2 = t_peak + t_max*(1 + model.get('z'))

					A=[]
					data_new = copy.copy(data)
					for i in range(len(data_new)):                    
						if data_new[i][0] <= t1 or data_new[i][0] >= t2:
							A.append(i)
					A=np.array(A)	
					for i in range(len(A)):
						#print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2)
						data_new.remove_row(A[i])
						A-=1
	
                    			res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, bounds={'t0':(t_peak-5, t_peak+5),'x0':(x0_guess-5, x0_guess+5),'x1':(x1_guess-5, x1_guess+5),'c':(c_guess-5, c_guess+5)})       
                    			chi2p=chi2
                    			resp=res
                    			fitted_modelp=fitted_model
                    			chi2=res.chisq
                    			print chi2p,chi2
                    			m+=1 

                		#final results
                		res=resp
				fitted_model=fitted_modelp

                		print 'number of iterations: ',m
				print res
				sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
				#plt.savefig('lc.eps')
				plt.show()
			except:
				fitfail.append(sn_name)
				print 'Error: fit fail for: ',sn_name
			print fitfail	

			outfile.write('\n %s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5],res.parameters[0],res.parameters[2],res.errors['x0'],res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c'],res.parameters[1],res.errors['t0'])) 
	outfile.close()


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

	#color = np.genfromtxt('/home/maria/Dropbox/Science/Supernovae/sncosmo/test_color.dat')
	#color2 = np.genfromtxt('/home/maria/Dropbox/Science/Supernovae/sncosmo/test_color2.dat')
	#c_salt = color[:,1]
	#c_snc = color[:,2]
	#plt.plot((c_snc-c_salt)/c_salt*100,c_salt,'+', (color2[:,2]-color2[:,1])/color2[:,1]*100,color2[:,1],'o')

def chi2(name=None):
	sn_name='lc-SDSS17240.list'
	head, data = read_lc_jla(sn_name)

	source = sncosmo.get_source('salt2', version='2.4')
	source.EBV_abs = head['@MWEBV']	
	source.z_abs = head['@Z_HELIO']
	source.Rv_abs = 3.1
	dust = sncosmo.CCM89Dust()

	model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
	model.set(mwebv=head['@MWEBV'], z=head['@Z_HELIO'], t0=54300, x0=7.0e-05,x1=-1,c=-0.04) #

	print model.parameters
	print sncosmo.chisq(data,model,modelcov=True)
	print model.bandflux('jla_SDSS::g', [54277.5, 54285.0, 54292.5, 54300., 54307.5, 54315., 54322.5], zp=25, zpsys='jla_AB_B12')
	#for i in range(13):
	#	print model.bandflux(data[i][1], float(data[i][0]), zp=float(data[i][4]), zpsys='AB_jla')
	sncosmo.plot_lc(data, model=model)
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
		if abs(dif) > np.sqrt(results('res_salt2.txt')[key][er_par]**2. + results('jla_data/jla_lcparams.txt')[key][er_par]**2.):
			print key
		plt.errorbar(results('res_salt2.txt')[key][par],results('jla_data/jla_lcparams.txt')[key][par],xerr=results('res_salt2.txt')[key][er_par],yerr=results('jla_data/jla_lcparams.txt')[key][er_par],marker = 'o',color='red')

		x = x + 1
	#plt.plot(range(x),np.zeros(x),'k')

	ax = plt.subplot(111)
	ax.xaxis.set_minor_locator(AutoMinorLocator(2))
	ax.yaxis.set_minor_locator(AutoMinorLocator(2))
	ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls = '--', color = 'black', alpha = 0.8) 
	plt.ylabel('snfit',fontsize=25)
	plt.xlabel('sncosmo',fontsize=25)
	plt.figtext(0.2, 0.8, par, fontsize=25)
	#plt.savefig('plot.png')
	plt.show()

def comparison_hist(par='tmax', er_par='dtmax'):
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
		#	print key
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
	#plt.savefig('plot.png')	
	plt.show()
