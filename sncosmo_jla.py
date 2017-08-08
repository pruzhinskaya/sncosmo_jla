#! /usr/bin/env python
import sncosmo
import builtins_jla
import os
from matplotlib import pyplot as plt
from read_jla import read_lc_jla
import numpy as np
import copy

t_min = -15
t_max = 45

t_min_sug = -12
t_max_sug = 42

def fit_salt2():
	#outfile = open('res_salt2.txt', 'w')
	#outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor')
	list_SDSS = ['lc-sn1990af.list']
	#list_SDSS = [raw_input('Name of SN: ')]		
	#for filename in os.listdir('/home/maria/Dropbox/Science/Supernovae/JLA_SALT2/JLA_fit/jla_data'):
	for filename in list_SDSS:
		#if 'lc-SDSS' == filename[:7]:
			sn_name = filename
			print sn_name
			head, data = read_lc_jla(sn_name, model = 'salt2')

			source = sncosmo.get_source('salt2', version='2.4')
			source.EBV_abs = head['@MWEBV']	
			source.z_abs = head['@Z_HELIO']
			source.Rv_abs = 3.1
	
			dust = sncosmo.CCM89Dust()
			#model = sncosmo.Model(source=source)
			model = sncosmo.Model(source=source,effects=[dust],effect_names=['mw'],effect_frames=['obs'])
			model.set(mwebv=head['@MWEBV'])
			model.set(z=head['@Z_HELIO'])
	
			#res, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], modelcov = True,bounds={'x1':(-10, 10),'c':(-5, 5)})

			res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'x1', 'c'], modelcov = False, guess_t0=True)

			t_peak = fitted_model.parameters[1]
			#x0_guess = fitted_model.parameters[2]
			#x1_guess = fitted_model.parameters[3]
			#c_guess = fitted_model.parameters[4]
			#print t_peak,fitted_model.parameters[4]

			t1 = t_peak + t_min*(1 + model.get('z'))
			t2 = t_peak + t_max*(1 + model.get('z'))

			A=[]
			for i in range(len(data)):                    
				if data[i][0] <= t1 or data[i][0] >= t2:
					A.append(i)
			A=np.array(A)	
			for i in range(len(A)):
				#print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data[A[i]][0],t1,t2)
				data.remove_row(A[i])
				A-=1

			#model.set(t0=t_peak)
			#er_t_peak = res.errors['t0']

			res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'x1', 'c'], modelcov = True, guess_t0=True)
			#res, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, bounds={'t0':(t_peak-5, t_peak+5),'x0':(x0_guess-5, x0_guess+5),'x1':(x1_guess-5, x1_guess+5),'c':(c_guess-5, c_guess+5)})

			#print fitted_model.get('t0')
			#print fitted_model.maxtime(), fitted_model.mintime()
			print res
			sncosmo.plot_lc(data, model=fitted_model, errors=res.errors)
			plt.savefig('SDSS9032_sugar.eps')
			plt.show()
	
			#outfile.write('%s %s \n' %(sn_name,fitted_model.bandmag('standard::b','jla1',fitted_model.parameters[1]))) 
			#outfile.write('%s \n' %(res)) 
		

			#outfile.write('\n')
		
			#outfile.write('%s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5],res.parameters[0],res.parameters[2],res.errors['x0'],res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c'],res.parameters[1],res.errors['t0'])) 
			#outfile.write('%s 999 %f 999 %f %f %f %f %f %f 999 999 %f %f 999 999 999 999 999 999 999' %(sn_name[3:-5],res.parameters[0],res.parameters[2],res.errors['x0'],res.parameters[3],res.errors['x1'],res.parameters[4],res.errors['c'],t_peak,er_t_peak))  
	#outfile.close()

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

def chi2(name=None):
	sn_name='lc-SDSS18604.list'
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
