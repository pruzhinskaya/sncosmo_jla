#! /usr/bin/env python
# The program for
import matplotlib.pyplot as plt
import numpy as np
import cPickle
from astropy import (cosmology, units as u, constants as const)
#from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d

HC_ERG_AA = const.h.cgs.value * const.c.to(u.AA / u.s).value
F_ab = 3631.*10**-23 # erg/s/cm2/Hz
h = 6.626070040*10**-27

step = 1000
U = np.linspace(3360,4048.2,step)
B = np.linspace(4048.2,4877.3,step)
V = np.linspace(4877.3,5876.3,step)
R = np.linspace(5876.3,7079.9,step)
I = np.linspace(7078.9,8530.0,step)

dic_bands = {'U_sugar': np.array(zip(U,np.ones(len(U)))),'B_sugar': np.array(zip(B,np.ones(len(B)))),'V_sugar': np.array(zip(V,np.ones(len(V)))),'R_sugar': np.array(zip(R,np.ones(len(R)))),'I_sugar': np.array(zip(I,np.ones(len(I))))}

for key in dic_bands:
	outfile = open('PF_bands/'+key+'.txt', 'w')
	data = dic_bands[key]
	for line in data:
		trans = line[1]
		wave = line[0]
		outfile.write('%10.10f %10.10f' %(wave,trans))
		outfile.write('\n')
	outfile.close()


dic_spectra=cPickle.load(open('SNF_SN/all_CABALLO_data_binning_speed_for_maria.pkl'))

def light_curves(band_name,sn_name):
	band = dic_bands[band_name]
	wl_b = band[:,0]
        tran_b = band[:,1]
	M_ab = []

	sn = dic_spectra[sn_name]
	for key in sn.keys():
		wl_sn = sn[key]['X']
		flux_sn = sn[key]['Y_flux_without_cosmology']
		phase = sn[key]['days']
		error_sn = np.sqrt(sn[key]['V_flux'])

		for i in range(len(flux_sn)):
			if flux_sn[i] != flux_sn[i]:
				flux_sn[i] = 0

		y_band = np.interp(wl_sn, wl_b, tran_b, 0, 0)
		f = np.trapz(y_band*flux_sn*wl_sn,wl_sn)/HC_ERG_AA
		fband = np.trapz(y_band*F_ab/wl_sn/h,wl_sn)
		#value = -2.5* np.log10(f/fband)
		#value = f/fband/100000.
		value = f/fband

		er_flux = (np.trapz(y_band*error_sn*wl_sn,wl_sn)/HC_ERG_AA)/fband
		#er = []
		#for e1,e2 in zip(y_band,error_sn):
		#	if e1 != 0:
		#		er.append(error_sn)
		#er_flux = np.sum(er)/fband
		#print f,fband,value,er_flux
		M_ab.append([phase,value,er_flux])
	return M_ab

for name in dic_spectra:
	outfile = open('SNF_SN/'+name, 'w')
	outfile.write('time        band        flux          fluxerr      zp  zpsys\n')
	for key in dic_bands:
		lc = light_curves(key,name)
		for line in lc:
			time = line[0]
			flux = line[1]
			error = line[2]
			outfile.write('%10.10f %s %10.10f %10.10f 0 ab' %(time,key,flux,error))
			outfile.write('\n')
	outfile.close()
