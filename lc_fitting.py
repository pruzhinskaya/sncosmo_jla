#! /usr/bin/env python
# The program for light curve fitting with Sugar model
import matplotlib.pyplot as plt
import numpy as np
import cPickle
from astropy import (cosmology, units as u, constants as const)
from scipy.optimize import curve_fit
from scipy import optimize

#from . import registry
#def get_bandpass(name):
#    """Get a Bandpass from the registry by name."""
#    return registry.retrieve(Bandpass, name)
#band = raw_input("type name of band " ) + '.dat'
#SNF=cPickle.load(open('UBVRI_BEDELLv1.pkl'))

HC_ERG_AA = const.h.cgs.value * const.c.to(u.AA / u.s).value
F_ab = 3631.*10**-23
h = 6.626070040*10**-27


gband = np.genfromtxt('/home/maria/Python/sncosmo/sdss/sdss_g.dat') # Download g_sdss filter
rband = np.genfromtxt('/home/maria/Python/sncosmo/sdss/sdss_r.dat') # Download r_sdss filter
iband = np.genfromtxt('/home/maria/Python/sncosmo/sdss/sdss_i.dat') # Download i_sdss filter
zband = np.genfromtxt('/home/maria/Python/sncosmo/sdss/sdss_z.dat') # Download z_sdss filter
Bband = np.genfromtxt('/home/maria/Python/sncosmo/bessell/bessell_b.dat') # Download bessell B filter
Vband = np.genfromtxt('/home/maria/Python/sncosmo/bessell/bessell_v.dat') # Download bessell V filter

d = {'rband': rband, 'iband': iband, 'Bband': Bband, 'Vband': Vband}

sn_r = np.genfromtxt('r.dat', usecols=[0,1]) # Download the observational data r
sn_i = np.genfromtxt('i.dat', usecols=[0,1]) # Download the observational data i
sn_B = np.genfromtxt('B.dat', usecols=[0,1]) # Download the observational data B
sn_V = np.genfromtxt('V.dat', usecols=[0,1]) # Download the observational data V

sn_data = {'rband': sn_r, 'iband': sn_i, 'Bband': sn_B, 'Vband': sn_V}

sugarfile = np.genfromtxt('/home/maria/Python/sncosmo/SUGAR_model.asci')
t = []
lamb = None
lamb_table = []
for line in sugarfile:
    if line[1] != lamb:
        if lamb != None:
            t.append(lamb_table)
        lamb_table = []
        lamb = line[1]
    lamb_table.append(line)

t = np.array(t)

def f(sn_data,q1,q2,q3,Av,Mgr):
	m_ab = []
	M_new = []
	m_last = []
	ans = {}
	#ans = dict( (bandname,[]) for bandname in d.keys())
	d_m_ab = dict( (bandname,[]) for bandname in d.keys())
	for i in range(0,len(t[0])):
		F = 10**(-0.4*(t[:,i,2] + q1*t[:,i,3] + q2*t[:,i,4] + q3*t[:,i,5] + Av*t[:,i,6] + Mgr*t[:,i,7]))
		x = t[:,1,1]
		phase = t[0,i,0]+54279
		for bandname,table in d.items():
			y_band = np.interp(x, table[:,0], table[:,1], 0, 0)
			f = np.trapz(y_band*F*x,x)/HC_ERG_AA
			fband = np.trapz(y_band*F_ab/x/h,x)
			value = -2.5* np.log10(f/fband)
			d_m_ab[bandname].append( [phase, value] )
	#return d_m_ab['iband']
	for bandname, datatable in d_m_ab.items():
		sndatum = sn_data[bandname]
		d_m_ab_dat = d_m_ab[bandname]
		dd = np.array(d_m_ab_dat)
		m_ab_new = np.interp(sndatum[:,0], dd[:,0], dd[:,1], 0, 0)
		ans[bandname] = m_ab_new
		#ans[bandname].append( [m_ab_new] )
	return ans

def chi2(q1,q2,q3,Av,Mgr):
	s_new = {}
	model = f(sn_data,q1,q2,q3,Av,Mgr)
	for bandname, datatable in sn_data.items():
		s = np.sum((-2.5*np.log10(sn_data[bandname][:,1]) - model[bandname])**2)
		s_new[bandname] = s
	chi = np.sum(s_new.values())
	return chi

#optimize.minimize( chi2, [1,1,1,1,1], method='Nelder-Mead')
myfit = optimize.fmin( lambda x: chi2(*x), [5,-7,1,1,40] )
print 'myfit parameters =',myfit

e=f(sn_data,myfit[0],myfit[1],myfit[2],myfit[3],myfit[4])['rband']
ee=f(sn_data,myfit[0],myfit[1],myfit[2],myfit[3],myfit[4])['iband']
e3=f(sn_data,myfit[0],myfit[1],myfit[2],myfit[3],myfit[4])['Bband']
e4=f(sn_data,myfit[0],myfit[1],myfit[2],myfit[3],myfit[4])['Vband']


plt.plot(sn_data['rband'][:,0],e-1,'r',sn_data['rband'][:,0],-2.5*np.log10(sn_data['rband'][:,1])-1,'ro')
plt.plot(sn_data['iband'][:,0],ee,'m',sn_data['iband'][:,0],-2.5*np.log10(sn_data['iband'][:,1]),'mo')
plt.plot(sn_data['Bband'][:,0],e3,'b',sn_data['Bband'][:,0],-2.5*np.log10(sn_data['Bband'][:,1]),'bo')
plt.plot(sn_data['Vband'][:,0],e4-2,'g',sn_data['Vband'][:,0],-2.5*np.log10(sn_data['Vband'][:,1])-2,'go')

plt.gca().invert_yaxis()
plt.title('SN fit')
plt.ylabel('magnitude')
plt.xlabel('T, days')
plt.show()

   # if sigma == None  or  (sigma == 0).any():
    #    sigma = np.ones_like(y)
    #return (( (f(x,q1,q2,q3,Av,Mgr) - y) / sigma )**2).sum() / (y.shape[0]-2)



			#m_ab.append(-2.5* np.log10(f/fband))
	#m_ab_new = np.interp(ph, t[1,:,0], m_ab, 0, 0)
		#print m_ab
		#M_new.append(m_ab)
		#M_new = np.array(M_new)
		#m_ab_new = np.interp(ph, t[1,:,0], M_new, 0, 0)
		#l1 = [item[0] for item in M_new]
		#l2 = [item[1] for item in M_new]
		#cc = [l1,l2]
	#return M_new
	#M_new = np.array(M_new)

	#for j in range(0,len(sn_data)-1):
	#	m_last = np.interp(sn_data[j][:,0]-54292, t[1,:,0], cc[j], 0, 0)
	#	m_last = np.append(m_last,m_last)
	#return m_last




	##l1 = [item[0] for item in M_new]
	##l2 = [item[1] for item in M_new]
	##m1 = np.interp(ph, t[1,:,0], l1, 0, 0)
	#for j in range(0,len(M_new[0])):
	#	ff = [item[j] for item in M_new]
	#model = {'rmodel':  [item[0] for item in M_new], 'imodel':  [item[1] for item in M_new]}

	##for name,table in zip(sn_data.items(),model.items():
	##	m_new = np.interp(ph, t[1,:,0], M_new, 0, 0)
	##return m_new



#def chi2(f, x, y, sigma=None):
#    if sigma == None  or  (sigma == 0).any():
#        sigma = np.ones_like(y)
#    return (( (f(x,q1,q2,q3,Av,Mgr) - y) / sigma )**2).sum() / (y.shape[0]-2)

#popt, pcov = curve_fit(f, t_i,sn_i)
#plt.plot(t_i, magapp_fit(t_i,92.9,-64,-74.6,-32.4,57),t_i,sn_i,'ro')
