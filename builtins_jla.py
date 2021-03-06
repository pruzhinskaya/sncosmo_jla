# for JLA filters, mag.systems and SUGAR model
"""Importing this module registers loaders for built-in data structures:

- Bandpasses
- MagSystems
- Sources
"""
from __future__ import print_function
import numpy as np
from scipy.interpolate import RectBivariateSpline as Spline2d
import sncosmo
import astropy.units as u
from sncosmo.models import _SOURCES

p = 'jla_data/'
# =============================================================================
# Bandpasses

jla_meta = {
    'filterset': 'jla',
    'reference': ('Betoule14', 'http://supernovae.in2p3.fr/salt/doku.php?id=instruments'),
    'description': 'Representation of instruments and magnitude systems used in JLA'}

bands = [('jla_SDSS::u', 'SDSS/u.dat', jla_meta),
         ('jla_SDSS::g', 'SDSS/g.dat', jla_meta),
         ('jla_SDSS::r', 'SDSS/r.dat', jla_meta),
         ('jla_SDSS::i', 'SDSS/i.dat', jla_meta),
         ('jla_SDSS::z', 'SDSS/z.dat', jla_meta),
	 ('jla_4SHOOTER2::Us', '4Shooter2/Us_4Shooter2.txt', jla_meta),
	 ('jla_4SHOOTER2::B', '4Shooter2/B_4Shooter2.txt', jla_meta),
	 ('jla_4SHOOTER2::V', '4Shooter2/V_4Shooter2.txt', jla_meta),
	 ('jla_4SHOOTER2::R', '4Shooter2/R_4Shooter2.txt', jla_meta),
	 ('jla_4SHOOTER2::I', '4Shooter2/I_4Shooter2.txt', jla_meta),
	 ('jla_KEPLERCAM::Us', 'Keplercam/Us_Keplercam.txt', jla_meta),
	 ('jla_KEPLERCAM::B', 'Keplercam/B_Keplercam.txt', jla_meta),
	 ('jla_KEPLERCAM::V', 'Keplercam/V_Keplercam.txt', jla_meta),
	 ('jla_KEPLERCAM::i', 'Keplercam/i_Keplercam.txt', jla_meta),
	 ('jla_KEPLERCAM::r', 'Keplercam/r_Keplercam.txt', jla_meta),
	 ('jla_STANDARD::U', 'SNLS3-Landolt-model/sux-shifted.dat', jla_meta),
	 ('jla_STANDARD::B', 'SNLS3-Landolt-model/sb-shifted.dat', jla_meta),
	 ('jla_STANDARD::V', 'SNLS3-Landolt-model/sv-shifted.dat', jla_meta),
	 ('jla_STANDARD::R', 'SNLS3-Landolt-model/sr-shifted.dat', jla_meta),
	 ('jla_STANDARD::I', 'SNLS3-Landolt-model/si-shifted.dat', jla_meta),
	 ('jla_ACSWF::F606W', 'ACSWF/F606W.dat', jla_meta),
	 ('jla_ACSWF::F625W', 'ACSWF/F625W.dat', jla_meta),
	 ('jla_ACSWF::F775W', 'ACSWF/F775W.dat', jla_meta),
	 ('jla_ACSWF::F814W', 'ACSWF/F814W.dat', jla_meta),
	 ('jla_ACSWF::F850LP', 'ACSWF/F850LP.dat', jla_meta),
	 ('jla_NICMOS2::F110W', 'NICMOS2/F110W.dat', jla_meta),
	 ('jla_NICMOS2::F160W', 'NICMOS2/F160W.dat', jla_meta),
	 ('jla_MEGACAMPSF::u', 'Megacam-PSF/u0.list', jla_meta),
	 ('jla_MEGACAMPSF::g', 'Megacam-PSF/g0.list', jla_meta),
	 ('jla_MEGACAMPSF::r', 'Megacam-PSF/r0.list', jla_meta),
	 ('jla_MEGACAMPSF::i', 'Megacam-PSF/i0.list', jla_meta),
	 ('jla_MEGACAMPSF::z', 'Megacam-PSF/z0.list', jla_meta),
	 ('jla_SWOPE2::u', 'Swope2/u_texas_WLcorr_atm.txt', jla_meta),
	 ('jla_SWOPE2::B', 'Swope2/B_texas_WLcorr_atm.txt', jla_meta),
	 ('jla_SWOPE2::g', 'Swope2/g_texas_WLcorr_atm.txt', jla_meta),
	 ('jla_SWOPE2::V', 'Swope2/V_LC3014_texas_WLcorr_atm.txt', jla_meta),
	 ('jla_SWOPE2::V1', 'Swope2/V_LC3009_texas_WLcorr_atm.txt', jla_meta),
	 ('jla_SWOPE2::V2', 'Swope2/V_LC9844_texax_WLcorr_atm.txt', jla_meta),
	 ('jla_SWOPE2::r', 'Swope2/r_texas_WLcorr_atm.txt', jla_meta),
	 ('jla_SWOPE2::i', 'Swope2/i_texas_WLcorr_atm.txt', jla_meta)]

for name, fname, meta in bands:
	gen = (r.encode('utf-8') for r in open(p + 'Instruments/' +  fname) if not r[0] in ('@', '#'))
	data = np.genfromtxt(gen)
	band = sncosmo.Bandpass(data[:,0],data[:,1],wave_unit=u.AA,name=name)
	sncosmo.registry.register(band)

# =============================================================================
# MagSystems

def load_spectral_magsys_fits2(relpath, name=None, version=None):
    data = np.genfromtxt(relpath)
    dispersion = data[:,0]
    flux_density = data[:,1]
    refspectrum = sncosmo.spectrum.Spectrum(dispersion, flux_density,
                           unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    return sncosmo.magsystems.SpectralMagSystem(refspectrum, name=name)

website = 'http://supernovae.in2p3.fr/salt/doku.php?id=instruments'
subclass = '`~sncosmo.SpectralMagSystem`'

VEGAHST_desc = 'Vega0.dat'
AB_jla_desc = 'B12-AB-off.dat'
VEGA2_desc = 'BD17-jla1.dat'
VEGA2_mb_desc = 'BD17-snls3.dat'


for name, fn, desc in [('jla_VEGAHST', 'alpha_lyr_stis_005.ascii', VEGAHST_desc),
		       ('jla_AB_B12_0', 'ab-spec.dat', AB_jla_desc),
                       ('jla_VEGA2_0', 'bd_17d4708_stisnic_003.ascii', VEGA2_desc),
		       ('jla_VEGA2_mb_0', 'bd_17d4708_stisnic_002.ascii', VEGA2_mb_desc)]:
    sncosmo.registry.register_loader(sncosmo.MagSystem, name, load_spectral_magsys_fits2,
                             args=[p + 'MagSys/' + fn],
                             meta={'subclass': subclass, 'url': website,
                                   'description': desc})

# offsets are in the sense (mag_SDSS - mag_AB) = offset
# -> for example: a source with AB mag = 0. will have SDSS mag = 0.06791
bands_ab = {'jla_SDSS::u': ('jla_AB_B12_0',  0.06791),
        'jla_SDSS::g': ('jla_AB_B12_0', -0.02028),
        'jla_SDSS::r': ('jla_AB_B12_0', -0.00493),
        'jla_SDSS::i': ('jla_AB_B12_0', -0.01780),
        'jla_SDSS::z': ('jla_AB_B12_0', -0.01015),
        'jla_MEGACAMPSF::g': ('jla_AB_B12_0', 0),
        'jla_MEGACAMPSF::r': ('jla_AB_B12_0', 0),
        'jla_MEGACAMPSF::i': ('jla_AB_B12_0', 0),
        'jla_MEGACAMPSF::z': ('jla_AB_B12_0', 0)}

sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_ab),'jla_AB_B12')

bands_BD17 ={'jla_STANDARD::U': ('jla_VEGA2_0', 9.724),
	'jla_STANDARD::B': ('jla_VEGA2_0', 9.907),
	'jla_STANDARD::V': ('jla_VEGA2_0', 9.464),
	'jla_STANDARD::R': ('jla_VEGA2_0', 9.166),
	'jla_STANDARD::I': ('jla_VEGA2_0', 8.846),
	'jla_SDSS::u': ('jla_VEGA2_0', 10.56),
	'jla_SDSS::g': ('jla_VEGA2_0', 9.64),
	'jla_SDSS::r': ('jla_VEGA2_0', 9.35),
	'jla_SDSS::i': ('jla_VEGA2_0', 9.25),
	'jla_SDSS::z': ('jla_VEGA2_0', 9.23),
	'jla_KEPLERCAM::Us': ('jla_VEGA2_0', 9.724), # U standard (normalement)
	'jla_KEPLERCAM::B': ('jla_VEGA2_0', 9.8803),
	'jla_KEPLERCAM::V': ('jla_VEGA2_0', 9.4722),
	'jla_KEPLERCAM::r': ('jla_VEGA2_0', 9.3524),
	'jla_KEPLERCAM::i': ('jla_VEGA2_0', 9.2542),
	'jla_4SHOOTER2::Us': ('jla_VEGA2_0', 9.724), # U standard (normalement)
	'jla_4SHOOTER2::B': ('jla_VEGA2_0', 9.8744),
	'jla_4SHOOTER2::V': ('jla_VEGA2_0', 9.4789),
	'jla_4SHOOTER2::R': ('jla_VEGA2_0', 9.1554),
	'jla_4SHOOTER2::I': ('jla_VEGA2_0', 8.8506),
	'jla_SWOPE2::u': ('jla_VEGA2_0',   10.514),
	'jla_SWOPE2::g': ('jla_VEGA2_0',   9.64406),
	'jla_SWOPE2::r': ('jla_VEGA2_0',   9.3516),
	'jla_SWOPE2::i': ('jla_VEGA2_0',   9.25),
	'jla_SWOPE2::B': ('jla_VEGA2_0',   9.876433),
	'jla_SWOPE2::V': ('jla_VEGA2_0',   9.476626),
	'jla_SWOPE2::V1': ('jla_VEGA2_0',  9.471276),
	'jla_SWOPE2::V2': ('jla_VEGA2_0',  9.477482)}

sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_BD17),'jla_VEGA2')

bands_BD17_mb ={'jla_STANDARD::U': ('jla_VEGA2_mb_0', 9.724),
	'jla_STANDARD::B': ('jla_VEGA2_mb_0', 9.907),
	'jla_STANDARD::V': ('jla_VEGA2_mb_0', 9.464),
	'jla_STANDARD::R': ('jla_VEGA2_mb_0', 9.166),
	'jla_STANDARD::I': ('jla_VEGA2_mb_0', 8.846),
	'jla_SDSS::u': ('jla_VEGA2_mb_0', 10.56),
	'jla_SDSS::g': ('jla_VEGA2_mb_0', 9.64),
	'jla_SDSS::r': ('jla_VEGA2_mb_0', 9.35),
	'jla_SDSS::i': ('jla_VEGA2_mb_0', 9.25),
	'jla_SDSS::z': ('jla_VEGA2_mb_0', 9.23),
	'jla_KEPLERCAM::B': ('jla_VEGA2_mb_0', 9.8803),
	'jla_KEPLERCAM::V': ('jla_VEGA2_mb_0', 9.4722),
	'jla_KEPLERCAM::r': ('jla_VEGA2_mb_0', 9.3524),
	'jla_KEPLERCAM::i': ('jla_VEGA2_mb_0', 9.2542),
	'jla_4SHOOTER2::B': ('jla_VEGA2_mb_0', 9.8744),
	'jla_4SHOOTER2::V': ('jla_VEGA2_mb_0', 9.4789),
	'jla_4SHOOTER2::R': ('jla_VEGA2_mb_0', 9.1554),
	'jla_4SHOOTER2::I': ('jla_VEGA2_mb_0', 8.8506)}

sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_BD17_mb),'jla_VEGA2_mb')
# =============================================================================
# Sources
p2 = '/Users/Maria/Dropbox/Science/Supernovae/sncosmo/SUGAR/sugar_model/'

from operator import itemgetter
from sncosmo.salt2utils import BicubicInterpolator

infile = open(p2 + 'SUGAR_model_v1.asci', 'r')
lines = infile.readlines()
infile.close()

listik = []
for line in lines:
    new_line = [float(i) for i in line.split()]
    listik.append(new_line)

s = sorted(listik, key=itemgetter(0))

names = ['0','1','2','3','4']
for i in names:
    outfile = open(p2 + 'sugar_template_' + i + '.dat', 'w')
    for line in s:
        j = 2+int(i)
        outfile.write('%4.4f %8.8f %8.8f' %(line[0],line[1],line[j]))
        outfile.write('\n')
    outfile.close()


class SUGARSource(sncosmo.Source):

    _param_names = ['q1', 'q2', 'q3', 'A', 'Mgr']
    param_names_latex = ['q_1', 'q_2', 'q_3', 'A', 'M_g']
    _SCALE_FACTOR = 1

    def __init__(self, modeldir=None,
                 m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
		 m2file='sugar_template_2.dat',
		 m3file='sugar_template_3.dat',
		 m4file='sugar_template_4.dat',
                 name=None, version=None):
        self.name = name
        self.version = version
        self._model = {}
        self._parameters = np.array([1., 1., 1., 1., 40.])

        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            phase, wave, values = sncosmo.read_griddata_ascii(p2 + names_or_objs[key])
            values *= self._SCALE_FACTOR
            # self._model[key] = Spline2d(phase, wave, values, kx=2, ky=2)
            self._model[key] = BicubicInterpolator(phase, wave, values)

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                self._phase = phase
                self._wave = wave

    def _flux(self, phase, wave):
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        m2 = self._model['M2'](phase, wave)
        m3 = self._model['M3'](phase, wave)
        m4 = self._model['M4'](phase, wave)
        return (10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4] + 48.59)) / (wave ** 2 / 299792458. * 1.e-10))


from sncosmo.builtins import DATADIR

# Sugar model
def load_sugarmodel(relpath, name=None, version=None):
    abspath = DATADIR.abspath(relpath, isdir=True)
    return SUGARSource(modeldir=abspath, name=name, version=version)

website = 'http://no'
PF16ref = ('PF16', 'PF et al. 2016 '
          '<http://arxiv.org/>')
for topdir, ver, ref in [('sugar_model', '1.0', PF16ref)]:
    meta = {'type': 'SN Ia', 'subclass': '`~sncosmo.SUGARSource`',
            'url': website, 'reference': ref}
    _SOURCES.register_loader('sugar', load_sugarmodel,
                             args=('/Users/Maria/Dropbox/Science/Supernovae/sncosmo/SUGAR/'+topdir,),
                             version=ver, meta=meta)
