#!/usr/bin/env python
# coding: utf-8

# import necessary modules
import numpy as np
from classy import Class
from math import pi
import matplotlib
import matplotlib.pyplot as plt
from math import pi

#####################################################
#
# Cosmological parameters and other CLASS parameters
#
#####################################################
gomp_settings = {# LambdaCDM parameters
                   'h': 0.67810,
                   'omega_b': 0.02238280,
                   'omega_cdm': 0.1201075,
#                   'sigma8': 0.8159,
                   'ln_A_s_1e10': 3.04478383,
                   'n_s': 0.9660499,
                   'reio_parametrization': 'reio_gomp',
                   'tau_reio': 0.05430842
}

tanh_settings = {# LambdaCDM parameters
		   'h': 0.67810,
		   'omega_b': 0.02238280,
		   'omega_cdm': 0.1201075,
#		   'sigma8': 0.8159,
           'ln_A_s_1e10': 3.04478383,
		   'n_s': 0.9660499,
		   'reio_parametrization': 'reio_camb',
		   'tau_reio': 0.05430842
}

common = {
	'output':'tCl,pCl,lCl,mPk',
	'lensing':'yes',
	'P_k_max_1/Mpc':3.0}

# max l
l_max = 2500

# gomp classy run
gomp = Class()
gomp.set(gomp_settings)
gomp.set(common)
gomp.compute()

print('The optical depth for Gomp is ', gomp.tau_reio())

gomp_cls = gomp.lensed_cl(l_max)
ll = gomp_cls['ell'][2:] # same for both of them
gomp_clTT = gomp_cls['tt'][2:]
gomp_clTE = gomp_cls['te'][2:]
gomp_clEE = gomp_cls['ee'][2:]
gomp_clBB = gomp_cls['bb'][2:]

# tanh classy run
tanh = Class()
tanh.set(tanh_settings)
tanh.set(common)
tanh.compute()

print('The optical depth for Tanh is ', tanh.tau_reio())

tanh_cls = tanh.lensed_cl(l_max)
tanh_clTT = tanh_cls['tt'][2:]
tanh_clTE = tanh_cls['te'][2:]
tanh_clEE = tanh_cls['ee'][2:]
tanh_clBB = tanh_cls['bb'][2:]

# plotting time
font = {'size':16, 'family':'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)
matplotlib.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]

factor = 1.e10 * ll * (ll + 1.) / 2. / pi

plt.xlim([2,l_max])
plt.ylim([1.e-8,10])
plt.xlabel(r"$\ell$")
plt.ylabel(r"$\ell (\ell+1) C_{\ell}^{XY} / (2 \pi) \, \, \, [\times 10^{10}]$")
plt.grid()

plt.loglog(ll, factor * gomp_clTT, 'r-', label=r'Gomp: $\mathrm{TT}$')
plt.loglog(ll, factor * tanh_clTT, 'r:', label=r'Tanh: $\mathrm{TT}$')
plt.loglog(ll, factor * np.abs(gomp_clTE), '-', color='purple', label=r'Gomp: $|\mathrm{TE}|$')
plt.loglog(ll, factor * np.abs(tanh_clTE), ':', color='purple', label=r'Tanh: $|\mathrm{TE}|$')
plt.loglog(ll, factor * gomp_clEE, 'b-', label=r'Gomp: $\mathrm{EE}$')
plt.loglog(ll, factor * tanh_clEE, 'b:', label=r'Tanh: $\mathrm{EE}$')
plt.loglog(ll, factor * gomp_clBB, 'g-', label=r'Gomp: $\mathrm{BB}$')
plt.loglog(ll, factor * tanh_clBB, 'g:', label=r'Tanh: $\mathrm{BB}$')
plt.legend(loc='right',bbox_to_anchor=(1.4,0.5))
plt.savefig('cls_comp_tau.pdf',bbox_inches='tight')


print(gomp_cls.keys())
# working on this


