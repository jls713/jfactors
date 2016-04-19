# -*- coding: utf-8 -*-
### Generates mass, J and D profiles for a range of alpha,beta,gamma models
### Produces Fig. 7 of Evans, Sanders & Geringer-Sameth (2016)

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('/home/jls/work/data/jfactors/')
from spherical_Jfactors import wyns_formulaJ_NFW, wyns_formulaD_NFW
from matplotlib.ticker import MaxNLocator # added
from spherical_Jfactors import *

angs = np.deg2rad(np.logspace(np.log10(5e-3),np.log10(0.6),30))
Rhalf = 0.03 ## 30pc
sig = 3.     ## 3km/s
G = 4.300918e-6 ## in units solar mass, km/s kpc
Mhalf = 4.*sig**2*Rhalf/G
rs = 0.15 ## scale radius of NFW units kpc
D = 30.   ## distance kpc
gamma = [0.,0.2,0.4,0.6,0.8,1.,1.2]
beta = [3.,3.5,4.,4.5,5.,5.5,6.]
alpha = [1.,1.5,2.]
rt = 10.

angs_dimless = angs*D/Rhalf
max_M,min_M = np.zeros(len(angs)),np.ones(len(angs))*1e50
max_J,min_J = np.zeros(len(angs)),np.ones(len(angs))*1e50
max_D,min_D = np.zeros(len(angs)),np.ones(len(angs))*1e50

f,a=plt.subplots(3,1,figsize=[3.32,5.5])
plt.subplots_adjust(hspace=0.)
for b,c in zip(beta,sns.color_palette()):
	for g in gamma:
		for al in alpha:
			rho0 = 1.
			M = integrate_rho_spherical_alphabetagamma(Rhalf,rho0,rs,al,b,g,rt)
			rho0=Mhalf/M

			for n,x in enumerate(angs):
				mm = integrate_rho_spherical_alphabetagamma(x*D,rho0,rs,al,b,g,rt)
				if(mm>max_M[n]):
					max_M[n]=mm
				if(mm<min_M[n]):
					min_M[n]=mm
				jj = integrate_J_spherical_alphabetagamma(x,D,rho0,rs,al,b,g,rt)
				jj2 = integrate_J_farfield_spherical_alphabetagamma(x,D,rho0,rs,al,b,g,rt)
				print jj,jj2
				if(jj>max_J[n]):
					max_J[n]=jj
				if(jj<min_J[n]):
					min_J[n]=jj
				dd = integrate_D_farfield_spherical_alphabetagamma(x,D,rho0,rs,al,b,g,rt)
				if(dd>max_D[n]):
					max_D[n]=dd
				if(dd<min_D[n]):
					min_D[n]=dd

a[0].fill_between(angs_dimless,min_M,max_M,alpha=0.5,color=sns.color_palette()[0])
a[1].fill_between(angs_dimless,min_J,max_J,alpha=0.5,color=sns.color_palette()[0])
a[2].fill_between(angs_dimless,min_D,max_D,alpha=0.5,color=sns.color_palette()[0])
a[1].plot(angs_dimless,wyns_formulaJ_NFW(sig,Rhalf*1000.,D,np.rad2deg(angs),rs),color='k')
a[2].plot(angs_dimless,wyns_formulaD_NFW(sig,Rhalf*1000.,D,np.rad2deg(angs),rs),color='k')
a[0].semilogx()
a[1].semilogx()
a[2].semilogx()
a[0].semilogy()
a[0].set_xlim(0.1,10.)
a[1].set_xlim(0.1,10.)
a[2].set_xlim(0.1,10.)
# a[1].set_xlabel(r'$\alpha/^\circ$')
a[2].set_xlabel(r'$D\theta/R_h$')
a[0].set_xticklabels([])
a[1].set_xticklabels([])
a[1].yaxis.set_major_locator(MaxNLocator(prune='upper'))
a[2].yaxis.set_major_locator(MaxNLocator(prune='upper'))
a[0].annotate(r'$\sigma_{\mathrm{los}}=3\,\mathrm{km\,s}^{-1},\,R_{\mathrm{half}}=30\,\mathrm{pc}$',xy=(0.1,2e8),annotation_clip=False,fontsize=14)
l=a[0].axvline(1.,ls='dashed',color='k')
l.set_dashes((3,1))
l=a[1].axvline(2.,ls='dashed',color='k')
l.set_dashes((3,1))
l=a[1].axvline(1.,ls='dashed',color='k')
l.set_dashes((3,1))
a[1].annotate('Walker et al. (2011)', xy=(2.1,18.5),rotation=90.,annotation_clip=False)
l=a[2].axvline(2.,ls='dashed',color='k')
l.set_dashes((3,1))
l=a[2].axvline(1.,ls='dashed',color='k')
l.set_dashes((3,1))
l=a[1].axvline(8.72,ls='dashed',color='r')
l.set_dashes((3,1))
l=a[2].axvline(8.72,ls='dashed',color='r')
l.set_dashes((3,1))
a[2].annotate(r'$\theta=0.5^\circ$', xy=(7.,17.),rotation=90.,annotation_clip=False)
a[2].set_ylim(15.5,19.)
a[0].set_ylabel(r'$M(D\theta)/\mathrm{M}_\odot$')
a[1].set_ylabel(r'$\log_{10} [J(\theta)/\mathrm{GeV\,cm}^{-5}]$')
a[2].set_ylabel(r'$\log_{10} [D(\theta)/\mathrm{GeV\,cm}^{-2}]$')
plt.savefig('spherical_comparison.pdf',bbox_inches='tight')
