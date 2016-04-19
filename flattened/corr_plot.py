# -*- coding: utf-8 -*-
### Generates Fig 4 of SEG 2016
### First need to generate the table using compute_corrections.py
## ============================================================================
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
sys.path.append('../spherical/')
from J_D_table import posh_names
## ============================================================================
## 1. Load in the data
g = np.genfromtxt("../spherical/dwarfs_Jfactors_ascii.dat")
g2 = np.genfromtxt("dwarfs_Jfactors_corr_ascii.dat")[1:]
nam2 = pd.read_csv("dwarfs_Jfactors_corr_ascii.dat",sep=' ')
nam = pd.read_csv("../spherical/dwarfs_Jfactors_ascii.dat",sep=' ')
nam = nam.rename(columns={'#Name': 'Name'})

## ============================================================================
## 2. Sort data by ellipticity
temp = [0]*len(g)
for j in range(len(g2)):
	for i in range(len(g)):
		if(nam2.Name[j]==nam.Name[i]):
			temp[j]=i
g = g[temp]
nam = nam.Name.values[temp]
nam2 = nam2.Name.values
temp = (-g.T[7]).argsort()
g = g[temp]
g2 = g2[temp]
nam = nam[temp]

## ============================================================================
## 3. Just prolate corrections against J-factor
f=plt.figure(figsize=[3.32,2.5])
plt.errorbar(g[:8].T[7],g2[:8].T[2],yerr=None,xerr=[g[:8].T[8],g[:8].T[9]],fmt='d',color=sns.color_palette()[2],ms=4,label='Classical dwarfs')
plt.errorbar(g[8:].T[7],g2[8:].T[2],yerr=None,xerr=[g[8:].T[8],g[8:].T[9]],fmt='o',color='k',ms=3,label='Ultrafaints')
plt.legend(loc=2,frameon=False)
plt.ylabel(r'$\mathcal{F}$')
plt.xlabel(r'$\log_{10}(J(0.5^\circ)/\,\mathrm{GeV^2\,cm}^{-5})$')
plt.savefig('corr_plot_prol.pdf',bbox_inches='tight')
plt.clf()

## ============================================================================
## 4. J factor against distance showing prolate corrections
plt.errorbar(g[:8].T[1],g[:8].T[7],xerr=g[:8].T[2],yerr=[g[:8].T[8],g[:8].T[9]],fmt='d',color=sns.color_palette()[2],ms=4,label='Classical dwarfs')
plt.errorbar(g[8:].T[1],g[8:].T[7],xerr=g[8:].T[2],yerr=[g[8:].T[8],g[8:].T[9]],fmt='o',color='k',ms=3,label='Ultrafaints')
for i in range(len(g)):
	plt.arrow(g[i][1],g[i][7],0.,g2[i][2],color=sns.color_palette()[0],head_width=7.,head_length=0.15,lw=1)
plt.ylim(15.,22.)
plt.axhline(g[-6][7]+g2[-6][2],color='k',ls='dashed')
plt.legend(loc=1,frameon=False)
plt.xlabel(r'$\mathrm{Distance/\,kpc}$')
plt.ylabel(r'$\log_{10}(J(0.5^\circ)/\,\mathrm{GeV^2\,cm}^{-5})$')
plt.savefig('corr_plot_prol_JD.pdf',bbox_inches='tight')
plt.clf()

## ============================================================================
## 5. J factor against distance showing oblate corrections
plt.errorbar(g[:8].T[1],g[:8].T[7],xerr=g[:8].T[2],yerr=[g[:8].T[8],g[:8].T[9]],fmt='d',color=sns.color_palette()[2],ms=4,label='Classical dwarfs')
plt.errorbar(g[8:].T[1],g[8:].T[7],xerr=g[8:].T[2],yerr=[g[8:].T[8],g[8:].T[9]],fmt='o',color='k',ms=3,label='Ultrafaints')
for i in range(len(g)):
	plt.arrow(g[i][1],g[i][7],0.,g2[i][1],color=sns.color_palette()[0],head_width=7.,head_length=0.15,lw=1)
plt.ylim(15.,22.)
plt.axhline(g[-6][7]+g2[-6][1],color='k',ls='dashed')
plt.legend(loc=1,frameon=False)
plt.xlabel(r'$\mathrm{Distance/\,kpc}$')
plt.ylabel(r'$\log_{10}(J(0.5^\circ)/\,\mathrm{GeV^2\,cm}^{-5})$')
plt.savefig('corr_plot_obl_JD.pdf',bbox_inches='tight')
plt.clf()

## ============================================================================
## 5. J factor against dwarf ranked by ellipticity with error-bars for
##	  spherical, oblate and prolate cases
f=plt.figure(figsize=[6.6,2.5])
plt.errorbar(np.arange(len(g)),g.T[7],xerr=None,yerr=[g.T[8],g.T[9]],fmt='d',color=sns.color_palette()[2],ms=3,label='Spherical',capsize=0)
plt.errorbar(np.arange(len(g))-0.15,g.T[7]+g2.T[1],xerr=None,yerr=[np.sqrt(g.T[8]**2+g2.T[2]**2),np.sqrt(g.T[9]**2+g2.T[3]**2)],fmt='o',color=sns.color_palette()[0],ms=2,label='Oblate',capsize=0)
plt.errorbar(np.arange(len(g))+0.15,g.T[7]+g2.T[4],xerr=None,yerr=[np.sqrt(g.T[8]**2+g2.T[5]**2),np.sqrt(g.T[9]**2+g2.T[6]**2)],fmt='s',color='k',ms=2,label='Prolate',capsize=0)
plt.ylim(14.,23.)
plt.xticks(np.arange(len(g)), map(lambda i:posh_names[i],nam), rotation='vertical')
[i.set_color(sns.color_palette()[2]) for i in plt.gca().get_xticklabels()[:3]]
## Add horizontal line for Ret II errorbar
l=plt.axhline(g[6][7]+np.sqrt(g[6][9]**2+g2[6][6]**2)+g2[6][4],color='k',ls='dashed',alpha=0.5)
plt.annotate(r'Reticulum II $+1\sigma$ prolate',xy=(len(g)-6,g[6][7]+np.sqrt(g[6][9]**2+g2[6][6]**2)+g2[6][4]),horizontalalignment='left',verticalalignment='bottom',color='k')
l.set_dashes((2,1))
l=plt.axhline(g[6][7]+np.sqrt(g[6][9]**2+g2[6][3]**2)+g2[6][1],color='k',ls='dashed',alpha=0.5)
plt.annotate(r'Reticulum II $+1\sigma$ oblate',xy=(len(g)-6,g[6][7]+np.sqrt(g[6][9]**2+g2[6][3]**2)+g2[6][1]-0.1),horizontalalignment='left',verticalalignment='top',color='k')
l.set_dashes((2,1))
plt.legend(loc=1,frameon=False,ncol=3)
plt.xlim(-1.,len(g))

## Display their ranks
temp = (-g.T[7]-g.T[9]).argsort()
ranks = np.empty(len(g), int)
ranks[temp] = np.arange(len(g))
print ranks,-g.T[7]-g.T[9]
for n,i in enumerate(ranks):
	plt.annotate(i+1,xy=(n,14.9),horizontalalignment='center',verticalalignment='center',color=sns.color_palette()[2])

## Display their ranks in prolate case
temp = (-g.T[7]-g2.T[4]-np.sqrt(g.T[9]**2+g2.T[6]**2)).argsort()
ranks = np.empty(len(g), int)
ranks[temp] = np.arange(len(g))
print ranks,-g.T[7]-g2.T[4]-np.sqrt(g.T[9]**2+g2.T[6]**2)
for n,i in enumerate(ranks):
	plt.annotate(i+1,xy=(n,14.5),horizontalalignment='center',verticalalignment='center',color='k')

plt.ylabel(r'$\log_{10}(J(0.5^\circ)/\,\mathrm{GeV^2\,cm}^{-5})$')
plt.savefig('flattening_corr.pdf',bbox_inches='tight')

## ============================================================================
## 5. D factor against dwarf ranked by ellipticity with error-bars for
##	  spherical, oblate and prolate cases
plt.clf()
f=plt.figure(figsize=[6.6,2.5])
plt.errorbar(np.arange(len(g)),g.T[13],xerr=None,yerr=[g.T[14],g.T[15]],fmt='d',color=sns.color_palette()[2],ms=3,label='Spherical',capsize=0)
plt.errorbar(np.arange(len(g))-0.15,g.T[13]+g2.T[7],xerr=None,yerr=[np.sqrt(g.T[14]**2+g2.T[8]**2),np.sqrt(g.T[15]**2+g2.T[9]**2)],fmt='o',color=sns.color_palette()[0],ms=2,label='Oblate',capsize=0)
plt.errorbar(np.arange(len(g))+0.15,g.T[13]+g2.T[10],xerr=None,yerr=[np.sqrt(g.T[14]**2+g2.T[11]**2),np.sqrt(g.T[15]**2+g2.T[12]**2)],fmt='s',color='k',ms=2,label='Prolate',capsize=0)
plt.ylim(15.,20.)
plt.xticks(np.arange(len(g)), map(lambda i:posh_names[i],nam), rotation='vertical')
[i.set_color(sns.color_palette()[2]) for i in plt.gca().get_xticklabels()[:3]]
plt.legend(loc=1,frameon=False,ncol=3)
plt.xlim(-1.,len(g))
plt.ylabel(r'$\log_{10}(D(0.5^\circ)/\,\mathrm{GeV\,cm}^{-2})$')
plt.savefig('flattening_corr_D.pdf',bbox_inches='tight')

## ============================================================================
