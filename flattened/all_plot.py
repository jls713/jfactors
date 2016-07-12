# -*- coding: utf-8 -*-
## Generates Figure 10 of SEG (2016)
## ============================================================================
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
sys.path.append('/home/jls/work/data/jfactors/spherical/')
from J_D_table import posh_names
## ============================================================================
## 1. Read in data and sort by ellipticity
g = np.genfromtxt("corr_triax_table_ascii.dat")
nam = pd.read_csv("corr_triax_table_ascii.dat",sep=' ',header=None).values.T[0]
temp = (g.T[1]).argsort()
nam = nam[temp]
g = g[temp]
## 2. Plot error bars for each assumption
plt.clf()
f=plt.figure(figsize=[6.6,2.5])
plt.errorbar(np.arange(len(g))-0.15,g.T[4],xerr=None,yerr=[g.T[6],g.T[5]],fmt='o',color=sns.color_palette()[0],ms=3,label='Uniform (U)',capsize=0)
plt.errorbar(np.arange(len(g)),g.T[7],xerr=None,yerr=[g.T[9],g.T[8]],fmt='s',color=sns.color_palette()[1],ms=3,label='Major axis (R)',capsize=0)
plt.errorbar(np.arange(len(g))+0.15,g.T[10],xerr=None,yerr=[g.T[12],g.T[11]],fmt='d',color=sns.color_palette()[2],ms=3,label='SJ 2016 (T)',capsize=0)

plt.xlim(-1.,len(g))
plt.xticks(np.arange(-1.,len(g))+1., map(lambda i:posh_names[i],nam), rotation='vertical')
plt.legend(ncol=3)
l=plt.axhline(0.,color='gray')
l.set_dashes((2,1))
plt.ylim(-.5,1.)
plt.ylabel(r'$\mathcal{F}_\mathrm{J}$')
plt.savefig('corr_triax_all.pdf',bbox_inches='tight')

## ============================================================================
