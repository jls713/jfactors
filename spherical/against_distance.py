# -*- coding: utf-8 -*-
### Generates J against distance for Fig. 1 of Evans, Sanders & Geringer-Sameth (2016)

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

f=plt.figure(figsize=[3.32,2.5])
g = np.genfromtxt("dwarfs_Jfactors_ascii.dat")
plt.errorbar(g[:8].T[1],g[:8].T[7],xerr=g[:8].T[2],yerr=[g[:8].T[8],g[:8].T[9]],fmt='d',color=sns.color_palette()[2],ms=4,label='Classical dwarfs')
plt.errorbar(g[8:23].T[1],g[8:23].T[7],xerr=g[8:23].T[2],yerr=[g[8:23].T[8],g[8:23].T[9]],fmt='o',color='k',ms=3,label='Ultrafaints')
plt.errorbar(g[23:].T[1],g[23:].T[7],xerr=g[23:].T[2],yerr=[g[23:].T[8],g[23:].T[9]],fmt='s',color=sns.color_palette()[0],ms=4,label='New')
plt.legend(frameon=False)
plt.xlabel(r'$\mathrm{Distance/\,kpc}$')
plt.ylabel(r'$\log_{10}(\mathrm{J}(0.5^\circ)/\,\mathrm{GeV^2\,cm}^{-5})$')
plt.savefig('J_against_distance.pdf',bbox_inches='tight')
