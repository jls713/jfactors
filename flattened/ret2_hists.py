## Generates Fig. 9 of SEG (2016)
## ============================================================================
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
## ============================================================================
## 1. load in data from the three assumptions
data,datama,datasj=np.genfromtxt('triaxial_results/ReticulumII_nop_hr'),np.genfromtxt('triaxial_results/ReticulumII_ma_hr'),np.genfromtxt('triaxial_results/ReticulumII_sj_hr')
## 2. Use seaborn KDE plot
sns.kdeplot(data.T[4],shade=True,label='Uniform (U)')
sns.kdeplot(datama.T[4],shade=True,ls='dotted',label='Major axis (R)')
sns.kdeplot(datasj.T[4],shade=True,ls='dashed',label='SJ 2016 (T)')
plt.xlabel(r'$\mathcal{F}_\mathrm{J}$')
plt.ylabel(r'$\mathrm{d}N/\mathrm{d}\mathcal{F}_\mathrm{J}$')
plt.xlim(-.5,.5)
plt.legend()
plt.savefig('RetII_hr.pdf',bbox_inches='tight')
## ============================================================================
