## Produce Table 5 of SEG (2016)
## ============================================================================
import numpy as np
import pandas as pd
import sys
import os.path
sys.path.append('/home/jls/work/data/jfactors/spherical/')
from J_D_table import posh_latex_names
## ============================================================================

def load_files(name):
	''' Load in three sample files for dwarf <name> '''
	name='triaxial_results/'+name
	if os.path.isfile(name+'_nop') and os.path.isfile(name+'_ma') and os.path.isfile(name+'_sj'):
		return np.genfromtxt(name+'_nop'),np.genfromtxt(name+'_ma'),np.genfromtxt(name+'_sj')
	else:
		return None,None,None
def write(l):
	''' Output median and \pm 1\sigma errors for correction factors in ascii
		form '''
	l = l.T[4]
	return r'$%0.3f^{+%0.3f}_{-%0.3f}$'%(np.median(l),np.percentile(l,84.1)-np.median(l),np.median(l)-np.percentile(l,15.9))
def write_ascii(l):
	''' Output median and \pm 1\sigma errors for correction factors in latex
		form '''
	l = l.T[4]
	return '%0.3f %0.3f %0.3f '%(np.median(l),np.percentile(l,84.1)-np.median(l),np.median(l)-np.percentile(l,15.9))

## ============================================================================
## 1. Read in data file and write headers to tables
data = pd.read_csv('../../../data/jfactors/data.dat',sep=' ')
ff = open('corr_triax_table.dat','w')
ffa = open('corr_triax_table_ascii.dat','w')
ff.write('\\begin{tabular}{lcccc}\n')
ff.write('\\hline\n\\hline\n')
ff.write('Name & Ellipticity & $\mathcal{F}_{\mathrm{J},U}$& $\mathcal{F}_{\mathrm{J},R}$& $\mathcal{F}_{\mathrm{J},T}$\\\\ \n')
ff.write('\\hline\n')
## 2. Loop over dwarfs and compute median and \pm 1 \sigma for correction factors
for i in data.ellip.argsort():
	d,e,f=load_files(data.Name[i])
	ellip_string='&$%0.2f^{+%0.2f}_{-%0.2f}$&'%(data.ellip[i],data.ellip_e1[i],data.ellip_e2[i])
	if(data.ellip_e1[i]!=data.ellip_e1[i]):
		ellip_string='&$<%0.2f$&'%(data.ellip[i])
	if(d==None):
		ff.write(posh_latex_names[data['Name'][i]]+ellip_string+'NaN&NaN&NaN\\\\\n')
		ffa.write(data['Name'][i]+' %0.2f %0.2f %0.2f '%(data.ellip[i],data.ellip_e1[i],data.ellip_e2[i])+'NaN '*9+'\n')
	else:
		ff.write(posh_latex_names[data['Name'][i]]+ellip_string
		         +write(d)+'&'+write(e)+'&'+write(f)+'\\\\\n')
		ffa.write(data['Name'][i]+' %0.2f %0.2f %0.2f '%(data.ellip[i],data.ellip_e1[i],data.ellip_e2[i])+write_ascii(d)+write_ascii(e)+write_ascii(f)+'\n')
ff.write('\\hline\n\\end{tabular}\n')
ff.close()
ffa.close()
## ============================================================================
