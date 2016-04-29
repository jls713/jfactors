## Generate tables 3 & 4 of Sanders, Evans & Geringer-Sameth (2016)
## Must run flattened.py first to generate fit results
## ============================================================================
import numpy as np
import pandas as pd
import sys
sys.path.append('../spherical/')
from spherical_Jfactors import asymmetric_gaussian_samples
from J_D_table import posh_latex_names
## ============================================================================
def compute_corrections(datafile):
	''' Generate tables with columns Name Oblate eO1 eO2 Prolate eP1 eP2 OblateD eOD1 eOD2 ProlateD ePD1 ePD2 -- also output in ascii format'''
	## ========================================================================
	## 1. Load in data and write table headers
	data = pd.read_csv(datafile,sep=' ')
	outfile=open('dwarfs_Jfactors_corr.dat','w')
	outfile2=open('dwarfs_Dfactors_corr.dat','w')
	outfile.write('\\begin{tabular}{lc|ccccc}\n')
	outfile2.write('\\begin{tabular}{lc|ccccc}\n')
	outfile.write('\\hline\n\\hline\n')
	outfile2.write('\\hline\n\\hline\n')
	outfile.write('Name & $\epsilon$ & $\log_{10}(\mathrm{J}_\mathrm{sph}(0.5^\circ))$ & Oblate $\mathcal{F}_\mathrm{J}$ & Prolate $\mathcal{F}_\mathrm{J}$ & $\log_{10}(\mathrm{J}_\mathrm{obl}(0.5^\circ))$ & $\log_{10}(\mathrm{J}_\mathrm{pro}(0.5^\circ))$\\\\ \n')
	outfile2.write('Name & $\epsilon$ & $\log_{10}(\mathrm{D}_\mathrm{sph}(0.5^\circ))$ & Oblate $\mathcal{F}_\mathrm{D}$ & Prolate $\mathcal{F}_\mathrm{D}$ & $\log_{10}(\mathrm{D}_\mathrm{obl}(0.5^\circ))$ & $\log_{10}(\mathrm{D}_\mathrm{pro}(0.5^\circ))$\\\\ \n')
	outfile.write('\\hline\n')
	outfile2.write('\\hline\n')
	dd = np.genfromtxt("../spherical/dwarfs_Jfactors_ascii.dat")

	outfile3=open('dwarfs_Jfactors_corr_ascii.dat','w')
	outfile3.write('Name Oblate eO1 eO2 Prolate eP1 eP2 OblateD eOD1 eOD2 ProlateD ePD1 ePD2\n')

	## ========================================================================
	## 2. Load in fit results from running flattened.py
	dat = np.genfromtxt("fit_results_ascii.dat")[1:]

	size = 5000 ## Number of samples used
	## ========================================================================
	## 3. Compute correction factors for series of samples
	for i in (data.ellip.values).argsort(): #range(len(data)):

		sample=np.exp(asymmetric_gaussian_samples(np.log(data.ellip.values[i]),[data.ellip_e2.values[i]/data.ellip.values[i],data.ellip_e1.values[i]/data.ellip.values[i]],size))
		sample=sample[sample<1.] ## Remove unphysical

		if(data.ellip_e1.values[i]>data.ellip.values[i] or data.ellip_e2.values[i]>data.ellip.values[i]):
			updown=np.count_nonzero(np.random.uniform(size=size)>0.5)
			sample1 = np.random.uniform(high=data.ellip.values[i],size=updown)
			sample2 = data.ellip.values[i]+np.fabs(np.random.normal(loc=0.,scale=data.ellip_e1.values[i],size=size-updown))
			sample = np.concatenate((sample1,sample2))

		## If only upper-bound sample uniformly from 0 to bound
		if(data.ellip_e1.values[i]!=data.ellip_e1.values[i]):
			sample = np.random.uniform(high=data.ellip.values[i],size=size)

		obs_corr = np.log10(1.-sample)*dat[0][0]
		pro_corr = np.log10(1./(1.-sample))*dat[0][1]
		obs_corrD = np.log10(1.-sample)*dat[2][0]
		pro_corrD = np.log10(1./(1.-sample))*dat[2][1]

		## Compute medians and standard deviations
		ocm,ocs1,ocs2 = np.median(obs_corr),np.percentile(obs_corr,15.9),np.percentile(obs_corr,84.1)
		ocs1=ocm-ocs1
		ocs2=ocs2-ocm
		ocmD,ocs1D,ocs2D = np.median(obs_corrD),np.percentile(obs_corrD,15.9),np.percentile(obs_corrD,84.1)
		ocs1D=ocmD-ocs1D
		ocs2D=ocs2D-ocmD
		pcm,pcs1,pcs2 = np.median(pro_corr),np.percentile(pro_corr,15.9),np.percentile(pro_corr,84.1)
		pcs1=pcm-pcs1
		pcs2=pcs2-pcm
		pcmD,pcs1D,pcs2D = np.median(pro_corrD),np.percentile(pro_corrD,15.9),np.percentile(pro_corrD,84.1)
		pcs1D=pcmD-pcs1D
		pcs2D=pcs2D-pcmD

		## Write to files
		string= str(posh_latex_names[data['Name'][i]])+"&"
		if(data.ellip_e1.values[i]==data.ellip_e1.values[i]):
			string+="$%0.2f_{-%0.2f}^{+%0.2f}$"%(data.ellip.values[i],data.ellip_e2.values[i],data.ellip_e1.values[i])+"&"
		else:
			string+="$<%0.2f$"%(data.ellip.values[i])+"&"
		outfile.write(string+
		        "$%0.2f_{-%0.2f}^{+%0.2f}$"%(dd[i][7],dd[i][8],dd[i][9])+\
				"&$%0.3f_{-%0.3f}^{+%0.3f}$"%(ocm,ocs1,ocs2)+\
				"&$%0.3f_{-%0.3f}^{+%0.3f}$"%(pcm,pcs1,pcs2)+\
				"&$%0.2f_{-%0.2f}^{+%0.2f}$"%(dd[i][7]+ocm,np.sqrt(dd[i][8]**2+ocs1**2),np.sqrt(dd[i][9]**2+ocs2**2))+\
				"&$%0.2f_{-%0.2f}^{+%0.2f}$"%(dd[i][7]+pcm,np.sqrt(dd[i][8]**2+pcs1**2),np.sqrt(dd[i][9]**2+pcs2**2))+"\\\\\n")
		outfile2.write(string+
		        "$%0.2f_{-%0.2f}^{+%0.2f}$"%(dd[i][13],dd[i][14],dd[i][15])+\
				"&$%0.3f_{-%0.3f}^{+%0.3f}$"%(ocmD,ocs1D,ocs2D)+\
				"&$%0.3f_{-%0.3f}^{+%0.3f}$"%(pcmD,pcs1D,pcs2D)+\
				"&$%0.2f_{-%0.2f}^{+%0.2f}$"%(dd[i][13]+ocmD,np.sqrt(dd[i][14]**2+ocs1D**2),np.sqrt(dd[i][15]**2+ocs2D**2))+\
				"&$%0.2f_{-%0.2f}^{+%0.2f}$"%(dd[i][13]+pcmD,np.sqrt(dd[i][14]**2+pcs1D**2),np.sqrt(dd[i][15]**2+pcs2D**2))+"\\\\\n")
		string= str(data['Name'][i])+" "+\
				"%0.5f %0.5f %0.5f"%(ocm,ocs1,ocs2)+" "+\
				"%0.5f %0.5f %0.5f"%(pcm,pcs1,pcs2)+" "+\
				"%0.5f %0.5f %0.5f"%(ocmD,ocs1D,ocs2D)+" "+\
				"%0.5f %0.5f %0.5f\n"%(pcmD,pcs1D,pcs2D)
		outfile3.write(string)
	outfile3.close()
	outfile.write('\\hline\n')
	outfile.write('\end{tabular}\n')
	outfile.close()
	outfile2.write('\\hline\n')
	outfile2.write('\end{tabular}\n')
	outfile2.close()

## ============================================================================
if __name__ == '__main__':
	compute_corrections('../data/data.dat')
## ============================================================================
