## Generate Figs 2 & 3 of SEG 2016 -- compute_J_ret2_sims.py must have been
## run before this
## ============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
sys.path.append('../../../code/jfactors/')
sys.path.append('../../../code/m2m/')
sys.path.append('../spherical/')
import flattened as fJ
import spherical_Jfactors as sJ
import jfactors_py as cJ
import seaborn as sns
## ============================================================================
## Load Ret 2 data
import ret_2
RetII = ret_2.RetII
r_maj_exp = RetII.r_maj_exp ## arcmin
r_maj = RetII.r_maj #arcmin
Distance = RetII.Distance #kpc
Velocity_dispersion = RetII.Velocity_dispersion #km/s
rh= (r_maj/60./180.*np.pi)*Distance ## in units kpc
## ============================================================================
q = 0.4
alpha = 1.
beta = 3.
gamma = 1.
gf = 1. ## np.sqrt(1.-RetII.e) We don't want the ellipticity correction here
		## as we are comparing to the spherical model
## ============================================================================
## Make figures
def make_plot(data,cols=[11,12,13],label='Face-on',
              output_file='J_profiles_sim_round.pdf'):
	''' Takes the results from ret_results.dat and plots a series of J-factor
		and D-factor profiles for each model '''
	f,a=plt.subplots(2,1,figsize=[3.32,4.])
	plt.subplots_adjust(hspace=0.)
	cm = sns.cubehelix_palette(8,start=.5,rot=-.75,as_cmap=True)
	cNorm = colors.Normalize(vmin=0.4,vmax=2.5)
	sM = cmx.ScalarMappable(norm=cNorm,cmap=cm)
	## 1. For each model in ret_results.dat we plot a profile of the J-factor
	##    against beam angle. ret_results.dat contains the central density,
	##	  scale-radii and tidal radii for each model
	tmax = np.logspace(-2.,np.log10(0.5),20)  ## angle grid
	for d in data[1:]:
		Jvals_r = np.zeros(len(tmax))
		Dvals_r = np.zeros(len(tmax))
		for n,i in enumerate(tmax):
			if(d[0]>1.): ## Prolate case
				ba,ca = 1./d[0],1./d[0]
				print ba,ca,i
				model = cJ.AlphaBetaGammaDensityProfile(np.array([alpha,beta,gamma]),d[cols[0]],d[cols[1]],d[cols[2]],np.array([1.,ba,ca]),True)
				Jvals_r[n] = np.log10(model.J_far_factor(Distance,i,"x")/sJ.GEV2cm5toMsol2kpc5)
				Dvals_r[n] = np.log10(model.D_far_factor(Distance,i,"x")/sJ.GEVcm2toMsolkpc2)
			else: ## Oblate case
				ba,ca = 1.,d[0]
				print ba,ca,i
				model = cJ.AlphaBetaGammaDensityProfile(np.array([alpha,beta,gamma]),d[cols[0]],d[cols[1]],d[cols[2]],np.array([1.,ba,ca]),True)
				Jvals_r[n] = np.log10(model.J_far_factor(Distance,i,"z")/sJ.GEV2cm5toMsol2kpc5)
				Dvals_r[n] = np.log10(model.D_far_factor(Distance,i,"z")/sJ.GEVcm2toMsolkpc2)
		l,=a[0].plot(tmax,Jvals_r,color=sM.to_rgba(d[0]))
		l2,=a[1].plot(tmax,Dvals_r,color=sM.to_rgba(d[0]))
		if(d[0]==1.):
			l.set_dashes((2,1))
			l2.set_dashes((2,1))
	## 2. Also plot the spherical formulae from Paper I
	l,=a[0].plot(tmax,sJ.wyns_formulaJ_NFW_data(Velocity_dispersion,rh*1000.*gf
	             ,Distance,tmax,2.*rh*gf),ls='dashed',color='k')
	l.set_dashes((4,1))
	l,=a[1].plot(tmax,sJ.wyns_formulaD_NFW_data(Velocity_dispersion,rh*1000.*gf,Distance,tmax,2.*rh*gf),ls='dashed',color='k')
	l.set_dashes((4,1))
	## 3. Add the colorbar
	divider = make_axes_locatable(a[0])
	cba = divider.append_axes("top",size="5%",pad=0.)
	cbl =matplotlib.colorbar.ColorbarBase(cba,cmap=cm,norm=cNorm,orientation='horizontal')
	cbl.set_label(r'$q$',labelpad=-30.4)
	cbl.ax.xaxis.set_ticks_position('top')
	a[0].yaxis.get_major_ticks()[0].label1.set_visible(False)
	a[0].set_xticklabels([])
	a[0].set_ylabel(r'$\log_{10}(J(\theta)/\,\mathrm{GeV^2\,cm}^{-5})$')
	a[1].set_xlabel(r'$\theta/\,\mathrm{deg}$')
	a[1].set_ylabel(r'$\log_{10}(D(\theta)/\,\mathrm{GeV\,cm}^{-2})$')
	a[0].text(0.9,0.1,label,horizontalalignment='right',verticalalignment='bottom',transform=a[0].transAxes,fontsize=14)
	plt.savefig(output_file,bbox_inches='tight')

if __name__ == '__main__':
	data = np.genfromtxt('ret_results.dat')
	make_plot(data) ## Face-on case (default)
	make_plot(data,cols=[8,9,10],label='Edge-on',
	          output_file='J_profiles_sim_edge.pdf') ## Edge-on case

## ============================================================================
