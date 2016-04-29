import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.special import gamma as Gamma
from scipy.integrate import quad

G = 4.300918e-6 ## in units solar mass, km/s kpc

GEV2cm5toMsol2kpc5 = 2.2482330e-07
GEVcm2toMsolkpc2 = 8.5358230e-15

def integrate_J_spherical_alphabetagamma(thetamax,D,rho0,rs,alpha,beta,gamma,rt):
	''' J for spherical alpha,beta,gamma model '''
	def rho(r):
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*np.sqrt(1-np.tanh(r/rt)**2)
	def J(ll,th):
		z = ll
		b = np.tan(th)*D
		x = np.sqrt(b*b+z*z)
		return np.sin(th)*(rho(x)**2)
	return np.log10(rho0*rho0*2.*np.pi*quad(lambda y: quad(lambda z: J(y,z), 0., thetamax)[0],-np.inf,np.inf)[0]/GEV2cm5toMsol2kpc5)

def integrate_J_farfield_spherical_alphabetagamma(thetamax,D,rho0,rs,alpha,beta,gamma,rt):
	''' J for spherical alpha,beta,gamma model in far-field limit'''
	def rho(r):
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*np.sqrt(1-np.tanh(r/rt)**2)
	def J(ll,b):
		z = ll
		x = np.sqrt(b*b+z*z)
		return b*(rho(x)**2)
	return np.log10(rho0*rho0*2.*np.pi*quad(lambda y: quad(lambda z: J(y,z), 0., thetamax*D)[0],-np.inf,np.inf)[0]/D/D/GEV2cm5toMsol2kpc5)


def integrate_D_spherical_alphabetagamma(thetamax,D,rho0,rs,alpha,beta,gamma,rt):
	''' D for spherical alpha,beta,gamma model'''
	def rho(r):
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*np.sqrt(1-np.tanh(r/rt)**2)
	def J(ll,th):
		z = ll
		b = np.tan(th)*D
		x = np.sqrt(b*b+z*z)
		return np.sin(th)*rho(x)
	return np.log10(rho0*2.*np.pi*quad(lambda y: quad(lambda z: J(y,z), 0., thetamax)[0],-np.inf,np.inf)[0]/GEVcm2toMsolkpc2)

def integrate_D_farfield_spherical_alphabetagamma(thetamax,D,rho0,rs,alpha,beta,gamma,rt):
	''' D for spherical alpha,beta,gamma model in far-field limit'''
	def rho(r):
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*np.sqrt(1-np.tanh(r/rt)**2)
	def J(ll,b):
		z = ll
		x = np.sqrt(b*b+z*z)
		return b*rho(x)
	return np.log10(rho0*2.*np.pi*quad(lambda y: quad(lambda z: J(y,z), 0., thetamax*D)[0],-np.inf,np.inf)[0]/D/D/GEVcm2toMsolkpc2)


def integrate_rho_spherical_alphabetagamma(R,rho0,rs,alpha,beta,gamma,rt):
	def rho(r):
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*np.sqrt(1-np.tanh(r/rt)**2)
	def J(x):
		return x*x*rho(x)
	return 4.*np.pi*rho0*quad(J, 0., R)[0]

def asymmetric_gaussian_samples(mean,sigma,N=1):
	''' sigmas = [lower error bar, upper error bar] '''
	updown=(np.random.uniform(size=N)>0.5)
	sigmas=[sigma[i] for i in updown]
	return mean+(2*updown-1.)*np.fabs(np.random.normal(loc=0.,scale=sigmas))

def barlow_asymmetric_gaussian_samples(mean,sigma,N=1):
	## Taken from Barlow (2003)
	## This produces very asymmetric looking distributions with sharp cut-offs
	## on the smaller error side
	## The bifurcated (or dimidated as Barlow corrects us) Gaussian
	## (implemented above) is in my mind better.
    alpha = .5*(sigma[1]-sigma[0])
    sig = .5*(sigma[1]+sigma[0])
    u = np.random.normal(loc=0.,scale=1.,size=N)
    return sig*u+mean+alpha*u*u

def HernquistX(s):
    """
      Computes X function from equations (33) & (34) of Hernquist (1990)
    """
    if(s<0.):
        raise ValueError("s must be positive in Hernquist X function")
    elif(s<1.):
        return np.log((1+np.sqrt(1-s*s))/s)/np.sqrt(1-s*s)
    elif(s==1.):
        return 1.
    else:
        return np.arccos(1./s)/np.sqrt(s*s-1)

def wyns_formulaJ_NFW(rho0,r_s,distance,angle):
	''' Analytic integration of J factor for NFW '''
	Delta2 = r_s**2-distance**2*angle**2
	X = distance*angle/r_s
	J = 2.*distance*angle*(7.*distance*r_s**3*angle-4.*distance**3*r_s*angle**3+3.*np.pi*Delta2**2)+6./r_s*(2*Delta2**3-2*r_s**4*Delta2-distance**4*r_s**2*angle**4)*np.array(map(lambda s:HernquistX(s),X))
	J *= np.pi*rho0**2*r_s**2/(3.*distance**2*Delta2**2)
	return np.log10(J/GEV2cm5toMsol2kpc5)

def wyns_formulaJ_NFW_data(sigma_los,r_half,distance,angle,r_s,walker_or_wolf="wolf"):
	'''
	J factor from M_half for NFW profile
	sigma_los in km/s, r_half in pc, distance in kpc, angle in deg, r_s in kpc
	'''
	r_half=0.001*r_half
	angle=np.deg2rad(angle)
	delta_Omega = 2.*np.pi*(1-np.cos(angle))

	if(walker_or_wolf=="wolf"):
		Mhalf = 4.*sigma_los**2*r_half/G
		r_half=4./3.*r_half
		rho0 = Mhalf/4./np.pi/r_s**3/(np.log((r_s+r_half)/r_s)-r_half/(r_s+r_half))
		return wyns_formulaJ_NFW(rho0,r_s,distance,angle)
	else:
		Mhalf = 2.5*sigma_los**2*r_half/G
		rho0 = Mhalf/4./np.pi/r_s**3/(np.log((r_s+r_half)/r_s)-r_half/(r_s+r_half))
		return wyns_formulaJ_NFW(rho0,r_s,distance,angle)

def wyns_formulaD_NFW(rho0,r_s,distance,angle):
	''' Analytic integration of J factor for NFW '''
	X = distance*angle/r_s
	D = np.log(X/2.)+np.array(map(lambda s:HernquistX(s),X))
	D *= 4.*np.pi*rho0*r_s**3/distance**2
	return np.log10(D/GEVcm2toMsolkpc2)

def wyns_formulaD_NFW_data(sigma_los,r_half,distance,angle,r_s,walker_or_wolf="wolf"):
	'''
	D factor from M_half for NFW profile
	sigma_los in km/s, r_half in pc, distance in kpc, angle in deg, r_s in kpc
	'''
	r_half=0.001*r_half
	angle=np.deg2rad(angle)
	delta_Omega = 2.*np.pi*(1-np.cos(angle))
	if(walker_or_wolf=="wolf"):
		Mhalf = 4.*sigma_los**2*r_half/G
		r_half=4./3.*r_half
		rho0 = Mhalf/4./np.pi/r_s**3/(np.log((r_s+r_half)/r_s)-r_half/(r_s+r_half))
		return wyns_formulaD_NFW(rho0,r_s,distance,angle)
	else:
		Mhalf = 2.5*sigma_los**2*r_half/G
		rho0 = Mhalf/4./np.pi/r_s**3/(np.log((r_s+r_half)/r_s)-r_half/(r_s+r_half))

		return wyns_formulaD_NFW(rho0,r_s,distance,angle)

def wyns_formulaJ(sigma_los,r_half,distance,angle,gamma=1.,walker_or_wolf="wolf"):
	'''
	J factor from M_half for power-law profile (slope = gamma)
	sigma_los in km/s, r_half in pc, distance in kpc, angle in deg, r_s in kpc
	'''
	r_half=0.001*r_half
	angle=np.deg2rad(angle)
	delta_Omega = 2.*np.pi*(1-np.cos(angle))

	fac = (2.5/4.)**2

	if(walker_or_wolf=="wolf"):
		fac=(0.25*(27./16.)*(4./3.)**gamma)**2

	if(gamma!=1. and gamma>.5 and gamma<1.5):
		factor = 2.*(3.-gamma)**2*Gamma(gamma-0.5)/(np.pi**(2-gamma)*(3.-2*gamma)*Gamma(gamma))
		return np.log10(factor*sigma_los**4*delta_Omega**(1.5-gamma)/G**2*distance**(1-2.*gamma)*r_half**(2*gamma-4.)/GEV2cm5toMsol2kpc5)+np.log10(fac)
	else:
		return np.log10(8./np.sqrt(np.pi)*sigma_los**4*np.sqrt(delta_Omega)/G**2/distance/(r_half**2)/GEV2cm5toMsol2kpc5)+np.log10(fac)


def wyns_formulaD(sigma_los,r_half,distance,angle,gamma=1.,walker_or_wolf="wolf"):
	'''
	D factor from M_half for power-law profile (slope = gamma)
	sigma_los in km/s, r_half in pc, distance in kpc, angle in deg, r_s in kpc
	'''
	r_half=0.001*r_half
	angle=np.deg2rad(angle)
	delta_Omega = 2.*np.pi*(1-np.cos(angle))

	fac = (2.5/4.)

	if(walker_or_wolf=="wolf"):
		fac=(0.25*(27./16.)*(4./3.)**gamma)

	if(gamma>1. and gamma<3.):
		factor = 2.*Gamma(gamma*0.5-0.5)/(np.pi**(1-0.5*gamma)*Gamma(gamma*0.5))
		return np.log10(factor*sigma_los**2*delta_Omega**(1.5-gamma*0.5)/G*distance**(1.-gamma)*r_half**(gamma-2.)/GEVcm2toMsolkpc2)+np.log10(fac)

def sample_errorsJ(sigma_los,esigma_los,r_half,er_half,distance,edistance,angle,eangle,gamma=1.,N=1000,nfw=-1.,walker_or_wolf="wolf"):
	''' Samples from sigma_los (km/s), r_half (pc), distance (kpc) and angle (deg) pdfs (gaussians) and returns median J value and pm 1 sigma '''
	if(esigma_los[0]==0.):
		## In this case sigma_los is the 95% upper limit. We sample from a uniform distribution from 0.1 km/s to sigma_los/0.95
		s=np.random.uniform(0.1,sigma_los/0.95,N)
	else:
		# s=asymmetric_gaussian_samples(sigma_los,esigma_los,N)
		s=np.exp(asymmetric_gaussian_samples(np.log(sigma_los),esigma_los/sigma_los,N))
	# r=asymmetric_gaussian_samples(r_half,er_half,N)
	r=np.exp(asymmetric_gaussian_samples(np.log(r_half),er_half/r_half,N))
	a=np.exp(asymmetric_gaussian_samples(np.log(angle),eangle/angle,N))
	d=np.random.normal(loc=distance,scale=edistance,size=N)
	if(nfw>0.):
		wf = wyns_formulaJ_NFW_data(s,r,d,a,nfw,walker_or_wolf)
	else:
		wf = wyns_formulaJ(s,r,d,a,gamma,walker_or_wolf)
	mean=np.nanmedian(wf)
	return np.array([mean,mean-np.nanpercentile(wf,15.87),np.nanpercentile(wf,84.13)-mean])

def wyns_formulaJ_error_sample(data,gamma=1.,gammaarray=None,angle='Max',nfw=[0.],N=1000,geo_factor=True,walker_or_wolf="wolf"):
	''' Performs J sampling for a set of data '''
	if(len(nfw)<len(data)):
		nfw=-1.*np.ones(len(data))
	angles=data.theta_max
	angerrs=[[1e-15,1e-15] for i in range(len(data))]
	if(angle=='Half'):
		angles=data.theta_half
		angerrs=[[data.etheta_half2[i],data.etheta_half1[i]] for i in range(len(data))]
	if(angle=='Half_05'):
		angles=0.5*np.ones(len(data))
		angerrs=[[1e-15,1e-15] for i in range(len(data))]
	geof=np.ones(len(data))
	if geo_factor:
		geof = np.sqrt(1.-data.ellip)
	if(isinstance(gammaarray,np.ndarray)):
		return np.array([
		            sample_errorsJ(data.sigma_los[i],
		                           [data.esigma_los2[i],data.esigma_los1[i]],
		                           data.R_half[i]*geof[i],
		                           [data.eR_half2[i]*geof[i],
		                           	data.eR_half1[i]*geof[i]],
		                           data.D[i],
		                           data.eD[i],
		                           angles[i],
		                           angerrs[i],
		                           gammaarray[i],
		                           N=N,
		                           nfw=nfw[i],
		                           walker_or_wolf=walker_or_wolf) for i in range(len(data))])
	return np.array([sample_errorsJ(data.sigma_los[i],
	                				[data.esigma_los2[i],data.esigma_los1[i]],
	                				data.R_half[i]*geof[i],
	                			   [data.eR_half2[i]*geof[i],
	                			   	data.eR_half1[i]*geof[i]],
	                				data.D[i],
	                				data.eD[i],
	                				angles[i],
	                				angerrs[i],
	                				gamma,
	                				N=N,
	                				nfw=nfw[i],
	                				walker_or_wolf=walker_or_wolf) for i in range(len(data))])

def sample_errorsD(sigma_los,esigma_los,r_half,er_half,distance,edistance,angle,eangle,gamma=1.,N=1000,nfw=-1.,walker_or_wolf="wolf"):
	''' Samples from sigma_los (km/s), r_half (pc), distance (kpc) and angle (deg) pdfs (gaussians) and returns median D value and pm 1 sigma '''
	if(esigma_los[0]==0.):
		## In this case sigma_los is the 95% upper limit. We sample from a uniform distribution from 0.1 km/s to sigma_los/0.95
		s=np.random.uniform(0.1,sigma_los,N)
	else:
		# s=asymmetric_gaussian_samples(sigma_los,esigma_los,N)
		s=np.exp(asymmetric_gaussian_samples(np.log(sigma_los),esigma_los/sigma_los,N))
	# r=asymmetric_gaussian_samples(r_half,er_half,N)
	r=np.exp(asymmetric_gaussian_samples(np.log(r_half),er_half/r_half,N))
	a=np.exp(asymmetric_gaussian_samples(np.log(angle),eangle/angle,N))
	d=np.random.normal(loc=distance,scale=edistance,size=N)
	if(nfw>0.):
		wf = wyns_formulaD_NFW_data(s,r,d,a,nfw,walker_or_wolf)
	else:
		wf = wyns_formulaD(s,r,d,a,gamma,walker_or_wolf)
	mean=np.nanmedian(wf)
	return np.array([mean,mean-np.nanpercentile(wf,15.87),np.nanpercentile(wf,84.13)-mean])

def wyns_formulaD_error_sample(data,gamma=1.,gammaarray=None,angle='Max',nfw=[0.],N=1000,geo_factor=True,walker_or_wolf="wolf"):
	''' Performs D sampling for a set of data '''
	if(len(nfw)<len(data)):
		nfw=-1.*np.ones(len(data))
	angles=data.theta_max
	angerrs=[[1e-15,1e-15] for i in range(len(data))]
	if(angle=='Half'):
		angles=data.dtheta_half
		angerrs=[[data.edtheta_half2[i],data.edtheta_half1[i]] for i in range(len(data))]
	if(angle=='Half_05'):
		angles=0.5*np.ones(len(data))
		angerrs=[[1e-15,1e-15] for i in range(len(data))]
	geof=np.ones(len(data))
	if geo_factor:
		geof = np.sqrt(1.-data.ellip)
	if(isinstance(gammaarray,np.ndarray)):
		return np.array([
		         sample_errorsD(data.sigma_los[i],
		                        [data.esigma_los2[i],data.esigma_los1[i]],
		                        data.R_half[i]*geof[i],
		                        [data.eR_half2[i]*geof[i],
		                         data.eR_half1[i]*geof[i]],
		                        data.D[i],
		                        data.eD[i],
		                        angles[i],
		                        angerrs[i],
		                        gammaarray[i],
		                        N=N,
		                        nfw=nfw[i],
		                        walker_or_wolf=walker_or_wolf) for i in range(len(data))])
	return np.array([sample_errorsD(data.sigma_los[i],
	                				[data.esigma_los2[i],data.esigma_los1[i]],
	                				data.R_half[i]*geof[i],
	                			   [data.eR_half2[i]*geof[i],
	                			    data.eR_half1[i]*geof[i]],
	                				data.D[i],
	                				data.eD[i],
	                				angles[i],
	                				angerrs[i],
	                				gamma,
	                				N=N,
	                				nfw=nfw[i],
	                				walker_or_wolf=walker_or_wolf) for i in range(len(data))])
# def add_thetas(ax,xrang,thetalist):
# 	ylim=ax.get_ylim()
# 	ax.set_ylim(ylim[0]-0.5,ylim[1])
# 	for x,t in zip(xrang,thetalist):
# 		ax.annotate(str(t)+r'$^\circ$',xy=(x,ylim[0]),horizontalalignment='center',verticalalignment='bottom',rotation=90)


# def make_table(data):
# 	WEJ2 = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half_05',N=10000,nfw=5.*data['R_half']/1000.)
# 	WED2 = wyns_formulaD_error_sample(data,gamma=1.,angle='Half_05',N=10000,nfw=5.*data['R_half']/1000.)
# 	WEJ3 = wyns_formulaJ_error_sample(data,gamma=1.,angle='Max',N=10000,nfw=5.*data['R_half']/1000.)
# 	WED3 = wyns_formulaD_error_sample(data,gamma=1.,angle='Max',N=10000,nfw=5.*data['R_half']/1000.)
# 	# outfile=open('dwarfs_Jfactors.dat','w')
# 	# outfile.write('\\begin{tabular}{lccccccc}\n')
# 	# outfile.write('\\hline\n\\hline\n')
# 	# outfile.write('Name & $\\theta_\mathrm{max}$ & $\\theta_{0.5}$ & $\\theta_{0.5, \mathrm{decay}}$ & $\log_{10} J(\\theta_\mathrm{max})$ & $\log_{10} J(0.5^\circ)$ & $\log_{10} D(\\theta_\mathrm{max})$ & $\log_{10} D(0.5^\circ)$\\\\ \n')
# 	# outfile.write('& [$^\circ$] & [$^\circ$] & [$^\circ$] & [$\mathrm{GeV^2\,cm}^{-5}$] & [$\mathrm{GeV^2\,cm}^{-5}$] & [$\mathrm{GeV\,cm}^{-2}$] & [$\mathrm{GeV\,cm}^{-2}$]\\\\\n')
# 	# outfile.write('\\hline\n')
# 	# for i in range(len(WEJ)):
# 	# 	string= str(data['Name'][i])+" & $"+\
# 	# 			str(data['theta_max'][i])+"$&$"+\
# 	# 			str(data['theta_half'][i])+"$&$"+\
# 	# 			str(data['dtheta_half'][i])+"$&"+\
# 	# 			"$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WEJ3[i][0],WEJ3[i][1],WEJ3[i][2])
# 	# 	if(i>21):
# 	# 		string+="-&"
# 	# 	else:
# 	# 		string+="$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WEJ2[i][0],WEJ2[i][1],WEJ2[i][2])
# 	# 	string+="$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WED3[i][0],WED3[i][1],WED3[i][2])
# 	# 	if(i>21):
# 	# 		string+="-&"
# 	# 	else:
# 	# 		string+="$%0.2f_{-%0.2f}^{+%0.2f}$"%(WED2[i][0],WEJ3[i][1],WEJ3[i][2])+"\\\\\n"
# 	# 	if(i==7 or i==22):
# 	# 		outfile.write('\\hline\n')
# 	# 	outfile.write(string)
# 	# outfile.write('\\hline\n')
# 	# outfile.write('\end{tabular}\n')
# 	# outfile.close()

# 	outfile=open('dwarfs_Jfactors.dat','w')
# 	outfile.write('\\begin{tabular}{lccccc}\n')
# 	outfile.write('\\hline\n\\hline\n')
# 	outfile.write('Name & $\\theta_\mathrm{max}$ & $\log_{10} J(\\theta_\mathrm{max})$ & $\log_{10} J(0.5^\circ)$ & $\log_{10} D(\\theta_\mathrm{max})$ & $\log_{10} D(0.5^\circ)$\\\\ \n')
# 	outfile.write('& [$^\circ$] & [$\mathrm{GeV^2\,cm}^{-5}$] & [$\mathrm{GeV^2\,cm}^{-5}$] & [$\mathrm{GeV\,cm}^{-2}$] & [$\mathrm{GeV\,cm}^{-2}$]\\\\\n')
# 	outfile.write('\\hline\n')
# 	for i in range(len(WEJ2)):
# 		string= str(data['Name'][i])+" & $"+\
# 				str(data['theta_max'][i])+"$&"+\
# 				"$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WEJ3[i][0],WEJ3[i][1],WEJ3[i][2])
# 		if(i>21):
# 			string+="-&"
# 		else:
# 			string+="$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WEJ2[i][0],WEJ2[i][1],WEJ2[i][2])
# 		string+="$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WED3[i][0],WED3[i][1],WED3[i][2])
# 		if(i>21):
# 			string+="-"+"\\\\\n"
# 		else:
# 			string+="$%0.2f_{-%0.2f}^{+%0.2f}$"%(WED2[i][0],WEJ3[i][1],WEJ3[i][2])+"\\\\\n"
# 		if(i==7 or i==22):
# 			outfile.write('\\hline\n')
# 		outfile.write(string)
# 	outfile.write('\\hline\n')
# 	outfile.write('\end{tabular}\n')
# 	outfile.close()

# if __name__ == '__main__':
# 	data = pd.read_csv('data.dat',sep=' ')
# 	make_table(data)

# 	gs_gammas=np.genfromtxt('geringer_sameth_gamma.dat',skip_header=49)
# 	# for i in range(len(gs_gammas)):
# 	# 	if(gs_gammas[i][23]<0.5):
# 	# 		gs_gammas[i][23]=0.50005
# 	cd=data[data.Class=='CD']
# 	uf=data[data.Class=='UF']
# 	labelrange=np.linspace(0.,len(data),len(data))
# 	labelscd=labelrange[:len(cd)]
# 	labelsuf=labelrange[len(cd):]
# 	f,a=plt.subplots(2,4,figsize=(16,8))
# 	plt.subplots_adjust(hspace=0.5)
# 	for ai in a:
# 		for aj in ai:
# 			aj.set_xticks(labelrange)
# 			aj.set_xticklabels(data.Name.values,rotation=90)
# 			aj.set_xlim(labelrange[0]-1,labelrange[-1]+1)

# 	for i in a[1]:
# 		ls=i.axvline(labelscd[-1]+.5,c='k',ls='dashed')
# 		ls.set_dashes((2,1))
# 		ls=i.axvline(labelsuf[13]+.5,c='k',ls='dashed')
# 		ls.set_dashes((2,1))

# 	a[0][0].errorbar(labelscd,cd.D,yerr=cd.eD,fmt='.')
# 	a[0][0].errorbar(labelsuf,uf.D.values,yerr=uf.eD.values,fmt='.')
# 	a[0][0].set_ylabel(r'Distance/kpc')
# 	a[0][1].errorbar(labelscd,cd.R_half,yerr=[cd.eR_half2,cd.eR_half1],fmt='.')
# 	a[0][1].errorbar(labelsuf,uf.R_half,yerr=[uf.eR_half2,uf.eR_half1],fmt='.')
# 	a[0][1].set_ylabel(r'$R_{\mathrm{half}}/\mathrm{pc}$')
# 	a[0][2].errorbar(labelscd,cd.sigma_los,yerr=[cd.esigma_los2,cd.esigma_los1],fmt='.')
# 	a[0][2].errorbar(labelsuf,uf.sigma_los,yerr=[uf.esigma_los2,uf.esigma_los1],fmt='.')
# 	a[0][2].arrow(labelsuf[9],uf.sigma_los.values[9],0.,-0.5,fc=sns.color_palette()[1],ec=sns.color_palette()[1],head_length=0.2,head_width=0.3)
# 	a[0][2].arrow(labelsuf[15],uf.sigma_los.values[15],0.,-0.5,fc=sns.color_palette()[1],ec=sns.color_palette()[1],head_length=0.2,head_width=0.3)
# 	a[0][2].arrow(labelsuf[17],uf.sigma_los.values[17],0.,-0.5,fc=sns.color_palette()[1],ec=sns.color_palette()[1],head_length=0.2,head_width=0.3)
# 	a[0][2].set_ylabel(r'$\sigma_{\mathrm{los}}/\mathrm{km\,s}^{-1}$')

# 	a[1][0].errorbar(labelscd,cd.Jmax,yerr=[cd.eJmax2,cd.eJmax1],fmt='.',color='k')
# 	a[1][0].errorbar(labelsuf,uf.Jmax,yerr=[uf.eJmax2,uf.eJmax1],fmt='.',color='k')
# 	WE = wyns_formulaJ_error_sample(data,gamma=1.)
# 	for i in range(len(data)):
# 		a[1][0].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2])
# 	# WE = wyns_formulaJ_error_sample(data,gamma=0.75)
# 	# for i in range(len(data)):
# 	# 	a[1][0].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[3])
# 	WE = wyns_formulaJ_error_sample(data,gamma=0.51)
# 	for i in range(len(data)):
# 		a[1][0].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4])

# 	WE = wyns_formulaJ_error_sample(data,gamma=1.,nfw=5.*data['R_half']/1000.)
# 	for i in range(len(data)):
# 		a[1][0].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0])

# 	add_thetas(a[1][0],labelrange,data.theta_max)
# 	a[1][0].set_ylabel(r'$\log_{10}(J_\mathrm{max}/\,\mathrm{GeV^2\,cm}^{-5})$')

# 	a[1][1].errorbar(labelscd,cd.Jmax.values-np.log10(2.),yerr=[cd.eJmax2,cd.eJmax1],fmt='.',label="",color='k')
# 	a[1][1].errorbar(labelsuf,uf.Jmax.values-np.log10(2.),yerr=[uf.eJmax2,uf.eJmax1],fmt='.',label="",color='k')
# 	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half')
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma=1$'
# 		a[1][1].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2],label=label)
# 	# WE = wyns_formulaJ_error_sample(data,gamma=0.75,angle='Half')
# 	# for i in range(len(data)):
# 	# 	label=None
# 	# 	if(i==0):
# 	# 		label=r'$\gamma=0.75$'
# 	# 	a[1][1].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[3],label=label)
# 	WE = wyns_formulaJ_error_sample(data,gamma=0.51,angle='Half')
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma=0.51$'
# 		a[1][1].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4],label=label)

# 	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half',nfw=5.*data['R_half']/1000.)
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'NFW'
# 		a[1][1].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0],label=label)

# 	gammas = gs_gammas.T[23]
# 	while(len(gammas)<len(data)):
# 		gammas = np.append(gammas,0.8)

# 	WE = wyns_formulaJ_error_sample(data,gammaarray=gammas,angle='Half')
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma_\mathrm{GS}$'
# 		a[1][1].fill_between([labelrange[i]-0.3,labelrange[i]+0.3], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=1.,facecolor="None",label=label)
# 	add_thetas(a[1][1],labelrange,data.theta_half)
# 	a[1][1].legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.0))
# 	a[1][1].set_ylabel(r'$\log_{10}(J_\mathrm{half}/\,\mathrm{GeV^2\,cm}^{-5})$')

# 	a[1][2].errorbar(labelscd,cd.dJmax.values-np.log10(2.),yerr=[cd.eJmax2,cd.edJmax1],fmt='.',color='k')
# 	a[1][2].errorbar(labelsuf,uf.dJmax.values-np.log10(2.),yerr=[uf.edJmax2,uf.edJmax1],fmt='.',color='k')

# 	WE = wyns_formulaD_error_sample(data,gamma=1.)
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma=1.$'
# 		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2],label=label)
# 	# WE = wyns_formulaD_error_sample(data,gamma=1.25)
# 	# for i in range(len(data)):
# 	# 	label=None
# 	# 	if(i==0):
# 	# 		label=r'$\gamma=1.25$'
# 	# 	a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[3],label=label)
# 	WE = wyns_formulaD_error_sample(data,gamma=1.49)
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma=1.49$'
# 		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4],label=label)

# 	WE = wyns_formulaD_error_sample(data,gamma=1.,nfw=5.*data['R_half']/1000.)
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'NFW'
# 		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0],label=label)

# 	WE = wyns_formulaD_error_sample(data,gammaarray=gammas,angle='Half')
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma_\mathrm{GS}$'
# 		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=1.,facecolor="None",label=label)
# 	add_thetas(a[1][2],labelrange,data.dtheta_half)
# 	a[1][2].legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.0))
# 	a[1][2].set_ylabel(r'$\log_{10}(D_\mathrm{half}/\,\mathrm{GeV\,cm}^{-2})$')

# 	a[1][3].errorbar(labelscd,cd.Jhalf.values,yerr=[cd.eJhalf2,cd.eJhalf1],fmt='.',label="",color='k')
# 	a[1][3].errorbar(labelsuf,uf.Jhalf.values,yerr=[uf.eJhalf2,uf.eJhalf1],fmt='.',label="",color='k')
# 	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half'	)
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma=1$'
# 		a[1][3].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2],label=label)
# 	# WE = wyns_formulaJ_error_sample(data,gamma=0.75,angle='Half_05')
# 	# for i in range(len(data)):
# 	# 	label=None
# 	# 	if(i==0):
# 	# 		label=r'$\gamma=0.75$'
# 	# 	a[1][3].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[3],label=label)


# 	WE = wyns_formulaJ_error_sample(data,gamma=0.51,angle='Half_05')
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma=0.51$'
# 		a[1][3].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4],label=label)

# 	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half_05',nfw=5.*data['R_half']/1000.)
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'NFW'
# 		a[1][3].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0],label=label)

# 	gammas = gs_gammas.T[23]
# 	while(len(gammas)<len(data)):
# 		gammas = np.append(gammas,0.8)

# 	WE = wyns_formulaJ_error_sample(data,gammaarray=gammas,angle='Half_05')
# 	for i in range(len(data)):
# 		label=None
# 		if(i==0):
# 			label=r'$\gamma_\mathrm{GS}$'
# 		a[1][3].fill_between([labelrange[i]-0.3,labelrange[i]+0.3], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=1.,facecolor="None",label=label)
# 	add_thetas(a[1][3],labelrange,np.ones(0.5)*len(data))
# 	a[1][3].legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.0))
# 	a[1][3].set_ylabel(r'$\log_{10}(J(0.5^\circ)/\,\mathrm{GeV^2\,cm}^{-5})$')

# 	plt.savefig('dwarfs_data.pdf',bbox_inches='tight')
