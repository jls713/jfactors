## Generates Figs 5, 6 & 7 of SEG 2016
## ============================================================================
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib.patches import Ellipse
import flattened as fJ
from scipy.optimize import curve_fit
import seaborn as sns
## ============================================================================
def read_wyn_table(table):
	''' For reading Wyn's Mathematica tables '''
	f = open(table,'r')
	l = f.read().replace("{","").replace('\t','').replace('*^','e').split('}')
	L = np.array(map(lambda i: np.float64(i.split(', ')),l[:-1]))
	return L
def load_wyn_data():
	def load_wyn_file(name):
		wyndata = np.genfromtxt(name)
		return wyndata[np.argsort(wyndata.T[0])]
	files = ['../../../code/jfactors/data/wyn_factor_n5_rr2.cf.dat',
			'../../../code/jfactors/data/wyn_factor_n5_rr20.cf.dat',
			'../../../code/jfactors/data/wyn_factor_n5_rr200.cf.dat']
	return [load_wyn_file(f) for f in files]
## ============================================================================
## Fitting functions
def func_qp(q, a, d):
	d=0.
	return a*np.log10(q)+d

def func_qm(q, a, d):
	d=0.
	return a*np.log10(q)+d

## ============================================================================
def qpot_from_q(q):
	''' Flattening in potential from density flattening for logarithmic '''
	return .5*np.sqrt(1.+np.sqrt(1.+8.*q*q))
## ============================================================================
## Results for axisymmetric cusps
def geometric_factor(q,gamma):
	''' Geometric factor for infinite axisymmetric cusp '''
	## \int dt (cos^2(t)+sin^2(t)/q^2)^{1/2-gamma}
	Fac = 1./(2.*np.pi*q*q)
	return np.array(
	      map(lambda Q:
	      1./(2.*np.pi*Q*Q)*quad(lambda t:
	      np.power(np.power(np.cos(t),2.)+np.power(np.sin(t),2.)/Q/Q,0.5-gamma)
	      ,0.,2.*np.pi)[0],q))

def Jkin(q,gamma,axis='short'):
	''' Kinematic factor for infinite axisymmetric cusp '''
	F=np.array(
	  map(lambda Q:(.5*
	  quad(lambda t: np.sin(t)**3*(np.sin(t)**2+np.cos(t)**2./Q/Q)**(-gamma/2.)
	       ,0.,np.pi)[0]
	 /quad(lambda t: np.sin(t)*np.cos(t)**2*(np.sin(t)**2+np.cos(t)**2./Q/Q)**(-gamma/2.),0.,np.pi)[0]),q))
	if(axis=='minor'):
		return ((2.*F+1.)/3.)**2
	else:
		return ((2.+1./F)/3.)**2

def check_equation_20(q,gamma=3.):
	'''
		Kinematic factor difference between analytic and numerical for
		gamma = 3
	'''
	Q = q/qpot_from_q(q)
	F=.5*quad(lambda t: np.sin(t)**3*(np.sin(t)**2+np.cos(t)**2./Q/Q)**(-gamma/2.),0.,np.pi)[0]/quad(lambda t: np.sin(t)*np.cos(t)**2*(np.sin(t)**2+np.cos(t)**2./Q/Q)**(-gamma/2.),0.,np.pi)[0]
	F2 = binney_tremaine_virial_ratio(q)
	Qst = 1./Q/Q
	if(gamma==3.):
		if(Qst>1.):
			Q = np.sqrt(Qst-1.)
			G = .5*(Qst*Q-np.sqrt(Qst)*np.arcsinh(Q))/(np.sqrt(Qst)*np.arcsinh(Q)-Q)
		else:
			Q = np.sqrt(1.-Qst)
			G = .5*(Qst*Q*Q-np.sqrt(Qst)*Q*np.arccos(np.sqrt(Qst)))/(np.sqrt(Qst)*Q*np.arccos(np.sqrt(Qst))-Q*Q)
	if(gamma==4.):
		if(Qst>1.):
			Q = np.sqrt(Qst-1.)
			G = (Qst*(Qst-2.)*np.arctan(Q)+Qst*Q)/(Qst*np.arctan(Q)-Q)/2.
		else:
			Q = np.sqrt(1.-Qst)
			G = (Qst*(Qst-2.)*np.arctanh(Q)+Qst*Q)/(Qst*np.arctanh(Q)-Q)/2.
	if(gamma==3.):
		q = np.linspace(0.2,4.)
		FF = np.zeros(len(q))
		FF2 = np.zeros(len(q))
		for n,i in enumerate(q):
			Q = i/qpot_from_q(i)
			FF[n]=.5*quad(lambda t: np.sin(t)**3*(np.sin(t)**2+np.cos(t)**2./Q/Q)**(-gamma/2.),0.,np.pi)[0]/quad(lambda t: np.sin(t)*np.cos(t)**2*(np.sin(t)**2+np.cos(t)**2./Q/Q)**(-gamma/2.),0.,np.pi)[0]
			FF2[n] = binney_tremaine_virial_ratio(i)
		plt.plot(q,np.log10(FF))
		plt.plot(q,np.log10(FF2))
		plt.savefig('tmp.pdf')
		plt.clf()
	if(gamma==3.):
		return F-G,F2-G
	else:
		return F-G

def binney_tremaine_virial_ratio(q):
	''' Table 2.1 of Binney & Tremaine (2008) -- <sigma_xx^2>/<sigma_zz^2>
		for stellar density stratified on spheroidal shells and dark matter
		density stratified on spheroidal shells of the same ellipticity '''
	if(q<1):
		e = np.sqrt(1.-q*q)
		return .5*(np.arcsin(e)/e-q)/(1./q-np.arcsin(e)/e)/q/q
	else:
		e = np.sqrt(1.-1./q/q)
		return .5*(q*q-.5*np.log((1.+e)/(1.-e))/e)/(.5*np.log((1.+e)/(1.-e))/e-1.)/q/q

def round_GF(q,gamma):
	''' Face-on geometric correction factor '''
	return 1./q

def edge_GF(q,gamma):
	''' Edge-on geometric correction factor '''
	return geometric_factor(q,gamma)

def WR(q_star,q_DM,gamma_star):
	''' Ratio of sigma_xx^2 to sigma_zz^2 '''
	q_eff = q_star/qpot_from_q(q_DM)
	return 1.5*np.sqrt(Jkin(q_eff,gamma_star,axis='minor'))-0.5

def edge_factor(q_star,q_DM,gamma_star,gamma_DM,geo_factor=False, finite=True):
	''' Combination of kinematic and geometric factors for J-factor edge-on'''
	q_eff = q_star/qpot_from_q(q_DM)
	gf = 1.
        if(geo_factor):
		gf = (q_star<1.)*q_star**(0.5*(4.-2.*gamma_DM))+(q_star>1.)*q_star**(-0.5*(4.-2.*gamma_DM))
	if not finite:
		return gf*((q_star>1.)*q_DM**(4-2.*gamma_DM)+(q_star<1.)*1.)*geometric_factor(q_DM,gamma_DM)*Jkin(q_eff,gamma_star,axis='major')
	else:
		return gf*((q_star>1.)*q_DM**(4-2.*gamma_DM)+(q_star<1.)*1.)/q_DM*Jkin(q_eff,gamma_star,axis='major')

def round_factor(q_star,q_DM,gamma_star,gamma_DM):
	''' Combination of kinematic and geometric factors for J-factor face-on'''
	q_eff = q_star/qpot_from_q(q_DM)
	return 1./q_DM*Jkin(q_eff,gamma_star,axis='minor')

def edge_factor_D(q_star,q_DM,gamma_star,gamma_DM,geo_factor=False,finite=True):
	''' Combination of kinematic and geometric factors for D-factor edge-on'''
	q_eff = q_star/qpot_from_q(q_DM)
	gf = 1.
	if(geo_factor):
		gf = (q_star<1.)*q_star**(0.5*(2.-gamma_DM))+(q_star>1.)*q_star**(-0.5*(2.-gamma_DM))
	if not finite:
		return gf*((q_star>1.)*q_DM**(2.-gamma_DM)+(q_star<1.)*1.)*geometric_factor(q_DM,.5*gamma_DM)*np.sqrt(Jkin(q_eff,gamma_star,axis='major'))
	else:
		return gf*((q_star>1.)*q_DM**(2.-gamma_DM)+(q_star<1.)*1.)*np.sqrt(Jkin(q_eff,gamma_star,axis='major'))

def round_factor_D(q_star,q_DM,gamma_star,gamma_DM):
	''' Combination of kinematic and geometric factors for D-factor face-on'''
	q_eff = q_star/qpot_from_q(q_DM)
	return np.sqrt(Jkin(q_eff,gamma_star,axis='minor'))

## ============================================================================
## Plotting functions
def both_round_pl(qrange,name='both_round.pdf',gamma_star = 2.,Dfactor=False):
	''' Generate Fig 7 '''
	if(Dfactor):
		fn = round_factor_D
	else:
		fn = round_factor

	## Plot M2M results
	datt = np.genfromtxt('ret_results.dat')
	plt.plot(datt.T[0],datt.T[2+2*Dfactor],'.',color='k',ms=5,zorder=10)

	## Fit curves with polynomial
	p, pcov = curve_fit(func_qm, datt.T[0][datt.T[0]<1.],datt.T[2+2*Dfactor][datt.T[0]<1.])
	p2, pcov = curve_fit(func_qp, datt.T[0][datt.T[0]>1.],datt.T[2+2*Dfactor][datt.T[0]>1.])
	ff = np.concatenate((func_qm(qrange[qrange<1.],p[0],p[1]),
	                    func_qp(qrange[qrange>1.],p2[0],p2[1])))
	l,=plt.plot(qrange,ff,label=r'$q^\eta$ fit',color='k')

	## Plot numerical results
	datt2 = np.genfromtxt("../../../code/jfactors/data/fig7_gf.cf.dat")
	datt2 = datt2[np.argsort(datt2.T[0])]
	l,=plt.plot(datt2.T[0],np.log10(datt2.T[1+Dfactor]),label='Virial Method',color=sns.color_palette()[1])
	l.set_dashes((2,1))

	gamma_DM = 1.

	## Plot single cusp model
	plt.fill_between(qrange,
	                 np.log10(fn(qrange,qrange,gamma_star-1,gamma_DM)),
	                 np.log10(fn(qrange,qrange,gamma_star+1,gamma_DM)),
	                 edgecolor="None",alpha=0.5,
	                 color=sns.color_palette()[0])
	loc = (3.,-0.72)
	rot = -16
	if(Dfactor):
		loc = (2.5,-0.1)
		rot = -6
	plt.annotate(r'$\gamma_\star=%i$'%(gamma_star-1),xy=loc,rotation=rot,color=sns.color_palette()[0])
	loc = (2.,-0.77)
	rot = -24
	if(Dfactor):
		loc = (2.5,-0.24)
		rot = -10
	plt.annotate(r'$\gamma_\star=%i$'%(gamma_star+1),xy=loc,rotation=rot,color=sns.color_palette()[0])
	l,=plt.plot(qrange,np.log10(fn(qrange,qrange,gamma_star,1.)),label=r'Cusp, $\gamma_\star=%i$'%gamma_star,color=sns.color_palette()[0])
	l.set_dashes((3,1))

	wyn = load_wyn_data()
	index = 1+Dfactor
	l,=plt.plot(wyn[0].T[0],np.log10(wyn[1].T[index]),label=r'Core, $\beta_\star=5$, $R_d/R_c=20$',color=sns.color_palette()[2])
	l.set_dashes((5,1))
	plt.fill_between(wyn[0].T[0],np.log10(wyn[0].T[index]),np.log10(wyn[2].T[index]),edgecolor="None",alpha=0.3,color=sns.color_palette()[2])

	loc = (2.0,-0.49)
	rot = -16
	if(Dfactor):
		loc = (3.2,-0.12)
		rot = -5
	plt.annotate(r'$R_d/R_c=2$',xy=loc,rotation=rot,color=sns.color_palette()[2])
	loc = (1.,-0.23)
	rot = -44
	if(Dfactor):
		loc = (1.2,-0.15)
		rot = -35
	plt.annotate(r'$R_d/R_c=200$',xy=loc,rotation=rot,color=sns.color_palette()[2])

	l=plt.axvline(1.,color='k',ls='dashed')
	l.set_dashes((2,1))

	## Plot little ellipses to show the viewing angle
	loc = (0.5,-0.5)
	arrow_length=0.1
	width, height = 0.5, 0.15
	hw = 0.05
	hl = 1.5*hw
	if(Dfactor):
		loc=(0.5,-0.25)
		height/=2
		arrow_length/=2
		hl/=2
	xloc = loc[1]+height/2.+arrow_length+hl*1.3
	ell = Ellipse(xy=loc, width=width, height=height, angle=0.,color='k')
	plt.arrow(loc[0],xloc,0.,-arrow_length,color='r',head_width=hw,head_length=hl)
	plt.gca().add_artist(ell)
	loc = (2.5,0.5)
	if(Dfactor):
		loc = (2.5,0.25)
		width,height=0.25,0.15
	ell = Ellipse(xy=loc, width=height, height=width, angle=0.,color='k')
	xloc = loc[1]+width/2.+arrow_length+hl*1.3
	plt.arrow(loc[0],xloc,0.,-arrow_length,color='r',head_width=hw,head_length=hl)
	plt.gca().add_artist(ell)

	plt.ylim(-1.,1.)
	if(Dfactor):
		plt.ylim(-0.5,0.5)
	plt.xlabel(r'$q$')
	plt.ylabel(r'$\mathcal{F}_\mathrm{J}$')
	if(Dfactor):
		plt.ylabel(r'$\mathcal{F}_\mathrm{D}$')
	plt.annotate(r'Viewed face-on', xy=(0.,1.01), xycoords='axes fraction',fontsize=8,
	 		     horizontalalignment='left', verticalalignment='bottom')
	plt.legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.05),handlelength=2.)

	plt.savefig(name,bbox_inches='tight')

	return [p[0],p2[0]]

def both_edge_pl(qrange,name='both_edge.pdf',gamma_star = 2.,Dfactor=False,geo_factor=False):
	''' Generate Fig 6 '''
	if(Dfactor):
		fn = edge_factor_D
	else:
		fn = edge_factor

	## M2M results
	datt = np.genfromtxt('ret_results.dat')
	plt.plot(datt.T[0],datt.T[1+2*Dfactor],'.',color='k',ms=5,zorder=10)

	## Fit polynomials
	p, pcov = curve_fit(func_qm, datt.T[0][datt.T[0]<1.],datt.T[1+2*Dfactor][datt.T[0]<1.])
	p2, pcov = curve_fit(func_qp, datt.T[0][datt.T[0]>1.],datt.T[1+2*Dfactor][datt.T[0]>1.])
	ff = np.concatenate((func_qm(qrange[qrange<1.],p[0],p[1]),
	                    func_qp(qrange[qrange>1.],p2[0],p2[1])))
	l,=plt.plot(qrange,ff,label=r'$q^\eta$ fit',color='k')

	## Plot numerical result
	datt2 = np.genfromtxt("../../../code/jfactors/data/fig7_gf.cf.dat")
	datt2 = datt2[np.argsort(datt2.T[0])]
	l,=plt.plot(datt2.T[0],np.log10(datt2.T[3+Dfactor]),label='Virial Method',color=sns.color_palette()[1])
	l.set_dashes((2,1))

	gamma_DM = 1.

	## Plot axisymmetric cusp model
	plt.fill_between(qrange,
	                 np.log10(fn(qrange,qrange,gamma_star,gamma_DM-.5,geo_factor)),
	                 np.log10(fn(qrange,qrange,gamma_star,gamma_DM+.5,geo_factor)),
	                 edgecolor="None",alpha=0.3,color=sns.color_palette()[0])
	loc = (2.28,0.55)
	rot = 12
	if(Dfactor):
		loc = (1.5,0.35)
		rot = 13
	plt.annotate(r'$\gamma_\mathrm{DM}=%0.1f$'%(gamma_DM-0.5),xy=loc,rotation=rot,color=sns.color_palette()[0])
	loc = (2.,-0.07)
	rot = 1
	if(Dfactor):
		loc = (2.3,0.25)
		rot = 5
	plt.annotate(r'$\gamma_\mathrm{DM}=%0.1f$'%(gamma_DM+0.5),xy=loc,rotation=rot,color=sns.color_palette()[0])

	l,=plt.plot(qrange,np.log10(fn(qrange,qrange,gamma_star,gamma_DM,geo_factor)),label=r'Cusp, $\gamma_{\rm DM}=1,\gamma_\star=%i$'%gamma_star,color=sns.color_palette()[0])
	l.set_dashes((3,1))

	wyn = load_wyn_data()
	index = 3+Dfactor
	l,=plt.plot(wyn[0].T[0],np.log10(wyn[1].T[index]),label=r'Core, $\beta_\star=5$, $R_d/R_c=20$',color=sns.color_palette()[2])
	l.set_dashes((5,1))
	plt.fill_between(wyn[0].T[0],np.log10(wyn[0].T[index]),np.log10(wyn[2].T[index]),edgecolor="None",alpha=0.3,color=sns.color_palette()[2])

	loc = (0.33,-0.16)
	rot = 25
	if(Dfactor):
		loc = (3.,-0.07)
		rot = 2
	plt.annotate(r'$R_d/R_c=2$',xy=loc,rotation=rot,color=sns.color_palette()[2])
	loc = (3.,0.21)
	rot = 6
	if(Dfactor):
		loc = (3.1,0.14)
		rot = 2
	plt.annotate(r'$R_d/R_c=200$',xy=loc,rotation=rot,color=sns.color_palette()[2])

	## Divide between prolate and oblate
	l=plt.axvline(1.,color='k',ls='dashed')
	l.set_dashes((2,1))

	## Plot little ellipses to show the viewing angle
	loc = (0.5,-0.75)
	arrow_length=0.1
	width, height = 0.5, 0.15
	hw = 0.05
	xloc = loc[0]+width/2.+arrow_length+hw*1.8
	ell = Ellipse(xy=loc, width=width, height=height, angle=0.,color='k')
	plt.arrow(xloc,loc[1],-arrow_length,0.,color='r',head_width=hw)
	plt.gca().add_artist(ell)
	loc = (3.,-0.22)
	if(Dfactor):
		loc = (2.5,0.7)
	ell = Ellipse(xy=loc, width=height, height=width, angle=0.,color='k')
	xloc = loc[0]+height/2.+arrow_length+hw*1.8
	plt.arrow(xloc,loc[1],-arrow_length,0.,color='r',head_width=hw)
	plt.gca().add_artist(ell)

	plt.ylim(-1.,1.)
	plt.xlabel(r'$q$')
	plt.ylabel(r'$\mathcal{F}_\mathrm{J}$')
	if(Dfactor):
		plt.ylabel(r'$\mathcal{F}_\mathrm{D}$')
	plt.legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.05),handlelength=2.)
	plt.annotate(r'Viewed edge-on', xy=(0.,1.01), xycoords='axes fraction',fontsize=8,
	 		     horizontalalignment='left', verticalalignment='bottom')

	plt.savefig(name,bbox_inches='tight')

	return [p[0],p2[0]]

def both_WR_pl(qrange,name='both_WR.pdf'):
	''' Generate Fig 5 '''
	plt.fill_between(qrange,np.log10(WR(qrange,qrange,2.)),np.log10(WR(qrange,qrange,4.)),edgecolor="None",alpha=0.3)
	plt.annotate(r'$\gamma_\star=2$',xy=(2.5,-0.17),rotation=-10)
	plt.annotate(r'$\gamma_\star=4$',xy=(2.5,-0.45),rotation=-17)
	l,=plt.plot(qrange,np.log10(WR(qrange,qrange,3.)),label=r'$\gamma_\star=3$')
	l.set_dashes((3,1))
	l=plt.axvline(1.,color='k',ls='dashed')

	plt.xlabel(r'$q$')
	plt.ylabel(r'$\log_{10}(\langle\sigma_{xx}^2\rangle/\langle\sigma_{zz}^2\rangle)$')

	datt2 = np.genfromtxt("../../../code/jfactors/data/fig6.wr.dat")
	datt2 = datt2[np.argsort(datt2.T[0])]
	l,=plt.plot(datt2.T[0],np.log10(datt2.T[1]),label='Tensor virial')
	l.set_dashes((2,1))

	vfunc = np.vectorize(binney_tremaine_virial_ratio)
	l,=plt.plot(qrange,np.log10(vfunc(qrange)),label='Spheroidal')
	l.set_dashes((4,1))

	datt = np.genfromtxt('ret_results.dat')
	plt.plot(datt.T[0],np.log10(datt.T[7]),'.',color='k',ms=5,lw=1)
	plt.legend(loc="lower center",ncol=3, bbox_to_anchor=(0.5, 1.0),handlelength=2.)
	plt.ylim(-0.6,0.6)
	plt.xlim(0.,4.)
	plt.savefig(name,bbox_inches='tight')


def both_round_GF_pl(qrange,name='both_round_GF.pdf'):
	''' Geometric factors for face-on case -- not in paper'''
	l,=plt.plot(qrange,np.log10(round_GF(qrange,0.)),label=r'$\gamma_{\rm DM}=0$')
	l,=plt.plot(qrange,np.log10(round_GF(qrange,0.5)),label=r'$\gamma_{\rm DM}=0.5$')
	l.set_dashes((1,1))
	l,=plt.plot(qrange,np.log10(round_GF(qrange,1.)),label=r'$\gamma_{\rm DM}=1$')
	l.set_dashes((3,1))
	l,=plt.plot(qrange,np.log10(round_GF(qrange,2.)),label=r'$\gamma_{\rm DM}=2$')
	l.set_dashes((5,2))
	l,=plt.plot(qrange,np.log10(round_GF(qrange,3.)),label=r'$\gamma_{\rm DM}=3$')
	l.set_dashes((2,2))
	l,=plt.plot(qrange,np.log10(round_GF(qrange,4.)),label=r'$\gamma_{\rm DM}=4$')
	l.set_dashes((7,2))

	plt.xlabel(r'$q$')
	plt.ylabel(r'$\log_{10}\mathrm{J}_{\rm geo,face}$')
	plt.annotate('Viewed face-on', xy=(0.,1.01), xycoords='axes fraction',fontsize=8,
	 		     horizontalalignment='left', verticalalignment='bottom')
	plt.legend(loc="lower center",ncol=3, bbox_to_anchor=(0.5, 1.05),handlelength=2.)
	plt.ylim(-0.4,0.4)
	plt.savefig(name,bbox_inches='tight')

def both_edge_GF_pl(qrange,name='both_edge_GF.pdf'):
	''' Geometric factors for edge-on case -- not in paper'''

	l,=plt.plot(qrange,np.log10(edge_GF(qrange,0.)),label=r'$\gamma_{\rm DM}=0$')
	l,=plt.plot(qrange,np.log10(edge_GF(qrange,0.5)),label=r'$\gamma_{\rm DM}=0.5$')
	l.set_dashes((1,1))
	l,=plt.plot(qrange,np.log10(edge_GF(qrange,1.)),label=r'$\gamma_{\rm DM}=1$')
	l.set_dashes((3,1))
	l,=plt.plot(qrange,np.log10(edge_GF(qrange,2.)),label=r'$\gamma_{\rm DM}=2$')
	l.set_dashes((5,2))
	l,=plt.plot(qrange,np.log10(edge_GF(qrange,3.)),label=r'$\gamma_{\rm DM}=3$')
	l.set_dashes((2,2))
	l,=plt.plot(qrange,np.log10(edge_GF(qrange,4.)),label=r'$\gamma_{\rm DM}=4$')
	l.set_dashes((7,2))

	plt.xlabel(r'$q$')
	plt.ylabel(r'$\log_{10}\mathrm{J}_{\rm geo,edge}$')
	plt.annotate('Viewed edge-on', xy=(0.,1.01), xycoords='axes fraction',fontsize=8,
	 		     horizontalalignment='left', verticalalignment='bottom')
	plt.legend(loc="lower center",ncol=3, bbox_to_anchor=(0.5, 1.05),handlelength=2.)
	plt.ylim(-0.4,0.4)
	plt.savefig(name,bbox_inches='tight')

## ============================================================================
## Generate plots and also fit the correction factors with polynomial functions
## of q and output in fit_results.dat
if __name__ == '__main__':

	print 'Diff between eq. (20) and numerical (q=1.4,gamma=3):',check_equation_20(1.4,3.)
	print 'Diff between eq. (20) and numerical (q=0.4,gamma=3):',check_equation_20(0.4,3.)
	print 'Diff between eq. (20) and numerical (q=0.4,gamma=4):',check_equation_20(0.4,4.)
	print 'Diff between eq. (20) and numerical (q=1.4,gamma=4):',check_equation_20(1.4,4.)
	print 'Diff between eq. (20) and numerical (q=1.1,gamma=3):',check_equation_20(1.1,3.)
	print 'Diff between eq. (20) and numerical (q=0.9,gamma=3):',check_equation_20(0.9,3.)

	name = 'ellcor'

	qrange=np.linspace(0.02,4.,100)
	g = open('fit_results_ascii.dat','w')
	with open('fit_results.dat','w') as f:
		f.write('\\begin{tabular}{lccc}\n')
		f.write('\\hline\n\\hline\n')
		f.write('& View & $\eta$ Oblate ($q<1$) & $\eta$ Prolate ($q>1$)\\\\\n')
		g.write('q<1 q>1\n')
		f.write('\\hline\n')
		X = both_edge_pl(qrange,name='both_edge_'+name+'.pdf',gamma_star=3.,geo_factor=True);plt.clf()
		g.write("%0.3f %0.3f\n"%(X[0],X[1]))
		f.write("$\mathcal{F}_\mathrm{J}$ & Edge-on&%0.3f&%0.3f\\\\\n"%(X[0],X[1]))
		X = both_round_pl(qrange,name='both_round_'+name+'.pdf',gamma_star=3.);plt.clf()
		g.write("%0.3f %0.3f\n"%(X[0],X[1]))
		f.write("&Face-on&%0.3f&%0.3f\\\\\n"%(X[0],X[1]))
		X = both_edge_pl(qrange,name='both_edge_D_'+name+'.pdf',gamma_star=3.,Dfactor=True,geo_factor=True);plt.clf()
		g.write("%0.3f %0.3f\n"%(X[0],X[1]))
		f.write("$\mathcal{F}_\mathrm{D}$ & Edge-on&%0.3f&%0.3f\\\\\n"%(X[0],X[1]))
		X = both_round_pl(qrange,name='both_round_D_'+name+'.pdf',gamma_star=3.,Dfactor=True);plt.clf()
		g.write("%0.3f %0.3f\n"%(X[0],X[1]))
		f.write("&Face-on&%0.3f&%0.3f\\\\\n"%(X[0],X[1]))
		f.write('\\hline\n')
		f.write('\\end{tabular}\n')
	g.close()
	both_WR_pl(qrange,name='both_WR_'+name+'.pdf');plt.clf()
	both_edge_GF_pl(qrange);plt.clf()
	both_round_GF_pl(qrange);plt.clf()

## ============================================================================
