## Checking triaxial cusps SEG(2016)
## ============================================================================
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib.patches import Ellipse
import flattened as fJ
from scipy.optimize import curve_fit
import seaborn as sns
## ============================================================================
def geometric_factor(q,gamma):
	''' Geometric factor for infinite axisymmetric cusp '''
	## \int dt (cos^2(t)+sin^2(t)/q^2)^{1/2-gamma}
	return quad(lambda t: np.power(np.power(np.cos(t),2.)+np.power(np.sin(t),2.)/q/q,0.5-gamma),0.,2.*np.pi)[0]
def f1(p,q,gamma):
	''' Virial ratio for infinite triaxial cusp '''
	return quad(lambda phi:quad(lambda t: np.cos(phi)**2*np.sin(t)**3*(np.sin(t)**2*np.cos(phi)**2+np.sin(t)**2*np.sin(phi)**2/p/p+np.cos(t)**2./q/q)**(-gamma/2.),0.,np.pi)[0],0.,2.*np.pi)[0]/quad(lambda phi:quad(lambda t: np.sin(t)*np.cos(t)**2*(np.sin(t)**2*np.cos(phi)**2+np.sin(t)**2*np.sin(phi)**2/p/p+np.cos(t)**2./q/q)**(-gamma/2.),0.,np.pi)[0],0.,2.*np.pi)[0]
def f2(p,q,gamma):
	''' Virial ratio for infinite triaxial cusp '''
	return quad(lambda phi:quad(lambda t: np.sin(phi)**2*np.sin(t)**3*(np.sin(t)**2*np.cos(phi)**2+np.sin(t)**2*np.sin(phi)**2/p/p+np.cos(t)**2./q/q)**(-gamma/2.),0.,np.pi)[0],0.,2.*np.pi)[0]/quad(lambda phi:quad(lambda t: np.sin(t)*np.cos(t)**2*(np.sin(t)**2*np.cos(phi)**2+np.sin(t)**2*np.sin(phi)**2/p/p+np.cos(t)**2./q/q)**(-gamma/2.),0.,np.pi)[0],0.,2.*np.pi)[0]
def jkin(p,q,gamma,th,ph):
	''' Kinematic factor for infinite triaxial cusp '''
	P = p/fJ.qpot_from_q(p)
	Q = q/fJ.qpot_from_q(q)
	ff1 = f1(P,Q,gamma)
	ff2 = f2(P,Q,gamma)
	return ((1.+ff1+ff2)/(np.cos(th)**2+ff1*np.sin(th)**2*np.cos(ph)**2+ff2*np.sin(th)**2*np.sin(ph)**2)/3.)**2

## ============================================================================

def jgeo_x(p,q,gamma):
	return p/q/q*geometric_factor(q/p,gamma)
def jgeo_y(p,q,gamma):
	return 1./p/q/q*geometric_factor(q,gamma)
def jgeo_z(p,q,gamma):
	return 1./p/p/q*geometric_factor(p,gamma)

def jkin_x(p,q,gamma):
	return jkin(p,q,gamma,.5*np.pi,0.)
def jkin_y(p,q,gamma):
	return jkin(p,q,gamma,.5*np.pi,.5*np.pi)
def jkin_z(p,q,gamma):
	return jkin(p,q,gamma,0.,0.)

def jtot_x(p,q,gammaDM,gammaST):
	return jgeo_x(p,q,gammaDM)*jkin_x(p,q,gammaST)
def jtot_y(p,q,gammaDM,gammaST):
	return jgeo_y(p,q,gammaDM)*jkin_y(p,q,gammaST)
def jtot_z(p,q,gammaDM,gammaST):
	return jgeo_z(p,q,gammaDM)*jkin_z(p,q,gammaST)

if __name__ == '__main__':
	q = 0.7
	p = 0.8
	gg = np.linspace(0.,5.,10)
	ggst = 3.
	qq = np.linspace(0.1,p,10)
	# plt.plot(gg,map(lambda g:jgeo_x(p,q,g),gg))
	# plt.plot(gg,map(lambda g:jgeo_y(p,q,g),gg))
	# plt.plot(gg,map(lambda g:jgeo_z(p,q,g),gg))
	plt.plot(qq,map(lambda g:jgeo_x(p,g,1.),qq))
	plt.plot(qq,map(lambda g:jgeo_y(p,g,1.),qq))
	plt.plot(qq,map(lambda g:jgeo_z(p,g,1.),qq))
	# plt.plot(gg,map(lambda g:jkin_x(p,g,ggst),qq))
	# plt.plot(gg,map(lambda g:jkin_y(p,g,ggst),qq))
	# plt.plot(gg,map(lambda g:jkin_z(p,g,ggst),qq))
	# plt.plot(gg,map(lambda g:jtot_x(p,q,g,ggst),gg))
	# plt.plot(gg,map(lambda g:jtot_y(p,q,g,ggst),gg))
	# plt.plot(gg,map(lambda g:jtot_z(p,q,g,ggst),gg))
	plt.savefig('tmp.pdf',bbox_inches='tight')
