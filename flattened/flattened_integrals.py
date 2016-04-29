import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.special import ellipk, ellipe
from scipy.integrate import quad,fixed_quad,nquad
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import Ellipse
from spherical_Jfactors import asymmetric_gaussian_samples

def integrate_totalmass_alphabetagamma(rho0,rs,alpha,beta,gamma,ba,ca,rt):
	def rho(x,y,z):
		a = 1./(ba*ca)**(1./3.)
		b = ba*a
		c = ca*a
		r = np.sqrt(x**2/a**2+y**2/b**2+z**2/c**2)
		fac = 1.
		if(rt>0.):
			fac=np.sqrt(1-np.tanh(r/rt)**2)
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*fac
	return rho0*quad(lambda x:quad(lambda y: quad(lambda z: rho(x,y,z), -np.inf, np.inf)[0], -np.inf, np.inf)[0],-np.inf,np.inf)[0]

def integrate_totalmass_alphabetagamma_sph(rho0,rs,alpha,beta,gamma,ba,ca,rt):
	def rho(r,p,t):
		x = r*np.sin(t)*np.cos(p)
		y = r*np.sin(t)*np.sin(p)
		z = r*np.cos(t)
		a = 1./(ba*ca)**(1./3.)
		b = ba*a
		c = ca*a
		r = np.sqrt(x**2/a**2+y**2/b**2+z**2/c**2)
		fac = 1.
		if(rt>0.):
			fac=np.sqrt(1-np.tanh(r/rt)**2)
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*fac
	return rho0*quad(lambda r:quad(lambda p: quad(lambda t: r*r*rho(r,p,t)*np.sin(t), 0., np.pi)[0], 0., 2.*np.pi)[0],-0.,np.inf)[0]

def integrate_Jthetamax_alphabetagamma(thetamax,D,rho0,rs,alpha,beta,gamma,ba,ca,rt,los='z'):
	def rho(x,y,z):
		a = 1./(ba*ca)**(1./3.)
		b = ba*a
		c = ca*a
		r = np.sqrt(x**2/a**2+y**2/b**2+z**2/c**2)
		fac = 1.
		if(rt>0.):
			fac=np.sqrt(1-np.tanh(r/rt)**2)
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*fac
	def J(ll,phi,theta):
		b=D*np.tan(theta)
		z = ll
		x = np.cos(phi)*b
		y = np.sin(phi)*b
		if(los=='z'):
			return np.sin(theta)*rho(x,y,z)**2
		elif(los=='y'):
			return np.sin(theta)*rho(x,z,y)**2
		elif(los=='x'):
			return np.sin(theta)*rho(z,x,y)**2
	return rho0*rho0*quad(lambda x:quad(lambda y: quad(lambda z: J(x,y,z), 0., thetamax)[0], 0., 2.*np.pi)[0],-np.inf,np.inf)[0]

def integrate_Dthetamax_alphabetagamma(thetamax,D,rho0,rs,alpha,beta,gamma,ba,ca,rt,los='z'):
	def rho(x,y,z):
		a = 1./(ba*ca)**(1./3.)
		b = ba*a
		c = ca*a
		r = np.sqrt(x**2/a**2+y**2/b**2+z**2/c**2)
		fac = 1.
		if(rt>0.):
			fac=np.sqrt(1-np.tanh(r/rt)**2)
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*fac
	def J(ll,phi,theta):
		b=D*np.tan(theta)
		z = ll
		x = np.cos(phi)*b
		y = np.sin(phi)*b
		if(los=='z'):
			return np.sin(theta)*rho(x,y,z)
		elif(los=='y'):
			return np.sin(theta)*rho(x,z,y)
		elif(los=='x'):
			return np.sin(theta)*rho(z,x,y)
	return rho0*quad(lambda x:quad(lambda y: quad(lambda z: J(x,y,z), 0., thetamax)[0], 0., 2.*np.pi)[0],-np.inf,np.inf)[0]

def integrate_rho_alphabetagamma(thetamax,D,rho0,rs,alpha,beta,gamma,ba,ca,rt):
	## Mass in cylinder
	def rho(x,y,z):
		a = 1./(ba*ca)**(1./3.)
		b = ba*a
		c = ca*a
		r = np.sqrt(x**2/a**2+y**2/b**2+z**2/c**2)
		fac = 1.
		if(rt>0.):
			fac=np.sqrt(1-np.tanh(r/rt)**2)
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*fac
	def J(ll,phi,theta):
		b=D*np.tan(theta)
		z = ll
		x = np.cos(phi)*b
		y = np.sin(phi)*b
		return np.sin(theta)*rho(x,y,z)
	return rho0*quad(lambda x:quad(lambda y: quad(lambda z: J(x,y,z), 0., thetamax)[0], 0., 2.*np.pi)[0],-np.inf,np.inf)[0]

def integrate_rho_ellp_alphabetagamma(r,rho0,rs,alpha,beta,gamma,ba,ca,rt):
	def rho(r):
		fac = 1.
		if(rt>0.):
			fac=np.sqrt(1-np.tanh(r/rt)**2)
		return np.power(r/rs,-gamma)*np.power(1+np.power(r/rs,alpha),((gamma-beta)/alpha))*fac
	a = 1./(ba*ca)**(1./3.)
	return 4.*np.pi*rho0*quad(lambda x:x*x*rho(x/a),0.,r)[0]
