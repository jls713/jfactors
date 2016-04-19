## Generate Fig 1 of Sanders, Evans & Geringer-Sameth
## ============================================================================
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/jls/work/code/m2m/')
from m2m_tools import *
import pynbody
from scipy.optimize import fmin
from subprocess import call
import ret_2
## ============================================================================
## Load in properties of Reticulum II
## ============================================================================
RetII = ret_2.RetII
r_maj_exp = RetII.r_maj_exp ## arcmin
r_maj = RetII.r_maj #arcmin
Distance = RetII.Distance #kpc
Velocity_dispersion = RetII.Velocity_dispersion #km/s
r0_unit = (r_maj/60./180.*np.pi)*Distance ## in units kpc
MtoL = 500.
Width = 70.
G = 4.300918e-6 ## in units solar mass, km/s kpc
b_over_a = 0.4
## ============================================================================
## Fitting functions for projected profiles
def fexp(p,R,m):
	''' Fit exponential '''
	if(p[0]<0.):
		return 1e10
	dr = np.log(np.exp(-R/p[0])*m/p[0]/p[0])
	return -np.sum(dr)
def fplum(p,R,m):
	''' Fit Plummer '''
	if(p[0]<0.):
		return 1e10
	dr = np.log(m*R/p[0]/p[0]/(1+(R/p[0])**2)**2)
	return -np.sum(dr)
## ============================================================================
def extract_final_gadget(input_file):
	''' Extract the final snapshot from.snp file, export to gadget and load into pynbody format '''
	output_gadget=input_file+".gadget"
	csh_command='$FALCON/bin/./s2g in='+input_file+' out='+output_gadget+' times=last;'
	run_csh_command(csh_command)

	snap = load_snp_pynbody_snapshot(output_gadget)

	return snap

def import_data(data_file):
	''' Load data file -- extract final snapshot to pynbody and add extras '''
	snap = extract_final_gadget(data_file)
	snap['pos'].convert_units('kpc')

	extract_final_snapshot(data_file,data_file+".tmp")
	process_density(data_file+".tmp",data_file+".density")
	data = read_dehnen_density_file(data_file+".density")
	call(['rm',data_file+".density",data_file+'.tmp'])
	data2=add_to_array(data)
	return data2
## ============================================================================
## Make Figure 1
def make_plot():
	f,a=plt.subplots(2,2,figsize=[3.32,3.32])
	plt.subplots_adjust(wspace=0.,hspace=0.)
	a[1][0].semilogy()
	a[1][0].semilogx()
	a[1][0].set_xlabel(r'$R/\mathrm{arcmin}$')
	a[1][0].set_ylabel(r'$M_\odot/ \mathrm{arcmin}^{-2}$')
	a[1][1].yaxis.tick_right()
	a[1][1].yaxis.set_label_position("right")
	a[1][1].set_xlabel(r'$v_\mathrm{los}/{\rm km\,s}^{-1}$')
	a[1][1].set_ylabel(r'$p(v_\mathrm{los}) /({\rm km\,s}^{-1})^{-1}$')
	a[1][1].set_xlim(-15.,15.)
	a[0][1].set_xlabel(r'$x/\mathrm{arcmin}$')
	a[0][1].yaxis.tick_right()
	a[0][1].yaxis.set_label_position("right")
	a[0][1].xaxis.tick_top()
	a[0][1].xaxis.set_label_position("top")
	bin_max = 5. ## Maximum radial bin in units of r_half
	bin_no = 30  ## Number of bins
	## =======================================================================
	## Oblate case
	data =import_data('/data/jls/m2m/ret2_flatNFW_flatPlummer_obl_0.4_M2M.snp')
	## 1. Measure rmaj from simulation
	data['R']=np.sqrt(data['x']**2+data['z']**2/0.4**2)
	rmajsim=fmin(fexp,[5.],args=(data['R'][data['R']<5.],data['m'][data['R']<5.]))[0] ## exponential
	rmajsim=fmin(fplum,[5.],args=(data['R'][data['R']<5.],data['m'][data['R']<5.]))[0] ## Plummer
	Mp = np.sum(data['m']) ## Total stellar mass

	## 2. Measure sigma_los from simulation
	sigma_unit = Velocity_dispersion # in units km/s
	sigma_los = np.sqrt(np.sum(data['m']*data['vy']**2)/np.sum(data['m']))

	GM = (sigma_unit/sigma_los)**2*(r0_unit/rmajsim)
	M = GM/G ## Total halo mass

	runit_am = (r0_unit/rmajsim)*60.*180./np.pi/Distance ## R unit in arcmin

	## 3. Plot projected profiles
	n,b,p=a[1][0].hist(data['R']*runit_am,
	                   weights=(M/MtoL)*data['m']/data['R']/runit_am/2./np.pi/(bin_max/bin_no*runit_am)/Mp,
	                   range=[0.,bin_max*runit_am],
	                   bins=bin_no,histtype='step',color='k')
	bc=.5*(b[1:]+b[:-1])
	bc = bc/(60.*180./np.pi/Distance)
	## 3.1. Plot best-fit
	l,=a[1][0].plot(bc*60.*180./np.pi/Distance,n[0]/(1+(bc/r0_unit)**2)**2,color=sns.color_palette()[0])

	## 4. Plot velocity dispersion profiles
	bins = np.linspace(-10.,10.,30)
	N =np.histogram(data['vy']*(sigma_unit/sigma_los),weights=data['m'],range=[-0.5*(sigma_unit/sigma_los),.5*(sigma_unit/sigma_los)],bins=bins,normed=True)
	l,=a[1][1].step(.5*(N[1][1:]+N[1][:-1]),N[0],color=sns.color_palette()[0],where='mid')

	print 'Size=',r0_unit,', VelDisp=',sigma_unit
	print 'From sim: Size=',rmajsim,', VelDisp=',sigma_los

	## 5. Plot projected 2D distributions
	snap.rotate_x(-90)
	snap['pos'].convert_units(str(1./(r0_unit/rmajsim*60.*180./np.pi/Distance))+' kpc')
	pl = pynbody.plot.image(snap.d,subplot=a[0][0],cmap="Greys",show_cbar=False,labelx=None,labely=r'$y/\mathrm{arcmin}$',av_z=True,width=Width)
	N = len(pl)
	W = Width
	x = np.array([W/N*(i+.5)-.5*W for i in range(N)])
	y = x
	a[0][0].contour(x,y,pl,extent=[x.min(),x.max(),y.min(),y.max()],linewidths=0.5,locator=LogLocator(base=2.),colors=[sns.color_palette()[0]])
	a[0][0].text(0.5,0.9,r'Oblate',horizontalalignment='center',verticalalignment='top',transform=a[0][0].transAxes,fontsize=14)

	## =======================================================================
	## Prolate case
	data =import_data('/data/jls/m2m/ret2_flatNFW_flatPlummer_pro_0.4_M2M.snp')

	## 1. Measure rmaj from simulation
	data['R']=np.sqrt(data['x']**2+data['y']**2/0.4**2)
	rmajsim=fmin(fexp,[5.],args=(data['R'][data['R']<5.],data['m'][data['R']<5.]))[0]
	rmajsim=fmin(fplum,[5.],args=(data['R'][data['R']<5.],data['m'][data['R']<5.]))[0]
	Mp = np.sum(data['m'])

	## 2. Measure sigma_los from simulation
	sigma_unit = Velocity_dispersion # in units km/s
	sigma_los = np.sqrt(np.sum(data['m']*data['vz']**2)/np.sum(data['m']))

	GM = (sigma_unit/sigma_los)**2*(r0_unit/rmajsim)
	M = GM/G ## Total halo mass

	runit_am = (r0_unit/rmajsim)*60.*180./np.pi/Distance ## R unit in arcmin

	## 3. Plot projected profiles
	n,b,p=a[1][0].hist(data['R']*runit_am,weights=(M/MtoL)*data['m']/data['R']/runit_am/2./np.pi/(bin_max/bin_no*runit_am)/Mp,range=[0.,bin_max*runit_am],bins=bin_no,histtype='step',color='k')
	bc=.5*(b[1:]+b[:-1])
	bc = bc/(60.*180./np.pi/Distance)
	## 3.1. Plot best-fit
	l, =a[1][0].plot(bc*60.*180./np.pi/Distance,n[0]/(1+(bc/r0_unit)**2)**2,color=sns.color_palette()[1],ls='dashed')
	l.set_dashes((3,1))

	## 4. Plot velocity dispersion profiles
	N =np.histogram(data['vz']*(sigma_unit/sigma_los),weights=data['m'],range=[-0.5*(sigma_unit/sigma_los),.5*(sigma_unit/sigma_los)],bins=bins,normed=True)
	l,=a[1][1].step(.5*(N[1][1:]+N[1][:-1]),N[0],color=sns.color_palette()[1],ls='dashed',where='mid')
	l.set_dashes((3,1))
	print 'Size=',r0_unit,', VelDisp=',sigma_unit
	print 'From sim: Size=',rmajsim,', VelDisp=',sigma_los

	## 5. Plot projected 2D distributions
	snap2['pos'].convert_units(str(1./(r0_unit/rmajsim*60.*180./np.pi/Distance))+' kpc')
	pl = pynbody.plot.image(snap2.d, subplot=a[0][1],cmap="Greys",show_cbar=False,labelx=None,labely=r'$y/\mathrm{arcmin}$',av_z=True,width=Width)
	N = len(pl)
	W = Width
	x = np.array([W/N*(i+.5)-.5*W for i in range(N)])
	y = x
	CS = a[0][1].contour(x,y,pl,extent=[x.min(),x.max(),y.min(),y.max()],linewidths=0.5,locator=LogLocator(base=2.),colors=[sns.color_palette()[1]],linestyles='dashed')
	for c in CS.collections:
	    c.set_dashes([(0, (2.0, 1.0))])

	a[0][1].text(0.5,0.9,r'Prolate',horizontalalignment='center',verticalalignment='top',transform=a[0][1].transAxes,fontsize=14)
	plt.setp(a[1][1].get_xticklabels()[0],visible=False)

	plt.savefig('paper/FlattenedPaper/models.pdf',bbox_inches='tight')

if __name__ == '__main__':
	make_plot()

## ============================================================================
