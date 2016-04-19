## Generate Fig 8 of SEG (2016) -- the triaxial sphere
## ============================================================================
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys
sys.path.append('../../../code/jfactors/')
sys.path.append('../../../code/m2m/')
import jfactors_py as cJ
import seaborn as sns
from matplotlib import _cntr as cntr
from collections import namedtuple
## ============================================================================
## Ret 2 properties
## Need to repeat this here rather than importing from ret_2 because there
## is some problem with the 3D plotting otherwise...
DwarfProperties = namedtuple("DwarfProperties","r_maj_exp r_maj Distance Velocity_dispersion e")
RetII = DwarfProperties(3.37,3.37*1.67,30.,3.22,0.6)  ## in arcmin, kpc, km/s
r_maj_exp = RetII.r_maj_exp ## arcmin
r_maj = RetII.r_maj #arcmin
D = RetII.Distance #kpc
slos = RetII.Velocity_dispersion #km/s
rh= (r_maj/60./180.*np.pi)*D ## in units kpc
## ============================================================================
def unit_sphere(outfile,gf=True,ell=False):
	''' Compute stuff on a unit sphere '''
	''' If ell = True, computes observed ellipticity '''
	''' If ell = False, computes correction factors '''
	''' Here we have specialized to the triaxial case '''
	## =======================================================================
	## 1. Create maps of the correction factor and ellipticity over sphere
	## Using T = 0.55
	ba=0.733
	ca=0.4
	ang = 0.5

	(n, m) = (160,160)

	fac=1.
	# Meshing a unit sphere according to n, m
	phi = np.linspace(0, 2 * np.pi, num=n, endpoint=False)
	theta = np.linspace(np.pi * 1./(m+1), np.pi*(1-1./(m+1)), num=m, endpoint=False)
	phi, theta = np.meshgrid(phi, theta)
	phi, theta = phi.ravel(), theta.ravel()
	phi = np.append(phi, [0.]) # Adding the north pole...
	theta = np.append(theta, [0.])
	phi = np.append(phi, [0.]) # Adding the south pole...
	theta = np.append(theta, [np.pi])
	mesh_x, mesh_y = (theta*np.cos(phi), theta*np.sin(phi))
	triangles = mtri.Triangulation(mesh_x, mesh_y).triangles
	x, y, z = fac*np.sin(theta)*np.cos(phi), fac*np.sin(theta)*np.sin(phi), fac*np.cos(theta)

	# Defining a custom color scalar field
	sphM = cJ.PaperModel(1.,0.999,rh,slos,True)
	sph = sphM.J_factor(0.,0.,D,ang,False,False)[0]
	pM = cJ.PaperModel(ba,ca,rh,slos,True)

	## We have the option of colouring the sphere according to J-factor or
	## ellipticity -- if ell = True colour by ellipticity
	def fnn(M,Th,Ph):
		if(ell):
			return M.ellipticity(Th,Ph)
		else:
			sphl = sph
			if(gf):
				sm = cJ.PaperModel(1.,1.,rh*np.sqrt(M.ellipticity(Th,Ph)),slos,True)
				sphl = sm.J_factor(Th,Ph,D,ang,False,False)[0]
				print Th,Ph,sphl
			return M.J_factor(Th,Ph,D,ang,False,False)[0]/sphl

	## Create a mesh that only covers an octant of the sphere -- as triaxial
	phi2 = np.linspace(0, np.pi/2., num=n/4+1, endpoint=True)
	theta2 = np.linspace(np.pi * 1./(m+1), np.pi/2., num=m/2+1, endpoint=True)
	phi2, theta2 = np.meshgrid(phi2, theta2)
	phi2, theta2 = phi2.ravel(), theta2.ravel()
	phi2 = np.append(phi2, [0.]) # Adding the north pole...
	theta2 = np.append(theta2, [0.])
	vals = np.array(map(lambda t,p:fnn(pM,t,p),theta2,phi2))

	## Now stack the results to cover the whole sphere
	allvals= np.reshape(vals[:-1],(m/2+1,n/4+1))
	allvals = np.hstack((allvals,allvals[:,::-1][:,1:]))
	allvals = np.hstack((allvals,allvals[:,1:-1]))
	allvals = np.vstack((allvals,allvals[::-1,:][1:-1,:]))
	allvals = np.append(np.append(allvals,vals[-1]),vals[-1])
	allvals = allvals.ravel()
	## The colour is the average over the values on the triangle
	colors = np.mean(np.log10(allvals[triangles]), axis=1)
	if(ell):
		colors = np.mean(allvals[triangles], axis=1)

	## Create a map of the ellipticity on the sphere
	ellip = np.array(map(lambda t,p:pM.ellipticity(t,p),theta2,phi2))
	ellip_all= np.reshape(ellip[:-1],(m/2+1,n/4+1))
	ellip_all = np.hstack((ellip_all,ellip_all[:,::-1][:,1:]))
	ellip_all = np.hstack((ellip_all,ellip_all[:,1:-1]))
	ellip_all = np.vstack((ellip_all,ellip_all[::-1,:][1:-1,:]))
	ellip_all = np.append(np.append(ellip_all,ellip[-1]),ellip[-1])
	ellip_all = ellip_all.ravel()

	## =======================================================================
	## 2. Plot the colour map
	fig = plt.figure(figsize=[3.5,4.])
	ax = fig.gca(projection='3d')
	ax.set_aspect('equal')
	cmap = sns.cubehelix_palette(8,start=.5,rot=-.75,as_cmap=True)
	triang = mtri.Triangulation(x, y, triangles)
	collec = ax.plot_trisurf(triang, z, cmap=cmap, shade=False, linewidth=0.,zorder=0)

	collec.set_array(colors)
	collec.autoscale()
	lbl = r'$\mathcal{F}_\mathrm{J}$'
	if(ell):
		lbl = r'$1-e$'
	plt.colorbar(collec,orientation='horizontal',label=lbl)
	# ax = plt.gca()
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_zticklabels([])
	ax.dist=7.

	max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0

	mean_x = x.mean()
	mean_y = y.mean()
	mean_z = z.mean()
	ax.set_xlim(mean_x - max_range, mean_x + max_range)
	ax.set_ylim(mean_y - max_range, mean_y + max_range)
	ax.set_zlim(mean_z - max_range, mean_z + max_range)

	## =======================================================================
	## 3. Add axes
	#draw a vector
	from matplotlib.patches import FancyArrowPatch
	from mpl_toolkits.mplot3d import proj3d

	class Arrow3D(FancyArrowPatch):
	    def __init__(self, xs, ys, zs, *args, **kwargs):
	        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
	        self._verts3d = xs, ys, zs

	    def draw(self, renderer):
	        xs3d, ys3d, zs3d = self._verts3d
	        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
	        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
	        FancyArrowPatch.draw(self, renderer)

	a = Arrow3D([-.8,-.8],[-.85,-.85],[.85,1.2], mutation_scale=40, lw=1, arrowstyle="-|>", color="k")
	ax.add_artist(a)
	ax.text(-.83,-.95,1.125, r"$z$")
	a = Arrow3D([-.8,-.45],[-.85,-.85],[.85,.85], mutation_scale=40, lw=1, arrowstyle="-|>", color="k")
	ax.text(-.62,-.9,0.77, r"$x$")
	ax.add_artist(a)
	a = Arrow3D([-.8,-.8],[-.85,-.5],[.85,.85], mutation_scale=40, lw=1, arrowstyle="-|>", color="k")
	ax.text(-.83,-.62,0.9, r"$y$")
	ax.add_artist(a)

	## =======================================================================
	## 4. Make a small ellipsoidal inset that shows isodensity
	phis = np.linspace(0, 2 * np.pi, num=n/(n/40), endpoint=False)
	thetas = np.linspace(np.pi * 1./(m+1), np.pi*(1-1./(m+1)), num=m/(m/40), endpoint=False)
	phis, thetas = np.meshgrid(phis, thetas)
	phis, thetas = phis.ravel(), thetas.ravel()
	phis = np.append(phis, [0.]) # Adding the north pole...
	thetas = np.append(thetas, [0.])
	phis = np.append(phis, [0.]) # Adding the south pole...
	thetas = np.append(thetas, [np.pi])
	xs, ys, zs = fac*np.sin(thetas)*np.cos(phis), fac*np.sin(thetas)*np.sin(phis), fac*np.cos(thetas)

	mesh_xs, mesh_ys = (thetas*np.cos(phis)-.7, thetas*np.sin(phis)*ba-.7)
	triangles = mtri.Triangulation(mesh_xs, mesh_ys).triangles
	fac = 0.3
	triangs = mtri.Triangulation(fac*xs-.7, fac*ba*ys-.7, triangles)
	collec = ax.plot_trisurf(triangs, fac*ca*zs-.9, color=sns.color_palette()[0],shade=True, linewidth=0.1,zorder=1)

	ax.view_init(28,-62)
	ax.dist=7.

	## =======================================================================
	## 5. Plot contours of constant observed ellipticity -- note we need to
	##    mask out those behind the sphere
	ae = ax.azim,ax.elev
	ae = np.deg2rad(ae)
	ae[1] = np.pi/2.-ae[1]
	ae[0] = ae[0]-np.pi

	def dotpp(Th,Ph):
		r = np.array([np.sin(ae[1])*np.cos(ae[0]),np.sin(ae[1])*np.sin(ae[0]),np.cos(ae[1])])
		vv = np.array([np.sin(Th)*np.cos(Ph),np.sin(Th)*np.sin(Ph),np.cos(Th)])
		return np.dot(r,vv)

	phi = np.reshape(phi[:-2],(n,m))
	ellip_all = np.reshape(ellip_all[:-2],(n,m))
	theta = np.reshape(theta[:-2],(n,m))
	phi = phi[:,3*n/40:25*n/40]
	theta = theta[:,3*n/40:25*n/40]
	ellip_all = ellip_all[:,3*n/40:25*n/40]
	## add contours
	C=cntr.Cntr(phi,theta,ellip_all)
	contour_list = [0.5,0.6,0.7,0.8,0.9]

	def plot_contour(th,ph):
		doo = np.array(map(lambda t,p:dotpp(t,p),th,ph))
		th = th[doo>0.]
		ph = ph[doo>0.]
		fac = 1.
		x,y,z = fac*np.sin(th)*np.cos(ph-np.pi), fac*np.sin(th)*np.sin(ph-np.pi), fac*np.cos(th)
		ax.plot(x,y,z,color='k',zorder=1,lw=0.5)

	for cc in contour_list:
		res = C.trace(cc)
		nseg = len(res) // 2
		segs, codes = res[:nseg], res[nseg:]
		ss = segs[0]
		ss = res[0]
		for j in ss:
			th,ph=ss.T[1],ss.T[0]
			plot_contour(th,ph)
			plot_contour(np.pi-th,ph)
			plot_contour(th,np.pi-ph)
			plot_contour(np.pi-th,np.pi-ph)
			# plot_contour(th,ph+np.pi)
			plot_contour(np.pi-th,ph+np.pi)
			plot_contour(th,2.*np.pi-ph)
			plot_contour(np.pi-th,2.*np.pi-ph)
	th,ph = np.pi/2.-0.03,-np.pi/2.-0.3
	ax.text3D(np.sin(th)*np.cos(ph),np.sin(th)*np.cos(ph),np.cos(th),r'$e=0.5$',zdir=np.array([-np.sin(ph),np.cos(ph),0.1]),fontsize=6,color='k')
	th,ph = np.pi/2.-0.4,-np.pi/2.-0.3
	ax.text3D(np.sin(th)*np.cos(ph),np.sin(th)*np.cos(ph),np.cos(th),r'$e=0.4$',zdir=np.array([-np.sin(ph),np.cos(ph),0.]),fontsize=6,color='k')
	th,ph = np.pi/2.-0.68,-np.pi/2.-0.1
	ax.text3D(np.sin(th)*np.cos(ph),np.sin(th)*np.cos(ph),np.cos(th),r'$e=0.3$',zdir=np.array([-np.sin(ph),np.cos(ph),-0.65]),fontsize=6,color='k')
	th,ph = np.pi/2.-0.6,-np.pi/2.+0.22
	ax.text3D(np.sin(th)*np.cos(ph),np.sin(th)*np.cos(ph),np.cos(th),r'$e=0.2$',zdir=np.array([-np.sin(ph),np.cos(ph),-1.85]),fontsize=6,color='k')
	th,ph = np.pi/2.-0.5,-np.pi/2.+0.36
	ax.text3D(np.sin(th)*np.cos(ph),np.sin(th)*np.cos(ph),np.cos(th),r'$e=0.1$',zdir=np.array([-np.sin(ph),np.cos(ph),-1.8]),fontsize=6,color='k')

	plt.show()
	plt.savefig(outfile,bbox_inches='tight',dpi=1000)

## ===========================================================================

if __name__ == '__main__':
	unit_sphere('triax_figure_gf.png',gf=True,ell=False)
	unit_sphere('triax_figure_ell.png',ell=True)

## ===========================================================================
