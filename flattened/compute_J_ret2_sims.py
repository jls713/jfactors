## Generates the table of J-factors for the M2M points of Fig 6,7 of SEG16
## ============================================================================
import sys
sys.path.append('../../code/m2m/')
import ret_2
import spherical_Jfactors as sJ
import numpy as np
## ============================================================================
## Load Ret 2 properties
RetII = ret_2.RetII
rh= (RetII.r_maj/60./180.*np.pi)*RetII.Distance ## in units kpc

## ============================================================================
## Grid of models
qrange = [0.4,0.5,0.6,0.7]

models = [ret_2.stellar_halo_model('ret2_ba04cusp_NFW_0.5_0.4',r'NFW $(b/a)=0.5$ $(c/a)=0.4$, Plummer $\gamma=1$ $(b/a)=0.4$ $(c/a)=0.38$',baN=0.5,caN=0.4,alpha=2.,beta=5.,gamma=1.,rnfwrs=2.,rtrunc=10.,bovera=0.4)]

sphr_model = ret_2.stellar_halo_model('ret2_roundNFW_roundPlummer',r'Spherical NFW, Spherical Plummer',baN=1.,caN=1.,alpha=2.,beta=5.,gamma=0.,rnfwrs=2.,rtrunc=10.,bovera=1.)

pro_models = [ret_2.stellar_halo_model('ret2_flatNFW_flatPlummer_pro'+'_%0.1f'%(q),
	                            r'NFW $(b/a)=%0.1f$ $(c/a)=%0.1f$, Plummer $(b/a)=%0.1f$ $(c/a)=%0.1f$' % (q,q,q,q),
	                            baN=q,
	                            caN=q,
	                            alpha=2.,
	                            beta=5.,
	                            gamma=0.,
	                            rnfwrs=2.,
	                            rtrunc=10.,
	                            bovera=q) for q in [0.4,0.5,0.6,0.7]]

obl_models = [ret_2.stellar_halo_model('ret2_flatNFW_flatPlummer_obl'+'_%0.1f'%(q),
	                            r'NFW $(c/a)=%0.1f$, Plummer $(c/a)=%0.1f$' % (q,q),
	                            baN=1.,
	                            caN=q,
	                            alpha=2.,
	                            beta=5.,
	                            gamma=0.,
	                            rnfwrs=2.,
	                            rtrunc=10.,
	                            bovera=1.) for q in qrange]

extra_models = [ret_2.stellar_halo_model('ret2_roundNFW_roundPlummer',r'Spherical NFW, Spherical Plummer',baN=1.,caN=1.,alpha=2.,beta=5.,gamma=0.,rnfwrs=2.,rtrunc=10.,bovera=1.),]

## ============================================================================

def compute_table(filename='ret_results.dat',angle=0.5,geo_factor=True):
	''' Generate the table of points -- columns of table are
		q CF_edge CF_round CFD_edge CFD_round KF_edge KF_round WR rho0z rsz rtz rho0x rsx rtx '''
	## ========================================================================
	## 1. Spherical analytic formulae -- note the ellipticity correction
	##    sqrt(1-e)
	gf = np.sqrt(1.-RetII.e)
	sph_J = sJ.wyns_formulaJ_NFW_data(RetII.Velocity_dispersion,rh*1000.*gf,RetII.Distance,[0.5],rh*2.*gf)
	sph_D = sJ.wyns_formulaD_NFW_data(RetII.Velocity_dispersion,rh*1000.*gf,RetII.Distance,[0.5],rh*2.*gf)
	print sph_J,sph_D
	## ========================================================================
	## 2. Spherical model
	FF = ret_2.find_J_value(extra_models[0],angle=angle)
	sph_J = np.log10(FF[-2])
	sph_D = np.log10(FF[-1])
	print sph_J,sph_D
	## ========================================================================
	## 3. Spherical models with ellipticity correction sqrt(1-e)
	FF_q = np.zeros((len(qrange),2))
	for n,q in enumerate(qrange):
		if(geo_factor==True):
			FF_q[n]=ret_2.find_J_value(extra_models[0],angle=angle,gf=np.sqrt(q))[-2:]
		else:
			FF_q[n]=FF[-2:]
	## ========================================================================
	## 4. Compute table
	data_file = open(filename,'w')
	data_file.write('q CF_edge CF_round CFD_edge CFD_round KF_edge KF_round WR rho0z rsz rtz rho0x rsx rtx\n')
	## 4.1. Prolate models
	for q,m,rq in zip(qrange,pro_models,FF_q):
		X = ret_2.find_J_value(m,angle=angle,los='y',flattening=q)
		Y = ret_2.find_J_value(m,angle=angle,los='x',flattening=1.)
		WR = 1./X[-3] ## Kinematic ratio
		KFr = ((2.*WR+1.)/3.)**2 ## Kinematic correction factor round
		KFe = ((2.+1./WR)/3.)**2
		data_file.write("%0.2f "+"%0.4f "*13+"\n"%(1./q,
		                np.log10(X[-2])-np.log10(rq[0]),
		                np.log10(Y[-2])-sph_J,
		                np.log10(X[-1])-np.log10(rq[1]),
		                np.log10(Y[-1])-sph_D,
		                KFe,KFr,WR,X[0],X[1],X[2],Y[0],Y[1],Y[2]))
	## 4.2. Spherical model
	for l,m in zip([1.,1.],extra_models):
		X = ret_2.find_J_value(m,angle=angle,los='z',flattening=1.)
		Y = ret_2.find_J_value(m,angle=angle,los='x',flattening=1.)
		WR = 1./X[-3]
		KFr = ((2.*WR+1.)/3.)**2
		KFe = ((2.+1./WR)/3.)**2
		data_file.write("%0.2f "+"%0.4f "*13+"\n"%(l,
		                np.log10(X[-2])-sph_J,
		                np.log10(Y[-2])-sph_J,
		                np.log10(X[-1])-sph_D,
		                np.log10(Y[-1])-sph_D,
		                KFe,KFr,WR,X[0],X[1],X[2],Y[0],Y[1],Y[2]))
	## 4.3. Oblate models
	for q,m,rq in zip(qrange,obl_models,FF_q):
		X = ret_2.find_J_value(m,angle=angle,los='x',flattening=q)
		Y = ret_2.find_J_value(m,angle=angle,los='z',flattening=1.)
		WR = X[-3]
		KFr = ((2.*WR+1.)/3.)**2
		KFe = ((2.+1./WR)/3.)**2
		data_file.write("%0.2f "+"%0.4f "*13+"\n"%(q,
		                np.log10(X[-2])-np.log10(rq[0]),
		                np.log10(Y[-2])-sph_J,
		                np.log10(X[-1])-np.log10(rq[1]),
		                np.log10(Y[-1])-sph_D,
		                KFe,KFr,WR,X[0],X[1],X[2],Y[0],Y[1],Y[2]))
	data_file.close()
	## ========================================================================

if __name__ == '__main__':
	compute_table()
