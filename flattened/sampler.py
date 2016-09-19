## Generate samples from triaxiality distributions for Figures 9 & 10 and Table 5  of Sanders, Evans & Geringer-Sameth
## ============================================================================
import numpy as np
from numpy import sqrt,cos,sin
import emcee
# import corner
import sys
sys.path.append('/home/jls/work/code/jfactors/')
import jfactors_py as jf
import pandas as pd
## ============================================================================
def observed_ellipticity(ba,ca,theta,phi):
    '''
        Observed ellipticity given intrinsic axis ratios and spherical polar line-of-sight (theta,phi)
        --- Contopoulos 1956, Weijmans 2014 Appendix
        This gives 1-\epsilon
    '''
    if(ca==1. and ba==1.):
    	return 1.
    st,ct=sin(theta),cos(theta)
    sp,cp=sin(phi),cos(phi)
    ba2=ba*ba
    ca2=ca*ca
    m1ba2 = 1.-ba2
    m1ca2 = 1.-ca2
    A = m1ca2*ct*ct+m1ba2*st*st*sp*sp+ba2+ca2
    B = m1ca2*ct*ct-m1ba2*st*st*sp*sp-ba2+ca2
    B*=B
    B+=4.*m1ca2*m1ba2*st*st*ct*ct*sp*sp
    B = sqrt(B)
    return sqrt((A-B)/(A+B))

def ba_from_Tca(T,ca):
    ''' b/a given triaxiality T and c/a '''
    ba = sqrt(1.-T*(1.-ca*ca))
    return ba

def axis_ratios(T,E):
    ''' (b/a,c/a) given triaxiality T and intrinsic ellipticity E=1-c/a '''
    return ba_from_Tca(T,1.-E),1.-E
## ============================================================================
## Priors
def logp(x):
    '''
        Prior distribution - limiting 0<T<1, 0<E<1, 0<theta<pi/2, 0<phi<pi/2
        and a factor sin(theta) for uniform sampling over sphere
    '''
    if(x[0]<0. or x[0]>1.):
        return -np.inf
    if(x[1]<0. or x[1]>1.):
        return -np.inf
    if(x[2]<0. or x[2]>.5*np.pi):
        return -np.inf
    if(x[3]<0. or x[3]>.5*np.pi):
        return -np.inf
    else:
        return np.log(np.sin(x[2]))

def major_axis_prior(theta,phi):
    ''' Major axis within 0.1 rad of line-of-sight '''
    sigma_phi=0.1
    sigma_theta=0.1
    return -phi**2/2./sigma_phi**2-(theta-.5*np.pi)**2/2./sigma_theta**2

def sanchez_janssen_prior(T,E):
    ''' Triaxiality = N(0.55, 0.04) and Ellipticity  (0.51,0.12) from
        Sanchez-Janssen et al. (2016) '''
    T0=0.55
    sigma_T=0.04
    E0=0.51
    sigma_E=0.12
    return -(T-T0)**2/2./sigma_T**2-(E-E0)**2/2./sigma_E**2

## ============================================================================
## Likelihood
def logl(x,e_mean,e_err,ma_prior,sj_prior):
    '''
        Evaluates likelihood of observed ellipticity x given ellipticity
        distributed normally with mean e_mean and s.d. e_err
        (can be asymmetric error-bars).
        ma_prior is a flag for the Major-Axis prior
        sj_prior is a flag for the Sanchez-Janssen prior
    '''
    ## Evaluate priors
    p = logp(x)
    if p==-np.inf:
    	return p
    ba,ca=axis_ratios(x[0],x[1])
    if(ca<0.05): ## a prior
        return -np.inf
    if(ma_prior):
        p+=major_axis_prior(x[2],x[3])
    if(sj_prior):
        p+=sanchez_janssen_prior(x[0],x[1])
    ## Evaluate full posterior
    oe = observed_ellipticity(ba,ca,x[2],x[3])
    if isinstance(e_err,list):
        if(oe>e_mean):
            return p-(e_mean-oe)**2/2./e_err[0]**2 - .5*np.log(2.*np.pi*e_err[0]**2)
        if(oe<e_mean):
            return p-(e_mean-oe)**2/2./e_err[1]**2 - .5*np.log(2.*np.pi*e_err[1]**2)
    return p-(e_mean-oe)**2/2./e_err**2 - .5*np.log(2.*np.pi*e_err**2)

## ============================================================================
## Main function
def compute_samples(nwalkers,e_m,e_s,ma_prior=False,sj_prior=False,withplots=None,nsteps=5000):
    ''' Generates a set of samples  of (T,E,theta,phi) that produce the
        distribution of observed ellipticity for each dwarf.
        ma_prior is a flag for the Major-Axis prior
        sj_prior is a flag for the Sanchez-Janssen prior
        e_m = (1-epsilon), e_s is associated error -- can be list
        (upper, lower errorbar)
    '''
    ndim=4
    ## Construct random initial samples uniformly distributed between lo and hi
    lo=[0.,0.,0.,0.]
    hi=[1.,1.,.5*np.pi,.5*np.pi]
    p0=np.array([np.random.uniform(low=lo[k],high=hi[k],size=nwalkers) for k in range(ndim)]).T

    sampler = emcee.EnsembleSampler(nwalkers,ndim,logl,
                                    args=[e_m,e_s,ma_prior,sj_prior])
    ## We do a burn-in of 3 n_steps
    pos,prob,state=sampler.run_mcmc(p0,3*nsteps)
    sampler.reset()
    ## Now production run
    pos,prob,state=sampler.run_mcmc(pos,nsteps)
    samples = sampler.chain.reshape((-1,ndim))
    if(withplots):
        fig=corner.corner(samples,labels=[r'$T$',r'$E$',r'$\theta$',r'$\phi$'])
        fig.savefig(withplots)
    print np.median(pos,axis=0),np.std(pos,axis=0)
    print np.median(sampler.acceptance_fraction)
    return samples

def samples(e_m,e_s,size,ffile,ma_prior=False,sj_prior=False,geo_factor=True,withplots=False):
    ''' Generates a set of <size> samples  of (T,E,theta,phi) that produce the
        distribution of observed ellipticity for each dwarf. For each sample
        calculates the correction factor and outputs the results to ffile
        ma_prior is a flag for the Major-Axis prior
        sj_prior is a flag for the Sanchez-Janssen prior
        geo_factor uses an additional factor sqrt(1-e) in spherical model
        '''
    ### 1. Compute samples from emcee
    nwalkers=50
    ellip = 1.-e_m ## We work with (1-\epsilon) as observed_ellipticity computes 1-\epsilon
    ## The errors can either be a list (upper, lower errorbar) or a single value
    ellip_errors = e_s
    if(isinstance(e_s,list)):
        ellip_errors=[e_s[1],e_s[0]]

    ## If e_s[0] or e_s = nan then e_m is a 90% upper limit and we use
    ## a Gaussian centered on ellip=1 with width e_m/2
    if(isinstance(e_s,list)):
        if(e_s[0]!=e_s[0]):
            ellip_errors=[0.,e_m/2.]
            ellip=1.
    else:
        if(e_s!=e_s):
            ellip_errors=e_m/2.
            ellip=1.

    samples = compute_samples(nwalkers,ellip,ellip_errors,ma_prior,sj_prior,withplots)
    ### 2. Take <size> samples and compute correction factors
    samples_redux=samples[np.random.randint(len(samples),size=size)]
    rh=0.05
    sig=3.22
    Dist=30.
    ang=0.5
    with_multipole=True
    sph_shape = np.array([0.999,0.99])
    ssp=jf.PaperModel(sph_shape[0],sph_shape[1],rh,sig,with_multipole)
    print_messages=False
    withD = False
    sph_viewing = np.array([0.01,0.01])
    Jt = ssp.J_factor(sph_viewing[0],sph_viewing[1],Dist,ang,print_messages,withD,-1.)[0]

    ff = open(ffile,'w')
    for k,i in enumerate(samples_redux):
        ba,ca=axis_ratios(i[0],i[1])
        oe = observed_ellipticity(ba,ca,i[2],i[3])
        print i,ba,ca,oe
        ss=jf.PaperModel(ba,ca,rh,sig,with_multipole)
        if(geo_factor):
            ssp2=jf.PaperModel(sph_shape[0],sph_shape[1],rh*np.sqrt(ss.ellipticity(i[2],i[3])),sig,with_multipole)
            Jt = ssp2.J_factor(sph_viewing[0],sph_viewing[1],Dist,ang,print_messages,withD,-1.)[0]
        rr=np.log10(ss.J_factor(i[2],i[3],Dist,ang,print_messages,withD,-1.)[0]/Jt)
        ll = logl(i,ellip,ellip_errors,ma_prior,sj_prior)
        ff.write('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n'%(i[0],i[1],i[2],i[3],rr,oe,ba,ca,ll))
    ff.close()

## ============================================================================
def run_grid(geo_factor=True):
    ''' For each dwarf compute sample of correction factors under the
        three assumptions and output to file '''
    data = pd.read_csv('../data/data.dat',sep=' ')
    N=500
    for i in range(19,len(data)):
        samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],N,
                'triaxial_results/'+data.Name[i]+'_nop',False,False,
                geo_factor=geo_factor,withplots=None)#'tmp.png')
        samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],N,
                'triaxial_results/'+data.Name[i]+'_ma',True,False,
                geo_factor=geo_factor,withplots=None)#'tmp.png')
        samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],N,
                'triaxial_results/'+data.Name[i]+'_sj',False,True,
                geo_factor=geo_factor,withplots=None)#'tmp.png')

def ret2(geo_factor=True):
    ''' For RetII compute sample of correction factors under the
        three assumptions and output to file (higher res than above) '''
    data = pd.read_csv('../data/data.dat',sep=' ')
    i=21
    samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],400,
            'triaxial_results/'+data.Name[i]+'_nop_hr',False,False,
            geo_factor=geo_factor,withplots=None)#'ret2_dist.png')
    samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],400,
            'triaxial_results/'+data.Name[i]+'_ma_hr',True,False,
            geo_factor=geo_factor,withplots=None)#'ret2_dist.png')
    samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],400,
            'triaxial_results/'+data.Name[i]+'_sj_hr',False,True,
            geo_factor=geo_factor,withplots=None)#'ret2_dist.png')
## ============================================================================
if __name__=="__main__":
    run_grid()
    # ret2()
## ============================================================================
