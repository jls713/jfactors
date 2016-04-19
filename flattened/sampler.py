## Generate samples from triaxiality distributions for Figures 9 & 10 and Table 5  of Sanders, Evans & Geringer-Sameth
## ============================================================================
import numpy as np
from numpy import sqrt,cos,sin
import emcee
import corner
import sys
sys.path.append('../../../code/jfactors/')
import jfactors_py as jf
import pandas as pd
## ============================================================================
def observed_ellipticity(ba,ca,theta,phi):
    ''' Observed ellipticity given intrinsic axis ratios and spherical polar line-of-sight (theta,phi)
        --- Contopoulos 1956, Weijmans 2014 Appendix '''
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
    ''' Prior distribution - limiting 0<T<1, 0<E<1, 0<theta<pi/2, 0<phi<pi/2'''
    if(x[0]<0. or x[0]>1.):
        return -np.inf
    if(x[1]<0. or x[1]>1.):
        return -np.inf
    if(x[2]<0. or x[2]>.5*np.pi):
        return -np.inf
    if(x[3]<0. or x[3]>.5*np.pi):
        return -np.inf
    else:
        return 0.

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
    p = logp(x)
    if p!=0.:
    	return p
    ba,ca=axis_ratios(x[0],x[1])
    oe = observed_ellipticity(ba,ca,x[2],x[3])
    if(ma_prior):
        p+=major_axis_prior(x[2],x[3])
    if(sj_prior):
        p+=sanchez_janssen_prior(x[0],x[1])
    if isinstance(e_err,list):
        if(oe>e_mean):
            return p-(e_mean-oe)**2/2./e_err[0]**2 - .5*np.log(2.*np.pi*e_err[0]**2)
        if(oe<e_mean):
            return p-(e_mean-oe)**2/2./e_err[1]**2 - .5*np.log(2.*np.pi*e_err[1]**2)
    return p-(e_mean-oe)**2/2./e_err**2 - .5*np.log(2.*np.pi*e_err**2)

## ============================================================================
## Main function
def samples(e_m,e_s,size,ffile,ma_prior=False,sj_prior=False,geo_factor=True):
    ''' Generates a set of <size> samples  of (T,E,theta,phi) that produce the
        distribution of observed ellipticity for each dwarf. For each sample
        calculates the correction factor and outputs the results to ffile
        ma_prior is a flag for the Major-Axis prior
        sj_prior is a flag for the Sanchez-Janssen prior
        geo_factor uses an additional factor sqrt(1-e) in spherical model
        '''
    ndim,nwalkers=4,50
    lo=[0.,0.,0.,0.]
    hi=[1.,1.,.5*np.pi,.5*np.pi]
    if(isinstance(e_s,list)):
        if(e_s[0]!=e_s[0]):
            e_s=[e_m/2.,0.]
            e_m=0.
    else:
        if(e_s!=e_s):
            e_s=ett/2.
            e_m=0.
    ### 1. Compute samples from emcee
    p0=np.array([np.random.uniform(low=lo[k],high=hi[k],size=nwalkers) for k in range(ndim)]).T
    sampler = emcee.EnsembleSampler(nwalkers,ndim,logl,args=[e_m,e_s,ma_prior,sj_prior])
    pos,prob,state=sampler.run_mcmc(p0,15000)
    sampler.reset()
    pos,prob,state=sampler.run_mcmc(pos,5000)
    samples = sampler.chain.reshape((-1,ndim))
    fig=corner.corner(samples,labels=[r'$T$',r'$E$',r'$\theta$',r'$\phi$'])
    fig.savefig("triangle.png")
    print np.median(pos,axis=0),np.std(pos,axis=0)

    ### 2. Take <size> samples and compute correction factors
    samples_redux=samples[np.random.randint(len(samples),size=size)]
    ssp=jf.PaperModel(0.999,0.99,0.05,3.22,True)
    Jt = ssp.J_factor(0.01,0.01,30.,0.5,False,False)[0]

    ff = open(ffile,'w')
    for k,i in enumerate(samples_redux):
        print i
        ss=jf.PaperModel(ba_from_Tca(i[0],i[1]),i[1],0.05,3.22,True)
        if(geo_factor):
            ssp=jf.PaperModel(0.999,0.99,0.05*np.sqrt(ss.ellipticity(i[2],i[3])),3.22,True)
            Jt = ssp.J_factor(0.01,0.01,30.,0.5,False,False)[0]
        rr=np.log10(ss.J_factor(i[2],i[3],30.,0.5,False,False)[0]/Jt)
        ba,ca=axis_ratios(i[0],i[1])
        oe = observed_ellipticity(ba,ca,i[2],i[3])
        ll = logl(i,e_m,e_s,ma_prior,sj_prior)
        ff.write('%0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n'%(i[0],i[1],i[2],i[3],rr,oe,ll))
    ff.close()

## ============================================================================
def run_grid(geo_factor=True):
    ''' For each dwarf compute sample of correction factors under the
        three assumptions and output to file '''
    data = pd.read_csv('../data.dat',sep=' ')
    N=100
    for i in range(len(data)):
        samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],N,
                'triaxial_results/'+data.Name[i]+'_nop',False,False,
                geo_factor=geo_factor)
        samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],N,
                'triaxial_results/'+data.Name[i]+'_ma',True,False,
                geo_factor=geo_factor)
        samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],N,
                'triaxial_results/'+data.Name[i]+'_sj',False,True,
                geo_factor=geo_factor)

def ret2(geo_factor=True):
    ''' For RetII compute sample of correction factors under the
        three assumptions and output to file (higher res than above) '''
    data = pd.read_csv('../data.dat',sep=' ')
    i=21
    samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],400,
            'triaxial_results/'+data.Name[i]+'_nop_hr',False,False,
            geo_factor=geo_factor)
    samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],400,
            'triaxial_results/'+data.Name[i]+'_ma_hr',True,False,
            geo_factor=geo_factor)
    samples(data.ellip[i],[data.ellip_e1[i],data.ellip_e2[i]],400,
            'triaxial_results/'+data.Name[i]+'_sj_hr',False,True,
            geo_factor=geo_factor)
## ============================================================================
if __name__=="__main__":
    run_grid()
    ret2()
## ============================================================================
