## Fit of spread in F_J due to triaxiality for SEG (2016)
import numpy as np
import matplotlib.pyplot as plt
import emcee
import pandas as pd
import corner
## ============================================================================
## Prior
def logp(x):
    if(x[0]<0.):
        return -np.inf
    return 0.
## Model
def model(fj,x):
    return 1.-np.power(10.,-fj/x[0])
    #return x[0]*fj**x[1]
## Likelihood
def logl_single(x,fj,e_mean,e_err):
    '''
        Evaluates single likelihood of error in F_J given ellipticity
        distributed normally with mean e_mean and s.d. e_err
        (can be asymmetric error-bars).
    '''
    oe = model(fj,x)
    if(e_err[0]!=e_err[0] or e_err[0]==np.nan):
        e_err=np.array([e_mean/2.,0.])
        e_mean=0.
    if isinstance(e_err,list) or isinstance(e_err,np.ndarray):
        if(oe>e_mean):
            return -(e_mean-oe)**2/2./e_err[0]**2 - .5*np.log(2.*np.pi*e_err[0]**2)
        if(oe<e_mean):
            return -(e_mean-oe)**2/2./e_err[1]**2 - .5*np.log(2.*np.pi*e_err[1]**2)
    else:
        return -(e_mean-oe)**2/2./e_err**2 - .5*np.log(2.*np.pi*e_err**2)
## Full likelihood
def logl(x,FJ,e_mean,e_err):
    '''
        Evaluates total likelihood of error in F_J given ellipticity
        distributed normally with mean e_mean and s.d. e_err
        (can be asymmetric error-bars).
    '''
    p = logp(x)
    if p!=0.:
        return p
    return p+np.sum([logl_single(x,fj,e_m,e_s) for fj,e_m,e_s in zip(FJ,e_mean,e_err)])
## ============================================================================
## Main function
def samples(ffile,cut=1.):
    ''' Generates a set of samples  of (a,b) that produce the ellipticities
        given the spread in the correction factors due to triaxiality for a
        model of the form e = a*fJ**b
        '''
    data = pd.read_csv(ffile,sep=' ',
                       names=['Name','e','e_e1','e_e2',
                                'FJU','FJU_e1','FJU_e2',
                                'FJR','FJR_e1','FJR_e2',
                                'FJT','FJT_e1','FJT_e2'],na_values='nan',index_col=False)
    data = data[data.e<cut].reset_index(drop=True)
    ndim,nwalkers=2,50
    lo=[0.,0.]
    hi=[1.,1.]
    ### 1. Compute samples from emcee
    p0=np.array([np.random.uniform(low=lo[k],high=hi[k],size=nwalkers) for k in range(ndim)]).T
    meanfJ = .5*(data['FJU_e1'].values+data['FJU_e1'].values)
    e_m = data['e'].values
    e_s = np.vstack((data['e_e1'].values,data['e_e2'].values)).T
    sampler = emcee.EnsembleSampler(nwalkers,ndim,logl,args=[meanfJ,e_m,e_s])
    pos,prob,state=sampler.run_mcmc(p0,1500)
    sampler.reset()
    pos,prob,state=sampler.run_mcmc(pos,500)
    samples = sampler.chain.reshape((-1,ndim))
    fig=corner.corner(samples,labels=[r'$a$',r'$b$'])
    fig.savefig("triangle.png")
    print np.median(samples,axis=0)
    print np.std(samples,axis=0)

if __name__ == '__main__':
    samples('corr_triax_table_ascii.dat',cut=0.6)
