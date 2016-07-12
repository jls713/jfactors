## Fit of spread in F_J due to triaxiality for SEG (2016)
import numpy as np
import matplotlib.pyplot as plt
import emcee
import pandas as pd
import corner
## ============================================================================
## Prior
def logp(x):
    if(x[0]<0.1 or x[0]>10.):
        return -np.inf
    return 0.
## Model
def model(fj,x):
    mm=1.-np.power(10.,-(fj-x[1])/x[0])
    return mm

def model2(fj,x):
    mm = 1.+(fj-x[1])/x[0]
    return mm

## Likelihood
def logl_single(x,fj,e_mean,e_err,fnn):
    '''
        Evaluates single likelihood of error in F_J given ellipticity
        distributed normally with mean e_mean and s.d. e_err
        (can be asymmetric error-bars).
    '''
    oe = fnn(fj,x)
    if(oe<0.):
        return -np.inf
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
def logl(x,FJ,e_mean,e_err,fnn):
    '''
        Evaluates total likelihood of error in F_J given ellipticity
        distributed normally with mean e_mean and s.d. e_err
        (can be asymmetric error-bars).
    '''
    p = logp(x)
    if p!=0.:
        return p
    return p+np.sum([logl_single(x,fj,e_m,e_s,fnn) for fj,e_m,e_s in zip(FJ,e_mean,e_err)])
## ============================================================================
## Main function
def samples(ffile,cut=1.):
    ''' Generates a set of samples  of (a,b) that produce the ellipticities
        given the spread in the correction factors due to triaxiality for a
        model of the form fnn
        '''
    data = pd.read_csv(ffile,sep=' ',
                       names=['Name','e','e_e1','e_e2',
                                'FJU','FJU_e1','FJU_e2',
                                'FJR','FJR_e1','FJR_e2',
                                'FJT','FJT_e1','FJT_e2'],na_values='nan',index_col=False)
    data = data[data.e<cut].reset_index(drop=True)
    ndim,nwalkers=2,50
    lo=[0.4,0.]
    hi=[1.,0.2]
    ### 1. Compute samples from emcee
    p0=np.array([np.random.uniform(low=lo[k],high=hi[k],size=nwalkers) for k in range(ndim)]).T
    meanfJ = .5*(data['FJU_e1'].values+data['FJU_e2'].values)
    e_m = data['e'].values
    e_s = np.vstack((data['e_e1'].values,data['e_e2'].values)).T
    fnn = model
    sampler = emcee.EnsembleSampler(nwalkers,ndim,logl,args=[meanfJ,e_m,e_s,fnn])
    pos,prob,state=sampler.run_mcmc(p0,1500)
    sampler.reset()
    pos,prob,state=sampler.run_mcmc(pos,500)
    samples = sampler.chain.reshape((-1,ndim))
    fig=corner.corner(samples,labels=[r'$a$',r'$b$'])
    fig.savefig("triangle.png")
    print np.median(samples,axis=0)
    print np.std(samples,axis=0)
    a1,b1=np.median(samples,axis=0)

    fnn = model2
    p0=np.array([np.random.uniform(low=lo[k],high=hi[k],size=nwalkers) for k in range(ndim)]).T
    sampler = emcee.EnsembleSampler(nwalkers,ndim,logl,args=[meanfJ,e_m,e_s,fnn])
    pos,prob,state=sampler.run_mcmc(p0,1500)
    sampler.reset()
    pos,prob,state=sampler.run_mcmc(pos,500)
    samples = sampler.chain.reshape((-1,ndim))
    fig=corner.corner(samples,labels=[r'$a$',r'$b$'])
    fig.savefig("triangle.png")
    print np.median(samples,axis=0)
    print np.std(samples,axis=0)
    a2,b2=np.median(samples,axis=0)

    plt.figure(figsize=[3.,2.5])
    plt.errorbar(data['e'].values,meanfJ,xerr=[data['e_e2'].values,data['e_e1'].values],fmt='o')
    ee = np.linspace(0.,1.)
    plt.plot(ee,b1-a1*np.log10(1.-ee))
    plt.plot(ee,-a2*(1.-ee)+b2)
    plt.savefig('fit_result.pdf')

    for ee in [0.2,0.4,0.6]:
        print ee,np.power(10.,b1-a1*np.log10(1.-ee)),np.power(10.,b2-a2*(1.-ee))

if __name__ == '__main__':
    samples('corr_triax_table_ascii.dat',cut=1.)
