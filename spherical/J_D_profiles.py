# -*- coding: utf-8 -*-
### Generates J and D profiles for Fig. 2,3,4,5,6,7 of Evans, Sanders & Geringer-Sameth (2016)
from spherical_Jfactors import *
from J_D_table import *


num_samples=50000 ## Number of samples used for MCing over errors

### Load in required data
data = pd.read_csv('../data/data.dat',sep=' ')
scale_radii=np.genfromtxt('../data/geringer_sameth_gamma.dat',skip_header=49).T[8]
scale_radii=np.power(10.,scale_radii)/1000.
ackermann_data = read_ackermann_data()

def plot_grid(ranges,f,a,plot_out,Name,geo_factor=True):
    ''' Plots a series of J profiles for the dwarfs in ranges from the data
        table data.dat. If geo_factor==True the half-light radii are
        ellipticity-corrected by a factor sqrt(1-e). '''
    plt.subplots_adjust(wspace=0.,hspace=0.)

    a[0][len(a[0])-1].annotate(Name, xy=(0.95,1.03),
                     xycoords='axes fraction', fontsize=16,
                     horizontalalignment='right', verticalalignment='bottom')
    k,j=0,0

    rnfwrs = 5.

    for i in ranges:

        print data['Name'][i]

        angles = np.linspace(1e-5,0.55,100)
        a[k][j].set_xlim(0.,0.55)
        a[k][j].set_ylim(14.,23.4)

        geof = 1.
        if(geo_factor):
            geof=np.sqrt(1.-data['ellip'][i])

        ## NFW
        J = np.array([
             map(lambda a:sample_errorsJ(
                            data['sigma_los'][i],
                            [data['esigma_los2'][i],data['esigma_los1'][i]],
                            data['R_half'][i]*geof,
                           [data['eR_half2'][i]*geof,data['eR_half1'][i]*geof],
                            data['D'][i],data['eD'][i],
                            a,[1e-10,1e-10], ## angle with no associated errors
                            gamma=1.,N=num_samples,
                            nfw=rnfwrs*data['R_half'][i]*geof/1000.,
                          walker_or_wolf="walker"),
                          angles)])[0]

        a[k][j].fill_between(angles,J.T[0]-J.T[1],J.T[0]+J.T[2],
                            alpha=0.5,color=sns.color_palette()[0],label='NFW')
        a[k][j].plot(angles,J.T[0],color=sns.color_palette()[0])

        if(data['theta_max'][i]!=0.5):
            a[k][j].errorbar([0.5],[data['Jhalf'][i]],
                             yerr=[[data['eJhalf2'][i]],[data['eJhalf1'][i]]],
                             fmt='.',color='k',capsize=0.)
        if(data['theta_half'][i]!=0.5):
            a[k][j].errorbar([data['theta_half'][i]],
                             [data['Jmax'][i]-np.log10(2.)],
                             yerr=[[data['eJmax2'][i]],[data['eJmax1'][i]]],
                             fmt='.',color='k',capsize=0.)

        ## Gamma=1 cusp
        J = np.array([
             map(lambda a:sample_errorsJ(
                            data['sigma_los'][i],
                            [data['esigma_los2'][i],data['esigma_los1'][i]],
                            data['R_half'][i]*geof,
                           [data['eR_half2'][i]*geof,data['eR_half1'][i]*geof],
                            data['D'][i],data['eD'][i],
                            a,[1e-10,1e-10], ## angle with no associated errors
                            gamma=1.,N=num_samples,
                          walker_or_wolf="walker"), angles)])[0]

        a[k][j].plot(angles,J.T[0],color=sns.color_palette()[2])
        a[k][j].fill_between(angles,J.T[0]-J.T[1],J.T[0]+J.T[2],
                    alpha=0.5,color=sns.color_palette()[2],label=r'$\gamma=1$')

        ## Gamma=0.51 cusp
        J = np.array([
             map(lambda a:sample_errorsJ(
                            data['sigma_los'][i],
                            [data['esigma_los2'][i],data['esigma_los1'][i]],
                            data['R_half'][i]*geof,
                           [data['eR_half2'][i]*geof,data['eR_half1'][i]*geof],
                            data['D'][i],data['eD'][i],
                            a,[1e-10,1e-10], ## angle with no associated errors
                            gamma=0.51,N=num_samples,
                          walker_or_wolf="walker"), angles)])[0]

        a[k][j].plot(angles,J.T[0],color=sns.color_palette()[4])
        a[k][j].fill_between(angles,J.T[0]-J.T[1],J.T[0]+J.T[2],
                 alpha=0.5,color=sns.color_palette()[4],label=r'$\gamma=0.51$')


        if(k==len(a)-1):
            a[k][j].set_xlabel(r'$\theta/^\circ$')
        else:
            a[k][j].set_xticklabels([])
        if(j==0):
            a[k][j].set_ylabel(r'$\log_{10}(\mathrm{J}(\theta)/\,\mathrm{GeV^2\,cm}^{-5})$')
        else:
            a[k][j].set_yticklabels([])
        a[k][j].annotate(posh_names[data['Name'][i]], xy=(0.95,0.05),
                     xycoords='axes fraction', fontsize=16,
                     horizontalalignment='right', verticalalignment='bottom')

        if(j==1 and k==0):
            a[k][j].legend(loc="lower left",ncol=4, bbox_to_anchor=(-0.85, 1.0))

        legendEntries=[]# list of plots that are going to be in the legend
        legendText=[]   # list of text messages for the legend

        if(data['Jmax'][i]==data['Jmax'][i]):
            GS = a[k][j].errorbar([data['theta_max'][i]],[data['Jmax'][i]],yerr=[[data['eJmax2'][i]],[data['eJmax1'][i]]],fmt='.',color='k',capsize=0.)
            legendEntries.append(GS)

            if(data['Name'][i]=="ReticulumII"):
                legendText.append('Bonnivard et al. (2015b)')
            elif(data['Name'][i]=="TucanaII"):
                legendText.append('Walker et al. (2016)')
            else:
                legendText.append('Geringer-Sameth et al. (2015)')

        bonni_data = read_bonnivard_table(data['Name'][i])

        if not bonni_data.empty:
            B, = a[k][j].plot(bonni_data['alpha'],bonni_data['J'],color='k')
            if data['Name'][i] not in ["ReticulumII","TucanaII"]:
                legendEntries.append(B)
                legendText.append('Bonnivard et al. (2015)')
            l,=a[k][j].plot(bonni_data['alpha'],bonni_data['eJm68'],ls='dashed',color='k')
            l.set_dashes((2,1))
            l,=a[k][j].plot(bonni_data['alpha'],bonni_data['eJp68'],ls='dashed',color='k')
            l.set_dashes((2,1))

        ackermann = ackermann_data[ackermann_data.name==data['Name'][i]]
        if(len(ackermann)>0):
            A = a[k][j].errorbar(0.5,ackermann['J'],yerr=ackermann['eJ'],ms=3,color='r',fmt='d',capsize=0.)
            if data['Name'][i] not in ["ReticulumII","TucanaII"]:
                legendEntries.append(A)
                legendText.append('Ackermann et al. (2014)')

        if((j==0 and k==0) or (data['Name'][i]=="ReticulumII" or data['Name'][i]=="TucanaII")):
            lgd = a[k][j].legend(legendEntries,legendText,numpoints=1,loc='upper left',fontsize=8)

        j+=1
        if(j==len(a[0])):
            j=0
            k+=1

    while(j<len(a[0]) and j!=0):
        a[k][j].axis('off')
        j+=1
    plt.savefig(plot_out,bbox_inches='tight')
    plt.clf()


def plot_decay_grid(ranges,f,a,plot_out,Name,geo_factor=True):
    ''' Plots a series of D profiles for the dwarfs in ranges from the data
        table data.dat. If geo_factor==True the half-light radii are
        ellipticity-corrected by a factor sqrt(1-e). '''
    plt.subplots_adjust(wspace=0.,hspace=0.)

    a[0][len(a[0])-1].annotate(Name, xy=(0.95,1.03),
                     xycoords='axes fraction', fontsize=16,
                     horizontalalignment='right', verticalalignment='bottom')
    k,j=0,0

    rnfwrs = 5.

    for i in ranges:

        print data['Name'][i]

        angles = np.linspace(1e-5,0.55,100)
        a[k][j].set_xlim(0.,0.55)
        a[k][j].set_ylim(14.,20.4)

        geof = 1.
        if(geo_factor):
            geof=np.sqrt(1.-data['ellip'][i])

        ## NFW profile
        J = np.array([
             map(lambda a:sample_errorsD(
                           data['sigma_los'][i],
                           [data['esigma_los2'][i],data['esigma_los1'][i]],
                           data['R_half'][i]*geof,
                           [data['eR_half2'][i]*geof,data['eR_half1'][i]*geof],
                           data['D'][i],data['eD'][i],
                           a,[1e-10,1e-10], ## angle with no associated errors
                           gamma=1.,N=num_samples,
                           nfw=rnfwrs*data['R_half'][i]*geof/1000.,
                          walker_or_wolf="walker"),
                        angles)])[0]


        a[k][j].fill_between(angles,J.T[0]-J.T[1],J.T[0]+J.T[2],
                         alpha=0.5,color=sns.color_palette()[0],label='NFW')
        a[k][j].plot(angles,J.T[0],color=sns.color_palette()[0])
        if(data['theta_max'][i]!=0.5):
            a[k][j].errorbar([0.5],[data['dJhalf'][i]],
                             yerr=[[data['edJhalf2'][i]],[data['edJhalf1'][i]]],
                             fmt='.',color='k',capsize=0.)
        if(data['theta_half'][i]!=0.5):
            a[k][j].errorbar([data['dtheta_half'][i]],
                             [data['dJmax'][i]-np.log10(2.)],
                             yerr=[[data['edJmax2'][i]],[data['edJmax1'][i]]],
                             fmt='.',color='k',capsize=0.)

        ## Gamma=1.01 cusp
        J = np.array([
             map(lambda a:sample_errorsD(
                           data['sigma_los'][i],
                           [data['esigma_los2'][i],data['esigma_los1'][i]],
                           data['R_half'][i]*geof,
                           [data['eR_half2'][i]*geof,data['eR_half1'][i]*geof],
                           data['D'][i],data['eD'][i],
                           a,[1e-10,1e-10], ## angle with no associated errors
                           gamma=1.01,N=num_samples,
                          walker_or_wolf="walker"), angles)])[0]
        a[k][j].plot(angles,J.T[0],color=sns.color_palette()[2])
        a[k][j].fill_between(angles,J.T[0]-J.T[1],J.T[0]+J.T[2],
                 alpha=0.5,color=sns.color_palette()[2],label=r'$\gamma=1.01$')

        J = np.array([
             map(lambda a:sample_errorsD(
                           data['sigma_los'][i],
                           [data['esigma_los2'][i],data['esigma_los1'][i]],
                           data['R_half'][i]*geof,
                           [data['eR_half2'][i]*geof,data['eR_half1'][i]*geof],
                           data['D'][i],data['eD'][i],
                           a,[1e-10,1e-10], ## angle with no associated errors
                           gamma=1.49,N=num_samples,
                          walker_or_wolf="walker"), angles)])[0]
        a[k][j].plot(angles,J.T[0],color=sns.color_palette()[4])
        a[k][j].fill_between(angles,J.T[0]-J.T[1],J.T[0]+J.T[2],
                 alpha=0.5,color=sns.color_palette()[4],label=r'$\gamma=1.49$')

        if(k==len(a)-1):
            a[k][j].set_xlabel(r'$\theta/^\circ$')
        else:
            a[k][j].set_xticklabels([])
        if(j==0):
            a[k][j].set_ylabel(r'$\log_{10}(\mathrm{D}(\theta)/\,\mathrm{GeV\,cm}^{-2})$')
        else:
            a[k][j].set_yticklabels([])

        a[k][j].annotate(posh_names[data['Name'][i]], xy=(0.95,0.05),
                     xycoords='axes fraction', fontsize=16,
                     horizontalalignment='right', verticalalignment='bottom')

        if(j==1 and k==0):
            a[k][j].legend(loc="lower left",ncol=4, bbox_to_anchor=(-0.85, 1.0))

        legendEntries=[]# list of plots that are going to be in the legend
        legendText=[]   # list of text messages for the legend

        if(data['Jmax'][i]==data['Jmax'][i]):
            GS = a[k][j].errorbar([data['theta_max'][i]],[data['dJmax'][i]],yerr=[[data['edJmax2'][i]],[data['edJmax1'][i]]],fmt='.',color='k',capsize=0.)
            legendEntries.append(GS)
            if(data['Name'][i]=="ReticulumII"):
                legendText.append('Bonnivard et al. (2015b)')
            elif(data['Name'][i]=="TucanaII"):
                legendText.append('Walker et al. (2016)')
            else:
                legendText.append('Geringer-Sameth et al. (2015)')

        bonni_data = read_bonnivard_table_decay(data['Name'][i])

        if not bonni_data.empty:
            B,=a[k][j].plot(bonni_data['alpha'],bonni_data['D'],color='k')
            if data['Name'][i] not in ["ReticulumII","TucanaII"]:
                legendEntries.append(B)
                legendText.append('Bonnivard et al. (2015)')
            l,=a[k][j].plot(bonni_data['alpha'],bonni_data['eDm68'],ls='dashed',color='k')
            l.set_dashes((2,1))
            l,=a[k][j].plot(bonni_data['alpha'],bonni_data['eDp68'],ls='dashed',color='k')
            l.set_dashes((2,1))

        if((j==0 and k==0) or (data['Name'][i] in ["ReticulumII","TucanaII"])):
            lgd = a[k][j].legend(legendEntries,legendText,numpoints=1,loc='upper left',fontsize=8)

        j+=1
        if(j==len(a[0])):
            j=0
            k+=1
    while(j<len(a[0]) and j!=0):
        a[k][j].axis('off')
        j+=1
    plt.savefig(plot_out,bbox_inches='tight')
    plt.clf()

def generate_paper_plots():
    ''' Generates the grids of plots for Paper I '''
    using_geo_factor=True

    f,a=plt.subplots(2,4,figsize=(9,4.5))
    plot_grid(range(8),f,a,'cd_grid.pdf','Classical Dwarfs',using_geo_factor)

    f,a=plt.subplots(4,4,figsize=(9,9))
    plot_grid(range(8,23),f,a,'uf_grid.pdf','Ultrafaints',using_geo_factor)

    f,a=plt.subplots(1,4,figsize=(9,2.5))
    plot_grid(range(23,27),f,[a],'predict_grid.pdf','Predictions',using_geo_factor)


    f,a=plt.subplots(2,4,figsize=(9,4.5))
    plot_decay_grid(range(8),f,a,'cd_grid_decay.pdf','Classical Dwarfs',using_geo_factor)

    f,a=plt.subplots(4,4,figsize=(9,9))
    plot_decay_grid(range(8,23),f,a,'uf_grid_decay.pdf','Ultrafaints',using_geo_factor)

    f,a=plt.subplots(1,4,figsize=(9,2.5))
    plot_decay_grid(range(23,27),f,[a],'predict_grid_decay.pdf','Predictions',using_geo_factor)

def generate_ret2_comparison():
    ''' Generates a plot to compare the estimates for the J-factor for Ret II
        from Paper I to Geringer-Sameth et al. (2015)
        NOT ELLIPTICITY-CORRECTED
        '''
    plt.clf()
    f=plt.figure(figsize=[3.32,2.5])
    i=20
    angles=np.array([0.5])
    xx = np.linspace(0.505,1.495)
    Jres = np.zeros(len(xx))
    Jup = np.zeros(len(xx))
    Jdo = np.zeros(len(xx))
    for j in range(len(Jres)):
        J = np.array([map(lambda a:sample_errorsJ(data['sigma_los'][i],[data['esigma_los2'][i],data['esigma_los1'][i]],data['R_half'][i],[data['eR_half2'][i],data['eR_half1'][i]],data['D'][i],data['eD'][i],a,[1e-10,1e-10],gamma=xx[j],N=num_samples), angles)])[0]
        Jres[j]=J[0][0]
        Jdo[j]=J[0][0]-J[0][1]
        Jup[j]=J[0][0]+J[0][2]
    plt.plot(xx,Jres,color=sns.color_palette()[0])
    plt.fill_between(xx,Jdo,Jup,alpha=0.5,color=sns.color_palette()[0])
    plt.errorbar([1.],[data['Jhalf'][i]],yerr=[[data['eJhalf2'][i]],[data['eJhalf1'][i]]],fmt='.',color='k')
    J = np.array([map(lambda a:sample_errorsJ(data['sigma_los'][i],[data['esigma_los2'][i],data['esigma_los1'][i]],data['R_half'][i],[data['eR_half2'][i],data['eR_half1'][i]],data['D'][i],data['eD'][i],a,[1e-10,1e-10],gamma=1.,N=num_samples,nfw=10.*data['R_half'][i]/1000.), angles)])[0]
    plt.errorbar([1.03],[J[0][0]],yerr=[[J[0][1]],[J[0][2]]],fmt='.',color='r')

    J = np.array([map(lambda a:sample_errorsJ(data['sigma_los'][i],[data['esigma_los2'][i],data['esigma_los1'][i]],data['R_half'][i],[data['eR_half2'][i],data['eR_half1'][i]],data['D'][i],data['eD'][i],a,[1e-10,1e-10],gamma=1.,N=num_samples,nfw=2.*data['R_half'][i]/1000.), angles)])[0]
    plt.errorbar([1.01],[J[0][0]],yerr=[[J[0][1]],[J[0][2]]],fmt='.',color='r')

    plt.annotate(data['Name'][i], xy=(0.95,0.95),
                     xycoords='axes fraction', fontsize=16,
                     horizontalalignment='right', verticalalignment='top')
    plt.annotate('Geringer-Sameth et al. (2015)', xy=(0.05,0.1),
                     xycoords='axes fraction', fontsize=9,
                     horizontalalignment='left', verticalalignment='top')
    plt.annotate('NFW, this work', xy=(0.5,0.45),
                     xycoords='axes fraction', fontsize=9,
                     horizontalalignment='left', verticalalignment='bottom',color='r')
    plt.xlabel(r'$\gamma$')
    plt.ylabel(r'$\log_{10}(J(\alpha)/\,\mathrm{GeV^2\,cm}^{-5})$')
    plt.savefig('ret_2.pdf',bbox_inches='tight')

if __name__ == '__main__':
    generate_paper_plots()
