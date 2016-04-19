# -*- coding: utf-8 -*-
### Generates J and D factor table for Evans, Sanders & Geringer-Sameth (2016)
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.special import gamma as Gamma
from spherical_Jfactors import *

### A set of strings to convert pandas table names into nicer display names

posh_names= {'BootesI':u'Boötes I',
'Carina':'Carina',
'Coma':'Coma Berenices',
'CVnI':'Canes Venatici I',
'CVnII':'Canes Venatici II',
'Draco':'Draco',
'Fornax':'Fornax',
'Hercules':'Hercules',
'LeoI':'Leo I',
'LeoII':'Leo II',
'LeoIV':'Leo IV',
'LeoV':'Leo V',
'LeoT':'Leo T',
'Sculptor':'Sculptor',
'Segue1':'Segue 1',
'Segue2':'Segue 2',
'Sextans':'Sextans',
'UrsaMajorI':'Ursa Major I',
'UrsaMajorII':'Ursa Major II',
'UrsaMinor':'Ursa Minor',
'ReticulumII': 'Reticulum II',
'TucanaII':'Tucana II',
'HydraII':'Hydra II',
'HorologiumI':'Horologium I',
'PiscesII':'Pisces II',
'GruI':'Grus I',
'Willman1':'Willman 1'}

posh_latex_names= {'BootesI':u'Bo\\"otes I',
'Carina':'Carina',
'Coma':'Coma Berenices',
'CVnI':'Canes Venatici I',
'CVnII':'Canes Venatici II',
'Draco':'Draco',
'Fornax':'Fornax',
'Hercules':'Hercules',
'LeoI':'Leo I',
'LeoII':'Leo II',
'LeoIV':'Leo IV',
'LeoV':'Leo V',
'LeoT':'Leo T',
'Sculptor':'Sculptor',
'Segue1':'Segue 1',
'Segue2':'Segue 2',
'Sextans':'Sextans',
'UrsaMajorI':'Ursa Major I',
'UrsaMajorII':'Ursa Major II',
'UrsaMinor':'Ursa Minor',
'ReticulumII': 'Reticulum II',
'TucanaII':'Tucana II',
'HydraII':'Hydra II',
'HorologiumI':'Horologium I',
'PiscesII':'Pisces II',
'GruI':'Grus I',
'Willman1':'Willman 1'}

bonnivard_names = {'BootesI':'boo1',
'Carina':'car',
'Coma':'coma',
'CVnI':'cvn1',
'CVnII':'cvn2',
'Draco':'dra',
'Fornax':'for',
'Hercules':'her',
'LeoI':'leo1',
'LeoII':'leo2',
'LeoIV':'leo4',
'LeoV':'leo5',
'LeoT':'leot',
'Sculptor':'scl',
'Segue1':'seg1',
'Segue2':'seg2',
'Sextans':'sex',
'UrsaMajorI':'uma1',
'UrsaMajorII':'uma2',
'UrsaMinor':'umi',
'Willman1':'wil1'}

def read_bonnivard_table(Name):
	''' Reads annihilation data from Bonnivard (2015) '''
	GEV2cm5toMsol2kpc5 = 2.2482330e-07
	if Name in bonnivard_names:
		data = np.genfromtxt('../bonnivard/'+bonnivard_names[Name]+'_Jalphaint_cls.output',
			     skip_header=5)
		data = np.delete(data,[2,5],1)
		df = pd.DataFrame(data,columns=['alpha','J','eJm68','eJp68','eJm95','eJp95'])
		df['J']=np.log10(df['J']/GEV2cm5toMsol2kpc5)
		df['eJm68']=np.log10(df['eJm68']/GEV2cm5toMsol2kpc5)
		df['eJp68']=np.log10(df['eJp68']/GEV2cm5toMsol2kpc5)
		df['eJm95']=np.log10(df['eJm95']/GEV2cm5toMsol2kpc5)
		df['eJp95']=np.log10(df['eJp95']/GEV2cm5toMsol2kpc5)
		return df
	else:
		return pd.DataFrame()

def read_bonnivard_table_decay(Name):
	''' Reads decay data from Bonnivard (2015)'''
	GEVcm2toMsolkpc2 = 8.5358230e-15
	if Name in bonnivard_names:
		data = np.genfromtxt('../bonnivard/'+bonnivard_names[Name]+'_Dalphaint_cls.output',
			     skip_header=5)
		data = np.delete(data,[2,5],1)
		df = pd.DataFrame(data,columns=['alpha','D','eDm68','eDp68','eDm95','eDp95'])
		df['D']=np.log10(df['D']/GEVcm2toMsolkpc2 )
		df['eDm68']=np.log10(df['eDm68']/GEVcm2toMsolkpc2)
		df['eDp68']=np.log10(df['eDp68']/GEVcm2toMsolkpc2)
		df['eDm95']=np.log10(df['eDm95']/GEVcm2toMsolkpc2)
		df['eDp95']=np.log10(df['eDp95']/GEVcm2toMsolkpc2)
		return df
	else:
		return pd.DataFrame()

def read_ackermann_data():
	''' Reads data from the Ackermann Fermi-LAT paper '''
	names = np.genfromtxt('../ackermann/ackermann_dwarfs.dat',skip_header=2,usecols=0,dtype=str)
	data = np.genfromtxt('../ackermann/ackermann_dwarfs.dat',skip_header=2)[:,4:6]
	df = pd.DataFrame(data,columns=['J','eJ'])
	df['name']=names
	return df

def make_table(data,geo_factor=True):
	''' Outputs two tables of J- and D-factors for the dwarfs using the NFW
		formula with rs = 5 R_half.
		geo_factor multiplies the half-light radii by a factor sqrt(1-e) to
		ellipticity correct them for use in the spherical formulae '''
	geof = np.ones(len(data))
	if(geo_factor):
		geof = np.sqrt(1.-data['ellip'])
	rnfwrs = 5.
	N=100000
	WEJ2 = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half_05',N=N,
	                                  nfw=rnfwrs*data['R_half']*geof/1000.,
	                                  geo_factor=geo_factor)
	WED2 = wyns_formulaD_error_sample(data,gamma=1.,angle='Half_05',N=N,
	                                  nfw=rnfwrs*data['R_half']*geof/1000.,
	                                  geo_factor=geo_factor)
	WEJ3 = wyns_formulaJ_error_sample(data,gamma=1.,angle='Max',N=N,
	                                  nfw=rnfwrs*data['R_half']*geof/1000.,
	                                  geo_factor=geo_factor)
	WED3 = wyns_formulaD_error_sample(data,gamma=1.,angle='Max',N=N,
	                                  nfw=rnfwrs*data['R_half']*geof/1000.,
	                                  geo_factor=geo_factor)

	outfile=open('dwarfs_Jfactors.dat','w')
	outfile.write('\\begin{tabular}{llccccc}\n')
	outfile.write('\\hline\n\\hline\n')
	outfile.write('Name & Distance & $\\theta_\mathrm{max}$ & $\log_{10} J(\\theta_\mathrm{max})$ & $\log_{10} J(0.5^\circ)$ & $\log_{10} D(\\theta_\mathrm{max})$ & $\log_{10} D(0.5^\circ)$\\\\ \n')
	outfile.write('&[$\mathrm{kpc}$]& [$^\circ$] & [$\mathrm{GeV^2\,cm}^{-5}$] & [$\mathrm{GeV^2\,cm}^{-5}$] & [$\mathrm{GeV\,cm}^{-2}$] & [$\mathrm{GeV\,cm}^{-2}$]\\\\\n')
	outfile.write('\\hline\n')
	for i in range(len(WEJ2)):
		string= posh_latex_names[data['Name'][i]]+\
				"&$%0.0f\pm%0.0f$"%(data['D'][i],data['eD'][i])+" & $"+\
				str(data['theta_max'][i])+"$&"+\
				"$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WEJ3[i][0],WEJ3[i][1],WEJ3[i][2])
		if(i>22):
			string+="-&"
		else:
			string+="$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WEJ2[i][0],WEJ2[i][1],WEJ2[i][2])
		string+="$%0.2f_{-%0.2f}^{+%0.2f}$&"%(WED3[i][0],WED3[i][1],WED3[i][2])
		if(i>22):
			string+="-"+"\\\\\n"
		else:
			string+="$%0.2f_{-%0.2f}^{+%0.2f}$"%(WED2[i][0],WEJ3[i][1],WEJ3[i][2])+"\\\\\n"
		if(i==8 or i==23):
			outfile.write('\\hline\n')
		outfile.write(string)
	outfile.write('\\hline\n')
	outfile.write('\end{tabular}\n')
	outfile.close()

	outfile=open('dwarfs_Jfactors_ascii.dat','w')
	outfile.write('#Name D eD thetamax Jmax eJmax1 eJmax2 J05 eJ051 eJ052 Dmax eDmax1 eDmax2 D05 eD051 eD052\n')
	for i in range(len(WEJ2)):
		string= data['Name'][i]+\
				" %0.0f %0.0f "%(data['D'][i],data['eD'][i])+\
				str(data['theta_max'][i])+" "+\
				"%0.2f %0.2f %0.2f "%(WEJ3[i][0],WEJ3[i][1],WEJ3[i][2])
		if(i>22):
			string+="%0.2f %0.2f %0.2f "%(WEJ3[i][0],WEJ3[i][1],WEJ3[i][2])
		else:
			string+="%0.2f %0.2f %0.2f "%(WEJ2[i][0],WEJ2[i][1],WEJ2[i][2])
		string+="%0.2f %0.2f %0.2f "%(WED3[i][0],WED3[i][1],WED3[i][2])
		if(i>22):
			string+="%0.2f %0.2f %0.2f\n"%(WED3[i][0],WED3[i][1],WED3[i][2])
		else:
			string+="%0.2f %0.2f %0.2f\n"%(WED2[i][0],WEJ3[i][1],WEJ3[i][2])
		outfile.write(string)
	outfile.close()

def add_thetas(ax,xrang,thetalist):
	''' Add theta values to a plot '''
	ylim=ax.get_ylim()
	ax.set_ylim(ylim[0]-0.5,ylim[1])
	for x,t in zip(xrang,thetalist):
		ax.annotate(str(t)+r'$^\circ$',xy=(x,ylim[0]),horizontalalignment='center',verticalalignment='bottom',rotation=90)

def summary_data_plot():
	''' Makes plots of data -- half-light radii, velocity dispersions and distances, along with J-factor estimates from the various methods '''

	gs_gammas=np.genfromtxt('geringer_sameth_gamma.dat',skip_header=49)

	cd=data[data.Class=='CD']
	uf=data[data.Class=='UF']
	labelrange=np.linspace(0.,len(data),len(data))
	labelscd=labelrange[:len(cd)]
	labelsuf=labelrange[len(cd):]
	f,a=plt.subplots(2,4,figsize=(16,8))
	plt.subplots_adjust(hspace=0.5)
	for ai in a:
		for aj in ai:
			aj.set_xticks(labelrange)
			aj.set_xticklabels(data.Name.values,rotation=90)
			aj.set_xlim(labelrange[0]-1,labelrange[-1]+1)

	for i in a[1]:
		ls=i.axvline(labelscd[-1]+.5,c='k',ls='dashed')
		ls.set_dashes((2,1))
		ls=i.axvline(labelsuf[13]+.5,c='k',ls='dashed')
		ls.set_dashes((2,1))

	a[0][0].errorbar(labelscd,cd.D,yerr=cd.eD,fmt='.')
	a[0][0].errorbar(labelsuf,uf.D.values,yerr=uf.eD.values,fmt='.')
	a[0][0].set_ylabel(r'Distance/kpc')
	a[0][1].errorbar(labelscd,cd.R_half,yerr=[cd.eR_half2,cd.eR_half1],fmt='.')
	a[0][1].errorbar(labelsuf,uf.R_half,yerr=[uf.eR_half2,uf.eR_half1],fmt='.')
	a[0][1].set_ylabel(r'$R_{\mathrm{half}}/\mathrm{pc}$')
	a[0][2].errorbar(labelscd,cd.sigma_los,yerr=[cd.esigma_los2,cd.esigma_los1],fmt='.')
	a[0][2].errorbar(labelsuf,uf.sigma_los,yerr=[uf.esigma_los2,uf.esigma_los1],fmt='.')
	a[0][2].arrow(labelsuf[9],uf.sigma_los.values[9],0.,-0.5,fc=sns.color_palette()[1],ec=sns.color_palette()[1],head_length=0.2,head_width=0.3)
	a[0][2].arrow(labelsuf[15],uf.sigma_los.values[15],0.,-0.5,fc=sns.color_palette()[1],ec=sns.color_palette()[1],head_length=0.2,head_width=0.3)
	a[0][2].arrow(labelsuf[17],uf.sigma_los.values[17],0.,-0.5,fc=sns.color_palette()[1],ec=sns.color_palette()[1],head_length=0.2,head_width=0.3)
	a[0][2].set_ylabel(r'$\sigma_{\mathrm{los}}/\mathrm{km\,s}^{-1}$')

	a[1][0].errorbar(labelscd,cd.Jmax,yerr=[cd.eJmax2,cd.eJmax1],fmt='.',color='k')
	a[1][0].errorbar(labelsuf,uf.Jmax,yerr=[uf.eJmax2,uf.eJmax1],fmt='.',color='k')
	WE = wyns_formulaJ_error_sample(data,gamma=1.)
	for i in range(len(data)):
		a[1][0].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2])

	WE = wyns_formulaJ_error_sample(data,gamma=0.51)
	for i in range(len(data)):
		a[1][0].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4])

	WE = wyns_formulaJ_error_sample(data,gamma=1.,nfw=5.*data['R_half']/1000.)
	for i in range(len(data)):
		a[1][0].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0])

	add_thetas(a[1][0],labelrange,data.theta_max)
	a[1][0].set_ylabel(r'$\log_{10}(J_\mathrm{max}/\,\mathrm{GeV^2\,cm}^{-5})$')

	a[1][1].errorbar(labelscd,cd.Jmax.values-np.log10(2.),yerr=[cd.eJmax2,cd.eJmax1],fmt='.',label="",color='k')
	a[1][1].errorbar(labelsuf,uf.Jmax.values-np.log10(2.),yerr=[uf.eJmax2,uf.eJmax1],fmt='.',label="",color='k')
	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half')
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma=1$'
		a[1][1].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2],label=label)

	WE = wyns_formulaJ_error_sample(data,gamma=0.51,angle='Half')
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma=0.51$'
		a[1][1].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4],label=label)

	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half',nfw=5.*data['R_half']/1000.)
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'NFW'
		a[1][1].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0],label=label)

	gammas = gs_gammas.T[23]
	while(len(gammas)<len(data)):
		gammas = np.append(gammas,0.8)

	WE = wyns_formulaJ_error_sample(data,gammaarray=gammas,angle='Half')
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma_\mathrm{GS}$'
		a[1][1].fill_between([labelrange[i]-0.3,labelrange[i]+0.3], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=1.,facecolor="None",label=label)
	add_thetas(a[1][1],labelrange,data.theta_half)
	a[1][1].legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.0))
	a[1][1].set_ylabel(r'$\log_{10}(J_\mathrm{half}/\,\mathrm{GeV^2\,cm}^{-5})$')

	a[1][2].errorbar(labelscd,cd.dJmax.values-np.log10(2.),yerr=[cd.eJmax2,cd.edJmax1],fmt='.',color='k')
	a[1][2].errorbar(labelsuf,uf.dJmax.values-np.log10(2.),yerr=[uf.edJmax2,uf.edJmax1],fmt='.',color='k')

	WE = wyns_formulaD_error_sample(data,gamma=1.)
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma=1.$'
		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2],label=label)

	WE = wyns_formulaD_error_sample(data,gamma=1.49)
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma=1.49$'
		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4],label=label)

	WE = wyns_formulaD_error_sample(data,gamma=1.,nfw=5.*data['R_half']/1000.)
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'NFW'
		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0],label=label)

	WE = wyns_formulaD_error_sample(data,gammaarray=gammas,angle='Half')
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma_\mathrm{GS}$'
		a[1][2].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=1.,facecolor="None",label=label)
	add_thetas(a[1][2],labelrange,data.dtheta_half)
	a[1][2].legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.0))
	a[1][2].set_ylabel(r'$\log_{10}(D_\mathrm{half}/\,\mathrm{GeV\,cm}^{-2})$')

	a[1][3].errorbar(labelscd,cd.Jhalf.values,yerr=[cd.eJhalf2,cd.eJhalf1],fmt='.',label="",color='k')
	a[1][3].errorbar(labelsuf,uf.Jhalf.values,yerr=[uf.eJhalf2,uf.eJhalf1],fmt='.',label="",color='k')
	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half'	)
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma=1$'
		a[1][3].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[2],label=label)

	WE = wyns_formulaJ_error_sample(data,gamma=0.51,angle='Half_05')
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma=0.51$'
		a[1][3].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[4],label=label)

	WE = wyns_formulaJ_error_sample(data,gamma=1.,angle='Half_05',nfw=5.*data['R_half']/1000.)
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'NFW'
		a[1][3].fill_between([labelrange[i]-0.2,labelrange[i]+0.2], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=0.5,edgecolor="None",color=sns.color_palette()[0],label=label)

	gammas = gs_gammas.T[23]
	while(len(gammas)<len(data)):
		gammas = np.append(gammas,0.8)

	WE = wyns_formulaJ_error_sample(data,gammaarray=gammas,angle='Half_05')
	for i in range(len(data)):
		label=None
		if(i==0):
			label=r'$\gamma_\mathrm{GS}$'
		a[1][3].fill_between([labelrange[i]-0.3,labelrange[i]+0.3], [WE[i][0]-WE[i][1],WE[i][0]-WE[i][1]], [WE[i][0]+WE[i][2],WE[i][0]+WE[i][2]],alpha=1.,facecolor="None",label=label)
	add_thetas(a[1][3],labelrange,np.ones(0.5)*len(data))
	a[1][3].legend(loc="lower center",ncol=2, bbox_to_anchor=(0.5, 1.0))
	a[1][3].set_ylabel(r'$\log_{10}(J(0.5^\circ)/\,\mathrm{GeV^2\,cm}^{-5})$')

	plt.savefig('dwarfs_data.pdf',bbox_inches='tight')



if __name__ == '__main__':
	data = pd.read_csv('data.dat',sep=' ')
	make_table(data)
	exit()
