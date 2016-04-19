## Make Table 1 of Sanders, Evans & Geringer-Sameth
## ============================================================================
import sys
sys.path.append('../../code/m2m/')
import ret_2
import spherical_Jfactors as sJ
import numpy as np
## ============================================================================
## Data on Reticulum II
## ============================================================================
import ret_22
RetII = ret_2.RetII
r_maj_exp = RetII.r_maj_exp ## arcmin
r_maj = RetII.r_maj #arcmin
Distance = RetII.Distance #kpc
Velocity_dispersion = RetII.Velocity_dispersion #km/s
rh= (r_maj/60./180.*np.pi)*Distance ## in units kpc
q = 1-.RetII.e
## ============================================================================
## List of models = (spherical model, oblate and three prolate models)
## ============================================================================
models = [ret_2.stellar_halo_model('ret2_roundNFW_roundPlummer',
                                   r'\makecell[l]{Spherical NFW\\ Spherical Plummer}',
                                   baN=1.,caN=1.,
                                   alpha=2.,beta=5.,gamma=0.,
                                   rnfwrs=2.,rtrunc=10.,bovera=1.),
		  ret_2.stellar_halo_model('ret2_flatNFW_flatPlummer_obl'+'_%0.1f'%(q),
	                            r'\makecell[l]{Oblate NFW, $p=1,\,q=%0.1f$\\ Oblate Plummer, $p=1,\,q=%0.1f$}' % (q,q),
		                            baN=q,caN=q,
		                            alpha=2.,beta=5.,gamma=0.,
		                            rnfwrs=2.,rtrunc=10.,bovera=1.),
		  ret_2.stellar_halo_model('ret2_flatNFW_flatPlummer_pro'+'_%0.1f'%(q),
	                        r'\makecell[l]{Prolate NFW, $p=1,\,q=%0.1f$\\ Prolate Plummer, $p=1,\,q=%0.1f$}' % (1./q,1./q),
		                            baN=q,caN=q,
		                            alpha=2.,beta=5.,gamma=0.,
		                            rnfwrs=2.,rtrunc=10.,bovera=q),
		  ret_2.stellar_halo_model('ret2_ba04cusp_NFW_0.5_0.4',
		                    r'\makecell[l]{Near-prolate NFW, $p=0.5$, $q=0.4$\\ Near-prolate cuspy Plummer\\\quad$\alpha_\star=2,\,\beta_\star=5,\,\gamma_\star=1,\,p=0.4,\,q=0.38$}',
		                           baN=0.5,caN=0.4,
		                           alpha=2.,beta=5.,gamma=1.,
		                           rnfwrs=2.,rtrunc=10.,bovera=0.4),
		  ret_2.stellar_halo_model('ret2_g0b4_Plummer',
		                    r'\makecell[l]{Prolate cored DM\\\quad$\alpha_\mathrm{DM}=1,\,\beta_\mathrm{DM}=4,\,\gamma_\mathrm{DM}=0,\,p=1,\,q=2.5$\\ Prolate Plummer, $p=1,\,q=2.5$}',
		                           baN=0.4,caN=0.4,
		                           alpha=2.,beta=5.,gamma=0.,
		                           inner_halo=0.,outer_halo=4.)]
## ============================================================================
## Analytic formulae -- note geometric correction factor sqrt(1-e)
## ============================================================================
sph_JW = sJ.wyns_formulaJ_NFW_data(Velocity_dispersion,rh*1000.*np.sqrt(q),Distance,[0.5],rh*2.*np.sqrt(q))
sph_DW = sJ.wyns_formulaD_NFW_data(Velocity_dispersion,rh*1000.*np.sqrt(q),Distance,[0.5],rh*2.*np.sqrt(q))
print sph_JW,sph_DW
## ============================================================================
## Spherical model -- note geometric correction factor sqrt(1-e)
## ============================================================================
FF = ret_2.find_J_value(models[0],los='z',gf=np.sqrt(q))
sph_J = np.log10(FF[-2])
sph_D = np.log10(FF[-1])
print sph_J,sph_D
## ============================================================================
## Generate table
## ============================================================================
outfile = open('ret2_table.dat','w')
## Header
outfile.write('\\begin{tabular}{lcccccc}\n')
outfile.write('\\hline\n\\hline\n')
outfile.write('Model & Paper I J & Paper I D & $\log_{10}(\mathrm{J}(0.5^\circ))$& $\log_{10}(\mathrm{D}(0.5^\circ))$& $\mathcal{F}_\mathrm{J}$& $\mathcal{F}_\mathrm{D}$\\\\ \n')
outfile.write('\\hline\n')
## Spherical model
string =models[0].descr+"&$%0.2f$&$%0.2f$"%(sph_JW,sph_DW)
string+="&$%0.2f$&$%0.2f$&$%0.2f$&$%0.2f$\\\\\n\\\\\n"%(sph_J,sph_D,0.,0.)
outfile.write(string)
## Oblate model
FF = ret_2.find_J_value(models[1],los='x',flattening=q,obl_or_prol='obl')
string =models[1].descr+"&-&-"
string+="&$%0.2f$&$%0.2f$&$%0.2f$&$%0.2f$\\\\\n\\\\\n"%(np.log10(FF[-2]),np.log10(FF[-1]),np.log10(FF[-2])-sph_J,np.log10(FF[-1])-sph_D)
outfile.write(string)
## Prolate models
for i in models[2:]:
	FF = ret_2.find_J_value(i,los='y',flattening=q,obl_or_prol='prol')
	string =i.descr+"&-&-"
	string+="&$%0.2f$&$%0.2f$&$%0.2f$&$%0.2f$\\\\\n\\\\\n"%(np.log10(FF[-2]),np.log10(FF[-1]),np.log10(FF[-2])-sph_J,np.log10(FF[-1])-sph_D)
	outfile.write(string)
## Footer
outfile.write('\\hline\n')
outfile.write('\\end{tabular}\n')
outfile.close()
## ============================================================================
