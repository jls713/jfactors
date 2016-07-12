#include "Python.h"
#include "utils.h"
#include "../general/utils.h"
#include "density.h"
#include "observed_density.h"
#include "potential.h"
#include "Multipole.h"
#include "two_component_model.h"
#include "analytic_results.h"

//=============================================================================
// Programs for paper plots
//=============================================================================

void kinematic_ratio_loop(void){
  std::cout<<"Round\n";
  VecDoub Q = {0.4,0.5,0.6,0.7};
  for(auto q: Q){
    PaperModel DD(1.,q);
    double K = DD.kinematic_ratio(0.,0.);
    std::cout<<q<<" "<<1./(K*K)<<std::endl;
  }
  for(auto q: Q){
    PaperModel DD(q,q);
    double K = DD.kinematic_ratio(PI/2.,0.);
    std::cout<<1./q<<" "<<1./(K*K)<<std::endl;
  }
  std::cout<<"Edge\n";
  for(auto q: Q){
    PaperModel DD(1.,q);
    double K = DD.kinematic_ratio(PI/2.,0.);
    std::cout<<q<<" "<<1./(K*K)<<std::endl;
  }
  for(auto q: Q){
    PaperModel DD(q,q);
    double K = DD.kinematic_ratio(0.,0.);
    std::cout<<1./q<<" "<<1./(K*K)<<std::endl;
  }
}

void sigratio_loop(std::string filename, bool with_multipole=false){
  std::ofstream outfile(filename+".wr.dat");
  VecDoub Q = create_range(0.2,1.,45);
  for(auto q: Q){
    PaperModel DD(1.,q, with_multipole);
    double K = DD.sigma_x_sigma_z();
    outfile<<q<<" "<<K<<std::endl;
  }
  for(auto q: Q){
    PaperModel DD(q,q, with_multipole);
    double K = DD.sigma_x_sigma_z();
    outfile<<1./q<<" "<<1./K<<std::endl;
  }
  outfile.close();
}

void rvratio_loop(std::string filename, bool with_multipole=false){

  double rsrdmratio=0.5; double rtdmrdmratio=10.; double rtsrdmratio=9.;
  VecDoub abg_st={2.,5.,0.};VecDoub abg_dm={1.,3.,1.};
  double Distance = 30.;
  double ang = 0.5;
  double rh = 0.049;
  double slos = 3.22;
  std::ofstream outfile(filename+".rv.dat");
  VecDoub Q = create_range(0.02,1.,45);
  for(auto q: Q){
    DoubleProfileModel DD(1.,q,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    double K = DD.velocity_ratio(PI/2.,0.);
    double K2 = DD.radius_ratio(PI/2.,0.);
    double K25 = DD.sigma_tot();
    double K3 = DD.J_factor(PI/2.,0.,30.,0.5,false)[0];
    outfile<<q<<" "<<K<<" "<<K2<<" "<<K25<<" "<<K3<<" ";
    DoubleProfileModel DD2(1.,1.,sqrt(q)*rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    K = DD2.velocity_ratio(PI/2.,0.);
    K2 = DD2.radius_ratio(PI/2.,0.);
    K25 = DD.sigma_tot();
    K3 = DD2.J_factor(PI/2.,0.,30.,0.5,false)[0];
    outfile<<K<<" "<<K2<<" "<<K25<<" "<<K3<<std::endl;
  }
  for(auto q: Q){
    DoubleProfileModel DD(q,q,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    double K = DD.velocity_ratio(0.,0.);
    double K2 = DD.radius_ratio(0.,0.);
    double K25 = DD.sigma_tot();
    double K3 = DD.J_factor(0.,0.,30.,0.5,false)[0];
    outfile<<1./q<<" "<<K<<" "<<K2<<" "<<K25<<" "<<K3<<" ";
    DoubleProfileModel DD2(1.,1.,sqrt(q)*rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    K = DD2.velocity_ratio(0.,0.);
    K2 = DD2.radius_ratio(0.,0.);
    K25 = DD.sigma_tot();
    K3 = DD2.J_factor(0.,0.,30.,0.5,false)[0];
    outfile<<K<<" "<<K2<<" "<<K25<<" "<<K3<<std::endl;
  }
  outfile.close();
}

void jfactor_loop(bool gobby=false){
  double Distance = 30.;
  double ang = 0.5;
  std::cout<<"Round\n";
  VecDoub Q = {0.4,0.5,0.6,0.7,0.99};
  for(auto q: Q){
    PaperModel DD(1.,q);
    double K = DD.J_factor(0.,0.,Distance,ang,gobby)[0];
    std::cout<<q<<" "<<K<<std::endl;
  }
  for(auto q: Q){
    PaperModel DD(q,q*0.99);
    double K = DD.J_factor(PI/2.,0.01,Distance,ang,gobby)[0];
    std::cout<<1./q<<" "<<K<<std::endl;
  }
  std::cout<<"Edge\n";
  for(auto q: Q){
    PaperModel DD(1.,q);
    double K = DD.J_factor(PI/2.,0.,Distance,ang,gobby)[0];
    std::cout<<q<<" "<<K<<std::endl;
  }
  for(auto q: Q){
    PaperModel DD(q,q*0.999);
    double K = DD.J_factor(0.,0.,Distance,ang,gobby)[0];
    std::cout<<1./q<<" "<<K<<std::endl;
  }
}
void cfactor_loop(std::string filename, bool gobby=false, bool with_multipole=false,bool geo_average=false,double rsrdmratio=0.5, double rtdmrdmratio=10., double rtsrdmratio=9., VecDoub abg_st={2.,5.,0.}, VecDoub abg_dm={1.,3.,1.}){
  double Distance = 30.;
  double ang = 0.5;
  double rh = 0.049;
  double slos = 3.22;
  std::ofstream outfile(filename+".cf.dat");
  VecDoub Q = create_range(0.25,0.99,20);

  DoubleProfileModel DR(1.,0.999,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
  DR.print();
  VecDoub Kr = DR.J_factor(0.,0.,Distance,ang,gobby,true);

  for(auto q: Q){
    DoubleProfileModel DD(1.,q,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    VecDoub K = DD.J_factor(0.,0.,Distance,ang,gobby,true);

    outfile<<q<<" "<<K[0]/Kr[0]<<" "<<K[1]/Kr[1];
    DoubleProfileModel DD2(1.,q,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    K = DD2.J_factor(PI/2.,0.,Distance,ang,gobby,true);

    if(geo_average){
      DoubleProfileModel DRs(1.,0.999,rh*sqrt(q),slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
      VecDoub Krgf = DRs.J_factor(0.,0.,Distance,ang,gobby,true);
      outfile<<" "<<K[0]/Krgf[0]<<" "<<K[1]/Krgf[1]<<std::endl;
    }
    else
      outfile<<" "<<K[0]/Kr[0]<<" "<<K[1]/Kr[1]<<std::endl;
  }
  for(auto q: Q){
    DoubleProfileModel DD(q,q*0.999,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    VecDoub K = DD.J_factor(PI/2.,0.0,Distance,ang,gobby,true);
    outfile<<1./q<<" "<<K[0]/Kr[0]<<" "<<K[1]/Kr[1];
    DoubleProfileModel DD2(q,q*0.999,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
    K = DD2.J_factor(0.,0.,Distance,ang,gobby,true);
    if(geo_average){
      DoubleProfileModel DRs(1.,0.999,rh*sqrt(q),slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
      VecDoub Krgf = DRs.J_factor(0.,0.,Distance,ang,gobby,true);
      outfile<<" "<<K[0]/Krgf[0]<<" "<<K[1]/Krgf[1]<<std::endl;
    }
    else
      outfile<<" "<<K[0]/Kr[0]<<" "<<K[1]/Kr[1]<<std::endl;
  }
  outfile.close();
}

void wyn_models_cfactor_loop(std::string filename, bool geo_average=false,double rsrdmratio=0.5,double n_star = 5.){

  double Distance = 30.;
  double ang = 0.5;
  double rh = 0.049;
  double slos = 3.22;
  std::ofstream outfile(filename+".cf.dat");
  VecDoub Q = create_range(0.25,0.99,20);

  CoredModel CC_sph(n_star,{rsrdmratio,1.},{1.,1.},{1.,1.});
  CC_sph.scale(rh,slos);
  double Jsph = CC_sph.J_factor_face_asymptote(Distance);
  //double Jsph = CC_sph.J_factor_arbitrary(Distance,ang,0.,0.);
  double Dsph = CC_sph.D_factor_arbitrary(Distance,ang,0.,0.);

  for(auto q:Q){
    CoredModel CC(n_star,{rsrdmratio,1.},{q,q},{1.,1.});
    CC.scale(rh,slos,"edge");
    double Je = CC.J_factor_edge(Distance,ang);
    double De = CC.D_factor_arbitrary(Distance,ang,.5*PI,0.);
    double Jes=Jsph, Des=Dsph;
    if(geo_average){
      CoredModel CCs(n_star,{rsrdmratio,1.},{1.,1.},{1.,1.});
      CCs.scale(rh*sqrt(q),slos,"edge");
      Jes = CCs.J_factor_edge(Distance,ang);
      Des = CCs.D_factor_arbitrary(Distance,ang,.5*PI,0.);
    }
    CoredModel CC2(n_star,{rsrdmratio,1.},{q,q},{1.,1.});
    CC2.scale(rh,slos,"round");
    double Jf = CC2.J_factor_face_asymptote(Distance);
    double Df = CC2.D_factor_arbitrary(Distance,ang,0.,0.);
    //double Jf = CC2.J_factor_arbitrary(Distance,ang,0.,0.);
    outfile<<q<<" "<<Jf/Jsph<<" "<<Df/Dsph<<" "<<Je/Jes<<" "<<De/Des<<std::endl;
  }
  for(auto q:Q){
    CoredModel CC(n_star,{rsrdmratio,1.},{1./q,1./q},{1.,1.});
    CC.scale(rh,slos,"edge");
    double Je = CC.J_factor_edge(Distance,ang);
    double De = CC.D_factor_arbitrary(Distance,ang,.5*PI,0.);
    double Jes=Jsph, Des=Dsph;
    if(geo_average){
      CoredModel CCs(n_star,{rsrdmratio,1.},{1.,1.},{1.,1.});
      CCs.scale(rh*sqrt(q),slos,"edge");
      Jes = CCs.J_factor_edge(Distance,ang);
      Des = CCs.D_factor_arbitrary(Distance,ang,.5*PI,0.);
    }
    CoredModel CC2(n_star,{rsrdmratio,1.},{1./q,1./q},{1.,1.});
    CC2.scale(rh,slos,"round");
    double Jf = CC2.J_factor_face_asymptote(Distance);
    double Df = CC2.D_factor_arbitrary(Distance,ang,0.,0.);
    outfile<<1./q<<" "<<Jf/Jsph<<" "<<Df/Dsph<<" "<<Je/Jes<<" "<<De/Des<<std::endl;
  }
  outfile.close();
}

int main(){
  // wyn_models_cfactor_loop("data/wyn_factor_n4_rr2",true,0.5,4.);
  // wyn_models_cfactor_loop("data/wyn_factor_n4_rr20",true,0.05,4.);
  // wyn_models_cfactor_loop("data/wyn_factor_n4_rr200",true,0.005,4.);

  // wyn_models_cfactor_loop("data/wyn_factor_n6_rr2",true,0.5,6.);
  // wyn_models_cfactor_loop("data/wyn_factor_n6_rr20",true,0.05,6.);
  // wyn_models_cfactor_loop("data/wyn_factor_n6_rr200",true,0.005,6.);

  wyn_models_cfactor_loop("data/wyn_factor_n5_rr2",true,0.5,5.);
  wyn_models_cfactor_loop("data/wyn_factor_n5_rr20",true,0.05,5.);
  wyn_models_cfactor_loop("data/wyn_factor_n5_rr200",true,0.005,5.);
  cfactor_loop("fig7_gf",false,true,true,.5,10.,9.,{2.,5.,0.},{1.,3.,1.});
  return 0;
  rvratio_loop("data/test",true); return 0;
  cfactor_loop("fig7_gf_betadm6",false,true,true,.5,10.,9.,{2.,5.,0.},{1.,6.,1.});
  // cfactor_loop("fig7_gf_cored",false,true,true,.5,10.,9.,{2.,5.,0.},{1.,3.,0.});
  return 0;
  cfactor_loop("fig7_gf_beta7",false,true,true,.5,10.,9.,{2.,7.,0.},{1.,3.,1.});return 0;

  cfactor_loop("fig7_ffffff2",false,true,true,2.,10.,9.,{2.,5.,0.},{1.,3.,1.});return 0.;
  // cfactor_loop("fig7_beta4",false,true,true,0.5,10.,9.,{2.,4.,0.},{1.,3.,0.});return 0.;
  // cfactor_loop("fig7_beta35",false,true,true,0.5,10.,9.,{2.,3.5,0.},{1.,3.,0.});
  cfactor_loop("fig7_betadm5",false,true,true,0.5,10.,9.,{2.,5.,0.},{1.,5.,1.});return 0.;
  // // cfactor_loop("fig7_dmcore",true,true,false,1.,9.,10.,{2.,5.,0.},{1.,5.,0.});
  // // cfactor_loop("fig7_dmbeta5",false,true,false,.5,9.,10.,{2.,5.,0.},{1.,5.,1.});
  // // cfactor_loop("fig7_rst5",false,true,false,.5,5.,10.,{2.,5.,0.},{1.,3.,1.});
  // // cfactor_loop("fig7_dmalpha2",false,true,false,.5,5.,10.,{2.,5.,0.},{2.,3.,1.});
  // // cfactor_loop("fig7_r02",false,true,false,.7);
  // return 0;
  cfactor_loop("fig7_gf", false,true,true); return 0;
  cfactor_loop("fig7_r01", false,true,false, 0.1);
  cfactor_loop("fig7_n4", false,true,false, 0.5, 4.);
}
