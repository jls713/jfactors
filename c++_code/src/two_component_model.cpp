// ============================================================================
#include <Python.h>
#include "utils_j.h"
#include <cmath>
#include <vector>
#include <string>
#include "../general/utils.h"
#include "density.h"
#include "observed_density.h"
#include "potential.h"
#include "Multipole.h"
#include "two_component_model.h"

//=============================================================================
// Utility functions
//=============================================================================

double qpot_from_q(double q){
  return .5*sqrt(1.+sqrt(1.+8.*q*q));
}

double wolf(double slos, double rh){
  return 4.*slos*slos*rh/Grav;
}

double wolf_sph(double slos, double rhs){
  return 3.*slos*slos*rhs/Grav;
}

//=============================================================================
// Specific model for the paper
//=============================================================================
DoubleProfileModel::DoubleProfileModel(double ba, double ca, double rh, double slos, bool use_multipole, double rsrdmratio, double rtdmrdmratio, double rtsrdmratio, VecDoub abg_st, VecDoub abg_dm)
  :rh(rh),slos(slos),use_multipole(use_multipole),rsrdmratio(rsrdmratio),rtdmrdmratio(rtdmrdmratio), rtsrdmratio(rtsrdmratio), abg_st(abg_st), abg_dm(abg_dm){
  VecDoub abc_st = {1.,ba,ca};
  VecDoub abc_dm = {1.,ba,ca};
  double rs_dm = 1., rs_st = rsrdmratio, rt_st = rtsrdmratio, rt_dm=rtdmrdmratio;
  double rho0_st = 1., rho0_dm = 1.;
  Stars = new AlphaBetaGammaDensityProfile(abg_st,rho0_st,rs_st,rt_st,abc_st,true);
  DM = new AlphaBetaGammaDensityProfile(abg_dm,rho0_dm,rs_dm,rt_dm,abc_dm,true);
}
void DoubleProfileModel::print(void){
  std::cout<<"=====\nStellar properties:\n";
  Stars->print();
  std::cout<<"=====\nDark-matter properties:\n";
  DM->print();
  std::cout<<"=====\n";
}
double DoubleProfileModel::observed_ellipticity(double theta, double phi){
  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);
  return ObsStars.ellipticity();
}
VecDoub DoubleProfileModel::J_factor_old(double theta, double phi, double D, double ang, bool gobby, bool with_D){

  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);
  ObservedTriaxialDensityProfile ObsDM(DM,theta,phi);

  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,DM->scale_radius(),0.001*DM->scale_radius(),10.*DM->tidal_radius());

  double Stars_Reff = ObsStars.half_light_radius();
  double KinematicRatio = (use_multipole?ObsStars.sigma_corr(&MEA):ObsStars.sigma_corr(&NFWP));;

  // As DM and stars flattened in same way we can cheat
  double Stars_rhalf = Stars->spherical_half_light_radius();
  double DM_ellipsoid = DM->mass_ellipsoid(Stars_rhalf);

  if(gobby){
    double DMmass = DM->mass();
    double StarMass = Stars->mass();
    double DMMass_sphere = DM->M_sphere(Stars_Reff,180./PI);
    double DMMass_ellipse = ObsDM.M_ellipse(Stars_Reff);
    double StarsMass_ellipse = ObsStars.M_ellipse(Stars_Reff);
    std::cout<<"Total mass DM:"<<DMmass<<std::endl;
    std::cout<<"Total mass Stars:"<<StarMass<<std::endl;
    std::cout<<"Stars half-light radius:"<<Stars_Reff<<std::endl;
    std::cout<<"Star mass in half-light radius:"<<StarsMass_ellipse<<std::endl;
    std::cout<<"DM mass in half-light radius:"<<DMMass_ellipse<<std::endl;
    std::cout<<"Kinematic Ratio:"<<KinematicRatio<<std::endl;
    std::cout<<"Sigma_tot:"<<slos/sqrt(KinematicRatio)<<std::endl;
  }
  // Now scale
  double RadiusRatio = rh/Stars_Reff;
  double Wolf = wolf(slos/sqrt(KinematicRatio),Stars_rhalf*RadiusRatio/DM->axes()[0]);
  double MassRatio = Wolf/DM_ellipsoid;
  if(gobby){
    std::cout<<"Wolf mass:"<<Wolf<<std::endl;
    std::cout<<"Mass ratio: "<<MassRatio<<", Radius ratio: "<<RadiusRatio<<std::endl;
  }
  double rho0_st = Stars->central_density()/pow(RadiusRatio,3.)*MassRatio;
  double rs_st = Stars->scale_radius()*RadiusRatio;
  double rt_st = Stars->tidal_radius()*RadiusRatio;
  VecDoub abg_st = Stars->alpha_beta_gamma();
  VecDoub abc_st = Stars->axes();
  double rho0_dm = DM->central_density()/pow(RadiusRatio,3.)*MassRatio;
  double rs_dm = DM->scale_radius()*RadiusRatio;
  double rt_dm = DM->tidal_radius()*RadiusRatio;
  VecDoub abg_dm = DM->alpha_beta_gamma();
  VecDoub abc_dm = DM->axes();

  AlphaBetaGammaDensityProfile StarsSc(abg_st,rho0_st,rs_st,rt_st,abc_st,false);
  AlphaBetaGammaDensityProfile DMSc(abg_dm,rho0_dm,rs_dm,rt_dm,abc_dm,false);
  ObservedTriaxialDensityProfile ObsStarsSc(&StarsSc,theta,phi);
  ObservedTriaxialDensityProfile ObsDMSc(&DMSc,theta,phi);

  if(gobby){
    // some checks
    double DMmass = DMSc.mass();
    double StarMass = StarsSc.mass();
    Stars_Reff = ObsStarsSc.half_light_radius();
    double StarsMass_ellipse = ObsStarsSc.M_ellipse(Stars_Reff);
    NFW NFWP2(1.,DMSc.scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
    MultipoleDensity MP2(&DMSc);
    MultipoleExpansion_Triaxial MEA2(&MP2,150,16,12,8,DMSc.scale_radius(),0.001*DMSc.scale_radius(),10.*DMSc.tidal_radius());
    KinematicRatio = (use_multipole?ObsStarsSc.sigma_corr(&MEA2):ObsStarsSc.sigma_corr(&NFWP2));

    // As DM and stars flattened in same way we can cheat
    double DMMass_ellipse = ObsDMSc.M_ellipse(Stars_Reff);
    std::cout<<"Total mass DM:"<<DMmass<<std::endl;
    std::cout<<"Total mass Stars:"<<StarMass<<std::endl;
    std::cout<<"Stars half-light radius:"<<Stars_Reff<<std::endl;
    std::cout<<"Star mass in half-light radius:"<<StarsMass_ellipse<<std::endl;
    std::cout<<"DM mass in half-light radius:"<<DMMass_ellipse<<std::endl;
    std::cout<<"Kinematic Ratio:"<<KinematicRatio<<std::endl;
  }
  double Jfactor = ObsDMSc.J_far_arbitrary_orientation(D,ang);
  if(gobby) std::cout<<"Jfactor: "<<Jfactor<<std::endl;
  if(with_D){
    double Dfactor = ObsDMSc.D_far_arbitrary_orientation(D,ang);
    return {Jfactor,Dfactor};
  }
  else return {Jfactor};
}

VecDoub DoubleProfileModel::J_factor(double theta, double phi, double D, double ang, bool gobby, bool with_D){

  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);

  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,DM->scale_radius(),0.001*DM->scale_radius(),10.*DM->tidal_radius());

  double Stars_Reff = ObsStars.half_light_radius();
  double VelocityDispersion = (use_multipole?ObsStars.sigma_los(&MEA):ObsStars.sigma_los(&NFWP));

  // Now scale
  double RadiusRatio = rh/Stars_Reff;
  double MassRatio = pow(slos/VelocityDispersion,2.)*RadiusRatio;

  double rho0_dm = DM->central_density()/pow(RadiusRatio,3.)*MassRatio;
  double rs_dm = DM->scale_radius()*RadiusRatio;
  double rt_dm = DM->tidal_radius()*RadiusRatio;
  VecDoub abg_dm = DM->alpha_beta_gamma();
  VecDoub abc_dm = DM->axes();

  AlphaBetaGammaDensityProfile DMSc(abg_dm,rho0_dm,rs_dm,rt_dm,abc_dm,false);
  ObservedTriaxialDensityProfile ObsDMSc(&DMSc,theta,phi);

  double Jfactor = ObsDMSc.J_far_arbitrary_orientation(D,ang);
  if(with_D){
    double Dfactor = ObsDMSc.D_far_arbitrary_orientation(D,ang);
    return {Jfactor,Dfactor};
  }
  else return {Jfactor};
}

VecDoub DoubleProfileModel::MassProfile(double theta, double phi, double D, VecDoub ang, bool gobby, std::string typ){

  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);

  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,DM->scale_radius(),0.001*DM->scale_radius(),DM->tidal_radius()>0.?10.*DM->tidal_radius():100.*DM->scale_radius());

  double Stars_Reff = ObsStars.half_light_radius();
  double VelocityDispersion = (use_multipole?ObsStars.sigma_los(&MEA):ObsStars.sigma_los(&NFWP));

  // Now scale
  double RadiusRatio = rh/Stars_Reff;
  double MassRatio = pow(slos/VelocityDispersion,2.)*RadiusRatio;

  double rho0_dm = DM->central_density()/pow(RadiusRatio,3.)*MassRatio;
  double rs_dm = DM->scale_radius()*RadiusRatio;
  double rt_dm = DM->tidal_radius()*RadiusRatio;
  VecDoub abg_dm = DM->alpha_beta_gamma();
  VecDoub abc_dm = DM->axes();

  AlphaBetaGammaDensityProfile DMSc(abg_dm,rho0_dm,rs_dm,rt_dm,abc_dm,false);
  ObservedTriaxialDensityProfile ObsDMSc(&DMSc,theta,phi);

  VecDoub masses(ang.size(),0.);
  if(typ=="cylinder")
    for(unsigned i=0;i<ang.size();++i)
      masses[i]=ObsDMSc.M_ellipse(D*ang[i]/180.*PI);
  else if(typ=="sphere")
    for(unsigned i=0;i<ang.size();++i)
      masses[i]=DMSc.M_sphere(D,ang[i]);
  else if(typ=="ellipsoid")
    for(unsigned i=0;i<ang.size();++i)
      masses[i]=DMSc.mass_ellipsoid(D*ang[i]/180.*PI);
  return masses;
}

double DoubleProfileModel::correction_factor(double theta, double phi, double D, double ang, bool gobby, bool geo_average){
  // Note that this function permanently changes the stars and dark matter
  // profiles to spherical

  double true_J = J_factor(theta,phi,D,ang,gobby)[0];

  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);
  double e = ObsStars.ellipticity();
  if(gobby) std::cout<<"Observed ellipticity: "<<e<<std::endl;
  if(geo_average) rh*=sqrt(e);

  VecDoub abg_st = Stars->alpha_beta_gamma(); // Plummer
  VecDoub abg_dm = DM->alpha_beta_gamma(); // NFW
  VecDoub abc_st = {1.,1.,0.999};
  VecDoub abc_dm = {1.,1.,0.999};
  double rs_dm = DM->scale_radius(), rs_st = Stars->scale_radius();
  double rt_st = DM->tidal_radius(), rt_dm=Stars->tidal_radius();
  double rho0_st = Stars->central_density(), rho0_dm = DM->central_density();
  delete Stars; delete DM;
  Stars = new AlphaBetaGammaDensityProfile(abg_st,rho0_st,rs_st,rt_st,abc_st,false);
  DM = new AlphaBetaGammaDensityProfile(abg_dm,rho0_dm,rs_dm,rt_dm,abc_dm,false);
  double sph_J = J_factor(theta,phi,D,ang,gobby)[0];
  std::cout<<true_J<<" "<<sph_J<<std::endl;
  return true_J/sph_J;
}
double DoubleProfileModel::kinematic_ratio(double theta, double phi){
  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);
  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,DM->scale_radius(),0.001*DM->scale_radius(),10.*DM->tidal_radius());
  double KinematicRatio = (use_multipole?ObsStars.sigma_corr(&MEA):ObsStars.sigma_corr(&NFWP));
  return KinematicRatio;
}
double DoubleProfileModel::sigma_x_sigma_z(void){
  ObservedTriaxialDensityProfile ObsStars(Stars,0.,0.);
  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion MEA(&MP,150,16,16,8,2.*DM->scale_radius(),0.001*DM->scale_radius(),10.*DM->tidal_radius(),false,true,true);
  double KinematicRatio = (use_multipole?ObsStars.sigma_x_sigma_z(&MEA):ObsStars.sigma_x_sigma_z(&NFWP));
  return KinematicRatio;
}
double DoubleProfileModel::radius_ratio(double theta, double phi){

  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);

  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,DM->scale_radius(),0.001*DM->scale_radius(),10.*DM->tidal_radius());

  double Stars_Reff = ObsStars.half_light_radius();
  return rh/Stars_Reff;
}

double DoubleProfileModel::velocity_ratio(double theta, double phi){

  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);

  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,DM->scale_radius(),0.001*DM->scale_radius(),10.*DM->tidal_radius());

  double VelocityDispersion = (use_multipole?ObsStars.sigma_los(&MEA):ObsStars.sigma_los(&NFWP));

  return slos/VelocityDispersion;
}
double DoubleProfileModel::sigma_tot(void){
  VecDoub ar = DM->axis_ratios();
  NFW NFWP(1.,DM->scale_radius(),qpot_from_q(ar[0]),qpot_from_q(ar[1]));
  MultipoleDensity MP(DM);
  MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,DM->scale_radius(),0.001*DM->scale_radius(),10.*DM->tidal_radius());
  return Stars->sigma_tot(&MEA);
}
double DoubleProfileModel::spherical_rh(void){
  return Stars->spherical_half_light_radius();
}
double DoubleProfileModel::projected_rh(double theta, double phi){
  ObservedTriaxialDensityProfile ObsStars(Stars,theta,phi);
  return ObsStars.half_light_radius();
}
