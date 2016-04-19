//=============================================================================
#include "observed_density.h"
//=============================================================================
// Observed density profiles
// -- Arbitrary density profiles viewed at (theta,phi) spherical polar viewing
//    angles
//=============================================================================

ObservedDensityProfile::ObservedDensityProfile(DensityProfile *DP, double theta, double phi):DP(DP),theta(theta),phi(phi){
  double cp = cos(phi), sp = sin(phi);
  double ct = cos(theta), st = sin(theta);
  n = {st*cp,st*sp,ct};
  phihat={-sp,cp,0.};
  thetahat={-cp*ct,-ct*sp,st};
}

ObservedTriaxialDensityProfile::ObservedTriaxialDensityProfile(FiniteMassTriaxialDensityProfile *TDP, double theta, double phi)
  :ObservedDensityProfile(TDP,theta,phi),DP(TDP){
    VecDoub ar = TDP->axis_ratios();
    double ba = ar[0], ca = ar[1];
    e = observed_ellipticity(ba,ca);
    theta_min=minor_axis_angle(ba,ca);
    stm=sin(theta_min);
    ctm=cos(theta_min);
  }

VecDoub ObservedDensityProfile::x_proj(VecDoub X){
  return phihat*X[0]+thetahat*X[1]+n*X[2];
}

double ObservedTriaxialDensityProfile::observed_ellipticity(double ba, double ca){
// Contopoulos 1956, Weijmans 2014 Appendix
    if(ca==1. and ba==1.) return 1.;
    double st= sin(theta),ct= cos(theta),sp=sin(phi),cp=cos(phi);
    double ba2=ba*ba, ca2=ca*ca, m1ba2 = 1.-ba2,m1ca2 = 1.-ca2;
    double A = m1ca2*ct*ct+m1ba2*st*st*sp*sp+ba2+ca2;
    double B = m1ca2*ct*ct-m1ba2*st*st*sp*sp-ba2+ca2;
    B*=B;
    B+=4.*m1ca2*m1ba2*st*st*ct*ct*sp*sp;
    B = sqrt(B);
    return sqrt((A-B)/(A+B));
}

double ObservedTriaxialDensityProfile::minor_axis_angle(double ba, double ca){
  // Weijmans 2014 Appendix
  double T = (1.-ba*ba)/(1.-ca*ca);
  if(ba==1. and ca==1.) T=0.;
  if(T==0.) return 0.;
  return .5*atan2((2.*T*sin(phi)*cos(phi)*cos(theta)),(sin(theta)*sin(theta)-T*cos(phi)*cos(phi)-sin(phi)*sin(phi)*cos(theta)*cos(theta)));
}

int projected_M_integrand(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    double y2[3];
    obs_mass_density_st *P = (obs_mass_density_st *) fdata;
    for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

    double rs = P->OTDP->scale_radius();
    y2[0]=rs*tan(y2[0]);

    double xp = y2[2]*cos(y2[1]);
    double yyp = y2[2]*sin(y2[1])*P->OTDP->ellipticity();

    // Rotate
    VecDoub ra = P->OTDP->rot_angles();
    double ctm=ra[0],stm=ra[1];
    double x = ctm*xp+stm*yyp;
    double yy =-stm*xp+ctm*yyp;

    // system coordinates
    VecDoub X = P->OTDP->x_proj({x,yy,y2[0]});
    double r = P->OTDP->rho(X);
    double jac = rs*(1.+y2[0]*y2[0]/rs/rs);
    fval[0]=y2[2]*r*jac;
    return 0;
}

double ObservedTriaxialDensityProfile::M_ellipse(double x){
  // mass in an elliptical cylinder of radius x
  VecDoub x2min = {-.5*PI,0.,0.};
  VecDoub x2max = {.5*PI,2.*PI,x};
  obs_mass_density_st P(x2min,x2max,this,0.);
  double err;
  return ellipticity()*integrate(&projected_M_integrand,&P,1e-3,0,INTEG,&err);
}

int projected_J_integrand(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){

    double y2[3];
    obs_density_st *P = (obs_density_st *) fdata;
    for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

    double rs = P->ODP->scale_radius();
    y2[0]=rs*tan(y2[0]);

    double x = y2[2]*cos(y2[1]);
    double yy = y2[2]*sin(y2[1]);

    // system coordinates
    VecDoub X = P->ODP->x_proj({x,yy,y2[0]});
  	double r = P->ODP->rho(X);
    double jac = rs*(1.+y2[0]*y2[0]/rs/rs);
    fval[0]=y2[2]*r*r*jac;
    return 0;
}

double ObservedTriaxialDensityProfile::J_far_arbitrary_orientation(double D, double ang){
	ang = ang/180.*PI;
	VecDoub x2min = {-.5*PI,0.,0.};
	VecDoub x2max = {.5*PI,2.*PI,ang*D};
	obs_density_st P(x2min,x2max,this,D);
	double err;
	return integrate(&projected_J_integrand,&P,1e-3,0,INTEG,&err)/D/D;
}

int projected_D_integrand(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){

    double y2[3];
    obs_density_st *P = (obs_density_st *) fdata;
    for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

    double rs = P->ODP->scale_radius();
    y2[0]=rs*tan(y2[0]);

    double x = y2[2]*cos(y2[1]);
    double yy = y2[2]*sin(y2[1]);

    // system coordinates
    VecDoub X = P->ODP->x_proj({x,yy,y2[0]});
    double r = P->ODP->rho(X);
    double jac = rs*(1.+y2[0]*y2[0]/rs/rs);
    fval[0]=y2[2]*r*jac;
    return 0;
}

double ObservedTriaxialDensityProfile::D_far_arbitrary_orientation(double D, double ang){
  ang = ang/180.*PI;
  VecDoub x2min = {-.5*PI,0.,0.};
  VecDoub x2max = {.5*PI,2.*PI,ang*D};
  obs_density_st P(x2min,x2max,this,D);
  double err;
  return integrate(&projected_D_integrand,&P,1e-3,0,INTEG,&err)/D/D;
}

struct hrf_st{
  ObservedTriaxialDensityProfile *OABGDP;
  double Mtot;
  hrf_st(ObservedTriaxialDensityProfile *OABGDP, double Mtot)
    :OABGDP(OABGDP), Mtot(Mtot){}
};

double hrf_fn(double x, void *p){
  hrf_st *P = (hrf_st *) p;
  double M = P->OABGDP->M_ellipse(x*x);
  return 0.5*P->Mtot-M;
}

double ObservedTriaxialDensityProfile::half_light_radius(void){
  double Mtot = DP->mass();
  root_find RF(1e-3,100);
  hrf_st HRF(this,Mtot);
  double rs = DP->scale_radius();
  VecDoub ar = DP->axis_ratios();
  double ca = ar[1];
  if(ca==1.) ca = 0.999;
  double lower = sqrt(rs*ca), upper = sqrt(rs/ca);
  while(hrf_fn(lower,&HRF)<0.) lower *= 0.8;
  while(hrf_fn(upper,&HRF)>0.) upper *= 1.2;
  double result = RF.findroot(&hrf_fn,lower,upper,&HRF);
  return result*result;
}
//=============================================================================
// Kinematic properties
// -- Ratios of velocity dispersions from tensor virial
//=============================================================================

int sig_integrand(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    sig_st *P = (sig_st *) fdata; VecDoub y2(3,0.);
    for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

    double rs = P->D->scale_radius();
    y2[2]=rs*tan(y2[2]);
    VecDoub X = {y2[2]*sin(y2[0])*cos(y2[1]),
                 y2[2]*sin(y2[0])*sin(y2[1]),
                 y2[2]*cos(y2[0])};
    VecDoub dP = P->Pot->Forces(X);
    fval[0]=sin(y2[0])*P->D->rho(X)*X[P->los]*(-dP[P->los])*y2[2]*y2[2]*rs*(1.+y2[2]*y2[2]/rs/rs);
    return 0;
}
int sig_los_integrand(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    sig_st *P = (sig_st *) fdata; VecDoub y2(3,0.);
    for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

    double rs = P->D->scale_radius();
    y2[2]=rs*tan(y2[2]);
    VecDoub X = {y2[2]*sin(y2[0])*cos(y2[1]),
                 y2[2]*sin(y2[0])*sin(y2[1]),
                 y2[2]*cos(y2[0])};
    VecDoub dP = P->Pot->Forces(X);
    VecDoub va = P->D->viewing_angles();
    double st= sin(va[0]),ct= cos(va[0]);
    double sp=sin(va[1]),cp=cos(va[1]);
    double f = ct*ct*X[2]*(-dP[2])+st*st*cp*cp*X[0]*(-dP[0])+st*st*sp*sp*X[1]*(-dP[1]);
    fval[0]=sin(y2[0])*P->D->rho(X)*f*y2[2]*y2[2]*rs*(1.+y2[2]*y2[2]/rs/rs);
    return 0;
}

double ObservedTriaxialDensityProfile::sigma_los(Potential_JS *Pot){
  VecDoub x2min = {0.,0.,0.};
  VecDoub x2max = {PI,2.*PI,.5*PI};
  sig_st P(x2min,x2max,this,Pot,0);
  double err;
  double xx = integrate(&sig_los_integrand,&P,1e-3,0,INTEG,&err);
  return sqrt(xx/DP->mass());
}

double ObservedTriaxialDensityProfile::sigma_corr(Potential_JS *Pot){
  VecDoub x2min = {0.,0.,0.};
  VecDoub x2max = {PI,2.*PI,.5*PI};
  sig_st P(x2min,x2max,this,Pot,0);
  double err;
  double xx = integrate(&sig_integrand,&P,1e-4,0,INTEG,&err);
  P = sig_st(x2min,x2max,this,Pot,1);
  double yy = integrate(&sig_integrand,&P,1e-4,0,INTEG,&err);
  P = sig_st(x2min,x2max,this,Pot,2);
  double zz = integrate(&sig_integrand,&P,1e-4,0,INTEG,&err);
  double f1 = xx/zz, f2 = yy/zz;
  double st= sin(theta),ct= cos(theta),sp=sin(phi),cp=cos(phi);
  double f = ct*ct+st*st*cp*cp*f1+st*st*sp*sp*f2;
  return 3.*f/(1.+f1+f2);
}

double ObservedTriaxialDensityProfile::sigma_x_sigma_z(Potential_JS *Pot){
  VecDoub x2min = {0.,0.,0.};
  VecDoub x2max = {PI,2.*PI,.5*PI};
  sig_st P(x2min,x2max,this,Pot,0);
  double err;
  double xx = integrate(&sig_integrand,&P,1e-4,0,INTEG,&err);
  P = sig_st(x2min,x2max,this,Pot,2);
  double zz = integrate(&sig_integrand,&P,1e-4,0,INTEG,&err);
  return xx/zz;
}
//=============================================================================
