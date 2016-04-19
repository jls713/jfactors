//=============================================================================
#include "utils_j.h"
#include <cmath>
#include <vector>
#include <string>
#include "../general/utils.h"
#include "density.h"
//=============================================================================
// Integrands for mass, J-factor and D-factor for general density profile
//=============================================================================

int M_integrand(const int ndim[],const double sphpol_scl[], const int*fdim, double fval[], void *fdata){
  double sphpol[3]; VecDoub X(3,0.);
  density_st *P = (density_st *) fdata;
  // scale to integration range
  for(int i=0;i<3;i++)
    sphpol[i]=(P->x2max[i]-P->x2min[i])*sphpol_scl[i]+P->x2min[i];
  X[0] = sphpol[2]*sin(sphpol[0])*cos(sphpol[1]);
  X[1] = sphpol[2]*sin(sphpol[0])*sin(sphpol[1]);
  X[2] = sphpol[2]*cos(sphpol[0]);
  double r = P->DP->rho(X);
  fval[0]=sin(sphpol[0])*r*sphpol[2]*sphpol[2];
  return 0;
}

int J_integrand(const int ndim[],const double cylpol_scl[], const int*fdim, double fval[], void *fdata){
  double cylpol[3]; VecDoub X(3,0.);
  density_st *P = (density_st *) fdata;
  // scale to integration range
  for(int i=0;i<3;i++)
    cylpol[i]=(P->x2max[i]-P->x2min[i])*cylpol_scl[i]+P->x2min[i];
  double rs = P->DP->scale_radius();
  cylpol[0]=rs*tan(cylpol[0]);
  double b=P->D*tan(cylpol[2]);
  X[0] = (P->los=="x"?cylpol[0]:cos(cylpol[1])*b);
  X[1] = (P->los=="y"?cylpol[0]:sin(cylpol[1])*b);
  X[2] = (P->los=="z"?cylpol[0]:(P->los=="x"?cos(cylpol[1])*b:sin(cylpol[1])*b));
  double r = P->DP->rho(X);
  double jac = rs*(1.+cylpol[0]*cylpol[0]/rs/rs);
  fval[0]=sin(cylpol[2])*r*r*jac;
  return 0;
}

int J_farfield_integrand(const int ndim[],const double cylpol_scl[], const int*fdim, double fval[], void *fdata){
  // small angle, distant
  double cylpol[3]; VecDoub X(3,0.);
  density_st *P = (density_st *) fdata;
  // scale to integration range
  for(int i=0;i<3;i++)
    cylpol[i]=(P->x2max[i]-P->x2min[i])*cylpol_scl[i]+P->x2min[i];
  double rs = P->DP->scale_radius();
  cylpol[0]=rs*tan(cylpol[0]);
  double b=cylpol[2];
  X[0] = (P->los=="x"?cylpol[0]:cos(cylpol[1])*b);
  X[1] = (P->los=="y"?cylpol[0]:sin(cylpol[1])*b);
  X[2] = (P->los=="z"?cylpol[0]:(P->los=="x"?cos(cylpol[1])*b:sin(cylpol[1])*b));
  double r = P->DP->rho(X);
  double jac = rs*(1.+cylpol[0]*cylpol[0]/rs/rs);
  fval[0]=b*r*r*jac;
  return 0;
}
int D_integrand(const int ndim[],const double cylpol_scl[], const int*fdim, double fval[], void *fdata){
  double cylpol[3]; VecDoub X(3,0.);
  density_st *P = (density_st *) fdata;
  // scale to integration range
  for(int i=0;i<3;i++)
    cylpol[i]=(P->x2max[i]-P->x2min[i])*cylpol_scl[i]+P->x2min[i];
  double rs = P->DP->scale_radius();
  cylpol[0]=rs*tan(cylpol[0]);
  double b=P->D*tan(cylpol[2]);
  X[0] = (P->los=="x"?cylpol[0]:cos(cylpol[1])*b);
  X[1] = (P->los=="y"?cylpol[0]:sin(cylpol[1])*b);
  X[2] = (P->los=="z"?cylpol[0]:(P->los=="x"?cos(cylpol[1])*b:sin(cylpol[1])*b));
  double r = P->DP->rho(X);
  double jac = rs*(1.+cylpol[0]*cylpol[0]/rs/rs);
  fval[0]=sin(cylpol[2])*r*jac;
  return 0;
}

int D_farfield_integrand(const int ndim[],const double cylpol_scl[], const int*fdim, double fval[], void *fdata){
  double cylpol[3]; VecDoub X(3,0.);
  density_st *P = (density_st *) fdata;
  // scale to integration range
  for(int i=0;i<3;i++)
    cylpol[i]=(P->x2max[i]-P->x2min[i])*cylpol_scl[i]+P->x2min[i];
  double rs = P->DP->scale_radius();
  cylpol[0]=rs*tan(cylpol[0]);
  double b=cylpol[2];
  X[0] = (P->los=="x"?cylpol[0]:cos(cylpol[1])*b);
  X[1] = (P->los=="y"?cylpol[0]:sin(cylpol[1])*b);
  X[2] = (P->los=="z"?cylpol[0]:(P->los=="x"?cos(cylpol[1])*b:sin(cylpol[1])*b));
  double r = P->DP->rho(X);
  double jac = rs*(1.+cylpol[0]*cylpol[0]/rs/rs);
  fval[0]=b*r*jac;
  return 0;
}

//=============================================================================
// Integrals for mass, J-factors and D-factors for general density profile
//=============================================================================

double DensityProfile::M_sphere(double D, double ang){
  // ang in degrees
  ang = ang/180.*PI;
  VecDoub x2min = {0.,0.,   0.};
  VecDoub x2max = {PI,2.*PI,D*ang};
  density_st P(x2min,x2max,this,D,"");
  double err;
  return integrate(&M_integrand,&P,1e-5,0,INTEG,&err);
}
double DensityProfile::J_factor(double D, double ang, std::string los){
  // ang in degrees
  ang = ang/180.*PI;
  VecDoub x2min = {-.5*PI,0.,   0.};
  VecDoub x2max = {.5*PI, 2.*PI,ang};
  density_st P(x2min,x2max,this,D,los);
  double err;
  return integrate(&J_integrand,&P,1e-4,0,INTEG,&err);
}
double DensityProfile::J_far_factor(double D, double ang, std::string los){
  // ang in degrees
  ang = ang/180.*PI;
  VecDoub x2min = {-.5*PI,0.,0.};
  VecDoub x2max = { .5*PI,2.*PI,ang*D};
  density_st P(x2min,x2max,this,D,los);
  double err;
  return integrate(&J_farfield_integrand,&P,5e-4,0,INTEG,&err)/D/D;
}
double DensityProfile::D_factor(double D, double ang, std::string los){
  // ang in degrees
  ang = ang/180.*PI;
  VecDoub x2min = {-.5*PI,0.,0.};
  VecDoub x2max = {.5*PI,2.*PI,ang};
  density_st P(x2min,x2max,this,D,los);
  double err;
  return integrate(&D_integrand,&P,1e-4,0,INTEG,&err);
}

double DensityProfile::D_far_factor(double D, double ang, std::string los){
  // ang in degrees
  ang = ang/180.*PI;
  VecDoub x2min = {-.5*PI,0.,0.};
  VecDoub x2max = {.5*PI,2.*PI,ang*D};
  density_st P(x2min,x2max,this,D,los);
  double err;
  return integrate(&D_farfield_integrand,&P,1e-4,0,INTEG,&err)/D/D;
}

//=============================================================================
// Properties of triaxial density profile
//=============================================================================

TriaxialDensityProfile::TriaxialDensityProfile(VecDoub ABC, bool normalize)
  : a(ABC[0]), b(ABC[1]), c(ABC[2]){
  ba = b/a; ca = c/a;
  if(normalize){
    a = pow(ba*ca,-1./3.);
    b = ba*a; c = ca*a;
  }
  abc=a*b*c; // 1 if normalized
  T = Triaxiality(ba,ca);
}
TriaxialDensityProfile::TriaxialDensityProfile(double a, double b, double c, bool normalize)
  :a(a),b(b),c(c){
  ba = b/a; ca = c/a;
  if(normalize){
    a = pow(ba*ca,-1./3.);
    b = ba*a; c = ca*a;
  }
  abc=a*b*c; // 1 if normalized
  T = Triaxiality(ba,ca);
}

double TriaxialDensityProfile::Triaxiality(double ba, double ca){
  if(ca==1. and ba!=1.)
    std::cerr<<"b/a<c/a\n";
  if(ca==1. and ba==1.) return 0.;
  return (1.-ba*ba)/(1.-ca*ca);
}

double TriaxialDensityProfile::b_over_a(double T, double ca){
  return sqrt(1.-(1.-ca*ca)*T);
}
double TriaxialDensityProfile::c_over_a(double T, double ba){
  return sqrt(1.-(1.-ba*ba)/T);
}

//=============================================================================
// Properties of finite mass triaxial density profile
//=============================================================================

int M_inf_integrand(const int ndim[],const double sphpol_scl[], const int*fdim, double fval[], void *fdata){
  double sphpol[3]; VecDoub X(3,0.);
  mass_density_st *P = (mass_density_st *) fdata;
  // scale to integration range
  for(int i=0;i<3;i++)
    sphpol[i]=(P->x2max[i]-P->x2min[i])*sphpol_scl[i]+P->x2min[i];
  double rs = P->DP->scale_radius();
  sphpol[2]=rs*tan(sphpol[2]);
  VecDoub ar = P->DP->axes();
  X[0] = ar[0]*sphpol[2]*sin(sphpol[0])*cos(sphpol[1]);
  X[1] = ar[1]*sphpol[2]*sin(sphpol[0])*sin(sphpol[1]);
  X[2] = ar[2]*sphpol[2]*cos(sphpol[0]);
  double r = P->DP->rho(X);
  double jac = rs*(1.+sphpol[2]*sphpol[2]/rs/rs);
  fval[0]=sin(sphpol[0])*r*jac*sphpol[2]*sphpol[2];
  return 0;
}


double FiniteMassTriaxialDensityProfile::mass(void){
  VecDoub x2min = {0.,0.,0.};
  VecDoub x2max = {PI,2.*PI,.5*PI};
  mass_density_st P(x2min,x2max,this,0.);
  double err;
  return abc*integrate(&M_inf_integrand,&P,1e-3,0,INTEG,&err);
}

int M_ell_integrand(const int ndim[],const double sphpol_scl[], const int*fdim, double fval[], void *fdata){
  double sphpol[3]; VecDoub X(3,0.);
  mass_density_st *P = (mass_density_st *) fdata;
  // scale to integration range
  for(int i=0;i<3;i++)
    sphpol[i]=(P->x2max[i]-P->x2min[i])*sphpol_scl[i]+P->x2min[i];
  ;
  VecDoub ar = P->DP->axis_ratios();
  X[0] = sphpol[2]*sin(sphpol[0])*cos(sphpol[1]);
  X[1] = ar[0]*sphpol[2]*sin(sphpol[0])*sin(sphpol[1]);
  X[2] = ar[1]*sphpol[2]*cos(sphpol[0]);
  double r = P->DP->rho(X);
  fval[0]=sin(sphpol[0])*r*sphpol[2]*sphpol[2];
  return 0;
}

double FiniteMassTriaxialDensityProfile::mass_ellipsoid(double r){
  VecDoub x2min = {0.,0.,0.};
  VecDoub x2max = {PI,2.*PI,r};
  mass_density_st P(x2min,x2max,this,0.);
  double err;
  return ba*ca*integrate(&M_ell_integrand,&P,1e-3,0,"Cuhre",&err);
}

int sig_tot_integrand(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
    sigt_st *P = (sigt_st *) fdata; VecDoub y2(3,0.);
    for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

    double rs = P->D->scale_radius();
    y2[2]=rs*tan(y2[2]);
    VecDoub X = {y2[2]*sin(y2[0])*cos(y2[1]),
                 y2[2]*sin(y2[0])*sin(y2[1]),
                 y2[2]*cos(y2[0])};
    VecDoub dP = P->Pot->Forces(X);
    double f = X[0]*(-dP[0])+X[1]*(-dP[1])+X[2]*(-dP[2]);
    fval[0]=sin(y2[0])*P->D->rho(X)*f*y2[2]*y2[2]*rs*(1.+y2[2]*y2[2]/rs/rs);
    return 0;
}

double FiniteMassTriaxialDensityProfile::sigma_tot(Potential_JS *Pot){
  VecDoub x2min = {0.,0.,0.};
  VecDoub x2max = {PI,2.*PI,.5*PI};
  sigt_st P(x2min,x2max,this,Pot);
  double err;
  double xx = integrate(&sig_tot_integrand,&P,1e-4,0,INTEG,&err);
  return sqrt(xx/mass());
}

struct shrf_st{
  FiniteMassTriaxialDensityProfile *D;
  double Mtot;
  shrf_st(FiniteMassTriaxialDensityProfile *D, double Mtot)
    :D(D), Mtot(Mtot){}
};
double shrf_fn(double x, void *p){
  shrf_st *P = (shrf_st *) p;
  double M = P->D->mass_ellipsoid(x*x);
  return 0.5*P->Mtot-M;
}

double FiniteMassTriaxialDensityProfile::spherical_half_light_radius(void){
  double Mtot = mass();
  root_find RF(1e-3,100);
  shrf_st HRF(this,Mtot);
  double rs = scale_radius();
  VecDoub ar = axis_ratios();
  double ca = ar[1];
  if(ca==1.)ca*=0.999;
  double lower = sqrt(rs*ca), upper = sqrt(rs/ca);
  while(shrf_fn(lower,&HRF)<0.) lower *= 0.8;
  while(shrf_fn(upper,&HRF)>0.) upper *= 1.2;
  double result = RF.findroot(&shrf_fn,sqrt(rs*ca),sqrt(rs/ca),&HRF);
  return result*result;
}
//=============================================================================
// Properties of (alpha,beta,gamma) triaxial density profile
//=============================================================================

AlphaBetaGammaDensityProfile::AlphaBetaGammaDensityProfile(VecDoub abg, double rho0, double rs, double rt, VecDoub ABC, bool normalize)
    :alpha(abg[0]),beta(abg[1]),gamma(abg[2]),rho0(rho0),rs(rs),rt(rt),FiniteMassTriaxialDensityProfile(ABC,normalize){
      if((gamma>=3. or beta<=3.) and rt<=0.)
        std::cerr<<"Infinite mass if gamma>3 or beta<3 and rt=0\n";
    }

AlphaBetaGammaDensityProfile::AlphaBetaGammaDensityProfile(double alpha, double beta, double gamma, double rho0, double rs, double rt, VecDoub ABC, bool normalize)
    :alpha(alpha),beta(beta),gamma(gamma),rho0(rho0),rs(rs),rt(rt),FiniteMassTriaxialDensityProfile(ABC,normalize){
      if((gamma>=3. or beta<=3.) and rt<=0.)
        std::cerr<<"Infinite mass if gamma>3 or beta<3 and rt=0\n";
    }

void AlphaBetaGammaDensityProfile::print(void){
  std::cout<<"Axis ratios: "<<a<<" "<<b<<" "<<c<<std::endl;
  std::cout<<"(alpha,beta,gamma): ("<<alpha<<","<<beta<<","<<gamma<<")\n";
  std::cout<<"rs="<<rs<<", rt="<<rt<<std::endl;
}

double AlphaBetaGammaDensityProfile::rho(VecDoub x){
  double r = sqrt(x[0]*x[0]/a/a+x[1]*x[1]/b/b+x[2]*x[2]/c/c);
  double fac = 1.;
  if(rt>0.){
      fac = tanh(r/rt);
      fac=sqrt(1.-fac*fac);
  }
  r = r/rs;
  return rho0*pow(r,-gamma)*pow(1.+pow(r,alpha),((gamma-beta)/alpha))*fac;
}

//=============================================================================
