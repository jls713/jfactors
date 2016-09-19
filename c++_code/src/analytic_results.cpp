//=============================================================================
#include "analytic_results.h"
#include "observed_density.h"
#include "Multipole.h"
#include "two_component_model.h"
#include <stdexcept>
/**
 * @brief Hernquist X function
 * @details Computes X function from equations (33) & (34) of Hernquist (1990)
 *
 * @param s argument
 * @return Hernquist X
 */
double HernquistX(double s){
    if(s<0.)
        throw std::invalid_argument( "HernquistX received negative value" );
    else if(s<1.)
        return log((1.+sqrt(1-s*s))/s)/sqrt(1-s*s);
    else if(s==1.)
        return 1.;
    else
        return acos(1./s)/sqrt(s*s-1);
}
//=============================================================================
// Spherical NFW profile
//=============================================================================
NFWDensityProfile::NFWDensityProfile(double rho0, double rs):rho0(rho0),rs(rs){}
double NFWDensityProfile::M_sphere(double D, double ang){
	ang = ang/180.*PI;
	double r = D*tan(ang);
	return 4.*PI*rho0*pow(rs,3.)*(log((rs+r)/rs)-r/(r+rs));
}
double NFWDensityProfile::J_factor(double D, double ang){
	ang = ang/180.*PI;
	double Da = D*ang, Da2 = Da*Da, Da3 = Da*Da2, Da4 = Da2*Da2;
	double Delta2 = rs*rs-Da*Da, Delta4 = Delta2*Delta2;
	double rs2 = rs*rs, rs3 = rs2*rs, rs4 = rs2*rs2;
	double X = Da/rs;
	double J = 2.*Da*(7.*Da*rs3-4.*rs*Da3+3.*PI*Delta4)+6./rs*(2.*Delta4*Delta2-2.*rs4*Delta2-rs2*Da4)*HernquistX(X);
	J *= PI*rho0*rho0*rs2/(3.*D*D*Delta4);
	return J;
}
double NFWDensityProfile::D_factor(double D, double ang){
	ang = ang/180.*PI;
	double Da = D*ang, Da2 = Da*Da, Da3 = Da*Da2, Da4 = Da2*Da2;
	double Delta2 = rs*rs-Da*Da, Delta4 = Delta2*Delta2;
	double rs2 = rs*rs, rs3 = rs2*rs, rs4 = rs2*rs2;
	double X = Da/rs;
	double DD = log(X/2.)+HernquistX(X);
	DD *= 4.*PI*rho0*rs3/D/D;
	return DD;
}
//=============================================================================
// Evans flattened cored model
//=============================================================================
CoredDMProfile::CoredDMProfile(double rs, double q, double v0):FiniteMassTriaxialDensityProfile({1.,1.,q},false),r_dm(rs),v0(v0),q(q){
	q_phi_dm = .5*sqrt(1.+sqrt(1.+8.*q*q));
}
double CoredDMProfile::phi_dm(double R, double z){
	return 0.5*v0*v0*log(r_dm*r_dm+R*R+z*z/q_phi_dm/q_phi_dm);
}

VecDoub CoredDMProfile::forces_dm(VecDoub x){
	double R = sqrt(x[0]*x[0]+x[1]*x[1]);
	double cnst = -v0*v0/(r_dm*r_dm+R*R+x[2]*x[2]/q_phi_dm/q_phi_dm);
	return {cnst*x[0],cnst*x[1],cnst*x[2]/q_phi_dm/q_phi_dm};
}

double CoredDMProfile::density_dm(double R, double z){
	double cnst = v0*v0/4./PI/Grav/q_phi_dm/q_phi_dm;
	double rr = r_dm*r_dm+R*R+z*z/q_phi_dm/q_phi_dm;
	return cnst*((2.*q_phi_dm*q_phi_dm+1.)*r_dm*r_dm+R*R+z*z*(2.-1./q_phi_dm/q_phi_dm))/rr/rr;
}

double CoredDMProfile::Phi(const VecDoub &X){
	double R = sqrt(X[0]*X[0]+X[1]*X[1]);
	double z = X[2];
	return phi_dm(R,z);
}
VecDoub CoredDMProfile::Forces(const VecDoub &X){
	return forces_dm(X);
}
double CoredDMProfile::rho(VecDoub X){
	double R = sqrt(X[0]*X[0]+X[1]*X[1]);
	double z = X[2];
	return density_dm(R,z);
}

void CoredDMProfile::scale(double rr, double vv){
	r_dm*=rr; v0*=vv;
}

CoredModel::CoredModel(double n_star, VecDoub radii, VecDoub flattening, VecDoub norms)
	:dm(radii[1],flattening[1],norms[1]),
         staaars({2.,n_star,0.},norms[0],radii[0],-1,{1.,1.,flattening[0]},false){}


// double CoredModel::density_star(double R, double z){
// 	return rho0*pow(r_star,n_star)/pow(r_star*r_star+R*R+z*z/q_star/q_star,n_star/2.);
// }

//double CoredModel::average_los_dispersion(void){
//	return v0/sqrt(n_star);
//}

double CoredModel::half_light_radius(void){
	VecDoub abg = staaars.alpha_beta_gamma();
	return staaars.scale_radius()*sqrt(pow(2.,2./(abg[1]-3.))-1.);
}
double CoredModel::mass_dm(double r){
	return dm.M_sphere(1.,r*180./PI);
}
// int mass_integrand_cored(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
//     sig_cored_st *P = (sig_cored_st *) fdata;
//     VecDoub y2(3,0.);
//     for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

//     double rs = P->C->scale_radius();
//     y2[2]=rs*tan(y2[2]);
//     VecDoub X = {y2[2]*sin(y2[0])*cos(y2[1]),
//                  y2[2]*sin(y2[0])*sin(y2[1]),
//                  y2[2]*cos(y2[0])};
//     double R = sqrt(X[0]*X[0]+X[1]*X[1]);
//     double z = X[2];

//     fval[0]=sin(y2[0])*P->C->density_star(R,z)*y2[2]*y2[2]*rs*(1.+y2[2]*y2[2]/rs/rs);
//     return 0;
// }

// double CoredModel::mass_star(void){
//   VecDoub x2min = {0.,0.,0.};
//   VecDoub x2max = {PI,2.*PI,.5*PI};
//   sig_cored_st P(x2min,x2max,this,0.,0.);
//   double err;
//   return integrate(&mass_integrand_cored,&P,1e-4,0,INTEG,&err);
// }

// int sig_los_integrand_cored(const int ndim[],const double y[], const int*fdim, double fval[], void *fdata){
//     sig_cored_st *P = (sig_cored_st *) fdata;
//     VecDoub y2(3,0.);
//     for(int i=0;i<3;i++) y2[i]=(P->x2max[i]-P->x2min[i])*y[i]+P->x2min[i];

//     double rs = P->C->scale_radius();
//     y2[2]=rs*tan(y2[2]);
//     VecDoub X = {y2[2]*sin(y2[0])*cos(y2[1]),
//                  y2[2]*sin(y2[0])*sin(y2[1]),
//                  y2[2]*cos(y2[0])};
//     double R = sqrt(X[0]*X[0]+X[1]*X[1]);
//     double z = X[2];
//     VecDoub dP = P->C->forces_dm(X);
//     double st= sin(P->theta),ct= cos(P->theta);
//     double sp=sin(P->phi),cp=cos(P->phi);
//     double f = ct*ct*X[2]*(-dP[2])+st*st*cp*cp*X[0]*(-dP[0])+st*st*sp*sp*X[1]*(-dP[1]);
//     fval[0]=sin(y2[0])*P->C->density_star(R,z)*f*y2[2]*y2[2]*rs*(1.+y2[2]*y2[2]/rs/rs);
//     return 0;
// }

// double CoredModel::sigma_los(double theta, double phi){
//   VecDoub x2min = {0.,0.,0.};
//   VecDoub x2max = {PI,2.*PI,.5*PI};
//   sig_cored_st P(x2min,x2max,this,theta,phi);
//   double err;
//   double xx = integrate(&sig_los_integrand_cored,&P,1e-4,0,INTEG,&err);
//   return sqrt(xx/staaars.mass());
// }

// double CoredModel::kinematic_ratio(double theta, double phi){
//   VecDoub x2min = {0.,0.,0.};
//   VecDoub x2max = {PI,2.*PI,.5*PI};
//   sig_cored_st P(x2min,x2max,this,.5*PI,0.);
//   double err;
//   double xx = integrate(&sig_los_integrand_cored,&P,1e-4,0,INTEG,&err);
//   P = sig_cored_st(x2min,x2max,this,.5*PI,.5*PI);
//   double yy = integrate(&sig_los_integrand_cored,&P,1e-4,0,INTEG,&err);
//   P = sig_cored_st(x2min,x2max,this,0.,0.);
//   double zz = integrate(&sig_los_integrand_cored,&P,1e-4,0,INTEG,&err);
//   double f1 = xx/zz, f2 = yy/zz;
//   double st= sin(theta),ct= cos(theta),sp=sin(phi),cp=cos(phi);
//   double f = ct*ct+st*st*cp*cp*f1+st*st*sp*sp*f2;
//   return 3.*f/(1.+f1+f2);
// }

double CoredModel::J_factor_face(double D, double ang){
	ang = ang*PI/180.;
        double q_phi_dm = dm.dm_flattening();
	double q2 = q_phi_dm*q_phi_dm, q4 = q2*q2;
        double r_dm = dm.scale_radius();
        double v0 = dm.amplitude();
	double y = r_dm/sqrt(r_dm*r_dm+D*ang*D*ang);
	double y3 = y*y*y, y5 = y*y*y3;
	return pow(v0,4.)/(96.*r_dm*D*D*Grav*Grav*q2*q_phi_dm)*(3.*(1.-y)-4.*q2*(y3-1.)+q4*(8.-3.*y-2.*y3-3.*y5));
}

double CoredModel::J_factor_face_asymptote(double D){
        double q_phi_dm = dm.dm_flattening();
	double q2 = q_phi_dm*q_phi_dm, q4 = q2*q2;
        double r_dm = dm.scale_radius();
        double v0 = dm.amplitude();
	return pow(v0,4.)/(96.*r_dm*D*D*Grav*Grav*q2*q_phi_dm)*(3.+4.*q2+8.*q4);
}

double CoredModel::D_factor_face(double D, double ang){
	ang = ang*PI/180.;
    double q_phi_dm = dm.dm_flattening();
	double q2 = q_phi_dm*q_phi_dm, q4 = q2*q2;
    double r_dm = dm.scale_radius();
    double v0 = dm.amplitude();
    double y = D*ang/r_dm;
    y = y*y/sqrt(1+y*y);
    return v0*v0*r_dm/Grav/D/D*q_phi_dm*y;
}

double cored_J_edgeon(double p, void *P){
	CoredModel *C = (CoredModel *) P;
	double cp = cos(p), sp = sin(p), s2p = sin(2.*p);
	double cp2 = cp*cp, cp4 = cp2*cp2;
	double sp2 = sp*sp, sp4 = sp2*sp2;
	double q = C->dm_flattening(), q2 = q*q, q4 = q2*q2;
	return (5.*(2.+q2)*s2p*s2p+6.*(2.+2.*q2+q4)*cp4+(30.-4./q2-4./q4)*sp4)/pow(cp2+sp2/q2,3.);
}

double cored_J_edgeon2(double p, void *P){
	CoredModel *C = (CoredModel *) P;
	double cp = cos(p), sp = sin(p);
	double cp2 = cp*cp, cp4 = cp2*cp2;
	double sp2 = sp*sp, sp4 = sp2*sp2;
	double s2p = sin(2.*p), c2p = cos(2.*p), c4p = cos(4.*p);
	double q = C->dm_flattening(), q2 = q*q, q4 = q2*q2, q6 = q2*q4, q8=q4*q4;
	return (6.-6.*q2+83.*q4+28.*q6+9.*q8+4.*(-1.+q)*(1.+q)*(2.+9.*q4+3.*q6)*c2p+(-1.+q2)*(-1.+q2)*(2.+2.*q2+3.*q4)*c4p)/4./q6/pow(cp2+sp2/q2,3.);
}

double CoredModel::J_factor_edge(double D, double ang){
	GaussLegendreIntegrator GL(50);
        double q_phi_dm = dm.dm_flattening();
	double q2 = q_phi_dm*q_phi_dm, q4 = q2*q2;
        double r_dm = dm.scale_radius();
        double v0 = dm.amplitude();
	return pow(v0,4.)/(384.*PI*r_dm*D*D*Grav*Grav*q_phi_dm*q_phi_dm)*GL.integrate(&cored_J_edgeon2,0.,2.*PI,this);
}

double CoredModel::J_factor_arbitrary(double D, double ang, double theta, double phi){
	ObservedTriaxialDensityProfile GG(&dm,theta,phi);
	return GG.J_far_arbitrary_orientation(D,ang);
}
double CoredModel::D_factor_arbitrary(double D, double ang, double theta, double phi){
	ObservedTriaxialDensityProfile GG(&dm,theta,phi);
	return GG.D_far_arbitrary_orientation(D,ang);
}

double CoredModel::sigma_los(double theta, double phi){
	ObservedTriaxialDensityProfile ObsStars(&staaars,theta,phi);
	return ObsStars.sigma_los(&dm);
}

double CoredModel::sigma_los_m(double theta, double phi){
	ObservedTriaxialDensityProfile ObsStars(&staaars,theta,phi);
  	MultipoleDensity MP(&dm);
  	MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,dm.scale_radius(),0.001*dm.scale_radius(),100.*dm.scale_radius());
	return ObsStars.sigma_los(&MEA);
}
void CoredModel::scale(double rh, double slos, std::string dir, double veldispradius){

	double theta=0.,phi=0.;
	if(dir=="round"){ theta=0.; phi=0.;}
	if(dir=="edge"){ theta=.5*PI; phi=0.;}

	ObservedTriaxialDensityProfile ObsStars(&staaars,theta,phi);
	double Rh = half_light_radius();
	double VD = ObsStars.sigma_los(&dm,veldispradius);
	// If viewing a prolate model edge-on, the measured major axis half-light
	// radius is the (half-light-radius of the spherical model)*q (where q>1)
	if(staaars.axis_ratios()[1]>1. and dir=="edge")
		{Rh=Rh*staaars.axis_ratios()[1];}
	double RadiusRatio = rh/Rh;
	double VelRatio = slos/VD;

	double MassRatio = pow(VelRatio,2.)*RadiusRatio;
    staaars = AlphaBetaGammaDensityProfile(staaars.alpha_beta_gamma(),
					staaars.central_density()/pow(RadiusRatio,3.)*MassRatio,
					staaars.scale_radius()*RadiusRatio,
					-1,
					staaars.axes(),false);
	dm = CoredDMProfile(RadiusRatio*dm.scale_radius(),
			    dm.dm_density_flattening(),
			    VelRatio*dm.amplitude());
}
double CoredModel::Mass(double theta, double phi, double r, std::string typ){
  if(typ=="sphere")
	return dm.M_sphere(1.,r*180./PI);
  else if(typ=="ellipsoid")
  	return dm.mass_ellipsoid(r);
}
