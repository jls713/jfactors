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
// Test programs
//=============================================================================
#include "gtest/gtest.h"

const double test_err = 1e-7;

namespace {

TEST(Profile,SphericalNFW){
	/** Test numerical integration of mass, J-factor and D-factor against
		analytic result for NFW **/
	double rho0 = 2.3, rs = 0.053, rt = 1000.;
	AlphaBetaGammaDensityProfile ABG({1.,3.,1.},rho0,rs,rt,{1.,1.,1.},false);
	NFWDensityProfile NFW(rho0,rs);
	EXPECT_NEAR(NFW.M_sphere(30.,0.5),ABG.M_sphere(30.,.5),test_err);
	EXPECT_NEAR(NFW.J_factor(30.,0.5),ABG.J_factor(30.,.5),test_err);
	EXPECT_NEAR(NFW.J_factor(30.,0.5),ABG.J_far_factor(30.,.5),test_err);
	EXPECT_NEAR(NFW.D_factor(30.,0.5),ABG.D_factor(30.,.5),test_err);
	EXPECT_NEAR(NFW.D_factor(30.,0.5),ABG.D_far_factor(30.,.5),test_err);
}

TEST(Profile,TriaxialUnNorm){
	/** Test unnormalized triaxial total mass, mass in an ellipsoid and
		ellipsoidal half-light radius. Here the mass distribution is stretched
		so the ellipsoidal half-light radius is a r_h  **/
	VecDoub abc = {1.,0.6,0.4};
	double rho0 = 2.3, rs = 0.053, rt = 1000.;
	AlphaBetaGammaDensityProfile ABG({2.,5.,0.},rho0,rs,rt,abc,false);
	EXPECT_NEAR((4*PI)/3.*abc[0]*abc[1]*abc[2]*rho0*rs*rs*rs,ABG.mass(),test_err);
	std::cout<<ABG.mass_ellipsoid(rs)<<" "<<ABG.M_sphere(1.,rs*180./PI)<<std::endl;
	EXPECT_NEAR((4*PI)/3.*abc[0]*abc[1]*abc[2]*rho0*rs*rs*rs,ABG.mass_ellipsoid(100.),test_err);
	EXPECT_NEAR(rs/sqrt(pow(2.,2./3.)-1.),ABG.spherical_half_light_radius(),1e-4);
	EXPECT_NEAR((4*PI)/3.*abc[0]*abc[1]*abc[2]*rho0*rs*rs*rs/2.,ABG.mass_ellipsoid(rs/sqrt(pow(2.,2./3.)-1.)),test_err);
	double T = (1.-abc[1]*abc[1]/abc[0]/abc[0])/(1.-abc[2]*abc[2]/abc[0]/abc[0]);
	EXPECT_DOUBLE_EQ(ABG.Triaxiality(),T);
	EXPECT_DOUBLE_EQ(ABG.b_over_a(),abc[1]/abc[0]);
	EXPECT_DOUBLE_EQ(ABG.c_over_a(),abc[2]/abc[0]);

}

TEST(Profile,TriaxialNorm){
	/** Test unnormalized triaxial total mass, mass in an ellipsoid and
		ellipsoidal half-light radius. **/
	VecDoub abc = {1.,.6,0.4};
	double rho0 = 2.3, rs = 0.053, rt = 1000.;
	AlphaBetaGammaDensityProfile ABG({2.,5.,0.},rho0,rs,rt,abc,true);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs,ABG.mass(),1e-6);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs,ABG.mass_ellipsoid(200.),test_err);
	VecDoub ar = ABG.axes();
	EXPECT_NEAR(rs/sqrt(pow(2.,2./3.)-1.)*ar[0],ABG.spherical_half_light_radius(),1e-4);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs/2.,ABG.mass_ellipsoid(rs/sqrt(pow(2.,2./3.)-1.)*ar[0]),test_err);
	double T = (1.-abc[1]*abc[1]/abc[0]/abc[0])/(1.-abc[2]*abc[2]/abc[0]/abc[0]);
	EXPECT_DOUBLE_EQ(ABG.Triaxiality(),T);
	EXPECT_DOUBLE_EQ(ABG.b_over_a(),abc[1]/abc[0]);
	EXPECT_DOUBLE_EQ(ABG.c_over_a(),abc[2]/abc[0]);
}
TEST(ObservedProfile,SphericalPlummer){
	/** Test total mass, contained mass, spherical half-light radius, mass in
		a column for Plummer sphere **/
	double rho0 = 2.3, rs = 0.053, rt = 1000.;
	AlphaBetaGammaDensityProfile ABG({2.,5.,0.},rho0,rs,rt,{1.,1.,1.},false);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs,ABG.M_sphere(1000.,10.),test_err);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs,ABG.mass(),1e-6);
	EXPECT_NEAR(rs/sqrt(pow(2.,2./3.)-1.),ABG.spherical_half_light_radius(),1e-4);
	ObservedTriaxialDensityProfile OP(&ABG,0.,0.);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs,OP.M_ellipse(100.),1e-6);
	EXPECT_NEAR(rs,OP.half_light_radius(),1e-4);
}
TEST(ObservedProfile,SphericalHalfLightRadii){
	/** Test half-light radius for a series of cored models **/
	double rho0 = 2.3, rs = 0.053, rt = -1.;
	VecDoub N = {3.5,4.,5.,6.,8.};
	for(auto n: N){
		AlphaBetaGammaDensityProfile ABG({2.,n,0.},rho0,rs,rt,{1.,1.,1.},false);
		EXPECT_NEAR(pow(2.,n-3.)*complete_beta(0.5*(n-3.),0.5*(n-3.))*PI*rho0*rs*rs*rs/(n-2.),ABG.mass(),1e-5);
		ObservedTriaxialDensityProfile OP(&ABG,0.,0.);
		EXPECT_NEAR(rs*sqrt(pow(2.,2./(n-3.))-1.),OP.half_light_radius(),5e-4);
	}
}

TEST(ObservedProfile,AxisymmetricHalfLightRadii){
	/** Test half-light radius for a series of cored models **/
	double rho0 = 2.3, rs = 0.053, rt = -1., q = 0.8;
	VecDoub N = {4.,5.,6.,8.};
	for(auto n: N){
		AlphaBetaGammaDensityProfile ABG({2.,n,0.},rho0,rs,rt,{1.,1.,q},false);
		EXPECT_NEAR(pow(2.,n-3.)*complete_beta(0.5*(n-3.),0.5*(n-3.))*PI*rho0*rs*rs*rs*q/(n-2.),ABG.mass(),1e-6);
		ObservedTriaxialDensityProfile OP(&ABG,PI/2.,0.);
		EXPECT_NEAR(rs*sqrt(pow(2.,2./(n-3.))-1.),OP.half_light_radius(),1e-4);
		ObservedTriaxialDensityProfile OP2(&ABG,0.,0.);
		EXPECT_NEAR(rs*sqrt(pow(2.,2./(n-3.))-1.),OP2.half_light_radius(),1e-4);
	}
}

TEST(ObservedProfile,AxisymmetricHalfLightRadiiProlate){
	/** Test half-light radius for a series of cored models **/
	double rho0 = 2.3, rs = 0.053, rt = -1., q = 0.8;
	VecDoub N = {4.,5.,6.,8.};
	for(auto n: N){
		AlphaBetaGammaDensityProfile ABG({2.,n,0.},rho0,rs,rt,{1.,q,q},false);
		EXPECT_NEAR(pow(2.,n-3.)*complete_beta(0.5*(n-3.),0.5*(n-3.))*PI*rho0*rs*rs*rs*q*q/(n-2.),ABG.mass(),1e-6);
		ObservedTriaxialDensityProfile OP(&ABG,0.,0.);
		EXPECT_NEAR(rs*sqrt(pow(2.,2./(n-3.))-1.),OP.half_light_radius(),1e-4);
		ObservedTriaxialDensityProfile OP2(&ABG,PI/2.,0.);
		EXPECT_NEAR(rs*sqrt(pow(2.,2./(n-3.))-1.)*q,OP2.half_light_radius(),1e-4);
	}
}

TEST(ObservedProfile, SphericalNFW){
	/** Test arbitrary orientation J-factors for spherical case **/
	double rho0 = 2.3, rs = 0.053, rt = 1000.;
	double theta = 0.33, phi = 0.78;
	AlphaBetaGammaDensityProfile ABG({1.,3.,1.},rho0,rs,rt,{1.,1.,1.},false);
	ObservedTriaxialDensityProfile OP(&ABG,theta,phi);
	EXPECT_NEAR(OP.J_far_arbitrary_orientation(30.,0.5),ABG.J_factor(30.,.5),test_err);
}
TEST(ObservedProfile,TriaxialPlummer){
	double rho0 = 2.3, rs = 1.053, rt = -1.;
	double theta = 0., phi = 0.;
	double ba=0.4001,ca=0.4;
	AlphaBetaGammaDensityProfile ABG({2.,5.,0.},rho0,rs,rt,{1.,ba,ca},false);
	ObservedTriaxialDensityProfile OP(&ABG,theta,phi);
	EXPECT_NEAR(ba,OP.ellipticity(),1e-5);
	EXPECT_NEAR(PI/2.,OP.minor_axis(),1e-3);
	double rr = OP.half_light_radius();
	EXPECT_NEAR(rs,rr,0.001);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs/2.*ba*ca,OP.M_ellipse(rr),1e-4);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs*ba*ca,ABG.mass(),1e-4);

	theta=PI/2.; phi=0.;
	ObservedTriaxialDensityProfile OP2(&ABG,theta,phi);
	EXPECT_NEAR(ca/ba,OP2.ellipticity(),1e-5);
	EXPECT_NEAR(0.,OP2.minor_axis(),1e-5);
	rr = OP2.half_light_radius();
	EXPECT_NEAR(rs*ba,rr,0.001);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs/2.*ba*ca,OP2.M_ellipse(rr),1e-3);

	theta=PI/2.; phi=PI/2.;
	ObservedTriaxialDensityProfile OP3(&ABG,theta,phi);
	EXPECT_NEAR(ca,OP3.ellipticity(),1e-5);
	EXPECT_NEAR(0.,OP3.minor_axis(),test_err);
	rr = OP3.half_light_radius();
	EXPECT_NEAR(rs,rr,0.001);
	EXPECT_NEAR((4*PI)/3.*rho0*rs*rs*rs/2.*ba*ca,OP3.M_ellipse(rr),1e-3);


	double T = ABG.Triaxiality();
	phi=0.; theta = atan(sqrt(T/(1.-T)));
	ObservedTriaxialDensityProfile OP4(&ABG,theta,phi);
	EXPECT_NEAR(1.,OP4.ellipticity(),test_err);
}

TEST(VelocityDispersions,Spherical){
	double rho0 = 2.3, rs = 1.053, rt = 10.;
	double ba = 1., ca = 1.;
	AlphaBetaGammaDensityProfile ABG({2.,5.,0.},rho0,rs,rt,{1.,ba,ca},false);
	AlphaBetaGammaDensityProfile NFW({1.,3.,1.},rho0,rs,rt,{1.,ba,ca},false);
	MultipoleDensity MP(&NFW);
	MultipoleExpansion_Triaxial MEA(&MP,150,16,12,8,NFW.scale_radius(),0.001*NFW.scale_radius(),10.*NFW.tidal_radius());

	double st = ABG.sigma_tot(&MEA);

	double theta = 0., phi = 0.;
	ObservedTriaxialDensityProfile OP(&ABG,theta,phi);
	double sl = OP.sigma_los(&MEA);

	EXPECT_NEAR(st/sqrt(3.),sl,1e-4);
	EXPECT_NEAR(OP.sigma_los(&MEA,10.*rt),sl,1e-4);
}
double binney_tremaine_virial_ratio(double q){
	/* Table 2.1 of Binney & Tremaine (2008) -- <sigma_xx^2>/<sigma_zz^2>
		for stellar density stratified on spheroidal shells and dark matter
		density stratified on spheroidal shells of the same ellipticity */
	if(q<1){
		double e = sqrt(1.-q*q);
		return .5*(asin(e)/e-q)/(1./q-asin(e)/e)/q/q;
	}
	else{
		double e = sqrt(1.-1./q/q);
		return .5*(q*q-.5*log((1.+e)/(1.-e))/e)/(.5*log((1.+e)/(1.-e))/e-1.)/q/q;
	}
}

TEST(Jfactors,DoubleModel){
	double q = 0.8, rh = 0.049, slos = 3.22;
	double rsrdmratio = 2., rtdmrdmratio = 10., rtsrdmratio = 10.;
	VecDoub abg_st = {2.,5.,0.}, abg_dm = {1.,3.,1.};
	bool with_multipole = true;
	DoubleProfileModel DD(1.,q,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
	EXPECT_NEAR(DD.sigma_x_sigma_z(),binney_tremaine_virial_ratio(q),1e-3);
	double theta = PI/2., phi = 0., D = 30., ang = 0.5;
	double JJ = DD.J_factor(theta,phi,D,ang,false,false)[0];
	double JJ2 = DD.J_factor_old(theta,phi,D,ang,false,false)[0];

	DoubleProfileModel DD_sph(1.,1.,sqrt(q)*rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
	double JJs = DD_sph.J_factor(theta,phi,D,ang,false,false)[0];
	double JJ2s = DD_sph.J_factor_old(theta,phi,D,ang,false,false)[0];
	EXPECT_NEAR(JJ/JJs,JJ2/JJ2s,0.05);

	double CC = DD.correction_factor(theta,phi,D,ang,false,true);
	EXPECT_NEAR(JJ/JJs,CC,0.001);

}


TEST(Jfactors,DoubleModelProlate){
	double q = .25, rh = 0.049, slos = 3.22;
	double rsrdmratio = 2., rtdmrdmratio = 10., rtsrdmratio = 10.;
	VecDoub abg_st = {2.,5.,0.}, abg_dm = {1.,3.,1.};
	bool with_multipole = true;
	DoubleProfileModel DD(q,q,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
	EXPECT_NEAR(DD.sigma_x_sigma_z(),1./binney_tremaine_virial_ratio(1./q),3e-3);
	double theta = 0., phi = 0., D = 30., ang = 0.5;
	double JJ = DD.J_factor(theta,phi,D,ang,false,false)[0];

	DoubleProfileModel DD_sph(1.,1.,sqrt(q)*rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
	double JJs = DD_sph.J_factor(theta,phi,D,ang,false,false)[0];

	double CC = DD.correction_factor(theta,phi,D,ang,false,true);
	EXPECT_NEAR(JJ/JJs,CC,0.005);

}

double BigF(double Rc, double Rd, double qp, double q){
	double D1 = sqrt(qp*qp-q*q);
	double D2 = sqrt(Rd*Rd-Rc*Rc);
	return (atan2(Rd*D1,q*D2)-atan2(Rc*D1,qp*D2))/D1/D2;
}

double sigma_R_wyn4(double Rc, double Rd, double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	double D12 = qp*qp-q*q;
	double D22 = Rd*Rd-Rc*Rc;
	double F = BigF(Rc,Rd,qp,q);
	return sqrt(qp*Rc/D12/D22*((qp*qp*Rd*Rd+q*q*(Rc*Rc-2.*Rd*Rd))*F-qp*Rc+q*Rd));
}

double sigma_z_wyn4(double Rc, double Rd, double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	double D12 = qp*qp-q*q;
	double D22 = Rd*Rd-Rc*Rc;
	double F = BigF(Rc,Rd,qp,q);
	return sqrt(q*q*Rc/D12*(qp*F-q/(q*Rc+qp*Rd)));
}

double sigma_R_wyn4_RcRd(double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	return sqrt(2.*qp*(2.*q+qp)/3./pow(q+qp,2.));
}
double sigma_z_wyn4_RcRd(double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	return sqrt(q*q/(q+qp)/(q+qp));
}

double sigma_R_wyn5(double Rc, double Rd, double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	double D2 = qp*qp*Rd*Rd-q*q*Rc*Rc;
	double lR = Rc*Rc-Rd*Rd;
	double lq = qp*qp-q*q;
	return sqrt(Rc*Rc/lR*(
	            2.*qp*Rd*Rd*Rd/lR/sqrt(D2)*acosh(qp*Rd/q/Rc)
	           -(qp*(2.*pow(qp*Rd,2.)+q*q*(Rc*Rc-3.*Rd*Rd))*acosh(qp/q))/lR/pow(lq,1.5)
	           +qp*qp/lq));
}

double sigma_z_wyn5(double Rc, double Rd, double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	double D2 = qp*qp*Rd*Rd-q*q*Rc*Rc;
	double lR = Rc*Rc-Rd*Rd;
	double lq = qp*qp-q*q;
	return sqrt(q*q*Rc*Rc/lR*(q*q*lR/lq/D2+qp*acosh(qp/q)/pow(lq,1.5)-qp*Rd*Rd*Rd/pow(D2,1.5)*acosh(qp*Rd/q/Rc)));
}

double sigma_R_wyn5_RcRd(double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	return sqrt(0.25*qp/pow(qp*qp-q*q,3.)*(5.*q*q*q*q*qp-7.*q*q*qp*qp*qp+2.*qp*qp*qp*qp*qp+3.*q*q*q*q*sqrt(qp*qp-q*q)*acosh(qp/q)));
}
double sigma_z_wyn5_RcRd(double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	return sqrt(.5*q*q/pow(qp*qp-q*q,3.)*(qp*qp*qp*qp+qp*qp*q*q-2.*q*q*q*q-3.*q*q*qp*sqrt(qp*qp-q*q)*acosh(qp/q)));
}

double sigma_R_wyn6(double Rc, double Rd, double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	double D12 = qp*qp-q*q;
	double D22 = Rd*Rd-Rc*Rc;
	double F = BigF(Rc,Rd,qp,q);
	return sqrt(qp*Rc*Rc/D12/D22/D22*((Rd*(Rc*Rc+2.*Rd*Rd)*D12-q*qp*Rc*D22)/(q*Rc+qp*Rd)-Rc*(3.*qp*qp*Rd*Rd+q*q*(Rc*Rc-4.*Rd*Rd))*F));
}

double sigma_z_wyn6(double Rc, double Rd, double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	double D12 = qp*qp-q*q;
	double D22 = Rd*Rd-Rc*Rc;
	double F = BigF(Rc,Rd,qp,q);
	return sqrt(q*q*Rc*Rc/D12/D22*(-qp*Rc*F+(q*qp*Rc*Rd+qp*qp*Rd*Rd-q*q*D22)/(pow(q*Rc+qp*Rd,2.))));
}

double sigma_R_wyn6_RcRd(double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	return sqrt(2.*qp*(8.*q*q+9.*q*qp+3.*qp*qp)/15./pow(q+qp,3.));
}
double sigma_z_wyn6_RcRd(double q, double qdm){
	double qp = .5*sqrt(1.+sqrt(1.+8.*qdm*qdm));
	return sqrt(q*q*(3.*q+qp)/3./pow(q+qp,3.));
}

TEST(AnalyticProfile,CoredModel5){

	double D = 3000., ang = 0.5, nstar=5.;
	double q = 0.3;
	CoredModel Csph(nstar,{1.,1.},{1.,1.},{1.,1.});
	EXPECT_NEAR(Csph.sigma_los(0.,0.),1./sqrt(nstar),1e-4);

	CoredModel C(nstar,{1.,2.},{1.,1.},{1.,1.});
	C.scale(0.05,3.);
	EXPECT_NEAR(C.J_factor_edge(D,ang),C.J_factor_face_asymptote(D),1e+5);
	EXPECT_NEAR(C.J_factor_arbitrary(D,ang,0.,0.),C.J_factor_face(D,ang),1e+5);
	EXPECT_NEAR(C.D_factor_arbitrary(D,ang,0.,0.),C.D_factor_face(D,ang),1e+5);

	CoredModel C2(nstar,{1.,2.},{q,q},{1.,1.});
	EXPECT_NEAR(sigma_z_wyn5(1.,1.00001,0.999,0.999),1./sqrt(5.),1e-3);
	EXPECT_NEAR(sigma_R_wyn5(1.,1.00001,0.999,0.999),sqrt(2./5.),1e-3);
	EXPECT_NEAR(C2.sigma_los(0.,0.),sigma_z_wyn5(1.,2.,q,q),1e-4);
	EXPECT_NEAR(C2.sigma_los(0.,0.),C2.sigma_los_m(0.,0.),1e-3);
	EXPECT_NEAR(C2.sigma_los(.5*PI,0.),sigma_R_wyn5(1.,2.,q,q)/sqrt(2.),1e-4);
	EXPECT_NEAR(C2.sigma_los(.5*PI,0.),C2.sigma_los_m(.5*PI,0.),1e-3);
	EXPECT_NEAR(C2.sigma_los(.5*PI,0.),C2.sigma_los(.5*PI,0.,100.),1e-3);
	C2.scale(0.05,3.);
	EXPECT_NEAR(C2.dm_density_flattening(),q,0.0001);
	EXPECT_NEAR(C2.J_factor_arbitrary(D,ang,0.,0.),C2.J_factor_face(D,ang),1e+6);
	EXPECT_NEAR(C2.J_factor_arbitrary(D,ang,PI/2.,0.),C2.J_factor_edge(D,ang),1e+7);

	CoredModel C3(nstar,{1.,2.},{1./q,1./q},{1.,1.});
	C3.scale(0.05,3.);
	EXPECT_NEAR(C3.J_factor_arbitrary(D,ang,0.,0.),C3.J_factor_face(D,ang),1e+4);
	EXPECT_NEAR(C3.J_factor_arbitrary(D,ang,PI/2.,0.),C3.J_factor_edge(D,ang),1e+5);
	CoredModel C4(nstar,{1.,1.},{q,q},{1.,1.});
	EXPECT_NEAR(C4.sigma_los(0.,0.),sigma_z_wyn5_RcRd(q,q),1e-4);
	EXPECT_NEAR(C4.sigma_los(.5*PI,0.),sigma_R_wyn5_RcRd(q,q)/sqrt(2.),1e-4);
	EXPECT_NEAR(C4.sigma_los(.5*PI,0.),C4.sigma_los(.5*PI,0.,100.),1e-3);
}

TEST(AnalyticProfile,CoredModel4){

	double D = 3000., ang = 0.5, nstar=4.;
	double q = 0.2;
	CoredModel Csph(nstar,{1.,1.},{1.,1.},{1.,1.});
	EXPECT_NEAR(Csph.sigma_los(0.,0.),1./sqrt(nstar),1e-4);

	CoredModel C(nstar,{1.,2.},{1.,1.},{1.,1.});
	C.scale(0.05,3.);
	EXPECT_NEAR(C.J_factor_edge(D,ang),C.J_factor_face_asymptote(D),1e+5);
	EXPECT_NEAR(C.J_factor_arbitrary(D,ang,0.,0.),C.J_factor_face(D,ang),1e+5);
	EXPECT_NEAR(C.D_factor_arbitrary(D,ang,0.,0.),C.D_factor_face(D,ang),1e+5);

	CoredModel C2(nstar,{1.,2.},{q,q},{1.,1.});
	EXPECT_NEAR(sigma_z_wyn4(1.,1.00001,0.999,0.999),1./sqrt(4.),1e-3);
	EXPECT_NEAR(C2.sigma_los(0.,0.),sigma_z_wyn4(1.,2.,q,q),1e-4);
	EXPECT_NEAR(C2.sigma_los(0.,0.),C2.sigma_los_m(0.,0.),1e-3);
	EXPECT_NEAR(C2.sigma_los(.5*PI,0.),sigma_R_wyn4(1.,2.,q,q)/sqrt(2.),1e-4);
	EXPECT_NEAR(C2.sigma_los(.5*PI,0.),C2.sigma_los_m(.5*PI,0.),1e-3);
	C2.scale(0.05,3.);
	EXPECT_NEAR(C2.J_factor_arbitrary(D,ang,0.,0.),C2.J_factor_face(D,ang),1e+7);
	EXPECT_NEAR(C2.J_factor_arbitrary(D,ang,PI/2.,0.),C2.J_factor_edge(D,ang),1e+7);

	CoredModel C3(nstar,{1.,2.},{1./q,1./q},{1.,1.});
	C3.scale(0.05,3.);
	EXPECT_NEAR(C3.J_factor_arbitrary(D,ang,0.,0.),C3.J_factor_face(D,ang),1e+7);
	EXPECT_NEAR(C3.J_factor_arbitrary(D,ang,PI/2.,0.),C3.J_factor_edge(D,ang),1e+7);
	CoredModel C4(nstar,{1.,1.},{q,q},{1.,1.});
	EXPECT_NEAR(C4.sigma_los(0.,0.),sigma_z_wyn4_RcRd(q,q),1e-4);
	EXPECT_NEAR(C4.sigma_los(.5*PI,0.),sigma_R_wyn4_RcRd(q,q)/sqrt(2.),1e-4);
	EXPECT_NEAR(C4.sigma_los(.5*PI,0.),C4.sigma_los(.5*PI,0.,100.),2e-3);
}

TEST(AnalyticProfile,CoredModel6){

	double D = 3000., ang = 0.5, nstar=6.;
	double q = 0.2;
	CoredModel Csph(nstar,{1.,1.},{1.,1.},{1.,1.});
	EXPECT_NEAR(Csph.sigma_los(0.,0.),1./sqrt(nstar),1e-4);

	CoredModel C(nstar,{1.,2.},{1.,1.},{1.,1.});
	auto rg = 2.;
	EXPECT_NEAR(C.mass_dm(rg),rg*rg*rg/Grav/(rg*rg+4.),1.);
	C.scale(0.05,3.);
	EXPECT_NEAR(C.J_factor_edge(D,ang),C.J_factor_face_asymptote(D),1e+5);
	EXPECT_NEAR(C.J_factor_arbitrary(D,ang,0.,0.),C.J_factor_face(D,ang),1e+5);
	EXPECT_NEAR(C.D_factor_arbitrary(D,ang,0.,0.),C.D_factor_face(D,ang),1e+5);

	CoredModel C2(nstar,{1.,2.},{q,q},{1.,1.});
	EXPECT_NEAR(sigma_z_wyn6(1.,1.00001,0.999,0.999),1./sqrt(6.),1e-3);
	EXPECT_NEAR(C2.sigma_los(0.,0.),sigma_z_wyn6(1.,2.,q,q),1e-4);
	EXPECT_NEAR(C2.sigma_los(0.,0.),C2.sigma_los_m(0.,0.),1e-3);
	EXPECT_NEAR(C2.sigma_los(.5*PI,0.),sigma_R_wyn6(1.,2.,q,q)/sqrt(2.),1e-4);
	EXPECT_NEAR(C2.sigma_los(.5*PI,0.),C2.sigma_los_m(.5*PI,0.),1e-3);
	C2.scale(0.05,3.);
	EXPECT_NEAR(C2.J_factor_arbitrary(D,ang,0.,0.),C2.J_factor_face(D,ang),1e+7);
	EXPECT_NEAR(C2.J_factor_arbitrary(D,ang,PI/2.,0.),C2.J_factor_edge(D,ang),1e+9);

	CoredModel C3(nstar,{1.,2.},{1./q,1./q},{1.,1.});
	C3.scale(0.05,3.);
	EXPECT_NEAR(C3.J_factor_arbitrary(D,ang,0.,0.),C3.J_factor_face(D,ang),1e+7);
	EXPECT_NEAR(C3.J_factor_arbitrary(D,ang,PI/2.,0.),C3.J_factor_edge(D,ang),1e+7);
	CoredModel C4(nstar,{1.,1.},{q,q},{1.,1.});
	EXPECT_NEAR(C4.sigma_los(0.,0.),sigma_z_wyn6_RcRd(q,q),1e-4);
	EXPECT_NEAR(C4.sigma_los(.5*PI,0.),sigma_R_wyn6_RcRd(q,q)/sqrt(2.),1e-4);
	EXPECT_NEAR(C4.sigma_los(.5*PI,0.),C4.sigma_los(.5*PI,0.,100.),1e-3);
}

TEST(Jfactors,TestModel){
	double q = 1., rh = 0.049, slos = 3.22;
	double rsrdmratio = 2., rtdmrdmratio = 10., rtsrdmratio = 10.;
	VecDoub abg_st = {2.,5.,0.}, abg_dm = {1.,3.,1.};
	bool with_multipole = true;
	DoubleProfileModel DD(1.,q,rh,slos,with_multipole,rsrdmratio,rtdmrdmratio,rtsrdmratio,abg_st,abg_dm);
	double theta = PI/2., phi = 0., D = 30., ang = 0.5;
	double JJ = DD.J_factor(theta,phi,D,ang,false,false)[0];
}

}
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
