// ============================================================================
#ifndef OBSDENS_H
#define OBSDENS_H
#include <Python.h>
#include "utils_j.h"
#include <cmath>
#include <vector>
#include <string>
#include "../general/utils.h"
#include "density.h"
#include "observed_density.h"
// ============================================================================

class ObservedDensityProfile{
  /**
   * @brief Observed density profile
   * -- density profile observed at arbitrary angles
   */
protected:
  DensityProfile *DP;           /* Base density profile pointer              */
  double theta, phi;            /* observation angles -- spherical polars    */
  double cp, sp, ct, st;        /* cos(phi), sin(phi), cos(theta), sin(theta)*/
  VecDoub phihat, thetahat;     /* phi unit vector, theta unit vector        */
  VecDoub n;                    /* radial unit vector */
public:
  /**
   * @brief Observed density profile
   * @details base class for observed density profiles -- takes density profile
   * and 'observes' it from a viewing angle (theta,phi) which are spherical
   * polar coordinates
   *
   * @param DP pointer to base density profile
   * @param theta spherical polar latitude angle -- viewing angle
   * @param phi spherical polar azimuthal angle  -- viewing angle
   */
  ObservedDensityProfile(DensityProfile *DP, double theta, double phi);
  /**
   * @brief J-factor for density profile
   * @details J-factor for density profile observed at (theta,phi) in the
   * limit of large distance from dwarf so beam = cylinder
   *
   * @param D distance to dwarf
   * @param ang beam angle in degrees
   *
   * @return J-factor
   */
  double J_far_arbitrary_orientation(double D, double ang);
  inline double rho(VecDoub x){return DP->rho(x);}
  inline double scale_radius(void){return DP->scale_radius();}
  inline VecDoub viewing_angles(void){return {theta,phi};}
  /**
   * @brief intrinsic coordinates from viewing coordinates
   * @details computes the intrinsic Cartesian (x,y,z) coordinates from a
   * a vector of observed coordinates (on sky location, line-of-sight distance)
   *
   * @param X Cartesian position in on-sky coordinates x,y,z
   * @return intrinsic coordinates
   */
  VecDoub x_proj(VecDoub X);
};
// ============================================================================

class ObservedTriaxialDensityProfile: public ObservedDensityProfile{
/**
 * @brief Observed triaxial density profile
 * -- triaxial density profile observed at arbitrary angles
*/
protected:
  FiniteMassTriaxialDensityProfile *DP; /*Base density profile pointer*/
  double e, theta_min;   /* observed ellipticity and minor axis on-sky angle */
  double stm,ctm;        /* sin(theta_min), cos(theta_min)                   */
  /**
   * @brief observed ellipticity
   * @details observed ellipticity from intrinsic ellipticity and viewing
   *  angles (equations A1 and A2 of Weijmans 2014)
   *
   * @param ba intermediate-to-major axis ratio
   * @param ca minor-to-major axis ratio
   *
   * @return on-sky ellipticity
   */
  double observed_ellipticity(double ba, double ca);
    /**
   * @brief minor axis orientation
   * @details on-sky minor axis orientation from intrinsic ellipticity
   * and viewing angles (equations A6 of Weijmans 2014)
   *
   * @param ba intermediate-to-major axis ratio
   * @param ca minor-to-major axis ratio
   *
   * @return on-sky minor axis orientation
   */
  double minor_axis_angle(double ba, double ca);
public:
  /**
   * @brief Observed triaxial density profile
   * @details observes a triaxial density profile from arbitrary viewing
   * angles (theta,phi)=spherical polar
   *
   * @param DP FiniteMassTriaxialDensityProfile pointer
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   */
  ObservedTriaxialDensityProfile(FiniteMassTriaxialDensityProfile *DP, double theta, double phi);
  VecDoub rot_angles(void){return {ctm,stm};}
  double ellipticity(void){return e;}
  double minor_axis(void){return theta_min;}
  /**
   * @brief J-factor for density profile
   * @details J-factor for density profile observed at (theta,phi) in the
   * limit of large distance from dwarf so beam = cylinder
   *
   * @param D distance to dwarf
   * @param ang beam angle in degrees
   *
   * @return J-factor
   */
  double J_far_arbitrary_orientation(double D, double ang);
  /**
   * @brief D-factor for density profile
   * @details D-factor for density profile observed at (theta,phi) in the
   * limit of large distance from dwarf so beam = cylinder
   *
   * @param D distance to dwarf
   * @param ang beam angle in degrees
   *
   * @return D-factor
   */
  double D_far_arbitrary_orientation(double D, double ang);
  /**
   * @brief mass inside elliptical cylinder
   * @details mass inside a cylinder with elliptical cross-section = observed
   * ellipticity of density profile
   *
   * @param major axis length of elliptical cross-section
   *
   * @return mass inside elliptical cylinder
   */
  double M_ellipse(double x);
  /**
   * @brief half-light elliptical radius
   * @details major axis length of cylinder that contains half the mass
   * @return half-light elliptical radius
   */
  double half_light_radius(void);
  /**
   * @brief l.o.s. velocity dispersion
   * @details velocity dispersion along the line-of-sight
   *
   * @param Pot pointer to a Potential_JS instance
   * @param radius within which to compute dispersion (if <0 over whole galaxy)
   * @return l.o.s. velocity dispersion
   */
  double sigma_los(Potential_JS *Pot, double radius=-1.);
  /**
   * @brief kinematic correction
   * @details ratio of sigma_los^2/sigma_tot^2
   *
   * @param Pot pointer to a Potential_JS instance
   * @return kinematic ratio
   */
  double sigma_corr(Potential_JS *Pot);
  /**
   * @brief ratio of sigma_xx^2 to sigma_zz^2
   *
   * @param Pot pointer to a Potential_JS instance
   * @return kinematic ratio of sigma_x^2/sigma_z^2
   */
  double sigma_x_sigma_z(Potential_JS *Pot);
};
// ============================================================================
struct obs_density_st{
  ObservedDensityProfile *ODP;
  VecDoub x2min,x2max;
  double D;
  /**
   * @brief structure for J-factor and D-factor at arbitrary angle integration
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param DP pointer to DensityProfile
   * @param D distance
   */
  obs_density_st(VecDoub x2min,VecDoub x2max,ObservedDensityProfile *ODP,double D)
    :x2min(x2min),x2max(x2max),ODP(ODP),D(D){}
};

struct obs_mass_density_st{
  ObservedTriaxialDensityProfile *OTDP;
  VecDoub x2min,x2max;
  double D;
    /**
   * @brief structure for mass inside ellipse integration
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param DP pointer to DensityProfile
   * @param D distance
   */
  obs_mass_density_st(VecDoub x2min,VecDoub x2max,ObservedTriaxialDensityProfile *OTDP,double D)
    :x2min(x2min),x2max(x2max),OTDP(OTDP),D(D){}
};
struct sig_st{
  ObservedTriaxialDensityProfile *D;
  VecDoub x2min,x2max;
  Potential_JS *Pot;
  int los; // x = 0, y = 1, z = 2
  /**
   * @brief structure for kinematic ratios
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param D pointer to ObservedTriaxialDensityProfile
   * @param Pot pointer to Potential_JS instance
   * @params los == component (x = 0, y = 1, z = 2)
   */
  sig_st(VecDoub x2min,VecDoub x2max,ObservedTriaxialDensityProfile *D,Potential_JS *Pot,int los)
    :x2min(x2min),x2max(x2max),D(D),Pot(Pot),los(los){}
};
#endif
// ============================================================================
