// ============================================================================
#ifndef TWOCOMP_H
#define TWOCOMP_H
#include <Python.h>
#include "utils_j.h"
#include <cmath>
#include <vector>
#include <string>
#include "../general/utils.h"
#include "density.h"
#include "observed_density.h"

//=============================================================================
// Interface to tact multipole code
//=============================================================================

class MultipoleDensity: public Density{
/**
 * @brief Wrapper for density profile to multipole code
 * @details wraps DensityProfile to pass to Multipole code
 */
protected:
  DensityProfile *DP;
public:
  MultipoleDensity(DensityProfile *DP):DP(DP){}
  double density(const VecDoub &X){return DP->rho(X)*Grav;}
};

//=============================================================================
class DoubleProfileModel{
  AlphaBetaGammaDensityProfile *Stars; /* Density profile for the stars */
  AlphaBetaGammaDensityProfile *DM;    /* Density profile for the DM    */
  double rh; 						               /* half-light radius in kpc      */
  double slos;			       /* line-of-sight velocity dispersion in km/s */
  bool use_multipole;      /* option to use multipole expansion         */
  double rsrdmratio; 			          		 /* scale radius of stars to DM */
  double rtdmrdmratio,rtsrdmratio;/* tidal radius of stars and DM to DM scale*/
  VecDoub abg_st, abg_dm;		    /* (alpha,beta,gamma) of stars and DM   */
public:

  DoubleProfileModel(double ba, double ca, double rh=0.049, double slos=3.22, bool use_multipole=true, double rsrdmratio=0.5, double rtdmrdmratio = 10., double rtsrdmratio=9., VecDoub abg_st={2.,5.,0.}, VecDoub abg_dm = {1.,3.,1.});
  DoubleProfileModel(VecDoub axis_ratios_st, VecDoub axis_ratios_dm, double rh, double slos, bool use_multipole, double rsrdmratio, double rtdmrdmratio, double rtsrdmratio, VecDoub abg_st, VecDoub abg_dm);
  ~DoubleProfileModel(){delete Stars; delete DM;}
  /**
   * @brief nice print for model
   */
  void print(void);
  double observed_ellipticity(double theta, double phi);
  /**
   * @brief compute J-factor
   * @details [long description]
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   * @param D distance
   * @param ang beam angle in deg
   * @param gobby print progress flag
   * @param with_D also compute D
   * @return vector of J (and D if required)
   */
  VecDoub J_factor(double theta, double phi, double D, double ang, bool gobby = true, bool with_D = false);
  VecDoub J_factor_old(double theta, double phi, double D, double ang, bool gobby = true, bool with_D = false);
    /**
   * @brief compute mass
   * @details compute mass in cylinder or sphere
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   * @param D distance
   * @param ang vector of beam angle in deg
   * @param gobby print progress flag
   * @param typ mass inside cylinder ('cylinder'), inside ellipsoid
   * ('ellipsoid') or inside sphere ('sphere')
   * @param radius within which to compute velocity dispersion (if <0 over whole galaxy)
   * @return mass
   */
  VecDoub MassProfile(double theta, double phi, double D, VecDoub ang, bool gobby = true, std::string typ="cylinder", double veldispradius=-1.);
    /**
   * @brief compute J-factor correction factor
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   * @param D distance
   * @param ang beam angle in deg
   * @param gobby print progress flag
   * @param geo_average (compare to geometric average spherical version)
   * @return correction factor
   */
  double correction_factor(double theta, double phi, double D, double ang, bool gobby = true, bool geo_average=true);
  /**
   * @brief kinematic ratio
   * @details ratio of sigma_los/sigma_tot
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   *
   * @return kinematic ratio
   */
  double kinematic_ratio(double theta, double phi);
	/**
   * @brief ratio of sigma_x^2 to sigma_z^2
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   *
   * @return kinematic ratio
   */
  double sigma_x_sigma_z(void);
  /**
   * @brief ratio of sigma_x^2 to sigma_z^2 and sigma_y^2 to sigma_z^2
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   *
   * @return kinematic ratios
   */
  VecDoub sigma_x_y_sigma_z(void);
    /**
   * @brief compute ratio of observed half-light major axis length to model
   * @details [long description]
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle

   * @return radius ratio
   */
  double radius_ratio(double theta, double phi);
  /**
   * @brief compute ratio of observed velocity dispersion length to model
   * @details [long description]
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle

   * @return velocity ratio
   */
   double velocity_ratio(double theta, double phi);
   /**
    * @brief total velocity dispersion
    */
   double sigma_tot(void);
   /**
    * @brief spherical half-light radius of stellar distribution
    * @return spherical rh
    */
   double spherical_rh(void);
   /**
    * @brief projected half-light radius of stellar distribution
    * @return projected rh
    */
   double projected_rh(double,double);
};

class PaperModel: public DoubleProfileModel{
/**
 * @brief Specific DoubleProfileModel for the paper
 */
public:
  PaperModel(double ba, double ca, double rh=0.049, double slos=3.22, bool use_multipole=true):DoubleProfileModel(ba,ca,rh,slos,use_multipole){}
};

#endif
