// ============================================================================
#ifndef DENSITY_H
#define DENSITY_H
#include <Python.h>
#include "utils_j.h"
// ============================================================================

class DensityProfile{
  /**
   * General density profile
   * -- base class for computing J-factors and D-factors of
   * -- general density profile
   */
public:
  /**
   * @brief evaluate density profile
   * @details density at the Cartesian coordinate X (to be overrid)
   *
   * @param X Cartesian coordinate = (x,y,z)
   * @return density
   */
  virtual double rho(VecDoub X) = 0;
  /**
   * @brief scale radius
   * @details scale radius
   * @return scale radius
   */
  virtual double scale_radius(void){return 1.;}
  /**
   * @brief mass inside sphere of angle ang (in degrees) at distance D
   * @param D distance
   * @param ang beam width in degrees
   * @return mass inside sphere
   */
  double M_sphere(double D, double ang);
  /**
   * @brief J-factor for beam with angle ang (in degrees) at distance D
   * @param D distance
   * @param ang beam width in degrees
   * @param los string giving viewing line-of-sight ('x','y' or 'z')
   * @return J-factor
   */
  double J_factor(double D, double ang, std::string los="z");
    /**
   * @brief J-factor for beam with angle ang (in degrees) at large distance D
   * @details differs from J_factor as it assumes distance D is large
   * compared to size of cluster so beam is cylinder
   * @param D distance
   * @param ang beam width in degrees
   * @param los string giving viewing line-of-sight ('x','y' or 'z')
   * @return J-factor
   */
  double J_far_factor(double D, double ang, std::string los="z");
  /**
   * @brief D-factor for beam with angle ang (in degrees) at distance D
   * @param D distance
   * @param ang beam width in degrees
   * @param los string giving viewing line-of-sight ('x','y' or 'z')
   * @return D-factor
   */
  double D_factor(double D, double ang, std::string los="z");
      /**
   * @brief D-factor for beam with angle ang (in degrees) at large distance D
   * @details differs from D_factor as it assumes distance D is large
   * compared to size of cluster so beam is cylinder
   * @param D distance
   * @param ang beam width in degrees
   * @param los string giving viewing line-of-sight ('x','y' or 'z')
   * @return D-factor
   */
  double D_far_factor(double D, double ang, std::string los="z");
};
// ============================================================================

class TriaxialDensityProfile: public DensityProfile{
   /**
   * General triaxial density profile
   * -- base class for constructing triaxial density profiles described by
   * rho(m^2) where m^2=x^2/a^2+y^2/b^2+z^2/c^2
   */
protected:
  double a, b, c;       /* axis lengths m^2=x^2/a^2+y^2/b^2+z^2/c^2  */
  double ba, ca;        /* axis ratios b/a, c/a                      */
  double T;             /* triaxiality T                             */
  double abc;           /* product of axis lengths (1 for normalized)*/
    /**
   * @brief Triaxiality
   * @details compute triaxiality parameter given general axis ratios
   *
   * @param ba intermediate-to-major axis ratios
   * @param ca minor-to-major axis ratios
   *
   * @return triaxiality parameter
   */
  double Triaxiality(double ba, double ca);
  /**
   * @brief b/a ratio
   * @details computes intermediate to major axis ratio given triaxiality
   * parameter T and minor-to-major axis ratio
   *
   * @param T triaxiality parameter
   * @param ca minor-to-major axis ratio
   *
   * @return intermediate-to-major axis ratio
   */
  double b_over_a(double T, double ca);
    /**
   * @brief c/a ratio
   * @details computes intermediate to major axis ratio given triaxiality
   * parameter T and intermediate-to-major axis ratio
   *
   * @param T triaxiality parameter
   * @param ba minor-to-major axis ratio
   *
   * @return minor-to-major axis ratio
   */
  double c_over_a(double T, double ba);
public:
  /**
   * @brief TriaxialDensityProfile constructor
   * @details initializes base class DensityProfile and sets axis ratios
   * (which are normalized if necessary)
   *
   * @param abc vector of axis lengths m^2=x^2/a^2+y^2/b^2+z^2/c^2
   * @param normalize bool to normalize such that abc=1
   */
  TriaxialDensityProfile(VecDoub abc, bool normalize=false);
  /**
   * @brief TriaxialDensityProfile constructor
   * @details initializes base class DensityProfile and sets axis ratios
   * (which are normalized if necessary)
   *
   * @param a ellipsoidal length in x direction
   * @param b ellipsoidal length in y direction
   * @param c ellipsoidal length in z direction
   * @param normalize bool to normalize such that abc=1
   */
  TriaxialDensityProfile(double a, double b, double c, bool normalize=false);
  /**
   * @brief evaluate density profile
   * @details density at the Cartesian coordinate X (to be overrid)
   *
   * @param X Cartesian coordinate = (x,y,z)
   * @return density
   */
  virtual double rho(VecDoub X) = 0;
  /**
   * @brief Triaxiality
   * @details compute triaxiality parameter given general axis ratios
   * @return triaxiality parameter
   */
  inline double Triaxiality(void) {return T;}
  /**
   * @brief b/a ratio
   * @details computes intermediate to major axis ratio given triaxiality
   * parameter T and minor-to-major axis ratio
   * @return intermediate-to-major axis ratio
   */
  inline double b_over_a(void) {return ba;}
    /**
   * @brief c/a ratio
   * @details computes intermediate to major axis ratio given triaxiality
   * parameter T and intermediate-to-major axis ratio
   * @return minor-to-major axis ratio
   */
  inline double c_over_a(void) {return ca;}
  /**
   * @brief axis ratios
   * @details gives intermediate-to-major axis ratios and minor-to-major
   * axis ratios
   * @return vector of (b/a,c/a)
   */
  inline VecDoub axis_ratios(){return {ba,ca};}
  /**
   * @brief axes lengths
   * @details ellipsoidal axes lengths
   * @return (a,b,c)
   */
  inline VecDoub axes(void){return {a,b,c};}
};
// ============================================================================
class FiniteMassTriaxialDensityProfile: public TriaxialDensityProfile{
  /**
   * @brief Triaxial finite mass density profiles
   * @details density profile for finite mass profiles.
   */
public:
  FiniteMassTriaxialDensityProfile(VecDoub abc, bool normalize)
    :TriaxialDensityProfile(abc,normalize){}
  /**
   * @brief mass
   * @details total mass of model
   * @return mass
   */
  double mass(void);
  /**
   * @brief ellipsoidal half-light radius
   * @details  ellipsoidal radius (with same axis ratios as model) that
   * contains half the total mass
   * @return ellipsoidal half-light radius
   */
  double spherical_half_light_radius(void);
  /**
   * @brief total velocity dispersion
   * @details compute total velocity dispersion of density profile in
   * given potential
   *
   * @param pointer to potential
   * @return total velocity dispersion
   */
  double sigma_tot(Potential_JS *Pot);
  /**
   * @brief mass inside an ellipsoid
   * @details mass in ellipsoid (with axis ratios b/a and c/a) of major axis
   * length r
   *
   * @param r major axis length
   * @return mass inside ellipsoid
   */
  double mass_ellipsoid(double r);
  /**
   * @brief mass inside an ellipsoid
   * @details mass in ellipsoid (with axis ratios ar) of major axis
   * length r
   *
   * @param r major axis length
   * @return mass inside ellipsoid
   */
  double mass_ellipsoid(double r, VecDoub ar);
};

// ============================================================================

class AlphaBetaGammaDensityProfile: public FiniteMassTriaxialDensityProfile{
  /**
   * @brief Triaxial double power-law density profiles
   * @details density profile for a broken double power-law profile.
   * The density is a function of the ellipsoidal radius
   * m^2=x^2/a^2+y^2/b^2+z^2/c^2 and of the form
   * rho(m) = rho0 (m/rs)^-gamma (1+(m/rs)^alpha)^-(beta-gamma)/alpha
   *          sech^2(m/rt)
   */
protected:
  double alpha, beta, gamma; /* slopes of density profile              */
  double rho0;               /* central density, density normalization */
  double rs, rt;             /* scale radius and tidal radius          */
public:
  /**
   * @brief Triaxial alpha,beta,gamma density profile constructor
   *
   * @param abg (alpha,beta,gamma)
   * @param rho0 'central density', overall normalization
   * @param rs scale radius
   * @param rt tidal radius
   * @param abc axis lengths m^2=x^2/a^2+y^2/b^2+z^2/c^2
   * @param normalize normalize model such that abc=1
   */
  AlphaBetaGammaDensityProfile(VecDoub abg, double rho0, double rs, double rt, VecDoub abc, bool normalize=false);
  /**
   * @brief Triaxial alpha,beta,gamma density profile constructor
   *
   * @param alpha
   * @param beta
   * @param gamma
   * @param rho0 'central density', overall normalization
   * @param rs scale radius
   * @param rt tidal radius
   * @param abc axis lengths m^2=x^2/a^2+y^2/b^2+z^2/c^2
   * @param normalize normalize model such that abc=1
   */
  AlphaBetaGammaDensityProfile(double alpha, double beta, double gamma, double rho0, double rs, double rt, VecDoub abc, bool normalize=false);
  /**
   * @brief density
   * @details rho(m) = rho0 (m/rs)^-gamma (1+(m/rs)^alpha)^-(beta-gamma)/alpha
   *          sech^2(m/rt)
   *
   * @param x Cartesian vector = (x,y,z)
   * @return density
   */
  double rho(VecDoub x);
  /** Model properties **/
  virtual double scale_radius(void){return rs;}
  virtual double tidal_radius(void){return rt;}
  virtual double central_density(void){return rho0;}
  VecDoub alpha_beta_gamma(void){return {alpha,beta,gamma};}
  /**
   * @brief print properties of model
   */
  void print(void);
};

// ============================================================================

struct densityJF_st{
  DensityProfile *DP;
  VecDoub x2min,x2max;
  double D;
  std::string los;
  /**
   * @brief structure for mass, J-factor and D-factor integration for
   * DensityProfile
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param DP pointer to DensityProfile
   * @param D distance
   * @param los line-of-sight string (='x','y','z')
   */
  densityJF_st(VecDoub x2min,VecDoub x2max,
             DensityProfile *DP,double D,std::string los)
    :x2min(x2min),x2max(x2max),DP(DP),D(D),los(los){}
};

struct mass_density_st{
  FiniteMassTriaxialDensityProfile *DP;
  VecDoub x2min,x2max;
  double D;
    /**
   * @brief structure for total mass for FiniteMassTriaxialDensityProfile
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param DP pointer to FiniteMassDensityProfile
   * @param D distance
   */
  mass_density_st(VecDoub x2min,VecDoub x2max,
                  FiniteMassTriaxialDensityProfile *DP,double D)
    :x2min(x2min),x2max(x2max),DP(DP),D(D){}
};
struct mass_density_diff_shape_st: mass_density_st{
  VecDoub ar;
  /**
   * @brief structure for mass inside ellipsoid of arbitrary shape for FiniteMassTriaxialDensityProfile
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param DP pointer to FiniteMassDensityProfile
   * @param D distance
   * @param ar axis ratios of ellipsoid
   */
  mass_density_diff_shape_st(VecDoub x2min,VecDoub x2max,
                  FiniteMassTriaxialDensityProfile *DP,double D,VecDoub ar)
    :mass_density_st(x2min,x2max,DP,D),ar(ar){}
};
struct sigt_st{
  FiniteMassTriaxialDensityProfile *D;
  VecDoub x2min,x2max;
  Potential_JS *Pot;
  /**
   * @brief structure for velocity dispersion
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param D pointer to FiniteMassTriaxialDensityProfile
   * @param Pot pointer to Potential_JS instance
   */
  sigt_st(VecDoub x2min,VecDoub x2max,
          FiniteMassTriaxialDensityProfile *D,Potential_JS *Pot)
    :x2min(x2min),x2max(x2max),D(D),Pot(Pot){}
};
// ============================================================================


#endif
