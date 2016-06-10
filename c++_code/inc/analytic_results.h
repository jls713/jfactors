#ifndef ANALYTIC_H
#define ANALYTIC_H
#include <Python.h>
#include "utils_j.h"
#include <cmath>
#include <vector>
#include <string>
#include "../general/utils.h"
#include "density.h"
#include "observed_density.h"

class NFWDensityProfile{
	/**
	 * Class for analytic results for spherical NFW profile
	 */
private:
	double rho0; /** central density **/
	double   rs; /** scale radius    **/
public:
	NFWDensityProfile(double rho0, double rs);
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
	* @return J-factor
	*/
	double J_factor(double D, double ang);
	/**
	* @brief D-factor for beam with angle ang (in degrees) at distance D
	* @param D distance
	* @param ang beam width in degrees
	* @return D-factor
	*/
	double D_factor(double D, double ang);
};

class CoredDMProfile: public FiniteMassTriaxialDensityProfile, public Potential_JS{
private:
	double r_dm;   /** scale radius of dm       **/
    double q;/** density flattening dm **/
    double q_phi_dm;/** potential flattening dm **/
	double v0;     /** amplitude of dm potential**/
public:
	/**
	 * @brief Cored DM model
	 * @details Wyn Evans cored models -- based on Evans (1993)
	 *			The dark matter potential is flattened logarithmic.
	 *
	 * @param radii r_dm
	 * @param flattening q_dm -- in the density
	 * @param norms v0
	 */
	CoredDMProfile(double rs, double q, double v0);
	inline double amplitude(void){return v0;}
	inline double scale_radius(void){return r_dm;}
	inline double dm_density_flattening(void){return q;}
	inline double dm_flattening(void){return q_phi_dm;}
	/**
	 * @brief DM potential
	 *
	 * @param R cylindrical R
	 * @param z cylindrical z
	 *
	 * @return potential of DM
	 */
	double phi_dm(double R, double z);
	/**
	 * @brief DM forces
	 *
	 * @param X -- Cartesian vector
	 *
	 * @return forces due to DM
	 */
	VecDoub forces_dm(VecDoub X);
	/**
	 * @brief DM density
	 *
	 * @param R cylindrical R
	 * @param z cylindrical z
	 *
	 * @return density of DM
	 */
	double density_dm(double R, double z);
	/**
	 * @brief DM potential
	 *
	 * @param X -- Cartesian vector
	 *
	 * @return potential of DM
	 */
	double Phi(const VecDoub &X);
	/**
	 * @brief DM forces
	 *
	 * @param X -- Cartesian vector
	 *
	 * @return forces due to DM
	 */
	VecDoub Forces(const VecDoub &X);
	/**
	 * @brief DM density
	 *
	 * @param X Cartesian vector
	 *
	 * @return density of DM
	 */
	double rho(VecDoub X);
	/**
	 * @brief scale model
	 * @details scale radii by rr and velocities by vv
	 *
	 * @param rr radial scaling
	 * @param vv velocity scaling
	 */
	void scale(double rr, double vv);
};

class CoredModel{
private:
//	double n_star; /** slope of density profile **/
//	double r_star; /** scale radius of stars    **/
//	double r_dm;   /** scale radius of dm       **/
//	double q_star; /** flattening of stars      **/
  //  double q_phi_dm;/** potential flattening dm **/
//	double rho0;   /** density norm for stars   **/
//	double v0;     /** amplitude of dm potential**/

	CoredDMProfile dm;
	AlphaBetaGammaDensityProfile staaars;

public:
	/**
	 * @brief Cored model
	 * @details Wyn Evans cored models -- based on Evans (1993)
	 *			The dark matter potential is flattened logarithmic.
	 *			The stellar profile is cored.
	 *
	 * @param n_star outer slope of the stellar profile
	 * @param radii (r_stellar, r_dm)
	 * @param flattening (q_star, q_dm)
	 * @param norms (rho0, v0)
	 */
	CoredModel(double n_star, VecDoub radii, VecDoub flattening, VecDoub norms);
	inline double dm_flattening(void){return dm.dm_flattening();}
	inline double dm_density_flattening(void){return dm.dm_density_flattening();}
	inline double scale_radius(void){return staaars.scale_radius();}
	inline double dm_scale_radius(void){return dm.scale_radius();}
	inline double dm_norm(void){return dm.amplitude();}

	VecDoub forces_dm(VecDoub X){return dm.forces_dm(X);}
	double density_star(double R, double z);
	/**
	 * @brief mass within sphere of radius r
	 *
	 * @param r radius
	 * @return mass within sphere radius r
	 */
	double mass_dm(double r);
	/**
	 * @brief total stellar mass
	 *
	 * @return mass of stars
	 */
	double mass_star(void);
	/**
	 * @brief l.o.s velocity dispersion averaged over density
	 * @return l.o.s. velocity dispersion
	 */
	double average_los_dispersion(void);
	/**
	 * @brief half-light radius of light
	 * @return half-light radius
	 */
	double half_light_radius(void);
	/**
	 * @brief J factor for model viewed face-on
	 * @details J factor for model viewed face-on i.e. isophotes look circular
	 *
	 * @param D distance
	 * @param ang beam angle in degrees
	 *
	 * @return J-factor face-on
	 */
	double J_factor_face(double D, double ang);
	/**
	 * @brief D factor for model viewed face-on
	 * @details D factor for model viewed face-on i.e. isophotes look circular
	 *
	 * @param D distance
	 * @param ang beam angle in degrees
	 *
	 * @return D-factor face-on
	 */
	double D_factor_face(double D, double ang);
	/**
	 * @brief J factor for model viewed face-on very distant
	 * @details J factor for model viewed face-on i.e. isophotes look circular that is very distant
	 *
	 * @param D distance
	 *
	 * @return J-factor face-on v. distant
	 */
	double J_factor_face_asymptote(double Distance);
	/**
	 * @brief J factor for model viewed edge-on
	 * @details J factor for model viewed edge-on i.e. isophotes look elliptical
	 *
	 * @param D distance
	 * @param ang beam angle in degrees
	 *
	 * @return J-factor edge-on
	 */
	double J_factor_edge(double D, double ang);
	double J_factor_arbitrary(double D, double ang, double theta, double phi);
	double D_factor_arbitrary(double D, double ang, double theta, double phi);
	/**
	 * @brief Kinematic ratio
	 * @details sigma_tot^2/sigma_los^2
	 *
	 * @param theta latitudinal viewing angle (spherical polar)
	 * @param phi azimuthal viewing angle (spherical polar)
	 *
	 * @return kinematic ratio
	 */
	double kinematic_ratio(double theta, double phi);
	/**
	 * @brief l.o.s. velocity dispersion
	 * @details sigma_los
	 *
	 * @param theta latitudinal viewing angle (spherical polar)
	 * @param phi azimuthal viewing angle (spherical polar)
	 *
	 * @return l.o.s. velocity dispersion
	 */
	double sigma_los(double theta, double phi);
	double sigma_los_m(double theta, double phi);
	/**
	 * @brief scale model
	 * @details scale model to match observables half-light radius and
	 * l.o.s. velocity dispersion
	 *
	 * @param rh half-light radius
	 * @param slos l.o.s. velocity dispersion
	 * @param dir viewing direction -- face-on (round) or edge-on (edge)
	 *
	 */
	void scale(double rh, double slos, std::string dir = "round");
    /**
   * @brief compute mass
   * @details compute mass in cylinder or sphere
   *
   * @param theta spherical polar latitude viewing angle
   * @param phi spherical polar azimuthal viewing angle
   * @param radius
   * @param gobby print progress flag
   * @param typ mass inside cylinder ('cylinder'), inside ellipsoid
   * ('ellipsoid') or inside sphere ('sphere')
   * @return mass
   */
   double Mass(double theta, double phi, double r, std::string typ="cylinder");
};

struct sig_cored_st{
  CoredModel *C;
  VecDoub x2min,x2max;
  double theta, phi;
  /**
   * @brief structure for kinematic ratios of cored models
   *
   * @param x2min vector of lower integration limits
   * @param x2max vector of upper integration limits
   * @param C pointer to CoredModel
   * @param theta latitudinal spherical polar viewing angle
   * @param phi     azimuthal spherical polar viewing angle
   */
  sig_cored_st(VecDoub x2min,VecDoub x2max,CoredModel *C,double theta,double phi)
    :x2min(x2min),x2max(x2max),C(C),theta(theta),phi(phi){}
};
#endif
