#include "Python.h"
#include "utils.h"
#include "../general/utils.h"
#include "density.h"
#include "observed_density.h"
#include "potential.h"
#include "Multipole.h"
#include "analytic_results.h"
#include "two_component_model.h"

void translate(std::exception const& e)
{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_RuntimeError, e.what());
}

BOOST_PYTHON_MODULE_INIT(jfactors_py) {
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  docstring_options doc_options(true);
  register_exception_translator<std::invalid_argument>(&translate);
  register_exception_translator<std::runtime_error>(&translate);
  class_<DoubleProfileModel, boost::noncopyable>("DoubleProfileModel",
          "Model for Sanders, Evans & Geringer-Sameth (2016)\n"
          "\n"
          "Args:\n"
          "    param1: ba axis ratio.\n"
          "    param2: ca axis ratio.\n"
          "    param3: rh.\n"
          "    param4: sigma_los.\n"
          "    param5: use multipole.\n"
          "    param6: rs/rdm.\n"
          "    param7: rtdm/rdm.\n"
          "    param8: rts/rdm.\n"
          "    param9: abg_st.\n"
          "    param10: abg_dm.\n"
          ,init<double, double, double, double, bool, double, double,double,VecDoub,VecDoub>())
    .def("J_factor", &DoubleProfileModel::J_factor,
         "Compute J factor\n"
          "\n"
          "Args:\n"
          "    param1: theta, viewing angle.\n"
          "    param2: phi, viewing angle.\n"
          "    param3: distance.\n"
          "    param4: ang, in degrees.\n"
          "    param5: print errors.\n"
          "    param6: return D-factor as well.\n"
          "\n"
          "Returns:\n"
          "    vector of J-factor (and D-factor if required)\n"
      "")
    .def("Mass", &DoubleProfileModel::MassProfile,
         "Compute mass\n"
          "\n"
          "Args:\n"
          "    param1: theta, viewing angle.\n"
          "    param2: phi, viewing angle.\n"
          "    param3: distance.\n"
          "    param4: ang, vector of angles in degrees.\n"
          "    param5: print errors.\n"
          "    param6: return inside cylinder ('cylinder'), spheroid ('spheroid') or sphere ('sphere').\n"
          "    param7: radius within which to compute dispersion"
          "\n"
          "Returns:\n"
          "    mass\n"
      "")
    .def("r_h", &DoubleProfileModel::spherical_rh,
         "Compute spherical half-light radius of stars\n"
          "\n"
          "Returns:\n"
          "    r_h\n"
      "")
    .def("R_h", &DoubleProfileModel::projected_rh,
         "Compute projected half-light radius of stars\n"
          "\n"
          "Args:\n"
          "    param1: theta, viewing angle.\n"
          "    param2: phi, viewing angle.\n"
          "Returns:\n"
          "    R_h\n"
      "")
    .def("ellipticity", &DoubleProfileModel::observed_ellipticity,
         "Compute ellipticity\n"
          "\n"
          "Args:\n"
          "    param1: theta, viewing angle.\n"
          "    param2: phi, viewing angle.\n"
          "\n"
          "Returns:\n"
          "    ellipticity\n"
      "");
  class_<PaperModel, bases<DoubleProfileModel>>("PaperModel","Model for paper", init<double,double,double,double,bool>());

  class_<DensityProfile, boost::noncopyable>("DensityProfile",no_init);
  class_<TriaxialDensityProfile, boost::noncopyable>("TriaxialDensityProfile",no_init);
  class_<AlphaBetaGammaDensityProfile, bases<TriaxialDensityProfile>>("AlphaBetaGammaDensityProfile",
    "Model for Sanders, Evans & Geringer-Sameth (2016)\n"
          "\n"
          "Args:\n"
          "    param1: (alpha,beta,gamma).\n"
          "    param2: rho0.\n"
          "    param3: rs.\n"
          "    param4: rt.\n"
          "    param5: abc axis lengths.\n"
          "    param6: normalize?.\n",
    init<VecDoub,double,double,double,VecDoub,double>())
   .def("J_far_factor", &AlphaBetaGammaDensityProfile::J_far_factor,
         "Compute J factor in far-field approx\n"
          "\n"
          "Args:\n"
          "    param1: D, distance.\n"
          "    param2: ang, beam angle.\n"
          "    param3: los, line-of-sight.\n"
          "\n"
          "Returns:\n"
          "    J factor\n"
      "")
   .def("D_far_factor", &AlphaBetaGammaDensityProfile::D_far_factor,
         "Compute D_factor in far-field approx\n"
          "\n"
          "Args:\n"
          "    param1: D, distance.\n"
          "    param2: ang, beam angle.\n"
          "    param3: los, line-of-sight.\n"
          "\n"
          "Returns:\n"
          "    D factor\n"
      "");
  class_<CoredModel, boost::noncopyable>("CoredModel",
          "Evans isothermal cored model (1993)\n"
          "\n"
          "Args:\n"
          "    param1: n_star.\n"
          "    param2: radii (r_s, r_dm).\n"
          "    param3: flattening (q_s,q_dm).\n"
          "    param4: norms (rho0,v0).\n"
          ,init<double,VecDoub,VecDoub,VecDoub>())
  .def("scale",&CoredModel::scale,
       "Scale model to half-light radius and velocity dispersion"
       "\n"
       "Args:\n"
       "    param1: rh, half-light radius\n"
       "    param2: slos, velocity dispersion\n"
       "    param3: dir, viewing direction -- face-on (round) or edge-on(edge)\n"
       "    param4: vdrad, radius within which velocity dispersion measured\n")
  .def("mass_dm",&CoredModel::mass_dm,
       "Mass within sphere of radius r"
       "\n"
       "Args:\n"
       "    param1: r, spherical radius\n")
  .def("mass",&CoredModel::Mass,
       "Compute mass\n"
          "\n"
          "Args:\n"
          "    param1: theta, viewing angle.\n"
          "    param2: phi, viewing angle.\n"
          "    param3: r, radius.\n"
          "    param4: return inside ellipsoid ('ellipsoid') or sphere ('sphere').\n"
          "\n"
          "Returns:\n"
          "    mass\n"
      "");
  to_python_converter<VecDoub, vector_to_ndarray<double>>();
  vector_from_ndarray<double>();
  import_array();
}
