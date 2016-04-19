#include "Python.h"
#include "utils.h"
#include "../general/utils.h"
#include "density.h"
#include "observed_density.h"
#include "potential.h"
#include "Multipole.h"
#include "two_component_model.h"

BOOST_PYTHON_MODULE_INIT(jfactors_py) {
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  docstring_options doc_options(true);
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

  to_python_converter<VecDoub, vector_to_ndarray<double>>();
  vector_from_ndarray<double>();
  import_array();
}
