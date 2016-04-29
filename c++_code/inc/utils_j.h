//=============================================================================
#ifndef UTILSJ_H
#define UTILSJ_H

#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/exception_translator.hpp>
#include <exception>
#include <boost/python/numeric.hpp>
#include <boost/python/list.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/docstring_options.hpp>
#include <numpy/arrayobject.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
using namespace boost::python;

#include "../general/cuba/cuba.h"
#include <iostream>
#include "../general/cuba/cuba.h"
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include "../general/utils.h"
#include "potential.h"

const double Grav = 4.300918e-6; // in units solar mass, km/s kpc

//=============================================================================
// General integration routine
//=============================================================================

const std::string INTEG = "Divonne";
const int nproc = 1;
const int SEED = time(0);

template<class c>
double integrate(integrand_t integrand, c *P, double IE, double AE, std::string type, double *err){

    int neval,fail,nregions;
    double integral[1],error[1],prob[1];
    int NSIZE = P->x2min.size();
    double prod = 1.;
    for(int i=0;i<NSIZE;i++)prod*=(P->x2max[i]-P->x2min[i]);

    if(type=="Vegas")
        Vegas(NSIZE,nproc,integrand,P,1,IE,AE,0,SEED,
        MINEVAL,MAXEVAL,NSTART,NINCREASE,NBATCH,GRIDNO,STATEFILE,SPIN,
        &neval,&fail,integral,error,prob);

    else if (type=="Suave")
        Suave(NSIZE,nproc,integrand,P,1,IE,AE,0,SEED,
        MINEVAL,MAXEVAL,NNEW,FLATNESS,STATEFILE,SPIN,&nregions,
        &neval,&fail,integral,error,prob);

    else if (type=="Cuhre")
        Cuhre(NSIZE,nproc,integrand,P,1,IE,AE,0,
        MINEVAL, MAXEVAL, 0, STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);

    else
        Divonne(NSIZE,nproc,integrand,P,1,IE,AE,0,SEED,
        MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
        BORDER, MAXCHISQ, MINDEVIATION,
        NGIVEN, LDXGIVEN, nullptr, NEXTRA, nullptr,STATEFILE,SPIN,
        &nregions, &neval, &fail, integral, error, prob);
    if(err)*err=prod*error[0];
    if(fail!=0)
      throw std::runtime_error("Error: Required accuracy not reached.");
    return prod*integral[0];
}

//=============================================================================
// Vector Python interface
//=============================================================================

template <class T>
struct vector_to_ndarray {
  static PyObject* convert(const std::vector<T>& v) {
    list l;
    typename std::vector<T>::const_iterator p;
    for(p=v.begin();p!=v.end();++p){
      l.append(object(*p));
    }
    return incref(numeric::array(l).ptr());
  }
};

template <typename T>
struct vector_from_ndarray {
  vector_from_ndarray() {
    converter::registry::push_back(&vector_from_ndarray<T>::convertible,
                                   &vector_from_ndarray<T>::construct,
                                   type_id<std::vector<T>>());
  }

  // Determine if obj_ptr can be converted in a std::vector<T>
  static void* convertible(PyObject* obj_ptr) {
    if (!PyArray_Check(obj_ptr)) {
      throw std::invalid_argument("You have passed a non-numpy array");
      return 0;
    }
    return obj_ptr;
  }

  // Convert obj_ptr into a std::vector<T>
  static void construct(PyObject* obj_ptr,
                        converter::rvalue_from_python_stage1_data* data) {

    list l(handle<>(borrowed(obj_ptr)));
    // Grab pointer to memory into which to construct the new std::vector<T>
    void* storage = ((converter::rvalue_from_python_storage<std::vector<T>>*)
                     data)->storage.bytes;
    // in-place construct the new std::vector<T> using the character data
    // extraced from the python object
    std::vector<T>& v = *(new (storage) std::vector<T>());
    // populate the vector from list contains !!!
    int le = len(l);
    v.resize(le);
    for (int i = 0; i != le; ++i) {
      v[i] = extract<T>(l[i]);
    }
    // Stash the memory chunk pointer for later use by boost.python
    data->convertible = storage;
  }
};

#endif
