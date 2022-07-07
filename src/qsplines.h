#ifndef _QSHEADER_
#define _QSHEADER_

#include <Rcpp.h>
#include "quaternion.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
using namespace boost::math::quadrature;
#include <boost/math/tools/roots.hpp>

typedef quaternion::Quaternion<double> qtrn;



std::size_t _check_time(double, Rcpp::NumericVector, bool);

Rcpp::NumericVector _seq_len(std::size_t);

qtrn slerp(qtrn, qtrn, double);

Rcpp::NumericVector _interpolateTimes(Rcpp::NumericVector, std::size_t, bool);

#endif
