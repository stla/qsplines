#ifndef _QSHEADER_
#define _QSHEADER_

#include <Rcpp.h>
#include "quaternion.h"

typedef quaternion::Quaternion<double> qtrn;



std::size_t _check_time(double, Rcpp::NumericVector, bool);

Rcpp::NumericVector _seq_len(std::size_t);

qtrn slerp(qtrn, qtrn, double);

#endif
