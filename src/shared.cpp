#include "qsplines.h"

std::vector<qtrn> _canonicalized(std::vector<qtrn> quaternions) {
  const std::size_t n = quaternions.size();
  std::vector<qtrn> out(n);
  qtrn p(1.0, 0.0, 0.0, 0.0);
  qtrn zero(0.0, 0.0, 0.0, 0.0);
  for(std::size_t i = 0; i < n; i++) {
    qtrn q = quaternions[i];
    if(quaternion::dot(p, q) < 0.0) {
      q = -q;
    }
    out[i] = q;
    p = q;
  }
  return out;
}

Rcpp::NumericVector _seq_len(std::size_t n) {
  Rcpp::NumericVector seq(n);
  for(std::size_t i = 1; i <= n; i++) {
    seq(i-1) = (double)(i);
  }
  return seq;
}

std::size_t _findInterval(double x, Rcpp::NumericVector vec) {
  size_t n = vec.size();
  if(x > vec(n - 1)) {
    return n;
  }
  std::size_t idx = 0;
  for(std::size_t i = 0; i < n - 1; i++) {
    if(x >= vec(i)) {
      idx++;
    } else {
      break;
    }
  }
  return idx;
}

std::size_t _check_time(double t, Rcpp::NumericVector keyTimes, bool special) {
  std::size_t n_keyTimes = keyTimes.size();
  double lastKeyTime = keyTimes(n_keyTimes - 1);
  if(t < keyTimes(0) || t > lastKeyTime) {
    Rcpp::stop(
        "The interpolating times must be within the range of `keyTimes`.");
  }
  std::size_t idx;
  if(t < lastKeyTime) {
    idx = _findInterval(t, keyTimes) - 1;
  } else {  // t = lastKeyTime
    if(special) {
      idx = n_keyTimes - 3;  // 2?
    } else {
      idx = n_keyTimes - 2;  // 1?
    }
  }
  return idx;
}

qtrn slerp(qtrn q0, qtrn q1, double t) {
  qtrn q2 = q1 * quaternion::inverse(q0);
  qtrn q2powt = quaternion::pow(q2, t);
  return q2powt * q0;
}

 
