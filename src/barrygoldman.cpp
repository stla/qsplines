#include "qsplines.h"

qtrn _eval_bg_single(
    double t,
    std::vector<qtrn> keyRotors,
    Rcpp::NumericVector keyTimes) {
  std::size_t idx = _check_time(t, keyTimes, true);
  qtrn q0 = keyRotors[idx];
  qtrn q1 = keyRotors[idx+1];
  double t0 = keyTimes[idx];
  double t1 = keyTimes[idx+1];
  const std::size_t nkeyRotors = keyRotors.size();
  const std::size_t nkeyTimes = keyTimes.size();
  qtrn q_1, q2;
  double t_1, t2;
  if(idx == 0) {
    q_1 = keyRotors[nkeyRotors - 2];
    if(quaternion::dot(q_1, q0) < 0){
      q_1 = -q_1;
    }
    t_1 = t0 - (keyTimes[nkeyTimes-1] - keyTimes[nkeyTimes - 2]);
  } else {
    q_1 = keyRotors[idx-1];
    t_1 = keyTimes[idx-1];
  }
  if(idx + 2 == nkeyRotors) {
    q2 = keyRotors[1];
    if(quaternion::dot(q1, q2) < 0){
      q2 = -q2;
    }
    t2 = t1 + (keyTimes[1] - keyTimes[0]);
  }else{
    q2 = keyRotors[idx+2];
    t2 = keyTimes[idx+2];
  }
  qtrn slerp_0_1 = slerp(q0, q1, (t - t0) / (t1 - t0));
  qtrn quat = slerp(
    slerp(
      slerp(q_1, q0, (t - t_1) / (t0 - t_1)),
      slerp_0_1,
      (t - t_1) / (t1 - t_1)
    ),
    slerp(
      slerp_0_1,
      slerp(q1, q2, (t - t1) / (t2 - t1)),
      (t - t0) / (t2 - t0)
    ),
    (t - t0) / (t1 - t0)
  );
  return quat;
}

std::vector<qtrn> _eval_bg_vector(
    Rcpp::NumericVector times,
    std::vector<qtrn> keyRotors,
    Rcpp::NumericVector keyTimes) {
  std::size_t n = times.size();
  std::vector<qtrn> quats(n);
  for(std::size_t i = 0; i < n; i++) {
    quats[i] = _eval_bg_single(times(i), keyRotors, keyTimes);
  }
  return quats;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix BarryGoldman_cpp(
  Rcpp::NumericMatrix keyRotorsR, Rcpp::NumericVector keyTimes, Rcpp::NumericVector times
) {
  std::vector<qtrn> keyRotors = _getRQuaternions(keyRotorsR);
  std::vector<qtrn> quats = _eval_bg_vector(times, keyRotors, keyTimes);
  return _getCQuaternions(quats);
} 
