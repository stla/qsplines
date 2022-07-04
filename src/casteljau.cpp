#include "qsplines.h"

std::vector<qtrn> _select_segment_and_normalize_t(
    std::vector<std::vector<qtrn>> segments,
    Rcpp::NumericVector keyTimes,
    double t,
    double* time,
    double* difftime) {
  const std::size_t idx = _check_time(t, keyTimes, false);
  const double t0 = keyTimes[idx];
  const double t1 = keyTimes[idx + 1];
  const double delta_t = t1 - t0;
  *difftime = delta_t;
  *time = (t - t0) / delta_t;
  return segments[idx];
}

qtrn _getRQuaternion(Rcpp::NumericVector qR) {
  qtrn quat(qR(0), qR(1), qR(2), qR(3));
  return quat;
}

std::vector<qtrn> _getRQuaternions(Rcpp::NumericMatrix Q) {
  std::size_t n = Q.ncol();
  std::vector<qtrn> quats(n);
  for(std::size_t j = 0; j < n; j++) {
    quats[j] = _getRQuaternion(Q(Rcpp::_, j));
  }
  return quats;
}

std::vector<std::vector<qtrn>> _getRSegments(Rcpp::List rsegments) {
  std::size_t nsegments = rsegments.size();
  std::vector<std::vector<qtrn>> segments(nsegments);
  for(std::size_t i = 0; i < nsegments; i++) {
    Rcpp::NumericMatrix segment = Rcpp::as<Rcpp::NumericMatrix>(rsegments(i));
    segments[i] = _getRQuaternions(segment);
  }
  return segments;
}

Rcpp::NumericVector _getCQuaternion(qtrn quat) { 
  return Rcpp::NumericVector::create(quat.a(), quat.b(), quat.c(), quat.d());
} 

Rcpp::NumericMatrix _getCQuaternions(std::vector<qtrn> quats) { 
  std::size_t n = quats.size();
  Rcpp::NumericMatrix Q(4, n);
  for(std::size_t j = 0; j < n; j++) {
    Rcpp::NumericVector qR = _getCQuaternion(quats[j]);
    Q(Rcpp::_, j) = qR;
  }
  return Q;
} 

std::vector<qtrn> _reduce_de_casteljau(std::vector<qtrn> segment, double t) {
  size_t l = segment.size();
  if(l < 2) {
    Rcpp::stop("Segment must have at least two quaternions.");
  }
  while(l > 2) {
    std::vector<qtrn> newsegment(l - 1);
    for(std::size_t i = 0; i < l - 1; i++) {
      qtrn one = segment[i];
      qtrn two = segment[i + 1];
      newsegment[i] = slerp(one, two, t);
    }
    segment = newsegment;
    l--;
  }
  return segment;
}

qtrn _eval_casteljau_single(double t,
                            std::vector<std::vector<qtrn>> segments,
                            Rcpp::NumericVector keyTimes) {
  double time, difftime;  // difftime not used here
  std::vector<qtrn> segment =
      _select_segment_and_normalize_t(segments, keyTimes, t, &time, &difftime);
  std::vector<qtrn> quats = _reduce_de_casteljau(segment, time);
  return slerp(quats[0], quats[1], time);
}

std::vector<qtrn> _eval_casteljau_vector(
    Rcpp::NumericVector times,
    std::vector<std::vector<qtrn>> segments,
    Rcpp::NumericVector keyTimes) {
  std::size_t n = times.size();
  std::vector<qtrn> quats(n);
  for(std::size_t i = 0; i < n; i++) {
    quats[i] = _eval_casteljau_single(times(i), segments, keyTimes);
  }
  return quats;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix DeCasteljau_cpp(
  Rcpp::List rsegments, Rcpp::NumericVector keyTimes, Rcpp::NumericVector times
){
  std::size_t nsegments = rsegments.size();
  std::size_t nkeyTimes = keyTimes.size();
  if(nkeyTimes == 0){
    keyTimes = _seq_len(nsegments + 1);
  }else if(nkeyTimes != nsegments + 1){
    Rcpp::stop("Number of key times must be one more than number of segments.");
  }
  std::vector<std::vector<qtrn>> segments = _getRSegments(rsegments);
  std::vector<qtrn> quats = _eval_casteljau_vector(times, segments, keyTimes);
  return _getCQuaternions(quats);
}

// {} []
// // [[Rcpp::export]]


// // [[Rcpp::export]]
// void test(){
//   qtrn q1(1.0/sqrt(2),1.0/sqrt(2), 0.0, 0.0);
//   qtrn q2(0.0, 0.0, 0.0, 1.0);
//   qtrn q3(0.5, 0.5, 0.5, 0.5);
//   std::vector<qtrn> quats = {q1, q2, q3};
//   std::array<qtrn, 2> arr = _calculate_control_quaternions(quats, 0.2, 0.4, 0.6, 0.1, 0.1, 0.1);
//   qtrn qq1 = arr[0];
//   qtrn qq2 = arr[1];
//   Rcpp::Rcout << qq1.w() << " " << qq1.x() << " " << qq1.y() << " " << qq1.z() << " --- ";
//   Rcpp::Rcout << qq2.w() << " " << qq2.x() << " " << qq2.y() << " " << qq2.z() << " --- ";
// }

// {} []
// // [[Rcpp::export]]