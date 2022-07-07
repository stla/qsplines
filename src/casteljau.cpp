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

std::array<qtrn, 2> _calculate_control_quaternions(
  std::array<qtrn, 3> quaternions, std::array<double, 3> times, 
  double t, double c, double b
){
  qtrn q_1 = quaternions[0];
  qtrn q0  = quaternions[1];
  qtrn q1  = quaternions[2];
  double t_1 = times[0];
  double t0  = times[1];
  double t1  = times[2];
  // if(t_1 == t0 || t0 == t1){ // 
  //   return {q0, q0};
  // }
  double A = (1 - t) * (1 + c) * (1 + b);
  double B = (1 - t) * (1 - c) * (1 - b);
  double C = (1 - t) * (1 - c) * (1 + b);
  double D = (1 - t) * (1 + c) * (1 - b);
  qtrn lq_in = quaternion::log(q0 * quaternion::inverse(q_1));
  qtrn lq_out = quaternion::log(q1 * quaternion::inverse(q0)); 
  Rcpp::NumericVector v_in = 
    Rcpp::NumericVector::create(0.0, lq_in.b(), lq_in.c(), lq_in.d());
  Rcpp::NumericVector v_out = 
    Rcpp::NumericVector::create(0.0, lq_out.b(), lq_out.c(), lq_out.d());
  v_in  = v_in / (t0 - t_1);
  v_out = v_out / (t1 - t0);
  Rcpp::NumericVector v0CD = 
    (C * (t1 - t0) * v_in + D * (t0 - t_1) * v_out) / (t1 - t_1);
  Rcpp::NumericVector v0AB = 
    (A * (t1 - t0) * v_in + B * (t0 - t_1) * v_out) / (t1 - t_1);
  std::array<qtrn, 2> out = {
    quaternion::exp(_getRQuaternion((t_1 - t0) * v0CD / 3.0)) * q0,
    quaternion::exp(_getRQuaternion((t1 - t0) * v0AB / 3.0)) * q0
  };
  return out;
}

template<typename T>
std::vector<std::array<T, 3>> makeTriplets(std::vector<T> vec) {
  const std::size_t n = vec.size();
  std::vector<std::array<T, 3>> triplets(n-2);
  for(size_t i = 0; i < n - 2; i ++) {
    triplets[i] = {vec[i], vec[i+1], vec[i+2]};
  }
  return triplets;
}

std::vector<std::array<double, 3>> makeTriplets_times(std::vector<double> times, bool closed) {
  if(closed){
    const std::size_t ntimes = times.size();
    double t1 = times[ntimes-1] + (times[1] - times[0]);
    times.insert(times.begin(), times[0] - (times[ntimes-1] - times[ntimes-2]));
    times.push_back(t1);
  }
  return makeTriplets<double>(times);
}

std::vector<std::array<qtrn, 3>> makeTriplets_rotors(std::vector<qtrn> rotors, bool closed) {
  if(closed){
    const std::size_t nrotors = rotors.size();
    qtrn prefix = rotors[nrotors - 2];
    if(quaternion::dot(prefix,rotors[0]) < 0.0){
      prefix = -prefix;
    }
    qtrn suffix = rotors[1];
    if(quaternion::dot(suffix, rotors[nrotors-1]) < 0.0){
      suffix = -suffix;
    }
    rotors.insert(rotors.begin(), prefix);
    rotors.push_back(suffix);
  }
  return makeTriplets<qtrn>(rotors);
}

qtrn _natural_control_quaternion(qtrn outer, qtrn inner_control){
  return quaternion::pow(
    inner_control * quaternion::inverse(outer), 0.5
  ) * outer;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix control_points_cpp(
  Rcpp::NumericVector keyTimesR, Rcpp::NumericMatrix keyRotorsR, bool closed, 
  double t, double c, double b
){
  std::vector<double> keyTimes(keyTimesR.begin(), keyTimesR.end());
  std::vector<qtrn> keyRotors = _getRQuaternions(keyRotorsR);
  std::vector<std::array<double, 3>> triplets_times = 
    makeTriplets_times(keyTimes, closed);
  std::vector<std::array<qtrn, 3>> triplets_rotors = 
    makeTriplets_rotors(keyRotors, closed);
  std::vector<qtrn> control_points(0);
  const std::size_t N = triplets_rotors.size(); 
  for(std::size_t i = 0; i < N; i++){ 
    std::array<qtrn, 3> qs = triplets_rotors[i];
    std::array<qtrn, 2> qb_qa = _calculate_control_quaternions(
      qs, triplets_times[i], t, c, b
    );
    qtrn q_before = qb_qa[0];
    qtrn q_after  = qb_qa[1];
    control_points.push_back(q_before);
    control_points.push_back(qs[1]);
    control_points.push_back(qs[1]);
    control_points.push_back(q_after);
  }
  return _getCQuaternions(control_points);
}

////////////////////////////////////////////////////////////////////////////////

double _eval2_casteljau_single(double t,
                               std::vector<std::vector<qtrn>> segments,
                               Rcpp::NumericVector keyTimes) {
  double time, difftime;  
  std::vector<qtrn> segment =
    _select_segment_and_normalize_t(segments, keyTimes, t, &time, &difftime);
  std::vector<qtrn> quats = _reduce_de_casteljau(segment, time);
  qtrn tangent = 
    quaternion::log(quats[1] * quaternion::inverse(quats[0]));
  std::size_t degree = segment.size() - 1;
  double lambda = (double)(2 * degree) / difftime;
  return sqrt(quaternion::unreal_norm_squared(lambda * tangent));
}

std::vector<qtrn> _eval2_casteljau_vector(
  std::vector<std::vector<qtrn>> segments, Rcpp::NumericVector keyTimes,
  std::size_t nintertimes
) {
  auto speed = [segments, keyTimes](double t) {
    return _eval2_casteljau_single(t, segments, keyTimes);
  };
  std::size_t nintervals = keyTimes.size() - 1;
  std::vector<double> integrated_speed(nintervals);
  for(std::size_t i = 0; i < nintervals; i++){
    double t0 = keyTimes[i];
    double t1 = keyTimes[i+1];
    double error;
    integrated_speed[i] = gauss_kronrod<double, 61>::integrate(speed, t0, t1, 0, 0, &error);
  }
  Rcpp::NumericVector newTimes(nintervals+1);
  newTimes[0] = 0;
  for(std::size_t i = 1; i <= nintervals; i++){
    newTimes[i] = newTimes[i-1] + integrated_speed[i-1];
  }
  Rcpp::Rcout << "LAST NEWTIME: " << newTimes(nintervals) << "  - ";
  Rcpp::NumericVector stimes = _interpolateTimes(newTimes, nintertimes, false);
  const std::size_t S = stimes.size();
  Rcpp::NumericVector times(S);
  for(std::size_t i = 0; i < S; i++){
    double s = stimes(i);
    if(s == newTimes(nintervals)){ // i.e dernier i ?
      Rcpp::Rcout << "FIRST IF; i = " << i << " --- ";
      times(i) = keyTimes(nintervals); // -1 ?
    }else if(s <= newTimes(0)){ // i.e. i=0
      Rcpp::Rcout << "SECOND IF; i = " << i << " --- ";
      times(i) = keyTimes(0);
    }else{
      const std::size_t idx = _check_time(s, newTimes, false);
      s = s - newTimes(idx);
      double t0 = keyTimes[idx];
      double t1 = keyTimes[idx+1];
      auto igspeed = [speed, t0, s](double t) {
        double error;
        return gauss_kronrod<double, 61>::integrate(speed, t0, t, 0, 0, &error) - s;
      };
      std::pair<double, double> root = boost::math::tools::bisect(
        igspeed,
        t0, t1,
        [](double l, double r){return abs(l-r) < 1e-12;}
      );
      times(i) = (root.first + root.second) / 2.0;
      Rcpp::Rcout << "TIMES(i):::::::::::::::: " << times(i) << " --- ";
    }
  }
  return _eval_casteljau_vector(times, segments, keyTimes);
} 

// [[Rcpp::export]]
Rcpp::NumericMatrix DeCasteljau_cpp2(
  Rcpp::List rsegments, Rcpp::NumericVector keyTimes, std::size_t nintertimes
){
  std::size_t nsegments = rsegments.size();
  std::size_t nkeyTimes = keyTimes.size();
  if(nkeyTimes == 0){
    keyTimes = _seq_len(nsegments + 1);
  }else if(nkeyTimes != nsegments + 1){
    Rcpp::stop("Number of key times must be one more than number of segments.");
  }
  std::vector<std::vector<qtrn>> segments = _getRSegments(rsegments);
  std::vector<qtrn> quats = _eval2_casteljau_vector(segments, keyTimes, nintertimes);
  return _getCQuaternions(quats);
}

// std::vector<qtrn> _eval2_casteljau_vector(
//     Rcpp::NumericVector times,
//     std::vector<std::vector<qtrn>> segments,
//     Rcpp::NumericVector keyTimes) {
//   std::size_t n = times.size();
//   std::vector<qtrn> quats(n);
//   for(std::size_t i = 0; i < n; i++) {
//     quats[i] = _eval2_casteljau_single(times(i), segments, keyTimes);
//   }
//   return quats;
// }

// {} []
// // [[Rcpp::export]]