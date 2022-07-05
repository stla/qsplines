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


// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_calculate_control_quaternions(
  Rcpp::NumericMatrix Rquaternions, Rcpp::NumericVector times, 
  double t, double c, double b
){
  std::vector<qtrn> quaternions = _getRQuaternions(Rquaternions);
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
  std::vector<qtrn> out = {
    quaternion::exp(_getRQuaternion((t_1 - t0) * v0CD / 3.0)) * q0,
    quaternion::exp(_getRQuaternion((t1 - t0) * v0AB / 3.0)) * q0
  };
  return _getCQuaternions(out);
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

// template<typename T>
// std::vector<std::vector<T>> makeTriplets(std::vector<T> vec) {
//   const std::size_t n = vec.size();
//   std::vector<std::vector<T>> triplets(n-2);
//   for(size_t i = 0; i < n - 2; i ++) {
//     triplets[i] = {vec[i], vec[i+1], vec[i+2]};
//   }
//   return triplets;
// }

// // [[Rcpp::export]]
// Rcpp::NumericMatrix makeTriplets_times(Rcpp::NumericVector Rtimes, bool closed) {
//   std::vector<double> times(Rtimes.begin(), Rtimes.end());
//   if(closed){
//     const std::size_t ntimes = times.size();
//     double t1 = times[ntimes-1] + (times[1] - times[0]);
//     times.insert(times.begin(), times[0] - (times[ntimes-1] - times[ntimes-2]));
//     times.push_back(t1);
//   }
//   std::vector<std::vector<double>> triplets = makeTriplets<double>(times); 
//   const std::size_t n = triplets.size();
//   Rcpp::NumericMatrix out(n, 3);
//   for(size_t i = 0; i < n; i++){
//     std::vector<double> vi = triplets[i];
//     out(i, Rcpp::_) = Rcpp::NumericVector::create(vi[0], vi[1], vi[2]);
//   }
//   return out;
// }

// // [[Rcpp::export]]
// Rcpp::List makeTriplets_rotors(Rcpp::NumericMatrix Rrotors, bool closed) {
//   std::vector<qtrn> rotors = _getRQuaternions(Rrotors);
//   if(closed){
//     const std::size_t nrotors = rotors.size();
//     qtrn prefix = rotors[nrotors - 2];
//     if(quaternion::dot(prefix,rotors[0]) < 0.0){
//       prefix = -prefix;
//     }
//     qtrn suffix = rotors[1];
//     if(quaternion::dot(suffix, rotors[nrotors-1]) < 0.0){
//       suffix = -suffix;
//     }
//     rotors.insert(rotors.begin(), prefix);
//     rotors.push_back(suffix);
//   }
//   std::vector<std::vector<qtrn>> triplets = makeTriplets<qtrn>(rotors);
//   const std::size_t n = triplets.size();
//   Rcpp::List out(n);
//   for(size_t i = 0; i < n; i++){
//     out(i) = _getCQuaternions(triplets[i]);
//   }
//   return out;
// }


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