# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

BarryGoldman_cpp <- function(keyRotorsR, keyTimes, times) {
    .Call(`_qsplines_BarryGoldman_cpp`, keyRotorsR, keyTimes, times)
}

DeCasteljau_cpp <- function(rsegments, keyTimes, times) {
    .Call(`_qsplines_DeCasteljau_cpp`, rsegments, keyTimes, times)
}

control_points_cpp <- function(keyTimesR, keyRotorsR, closed, t, c, b) {
    .Call(`_qsplines_control_points_cpp`, keyTimesR, keyRotorsR, closed, t, c, b)
}

DeCasteljau_cpp2 <- function(rsegments, keyTimes, nintertimes) {
    .Call(`_qsplines_DeCasteljau_cpp2`, rsegments, keyTimes, nintertimes)
}

