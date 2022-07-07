.getQMatrix <- function(quaternions){
  stopifnot(is.quaternion(quaternions))
  as.matrix(quaternions)
}

.isNumericVector <- function(x){
  is.numeric(x) && !anyNA(x)
}

.isBoolean <- function(x){
  is.logical(x) && length(x) == 1L && !is.na(x)
}

#' @title Spline using the De Casteljau algorithm
#' @description Constructs a quaternions spline using the De Casteljau
#'   algorithm.
#'
#' @param segments a list of vectors of unit quaternions; each segment must
#'   contain at least two quaternions
#' @param keyTimes the times corresponding to the segment boundaries, an
#'   increasing vector of length \code{length(segments)+1}; if \code{NULL},
#'   it is set to \code{1, 2, ..., length(segments)+1}
#' @param n_intertimes a positive integer used to linearly interpolate the 
#'   times given in \code{keyTimes} in order that there are
#'   \code{n_intertimes - 1} between two key times (so one gets the key
#'   times if \code{n_intertimes = 1}); this parameter must be given if 
#'   \code{constantSpeed=TRUE} and if it is given when 
#'   \code{constantSpeed=FALSE}, then it has precedence over \code{times}
#' @param times the interpolating times, they must lie within the range of
#'   \code{keyTimes}; ignored if \code{constantSpeed=TRUE} or if 
#'   \code{n_intertimes} is given
#' @param constantSpeed Boolean, whether to re-parameterize the spline to
#'   have constant speed; in this case, \code{"times"} is ignored and a
#'   function is returned, with an attribute \code{"times"}, the vector of
#'   new times corresponding to the key rotors
#'
#' @return A vector of quaternions whose length is the number of interpolating 
#'   times.
#' @export
#' 
#' @note This algorithm is rather for internal purpose. It serves for example
#'   as a base for the \link[=KochanekBartels]{Konachek-Bartels} algorithm.
DeCasteljau <- function(
  segments, keyTimes = NULL, n_intertimes, times, constantSpeed = FALSE
){
  stopifnot(is.list(segments))
  stopifnot(is.null(keyTimes) || .isNumericVector(keyTimes))
  stopifnot(missing(n_intertimes) || .isPositiveInteger(n_intertimes))
  stopifnot(missing(times) || .isNumericVector(times))
  stopifnot(.isBoolean(constantSpeed))
  n_segments <- length(segments)
  if(is.null(keyTimes) && missing(times)){
    keyTimes <- seq_len(n_segments + 1L)
  }else if(length(keyTimes) != n_segments + 1L){
    stop("Number of key times must be one more than number of segments.")
  }
  if(constantSpeed){
    if(missing(n_intertimes)){
      stop("With `constantSpeed=TRUE`, you must supply `n_intertimes`.")
    }
  }else{
    if(missing(times) && missing(n_intertimes)){
      stop(
        "With `constantSpeed=FALSE`, you must supply `n_intertimes` or `times`."
      )
    }
  }
  if(!constantSpeed && !missing(n_intertimes)){
    times <- interpolateTimes(keyTimes, n_intertimes, FALSE)
  }
  segments <- lapply(segments, .getQMatrix)
  if(constantSpeed){
    Q <- DeCasteljau_cpp2(segments, keyTimes, n_intertimes)
  }else{
    Q <- DeCasteljau_cpp(segments, keyTimes, times) 
  }
  as.quaternion(Q)
}

.natural_control_quaternion <- function(outer, inner_control){
  (inner_control * onion_inverse(outer))^(1 / 2) * outer
}

#' @title Kochanek-Bartels quaternions spline
#' @description Constructs a quaternions spline by the Kochanek-Bartels
#'   algorithm.
#'
#' @param keyRotors a vector of unit quaternions (rotors) to be interpolated
#' @param keyTimes the times corresponding to the key rotors; must be an
#'   increasing vector of the same length a \code{keyRotors} if
#'   \code{endcondition = "natural"} or of length one more than number of key
#'   rotors if \code{endcondition = "closed"}
#' @param tcb a vector of three numbers respectively corresponding to tension,
#'   continuity and bias
#' @param times the times of interpolation; each time must lie within the range
#'   of the key times; this parameter can be missing if \code{keyTimes} is
#'   \code{NULL} and \code{n_intertimes} is not missing, and it is ignored if
#'   \code{constantSpeed=TRUE}
#' @param n_intertimes if given, this argument has precedence over \code{times};
#'   \code{keyTimes} can be \code{NULL} and \code{times} is constructed by
#'   linearly interpolating the key times such that there are
#'   \code{n_intertimes - 1} between two key times (so the times are the key
#'   times if \code{n_intertimes = 1})
#' @param endcondition start/end conditions, can be \code{"closed"} or
#'   \code{"natural"}
#' @param constantSpeed Boolean, whether to re-parameterize the spline to
#'   have constant speed; in this case, \code{"times"} is ignored and you 
#'   must set the interpolating times with the help of \code{n_intertimes}
#'
#' @return A vector of quaternions having the same length as the \code{times}
#'   vector.
#' @export
#'
#' @note The algorithm with constant speed is very slow.
#' 
#' @examples 
#' library(qsplines)
#' # Using a Kochanek-Bartels quaternions spline to construct 
#' #   a spherical curve interpolating some key points on the 
#' #     sphere of radius 5
#'     
#' # helper function: spherical to Cartesian coordinates
#' sph2cart <- function(rho, theta, phi){
#'   return(c(
#'     rho * cos(theta) * sin(phi),
#'     rho * sin(theta) * sin(phi),
#'     rho * cos(phi)
#'   ))
#' }
#'
#' # construction of the key points on the sphere
#' keyPoints <- matrix(nrow = 0L, ncol = 3L)
#' theta_ <- seq(0, 2*pi, length.out = 9L)[-1L]
#' phi <- 1.3
#' for(theta in theta_){
#'   keyPoints <- rbind(keyPoints, sph2cart(5, theta, phi))
#'   phi = pi - phi
#' }
#' n_keyPoints <- nrow(keyPoints)
#'
#' # construction of the key rotors; the first key rotor 
#' #   is the identity quaternion and rotor i sends the 
#' #     key point i-1 to the key point i
#' keyRotors <- quaternion(length.out = n_keyPoints)
#' rotor <- keyRotors[1L] <- H1
#' for(i in seq_len(n_keyPoints - 1L)){
#'   keyRotors[i+1L] <- rotor <-
#'     quaternionFromTo(
#'       keyPoints[i, ]/5, keyPoints[i+1L, ]/5
#'     ) * rotor
#' }
#'
#' # Kochanek-Bartels quaternions spline
#' \donttest{rotors <- KochanekBartels(
#'   keyRotors, n_intertimes = 25L, 
#'   endcondition = "closed", tcb = c(-1, 5, 0)
#' )
#'
#' # construction of the interpolating points on the sphere
#' points <- matrix(nrow = 0L, ncol = 3L)
#' keyPoint1 <- rbind(keyPoints[1L, ])
#' for(i in seq_along(rotors)){
#'   points <- rbind(points, rotate(keyPoint1, rotors[i]))
#' }
#'
#' # visualize the result with the 'rgl' package
#' library(rgl)
#' spheres3d(0, 0, 0, radius = 5, color = "lightgreen")
#' spheres3d(points, radius = 0.2, color = "midnightblue")
#' spheres3d(keyPoints, radius = 0.25, color = "red")}
KochanekBartels <- function(
  keyRotors, keyTimes = NULL, tcb = c(0, 0, 0),
  times, n_intertimes, endcondition = "natural",
  constantSpeed = FALSE
){
  endcondition <- match.arg(endcondition, c("closed", "natural"))
  stopifnot(.isBoolean(constantSpeed))
  closed <- endcondition == "closed"
  keyRotors <- .check_keyRotors(keyRotors, closed)
  n_keyRotors <- length(keyRotors)
  keyTimes <- .check_keyTimes(keyTimes, n_keyRotors)
  if(!constantSpeed && !missing(n_intertimes)){
    stopifnot(.isPositiveInteger(n_intertimes))
    times <- interpolateTimes(keyTimes, n_intertimes, !closed)
    if(closed){
      # times <- head(times, -1L)
    }
  }else{
    times <- numeric(0L)
  }
  # if(is.null(keyTimes) && !missing(n_intertimes)){
  #   stopifnot(.isPositiveInteger(n_intertimes))
  #   times <- seq(
  #     1, n_keyRotors, length.out = n_intertimes * (n_keyRotors - 1L) + 1L
  #   )
  #   if(closed){
  #     times <- head(times, -1L)
  #   }
  # }
  # keyTimes <- .check_keyTimes(keyTimes, n_keyRotors)
  control_points <- as.quaternion(control_points_cpp(
    keyTimes, as.matrix(keyRotors), closed, tcb[1L], tcb[2L], tcb[3L]
  ))
  n_control_points <- length(control_points)
  if(closed){
    stopifnot(4*length(keyTimes) == n_control_points)
    control_points <- control_points[3L:(n_control_points-2L)]
  }else if(n_control_points == 0L){
    # two quaternions -> slerp
    stopifnot(n_keyRotors == 2L)
    stopifnot(length(keyTimes) == 2L)
    q0 <- keyRotors[1L]
    q1 <- keyRotors[2L]
    offset <- (q1 * onion_inverse(q0))^(1/3)
    control_points <- c(q0, offset*q0, onion_inverse(offset)*q1, q1)
  }else{ # natural
    control_points <- c(
      keyRotors[1L],
      .natural_control_quaternion(
        keyRotors[1L], control_points[1L]
      ),
      control_points,
      .natural_control_quaternion(
        keyRotors[n_keyRotors], control_points[n_control_points]
      ),
      keyRotors[n_keyRotors]
    )
  }
  segments <- lapply(4L*seq_len(length(control_points) %/% 4L)-3L, function(i){
    control_points[c(i, i+1L, i+2L, i+3L)]
  })
  DeCasteljau(
    segments, keyTimes, n_intertimes, times, constantSpeed
  )
}
