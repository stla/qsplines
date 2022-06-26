library(onion)

sph2cart <- function(rho, theta, phi){
  return(c(
    rho * cos(theta) * sin(phi),
    rho * sin(theta) * sin(phi),
    rho * cos(phi)
  ))
}

cart2latlong <- function(x, y, z){
  cbind(asin(z / sqrt(x*x+y*y+z*z)), atan2(y, x)) * 180/pi
}
# construction of the key points on the sphere
keyPoints <- matrix(nrow = 0L, ncol = 3L)
theta_ <- seq(0, 2*pi, length.out = 9L)[-1L]
phi <- 1.3
for(theta in theta_){
  keyPoints <- rbind(keyPoints, sph2cart(5, theta, phi))
  phi = pi - phi
}
n_keyPoints <- nrow(keyPoints)
# construction of the key rotors; the first key rotor is the identity
#   quaternion and rotor i sends the key point i-1 to the key point i
keyRotors <- quaternion(length.out = n_keyPoints)
rotor <- keyRotors[1L] <- H1
for(i in seq_len(n_keyPoints - 1L)){
  keyRotors[i+1L] <- rotor <-
    quaternionFromTo(keyPoints[i, ]/5, keyPoints[i+1L, ]/5) * rotor
}
Spline <- function(tcb){
  # Kochanek-Bartels quaternions spline
  rotors <- KochanekBartels(
    keyRotors, n_intertimes = 15L, endcondition = "closed", tcb=tcb
  )
  # construction of the interpolating points on the sphere
  points <- matrix(nrow = 0L, ncol = 3L)
  keyPoint1 <- rbind(keyPoints[1L, ])
  for(i in seq_along(rotors)){
    points <- rbind(points, rotate(keyPoint1, rotors[i]))
  }
  points
}

shinyServer(
  function(input, output, session){
    SPLINE <- reactiveVal(Spline(c(-1, 5, 0)))
    
    observeEvent(list(input[["numt"]], input[["numc"]], input[["numb"]]),{
      SPLINE(Spline(c(input[["numt"]], input[["numc"]], input[["numb"]])))
    }, ignoreInit = TRUE, priority = 1)
    
    output[["sphere"]] <- renderGlobe({
      # visualize the result with the 'rgl' package
      # trigger()
      #par3d(userMatrix = input$par3d$userMatrix)
      points <- SPLINE()
      # trigger(trigger()+1)
      allpoints <- rbind(points, keyPoints)
      allpoints <- cart2latlong(allpoints[,1], allpoints[, 2], allpoints[, 3])
      lat <- allpoints[, 1]
      long <- allpoints[, 2]
      
      # spheres3d(points, radius = 0.2, color = "midnightblue")
      # spheres3d(keyPoints, radius = 0.25, color = "red")
      globejs(
        lat = lat,
        long = long,
        value = 5,
        color = c(rep("midnightblue", length(points)), rep("red", length(keyPoints))),
        pointsize = 10
      )
      
    })
    
  }
)



