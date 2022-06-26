library(onion)

getid <- function(scene){
  for(obj in scene$objects){
    if(obj$type != "spheres"){
      next
    }
    if(nrow(obj$vertices) == 120){
      print("@@@@@@@@@@@@@@@@@120@@@@@@@@@@@@@@@")
      return(obj$id)
    }
  }
}


# helper function: spherical to Cartesian coordinates
sph2cart <- function(rho, theta, phi){
  return(c(
    rho * cos(theta) * sin(phi),
    rho * sin(theta) * sin(phi),
    rho * cos(phi)
  ))
}

try(close3d())
try(clear3d())
#open3d(c(0,0,500,500))
open3d(useNULL = TRUE)
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
# SplineStatic <- function(){
#   # Kochanek-Bartels quaternions spline
#   rotors <- KochanekBartels(
#     keyRotors, n_intertimes = 15L, endcondition = "closed", 
#     tcb = c(0, 0, 0)
#   )
#   # construction of the interpolating points on the sphere
#   points <- matrix(nrow = 0L, ncol = 3L)
#   keyPoint1 <- rbind(keyPoints[1L, ])
#   for(i in seq_along(rotors)){
#     points <- rbind(points, rotate(keyPoint1, rotors[i]))
#   }
#   list(points = points, keyPoints=keyPoints)
# }
# points <- SplineStatic()[["points"]]
# keyPoints <- SplineStatic()[["keyPoints"]]

obj <- drawSphere(radius = 5, color = "lightgreen")
scene <- scene3d()
#close3d()
#root <- scene$rootSubscene
#rgl.close()

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
    options(rgl.useNULL = TRUE)
    #plot3d(scene)
    #root <- currentSubscene3d()
    #dev <- cur3d()
    save <- options(rgl.inShiny = TRUE)
    on.exit(options(save))

    trigger <- reactiveVal(0)
    SPLINE <- reactiveVal(Spline(c(-1, 5, 0)))
    OBJ <- NULL
    observeEvent(list(input[["numt"]], input[["numc"]], input[["numb"]]), {
      print("DDDDEEEELLLLLLLLLLLLLLL")
      delFromSubscene3d(OBJ)
    }, ignoreInit = TRUE, priority = 10)
    
    observeEvent(list(input[["numt"]], input[["numc"]], input[["numb"]]),{
      SPLINE(Spline(c(input[["numt"]], input[["numc"]], input[["numb"]])))
    }, ignoreInit = TRUE, priority = 1)
    
    id <- w <- NULL
    subid=NULL

    output[["sphere"]] <- renderRglwidget({
      # visualize the result with the 'rgl' package
      # trigger()
      #par3d(userMatrix = input$par3d$userMatrix)
      points <- SPLINE()
      # trigger(trigger()+1)
      OBJ <<- spheres3d(points, radius = 0.2, color = "midnightblue")
      scene <<- scene3d()
      # delFromSubscene3d(ids)#getid(scene))
      spheres3d(keyPoints, radius = 0.25, color = "red")
      rglwidget()

    })
    
    observeEvent(input$test, {
      #root <- currentSubscene3d()
      #useSubscene3d(root)
      #useSubscene3d(root)
      #set3d(dev)
      delFromSubscene3d(getid(scene))
      trigger(trigger()+1)
      # M <- diag(4)
      # M[1:2,1:2] = rbind(c(0,1),c(1,0))
      # shinyGetPar3d("userMatrix", session=session)
    })
    
#     output$thecontroller <-
#       renderPlaywidget({
#         M <- diag(4)#r3dDefaults$userMatrix
#         print("MMMMMMMMMMMMMMMMMmm")
#         print(M)
#         fn <- par3dinterp(times = (0:2)*0.75, userMatrix = list(M,
#                                                                 rotate3d(M, pi/2, 1, 0, 0),
#                                                                 rotate3d(M, pi/2, 0, 1, 0)),
#                           scale = c(0.5, 1, 2))
#         control <- par3dinterpControl(fn, 0, 3, steps = 15)
#         
# #        rglwidget %>%
#         print("XXXXXXXXXX")
#         drawSphere()
#         rglwidget() %>% playwidget(controls=control, components = c("Play", "Reset"))#propertyControl(value=M, tags="spline"),
#       })
    
      observe({
        print(names(input))
      })
    
    
  })