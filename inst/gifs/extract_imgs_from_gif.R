# setwd("C:/SL/MyPackages/qsplines/inst/gifs")
# 
# gif <- "KochanekBartels.gif"
# file.copy(gif, to = "gif.gif")
# command <- "magick convert gif.gif -coalesce gif-%05d.png"
# system(command)
# 
# command <- "gifski --frames=gif-0*.png --fps=10 -o KochanekBartels.gif"
# system(command)
# 
# appendGIFs <- function(gif1, gif2, gifout="result.gif", vertically=FALSE, 
#                        fps = 10, convert = "magick convert"){
#   currentdir <- getwd()
#   on.exit(setwd(currentdir))
#   tmpdir <- tempdir()
#   invisible(file.remove(list.files(tmpdir, full.names = TRUE, pattern = "*.gif$")))
#   file.copy(gif1, to = file.path(tmpdir, "gif1.gif"))
#   file.copy(gif2, to = file.path(tmpdir, "gif2.gif"))
#   setwd(tmpdir)
#   command <- sprintf("%s gif1.gif -coalesce gif1-%%05d.png", convert)
#   system(command)
#   command <- sprintf("%s gif2.gif -coalesce gif2-%%05d.png", convert)
#   system(command)
#   nframes <- length(list.files(tmpdir, pattern = "^gif1-.*png$"))
#   for(i in 1:nframes){
#     command <- sprintf("%s gif1-%05d.png gif2-%05d.png %sappend gif-%05d.png", 
#                        convert, i-1, i-1, ifelse(vertically, "-", "+"), i)
#     system(command)
#   }
#   command <- sprintf("gifski --frames=gif-*.png --fps=%d -s 512x256 -o result.gif", fps)
#   system(command)
#   file.copy("result.gif", file.path(currentdir, gifout), overwrite=TRUE)
# }
# 
# appendGIFs("BarryGoldman.gif", "KochanekBartels.gif")
