#plot FANC neurons
library(fancr)
#library(RColorBrewer)
library(googlesheets4)
#library(webshot2)
library(malevnc)

dn_lut = c(DNg02 = "#3A9E79", DNp54 = "#D15E13", DNp08 = "#7670B1", DNg32 = "#DE2A88", DNp38 = "#6DA62D", "26730" = "#E1AB24", DNg29 = "#A0711C", "17526" = "#666666", "DN (summed)" = "#CBCBCB")

#just plot one of each side
fanc_ahns =c(msahn_r = "648518346489573207", mtahn_l = "648518346488561230")
manc_ahns = c(msahn_r = 12536, mtahn_l = 42819)
fanc_color_lut = c("limegreen", "magenta4")
manc_color_lut = c("steelblue4", "yellow4")

#fanc_set_token()
#should read from a spreadsheet instead of hardcoding

neuron_meshes = read_fanc_meshes(fanc_ahns)

FANC.hxsurf = as.hxsurf(FANC.surf)

#open3d(antialias=4)
open3d()
userMatrix = matrix(c(0.05955527,0.98893923,-0.13578475,0,-0.99598593,0.04977947,-0.07428607,0,-0.06670625,0.13966437,0.98794132,0,0.00000000,0.00000000,0.00000000,1), byrow=T, nrow = 4, ncol = 4)
par3d(windowRect=c(0,0,1920,600),zoom=0.5,userMatrix=userMatrix)
plot3d(FANC.hxsurf,alpha=0.1, color="grey")
par3d(ignoreExtent = TRUE)
for (i in 1:length(neuron_meshes)) wire3d(neuron_meshes[[i]], col = fanc_color_lut[i], lit = FALSE)
snapshot3d(paste0("ahns_fanc_", Sys.Date(), ".png"),fmt="png",webshot = FALSE)
#rgl.postscript("msahn_dns_fanc.svg",fmt="svg")
rgl.close()

dorsal = matrix(c(0.97966915,0.1877270,-0.07075554,0,-0.01514234,0.4208734,0.90699244,0,0.20004621,-0.8874811,0.41515934,0,0.00000000,0.0000000,0.00000000,1), nrow = 4, ncol = 4, byrow = TRUE)
ventral = matrix(c(-0.984182000,-0.1682267,0.05554246,0,-0.002579868,0.3270945,0.94498819,0,-0.177139848,0.9298971,-0.32235441,0,0.000000000,0.0000000,0.00000000,1), nrow = 4, ncol = 4, byrow = TRUE)

manc_neuron_meshes = read_manc_meshes(manc_ahns)
for (i in 1:length(manc_neuron_meshes)) manc_neuron_meshes[[i]]$vb[1:3,] = manc_neuron_meshes[[i]]$vb[1:3,]/1000
open3d()
par3d(userMatrix = ventral, zoom = 0.613913536071777, windowRect = c(1792L, 45L, 2445L, 1079L), FOV = 0)
plot3d(MANC.surf,alpha=0.1, color="grey")
par3d(ignoreExtent = TRUE)
for (i in 1:length(manc_neuron_meshes)) wire3d(manc_neuron_meshes[[i]], col = manc_color_lut[i], lit = FALSE)
snapshot3d(paste0("ahns_manc_", Sys.Date(), ".png"),fmt="png",webshot = FALSE)
#rgl.postscript("mtahn_dns_fanc.svg",fmt="svg")
rgl.close()

manc_meshes_xform = transform_fanc2manc(read_manc_meshes(manc_ahns), inverse = T)
open3d()
userMatrix = matrix(c(0.05955527,0.98893923,-0.13578475,0,-0.99598593,0.04977947,-0.07428607,0,-0.06670625,0.13966437,0.98794132,0,0.00000000,0.00000000,0.00000000,1), byrow=T, nrow = 4, ncol = 4)
par3d(windowRect=c(0,0,1920,600),zoom=0.5,userMatrix=userMatrix)
plot3d(FANC.hxsurf,alpha=0.1, color="grey")
par3d(ignoreExtent = TRUE)
wire3d(neuron_meshes[[1]], col = fanc_color_lut[1], lit = FALSE)
wire3d(manc_meshes_xform[[1]], col = manc_color_lut[1], lit = FALSE)
snapshot3d(paste0("msahns_manc_to_fanc_xform_", Sys.Date(), ".png"),fmt="png",webshot = FALSE)
#rgl.postscript("msahn_dns_fanc.svg",fmt="svg")
rgl.close()

open3d()
userMatrix = matrix(c(0.05955527,0.98893923,-0.13578475,0,-0.99598593,0.04977947,-0.07428607,0,-0.06670625,0.13966437,0.98794132,0,0.00000000,0.00000000,0.00000000,1), byrow=T, nrow = 4, ncol = 4)
par3d(windowRect=c(0,0,1920,600),zoom=0.5,userMatrix=userMatrix)
plot3d(FANC.hxsurf,alpha=0.1, color="grey")
par3d(ignoreExtent = TRUE)
wire3d(neuron_meshes[[2]], col = fanc_color_lut[2], lit = FALSE)
wire3d(manc_meshes_xform[[2]], col = manc_color_lut[2], lit = FALSE)
snapshot3d(paste0("mtahns_manc_to_fanc_xform_", Sys.Date(), ".png"),fmt="png",webshot = FALSE)
rgl.close()