###~~~~~~~~~~~~~loading rootfile at beginning of session and read in correct file

library(geomorph) # version 3.0.6

library(abind) # version 1.4-3

library(caper)

root.file.path=getwd() # or relplace with root file e.g. "C:/Users/Kaylene/Documents/Git Hub/"

# Load in one of the following datasets:

  #full crania dataset
  load(paste(root.file.path, "/Processed_data/crania_with_fossils.rda", sep="")) 

  #full crania dataset without fossils
  load(paste(root.file.path, "/Processed_data/crania_without_fossils.rda", sep="")) 

  #crania dataset without curve of the molar row
  load(paste(root.file.path, "/Processed_data/crania_without_molar_curve.rda", sep="")) 

  #crania dataset without curve of molar row OR fossils
  load(paste(root.file.path, "/Processed_data/crania_without_molars_and_fossils.rda", sep="")) 

  #Dentary dataset
  load(paste(root.file.path, "/Processed_data/dentary_with_fossils.rda", sep="")) 

  #Dentary dataset without fossils
  load(paste(root.file.path, "/Processed_data/dentary_without_fossils.rda", sep=""))

####~~~~~~~~~~~~PCA PLOTTING for Butler et al. - Kangaroo 3D geometric Morpometrics~~~~~~~~~~~~~~~~~###

  
###Setting up
  
  cols_diet <- c("blue", "black", "dark green", "red", "yellow")
  names(cols_diet) <- levels(classifiers$Diet)
  cols_diet <- cols_diet[match(classifiers$Diet, names(cols_diet))] 
  
  cols_loco <- c("black", "blue", "red", "dark green")
  names(cols_loco) <- levels(classifiers$Movement)
  cols_loco <- cols_loco[match(classifiers$Movement, names(cols_loco))] 
  
  
  
  
  point_diet <-c(21,23,24,22,25)
  point_loco <-c(23,21,24,22)
  names(point_diet)=levels(classifiers$Diet)
  names(point_loco)=levels(classifiers$Movement)
  points_diet <- point_diet[match(classifiers$Diet, names(point_diet))] 
  points_loco <-point_loco[match(classifiers$Movement , names(point_loco))]
  
  
  
### PCA plot for diet with centroid info

plot.new()
  

plot(PCA$pc.scores[,1], PCA$pc.scores[,2], asp=T, 
     xlab= paste("PC 1 ","(", round(100*Roo.PCA$pc.summary$importance[2,1],1),"%)",sep=""), 
     ylab= paste("PC 2 ","(", round(100*Roo.PCA$pc.summary$importance[2,2],1),"%)",sep=""), pch = points_diet,cex = (CSize/10), bg = cols_diet) # plots PCA with specimens scales to CScize
text(PCA$pc.scores[,1], PCA$pc.scores[,2], labels = classifiers$Number, cex=0.5, pos=1)
#par(mar=c(1, 5, 1, 1))  # sets the margins


legend("bottomleft", legend=as.vector(levels(classifiers$Diet)), 
                  pt.bg= col.gp[levels(classifiers$Diet)],
                  pch=point_diet, cex=0.9, pt.cex=1.3, ncol = 2, bty="n") # plots legend below

plot.new()

plot(PCA$pc.scores[,1], PCA$pc.scores[,3], asp=T, 
     xlab= paste("PC 1 ","(", round(100*Roo.PCA$pc.summary$importance[2,1],1),"%)",sep=""), 
     ylab= paste("PC 3 ","(", round(100*Roo.PCA$pc.summary$importance[2,3],1),"%)",sep=""), pch = point.gp,cex = (CSize/10), bg = cols_diet) # plots PCA with specimens scales to CScize
text(PCA$pc.scores[,1], PCA$pc.scores[,3], labels = classifiers$Number, cex=0.5, pos=1)




####To get PCA importance
Roo.PCA$pc.summary$importance[,1:4]

##############PCA Plot - Locomotion################
#cols_loco <- c("black", "blue", "red", "dark green") #remove black colour if fossils removed
#names(cols)=levels(classifiers$Movement) #Movemnet = locmotor mode
#col_loco <- cols[match(classifiers$Movement, names(cols))]

#point2 <-c(23,21,24,22)
#names(point2)=levels(classifiers$Movement)
#point2.gp <- point2[match(classifiers$Movement, names(point2))]

#Create PCA lables
#pc1lab <-paste("Principal Component 1 ","(",
               #round(100*Roo.PCA$pc.summary$importance[2,1],1),"%)",sep="")
#pc2lab <-paste("Principal Component 2 ","(",
              # round(100*Roo.PCA$pc.summary$importance[2,2],1),"%)",sep="")
#pc3lab <-paste("Principal Component 3 ","(",
               #round(100*Roo.PCA$pc.summary$importance[2,3],1),"%)",sep="")



### PCA plot for diet with centroid info

plot.new()


plot(PCA$pc.scores[,1], PCA$pc.scores[,2], asp=T, 
     xlab= paste("PC 1 ","(", round(100*Roo.PCA$pc.summary$importance[2,1],1),"%)",sep=""), 
     ylab= paste("PC 2 ","(", round(100*Roo.PCA$pc.summary$importance[2,2],1),"%)",sep=""), pch = points_loco,cex = (CSize/10), bg = cols_loco) # plots PCA with specimens scales to CScize
text(PCA$pc.scores[,1], PCA$pc.scores[,2], labels = classifiers$Number, cex=0.5, pos=1)
#par(mar=c(1, 5, 1, 1))  # sets the margins


legend(-0.17, -0.035, legend=as.vector(levels(classifiers$Movement)), 
       pt.bg= cols_loco[levels(classifiers$Movement)],
       pch=point_loco, cex=0.9, pt.cex=1.3, ncol = 2, bty="n") # plots legend below

plot.new()

plot(PCA$pc.scores[,1], PCA$pc.scores[,3], asp=T, 
     xlab= paste("PC 1 ","(", round(100*Roo.PCA$pc.summary$importance[2,1],1),"%)",sep=""), 
     ylab= paste("PC 3 ","(", round(100*Roo.PCA$pc.summary$importance[2,3],1),"%)",sep=""), pch = point.gp,cex = (CSize/10), bg = cols_loco) # plots PCA with specimens scales to CScize
text(PCA$pc.scores[,1], PCA$pc.scores[,3], labels = classifiers$Number, cex=0.5, pos=1)


####~~~~~~~~~~~~~~~~~~~~~~~Code to import 3D mesh~~~~~~~~~~~~~~~~~~~~~~~~~~###  

############Code for analysis of crania with curve of molar row#################

#next define file pathway depending on the dataset

#Full crania dataset pathways
    model_file_path="Processed_data/Crania_Models/"
    figure_file_path="Processed_data/Crania_Figures/"

#Curve of molar row pathways
    #model.file.path="Processed_data/No_molar_curve_Models/"
    #figure.file.path="Processed_data/No_molar_curve_Figures/"

#Dentary pathways
    #model.file.path="Processed_data/Dentary_Models/"
    #figure.file.path="Processed_data/Dentary_Figures/"

####################Import ref mesh#################################
    
    #First, warp the skull closest to the mean to the mean shape
ref <- mshape(coords) #reference coordinates are the mean shape
    
# Load in the PLY file and prepare for plotting
mesh <- read.ply("Raw_data/Crania/Macropus_parma-MZRC6451-cranium.ply") #model for mean shape #change mean specimen depending on dataset
mesh$material <- list(color = "pink") ; open3d(); shade3d(mesh) 
mesh.coords <- read.morphologika("Raw_data/Crania/Morphologika files/Macropus_parma-MZRC6451.txt") #change mean specimen depending on dataset
  
  #For Crnaia without molar row 
     #mesh.coords <- read.morphologika("Raw_data/Crania/Macropus_parma_reduced.txt")

refmesh <- warpRefMesh(mesh, mesh.coords[,,1]/1000, ref, centered=FALSE)

#Code for analysis of dentary

      #     mesh <- read.ply("Raw_data/Dentary/Petrogale_herberti.ply")
      #     mesh$material <- "pink" ; open3d(); shade3d(mesh)
      # 
      #     mesh.coords <- read.morphologika("Raw_data/Dentary/Petrogale_herberti_reduced.txt")
      #     refmesh <- warpRefMesh(mesh, mesh.coords[,,1]/1000, ref, centered=FALSE)

open3d(); shade3d(refmesh);writePLY(paste(model.file.path, "/refmesh.ply", sep=""),withColors=T)


####################Create PC renderings and export images #######################################

# Warp the refmesh to the min and max of PC 1 and PC 2 and save mesh to WD/Models

mPC1min <- plotRefToTarget(ref, PCA$pc.shapes$PC1min, method="surface", 
                           mesh = refmesh, mag=1); writePLY(paste(model.file.path, "/PC1min.ply", sep=""),withColors=T)
mPC1max <- plotRefToTarget(ref, PCA$pc.shapes$PC1max, method="surface", 
                           mesh = refmesh, mag=1); writePLY(paste(model.file.path, "/PC1max.ply", sep=""), withColors=T)
mPC2min <- plotRefToTarget(ref, PCA$pc.shapes$PC2min, method="surface", 
                           mesh = refmesh, mag=1); writePLY(paste(model.file.path, "/PC2min.ply", sep=""),withColors=T)
mPC2max <- plotRefToTarget(ref, PCA$pc.shapes$PC2max, method="surface", 
                           mesh = refmesh, mag=1); writePLY(paste(model.file.path, "/PC2max.ply", sep=""), withColors=T)


## For PCs 3 and 4, need a new PCA:
PCA2 <- plotTangentSpace(coords, axis1 =3, axis2 =4)
PCA$pc.shapes <- append(PCA$pc.shapes, PCA2$pc.shapes) # here I am putting the PC3 and PC4 into the orginal PCA$pc.shapes list, just to make it simpler for calling later
remove(PCA2)

mPC3min <- plotRefToTarget(ref, PCA$pc.shapes$PC3min, method="surface", 
                           mesh = refmesh, mag=1); writePLY(paste(model.file.path, "/PC3min.ply", sep=""), withColors=T)
mPC3max <- plotRefToTarget(ref, PCA$pc.shapes$PC3max, method="surface", 
                           mesh = refmesh, mag=1); writePLY(paste(model.file.path, "/PC3max.ply", sep=""), withColors=T)



save.image(paste(root.file.path, "/Processed_data/crania_3D_rendering.rda", sep="")) #change depending on dataset
#save.image(paste(root.file.path, "/Processed_data/crania_less_molar_3D_rendering.rda", sep="")) #change depending on dataset
#save.image(paste(root.file.path, "/Processed_data/dentary_3D_rendering.rda", sep="")) #change depending on dataset

### Create TPS warps for figures

open3d(); view3d(fov=0);shade3d(mPC1min) # open a mesh on rgl
usrMat.dorsal <- par3d()$userMatrix # set by hand, adjust specimen into positon
usrMat.lateral <- par3d()$userMatrix # set by hand, adjust specimen into positon
usrMat.posterior <- par3d()$userMatrix # set by hand, adjust specimen into positon
usrMat.ventral <- par3d()$userMatrix # set by hand, adjust specimen into positon
windowRect <- c(0,0,600,600) # sets size of window to 600 by 600 pixels

## Load the specimen in the required positions 
open3d(FOV=0, userMatrix = usrMat.dorsal, windowRect=windowRect) # for a dorsal view
open3d(FOV=0, userMatrix = usrMat.lateral, windowRect=windowRect) # for a lateral view
open3d(FOV=0, userMatrix = usrMat.posterior, windowRect=windowRect) # for a posterior view
open3d(FOV=0, userMatrix = usrMat.ventral, windowRect=windowRect) # for a lateral view

## Plotting meshes and saving a screenshot of each; run all three lines together
# Plots by hand see alterante code below
open3d(FOV=0, userMatrix = usrMat.dorsal, windowRect=windowRect)
shade3d(refmesh, alpha=1, lit = T, shininess = 128.0, specular=49)
rgl.snapshot(paste(figure.file.path, "/Mesh_refmesh_dorsal.png", sep=""), fmt="png", top=TRUE)

open3d(FOV=0, userMatrix = usrMat.lateral, windowRect=windowRect)
shade3d(mPC1min, alpha=1, lit = T, shininess = 128.0, specular=49)
rgl.snapshot(paste(figure.file.path, "/Mesh_refmesh_lateral.png", sep=""), fmt="png", top=TRUE)

open3d(FOV=0, userMatrix = usrMat.posterior, windowRect=windowRect)
shade3d(mPC1min, alpha=1, lit = T, shininess = 128.0, specular=49)
rgl.snapshot(paste(figure.file.path,"/Mesh_refmesh_posterior.png", sep=""), fmt="png", top=TRUE)

open3d(FOV=0, userMatrix = usrMat.ventral, windowRect=windowRect)
shade3d(mPC1min, alpha=1, lit = T, shininess = 128.0, specular=49)
rgl.snapshot(paste(figure.file.path, "/Mesh_refmesh_ventral.png", sep=""), fmt="png", top=TRUE)

## This loop for the above

views <- c("dorsal", "lateral", "posterior", "ventral")
PCs <- c("PC1min", "PC1max", "PC2min", "PC2max", "PC3min", "PC3max")
for(i in PCs){
  for(j in views){
    open3d(FOV=0, userMatrix=as.name(paste("usrMat", j, sep=".")), windowRect=windowRect)
    shade3d(get(paste("m",i,sep="")), alpha=1, lit = T, shininess = 128.0, specular=49) # plot mesh
    options(warn = 0) # supresses warnings in case they appear in console
    filename <- paste(figure.file.path,"Mesh_",i,"_",j,".png", sep="")
    rgl.snapshot(filename, fmt="png", top=TRUE)
    rgl.close() # closes the window!
  }} 

#This to export lolipop diagrams by hand
plotRefToTarget(PCA$pc.shapes$PC1min, PCA$pc.shapes$PC1max, mag=1, method="vector", gridPars=gridPar(pt.bg = "red", pt.size = 1), label = F)
plotRefToTarget(PCA$pc.shapes$PC2min, PCA$pc.shapes$PC2max, mag=1, method="vector", gridPars=gridPar(pt.bg = "red", pt.size = 1), label = F)
plotRefToTarget(PCA$pc.shapes$PC3min, PCA$pc.shapes$PC3max, mag=1, method="vector", gridPars=gridPar(pt.bg = "red", pt.size = 1), label = F)


views <- c("dorsal", "lateral", "posterior", "ventral")
PCs <- c("PC1min", "PC1max", "PC2min", "PC2max", "PC3min", "PC3max")
####Lollipop loop
attach(PCA$pc.shapes) # so we can call each PC shape from this list
for(i in PCs){
  for(j in views){
    tmp <- usrMat.lateral # set which view to use
    open3d(FOV=0, userMatrix=as.name(paste("usrMat", j, sep=".")), windowRect=windowRect)
    plotRefToTarget(ref, PC1min, method="vector", mag=1, axes=F)
    #shade3d(ref, alpha=0.5)
    options(warn = 0) # supresses warnings in case they appear in console
    filename <- paste(figure.file.path,"Lollipop_",i,"_",j,".png", sep="")
    rgl.snapshot(filename, fmt="png", top=TRUE)
    rgl.close() # closes the window!
  }}
detach(PCA$pc.shapes)  # clean up environment

