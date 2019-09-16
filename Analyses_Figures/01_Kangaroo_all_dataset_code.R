####~~~~~~~~~~~~~~~~~~Code for Butler et al. - Kangaroo 3D geometric Morpometrics - crania with all landmarks~~~~~~~~~~~~~~~~~###

# Load first libraries
library(geomorph) # version 3.0.6

library(abind) # version 1.4-3

library(caper)

#NOTE: The following code is designed to run all at once to export read in files for each dataset

##Define root directory - only when first opening R
root.file.path=getwd() # or relplace with root file pathway e.g. "C:/Users/Kaylene/Documents/Git Hub/"

############################################Crania with fossils ##################################################


# set Working directory 
setwd(paste(root.file.path, "/Raw_data/Crania", sep=""))

#### loading 3d coordinate files

###FIRST REPLICATE - Using the 48 landmark dataset
####the following reads in all morphologica files### 
filelist <- list.files(path= "Morphologika files", pattern = "*.txt")
names <- gsub (".txt", "", filelist) # extracts names of specimens from the file name
for (i in 1:length(filelist)){                 #code to cut spemen numbers from names
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Morphologika files/", filelist, sep="") # rename with path
coords <- NULL # make empty object that will be filled with 3D array of coordinate data
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  #If loop ONLY for specimens that are not scaled
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
Roo <- arrayspecs(coords, k, 3)
dimnames(Roo)[[3]] <- names
remove(i, filelist, names, k, coords, temp) #remove unnecessary files

##SECOND REPLICATE
filelist <- list.files(path= "Rep2", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Rep2/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
rep.2 <- arrayspecs(coords, k, 3)
dimnames(rep.2)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 


#########Estimate missing data################
for (i in 1:length(Roo)){
  if(Roo[i]==9999){
    Roo[i]<-NA
  }
} #Replaces 9999 with NA

Roo<-estimate.missing(Roo,method="Reg")

####Rep,2 Missing
for (i in 1:length(rep.2)){
  if(rep.2[i]==9999){
    rep.2[i]<-NA
  }
}

rep.2<-estimate.missing(rep.2,method="Reg")

###Check if file lists are similar
filelist1=list.files(path= "Morphologika files", pattern = "*.txt")
filelist2=list.files(path= "Rep2", pattern = "*.txt")
match(filelist1, filelist2)

Roo<- abind(Roo,rep.2) # abind concatenates (stacks) two arrays
remove(rep.2) # remove unneeded file

dimnames(Roo)
attributes(Roo)  # shows the specimens loaded with file names
class(Roo)  #check that this is an array
Roo

plot3d(Roo[,,14], asp=F)#plot the 5th species of this array
text3d(Roo[,,2],text=c(1:71))#Gives the landmarks numbers 1 for 1st coord row etc...

### what lms are semilandmarks on curves (curveslide), patches(surfslide) 

curveslide=as.matrix(read.csv("curveslidenew.csv", header=T))

####Procrustes superimposition

Roo.GPA <- gpagen(Roo,curves=curveslide)#if using surface/curve semilandmarks, add curves=curveslide, surfaces=surfslide 

#####Now run PCA
Roo.PCA=plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]])

#for checking, identify the numbers of the specimens that are OK and that are outliers, then run the following:
plotRefToTarget(Roo.GPA$coords[,,7], Roo.GPA$coords[,,49],  method="vector", label = T)
dimnames(Roo)

### Next, average REPLICATES

individualsBase=dimnames(Roo.GPA$coords)[[3]] # make a vector of the names of the specimens (with two replicates, this will mean each name given twice)
means <- (aggregate(two.d.array(Roo.GPA$coords) ~ individualsBase, FUN=mean))[,-1] #Here the 3D array converts into a 2D array and then is averaged by specimen
rownames(means) <- unique(dimnames(Roo.GPA$coords)[[3]]) 
Roo.GPA_Aver <- arrayspecs(means, p=ncol(means)/3, k=3)# p equals number of landmarks; arrayspecs converts landmark matrix back into 3D array
Roo.PCA_Aver=plotTangentSpace(Roo.GPA_Aver, label=rownames(means) )

findMeanSpec(Roo.GPA_Aver) #specimen closest to the mean shape 

plotRefToTarget(Roo.GPA_Aver[,,14], Roo.GPA_Aver[,,15],  method="vector", label = T)

#Aggregating also for Csize
CSize <- (aggregate(Roo.GPA$Csize ~ individualsBase, FUN=mean))[-1]
rownames(CSize) <- unique(dimnames(Roo.GPA$coords)[[3]])

landpairs=as.matrix(read.csv("Landpairs.csv", header = F))#make a landmark pair file to adjust for symmetry below
Roo.GPASymm<-bilat.symmetry(Roo.GPA_Aver, ind=dimnames(Roo.GPA_Aver)[[3]], side=NULL, replicate=NULL, object.sym=TRUE, land.pairs =landpairs)

######################################################################################
# Comparing PCAs
RooPCASymm=plotTangentSpace(Roo.GPASymm$symm.shape,label=dimnames(Roo.GPASymm$symm.shape)[[3]], warpgrids = T)

# for comparison, without symmetry: 
plotTangentSpace(Roo.GPA_Aver, label=dimnames(Roo.GPA_Aver)[[3]], warpgrids = T)
plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]], warpgrids = T)

findMeanSpec(Roo.GPASymm$symm.shape)#specimen closest to the mean shape for finding the correct mesh for warping

######################################################################################

###########Looking at differences between minimum and maximum PCAS
ref <- mshape(Roo.GPA$coords)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC1min)
plotRefToTarget(Roo[,,15], Roo[,,59], method="vector", label = T)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC2min)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC3min)

####~~~~~~~~~~~~Define coords and Csize~~~~~~~~~~~~~~~~~###
###Define coords
ibase=unique(individualsBase)
SpeciesGroupsSingle <- ibase #redunant line now

SpeciesGroupsSingle<-as.factor(SpeciesGroupsSingle)
is.factor(SpeciesGroupsSingle) # TRUE; the function read.table automatically makes text factors
coords <- (aggregate(two.d.array(Roo.GPASymm$symm.shape) ~ SpeciesGroupsSingle, FUN=mean))[,-1]
rownames(coords) <- unique(SpeciesGroupsSingle)
coords <- arrayspecs(coords,p=dim(Roo.GPASymm$symm.shape)[1],k=3) # this is the 3d array for plotTangentSpace

CSize=as.matrix(CSize)
CSize<- (aggregate(CSize ~ SpeciesGroupsSingle, FUN=mean))[,-1] #new centroid size vector
names(CSize)<-unique(SpeciesGroupsSingle) # object is a vector so use names() not rownames()

coords 
CSize #Associated centroid sizes

###Define 'PCA'

PCA <- plotTangentSpace(coords, warpgrids=F) 

summary(PCA)

####read in classifiers
classifiers=read.table("classifiers.txt", header=T)
classifiers$Fossil <- c(rep(1,length(rownames(classifiers))));classifiers$Fossil [ c(which(classifiers$Diet=="fossil"))] <-2
rownames(classifiers) <- classifiers$Species

trees <- read.nexus("Roo.con.tre") 

save.image(paste(root.file.path, "/Processed_data/crania_with_fossils.rda", sep=""))

setwd(root.file.path)



#########################################Crania without fossils#######################################################Clear environment
rm(list = ls(all.names = TRUE))

root.file.path=getwd() 
setwd(paste(root.file.path, "/Raw_data/Without fossils/Crania", sep=""))

filelist <- list.files(path= "Morphologika files", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Morphologika files/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
Roo <- arrayspecs(coords, k, 3)
dimnames(Roo)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

##SECOND REPLICATE
filelist <- list.files(path= "Rep2", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Rep2/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
rep.2 <- arrayspecs(coords, k, 3)
dimnames(rep.2)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

filelist1=list.files(path= "Morphologika files", pattern = "*.txt")
filelist2=list.files(path= "Rep2", pattern = "*.txt")
match(filelist1, filelist2)

Roo<- abind(Roo,rep.2) 
remove(rep.2) 

dimnames(Roo)
attributes(Roo)  
class(Roo) 
Roo

plot3d(Roo[,,14], asp=F)
text3d(Roo[,,2],text=c(1:71))

curveslide=as.matrix(read.csv("curveslidenew.csv", header=T))


Roo.GPA <- gpagen(Roo,curves=curveslide)

Roo.PCA=plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]])

plotRefToTarget(Roo.GPA$coords[,,7], Roo.GPA$coords[,,49], method="vector", label = T)
dimnames(Roo)

individualsBase=dimnames(Roo.GPA$coords)[[3]] 
means <- (aggregate(two.d.array(Roo.GPA$coords) ~ individualsBase, FUN=mean))[,-1] 
rownames(means) <- unique(dimnames(Roo.GPA$coords)[[3]]) 
Roo.GPA_Aver <- arrayspecs(means, p=ncol(means)/3, k=3)
Roo.PCA_Aver=plotTangentSpace(Roo.GPA_Aver, label=rownames(means) )

findMeanSpec(Roo.GPA_Aver)

plotRefToTarget(Roo.GPA_Aver[,,14], Roo.GPA_Aver[,,15], method="vector", label = T)

#same needs doing for Csize:
CSize <- (aggregate(Roo.GPA$Csize ~ individualsBase, FUN=mean))[-1]
rownames(CSize) <- unique(dimnames(Roo.GPA$coords)[[3]])

landpairs=as.matrix(read.csv("Landpairs.csv", header = F))
Roo.GPASymm<-bilat.symmetry(Roo.GPA_Aver, ind=dimnames(Roo.GPA_Aver)[[3]], side=NULL, replicate=NULL, object.sym=TRUE, land.pairs =landpairs)

######################################################################################

RooPCASymm=plotTangentSpace(Roo.GPASymm$symm.shape,label=dimnames(Roo.GPASymm$symm.shape)[[3]], warpgrids = T) 

plotTangentSpace(Roo.GPA_Aver, label=dimnames(Roo.GPA_Aver)[[3]], warpgrids = T)
plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]], warpgrids = T)

findMeanSpec(Roo.GPASymm$symm.shape)

######################################################################################

ref <- mshape(Roo.GPA$coords)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC1min)
plotRefToTarget(Roo[,,15], Roo[,,59], method="vector", label = T)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC2min)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC3min)

####~~~~~~~~~~~~Define coords and Csize~~~~~~~~~~~~~~~~~###
ibase=unique(individualsBase)
SpeciesGroupsSingle <- ibase 

SpeciesGroupsSingle<-as.factor(SpeciesGroupsSingle)
is.factor(SpeciesGroupsSingle)
coords <- (aggregate(two.d.array(Roo.GPASymm$symm.shape) ~ SpeciesGroupsSingle, FUN=mean))[,-1]
rownames(coords) <- unique(SpeciesGroupsSingle)
coords <- arrayspecs(coords,p=dim(Roo.GPASymm$symm.shape)[1],k=3)

CSize=as.matrix(CSize)
CSize<- (aggregate(CSize ~ SpeciesGroupsSingle, FUN=mean))[,-1]
names(CSize)<-unique(SpeciesGroupsSingle)

coords 
CSize 

PCA <- plotTangentSpace(coords, warpgrids=F) 

summary(PCA)

classifiers=read.table("classifiers_no_fossils.txt", header=T)
classifiers$Fossil <- c(rep(1,length(rownames(classifiers))));classifiers$Fossil [ c(which(classifiers$Diet=="fossil"))] <-2
rownames(classifiers) <- classifiers$Species

cols <- c("blue", "dark green", "red", "yellow")
names(cols)=levels(classifiers$Diet)
col.gp <- cols[match(classifiers$Diet, names(cols))]

point <-c(21,23,24,22,25)
names(point)=levels(classifiers$Diet)
point.gp <- point[match(classifiers$Diet, names(point))]

trees <- read.nexus("Roo.con.tre")

save.image(paste(root.file.path, "/Processed_data/crania_without_fossils.rda", sep=""))

setwd(root.file.path)



################################crania wihtout molar row###################################################

rm(list = ls(all.names = TRUE))

root.file.path=getwd() 
setwd(paste(root.file.path, "/Raw_data/Crania", sep=""))
 
filelist <- list.files(path= "Morphologika files", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Morphologika files/", filelist, sep="")
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
Roo <- arrayspecs(coords, k, 3)
dimnames(Roo)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

filelist <- list.files(path= "Rep2", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Rep2/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
rep.2 <- arrayspecs(coords, k, 3)
dimnames(rep.2)[[3]] <- names
remove(i, filelist, names, k, coords, temp)


for (i in 1:length(Roo)){
  if(Roo[i]==9999){
    Roo[i]<-NA
  }
}


Roo<-estimate.missing(Roo,method="Reg")

for (i in 1:length(rep.2)){
  if(rep.2[i]==9999){
    rep.2[i]<-NA
  }
}

rep.2<-estimate.missing(rep.2,method="Reg")


############REMOVE CURVE OF MOLAR ROWS###########

Roo<-Roo[-c(6,7,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47),,] #fixed points no curve slide molars,no p3 points - curveslide 17-24
rep.2<-rep.2[-c(6,7,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47),,]

filelist1=list.files(path= "Morphologika files", pattern = "*.txt")
filelist2=list.files(path= "Rep2", pattern = "*.txt")
match(filelist1, filelist2)

Roo<- abind(Roo,rep.2) 
remove(rep.2)

dimnames(Roo)
attributes(Roo)  
class(Roo)  
Roo

plot3d(Roo[,,3], asp=F) 
text3d(Roo[,,2],text=c(1:71)) 

classifiers=read.table("classifiers.txt", header=T)
rownames(classifiers) <- classifiers$Species


curveslide=as.matrix(read.csv("curveslidenew19_26.csv", header=T))


Roo.GPA <- gpagen(Roo,curves=curveslide) 


Roo.PCA=plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]])

plotRefToTarget(Roo.GPA$coords[,,7], Roo.GPA$coords[,,49], method="vector", label = T)
dimnames(Roo)


individualsBase=dimnames(Roo.GPA$coords)[[3]] 
means <- (aggregate(two.d.array(Roo.GPA$coords) ~ individualsBase, FUN=mean))[,-1] 
rownames(means) <- unique(dimnames(Roo.GPA$coords)[[3]]) 
Roo.GPA_Aver <- arrayspecs(means, p=ncol(means)/3, k=3)
Roo.PCA_Aver=plotTangentSpace(Roo.GPA_Aver, label=rownames(means) )

findMeanSpec(Roo.GPA_Aver) 

plotRefToTarget(Roo.GPA_Aver[,,14], Roo.GPA_Aver[,,15], method="vector", label = T)

CSize <- (aggregate(Roo.GPA$Csize ~ individualsBase, FUN=mean))[-1]
rownames(CSize) <- unique(dimnames(Roo.GPA$coords)[[3]])

landpairs=as.matrix(read.csv("Landpairs_minus_molars.csv", header = F))
Roo.GPASymm<-bilat.symmetry(Roo.GPA_Aver, ind=dimnames(Roo.GPA_Aver)[[3]], side=NULL, replicate=NULL, object.sym=TRUE, land.pairs =landpairs)

######################################################################################

RooPCASymm=plotTangentSpace(Roo.GPASymm$symm.shape,label=dimnames(Roo.GPASymm$symm.shape)[[3]], warpgrids = T)#for some reason, this works and is consistent with the previous PCA. Note the way the label is written out - you have to specify the dimname of Brain.GPASymm$symm.shape (rather than Brain.GPASymm[[3]] and such!! ) 


plotTangentSpace(Roo.GPA_Aver, label=dimnames(Roo.GPA_Aver)[[3]], warpgrids = T)
plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]], warpgrids = T)

findMeanSpec(Roo.GPASymm$symm.shape)

######################################################################################

ref <- mshape(Roo.GPA$coords)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC1min)
plotRefToTarget(Roo[,,15], Roo[,,59], , method="vector", label = T)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC2min)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC3min)

######################Define coords and Csize#####################################################

ibase=unique(individualsBase)
SpeciesGroupsSingle <- ibase 

SpeciesGroupsSingle<-as.factor(SpeciesGroupsSingle)
is.factor(SpeciesGroupsSingle) 
coords <- (aggregate(two.d.array(Roo.GPASymm$symm.shape) ~ SpeciesGroupsSingle, FUN=mean))[,-1]
rownames(coords) <- unique(SpeciesGroupsSingle)
coords <- arrayspecs(coords,p=dim(Roo.GPASymm$symm.shape)[1],k=3) 

CSize=as.matrix(CSize)
CSize<- (aggregate(CSize ~ SpeciesGroupsSingle, FUN=mean))[,-1]
names(CSize)<-unique(SpeciesGroupsSingle) 

coords
CSize

PCA <- plotTangentSpace(coords, warpgrids=F) 
write.csv(PCA$pc.scores, file="PCScores.csv")
summary(PCA)

classifiers=read.table("classifiers.txt", header=T)
classifiers$Fossil <- c(rep(1,length(rownames(classifiers))));classifiers$Fossil [ c(which(classifiers$Diet=="fossil"))] <-2
rownames(classifiers) <- classifiers$Species

cols <- c("blue", "black", "dark green", "red", "yellow")
names(cols)=levels(classifiers$Diet)
col.gp <- cols[match(classifiers$Diet, names(cols))]

point <-c(21,23,24,22,25)
names(point)=levels(classifiers$Diet)
point.gp <- point[match(classifiers$Diet, names(point))]

trees <- read.nexus("Roo.con.tre") 

save.image(paste(root.file.path, "/Processed_data/crania_without_molar_curve.rda", sep=""))

setwd(root.file.path)


###########################################Crania without curve of molar row AND without fossils###########################

rm(list = ls(all.names = TRUE))

root.file.path=getwd() 

setwd(paste(root.file.path, "/Raw_data/Without fossils/Crania", sep=""))

filelist <- list.files(path= "Morphologika files", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Morphologika files/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
Roo <- arrayspecs(coords, k, 3)
dimnames(Roo)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

##SECOND REPLICATE if needed
filelist <- list.files(path= "Rep2", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){               
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Rep2/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
rep.2 <- arrayspecs(coords, k, 3)
dimnames(rep.2)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

filelist <- list.files(path= "Rep2", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Rep2/", filelist, sep="") 
coords <- NULL
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
rep.2 <- arrayspecs(coords, k, 3)
dimnames(rep.2)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

############REMOVE CURVE OF MOLAR ROWS###########

Roo<-Roo[-c(6,7,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47),,]
rep.2<-rep.2[-c(6,7,30,31,32,33,34,35,36,37,40,41,42,43,44,45,46,47),,]

filelist1=list.files(path= "Morphologika files", pattern = "*.txt")
filelist2=list.files(path= "Rep2", pattern = "*.txt")
match(filelist1, filelist2)

Roo<- abind(Roo,rep.2) 
remove(rep.2)

dimnames(Roo)
attributes(Roo)  
class(Roo)  
Roo

plot3d(Roo[,,3], asp=F) 
text3d(Roo[,,2],text=c(1:71)) 

curveslide=as.matrix(read.csv("curveslidenew19_26.csv", header=T))


Roo.GPA <- gpagen(Roo,curves=curveslide) 

Roo.PCA=plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]])

plotRefToTarget(Roo.GPA$coords[,,7], Roo.GPA$coords[,,49], method="vector", label = T)
dimnames(Roo)

individualsBase=dimnames(Roo.GPA$coords)[[3]] 
means <- (aggregate(two.d.array(Roo.GPA$coords) ~ individualsBase, FUN=mean))[,-1] 
rownames(means) <- unique(dimnames(Roo.GPA$coords)[[3]]) 
Roo.GPA_Aver <- arrayspecs(means, p=ncol(means)/3, k=3)
Roo.PCA_Aver=plotTangentSpace(Roo.GPA_Aver, label=rownames(means) )

findMeanSpec(Roo.GPA_Aver) 

plotRefToTarget(Roo.GPA_Aver[,,14], Roo.GPA_Aver[,,15], method="vector", label = T)

CSize <- (aggregate(Roo.GPA$Csize ~ individualsBase, FUN=mean))[-1]
rownames(CSize) <- unique(dimnames(Roo.GPA$coords)[[3]])

landpairs=as.matrix(read.csv("Landpairs_minus_molars.csv", header = F))#make a landmark pair file
Roo.GPASymm<-bilat.symmetry(Roo.GPA_Aver, ind=dimnames(Roo.GPA_Aver)[[3]], side=NULL, replicate=NULL, object.sym=TRUE, land.pairs =landpairs)

######################################################################################

RooPCASymm=plotTangentSpace(Roo.GPASymm$symm.shape,label=dimnames(Roo.GPASymm$symm.shape)[[3]], warpgrids = T)#for some reason, this works and is consistent with the previous PCA. Note the way the label is written out - you have to specify the dimname of Brain.GPASymm$symm.shape (rather than Brain.GPASymm[[3]] and such!! ) 


plotTangentSpace(Roo.GPA_Aver, label=dimnames(Roo.GPA_Aver)[[3]], warpgrids = T)
plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]], warpgrids = T)

findMeanSpec(Roo.GPASymm$symm.shape)

######################################################################################

ref <- mshape(Roo.GPA$coords)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC1min)
plotRefToTarget(Roo[,,15], Roo[,,59], method="vector", label = T)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC2min)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC3min)

ibase=unique(individualsBase)
SpeciesGroupsSingle <- ibase 

SpeciesGroupsSingle<-as.factor(SpeciesGroupsSingle)
is.factor(SpeciesGroupsSingle) 
coords <- (aggregate(two.d.array(Roo.GPASymm$symm.shape) ~ SpeciesGroupsSingle, FUN=mean))[,-1]
rownames(coords) <- unique(SpeciesGroupsSingle)
coords <- arrayspecs(coords,p=dim(Roo.GPASymm$symm.shape)[1],k=3) 

CSize=as.matrix(CSize)
CSize<- (aggregate(CSize ~ SpeciesGroupsSingle, FUN=mean))[,-1]
names(CSize)<-unique(SpeciesGroupsSingle) 

coords
CSize

PCA <- plotTangentSpace(coords, warpgrids=F) 
write.csv(PCA$pc.scores, file="PCScores.csv")
summary(PCA)

classifiers=read.table("classifiers_no_fossils.txt", header=T)
classifiers$Fossil <- c(rep(1,length(rownames(classifiers))));classifiers$Fossil [ c(which(classifiers$Diet=="fossil"))] <-2
rownames(classifiers) <- classifiers$Species

cols <- c("blue", "dark green", "red", "yellow")
names(cols)=levels(classifiers$Diet)
col.gp <- cols[match(classifiers$Diet, names(cols))]

point <-c(21,23,24,22,25)
names(point)=levels(classifiers$Diet)
point.gp <- point[match(classifiers$Diet, names(point))]

trees <- read.nexus("Roo.con.tre")

save.image(paste(root.file.path, "/Processed_data/crania_without_molars_and_fossils.rda", sep=""))

setwd(root.file.path)


#########################################Dentary Code#############################################################

rm(list = ls(all.names = TRUE))

root.file.path=getwd() 

setwd(paste(root.file.path, "/Raw_data/Dentary", sep=""))

filelist <- list.files(path= "Morphologika files", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Morphologika files/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
Roo <- arrayspecs(coords, k, 3)
dimnames(Roo)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

filelist <- list.files(path= "Rep2", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Rep2/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
rep.2 <- arrayspecs(coords, k, 3)
dimnames(rep.2)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 


#########Estimate missing data################
for (i in 1:length(Roo)){
  if(Roo[i]==9999){
    Roo[i]<-NA
  }
}


write.morphologika(Roo, file="Dentary coordinates") 

Roo<-estimate.missing(Roo,method="Reg")

####Rep,2 Missing

for (i in 1:length(rep.2)){
  if(rep.2[i]==9999){
    rep.2[i]<-NA
  }
}

rep.2<-estimate.missing(rep.2,method="Reg")

plot3d(Roo[,,2], asp=F)
text3d(Roo[,,2],text=c(1:71))

Roo<-Roo[c(1,3,4,5,7,10,11,12,15,21,22,25,28,31,34,35,40),,] 
rep.2<-rep.2[c(1,3,4,5,7,10,11,12,15,21,22,25,28,31,34,35,40),,]

filelist1=list.files(path= "Morphologika files", pattern = "*.txt")
filelist2=list.files(path= "Rep2", pattern = "*.txt")
match(filelist1, filelist2)

Roo<- abind(Roo,rep.2) 
remove(rep.2) 

dimnames(Roo)
attributes(Roo) 
class(Roo)  
Roo

plot3d(Roo[,,23], asp=F) 
text3d(Roo[,,23],text=c(1:71))

classifiers=read.table("classifiers.txt", header=T)
rownames(classifiers) <- classifiers$Species

curveslide=as.matrix(read.csv("curveslidenew.csv", header=T))

Roo.GPA <- gpagen(Roo)  

Roo.PCA=plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]])

plotRefToTarget(Roo.GPA$coords[,,1], Roo.GPA$coords[,,2], method="vector", label = T)
dimnames(Roo)

individualsBase=dimnames(Roo.GPA$coords)[[3]] 
means <- (aggregate(two.d.array(Roo.GPA$coords) ~ individualsBase, FUN=mean))[,-1] 
rownames(means) <- unique(dimnames(Roo.GPA$coords)[[3]]) 
Roo.GPA_Aver <- arrayspecs(means, p=ncol(means)/3, k=3)
Roo.PCA_Aver=plotTangentSpace(Roo.GPA_Aver, label=rownames(means) )

findMeanSpec(Roo.GPA_Aver) 

plotRefToTarget(Roo.GPA_Aver[,,14], Roo.GPA_Aver[,,15], method="vector", label = T)


CSize <- (aggregate(Roo.GPA$Csize ~ individualsBase, FUN=mean))[-1]
rownames(CSize) <- unique(dimnames(Roo.GPA$coords)[[3]])#

plotTangentSpace(Roo.GPA_Aver, label=dimnames(Roo.GPA_Aver)[[3]], warpgrids = T)
plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]], warpgrids = T)

findMeanSpec(Roo.GPA$coords) #findMeanSpec(Roo.GPASymm$symm.shape) not used as no symm for a single side of a dentary
findMeanSpec(Roo.GPA_Aver)


######################################################################################

ref <- mshape(Roo.GPA$coords)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC1min)
plotRefToTarget(Roo[,,3], Roo[,,44], method="vector", label = T)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC2min)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC3min)

######################Define coords and Csize#####################################################

ibase=unique(individualsBase)
SpeciesGroupsSingle <- ibase 

SpeciesGroupsSingle<-as.factor(SpeciesGroupsSingle)
is.factor(SpeciesGroupsSingle)s
coords <- (aggregate(two.d.array(Roo.GPA_Aver) ~ SpeciesGroupsSingle, FUN=mean))[,-1]
rownames(coords) <- unique(SpeciesGroupsSingle)
coords <- arrayspecs(coords,p=dim(Roo.GPA_Aver)[1],k=3)

CSize=as.matrix(CSize)
CSize<- (aggregate(CSize ~ SpeciesGroupsSingle, FUN=mean))[,-1]
names(CSize)<-unique(SpeciesGroupsSingle) 

coords 
CSize


PCA <- plotTangentSpace(coords, warpgrids=F) 
write.csv(PCA$pc.scores, file="PCScores.csv")
summary(PCA)

####read in classifiers
classifiers=read.table("classifiers.txt", header=T)
classifiers$Fossil <- c(rep(1,length(rownames(classifiers))));classifiers$Fossil [ c(which(classifiers$Diet=="fossil"))] <-2
rownames(classifiers) <- classifiers$Species

cols <- c("blue", "black", "dark green", "red", "yellow")
names(cols)=levels(classifiers$Diet)
col.gp <- cols[match(classifiers$Diet, names(cols))]

point <-c(21,23,24,22,25)
names(point)=levels(classifiers$Diet)
point.gp <- point[match(classifiers$Diet, names(point))]

trees <- read.nexus("Roo.con.tre") 

save.image(paste(root.file.path, "/Processed_data/dentary_with_fossils.rda", sep=""))

setwd(root.file.path)






#############################################Dentary Code no fossils#################################################
library(geomorph) # version 3.0.6

library(abind) # version 1.4-3

library(caper)

root.file.path=getwd() 

setwd(paste(root.file.path, "/Raw_data/Without fossils/Dentary", sep=""))

filelist <- list.files(path= "Morphologika files", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                 
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Morphologika files/", filelist, sep="") # rename with path
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
Roo <- arrayspecs(coords, k, 3)
dimnames(Roo)[[3]] <- names
remove(i, filelist, names, k, coords, temp) 

filelist <- list.files(path= "Rep2", pattern = "*.txt")
names <- gsub (".txt", "", filelist) 
for (i in 1:length(filelist)){                
  hyphen = "-"
  nameI <- names[i]
  replicatename = grepl(hyphen, nameI, fixed=FALSE)
  if(replicatename=="TRUE"){
    nameI <- gsub("\\-.*", "", nameI)
  }
  names[i] <- nameI
}
filelist <- paste("Rep2/", filelist, sep="") 
coords <- NULL 
for (i in 1:length(filelist)){
  temp <- read.morphologika(filelist[i])
  if(grepl("-x1000", filelist[i], fixed=FALSE)==TRUE){  
    temp <- temp/1000
  }
  k <- dim(temp)[1]
  coords <- rbind(coords, two.d.array(temp)) }
rep.2 <- arrayspecs(coords, k, 3)
dimnames(rep.2)[[3]] <- names
remove(i, filelist, names, k, coords, temp)


#########Estimate missing data################
for (i in 1:length(Roo)){
  if(Roo[i]==9999){
    Roo[i]<-NA
  }
}

Roo<-estimate.missing(Roo,method="Reg")

for (i in 1:length(rep.2)){
  if(rep.2[i]==9999){
    rep.2[i]<-NA
  }
}

rep.2<-estimate.missing(rep.2,method="Reg")

plot3d(Roo[,,2], asp=F)
text3d(Roo[,,2],text=c(1:71))

#BEST

Roo<-Roo[c(1,3,4,5,7,10,11,12,15,21,22,25,28,31,34,35,40),,] #landmarks included in analysis
rep.2<-rep.2[c(1,3,4,5,7,10,11,12,15,21,22,25,28,31,34,35,40),,]

filelist1=list.files(path= "Morphologika files", pattern = "*.txt")
filelist2=list.files(path= "Rep2", pattern = "*.txt")
match(filelist1, filelist2)

Roo<- abind(Roo,rep.2) 
remove(rep.2) 

dimnames(Roo)
attributes(Roo) 
class(Roo) 
Roo

plot3d(Roo[,,23], asp=F) 
text3d(Roo[,,23],text=c(1:71))


curveslide=as.matrix(read.csv("curveslidenew.csv", header=T))


Roo.GPA <- gpagen(Roo)  


Roo.PCA=plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]])

plotRefToTarget(Roo.GPA$coords[,,1], Roo.GPA$coords[,,2], method="vector", label = T)
dimnames(Roo)


individualsBase=dimnames(Roo.GPA$coords)[[3]] 
means <- (aggregate(two.d.array(Roo.GPA$coords) ~ individualsBase, FUN=mean))[,-1] 
rownames(means) <- unique(dimnames(Roo.GPA$coords)[[3]]) 
Roo.GPA_Aver <- arrayspecs(means, p=ncol(means)/3, k=3)
Roo.PCA_Aver=plotTangentSpace(Roo.GPA_Aver, label=rownames(means) )

findMeanSpec(Roo.GPA_Aver) 

plotRefToTarget(Roo.GPA_Aver[,,14], Roo.GPA_Aver[,,15], method="vector", label = T)

CSize <- (aggregate(Roo.GPA$Csize ~ individualsBase, FUN=mean))[-1]
rownames(CSize) <- unique(dimnames(Roo.GPA$coords)[[3]])#

plotTangentSpace(Roo.GPA_Aver, label=dimnames(Roo.GPA_Aver)[[3]], warpgrids = T)
plotTangentSpace(Roo.GPA$coords, label=dimnames(Roo.GPA$coords)[[3]], warpgrids = T)


findMeanSpec(Roo.GPA$coords)
findMeanSpec(Roo.GPA_Aver)


######################################################################################

ref <- mshape(Roo.GPA$coords)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC1min)
plotRefToTarget(Roo[,,3], Roo[,,44], method="vector", label = T)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC2min)
plotRefToTarget(ref, Roo.PCA$pc.shapes$PC3min)

######################Define coords and Csize#####################################################

ibase=unique(individualsBase)
SpeciesGroupsSingle <- ibase 

SpeciesGroupsSingle<-as.factor(SpeciesGroupsSingle)
is.factor(SpeciesGroupsSingle) 
coords <- (aggregate(two.d.array(Roo.GPA_Aver) ~ SpeciesGroupsSingle, FUN=mean))[,-1]
rownames(coords) <- unique(SpeciesGroupsSingle)
coords <- arrayspecs(coords,p=dim(Roo.GPA_Aver)[1],k=3)

CSize=as.matrix(CSize)
CSize<- (aggregate(CSize ~ SpeciesGroupsSingle, FUN=mean))[,-1]
names(CSize)<-unique(SpeciesGroupsSingle)

coords 
CSize

PCA <- plotTangentSpace(coords, warpgrids=F) 
write.csv(PCA$pc.scores, file="PCScores.csv")
summary(PCA)

classifiers=read.table("classifiers_no_fossils.txt", header=T)
classifiers$Fossil <- c(rep(1,length(rownames(classifiers))));classifiers$Fossil [ c(which(classifiers$Diet=="fossil"))] <-2
rownames(classifiers) <- classifiers$Species

cols <- c("blue", "dark green", "red", "yellow")
names(cols)=levels(classifiers$Diet)
col.gp <- cols[match(classifiers$Diet, names(cols))]

point <-c(21,23,24,22,25)
names(point)=levels(classifiers$Diet)
point.gp <- point[match(classifiers$Diet, names(point))]

trees <- read.nexus("Roo.con.tre")

save.image(paste(root.file.path, "/Processed_data/dentary_without_fossils.rda", sep=""))

setwd(root.file.path)
