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

iter = 10000 # For all analyses using iterations

####~~~~~~~~~~~~PHYLOGENETIC ANALYSES~~~~~~~~~~~~~~~~~###

  #Making sure tree is set and has Grafen-ultrametricized branch lengths
trees <- compute.brlen(trees,  method="Grafen")
plot(trees)

###################Caper code########################
###We used mequite to add species to the Butler et al.(2018) matrix

library(geiger)
NamesRoo=as.character(classifiers$Species)
RooData=classifiers
RooData=cbind(NamesRoo, RooData) ##RooData will be body mass calculations, diet etc.
RooCD=comparative.data(trees, RooData, names.col =NamesRoo) 
name.check(RooCD$phy, RooCD$data)

#additional check
RooCD$phy$tip.label== rownames(RooCD$data) #all data points = true

phylo_coords <- coords[,,RooCD$phy$tip.label]
phylo_CSize <- as.matrix(CSize[RooCD$phy$tip.label]) 
phylo_classifiers <- classifiers[RooCD$phy$tip.label,] # reorders to tip labels
plot(RooCD$phy)

name.check(RooCD$phy,PCA$pc.scores)

### Evolutionary Allometry - using CSize ###

phylo_allom <- RooCD$phy
phylo_allom <- procD.pgls(phylo_coords ~ log(phylo_CSize), RooCD$phy,iter=iter) #Phylogenetic ANOVA/Regression For Shape Data
summary(phylo_allom)

### Evolutionary Allometry - Using Body Mass (g)
rpadataframe = geomorph.data.frame(phylo_shape = phylo_coords, phylo_mass = phylo_classifiers$Bodymass_g)
phylo_allo_mass <- RooCD$phy
phylo_allo_mass <- procD.pgls(phylo_shape ~ log(phylo_mass), RooCD$phy,iter=iter, data=rpadataframe)
summary(phylo_allo_mass)

#####Plotting body mass vs csize
plot((classifiers$Bodymass_g )~  (CSize))
text((classifiers$Bodymass_g )~ (CSize), labels=classifiers$Species)

####Plot tree using Diet classifer just for a brief check
plot(RooCD$phy, type="phylogram", cex=0.5, tip.color = cols[phylo_classifiers$Diet] )

######~~~~~~~~~~~~How does size relate to the main variation as per the first few PCs?~~~~~~~~~~~~~~~~~###

summary(lm(PCA$pc.scores[,1] ~ log(CSize))) # R2 =  0.3038 

summary(lm(PCA$pc.scores[,2] ~ log(CSize))) # R2 =  0.1706  
summary(lm(PCA$pc.scores[,3] ~ log(CSize))) # R2 =  0.1128 

#Body mass
summary(lm(PCA$pc.scores[,1] ~ classifiers$Bodymass_g)) # R2 =  0.2597  
summary(lm(PCA$pc.scores[,2] ~ classifiers$Bodymass_g)) # R2 = -0.01161  
summary(lm(PCA$pc.scores[,3] ~ classifiers$Bodymass_g)) # R2 =  0.1815   

###########~~~~~~~~~~~~~~~~~~~~~~~~~Analyses below to run witout including fossils to assess known dietary/locomotor shapes~~~~~~~~~~~~~~~~~~~~~~~~~~#########


####~~~~~~~~~~~~Procrustes analysis of Diet and locmotor mode~~~~~~~~~~~~~~~~~###

##make dataframe to analyses data below

geodataframe <- geomorph.data.frame(shape = coords, size = as.matrix(CSize), diet = classifiers$Diet, clades = classifiers$Clade, species = classifiers$Species, movement = classifiers$Movement) # make the geomorph data frame

# how much of the shape variation can be explained by diet or locomotion? 

summary(procD.lm(shape ~ log(size) + diet, data = geodataframe)) 
summary(procD.lm(log(size) ~ diet, data = geodataframe)) 

summary(procD.lm(shape ~ log(size) + movement, data = geodataframe))
summary(procD.lm(log(size) ~ movement, data = geodataframe)) 

####~~~~~~~~~~~~Diet and locmotion in phlyogenetic framework~~~~~~~~~~~~~~~~~###

#Run all code without inclusion of fossils for this section

#Diet in a phylogenetic framework
rpadataframe = geomorph.data.frame(phylo_shape = phylo_coords, phylo_CSize=phylo_CSize, phylo_Movement = phylo_classifiers$Diet)

phylo_allo_diet <- RooCD$phy
phylo_allo_diet <- procD.pgls(phylo_coords ~ log(phylo_CSize) + phylo_classifiers$Diet, RooCD$phy,iter=iter, data=rpadataframe)
summary(phylo_allo_diet)

### Locomotion in a phylogenetic framework
rpadataframe = geomorph.data.frame(phylo_shape = phylo_coords, phylo_CSize = phylo_CSize, phylo_Movement = phylo_classifiers$Movement)
phylo_allo_movement <- RooCD$phy
phylo_allo_movement <- procD.pgls(phylo_coords ~ log(phylo_CSize)+phylo_Movement, RooCD$phy,iter=iter, data=rpadataframe)
summary(phylo_allo_movement)

