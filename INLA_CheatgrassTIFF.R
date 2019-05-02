# INLA_Cheatgrass assessing the effect that fire count has on percent cheatgrass cover with spatial autocorrelation 

library("INLA")
library('raster')
library('INLA')
library('gridExtra')
library('plotly')
library('lattice')
library('rgdal')


setwd("H:/HierarchicalModel/Semester_Project") 

#Load data 
fire_shp <- st_read("Fire_count.shp") # Fire shapefile maybe  ... 
fire <- raster('Raster_fire9.tif') # Bring in fire count raster, same cell size as BRTE raster
BRTE <- raster('BRTE_clip6.tif') #Load BRTE model raster


# Rasterize the fire data 
#r <- raster(ncols= ncol(BRTE), nrows=nrow(BRTE))
#rast_ext <- extent(BRTE)
#fire_rast <- rasterize(fire_shp,r)

# Visualize data intially 
#plot(BRTE)
#plot(fire)


# Set extent of fire to match extent of BRTE 
compareRaster(fire, BRTE) # comares the extents of raster, must be TRUE
extent(BRTE) # shows extent of BRTE
extent(fire) # shows extent of fire 
b <- extent(BRTE) # sets b as extent of BRTE
Fire_2 <- setExtent(fire,b,keepres = TRUE) #modifies extent of fire data to that of BRTE, Keppers tell to keep columns as true 
extent(Fire_2) # checks the extent of fire_2
compareRaster(Fire_2, BRTE) #ensures that extent is true now 

#Stack the rasters on top of each other
BRTE_stack <- stack(Fire_2, BRTE)
# creates a dataframe of the stacked raster 
BRTE_data <- rasterToPoints(BRTE_stack) # takes a little bit to run ~ 5 min
BRTE_data2 <- as.data.frame(BRTE_data) # puts it into a workable data frame 
head(BRTE_data)

length(which(BRTE_data2$Raster_fire9>3))
hist(BRTE_data2$BRTE_clip6)
# Subset the data: 
set.seed(500)
I <- sample(1:nrow(BRTE_data2),100) 
BRTE_samp <- BRTE_data2[I,]

# Spatially subset the data: create folds in the data - get data from Juami ... 
head(BRTE_samp[,c(1,2)])

points(BRTE_samp[,c(2)],BRTE_samp[,c(1)], col = "red", pch = 16) # points are not correct ... 

# Start INLA Model: 
# Step one: create mesh 
mesh3 <- inla.mesh.2d(loc = BRTE_samp[,c(2,1)], max.edge = c(0.03,0.2), cutoff = .01)

# loc = information about your spatial domain,your data with the points to create the mesh , can be all data, C(1,2) call the col with coordinates 
# larger cutoff is less verticies 
# What is the difference between max.edge and cutoff smaller the cutoff the smaller the max.edge can be 

plot(mesh3)
points(BRTE_samp[,c(2)],BRTE_samp[,c(1)], col = "red", pch = 16)


#2 
#SPE is the sochastic part - partial differential equation 
spde2 <- inla.spde2.matern(mesh3)



# to assist with keeping track of which elements relate to what effects, I can create an index
s.index <- inla.spde.make.index("spatial.field", n.spde = spde2$n.spde)

# to create a projector matrix (links the m mesh vertices to the n responses)
A_matrix <- inla.spde.make.A(mesh3, loc = as.matrix(BRTE_samp[,c(2,1)]))


# Uses gaussian markov random fields 

# Step 6 - stack the data , need to define the intercept (b/c you need to provide the intercept)
# to create the stack - puts all the pieces together, includes the data, response, predictor, and effects 
stack.brte <- inla.stack(data = list(y = BRTE_samp$BRTE_clip6), A = list(A_matrix, 1), effects = list(s.index, list(n_burns = BRTE_samp$Raster_fire9, tag = "Burns"))) # What you would feed into stan or jags? 


# Y = response variable 


# Step 7. define your INLA model - like stanglm or something for inputs 
# run the model including spatial random effect
cheatgrass.model.inla <- inla(y ~ -1 + n_burns + f(spatial.field, model = spde2)
                              , data = inla.stack.data(stack.brte),
                              control.predictor = list(A = inla.stack.A(stack.brte),compute=TRUE),family = "nbinomial",verbose = T)

summary(cheatgrass.model.inla)
plot(cheatgrass.model.inla$marginals.fixed$n_burns,type="l")



# project the spatial random effect
gproj <- inla.mesh.projector(mesh3)
g.mean <- inla.mesh.project(gproj, cheatgrass.model.inla$summary.random$spatial.field$mean)
g.sd <- inla.mesh.project(gproj, cheatgrass.model.inla$summary.random$spatial.field$sd)
grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions=terrain.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd',col.regions=terrain.colors(16)),nrow=1)
par(mfrow = c(1, 2))
plot(lossYearMcCall, main = 'Deforestation 2001 - 2018 (Mcall window)')
rotate <- function(x) t(apply(x, 2, rev))
r<-raster(rotate(rotate(rotate(g.sd))))
plot(r, xlab='', ylab='', main='GF representation')
dev.off()

# plot model residuals
fitted <- deforest.model.inla$"summary.fitted.values"[,1][1:length(id)]
residuals <- (fitted - trainingData$deforestation)
plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')

# get the spatial parameters of the spatial random effect
spatial.parameters <- inla.spde2.result(inla = deforest.model.inla, name = "spatial.field", spde = spde, do.transform = T)
# nominal variance (the sill)
sill <- inla.emarginal(function(x) x, spatial.parameters$marginals.variance.nominal[[1]])
sill
# plot posterior distribution
plot(spatial.parameters$marginals.variance.nominal[[1]],type='l',main='Sill')

# range
range <- inla.emarginal(function(x) x, spatial.parameters$marginals.range.nominal[[1]])
range
# plot posterior distribution
plot(spatial.parameters$marginals.range.nominal[[1]],type='l',main='Range')

# nugget
nugget <- inla.emarginal(function(x) 1/x, deforest.model.inla$marginals.hyperpar$`Precision for the Gaussian observations`)
nugget
# plot posterior distribution
plot(inla.tmarginal(function(x) 1/x, deforest.model.inla$marginals.hyperpar$`Precision for the Gaussian observations`),type='l',main='Nugget')

# model validation with test data
nn <- setdiff(c(1:n),id)
sampleidTest <- sample(nn,5000)
coordSampleTest <- xyFromCell(lossYearMcCall, sampleidTest)
testData <- data.frame(coordSampleTest, deforestation = extract(lossYearMcCall, coordSampleTest))

# matrix for the predicitons
A_pred <- inla.spde.make.A(mesh = mesh1, loc = coordSampleTest)
## stack
stackP <- inla.stack(data = list(y = NA),
                     A = list(A_pred, 1), effects = list(s.index,list(b0=rep(1,length(testData$deforestation))
                                                                      ,remove.unused=TRUE,compress=T, tag = "spde")), tag = "predictions")
stackPRED<-inla.stack(stack.deforest,stackP)
# run the model
modINLAP<-inla(y ~ -1 + b0 + f(spatial.field, model = spde), data = inla.stack.data(stackPRED,spde=spde),family = "gaussian",
               control.predictor=list(A=inla.stack.A(stackPRED),compute =TRUE),verbose=T)
# r2
observed<-(testData$deforestation)
meanOb <- mean(observed)
idT<-inla.stack.index(stack=stackPRED , tag='predictions')
fitted<-modINLAP$summary.fitted.values$mean[idT$data]
residuals <- (fitted - testData$deforestation)
numerator<- sum((fitted - meanOb)^2)
denominator <- sum((observed - meanOb)^2)
r2<- (numerator) / (denominator)
print('R2');r2
plot(fitted,residuals,main=('Residuals vs Fitted (Including spatial term)')); abline(h=0,lty='dashed',col='red')

plot(cheatgrass.model.inla$marginals.fixed$n_burns,type="l")

# Get a burn index from landsat or sentinel - for raster laster, spatial resolution quesiton, compare MODIS and effect of fire with different types of data 
# MOD Is everywhere 
