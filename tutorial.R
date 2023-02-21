
if(!require(ggregplot)) devtools::install_github("gfalbery/ggregplot") # Installing Greg's package for plotting functions!

library(INLA); library(ggplot2); library(ggregplot)
library(tidyverse)
library(RColorBrewer)

Root <- getwd() # This should be the path to your working directory

Hosts <- read.csv(paste0(Root, "/HostCaptures.csv"), header = T)

head(Hosts)

substr(names(Hosts), 1, 1) <- toupper(substr(names(Hosts), 1, 1)) # Giving the host names capital letters

phen <- c("Grid", "ID", "Easting", "Northing") # Base columns with spatial information we'll need

resp <- "Parasite.count" # Response variable

covar <- c("Month", # Julian month of sampling
           "Sex", # Sex
           "Smi", # Body condition
           "Supp.corrected", # Nutrition supplementation
           "Treated") # Treatment

TestHosts <- na.omit(Hosts[, c(phen, resp, covar)]) # Getting rid of NA's, picking adults
# We are using the [] to subset and only extract specific columns

# Turning variables into factors
TestHosts$Month <- as.factor(TestHosts$Month)
TestHosts$Grid <- as.factor(TestHosts$Grid)

TestHosts$Parasite.count <- round(TestHosts$Parasite.count) # Parasite counts should be integers

table(table(TestHosts$ID)) # Enough repeat samples for a mixed model?

# Setting up a custom theme
THEME <- theme(axis.text.x = element_text(size = 12,colour = "black"),
               axis.text.y = element_text(size = 12, colour = "black"),
               axis.title.x = element_text(vjust = -0.35),
               axis.title.y = element_text(vjust = 1.2)) + theme_bw()

(samp_locations <- ggplot(TestHosts, aes(Easting, Northing)) + 
    geom_jitter(aes(colour = factor(Grid))) + coord_fixed() + 
    THEME + 
    labs(colour = "Grid"))

length(unique(TestHosts$ID))

table(with(TestHosts, tapply(Grid, ID, function(x) length(unique(x)))))

# First without random effects ####

# Specify the formula
f0.1 <- as.formula(paste0(resp, " ~ ", # Response first
                          paste(covar, collapse = " + ") # Collapse the vector of covariates
))

# Run the model
IM0.1  <- inla(Parasite.count ~ Month + Sex + Smi + Supp.corrected + Treated, 
               family = "nbinomial", # Specify the family. Can be a wide range (see r-inla.org).
               data = TestHosts) # Specify the data

# Run the model # (This is the same thing)
IM0.1  <- inla(f0.1, 
               family = "nbinomial", # Specify the family. Can be a wide range (see r-inla.org).
               data = TestHosts) # Specify the data

# Then with an ID random effect ####

f0.2 <- as.formula(paste0(resp, " ~ ", 
                          paste(covar, collapse = " + "), 
                          " +  f(ID, model = 'iid')")) # This is how you include  a typical random effect.

IM0.2  <- inla(f0.2, 
               family = "nbinomial",
               data = TestHosts) 

summary(IM0.1)
summary(IM0.2)

library('magrittr')
Efxplot(list(IM0.1, IM0.2))

# Let's try it on our data ####

HostModelSel <- INLAModelSel(resp, covar, "ID", "iid", "nbinomial", TestHosts)

Finalcovar <- HostModelSel$Removed[[length(HostModelSel$Removed)]]

f1 <- as.formula(paste0(resp, " ~ ", 
                        paste(Finalcovar, collapse = " + "), 
                        "+ f(ID, model = 'iid')")) 

IM1 <- inla(f1,
            family = "nbinomial",
            data = TestHosts,
            control.compute = list(dic = TRUE)) 

summary(IM1)

Locations = cbind(TestHosts$Easting, TestHosts$Northing) # using the sampling locations 

MeshA <- inla.mesh.2d(jitter(Locations), max.edge = c(20, 40))
MeshB <- inla.mesh.2d(Locations, max.edge = c(20, 40))
MeshC <- inla.mesh.2d(Locations, max.edge = c(10, 20))

Mesh <- MeshB

plot(MeshA)

plot(MeshB)

plot(MeshC)

points(Locations, col = "red", pch = 2)

# Making the A matrix

HostsA <- inla.spde.make.A(Mesh, loc = Locations) # Making A matrix
Hosts.spde = inla.spde2.pcmatern(mesh = Mesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.Host <- inla.spde.make.index('w', n.spde = Hosts.spde$n.spde) # making the w



# Making the model matrix #### 

X0 <- model.matrix(as.formula(paste0(" ~ -1 + ", paste(Finalcovar, collapse = " + "))), data = TestHosts) # make the model matrix using the final model selection formula without a response variable.

X <- as.data.frame(X0[,-which(colnames(X0)%in%c("Month7"))]) # convert to a data frame. Eliminate the base level of the first categorical variable if applicable (you will manually specify an intercept below) 

head(X)

# Making the stack ####

N <- nrow(TestHosts)

StackHost <- inla.stack(
  data = list(y = TestHosts[,resp]), # specify the response variable
  
  A = list(1, 1, 1, HostsA), # Vector of Multiplication factors for random and fixed effects              
  
  effects = list(
    
    Intercept = rep(1, N), # specify the manual intercept!
    
    X = X, # attach the model matrix
    
    ID = TestHosts$ID, # insert vectors of any random effects
    
    w = w.Host)) # attach the w 

f1 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + ")))
f2 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), " +  f(ID, model = 'iid')"))
f3 <- as.formula(paste0("y ~ -1 + Intercept + ", paste0(colnames(X), collapse = " + "), " +  f(ID, model = 'iid') + f(w, model = Hosts.spde)"))


IM1 <- inla(f1, # Base model (no random effects)
            family = "nbinomial",
            data = inla.stack.data(StackHost),
            control.compute = list(dic = TRUE),
            control.predictor = list(A = inla.stack.A(StackHost))
)

IM2 <- inla(f2, # f1 + Year and ID random effects
            family = "nbinomial",
            data = inla.stack.data(StackHost),
            control.compute = list(dic = TRUE),
            control.predictor = list(A = inla.stack.A(StackHost))
)

IM3 <- inla(f3, # f2 + SPDE random effect 
            family = "nbinomial",
            data = inla.stack.data(StackHost),
            control.compute = list(dic = TRUE),
            control.predictor = list(A = inla.stack.A(StackHost))
)

SpatialHostList <- list(IM1, IM2, IM3)

ggField(IM3, Mesh, Groups = 1) +
  scale_fill_brewer(palette = "Blues") 

# always use a single-dimension colour palette if you can! It's just easier on the eyes, 
# better for colourblind people, makes sense in black and white, etc.

# ignore the Groups part of the function for now. That'll come later.

# function takes (a list of) models and plots the decay of spatial autocorrelation across a user-defined range

# let's try it on our model ###

# Define the maximum range as something reasonable: the study area is 80 eastings wide, so lets go for:

Maxrange = 40

INLARange(list(IM3), maxrange = Maxrange)

sapply(SpatialHostList, function(f) f$dic$dic)


# Let's try it on our data ####

INLADICFig(SpatialHostList, ModelNames = c("Base", "IID", "SPDE"))
