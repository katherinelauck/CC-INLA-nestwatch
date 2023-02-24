##### Implement INLA with nestwatch data to check for sensitivity to spatial autocorrelation
##### Author: Katherine Lauck
##### Last updated: 20230221

library(INLA); library(ggplot2); library(ggregplot)
library(tidyverse)
library(RColorBrewer)
library(magrittr)

d <- read_rds("success-cleaned.rds")

m <- read_rds("mainv1_withregion.rds")

head(d)

phen <- c("Region", "species","UnCoor","year","lon", "lat") # Base columns with spatial information we'll need

resp <- "at_least_one_success" # Response variable

covar <- c("Tmax_std_gridmet", # Temp anomaly
           "NewLU1", # Local land use (Ag, Forest, Natural_open, Human)
           "Tmax_std_gridmet_sq", # Temp anomaly squared
           "pcpbefore_raw_gridmet", # precipitation
           "NLCD_p_forest", # Landscape forest % cover
           "NLCD_p_human", # landscape human % cover
           "NLCD_p_ag", # landscape ag % cover
           "substrate_binary", # binary indication of whether in nestbox or not
           "laydate_scaled") # Julian date laydate

TestHosts <- na.omit(d[, c(phen, resp, covar)]) # Getting rid of NA's, picking adults
# We are using the [] to subset and only extract specific columns

TestHosts$UnCoor <- as.factor(TestHosts$UnCoor)

# Setting up a custom theme
THEME <- theme(axis.text.x = element_text(size = 12,colour = "black"),
               axis.text.y = element_text(size = 12, colour = "black"),
               axis.title.x = element_text(vjust = -0.35),
               axis.title.y = element_text(vjust = 1.2)) + theme_bw()

(samp_locations <- ggplot(TestHosts, aes(lon, lat)) + 
    geom_point() + 
    THEME)

fixed.effects <- formula(m) %>% as.character() %>% str_replace(pattern = " \\+ \\(.*\\)","")

mainv1_formula_noRegion <- as.formula(paste0(resp, " ~ ",
                                    fixed.effects[3],
                                    " + f(species, model = 'iid')",
                                    " + f(year, model = 'iid')",
                                    " + f(UnCoor, model = 'iid')"))


mainv1_formula_withRegion <- as.formula(paste0(resp, " ~ ",
                                        fixed.effects[3],
                                        " + f(species, model = 'iid')",
                                        " + f(year, model = 'iid')",
                                        " + f(UnCoor, model = 'iid')",
                                        " + f(Region, model = 'iid')"))

TestHosts_samp <- TestHosts # slice_sample(TestHosts,prop = .1)

inla_noautocorr <- inla(mainv1_formula_withRegion,
                       family = "binomial",
                       data = TestHosts_samp,
                       control.compute = list(residuals = TRUE))
write_rds(inla_noautocorr,"results/revisions/inla_noautocorr.rds")

Efxplot(inla_noautocorr)

Locations = cbind(TestHosts_samp$lon,TestHosts_samp$lat)

MeshA <- inla.mesh.2d(Locations, max.edge = c(1,1.2))
MeshB <- inla.mesh.2d(Locations, max.edge = c(2,2.2))
MeshC <- inla.mesh.2d(Locations, max.edge = c(.5,.7))

plot(MeshC)
points(Locations,col = "red",pch = 2)

# Making the A matrix

HostsA <- inla.spde.make.A(MeshC, loc = Locations) # Making A matrix
Hosts.spde = inla.spde2.pcmatern(mesh = MeshC, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.Host <- inla.spde.make.index('w', n.spde = Hosts.spde$n.spde) # making the w


# Making the model matrix #### 



X0 <- model.matrix(as.formula(paste0(" ~ -1 + ", fixed.effects[3])), data = TestHosts_samp) # make the model matrix using the final model selection formula without a response variable.

X <- as.data.frame(X0[,-which(colnames(X0)%in%c("NewLU1Forest"))]) # convert to a data frame. Eliminate the base level of the first categorical variable if applicable (you will manually specify an intercept below) 

head(X)

names(X) <- str_replace(names(X),"\\:",".") # interaction character not allowed in inla.stack. How to match names with formula?

# Making the stack ####

N <- nrow(TestHosts_samp)

StackHost <- inla.stack(
  data = list(y = TestHosts_samp[,resp]), # specify the response variable
  
  A = list(1, 1, 1, 1, 1, HostsA), # Vector of Multiplication factors for random and fixed effects              
  
  effects = list(
    
    Intercept = rep(1, N), # specify the manual intercept!
    
    X = X, # attach the model matrix
    
    species = TestHosts_samp$species, # insert vectors of any random effects
    UnCoor = TestHosts_samp$UnCoor,
    year = TestHosts_samp$year,
    
    w = w.Host)) # attach the w 

f_autocorr <- as.formula(paste0("at_least_one_success ~ -1 + Intercept + ",
                                paste0(colnames(X), collapse = " + "),
                                " + f(species, model = 'iid')",
                                " + f(year, model = 'iid')",
                                " + f(UnCoor, model = 'iid')"))

inla_autocorr <- inla(f_autocorr,
                      family = "binomial",
                      data = inla.stack.data(StackHost),
                      control.compute = list(dic = TRUE, residuals = TRUE),
                      control.predictor = list(A = inla.stack.A(StackHost)))
write_rds(inla_autocorr,"results/revisions/inla_autocorr.rds")

effects <- Efxplot(list(inla_noautocorr,inla_autocorr))
write_rds(effects, "results/revisions/effects_plot.rds")

SpatialHostList <- list(inla_noautocorr,inla_autocorr)

dic <- INLADICFig(SpatialHostList)
write_rds(dic, "results/revisions/dic_plot.rds")
