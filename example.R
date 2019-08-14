#########################################################################################
##### FLEXIBLE MODELING OF RATIO OUTCOMES IN CLINLICAL AND EPIDEMIOLOGICAL RESEARCH ##### 
#####                               Electronic Supplement                           #####
#########################################################################################
#####                               Author: Moritz Berger                           #####
#########################################################################################


# Exemplary analysis fitting an eGB2 GAMLSS 

# Note: The parameter theta corresponds to gamma in the manuscript. 

# load library and functions
library("gamboostLSS")
source("gamboostLSSFamily_eGB2.R")

# load simulated data set 
load("sim_data.rda")
# the data set contains the outcome variable y simulated from the eGB2 distribution
# and four explanatory variables x2,...,x5 that were related to theta 


# build formula for gamboostLSS using linear base-learners 
form <- paste0("bols(", names(data[-1]), ", intercept=F)", collapse="+")
form <- formula(paste0("y~", form))

# build matrix for cross-validation 
folds      <- cv(rep(1,500), type="kfold", B=10)

##########

# fit eGB2 GAMLSS  
mod_eGB2 <- gamboostLSS(formula = list(theta = form,
                                       alpha = form,
                                       rho   = form),
                        families = eGB2LSS(stabilization = "none"), data = data, method = "noncyclic",
                        control = boost_control(mstop = 1000, trace = T, nu = c(theta = 0.1, alpha = 0.1, rho = 0.1)))

# select optimal stopping iteration 
sel_eGB2    <- cvrisk(mod_eGB2, folds=folds, trace=T)
plot(sel_eGB2)
opt_eGB2    <- mstop(sel_eGB2)
modopt_eGB2 <- mod_eGB2[opt_eGB2]

# coefficient estimates 
coef(modopt_eGB2)$theta
coef(modopt_eGB2)$alpha
coef(modopt_eGB2)$rho

# fitted values 
boxplot(predict(modopt_eGB2, type="response", parameter="theta"))
boxplot(predict(modopt_eGB2, type="response", parameter="alpha"))              
boxplot(predict(modopt_eGB2, type="response", parameter="rho"))

                
