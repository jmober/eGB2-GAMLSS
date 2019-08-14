#########################################################################################
##### FLEXIBLE MODELING OF RATIO OUTCOMES IN CLINLICAL AND EPIDEMIOLOGICAL RESEARCH ##### 
#####                               Electronic Supplement                           #####
#########################################################################################
#####                               Author: Moritz Berger                           #####
#########################################################################################

## Functions: 
# eGB2_theta()
# eGB2_alpha()
# eGB2_rho() 
# eGB2LSS() 

## Description: 
# The functions implement the GAMLSS family 
# for fitting the extended GB2 model. 

# Note: The parameter theta corresponds to gamma in the manuscript. 

eGB2_theta <- function(theta, alpha, rho, stabilization){
  
  loss <- function(alpha, rho, y, f, w = 1) {
    yy  <- y/(1+y)
    my  <- 1 + (exp(-f)-1) * yy
    dmy <- my^2 - 4 * rho * exp(-f) * yy * (1-yy)
    fy  <- lgamma(2*alpha) - 2*lgamma(alpha) - alpha * f + 
      alpha * log(1-rho) + log(my) + (alpha-1)*log(yy) +(alpha+1)*log(1-yy) -
      (alpha+0.5) * log(dmy)
    return(-fy)
  }
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, alpha=alpha, rho=rho))
  }
  ngradient <- function(y, f, w = 1) {
    yy  <- y/(1+y)
    my  <- 1 + (exp(-f)-1) * yy
    dmy <- my^2 - 4 * rho * exp(-f) * yy * (1-yy)
    ngr <- -alpha - (yy * exp(-f)) / my + (2 * alpha + 1) * exp(-f) * yy * (my - 2 * rho * (1-yy)) / dmy
    ngr <- gamboostLSS:::stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr)
  }
  offset <- function(y, w = 1) {
    0
  }
  response <- function(f) {
    exp(f)
  }
  Family(ngradient = ngradient, risk = risk, offset = offset, 
         check_y = function(y) {
           stopifnot(all(y > 0))
           y
         }, response = response, 
         rclass = function(f) f, name = "eGB2 theta")
  
}

eGB2_alpha <- function(theta, alpha, rho, stabilization){
  
  loss <- function(theta, rho, y, f, w = 1) {
    yy  <- y/(1+y)
    my  <- 1 + (1/theta-1) * yy
    dmy <- my^2 - 4 * rho / theta * yy * (1-yy)
    fy  <- lgamma(2*exp(f)) - 2*lgamma(exp(f)) - exp(f)*log(theta) +
      exp(f) * log(1-rho) + log(my) + (exp(f)-1)*log(yy) +(exp(f)+1)*log(1-yy) -
      (exp(f)+0.5) * log(dmy)
    return(-fy)
  }
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, theta=theta, rho=rho))
  }
  ngradient <- function(y, f, w = 1) {
    yy  <- y/(1+y)
    my  <- 1 + (1/theta-1) * yy
    dmy <- my^2 - 4 * rho / theta * yy * (1-yy)
    ngr <- exp(f)*(2*digamma(2*exp(f))-2*digamma(exp(f))-log(theta)+log(1-rho)+log(y)-2*log(1+y)-log(dmy))
    ngr <- gamboostLSS:::stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr)
  }
  offset <- function(y, w = 1) {
    0
  }
  response <- function(f) {
    exp(f)
  }
  Family(ngradient = ngradient, risk = risk, offset = offset, 
         check_y = function(y) {
           stopifnot(all(y > 0))
           y
         }, response = response, 
         rclass = function(f) f, name = "eGB2 alpha")
  
}

eGB2_rho <- function(theta, alpha, rho, stabilization){
  
  loss <- function(theta, alpha, y, f, w = 1) {
    e1e <- exp(f)/(1+exp(f))
    yy  <- y/(1+y)
    my  <- 1 + (1/theta-1) * yy
    dmy <- my^2 - 4 * e1e / theta * yy * (1-yy)
    fy  <- lgamma(2*alpha) - 2*lgamma(alpha) - alpha*log(theta) +
      alpha * log(1-e1e) + log(my) + (alpha-1)*log(yy) +(alpha+1)*log(1-yy) -
      (alpha+0.5) * log(dmy)
    return(-fy)
  }
  risk <- function(y, f, w = 1) {
    sum(w * loss(y = y, f = f, theta=theta, alpha=alpha))
  }
  ngradient <- function(y, f, w = 1) {
    e1e  <- exp(f)/(1+exp(f))
    e1e2 <- exp(f)/(1+exp(f))^2
    yy   <- y/(1+y)
    my   <- 1 + (1/theta-1) * yy
    dmy  <- my^2 - 4 * e1e / theta * yy * (1-yy)
    ngr  <- -alpha*e1e2/(1-e1e)+(alpha+0.5)/dmy/theta*4*yy*(1-yy)*e1e2
    ngr  <- gamboostLSS:::stabilize_ngradient(ngr, w = w, stabilization)
    return(ngr)
  }
  offset <- function(y, w = 1) {
    0
  }
  response <- function(f) {
    exp(f)/(1+exp(f))
  }
  Family(ngradient = ngradient, risk = risk, offset = offset, 
         check_y = function(y) {
           stopifnot(all(y > 0))
           y
         }, response = response, 
         rclass = function(f) f, name = "eGB2 rho")
  
}

eGB2LSS <- function (theta = NULL, alpha = NULL, rho = NULL,
                            stabilization = c("none", "MAD", "L2")){
                            stabilization <- gamboostLSS:::check_stabilization(stabilization)
                            Families(theta=eGB2_theta(theta, alpha, rho, stabilization = stabilization),
                                     alpha=eGB2_alpha(theta, alpha, rho, stabilization = stabilization),
                                     rho=eGB2_rho(theta, alpha, rho, stabilization = stabilization))
                            }


