% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selfmade.R
\name{mocasin}
\alias{mocasin}
\title{Function which computes selective p-values and intervals for \code{gamm4} and 
\code{merMod} objects}
\usage{
mocasin(mod, checkFun, this_y = NULL, nrSamples = 1000,
  bayesian = TRUE, varInTestvec = c("est", "minMod", "varY",
  "supplied"), varForSampling = c("est", "minMod", "varY", "supplied"),
  VCOV_vT = NULL, VCOV_sampling = NULL, conditional = TRUE,
  name = NULL, nrlocs = 7, which = NULL, vT = NULL, G = NULL,
  efficient = TRUE, trace = TRUE)
}
\arguments{
\item{mod}{an object of class \code{merMod} or result of \code{gamm4} function}

\item{checkFun}{a function of \code{y}, a vector of the same length as 
the original response vector which returns \code{TRUE} or \code{FALSE} depending
on whether the selection event for a given \code{y} corresponds to the 
original model selection. See the example for more details.}

\item{this_y}{original response vector (explicit reference may be necessary 
for certain model classes)}

\item{nrSamples}{integer; the number of Monte Carlo samples to be used for 
inference (defaults to 1000)}

\item{bayesian}{logical; whether or not to use a bayesian type covariance}

\item{varInTestvec}{for expert use only; variance used in the test vector definition}

\item{varForSampling}{variance used for inference; per default the estimated variance 
of \code{mod} is used. 
Other options are a conservative estimate based on the variance of the 
response is used ("varY") or to supply a numeric value to base inference 
on a customize variance}

\item{VCOV_vT}{for expert use only; VCOV used in the test vector definition}

\item{VCOV_sampling}{covariance matrix of dimension of the response used for inference; 
per default the estimated covariance of \code{mod} is used. 
Otherwise a matrix must be supplied on which basis inference is conducted. 
If the true 
covariance is unknown, an conservative alternative to plugging in the 
estimator is given 
by using the covariance of the refitted mixed model, for which all fixed 
effects but the intercept 
are excluded.}

\item{conditional}{logical; determines whether to use the conditional or 
marginal approach
when \code{mod} is of class \code{merMod}, i.e., inference is sought for a 
linear mixed model}

\item{name}{character; for the \code{gamm4}-case: the name of the covariate, 
for which inference is done}

\item{nrlocs}{integer; for the \code{gamm4}-case: the number of locations 
tested for non-linear effects}

\item{which}{integer; for the \code{merMod}-case: defining the effect for 
which inference is done}

\item{vT}{list of vectors (optional); if inference is sought for a customized 
test vector, this argument
can be used}

\item{G}{true random effect covariance (optional)}

\item{efficient}{logical; whether or not to compute the test statistic based on
an (efficient) weighted LS estimator instead of a OLS estimator for the marginal model}

\item{trace}{logical; if TRUE, a progress bar is printed in the console}
}
\description{
Function which computes selective p-values and intervals for \code{gamm4} and 
\code{merMod} objects
}
\details{
Note that the additive and conditional mixed model approach 
currently only works for a diagonal error covariance.
}
\examples{

library(lme4)
if(require(lmerTest)){

##### BASED ON lmerTest HELP PAGE #########
# define function to fit a model based on response
modFun <- function(y)
{
  ham$y <- y
  lmer(y ~ Gender + Information * Product + (1 | Consumer) +
  (1 | Product), data=ham)
  
  }

# define a function to select a model (based on a model)
selFun <- function(mod) step(mod)

# define a function which extracts the results of the selection procedure
extractSelFun <- function(this_mod){

this_mod <- attr(this_mod, "model")
if(class(this_mod)=="lm") 
  return(attr(this_mod$coefficients, "names")) else
    return(c(names(fixef(this_mod)), names(getME(this_mod, "theta"))))

}


## backward elimination of non-significant effects:
(step_result <- selFun(modFun(ham$Informed.liking)))
attr(step_result, "model")
## Elimination tables for random- and fixed-effect terms:
(sel <- extractSelFun(step_result))

## Now we can finally define the function checking the congruency
## with the original selection

checkFun <- function(yb){ 

 this_mod <- modFun(yb)
 setequal( extractSelFun(selFun(this_mod)), sel )
 
 }

# Now let's compute valid p-values conditional on the selection
res <- mocasin(attr(step_result, "model"), this_y = ham$Informed.liking,
          checkFun = checkFun, which = 1:4, nrSamples = 50, trace = FALSE)

# print(res)
} 

# gamm4 example similar to the one from gamm4 help page
if(require(gamm4)){
set.seed(0) 
dat <- gamSim(1,n=500,scale=2) ## simulate 4 term additive truth

dat$y <- 3 + dat$x0^2 + rnorm(n=500)
br <- gamm4(y~ s(x0) + s(x1), data = dat)
summary(br$gam) ## summary of gam

# do not use any selection
checkFun <- function(yb) TRUE

res <- mocasin(br, this_y = dat$y,
          checkFun = checkFun, 
          nrlocs = c(0.7,1), 
          nrSamples = 100)
}
}
