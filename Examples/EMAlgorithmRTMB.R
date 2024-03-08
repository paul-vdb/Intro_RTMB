library(RTMB)

## Test EM Algorithm in general for RTMB:

## Simulate some data:
mu <- c(2, 8, 15)
sd <- c(0.5, 1, 2)
wgts <- rgamma(3, 1, 1)
wgts <- wgts/sum(wgts)  ## dirichlet (1,1,1)

dat <- list()
dat$n <- 1000
dat$id <- sample(3, dat$n, replace = TRUE, prob = wgts)
dat$length <- rnorm(dat$n, mu[dat$id], sd[dat$id])
plot(density(dat$length))

## Create Expected value of the mixture weights function, given the previous data.
## All the previous data must be contained in pars.
pars <- list()
pars$mu <- c(0,0,0)
pars$logSD <- c(0,0,0)
pars$transP <- c(0,0)

transformP <- function (alpha){
	p <- 1/(1+exp(-alpha[1]))
	if(length(alpha) == 1) return(c(p, 1-p))
	for( i in 2:length(alpha) )	p <- c(p, 1/(1 + exp(-alpha[i])) * (1-sum(p[1:(i-1)])))
	return(c(p, 1-sum(p)))
}

iTransformP <- function(p){
  alpha <- log(p[1]/(1-p[1]))
  for( i in 2:(length(p)-1) ){
    pi <- p[i]/(1-sum(p[1:(i-1)]))
    alpha <- c(alpha, log(pi / (1-pi)) )
  }
  alpha
}

EMAlgorithm <- function(data, K, maxIter = 100, init = NULL){
  negLL <- function(pars){
    getAll(pars, data)
    sd <- exp(logSD)
    p <- transformP(transP)
    K <- length(mu)
    negLL <- 0
    ADREPORT(sd)
    ADREPORT(p)
    for( i in 1:n) {
      negLL <- negLL - log(sum(p*dnorm(length[i], mean = mu, sd = sd)))
    }
    negLL
  }

  Qfn <- function(pars){
    getAll(pars, data)
    sd <- exp(logSD)
    p <- transformP(transP)
    p0 <- transformP(pars0$transP)
    sd0 <- exp(pars0$logSD)
    mu0 <- pars0$mu
    K <- length(mu)
    Q <- 0
    wgts0 <- matrix(0,nrow=n,ncol=K)

    for( i in 1:n) {
      wgts0[i,] <- p0*dnorm(length[i], mean = mu0, sd = sd0)
      wgts0[i,] <- wgts0[i,]/sum(wgts0[i,])
      Q <- Q + 
        sum(wgts0[i,]*dnorm(length[i], mean = mu, sd = sd, log = TRUE) +
        wgts0[i,]*log(p))
    }
    -Q
  }

  ## smart init:
  if(is.null(init)){
    pars <- list()
    pars$mu <- quantile(data$length, cumsum(rep(1/K, K)) )
    pars$logSD <- log(rep(sd(data$length)/K, K))
    pars$transP <- iTransformP(rep(1/K, K))
    pars0 <- pars
  }else{
    pars <- init
    pars0 <- pars
  }
  Q0 <- -Inf
  Q1 <- Qfn(pars)
  iters <- 1
  
  while( abs(Q1-Q0) > 0.0000001 & iters < maxIter){
    obj <- MakeADFun(Qfn, pars0, silent = TRUE)
    opt <- nlminb (obj$par, obj$fn, obj$gr)
    pars0$mu <- obj$env$last.par.best[names(obj$env$last.par.best) == "mu"]
    pars0$logSD <- obj$env$last.par.best[names(obj$env$last.par.best) == "logSD"]
    pars0$transP <- obj$env$last.par.best[names(obj$env$last.par.best) == "transP"]
    Q0 <- Q1
    Q1 <- opt$objective
    iters <- iters + 1
  }

  obj <- MakeADFun(negLL, pars0, silent = TRUE)
  return(obj)
}

## Initialize at the real value for testing.
init <- list()
init$mu <- mu
init$logSD <- log(sd)
init$transP <- iTransformP(wgts)

objEM <- EMAlgorithm(data=dat, K=3, maxIter=100, init = init)
objEM$env$last.par.best
sdrep <- sdreport(objEM)
pl <- as.list(sdrep, "Est", report=TRUE)
plsd <- as.list(sdrep, "Std", report=TRUE)
objEM$fn(objEM$env$last.par.best)

## Compare with R
## Does not match. But this is a different Log Likelihood
library(mixdist)
Model.Constr <- mixconstr(conpi = "NONE", conmu = "NONE", consigma = "NONE")
init <- data.frame(pi = wgts, mu = mu, sigma= sd)
h.up <- hist(dat$length, plot = F, nclass = 20)
x.up <- h.up$breaks[2:length(h.up$breaks)]
y.up <- h.up$counts
dat2 <- as.mixdata(data.frame(x.up, y.up))
fitdat <- mix(mixdat = dat2, mixpar=init, dist = "norm", emsteps = 10)

## Does match! This is the same.
library(mixtools)
gm <- normalmixEM(dat$length, k=3, lambda=init$pi, mu=init$mu, sigma=init$sigma)
gm$lambda
gm$mu
gm$sigma

## Why we need EM algorithm:
negLL <- function(pars){
  getAll(pars, data)
  sd <- exp(logSD)
  p <- transformP(transP)
  K <- length(mu)
  negLL <- 0
  ADREPORT(sd)
  ADREPORT(p)
  for( i in 1:n) {
    negLL <- negLL - log(sum(p*dnorm(length[i], mean = mu, sd = sd)))
  }
  negLL
}

obj <- MakeADFun(negLL, pars, silent = TRUE)
opt <- nlminb (obj$par, obj$fn, obj$gr)
sdrep <- sdreport(objEM)
pl <- as.list(sdrep, "Est", report=TRUE)
plsd <- as.list(sdrep, "Std", report=TRUE)
