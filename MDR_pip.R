### Code for brief crash course on Bayesian TPC estimation.

# Note: Please change this to your working directory!
setwd("/Users/mauricio/git/bayesTPC")

set.seed(42) # Set seed for reproducibility.

# We will need these packages. Please make sure you install them beforehand!
library('R2jags')
library('mcmcplots')
library('MCMCvis')

# Read data.
data <- read.csv("MDR_pipiens.csv")
data

plot(data$T, 1 / data$value, pch=20, xlab="Temperature [°C]",
     ylab="Developmental Rate (eggs)", col="steelblue", xlim=c(10, 40),
     ylim=c(0,0.12))

# We know what temperatures mean. But sometimes our model has parameters that
# don't have a straightforward interpretation. To set a reasonable prior 
# distribution for these parameters, it's a good idea to plot our function.
briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax)  )
  return(result)
}

T <- seq(5, 40, 0.01)

# Pick some ballpark reasonable values of Tmin and Tmax,
Tmin = 10
Tmax = 37

# Try different values of c to get an idea of what this parameter does to the
# curve and what would be reasonable values.
plot(T, briere(T, Tmin, Tmax, 0.01), type="l")
plot(T, briere(T, Tmin, Tmax, 0.005), type="l")

plot(T, briere(T, Tmin, Tmax, 0.001), type="l")
plot(T, briere(T, Tmin, Tmax, 0.0001), type="l")

# If we assume that developmental rates are usually between zero and one, 
# reasonable values of c should probably be below 0.001 or so. To be safe, let's
# assume they are below 0.01.

# Plot Topt for these parameters to check calculation is correct. We will need
# this formula later to find optimal temperature from fitted model parameters.
abline(v=(4*Tmax + 3*Tmin)/10 + sqrt(4 * Tmax^2 + (9/4) * Tmin^2 - 4*Tmax*Tmin) / 5)
## Note: This is the formula for the Briere model. It will not work for other
## models.


# Simple model with uniform priors.
sink("briere_unif.txt")
cat("
    model{
    
    ## Priors
    c ~ dunif(0, 0.01)
    Tmin ~ dunif(0, 20)
    Tmax ~ dunif(28, 45)
    sigma ~ dunif(0, 0.2)
    tau <- 1 / (sigma * sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
      trait.mu[i] <- c * temp[i] * (temp[i] - Tmin) * sqrt((Tmax - temp[i]) * (Tmax > temp[i])) * (Tmin < temp[i])
      trait[i] ~ dnorm(trait.mu[i], tau)
    }
    
    ## Derived Quantities
    Topt = (4*Tmax + 3*Tmin)/10 + sqrt(4 * Tmax^2 + (9/4) * Tmin^2 - 4*Tmax*Tmin) / 5
    
    } # close model
    ")
sink()

# Simple model with informative priors for temperature parameters.
sink("briere_inf.txt")
cat("
    model{
    
    ## Priors
    c ~ dunif(0, 0.001)
    Tmin ~ dnorm(10, 1 / 25) # 95% likely to be in (0, 20)
    Tmax ~ dnorm(35, 1 / 25) # 95% likely to be in (25, 45)
    sigma ~ dunif(0, 0.2)
    tau <- 1 / (sigma * sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
      trait.mu[i] <- c * temp[i] * (temp[i] - Tmin) * sqrt((Tmax - temp[i]) * (Tmax > temp[i])) * (Tmin < temp[i])
      trait[i] ~ dnorm(trait.mu[i], tau)
    }
    
    ## Derived Quantities
    Topt = (4*Tmax + 3*Tmin)/10 + sqrt(4 * Tmax^2 + (9/4) * Tmin^2 - 4*Tmax*Tmin) / 5
    
    } # close model
    ")
sink()

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  c = runif(0.0001, 0.001),
  Tmax = runif(1, min=35, max=40),
  Tmin = runif(1, min=5, max=10),
  sigma = runif(1, min=0.001, max=0.01))}

##### Parameters to Estimate
parameters <- c("c", "Tmin", "Tmax","sigma", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Organize Data for JAGS
trait <- 1 / data$value
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
##### Run JAGS
unif.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_unif.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())


# Examine output.
unif.out
round(unif.out$BUGSoutput$summary, 6)
mcmcplot(unif.out)

# We can save output for further processing.
#save(unif.out, file = "MDR_Cpip_jags_unif.Rdata")


chains.unif <- MCMCchains(unif.out, params=c("Tmin", "Tmax", "c"))
temps <- seq(5, 45, 0.1) 
curves <- apply(chains.unif, 1, function(x) briere(temps, x[1], x[2], x[3]))

meancurve <- apply(curves, 1, mean)
CI <- apply(curves, 1, quantile, c(0.025, 0.975))


plot(data$T, 1 / data$value, xlab="Temperature [°C]",
     ylab="MDR", xlim=c(5, 45), ylim=c(0, 0.16), pch=20, col="steelblue",
     main="Uniform priors")

lines(temps, meancurve, col="black")
lines(temps, CI[1,], col="gray")
lines(temps, CI[2,], col="gray")

## Check prior-posterior overlap
# Generate samples from priors.
prior.Tmin <- runif(125000, 0, 20) # 125000: Number of samples needs to match MCMC chain.
prior.Tmax <- runif(125000, 28, 45)
prior.c <- runif(125000, 0, 0.01)
prior.sigma <- runif(125000, 0, 0.2)

# Plot prior along with posterior.
MCMCtrace(unif.out, params=c("Tmin", "Tmax", "c", "sigma"), 
          priors=cbind(prior.Tmin, prior.Tmax, prior.c, prior.sigma))


##### Run JAGS
inf.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())

inf.out
round(inf.out$BUGSoutput$summary, 6) # Add more digits so we can see results for c.
# Save results for later use.
#save(inf.out, file = "MDR_Cpip_jags_inf.Rdata")
mcmcplot(inf.out)


chains.inf <- MCMCchains(inf.out, params=c("Tmin", "Tmax", "c"))
temps <- seq(5, 45, 0.1) 
curves <- apply(chains.inf, 1, function(x) briere(temps, x[1], x[2], x[3]))

# Find mean curve and credible intervals.
meancurve <- apply(curves, 1, mean)
CI <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data$T, 1 / data$value, xlab="Temperature [°C]",
     ylab="MDR", xlim=c(5, 45), ylim=c(0, 0.16), pch=20, col="steelblue",
     main="Informative priors")

# Plot mean curve and 95% credible interval.
lines(temps, meancurve, col="black")
lines(temps, CI[1,], col="gray")
lines(temps, CI[2,], col="gray")

dim(curves)

## Check prior-posterior overlap
# Generate samples from priors.
prior.Tmin <- rnorm(125000, 10, 5) # 125000: Number of samples in MCMC chain.
prior.Tmax <- rnorm(125000, 35, 5)
prior.c <- runif(125000, 0, 0.001)
prior.sigma <- runif(125000, 0, 0.2)

MCMCtrace(inf.out, params=c("Tmin", "Tmax", "c", "sigma"), 
          priors=cbind(prior.Tmin, prior.Tmax, prior.c, prior.sigma))
