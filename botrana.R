### Code for brief crash course on Bayesian TPC estimation.

# Note: Please change this to your working directory!
setwd("/Users/mauricio/Dropbox/mordecailab/bayesTPC")

# We will need these packages. Please make sure you install them beforehand!
library('R2jags')
library('mcmcplots')
library('MCMCvis')

# Read data.
data <- read.csv("briere_data_l_botrana.csv")
data

plot(data$T, 1 / data$pupae, pch=20, xlab="Temperature [°C]",
     ylab="Developmental Rate (eggs)", col="steelblue")

# We know what temperatures mean. But sometimes our model has parameters that
# don't have a straightforward interpretation. To set a reasonable prior 
# distributions for these parameters, it's a good idea to plot our function.
briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax)  )
  return(result)
}

T <- seq(5, 40, 0.01)

# Pick some ballpark reasonable values of Tmin and Tmax,
Tmin = 8
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

# Simple model with uniform priors.
sink("briere_unif.txt")
cat("
    model{
    
    ## Priors
    c ~ dunif(0, 0.01)
    Tmin ~ dunif(0, 20)
    Tmax ~ dunif(28, 45)
    sigma ~ dunif(0, 0.5)
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
    Tmin ~ dnorm(10, 1 / 25) # Assumes Tmin has a 95% chance to be in interval (0, 20).
    Tmax ~ dnorm(35, 1 / 25) # Assumes Tmax has a 95% chance to be in interval (25, 45).
    sigma ~ dunif(0, 0.5)
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

##### inits Function
inits<-function(){list(
  c = runif(0.00001, 0.001),
  Tmax = runif(1, min=35, max=40),
  Tmin = runif(1, min=5, max=10),
  sigma = runif(1, min=0.01, max=0.15))}

##### Parameters to Estimate
parameters <- c("c", "Tmin", "Tmax","sigma", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Organize Data for JAGS
trait <- 1 / data$pupae
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
##### Run JAGS
pupae.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_unif.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())

round(pupae.out$BUGSoutput$summary, 4)
mcmcplot(pupae.out)

chains.unif <- MCMCchains(pupae.out, params=c("Tmin", "Tmax", "c"))
temps <- seq(5, 45, 0.1) 
curves <- apply(chains.unif, 1, function(x) briere(temps, x[1], x[2], x[3]))

# Get mean curve and 2.5% and 97.5% quantiles (for CI).
meancurve <- apply(curves, 1, mean)
CI <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data$T, 1 / data$pupae, xlab="Temperature [°C]",
     ylab="pupae development rate", xlim=c(5, 45), ylim=c(0, 0.16), pch=20, col="steelblue",
     main="Uniform priors")


# Plot mean prediction and 95% credible interval.
lines(temps, meancurve, col="black")
lines(temps, CI[1,], col="gray")
lines(temps, CI[2,], col="gray")

## Check prior-posterior overlap
# Generate samples from priors.
prior.Tmin <- runif(125000, 0, 20) # 125000: Number of samples needs to match MCMC chain.
prior.Tmax <- runif(125000, 28, 45)
prior.c <- runif(125000, 0, 0.01)
prior.sigma <- runif(125000, 0, 0.5)

# Plot prior along with posterior.
MCMCtrace(pupae.out, params=c("Tmin", "Tmax", "c", "sigma"), 
          priors=cbind(prior.Tmin, prior.Tmax, prior.c, prior.sigma))



##### Run JAGS
pupae.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())

round(pupae.out.inf$BUGSoutput$summary, 4)
mcmcplot(pupae.out.inf)


chains.inf <- MCMCchains(pupae.out.inf, params=c("Tmin", "Tmax", "c"))
temps <- seq(5, 45, 0.1) 
curves <- apply(chains.inf, 1, function(x) briere(temps, x[1], x[2], x[3]))



# Get mean curve and 2.5% and 97.5% quantiles (for CI).
meancurve <- apply(curves, 1, mean)
CI <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data$T, 1 / data$pupae, xlab="Temperature [°C]",
     ylab="pupae development rate", xlim=c(5, 45), ylim=c(0, 0.16), pch=20, col="steelblue",
     main="Informative priors")

# Plot mean prediction and credible interval.
lines(temps, meancurve, col="black")
lines(temps, CI[1,], col="gray")
lines(temps, CI[2,], col="gray")

## Check prior-posterior overlap
# Generate samples from priors.
prior.Tmin <- rnorm(125000, 10, 5) # 125000: Number of samples in MCMC chain.
prior.Tmax <- rnorm(125000, 35, 5)
prior.c <- runif(125000, 0, 0.001)
prior.sigma <- runif(125000, 0, 0.2)

MCMCtrace(pupae.out.inf, params=c("Tmin", "Tmax", "c", "sigma"), 
          priors=cbind(prior.Tmin, prior.Tmax, prior.c, prior.sigma))

