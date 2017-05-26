###################################################################################

#   Stage-based N-mixture model for white-tailed deer based on model by 
#   Zipkin 2014 Ecology 95:22-29. 3 stages include fawn, does, and bucks. Dynamics
#   are modled based on an autoregressive model (i.e. based on previous year 
#   abundance). I will use Poisson, zero-inflated Poisson, and Negative Binomial 
#   model for abundance and use RJMCMC to determine which model is best. Also, check
#   for collinearity among variables and either exlude or use PCA. 

###################################################################################

###################################################

#         Make fake data                          #

###################################################


data.fn <- function (nsites=20, nprim=5, nocc=20, nstages=3, b0.lam=c(0.9,1.6,1.6), 
                     b0.gam=c(0.1,-2.5,-2.5), b0.ome=c(.3,.9,.8), b0.p=c(.8,.8,.8)){
  
  # Definitions
  
    # nsites = number of sites
    # nprim = number of primary periods
    # nocc = number of secondary periods
    # nstages = number of stages; fawn, doe, buck
    # b0.lam = intercept for initial abundance for fawns, does, bucks
    # b0.gam = intercept for recruitment for fawns, does, bucks
    # b0.ome = intercept for apparent survival for fawns, does, bucks
    # b0.p = intercept for detection probability for fawns, does, bucks
  
  # Arrays
  
    N <- array(dim=c(nsites, nstages, nprim))
    y <- array(dim=c(nsites, nocc, nprim, nstages))
    S <- array(dim=c(nsites, nstages, nprim))
    G <- array(dim=c(nsites, nstages, nprim))
    p <- array(dim=c(nstages))
    
    
  # Initial abundance for primary period 1
    
    lam.f <- exp(b0.lam[1])
    lam.d <- exp(b0.lam[2])
    lam.b <- exp(b0.lam[3])
    
    for(i in 1:nsites){
      N[i,1,1] <- rpois(1,lam.f) 
      N[i,2,1] <- rpois(1,lam.d)
      N[i,3,1] <- rpois(1,lam.b)
    }
    
  # Apparent survival and recruitment
    
    ome.f <- 1 / (1 + exp(-b0.ome[1]))
    ome.d <- 1 / (1 + exp(-b0.ome[2]))
    ome.b <- 1 / (1 + exp(-b0.ome[3]))
    
    gam.f <-exp(b0.gam[1])
    gam.d <-exp(b0.gam[2])
    gam.b <-exp(b0.gam[3])
    
    for(i in 1:nsites){
      for(t in 2){
        S[i,1,t] <- rbinom(1, N[i,1,t-1], ome.f)
        S[i,2,t] <- rbinom(1, N[i,2,t-1], ome.d)
        S[i,3,t] <- rbinom(1, N[i,3,t-1], ome.b)
        
        G[i,1,t] <- rpois(1,gam.f * N[i,2,t-1])
        G[i,2,t] <- rpois(1,gam.d * sum(N[i,,t-1]))
        G[i,3,t] <- rpois(1,gam.b * sum(N[i,,t-1]))
        
        N[i,1,t] <- G[i,1,t]
        N[i,2,t] <- round(0.5 * S[i,1,t]) + S[i,2,t] + G[i,2,t]
        N[i,3,t] <- round(0.5 * S[i,1,t]) + S[i,3,t] + G[i,3,t]
      }
      
      for(t in 3){
        S[i,1,t] <- rbinom(1, N[i,1,t-1], ome.f)
        S[i,2,t] <- rbinom(1, N[i,2,t-1], ome.d)
        S[i,3,t] <- rbinom(1, N[i,3,t-1], ome.b)
        
        G[i,1,t] <- rpois(1,gam.f * N[i,2,t-1])
        G[i,2,t] <- rpois(1,gam.d * sum(N[i,,t-1]))
        G[i,3,t] <- rpois(1,gam.b * sum(N[i,,t-1]))
        
        N[i,1,t] <- S[i,1,t]
        N[i,2,t] <- S[i,2,t] + G[i,2,t]
        N[i,3,t] <- S[i,3,t] + G[i,3,t]
      }
      
      for(t in 4){
        S[i,1,t] <- rbinom(1, N[i,1,t-1], ome.f)
        S[i,2,t] <- rbinom(1, N[i,2,t-1], ome.d)
        S[i,3,t] <- rbinom(1, N[i,3,t-1], ome.b)
        
        G[i,1,t] <- rpois(1,gam.f * N[i,2,t-1])
        G[i,2,t] <- rpois(1,gam.d * sum(N[i,,t-1]))
        G[i,3,t] <- rpois(1,gam.b * sum(N[i,,t-1]))
        
        N[i,1,t] <- G[i,1,t]
        N[i,2,t] <- round(0.5 * S[i,1,t]) + S[i,2,t] + G[i,2,t]
        N[i,3,t] <- round(0.5 * S[i,1,t]) + S[i,3,t] + G[i,3,t]
      }
      
      for(t in 5){
        S[i,1,t] <- rbinom(1, N[i,1,t-1], ome.f)
        S[i,2,t] <- rbinom(1, N[i,2,t-1], ome.d)
        S[i,3,t] <- rbinom(1, N[i,3,t-1], ome.b)
        
        G[i,1,t] <- rpois(1,gam.f * N[i,2,t-1])
        G[i,2,t] <- rpois(1,gam.d * sum(N[i,,t-1]))
        G[i,3,t] <- rpois(1,gam.b * sum(N[i,,t-1]))
        
        N[i,1,t] <- S[i,1,t]
        N[i,2,t] <- S[i,2,t] + G[i,2,t]
        N[i,3,t] <- S[i,3,t] + G[i,3,t]
      }
    }
    
    for(s in 1:nstages){
      p[s] <- 1 / (1 + exp(-b0.p[s]))

      for(t in 1:nprim){
        for(j in 1:nocc){
          for(i in 1:nsites){
            y[i,j,t,s] <- rbinom(1, N[i,s,t], p[s])
          }
        }
      }
    }
    
    
    return(list(nsites=nsites, nstages=nstages, nocc=nocc, nprim=nprim, y=y,
                N=N, b0.gam=b0.gam, b0.p=b0.p, b0.lam=b0.lam, b0.ome=b0.ome, 
                G=G, S=S))
    
}


datum <- data.fn(nsites=20, nprim=5, nocc=20, nstages=3, b0.lam=c(0.9,1.6,1.6), 
                 b0.gam=c(0.1,-2.5,-2.5), b0.ome=c(.3,.9,.8), b0.p=c(.8,.8,.8))

attach(datum)

###################################################

#         Run the model                           #

###################################################


library(R2jags)
library(mcmcplots)


sink("Pmodel.txt")
cat("
    model{
    
      # Specify priors for all parameters in the model

        # Intercept for initial abundance (lam), per-capita recruitment (gam), 
        # apparent survival (ome), and detection (p) for stage s

        for(s in 1:nstages){
          b0.lam[s] ~ dnorm(0,0.001)
          b0.gam[s] ~ dnorm(0,0.001)
          b0.ome[s] ~ dnorm(0,0.001)
          b0.p[s] ~ dnorm(0,0.001)
        }
  
      # Ecological process model
        
        # Initial abundance for primary period 1
          
          # Apparent survival and recruitment
          
              logit(omega[1]) <- b0.ome[1]
              logit(omega[2]) <- b0.ome[2]
              logit(omega[3]) <- b0.ome[3]
              
              log(gamma[1]) <- b0.gam[1]
              log(gamma[2]) <- b0.gam[2]
              log(gamma[3]) <- b0.gam[3]

        for(i in 1:nsites){
          log(lambda[i,1]) <- b0.lam[1]
          log(lambda[i,2]) <- b0.lam[2]
          log(lambda[i,3]) <- b0.lam[3]

          N[i,1,1] ~ dpois(lambda[i,1])
          N[i,2,1] ~ dpois(lambda[i,2])
          N[i,3,1] ~ dpois(lambda[i,3])


          # Abundance for primary periods 2 through number of primary periods

          for(t in 2:nprim){
            S[i,1,t] ~ dbin(omega[1], N[i,1,t-1])
            S[i,2,t] ~ dbin(omega[2], N[i,2,t-1])
            S[i,3,t] ~ dbin(omega[3], N[i,3,t-1])

            G[i,1,t] ~ dpois(gamma[1] * N[i,2,t-1])
            G[i,2,t] ~ dpois(gamma[2] * sum(N[i,,t-1]))
            G[i,3,t] ~ dpois(gamma[3] * sum(N[i,,t-1]))

          }
        
          # Abundance for fall periods (2 and 4)
          for(t in c(2,4)){
            N[i,1,t] <- G[i,1,t]
            N[i,2,t] <- rbin0.5 * S[i,1,t] + S[i,2,t] + G[i,2,t]
            N[i,3,t] <- 0.5 * S[i,1,t] + S[i,3,t] + G[i,3,t]
          }
          
          # Abundance for spring periods (3 and 5)
          for(t in c(3,5)){
            N[i,1,t] <- S[i,1,t]
            N[i,2,t] <- S[i,2,t] + G[i,2,t]
            N[i,3,t] <- S[i,3,t] + G[i,3,t]
          }
       
        }

        # Observation model

        for(s in 1:nstages){
          logit(p[s]) <- b0.p[s]
          
          for(t in 1:nprim){
            for(j in 1:nocc){
              for(i in 1:nsites){
                y[i,j,t,s] ~ dbin(p[s], N[i,s,t])
              }
            }
          }
        }



    }",fill=TRUE)
sink()

win.data<-list(y=y, nsites=nsites, nstages=nstages, nocc=nocc, nprim=nprim)


Nst <- apply(y, c(1,4,3), max, na.rm=TRUE)
Nst[is.na(Nst)]<-1

inits<-function(){list(N=Nst[,,1], b0.lam=rnorm(1,0,1), b0.gam=rnorm(1,0,1), 
                       b0.p=rnorm(1,0,1), b0.ome=rnorm(1,0,1))}



params<-c("b0.lam", "b0.gam", "b0.p", "b0.ome", "N", "S", "G")

ni<-1000
nt<-1
nb<-50
nc<-3

out.P<-jags(data=win.data,inits=inits, params,model.file="Pmodel.txt",n.thin=nt, 
            n.chains=nc,n.burnin=nb, n.iter=ni)










































###################################################################################

# Bayesian code for the Dail-Madsen open population model. Before begininng, 
# we tested for collinearity between covariates, and covariates which were highly 
# correlated were not used in the same model. Due to collinearty between the habitat
# covariates, we used a PCA and used the principle components instead of the actual 
# variables. We then used the program Unmarked to determine which process model to 
# use (Poisson or zero-inflated Poisson) and which dynamics to use. 

###################################################################################

# Get the number of sites and number of survey occasions from the data
M <- nrow(Models) #Number of sites
T <- 5 #Number of primary periods. Must have equal number of secondary periods in each primary period. For missing values, put blanks
JT <- ncol(Models)
J <- JT/T #The number of secondary survey occasions within each primary period


# Start the exponential growth model with a Poisson distribution

sink("model.txt")
cat("
    model{

      # Define priors for all parameters
        B0.lam~dnorm(0, 0.1) #Intercept for poisson log-linear abundance model
        B0.r~dnorm(0, 0.1) #Intercept for intrinsic per capita rate of increase
        B0.p~dunif(0, 1) #Intercept for detection probabilty on the logit scale
      
      #The likelihood. Loop over M sites
        for(i in 1:M){
        loglam[i,1]<- B0.lam #The Poisson generalized linear model for initial abundance
        lam[i,1]<-exp(loglam[i,1]) #Expected initial abundance at each site
        N[i,1]~dpois(lam[i,1]) #Assumes actual counts are poisson distributed 
        
      #Now loop over secondary periods within primary period 1 for Binomial observation process
      for(j in 1:J){
      logitp[i,1,j]<-B0.p #The binomial detection model for year 1
      p[i,1,j]<-1/(1+exp(-logitp[i,1,j])) #The inverse logit transformation
      y[i,1,j]~dbin(p[i,1,j], N[i,1])
      }
      
      #Loop over the remaining years
      for(t in 2:T){
      logr[i,t-1]<-B0.r #The Poisson model for recruitment rate
      r[i,t-1]<-exp(logr[i,t-1]) 
      N[i,t]~dpois(r[i,t-1]*N[i,t-1])
      
      for(j in 1:J){
      logitp[i,t,j]<-B0.p #The binomial detection model for all years
      p[i,t,j]<-1/(1+exp(-logitp[i,t,j])) #The inverse logit transformation
      y[i,t,j]~dbin(p[i,t,j], N[i,t])
      eval[i,t,j]<-p[i,t,j]*N[i,t]
      E[i,t,j]<-pow((y[i,t,j]-eval[i,t,j]),2)/(eval[i,t,j]+0.5)
      
      y.new[i,t,j]~dbin(p[i,t,j],N[i,t])
      E.new[i,t,j]<-pow((y.new[i,t,j]-eval[i,t,j]),2)/(eval[i,t,j]+0.5)
      }
      }
      }
      for(t in 1:T){
      Ntot[t]<-sum(N[,t])
      }
      fit<-sum(E[,,])
      fit.new<-sum(E.new[,,])
    }",fill=TRUE)
sink()

win.data<-list(M=26, T=5, J=J, y=matrix(unlist((Models)),nrow=26))

Nst<-apply(Models, c(1,3), max)+1
Nst[is.na(Nst)]<-1
inits<-function(){list(N=rpois(20,25), B0.lam=rnorm(1,0,1), B0.r=rnorm(1,0,1), B0.p=rnorm(1,0,1))}

params<-c("B0.lam", "B0.r", "B0.p", "Ntot", "fit", "fit.new")

ni<-1000
nt<-1
nb<-50
nc<-3

out.P<-bugs(data=win.data,inits=inits, params,model.file="P_model.txt",n.thin=nt, n.chains=nc,n.burnin=nb, n.iter=ni,debug=TRUE,DIC=TRUE,working.directory=getwd())



#Start the exponential growth model with a zero-inflated Poisson distribution
sink("model.txt")
cat("
    model{
    #Define priors for all parameters
    B0.lam~dnorm(0, 0.1) #Intercept for poisson log-linear abundance model
    B0.psi~dunif(0,1) #Intercept for the zero-inflation parameter, i.e. occupancy
    B0.r~dnorm(0, 0.1) #Intercept for intrinsic per capita rate of increase
    B0.p~dunif(0, 1) #Intercept for detection probabilty on the logit scale
    
    #The likelihood. Loop over M sites
    for(i in 1:M){
    logitpsi[i]<-B0.psi #Binomial model for occupancy
    psi[i]<-1/(1+exp(-logitpsi[i]))
    occu[i]~dbern(1-psi[i])
    loglam[i,1]<- B0.lam #The Poisson generalized linear model for initial abundance
    lam[i,1]<-exp(loglam[i,1]) #Expected initial abundance at each site
    N[i,1]~dpois(lam[i,1]*occu[i]) #Assumes actual counts are poisson distributed 
    
    #Now loop over secondary periods within primary period 1 for Binomial observation process
    for(j in 1:J){
    logitp[i,1,j]<-B0.p #The binomial detection model for year 1
    p[i,1,j]<-1/(1+exp(-logitp[i,1,j])) #The inverse logit transformation
    y[i,1,j]~dbin(p[i,1,j], N[i,1])
    }
    
    #Loop over the remaining years
    for(t in 2:T){
    logr[i,t-1]<-B0.r #The Poisson model for recruitment rate
    r[i,t-1]<-exp(logr[i,t-1]) 
    N[i,t]~dpois(r[i,t-1]*N[i,t-1]*occu[i])
    
    for(j in 1:J){
    logitp[i,t,j]<-B0.p #The binomial detection model for all years
    p[i,t,j]<-1/(1+exp(-logitp[i,t,j])) #The inverse logit transformation
    y[i,t,j]~dbin(p[i,t,j], N[i,t])
    }
    }
    }
    for(t in 1:T){
    Ntot[t]<-sum(N[,t])
    }
    }",fill=TRUE)
sink()



#Start the exponential growth model with a negative binomial distribution for initial abundance
sink("model.txt")
cat("
    model{
    #Define priors for all parameters
    B0.lam~dnorm(0, 0.1) #Intercept for poisson log-linear abundance model
    alpha~dunif(0,400)
    B0.r~dnorm(0, 0.1) #Intercept for intrinsic per capita rate of increase
    B0.p~dunif(0, 1) #Intercept for detection probabilty on the logit scale
    
    #The likelihood. Loop over M sites
    for(i in 1:M){
    loglam[i,1]<- B0.lam #The Poisson generalized linear model for initial abundance
    lam[i,1]<-exp(loglam[i,1]) #Expected initial abundance at each site
    P.NB[i,1]<-alpha/(alpha+lam[i,1])
    N[i,1]~dnegbin(P[i,1], alpha) #Assumes actual counts are poisson distributed 
    
    #Now loop over secondary periods within primary period 1 for Binomial observation process
    for(j in 1:J){
    logitp[i,1,j]<-B0.p #The binomial detection model for year 1
    p[i,1,j]<-1/(1+exp(-logitp[i,1,j])) #The inverse logit transformation
    y[i,1,j]~dbin(p[i,1,j], N[i,1])
    }
    
    #Loop over the remaining years
    for(t in 2:T){
    logr[i,t-1]<-B0.r #The Poisson model for recruitment rate
    r[i,t-1]<-exp(logr[i,t-1]) 
    N[i,t]~dpois(r[i,t-1]*N[i,t-1])
    
    for(j in 1:J){
    logitp[i,t,j]<-B0.p #The binomial detection model for all years
    p[i,t,j]<-1/(1+exp(-logitp[i,t,j])) #The inverse logit transformation
    y[i,t,j]~dbin(p[i,t,j], N[i,t])
    }
    }
    }
    for(t in 1:T){
    Ntot[t]<-sum(N[,t])
    }
    }",fill=TRUE)
sink()





#Start the auto-regressive model 
sink("model.txt")
cat("
    model{
    #Define priors for all parameters
    B0.lam~dnorm(0, 0.1) #Intercept for poisson log-linear abundance model
    B0.gam~dnorm(0, 0.1) #Intercept for gamma
    B0.om~dunif(0, 1) #Intercept for omega
    B0.p~dunif(0, 1) #Intercept for detection probabilty on the logit scale
    
    #The likelihood. Loop over M sites
    for(i in 1:M){
    loglam[i,1]<- B0.lam #The Poisson generalized linear model for initial abundance
    lam[i,1]<-exp(loglam[i,1]) #Expected initial abundance at each site
    N[i,1]~dpois(lam[i,1]) #Assumes actual counts are poisson distributed 
    
    #Now loop over secondary periods within primary period 1 for Binomial observation process
    for(j in 1:J){
    logitp[i,1,j]<-B0.p #The binomial detection model for year 1
    p[i,1,j]<-1/(1+exp(-logitp[i,1,j])) #The inverse logit transformation
    y[i,1,j]~dbin(p[i,1,j], N[i,1])
    }
    #Loop over the remaining years
    for(t in 2:T){
    logitom[i,t-1]<-B0.om #The binomial apparent survival model between each primary period
    omega[i,t-1]<-1/(1+exp(-logitom[i,t-1])) #Inverse logit transformation
    S[i,t-1]~dbin(omega[i,t-1], N[i,t-1]) #Apparent survival
    
    loggam[i,t-1]<-B0.gam #The Poisson model for recruitment rate
    gamma[i,t-1]<-exp(loggam[i,t-1])
    G[i,t-1]~dpois(gamma[i,t-1]*N[i,t-1])
    N[i,t]<-S[i,t-1]+G[i,t-1]
    
    for(j in 1:J){
    logitp[i,t,j]<-B0.p #The binomial detection model for all years
    p[i,t,j]<-1/(1+exp(-logitp[i,t,j])) #The inverse logit transformation
    y[i,t,j]~dbin(p[i,t,j], N[i,t])
    }
    }
    }
    for(t in 1:T){
    Ntot[t]<-sum(N[,t])
    }
    }",fill=TRUE)
sink()



#Start the exponential growth model 
sink("model.txt")
cat("
    model{
    #Define priors for all parameters
    B0.lam~dnorm(0, 0.1) #Intercept for poisson log-linear abundance model
    B0.r~dnorm(0, 0.1) #Intercept for intrinsic per capita rate of increase
    B0.p~dunif(0, 1) #Intercept for detection probabilty on the logit scale
    
    #The likelihood. Loop over M sites
    for(i in 1:M){
    loglam[i,1]<- B0.lam #The Poisson generalized linear model for initial abundance
    lam[i,1]<-exp(loglam[i,1]) #Expected initial abundance at each site
    N[i,1]~dpois(lam[i,1]) #Assumes actual counts are poisson distributed 
    
    #Now loop over secondary periods within primary period 1 for Binomial observation process
    for(j in 1:J){
    logitp[i,1,j]<-B0.p #The binomial detection model for year 1
    p[i,1,j]<-1/(1+exp(-logitp[i,1,j])) #The inverse logit transformation
    y[i,1,j]~dbin(p[i,1,j], N[i,1])
    }
    
    #Loop over the remaining years
    for(t in 2:T){
    logr[i,t-1]<-B0.r #The Poisson model for recruitment rate
    r[i,t-1]<-exp(logr[i,t-1]) 
    N[i,t]~dpois(r[i,t-1]*N[i,t-1])
    
    for(j in 1:J){
    logitp[i,t,j]<-B0.p #The binomial detection model for all years
    p[i,t,j]<-1/(1+exp(-logitp[i,t,j])) #The inverse logit transformation
    y[i,t,j]~dbin(p[i,t,j], N[i,t])
    }
    }
    }
    for(t in 1:T){
    Ntot[t]<-sum(N[,t])
    }
    }",fill=TRUE)
sink()




#Start the Ricker-logistic model 
sink("model.txt")
cat("
    model{
    #Define priors for all parameters
    B0.lam~dnorm(0, 0.1) #Intercept for poisson log-linear abundance model
    B0.r~dnorm(0, 0.1) #Intercept for instantaneous population growth rate
    K~dunif(0, 500) #stable equilibrium
    B0.p~dunif(0, 1) #Intercept for detection probabilty on the logit scale
    
    #The likelihood. Loop over M sites
    for(i in 1:M){
    loglam[i,1]<- B0.lam #The Poisson generalized linear model for initial abundance
    lam[i,1]<-exp(loglam[i,1]) #Expected initial abundance at each site
    N[i,1]~dpois(lam[i,1]) #Assumes actual counts are poisson distributed 
    
    #Now loop over secondary periods within primary period 1 for Binomial observation process
    for(j in 1:J){
    logitp[i,1,j]<-B0.p #The binomial detection model for year 1
    p[i,1,j]<-1/(1+exp(-logitp[i,1,j])) #The inverse logit transformation
    y[i,1,j]~dbin(p[i,1,j], N[i,1])
    }
    
    #Loop over the remaining years
    for(t in 2:T){
    logr[i,t-1]<-B0.r #The Poisson model for Instantaneous growth rate
    r[i,t-1]<-exp(logr[i,t-1])
    muN[i,t-1]<-N[i,t-1]*exp(r[i,t-1]*(1-N[i,t-1]/K))
    N[i,t]~dpois(muN[i,t-1])
    
    for(j in 1:J){
    logitp[i,t,j]<-B0.p #The binomial detection model for all years
    p[i,t,j]<-1/(1+exp(-logitp[i,t,j])) #The inverse logit transformation
    y[i,t,j]~dbin(p[i,t,j], N[i,t])
    }
    }
    }
    for(t in 1:T){
    Ntot[t]<-sum(N[,t])
    }
    }",fill=TRUE)
sink()




#Start the Gompertz-logistic model 
sink("model.txt")
cat("
    model{
    #Define priors for all parameters
    B0.lam~dnorm(0, 0.1) #Intercept for poisson log-linear abundance model
    B0.r~dnorm(0, 0.1) #Intercept for instantaneous population growth rate
    K~dunif(0, 500) #stable equilibrium
    B0.p~dunif(0, 1) #Intercept for detection probabilty on the logit scale
    
    #The likelihood. Loop over M sites
    for(i in 1:M){
    loglam[i,1]<- B0.lam #The Poisson generalized linear model for initial abundance
    lam[i,1]<-exp(loglam[i,1]) #Expected initial abundance at each site
    N[i,1]~dpois(lam[i,1]) #Assumes actual counts are poisson distributed 
    
    #Now loop over secondary periods within primary period 1 for Binomial observation process
    for(j in 1:J){
    logitp[i,1,j]<-B0.p #The binomial detection model for year 1
    p[i,1,j]<-1/(1+exp(-logitp[i,1,j])) #The inverse logit transformation
    y[i,1,j]~dbin(p[i,1,j], N[i,1])
    }
    
    #Loop over the remaining years
    for(t in 2:T){
    logr[i,t-1]<-B0.r #The Poisson model for Instantaneous growth rate
    r[i,t-1]<-exp(logr[i,t-1])
    muN[i,t-1]<-N[i,t-1]*exp(r[i,t-1]*(1-log(N[i,t-1]+1)/log(K+1)))
    N[i,t]~dpois(muN[i,t-1])
    
    for(j in 1:J){
    logitp[i,t,j]<-B0.p #The binomial detection model for all years
    p[i,t,j]<-1/(1+exp(-logitp[i,t,j])) #The inverse logit transformation
    y[i,t,j]~dbin(p[i,t,j], N[i,t])
    }
    }
    }
    for(t in 1:T){
    Ntot[t]<-sum(N[,t])
    }
    }",fill=TRUE)
sink()
