library("popbio")

#projection matrices based on values obtained in the literature; lambdas calculated and, when possible, 
#compared to lambdas given in the literature to determine best values for survival and seed output
proj1 <- matrix(c(
  0,27.93,
  0.2,0
),nrow=2,ncol=2,byrow=TRUE)
#lambda1 <- lambda(proj1)

proj2 <- matrix(c(
  0,10.74,
  0.39,0
),nrow=2,ncol=2,byrow=TRUE)
#lambda2 <- lambda(proj2)

proj3 <- matrix(c(
  0,10,
  0.555,0
),nrow=2,ncol=2,byrow=TRUE)
#lambda3 <- lambda(proj3)

proj4 <- matrix(c(
  0,10.7,
  0.39,0
),nrow=2,ncol=2,byrow=TRUE)
#lambda4 <- lambda(proj4)

#obtains the mean and standard deviation of the projections matrices and lambdas
#alllambda <- c(lambda1,lambda2,lambda3,lambda4)
#sdlambda <- sd(alllambda)
#meanlambda <- mean(alllambda)

projections <- array(c(proj1,proj2,proj3,proj4), dim=c(2,2,4)) 
sdmatrix <- apply(projections,c(1,2),sd,na.rm=T)
meanmatrix <- apply(projections,c(1,2),mean,na.rm=T)

#numbers obtained from the literature
#-seed refers to the proportion of seeds left after treatment
#-pop refers to the proportion of population left after treatment
grazeseed <- c(0.328,0.12)
burnseed <- c(1.12,0.2)
grazeburnseed <- 0.1
grazepop <- c(0.55,0.41)
burnpop <- c(0.58,0.8)
grazeburnpop <- 0.1

x <- data.frame(x=c(660,4000,178,90,4850),y=c(10.7,1.75,10,28,1.61),x.log=log(c(660,4000,178,90,4850)))
reprodensity <- lm(x$y~log(x$x))
reprodensityform <- reprodensity$coefficients

initabund <- 903
initabundsd <- 527
K <- 4000
nreps <- 1500
e <- exp(1)

mylist <- list()

treatment <- function(treatment){
  
  matrixfake <- matrix(nrow = 2,ncol = 2)
  fakematrix <- matrix(nrow = 2,ncol = 2)
  #from a normal distribution (environmental stochasticity), determining a base seed output and survival 
  #values based on the matrix created from the literature
  for(rep in 1:nreps){
    mylist$reprohold[[rep]] <- rnorm(1,meanmatrix[3],sdmatrix[3]) 
    mylist$survhold[[rep]] <- rnorm(1,meanmatrix[2],sdmatrix[2])
    fakematrix <- matrix(c(0,rnorm(1,meanmatrix[2],sdmatrix[2]),rnorm(1,meanmatrix[3],sdmatrix[3]),0),
                         nrow=2,ncol=2)
    mylist$controllambda[[rep]] <- lambda(fakematrix)
  }
  
  repro <- mean(mylist$reprohold)
  surv <- mean(mylist$survhold)
  
  mylist$lambda[[1]] <- mean(mylist$controllambda)
  #initial abundance before treatment
  mylist$pop[[1]] <- initabund
  #calculated value of seed output before treatment
  mylist$repro[[1]] <- repro
  
  #depending on treatment specified, creating population values of abundance the first year after treatment
  #and seed output during treatment
  if(treatment=="graze"){
    for(rep in 2:nreps){
      mylist$pop[[rep]] <- rnorm(1,mean(grazepop),sd(grazepop))
      mylist$repro[[rep]] <- rnorm(1,mean(grazeseed),sd(grazeseed))
      matrixfake <- matrix(c(0,rnorm(1,meanmatrix[2],sdmatrix[2]),
                             mylist$repro[[1]]*rnorm(1,mean(grazeseed),sd(grazeseed)),0),nrow=2,ncol=2)
      mylist$lambda[[rep]] <- lambda(matrixfake)
    }
    return(mylist)
  }
  if(treatment=="burn"){
    for(rep in 2:nreps){
      mylist$pop[[rep]] <- rnorm(1,mean(burnpop),0.9)
      mylist$repro[[rep]] <- rnorm(1,mean(burnseed),0.9)
      matrixfake <- matrix(c(0,rnorm(1,meanmatrix[2],sdmatrix[2]),
                             mylist$repro[[1]]*rnorm(1,mean(burnseed),0.9),0),nrow=2,ncol=2)
      mylist$lambda[[rep]] <- lambda(matrixfake)
    }
    return(mylist)
  }
  if(treatment=="grazeburn"){
    for(rep in 2:nreps){
      mylist$pop[[rep]] <- rnorm(1,mean(grazeburnpop),0.108)
      mylist$repro[[rep]] <- rnorm(1,mean(grazeburnseed),0.108)
      matrixfake <- matrix(c(0,rnorm(1,meanmatrix[2],sdmatrix[2]),
                             mylist$repro[[1]]*rnorm(1,mean(grazeburnseed),0.15),0),nrow=2,ncol=2)
      mylist$lambda[[rep]] <- lambda(matrixfake)
    }
    return(mylist)
  }
}

#for each treatment, running the function, deriving the control matrix for the simulation and its lambda, 
#putting the seed output into the population matrix, and deriving the initial abundance after treatment
grazingfake <- treatment("graze")
grazingcontrol <- mean(grazingfake$controllambda)
grazing <- mean(grazingfake$repro[2:length(grazingfake$repro)])
sdgrazing <- sd(grazingfake$repro[2:length(grazingfake$repro)])
grazingabund <- mean(grazingfake$pop[2:length(grazingfake$pop)])
grazingabundsd <- sd(grazingfake$pop[2:length(grazingfake$pop)])
g <- mean(grazingfake$lambda)
sdg <- sd(grazingfake$lambda)

burningfake <- treatment("burn")
burningcontrol <- mean(burningfake$controllambda)
burning <- mean(burningfake$repro[2:length(burningfake$repro)])
sdburning <- sd(burningfake$repro[2:length(burningfake$repro)])
burningabund <- mean(burningfake$pop[2:length(burningfake$pop)])
burningabundsd <- sd(burningfake$pop[2:length(burningfake$pop)])
b <- mean(burningfake$lambda)
sdb <- sd(burningfake$lambda)

grazingburningfake <- treatment("grazeburn")
grazingburningcontrol <- mean(grazingburningfake$controllambda)
grazingburning <- mean(grazingburningfake$repro[2:length(grazingburningfake$repro)])
sdgrazingburning <- sd(grazingburningfake$repro[2:length(grazingburningfake$repro)])
grazingburningabund <- mean(grazingburningfake$pop[2:length(grazingburningfake$pop)])
grazingburningabundsd <- sd(grazingburningfake$pop[2:length(grazingburningfake$pop)])
gb <- mean(grazingburningfake$lambda)
sdgb <- sd(grazingburningfake$lambda)


#obtains the mean, standard deviation, and lambda of all control matrices from the simulations put together
allcontrol <- c(grazingcontrol,burningcontrol,grazingburningcontrol)
sdcontrol <- sd(allcontrol)
meancontrol <- mean(allcontrol)
seeds <- c(mean(grazingfake$reprohold),mean(burningfake$reprohold),mean(grazingburningfake$reprohold))
seedoutput <- mean(seeds)
seedoutputsd <- sd(seeds)

