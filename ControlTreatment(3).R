#add one for years of simulation run
nyears <- 5
#r derived from lambda
R_max <- log(meancontrol,e)
#Initial abundance based on mean values defined in the literature
Init_N <- initabund
Init_Nsd <- initabundsd
#the standard deviation of the matrix
SD_lambda <- sdcontrol    

alist <- list()

#the Ricker formula to determine population based on carrying capacity
Ricker <- function(prev_abund){  
  x <- prev_abund * e^((rnorm(1,R_max,SD_lambda))*(1-(prev_abund/K)))
  if(x==0){return(prev_abund)}
  return(x)
}

#to derive the seedbank each year
Seedbank <- function(prev_abund,lastyear){
  #the last year's population multiplyed by seed output per plant and the population's carrying capacity
  if(lastyear > 4000){
    x <- otherwise*lastyear
    return(x)
  }
  x <- (reprodensityform[[2]]*log(lastyear)+reprodensityform[[1]])*lastyear
  if(is.nan(x)|is.infinite(x)){x <- 0}
  if(x==0){
    x <- otherwise*lastyear
  }
  return(x)
}

#PVA model
PVA <- function(nreps,nyears,Init_N,R_max){
  
  PopArray2 <- array(0,dim=c((nyears+1),nreps))
  Seeds <- array(0,dim=c((nyears+1),nreps))
  
  for(rep in 1:nreps){
    
    PopArray2[1,rep] <- max(0,rnorm(1,Init_N,Init_Nsd))
    Seeds[1,rep] <-(reprodensityform[[2]]*log(Init_N)+reprodensityform[[1]])*Init_N
    
    for(y in 2:(nyears+1)){
      nextyear <- max(0,Ricker(PopArray2[y-1,rep]))
      PopArray2[y,rep] <- nextyear 
      nextyearseed <- max(0,Seedbank(Seeds[y-1,rep],PopArray2[y-1,rep]))
      Seeds[y,rep] <- nextyearseed
    }
  }
  alist$poparray <- PopArray2
  alist$seeds <- Seeds
  return(alist)
}

EndControl <- PVA(nreps,nyears,Init_N,R_max)


fakecontrol <- matrix(NA,nrow=6,ncol=2)
for(i in 1:nrow(EndControl$poparray)){
  mean <- mean(EndControl$poparray[i,])
  sd <- sd(EndControl$poparray[i,])
  n <- nreps
  error <- qnorm(0.975)*sd/sqrt(n)
  fakecontrol[i,1] <- mean-error
  fakecontrol[i,2] <- mean+error
}

controlnew <- NULL
for(i in 1:nrow(EndControl$poparray)){
  controlnew[i] <- mean(EndControl$poparray[i,])
}
time <- 1:6
plot(controlnew~time,ylim=c(700,4500),xaxp=c(1,6,5),xlab="Year",ylab=expression("Plants per"~m^{2}),pch=19,
     col="black",main="Control Population Projection",cex=1)
arrows(time,fakecontrol[,1],time,fakecontrol[,2],length=0)

endpop <- mean(EndControl$poparray[nrow(EndControl$poparray),])
endseed <- mean(EndControl$seeds[nrow(EndControl$seeds),])

acontrol <- endpop
scontrol <- sd(EndControl$poparray[nrow(EndControl$poparray),])
ncontrol <- nreps
errorcontrol <- qnorm(0.975)*scontrol/sqrt(ncontrol)
leftcontrol <- acontrol-errorcontrol
rightcontrol <- acontrol+errorcontrol
endci <- c(leftcontrol,rightcontrol)

acontrols <- endseed
scontrols <- sd(EndControl$seeds[nrow(EndControl$seeds),])
ncontrols <- nreps
errorcontrols <- qnorm(0.975)*scontrols/sqrt(ncontrols)
leftcontrols <- acontrols-errorcontrols
rightcontrols <- acontrols+errorcontrols
endcis <- c(leftcontrols,rightcontrols) 




