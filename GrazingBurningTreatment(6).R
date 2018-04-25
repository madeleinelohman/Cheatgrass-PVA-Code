nyears <- 5
R_max <- log(gb,e)
Init_N <- initabund
Init_Nsd <- initabundsd
SD_lambda <- sdgb
alist <- list()

Ricker <- function(prev_abund){  
  x <- prev_abund * e^((rnorm(1,R_max,SD_lambda))*(1-(prev_abund/K)))
  return(x)
}

GrazeBurnSeedbank <- function(prev_abund,lastyear){
  if(lastyear > 4000){
    x <- otherwise*lastyear
    return(x)
  }
  x <- max(0,rnorm(1,grazingburning,0.4))*(reprodensityform[[2]]*log(lastyear)+reprodensityform[[1]])*lastyear
  if(is.nan(x)|is.infinite(x)){x <- 0}
  if(x==0){
    x <- grazingburning*lastyear*rnorm(1,seedoutput,seedoutputsd)
  }
  return(x)
}

GrazeBurnPVA <- function(nreps,nyears,Init_N,R_max){
  
  PopArray2 <- array(0,dim=c((nyears+1),nreps))
  Seeds <- array(0,dim=c((nyears+1),nreps))
  
  for(rep in 1:nreps){
    
    PopArray2[1,rep] <- max(0,rnorm(1,grazingburningabund,0.4))*Init_N
    if(PopArray2[1,rep]<=0){
      PopArray2[1,rep] <- Init_N*grazingburningabund
    }
     Seeds[1,rep] <- max(0,rnorm(1,grazingburning,0.4))*(reprodensityform[[2]]*log(Init_N)+reprodensityform[[1]])*Init_N
     if(Seeds[1,rep]<=0){
       Seeds[1,rep] <- grazingburning*(reprodensityform[[2]]*log(Init_N)+reprodensityform[[1]])*Init_N
     }
     
    for(y in 2:(nyears+1)){
      nextyear <- max(0,trunc(Ricker(PopArray2[y-1,rep])))
      if(nextyear==0){
        nextyear <- PopArray2[y-1,rep]
      }
      PopArray2[y,rep] <- nextyear
      nextyearseed <- max(0,GrazeBurnSeedbank(Seeds[y-1,rep],PopArray2[y-1,rep]))
      Seeds[y,rep] <- nextyearseed
    }
  }
  alist$poparray <- PopArray2
  alist$seeds <- Seeds
  return(alist)
}

GrazeBurnDefault <- GrazeBurnPVA(nreps,nyears,Init_N,R_max)

fakegb <- matrix(NA,nrow=6,ncol=2)
for(i in 1:nrow(GrazeBurnDefault$poparray)){
  mean <- mean(GrazeBurnDefault$poparray[i,])
  sd <- sd(GrazeBurnDefault$poparray[i,])
  n<- nreps
  error <- qnorm(0.975)*sd/sqrt(n)
  fakegb[i,1] <- mean-error
  fakegb[i,2] <- mean+error
}

gbv <- NULL
for(i in 1:nrow(GrazeBurnDefault$poparray)){
  gbv[i] <- mean(GrazeBurnDefault$poparray[i,])
}
plot(gbv~time,type="p",xaxp=c(1,6,5),xlab="Year",ylab="Plants per meter^2",pch=19,
     col="black",main="Graze/Burn Population Projection",ylim=c(-4,300))
arrows(time,fakegb[,1],time,fakegb[,2],length=0)

new <- GrazeBurnDefault$seeds[6,]>480
newefake <- which(GrazeBurnDefault$poparray[6,]<480)

grazeburnendpop <- mean(GrazeBurnDefault$poparray[nrow(GrazeBurnDefault$poparray),])
grazeburnendseed <- mean(GrazeBurnDefault$seeds[nrow(GrazeBurnDefault$poparray),])

agrazeburn <- grazeburnendpop
sgrazeburn <- sd(GrazeBurnDefault$poparray[nrow(GrazeBurnDefault$poparray),])
ngrazeburn <- nreps
errorgrazeburn <- qnorm(0.975)*sgrazeburn/sqrt(ngrazeburn)
leftgrazeburn <- agrazeburn-errorgrazeburn
rightgrazeburn <- agrazeburn+errorgrazeburn
cigrazeburn <- c(leftgrazeburn,rightgrazeburn)

agrazeburns <- grazeburnendseed
sgrazeburns <- sd(GrazeBurnDefault$seeds[nrow(GrazeBurnDefault$seeds),])
ngrazeburns <- nreps
errorgrazeburns <- qnorm(0.975)*sgrazeburns/sqrt(ngrazeburns)
leftgrazeburns <- agrazeburns-errorgrazeburns
rightgrazeburns <- agrazeburns+errorgrazeburns
cigrazeburns <- c(leftgrazeburns,rightgrazeburns)


