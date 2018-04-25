nyears <- 5
R_max <- log(b,e)
Init_N <- initabund
Init_Nsd <- initabundsd
SD_lambda <- sdb
alist <- list()

Ricker <- function(prev_abund){ 
  x <- prev_abund * e^(rnorm(1,R_max,SD_lambda))*(1-(prev_abund/K))
  return(x)
}

BurnSeedbank <- function(prev_abund,lastyear){
  if(prev_abund > 4000){
    x <- otherwise*lastyear
    return(x)
  }
  x <- max(0,rnorm(1,burning,sdburning))*(reprodensityform[[2]]*log(lastyear)+reprodensityform[[1]])*lastyear
  if(is.nan(x)|is.infinite(x)){x <- 0}
  if(x==0){
    x <- rnorm(1,burning,sdburning)*lastyear*rnorm(1,seedoutput,seedoutputsd)
  }
  return(x)
}

BurnPVA <- function(nreps,nyears,Init_N,R_max){
  
  PopArray2 <- array(0,dim=c((nyears+1),nreps))
  Seeds <- array(0,dim=c((nyears+1),nreps))
  
  for(rep in 1:nreps){
    
    PopArray2[1,rep] <- max(0,rnorm(1,burningabund,burningabundsd))*Init_N
    if(PopArray2[1,rep]<=0){
      PopArray2[1,rep] <- Init_N*burningabund
    }
    Seeds[1,rep] <- rnorm(1,burning,sdburning)*(reprodensityform[[2]]*log(Init_N)+reprodensityform[[1]])*Init_N
    if(Seeds[1,rep]<=0){
      Seeds[1,rep] <- burning*(reprodensityform[[2]]*log(Init_N)+reprodensityform[[1]])*Init_N
    }
    
    for(y in 2:(nyears+1)){
      nextyear <- max(0,Ricker(PopArray2[y-1,rep]))
      PopArray2[y,rep] <- nextyear
      nextyearseed <- max(0,BurnSeedbank(Seeds[y-1,rep],PopArray2[y-1,rep]))
      Seeds[y,rep] <- nextyearseed
    }
  }
  alist$poparray <- PopArray2
  alist$seeds <- Seeds
  return(alist)
}

BurnDefault <- BurnPVA(nreps,nyears,Init_N,R_max)


burnendpop <- mean(BurnDefault$poparray[nrow(BurnDefault$poparray),])
burnendseed <- mean(BurnDefault$seeds[nrow(BurnDefault$seeds),])


aburn <- burnendpop
sburn <- sd(BurnDefault$poparray[nrow(BurnDefault$poparray),])
nBurn <- nreps
errorburn <- qnorm(0.975)*sburn/sqrt(nBurn)
leftburn <- aburn-errorburn
rightburn <- aburn+errorburn
ciburn <- c(leftburn,rightburn)

aburns <- burnendseed
sburns <- sd(BurnDefault$seed[nrow(BurnDefault$seed),])
nBurns <- nreps
errorburns <- qnorm(0.975)*sburns/sqrt(nBurns)
leftburns <- aburns-errorburns
rightburns <- aburns+errorburns
ciburns <- c(leftburns,rightburns)
