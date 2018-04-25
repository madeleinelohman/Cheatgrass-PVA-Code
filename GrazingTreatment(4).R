nyears <- 5
R_max <- log(g,e)     
Init_N <- 903
Init_Nsd <- initabundsd
SD_lambda <- sdg
alist <- list()

Ricker <- function(prev_abund){  
  prev_abund * e^((rnorm(1,R_max,SD_lambda))*(1-(prev_abund/K)))
}

GrazeSeedbank <- function(prev_abund,lastyear){
  if(lastyear > 4000){
    x <- otherwise*lastyear
    return(x)
  }
  x <- rnorm(1,grazing,sdgrazing)*(reprodensityform[[2]]*log(lastyear)+reprodensityform[[1]])*lastyear
  #to ensure 
  if(is.nan(x)|is.infinite(x)){x <- 0}
  if(x==0){
    x <- rnorm(1,grazing,sdgrazing)*lastyear*rnorm(1,seedoutput,seedoutputsd)
  }
  return(x)
}

GrazePVA <- function(nreps,nyears,Init_N,R_max){
  
  PopArray2 <- array(0,dim=c((nyears+1),nreps))
  Seeds <- array(0,dim=c((nyears+1),nreps))
  
  for(rep in 1:nreps){
    
    PopArray2[1,rep] <- rnorm(1,grazingabund,grazingabundsd)*Init_N
    Seeds[1,rep] <- max(0,rnorm(1,grazing,sdgrazing))*(reprodensityform[[2]]*log(Init_N)+reprodensityform[[1]])*Init_N
    
    
    for(y in 2:(nyears+1)){
      
      nextyear <- max(0,Ricker(PopArray2[y-1,rep]))
      if(nextyear==0){
        nextyear <- PopArray2[y-1,rep]
      }
      PopArray2[y,rep] <- nextyear
      nextyearseed <- max(0,GrazeSeedbank(Seeds[y-1,rep],PopArray2[y-1,rep]))
      Seeds[y,rep] <- nextyearseed
    }
  }
  alist$poparray <- PopArray2
  alist$seeds <- Seeds
  return(alist)
}

GrazeDefault <- GrazePVA(nreps,nyears,Init_N,R_max)

grazeendpop <- mean(GrazeDefault$poparray[nrow(GrazeDefault$poparray),])
grazeendseed <- mean(GrazeDefault$seeds[nrow(GrazeDefault$seeds),])


agraze <- grazeendpop
sgraze <- sd(GrazeDefault$poparray[nrow(GrazeDefault$poparray),])
ngraze <- nreps
errorgraze <- qnorm(0.975)*sgraze/sqrt(ngraze)
leftgraze <- agraze-errorgraze
rightgraze <- agraze+errorgraze
cigraze <- c(leftgraze,rightgraze)


agrazes <- grazeendseed
sgrazes <- sd(GrazeDefault$seeds[nrow(GrazeDefault$seeds),])
ngrazes <- nreps
errorgrazes <- qnorm(0.975)*sgrazes/sqrt(ngrazes)
leftgrazes <- agrazes-errorgrazes
rightgrazes <- agrazes+errorgrazes
cigrazes <- c(leftgrazes,rightgrazes)



