#population values for control experiments
x <- c(660,4000,178,90,4850)
#seed output per plant for control
y <- c(10.7,1.75,10,28,1.61)
mydf5 <- data.frame(x,y)

#mean of seed output at highest densities; used to calculate seed output at K
otherwise <- mean(c(1.75,1.61))

plot(mydf5,xlab=expression("Plants per"~m^{2}),ylab="Seeds per plant",pch=19,ylim=c(0,30),
     main="Density Dependent Seed Production")

#to produce logarithmic line, predict function was used
fit <- lm(y~log(x))
x=seq(from=1,to=5000,length.out=1000)
y=predict(fit,newdata=list(x=seq(from=1,to=5000,length.out=1000)))
mydf6 <- data.frame(x,y)
lines(y~x,data=mydf6,xpd=F)

