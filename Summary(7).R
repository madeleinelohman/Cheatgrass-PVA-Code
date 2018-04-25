mydfpop <- c(EndControl$poparray[nrow(EndControl$poparray),],
                  BurnDefault$poparray[nrow(BurnDefault$poparray),],
                  GrazeDefault$poparray[nrow(GrazeDefault$poparray),],
                  GrazeBurnDefault$poparray[nrow(GrazeBurnDefault$poparray),])
mylistpop <- c("Control","Graze","Burn","GrazeBurn")
mydf1pop <- data.frame(rep(NA,6000),rep(mylistpop,each=1500))
mydf1pop[,1] <- mydfpop
colnames(mydf1pop) <- c("values","treatment")
mynewpop <- aov(values ~ treatment, data = mydf1pop)
summary(mynewpop)


endpop
grazeendpop
burnendpop
grazeburnendpop

popci <- c(endci,cigraze,ciburn,cigrazeburn)
popcomparepoints <- c(endpop,grazeendpop,burnendpop,grazeburnendpop)
pop <- data.frame(c(popci,popcomparepoints))
popcomparetreatment <- c("Control","Grazing","Burning",paste("Grazing","and","Burning",sep="\n"))

barCenters <- barplot(popcomparepoints,names.arg = popcomparetreatment,beside = true, las = 2,cex.names = 1, 
                      main = "Population after 6 Years",ylab = expression("Plants per"~m^{2}),cex=0.9,border = "black",
                      ylim=c(0,max(popci)))

segments(barCenters, popci[c(2,4,6,8)], barCenters,
         popci[c(1,3,5,7)], lwd = 1.5)

arrows(barCenters, popci[c(2,4,6,8)], barCenters,
       popci[c(1,3,5,7)], lwd = 1.5, angle = 90,
       code = 3, length = 0.05)



mydfs <- c(EndControl$seeds[nrow(EndControl$seeds),],
          GrazeDefault$seeds[nrow(GrazeDefault$seeds),],
          BurnDefault$seeds[nrow(BurnDefault$seeds),],
          GrazeBurnDefault$seeds[nrow(GrazeBurnDefault$seeds),])
mylist <- c("Control","Graze","Burn","GrazeBurn")
mydf1s <- data.frame(rep(NA,6000),rep(mylist,each=1500))
mydf1s[,1] <- mydfs
colnames(mydf1s) <- c("values","treatment")
mynews <- aov(values ~ treatment, data = mydf1)
summary(mynews)


endseed
grazeendseed
burnendseed
grazeburnendseed

seedci <- c(endcis,cigrazes,ciburns,cigrazeburns)
seedcomparepoints <- c(endseed,grazeendseed,burnendseed,grazeburnendseed)
seedcomparetreatment <- c("Control","Grazing","Burning",paste("Grazing","and","Burning",sep="\n"))

barCenterss <- barplot(seedcomparepoints,names.arg=seedcomparetreatment,beside = true, las = 2,cex.names = 1, 
                      main = "Seedbank after 6 Years",ylab = expression("Seeds per"~m^{2}),cex=0.9,border = "black",
                      ylim=c(0,max(seedci)))

segments(barCenters, seedci[c(2,4,6,8)], barCenters,
         seedci[c(1,3,5,7)], lwd = 1.5)

arrows(barCenters, seedci[c(2,4,6,8)], barCenters,
       seedci[c(1,3,5,7)], lwd = 1.5, angle = 90,
       code = 3, length = 0.05)







