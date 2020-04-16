library(ggthemes)
library(tidyverse)
library(scales)
library(ggplot2)
library(ggridges)

source("../R/read-admb.R")
source("../R/plot_ssb.R")
m0 <- read_rep("m0/pop_R.rep")
m1 <- read_rep("m1/pop_R.rep")
m2 <- read_rep("m2/pop_R.rep")
m3 <- read_rep("m3/pop_R.rep")
#-------------------------------------------------------------------------------
#Selectivity
psel <- function(M,title="",styr=1961) {
  df  <- data.frame(M$Selectivity[,1:29] ); names(df) <- c("yr",2:29)
  sdf <- gather(df,age,sel,2:29) %>% filter(yr>styr) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
  ggplot(sdf,aes(x=age,y=as.factor(yr),height = sel)) + geom_density_ridges(stat = "identity",scale = 5.8, alpha = .2,color="blue",fill="yellow",size=.5)+ xlim(c(1,29))+ ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$yr)))) + theme_few() + ggtitle(title)
}
psel(m0,title="Hulson, blocks")
psel(m1,title="3pDL, blocks")
psel(m2,title="3pDL, 3yr")
psel(m3,title="3pDL, 2yr")

#-------------------------------------------------------------------------------
#SSB
ssb1<- data.frame(m0$SSB,model="Hulson")
ssb2<- data.frame(m3$SSB,model="3pDL")
ssb <- rbind(ssb1,ssb2)
ssb <- rbind(ssb.mle)
names(ssb) <- c("Year","SSB", "SD","LB","UB","Model")
ggplot(ssb,aes(x=Year,y=SSB, ymin=LB,ymax=UB,fill=Model)) + geom_ribbon( alpha = .3) + xlab("Year") + theme_few()
ssb <- data.frame(m$SSB)
names(ssb) <- c("Year","SSB", "SD","LB","UB")
ggplot(ssb,aes(x=Year,y=SSB, ymin=LB,ymax=UB),fill='salmon') + geom_ribbon( alpha = .3) + xlab("Year") + theme_few()

#-------------------------------------------------------------------------------
#Rec
rdf <- cbind(data.frame(m1$R,case="3pDL, blocks"))
rdf <- rbind(rdf,cbind(data.frame(m0$R,case="Hulson, blocks")))
names(rdf) <- c("yr","R","se","lb","ub","case")
#rdf  <- rdf %>% filter(yr>1970,yr<2020)
mnR <- mean(m0$R[,2])
dodge <- position_dodge(width=0.8)
ggplot(rdf,aes(x=yr-1,y=R,fill=case)) + xlab("Year class") + ylab("Age 2 recruits (thousands)") + 
       geom_bar(width=0.75,position="dodge",stat="identity",color="black") + 
       scale_x_continuous(breaks=seq(1970,2019,5)) +
       geom_errorbar(aes(ymin=lb,ymax=ub),width=.3,colour="blue",position=dodge) + theme_few() + geom_hline(aes(yintercept=mnR))


#-------------------------------------------------------------------------------
# Look at fits... NOTE this is read_admb which gets lots more than read_rep used above...
m0f <- read_admb("m0/pop")
m1f <- read_admb("m1/pop")
m2f <- read_admb("m2/pop")
(m0f$fit$nlogl)
(m1f$fit$nlogl)
(m2f$fit$nlogl)

#-------------------------------------------------------------------------------
#Plot fit to surve (NOT WORKING YET)
# source("R/plot_ind.R")
# lstOut1  <- list( "Base"= m0, "3PdL blocks"=m1)
# plot_ind(lstOut1)


