###Application of kriging when the rainfall intensities are not normally distributed ###

#It can be applied for any cases even when the data is normally distributed, but
#it's computationally demanding compare to kriging_nor.R

#v4: extrapolation method changed

install.packages("sp","rgdal","gstat","rgeos","zoo","ggplot2","calibrate","reshape","xts")
library(sp)
library(rgdal)
library(gstat)
library(rgeos)
library(zoo) 
library(ggplot2)
library(calibrate) # textxy
library(reshape)
library(xts)

setwd('..')
source("./2_Scripts/kriging_genutil.R")
source("./2_Scripts/weightRain.R")


num.pts=NULL
############################Part-1: Estimation of Variogram #######################################
#inputs
dat.year=c(113)
npoint=8 #number of stations
avr.int=30 #averaging interval of rainfall intensity (min)
lowthres=10
upthres=200
vgbound=c(25,50,100,150,200,250,300,350,400)
#vgbound=c(75,105,125,175,250,305,400)
#vgbound=c(25,50,100,200,300,400)
#vgbound=c(1:5*100,1200,1500,2000,2500,3000,4000)


#set threshold
rain.data=read.zoo(paste0("./1_Data/int_brad_",avr.int,".csv"),sep = "," ,
                   header = TRUE,tz="",format = "%d/%m/%Y %H:%M",row.names=NULL)
rain.data=xts(rain.data)
#rain.data=rain.data[.indexyear(rain.data)%in% dat.year]
rain.data=zoo.to.data.frame(rain.data, "time")
coord=read.csv("./1_Data/loc_brad.csv")
geo_data=geostat_data(rain.data,coord,lowthres,upthres,samp.per=80)
nrow(geo_data)

#Table for paper
num.pts=rbind(num.pts,nrow(geo_data))
write.csv(num.pts,"./3_Results/num.pts.csv")


# read and explore data
d <- data.frame(geo_data)
names(d); dim(d)
max(d$x) - min(d$x); max(d$y) - min(d$y)
hist(d$int, col="Lightblue")

# small shift of x-coordinate to avoid predicting at observation locations
d$x <- d$x + 0.01

# # for paper
# #start
# D1=d$int
# D2=d$int
# D3=d$int
# 
# #plot
# Dall=c(D1,D2,D3)
# names=c(rep("0.1-5.0 mm/hr",length(D1)),rep("5.1-10.0 mm/hr",length(D2)),
#         rep("> 10.1 mm/hr",length(D3)))
# D.df=data.frame(Dall,names)
# D.df$names = factor(D.df$names, levels=unique(D.df$names))
# theme_set(theme_bw())
# ggplot(data=D.df, aes(D.df$Dall)) + 
#   labs(x="Intensity [mm/hr]", y="Count [-]")+
#   geom_freqpoly()+
#   #geom_histogram(fill="transparent",border="black")+
#   facet_wrap(~names,ncol = 3,scales = "free")
# 
# #end

#1. standardisation steps
mean.tp=rollapply(d$int,width=npoint,by=npoint, FUN=mean)
sd.tp=rollapply(d$int,width=npoint,by=npoint, FUN=sd)
plot(mean.tp, type="b"); plot(sd.tp)

#2. copy data to calculate standardised set for which stationary variogram is calculated
d2 <- d #Becouse we need d later
d2$mean.tp=rep(mean.tp,each=npoint)
d2$sd.tp=rep(sd.tp,each=npoint)
d2$std.int=(d2$int-d2$mean.tp)/d2$sd.tp
hist(d2$std.int, col="Lightblue")

# # for paper
# #start
# D1=d2$std.int
# D2=d2$std.int
# D3=d2$std.int
# 
# #plot
# Dall=c(D1,D2,D3)
# names=c(rep("< 5.0 mm/h",length(D1)),rep("5.0-10.0 mm/h",length(D2)),
#         rep("> 10.0 mm/h",length(D3)))
# D.df=data.frame(Dall,names)
# D.df$names = factor(D.df$names, levels=unique(D.df$names))
# theme_set(theme_bw())
# ggplot(data=D.df, aes(D.df$Dall)) +
#   labs(x="Standardised intensity [-]", y="Count [-]")+
#   geom_freqpoly(binwidth=0.3)+
#   #geom_histogram(fill="transparent",border="black")+
#   facet_wrap(~names,ncol = 3,scales = "free_y")+
#   theme(axis.text=element_text(size=8),
#         axis.title=element_text(size=10,face="bold"),
#         strip.text.x = element_text(size = 8),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())
# #end

#3. NQT transformation (Chose a transformation method which makes the data to roughly normally distributed)
intNorm=qqnorm(d2$std.int, plot.it = T) #transformed to normal domain
hist(intNorm$x)
d2$nqt.int=intNorm$x

# for paper
#start
D1=d2$nqt.int
D2=d2$nqt.int
D3=d2$nqt.int

#plot
Dall=c(D1,D2,D3)
names=c(rep("< 5.0 mm/h",length(D1)),rep("5.0-10.0 mm/h",length(D2)),
        rep("> 10.0 mm/h",length(D3)))
D.df=data.frame(Dall,names)
D.df$names = factor(D.df$names, levels=unique(D.df$names))
theme_set(theme_bw())
ggplot(data=D.df, aes(D.df$Dall)) + 
  labs(x="Standardised intensity - after NST [-]", y="Count [-]")+
  geom_freqpoly(binwidth=0.3)+
  #geom_histogram()+
  facet_wrap(~names,ncol = 3,scales = "free_y")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"),
        strip.text.x = element_text(size = 8),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
#end
# #for extrapolation during backtransformation
# lm.head=lm(head(sort(intNorm$x),npoint)~head(sort(intNorm$y),npoint))
# lm.tail=lm(tail(sort(intNorm$x),npoint)~tail(sort(intNorm$y),npoint))
# 
# plot(head(sort(intNorm$y),npoint),head(sort(intNorm$x),npoint))
# abline(lm.head[[1]][[1]],lm.head[[1]][[2]])
# 
# plot(tail(sort(intNorm$y),npoint),tail(sort(intNorm$x),npoint))
# abline(lm.tail[[1]][[1]],lm.tail[[1]][[2]])

#4. variogram estimation
# large shift x-coordinate each time instant to enable variogram estimation using all data
# to achieve that only paired comparisons of observations at the same time are used
d2$x <- d2$x + 100000*(rep(1:(nrow(d2)/npoint),each=npoint)-1)

coordinates(d2) <- ~x+y
g <- gstat(formula=nqt.int~1, data=d2)
vg <- variogram(g, boundaries=vgbound)
vgm <- vgm(psill=max(vg$gamma), "Exp", range=100, nugget=vg$gamma[1])
vgm <- fit.variogram(vg, vgm)
plot(vg, vgm, plot.nu=T,main="Avr. Int : 15 min, Range :1-5 mm/hr")

#assign variograms according to range to plot (2min)
#range 1: 1-5 mm/hr
vgm1=vgm
vg1=vg

#range 2: 5-10 mm/hr
vgm2=vgm
vg2=vg

#range 3: >10 mm/hr
vgm3=vgm
vg3=vg

#assign variograms according to range to plot (5min)
#range 1: 1-5 mm/hr
vgm4=vgm
vg4=vg

#range 2: 5-10 mm/hr
vgm5=vgm
vg5=vg

#range 3: >10 mm/hr
vgm6=vgm
vg6=vg

#assign variograms according to range to plot(15min)
#range 1: 1-5 mm/hr
vgm7=vgm
vg7=vg

#range 2: 5-10 mm/hr
vgm8=vgm
vg8=vg

#range 3: >10 mm/hr
vgm9=vgm
vg9=vg

#assign variograms according to range to plot (30min)
#range 1: 1-5 mm/hr
vgm10=vgm
vg10=vg

#range 2: 5-10 mm/hr
vgm11=vgm
vg11=vg

#range 3: >10 mm/hr
vgm12=vgm
vg12=vg

#saving the results
vgm_all=rbind(vgm1,vgm2,vgm3,vgm4,vgm5,vgm6,vgm7,vgm8,vgm9,vgm10,vgm11,vgm12)
save(vgm_all, file = "./1_Data/vgm_all13.RData")
save.image()
write.csv(vgm_all,"./1_Data/vgm_all13.csv") #doesn't keep the format

vg_all=rbind(vg1,vg2,vg3,vg4,vg5,vg6,vg7,vg8,vg9,vg10,vg11,vg12)
vg_all$ts=rep(c("2 min","5 min","15 min" ,"30 min"),each=27)
vg_all$range=rep(c("< 5.0 mm/h","5.0-10.0 mm/h","> 10.0 mm/h"),each=9,times=4)
write.csv(vg_all,"./1_Data/vg_all.csv")

#Additional work for referee Comments (narrower intensity class)
vgm_all=rbind(vgm1,vgm2,vgm3,vgm4,vgm5,vgm6,vgm7,vgm8,vgm9)
save(vgm_all, file = "./1_Data/vgm_all_112416.RData")
save.image()
write.csv(vgm_all,"./1_Data/vgm_all_112416.csv") #doesn't keep the format

vg_all=rbind(vg1,vg2,vg3,vg4,vg5,vg6,vg7,vg8,vg9)
write.csv(vg_all,"./1_Data/vg_all_112416.csv")

#Additional work for referee Comments (80% data)
vgm_all=rbind(vgm1,vgm2,vgm3,vgm4,vgm5,vgm6,vgm7,vgm8,vgm9,vgm10,vgm11,vgm12)
save(vgm_all, file = "./1_Data/vgm_all_1124.RData")
save.image()
write.csv(vgm_all,"./1_Data/vgm_all_112416.csv") #doesn't keep the format

vg_all=rbind(vg1,vg2,vg3,vg4,vg5,vg6,vg7,vg8,vg9,vg10,vg11,vg12)
vg_all$ts=rep(c("2 min","5 min","15 min" ,"30 min"),each=27)
vg_all$range=rep(c("< 5.0 mm/h","5.0-10.0 mm/h","> 10.0 mm/h"),each=9,times=4)
write.csv(vg_all,"./1_Data/vg_all_112416.csv")

#saving the plot-for variograms of 1-5 mm/hr on a single plot
#start
vgm_1to5=rbind(vgm1,vgm2,vgm3)
save(vgm_1to5, file = "./1_Data/vgm_1to5.RData")
save.image()
write.csv(vgm_1to5,"./1_Data/vgm_1to5.csv") #doesn't keep the format

vg_1to5=rbind(vg1,vg2,vg3)
vg_1to5$ts=rep(c("30 min","60 min","120min"),each=9)
write.csv(vg_1to5,"./1_Data/vg_1to5.csv")
#end

#Plot all variograms in a single plot
plot(gamma~dist, vg1, ylim = c(0, 2), #1.05*max(vg1$gamma,vg2$gamma,vg3$gamma)),
     xlim=c(0,400), pty=1, col='red', main=paste0("Avr. Interval - ",avr.int,"min"),
     ylab = 'semivariance', xlab = 'distance (m)') 
lines(variogramLine(vgm1, 400), lty=1,col='red')
textxy(vg1$dist, vg1$gamma,  vg1$np,col="red")

points(gamma~dist, vg2, pch=2, col='blue') 
lines(variogramLine(vgm2, 400), lty=2,col='blue')
textxy(vg2$dist, vg2$gamma,  vg2$np,col="blue")

points(gamma~dist, vg3, pch=3, col='green')
lines(variogramLine(vgm3, 400), lty=3,col='green')
textxy(vg3$dist, vg3$gamma,  vg3$np,col="green")

legend("bottomright",title="Intensity Range (mm/hr)", legend=c("0.1-5", "5.1-10", ">10.1"), 
       pch=c(1,2,3),lty=c(1,2,3), col=c("red","blue","green"),bty="n")


#Plot all variograms in a single plot - linear model
plot(gamma~dist, vg1, ylim = c(0, 2), #1.05*max(vg1$gamma,vg2$gamma,vg3$gamma)),
     xlim=c(0,400), pty=1, col='red', main=paste0("Avr. Interval - ",avr.int,"min"),
     ylab = 'semivariance', xlab = 'distance (m)')
lin.mod=lm(gamma~dist,vg1)
abline(a=lin.mod[[1]][1],b=lin.mod[[1]][2], lty=1,col='red')
textxy(vg1$dist, vg1$gamma,  vg1$np,col="red")

points(gamma~dist, vg2, pch=2, col='blue') 
lin.mod=lm(gamma~dist,vg2)
abline(a=lin.mod[[1]][1],b=lin.mod[[1]][2], lty=2,col='blue')
textxy(vg2$dist, vg2$gamma,  vg2$np,col="blue")

points(gamma~dist, vg3, pch=3, col='green')
lin.mod=lm(gamma~dist,vg3)
abline(a=lin.mod[[1]][1],b=lin.mod[[1]][2], lty=3,col='green')
textxy(vg3$dist, vg3$gamma,  vg3$np,col="green")

legend("bottomright",title="Intensity Range (mm/hr)", legend=c("0.1-5", "5.1-10", ">10.1"), 
       pch=c(1,2,3),lty=c(1,2,3), col=c("red","blue","green"),bty="n")


############################Part-2: Application of Kriging#######################################

#select a event
npoint=8
timescale=c(30) #min
event=c(1)
coord=read.csv("./1_Data/loc_brad.csv")
#event=c(11) #for peaks
load("./1_Data/vgm_all.RData")

for(m in 1:12) 
{ 
  nam <- paste0("vgm", m)
  assign(nam, vgm_all[m:(m+1),])
}

kriging.all=NULL

for (k in event) {
  for (l in timescale) {
    
    #rain.data=read.csv(paste0("./3_Results/peaks_",l,".csv"))[,-c(1,11)] #only peaks
    rain.data=read.csv(paste0("./3_Results/allevents_",l,".csv")) #all 30 min
    #rain.data=read.csv(paste0("./3_Results/",k,"_",l,".csv")) # for each event
    geo_data=geostat_data(rain.data,coord,0,200) # doesn't include 0
    
    #1. Select vgm models
    if (l==2) {
      vgm_l=vgm1
      vgm_m=vgm2
      vgm_h=vgm3
    } else if (l==5) {
      vgm_l=vgm4
      vgm_m=vgm5
      vgm_h=vgm6
    } else if (l==15) {
      vgm_l=vgm7
      vgm_m=vgm8
      vgm_h=vgm9
    } else if (l==30) {
      vgm_l=vgm10
      vgm_m=vgm11
      vgm_h=vgm12
    } 
    
    d <- geo_data
    
    mean.tp=rollapply(d$int,width=npoint,by=npoint, FUN=mean)
    sd.tp=rollapply(d$int,width=npoint,by=npoint, FUN=sd)
    plot(mean.tp, type="b"); plot(sd.tp)
    
    #2. copy data to calculate standardised set 
    
    d$mean.tp=rep(mean.tp,each=npoint)
    d$sd.tp=rep(sd.tp,each=npoint)
    d$std.int=(d$int-d$mean.tp)/d$sd.tp
    hist(d$std.int, col="Lightblue")
    
    #3. NQT transformation (Chose a transformation method which makes the data to roughly normally distributed)
    intNorm=qqnorm(d$std.int, plot.it = T) #transformed to normal domain
    hist(intNorm$x)
    d$nqt.int=intNorm$x
    
    #for extrapolation during backtransformation
    lm.head=lm(head(sort(intNorm$x),npoint)~head(sort(intNorm$y),npoint))
    lm.tail=lm(tail(sort(intNorm$x),npoint)~tail(sort(intNorm$y),npoint))
    
    plot(head(sort(intNorm$y),npoint),head(sort(intNorm$x),npoint))
    
    if (is.na(lm.head[[1]][[2]])) {
      abline(v=lm.head[[1]][[1]])
    } else {
      abline(lm.head[[1]][[1]],lm.head[[1]][[2]])
    }
   
    
    plot(tail(sort(intNorm$y),npoint),tail(sort(intNorm$x),npoint))
    if (is.na(lm.tail[[1]][[2]])) {
      abline(v=max(intNorm$y))
    } else {
      abline(lm.tail[[1]][[1]],lm.tail[[1]][[2]])
    }
        
    act.pred.all=NULL
    act.var.all=NULL
    time.all=NULL
    
    row.index=(0:(nrow(geo_data)/npoint-1)*npoint+1)
    
    for (i in row.index) {
      #1. choose a time instant
      
      d2 <- d[i:(i+npoint-1),]
      time <- d$time[i]
      coordinates(d2) <- ~x+y
      
      #4. Calculation of point kriging spatial stochastic application 
      
      # define grid
      xy <- expand.grid(seq(415300,415735,25), seq(432715, 432890,25))
      names(xy) <- c('x','y')
      gridded(xy) <- ~x+y
      
      # using spatial stochastic simulation, first simulate at point support
      if (mean(d2$int)>10) {
        g <- gstat(formula=nqt.int~1, locations=d2, model=vgm_h, nmax=100)
        pksim <- predict(g, newdata=xy, nsim=500, debug.level=-1)
        
        
      } else if (mean(d2$int)>5) {
        g <- gstat(formula=nqt.int~1, locations=d2, model=vgm_m, nmax=100)
        pksim <- predict(g, newdata=xy, nsim=500, debug.level=-1)
        
        
      } else if (mean(d2$int)>0) {
        g <- gstat(formula=nqt.int~1, locations=d2, model=vgm_l, nmax=100)
        pksim <- predict(g, newdata=xy, nsim=500, debug.level=-1)
        
      }
      
      #5. Back transformation
      pksim.ori=apply(pksim@data,c(1,2),function (x) approx(intNorm$x,intNorm$y,x,rule=1)$y)
      
      #use extrapolation
      aaa=which(is.na(pksim.ori), TRUE)
      if (nrow(aaa)>0) {
        
        for (j in 1:nrow(aaa)){
          
          xxx=pksim@data[aaa[j,1],aaa[j,2]]
          
          if (xxx > max(intNorm$x)) {
            
            if (is.na(lm.tail[[1]][[2]])){
              pksim.ori[aaa[j,1],aaa[j,2]]=max(intNorm$y)
            } else {
              pksim.ori[aaa[j,1],aaa[j,2]]=(xxx-lm.tail[[1]][[1]])/lm.tail[[1]][[2]]
            } 
            
          } else {
            if (is.na(lm.head[[1]][[2]])){
              pksim.ori[aaa[j,1],aaa[j,2]]=min(intNorm$y)
            } else {
              pksim.ori[aaa[j,1],aaa[j,2]]=(xxx-lm.head[[1]][[1]])/lm.head[[1]][[2]]
            }
            
          }
        }
      }
      
      # Predicted block average equals the mean and standard deviation of spatial means over realisations
      plot(colMeans(pksim.ori))
      bksim.pred=mean(colMeans(pksim.ori))
      bksim.var=var(colMeans(pksim.ori))
      
      #8. transformation 2
      act.pred=bksim.pred*sd.tp[which(row.index==i)]+mean.tp[which(row.index==i)]
      act.var= bksim.var*(sd.tp[which(row.index==i)])^2
      
      act.pred.all=cbind(act.pred.all,act.pred)
      act.var.all=cbind(act.var.all,act.var)
      time.all=cbind(time.all,time)
      
    }
    
    
    #plotting the results
    
    times=geo_data[,3][0:(nrow(geo_data)/npoint-1)*npoint+1]
    event.kriging=data.frame(matrix(NA,nrow=length(act.var.all),
                                    ncol=3))
    event.kriging[,1]=times
    event.kriging[,2]=t(act.pred.all)
    event.kriging[,3]=t(sqrt(act.var.all))
    colnames(event.kriging)=c("time","pred","sd")
    
    event.kriging$CI95up=qnorm(0.975, mean=event.kriging$pred, sd=event.kriging$sd)
    event.kriging$CI95low=qnorm(0.025, mean=event.kriging$pred, sd=event.kriging$sd)
    
    event.kriging$CI50up=qnorm(0.75, mean=event.kriging$pred, sd=event.kriging$sd)
    event.kriging$CI50low=qnorm(0.25, mean=event.kriging$pred, sd=event.kriging$sd)
    
    # event.kriging$CI95up=event.kriging$pred+1.96*event.kriging$sd
    # event.kriging$CI95low=event.kriging$pred-1.96*event.kriging$sd
    
    #kriging.zoo=read.zoo(event.kriging, header = TRUE, index = 1,tz = "",format = "%d/%m/%Y %H:%M") # some times the below one doesn't work don't know the reason
    kriging.zoo=read.zoo(event.kriging)
    #meas.zoo=read.zoo(rain.data, header = TRUE, index = 1,tz = "",format = "%d/%m/%Y %H:%M") # some times the below one doesn't work don't know the reason
    meas.zoo=read.zoo(rain.data)
    kriging.data=merge.zoo(meas.zoo,kriging.zoo, fill=0)
    plot(kriging.data$pred, ylim=c(0,1.25*max(kriging.data$pred)),
         xlab=paste0("Time (",as.Date(index(kriging.data[1,])),")"),
         ylab="Rainfall Intensity [mm/hr]")
    
    polygon(c(index(kriging.data$CI95up),rev(index(kriging.data$CI95low))),
            c(coredata(kriging.data$CI95up),coredata(rev(kriging.data$CI95low))),
            col="gray86",border=NA)
    
    polygon(c(index(kriging.data$CI50up),rev(index(kriging.data$CI50low))),
            c(coredata(kriging.data$CI50up),coredata(rev(kriging.data$CI50low))),
            col="gray48",border=NA)
    
    lines(kriging.data$pred,type="b",lwd=2)
    legend("topright",bty="n", 
           legend=c("Kriging Prediction","50% CI", "95% CI"),
           col = c("black", "gray48","gray86"),
           pch=c(1,22,22),lty=c(1,NA,NA), pt.bg=c(NA, "gray48","gray86"),
           pt.cex=c(1,3,3),lwd=c(2,NA,NA),pt.lwd=c(NA,NA,NA))
        
    #using ggplot
    kriging.gg=kriging.data[,c(9,11:14)]
    autoplot(kriging.gg, facet = NULL)
    
    # #plotting kriging vs average and thiessen
    # prad.kri=kriging.data[,9]
    # pred.avr=zoo(rowMeans(meas.zoo),time(meas.zoo))
    # #pred.thi
    # thi.wts=read.csv("thiessenWeights.csv",header=F)[,1]
    # pred.thi=zoo(weightRain(coredata(meas.zoo),thi.wts,1:8),time(meas.zoo))
    # pred.all=merge.zoo(pred.avr,pred.thi,prad.kri)
    # autoplot(pred.all, facet = NULL)
  
    #confidant bands 
    z=fortify(kriging.gg)
    ggplot(z, aes(x=Index, y=pred, group = 1)) +
      geom_line() +
      geom_ribbon(aes(ymax = CI50up, ymin = CI50low), fill = 1, alpha = 0.2) +
      geom_ribbon(aes(ymax = CI95up, ymin = CI95low), fill = 2, alpha = 0.2)
    
    #kriging.all=rbind(kriging.all,coredata(kriging.data[,-c(1:8)])) #for peaks
  }
  kriging.all=rbind(kriging.all,kriging.data[,-c(1:8)])
}


#preparing data
kriging.2=kriging.data[,-c(1:8)]
kriging.2=data.frame(Date=time(kriging.2),coredata(kriging.2))
kriging.2$ts=rep("2 min", nrow(kriging.2))

kriging.5=kriging.data[,-c(1:8)]
kriging.5=data.frame(Date=time(kriging.5),coredata(kriging.5))
kriging.5$ts=rep("5 min", nrow(kriging.5))

kriging.15=kriging.data[,-c(1:8)]
kriging.15=data.frame(Date=time(kriging.15),coredata(kriging.15))
kriging.15$ts=rep("15 min", nrow(kriging.15))

kriging.30=kriging.data[,-c(1:8)]
kriging.30=data.frame(Date=time(kriging.30),coredata(kriging.30))
kriging.30$ts=rep("30 min", nrow(kriging.30))

kriging.comb=rbind(kriging.2,kriging.5,kriging.15,kriging.30)
kriging.comb$ts = factor(kriging.comb$ts, levels=unique(kriging.comb$ts)) #to adjust the order of facet_wrap in ggplot
write.csv(kriging.comb,"./3_Results/kriging.comb.csv")

#if read from csv file
kriging.comb=read.csv("./1_Data/kriging.comb.csv")
kriging.comb$Date = as.POSIXct(kriging.comb$Date,format='%Y-%m-%d %H:%M:%S')
kriging.comb$ts = factor(kriging.comb$ts, levels=unique(kriging.comb$ts)) #to adjust the order of facet_wrap in ggplot

