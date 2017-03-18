#Reference:
## http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/

library(gstat)
library(zoo) 
library(ggplot2)
library(grid) # to use "unit" in ggplot

setwd('..')
source("./2_Scripts/kriging_genutil.R")

#1. Variogram plots
vg_all=read.csv("./1_Data/vg_all_112416.csv")
load("./1_Data/vgm_all_112416.RData")

vg_all
vg_all = vg_all[-c(1,5,6,7)]
vg_all$log_np = log(vg_all$np)
vg_all$ts = factor(vg_all$ts, levels=unique(vg_all$ts)) #to adjust the order of facet_wrap in ggplot
vg_all$range = factor(vg_all$range, levels=unique(vg_all$range)) #to adjust the order of facet_wrap in ggplot

fitted = seq(0.01, max(vg_all$dist), length = 100)

gamma_all=NULL
abc=seq(1,23,2)
for ( i in abc) {
  gamma = variogramLine(vgm_all[i:(i+1),], dist_vector = fitted)$gamma
  gamma_all=cbind(gamma_all,gamma)
}
gamma_all=c(gamma_all)

fitted=cbind(rep(fitted,24),gamma_all)
fitted=as.data.frame(fitted)
colnames(fitted)=c("dist","gamma")
fitted$ts=rep(c("2 min","5 min","15 min" ,"30 min"),each=300)
fitted$range=rep(c("< 5.0 mm/h","5.0-10.0 mm/h","> 10.0 mm/h"),each=100,times=4)

fitted$ts = factor(fitted$ts, levels=unique(fitted$ts))  #to adjust the order of facet_wrap in ggplot
fitted$range=factor(fitted$range, levels=unique(fitted$range))  #to adjust the order of facet_wrap in ggplot

theme_set(theme_bw())

#use this fond/pointer size for larger plot
ggplot(vg_all, aes(x = dist, y = gamma, shape=range,colour=range)) + 
  geom_point(size=2.5)+
  geom_line(data=fitted,size=0.7)+
  labs(x="Distance [m]", y="Semivariance [-]")+
  facet_wrap(~ts)+
  scale_linetype_manual(values=c("solid","solid","solid"))+
  theme(legend.position="bottom",legend.title=element_blank(),
        legend.text=element_text(face="bold",size=10),legend.key=element_blank(),
        legend.key.width=unit(1,'cm'),legend.margin = unit(0.05, "cm"),legend.key.height=unit(0.05, "cm"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=10,face="bold"),
        axis.title=element_text(size=12,face="bold"))+
  theme(strip.text.x = element_text(size = 10))
  
#for paper 
ggplot(vg_all, aes(x = dist, y = gamma, shape=range,linetype=range)) + 
  geom_point(size=1.8)+
  geom_line(data=fitted,size=0.7)+
  #geom_line(data=fitted,size=0.7)+
  labs(x="Distance [m]", y="Semivariance [-]")+
  scale_x_continuous(breaks=c(0,100,200,300,400))+
  facet_wrap(~ts)+
  scale_linetype_manual(values=c("solid","dotted","longdash"))+
  theme(legend.position="bottom",legend.title=element_blank(),
        legend.text=element_text(face="bold",size=8),legend.key=element_blank(),
        legend.key.width=unit(1,'cm'),legend.margin = unit(0.05, "cm"),legend.key.height=unit(0.05, "cm"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"))+
  theme(strip.text.x = element_text(size = 8))

#Additional work for referee Comments
vg_all=read.csv("./1_Data/vg_all_112416.csv")
vg_all = vg_all[-c(64:81),-c(1,5,6,7)]
vg_all$range=rep(c("<2.0",">=2.0-<4.0",">=4.0-<6.0",">=6.0-<8.0",">=8.0-<10.0",">=10.0-<12.0",">=12.0-<14.0"),each=9)


load("./1_Data/vgm_all_112416.RData")
fitted = seq(0.01, max(vg_all$dist), length = 100)

gamma_all=NULL
abc=seq(1,13,2)
for ( i in abc) {
  gamma = variogramLine(vgm_all[i:(i+1),], dist_vector = fitted)$gamma
  gamma_all=cbind(gamma_all,gamma)
}
gamma_all=c(gamma_all)

fitted=cbind(rep(fitted,14),gamma_all)
fitted=as.data.frame(fitted)
colnames(fitted)=c("dist","gamma")
#fitted$range=rep(c("\u2265 2","\u2265 4","\u2265 6","\u2265 8","\u2265 10","\u2265 12","\u2265 14"),each=100) # doesn't work
fitted$range=rep(c("(0,2]","(2,4]","(4,6]","(6,8]","(8,10]","(10,12]","(12,14]"),each=100)
fitted$range=factor(fitted$range, levels=unique(fitted$range))

ggplot(vg_all, aes(x = dist, y = gamma, color=range)) + 
  #geom_point(size=1.8)+
  geom_line(data=fitted,size=0.8)+
  #geom_line(data=fitted,size=0.7)+
  labs(x="Distance [m]", y="Semivariance [-]")+
  scale_x_continuous(breaks=c(0,100,200,300,400))+
  scale_colour_discrete(name="Intensity range\n (mm/hr)")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"))+
  theme(strip.text.x = element_text(size = 8))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"))+
  theme(legend.position=c(.85,.4),legend.text=element_text(face="bold",size=8),
        legend.key=element_blank())
 
#1.1 Variogram plots for different timescales at 1-5 mm/hr
vg_1to5=read.csv("./1_Data/vg_1to5.csv")
load("./1_Data/vgm_1to5.RData")

vg_1to5
vg_1to5 = vg_1to5[-c(1,5,6,7)]
vg_1to5$log_np = log(vg_1to5$np)
vg_1to5$ts = factor(vg_1to5$ts, levels=unique(vg_1to5$ts)) #to adjust the order of facet_wrap in ggplot

fitted = seq(0.01, max(vg_1to5$dist), length = 100)

gamma_1to5=NULL
abc=seq(1,5,2)
for ( i in abc) {
  gamma = variogramLine(vgm_1to5[i:(i+1),], dist_vector = fitted)$gamma
  gamma_1to5=cbind(gamma_1to5,gamma)
}
gamma_1to5=c(gamma_1to5)

fitted=cbind(rep(fitted,6),gamma_1to5)
fitted=as.data.frame(fitted)
colnames(fitted)=c("dist","gamma")
fitted$ts=rep(c("30 min","60 min","120min"),each=100)


fitted$ts = factor(fitted$ts, levels=unique(fitted$ts))  #to adjust the order of facet_wrap in ggplot

theme_set(theme_bw())

ggplot(vg_1to5, aes(x = dist, y = gamma)) + 
  #geom_point(aes(shape=ts),size=3)+
  scale_y_continuous(limits = c(0, 2))+
  geom_line(data=fitted,aes(linetype=ts),size=1)+
  labs(x="Distance [m]", y="Semivariance [-]")+
  scale_linetype_manual(values=c("solid","dashed","dotted"))+
  theme(legend.position="none",legend.title=element_blank(),
        legend.text=element_text(face="bold",size=12),
        legend.key=element_blank(),legend.key.width=unit(3,'cm'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=12,face="bold"),
        axis.title=element_text(size=12,face="bold"))+
  annotate("text", label = "120 min", x = 400, y = 1.6, size=4)+
  annotate("text", label = "60 min", x = 400, y = 1.5, size=4)+
  annotate("text", label = "30 min", x = 400, y = 1.4, size=4)

#2. plotting kriging prediction with confidance bands

#reading and formatting data for ggplot
kriging.comb=read.csv("./3_Results/kriging.comb.csv")
kriging.comb$Date = as.POSIXct(kriging.comb$Date,format='%Y-%m-%d %H:%M:%S')
kriging.comb$ts = factor(kriging.comb$ts, levels=unique(kriging.comb$ts)) #to adjust the order of facet_wrap in ggplot

#plotting
theme_set(theme_bw())
ggplot(kriging.comb, aes(x = Date, y = pred)) + 
  geom_point(size=1.5)+
  geom_line(linetype="dotted") +
  #geom_ribbon(aes(ymax = CI50up, ymin = CI50low), fill = 504, alpha =0.2 )+
  geom_ribbon(aes(ymax = CI95up, ymin = CI95low), fill = 505, alpha =0.35 )+
  facet_wrap(~ts)+
  theme(legend.position="bottom",legend.title=element_blank())+
  scale_fill_manual("",values=c("grey","grey70"),labels=c("95% confidence interval","50% confidence interval"))+
  labs(x="Time [hh:mm]", y="Predicted rainfall [mm/h]")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10,face="bold"))+
  theme(strip.text.x = element_text(size = 8))


#3. plotting peaks pred and variance

#reading and formatting data for ggplot
event.peaks=read.csv("./3_Results/eventPeaks.csv")
event.peaks$ts=rep(c("2 min","5 min", "15 min","30 min"),each=19)
event.peaks$index=rep(1:19,4)
event.peaks$ts = factor(event.peaks$ts, levels=unique(event.peaks$ts)) #to adjust the order of facet_wrap in ggplot
event.peaks$cod=signif(event.peaks$sd/event.peaks$pred*100, 2)

#plotting
theme_set(theme_bw())
ggplot(event.peaks, aes(x = index, y = pred, ymax = CI95up, ymin = CI95low)) +
  geom_pointrange(size=0.35)+
  #geom_point(aes(colour=cod*100))+
  facet_wrap(~ts)+
  scale_x_continuous(breaks = 1:19)+
  scale_y_continuous(breaks = seq(0,100,20))+
  theme(panel.grid.minor.x=element_blank())+
  labs(x="Event ID", y="Peak rainfall [mm/h]")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10, face="bold"))+
  theme(strip.text.x = element_text(size = 8))


#plotting: Additional work for referee Comments
theme_set(theme_bw())
ggplot(event.peaks, aes(x = index, y = pred, label = cod)) +
  #geom_pointrange(size=0.35)+
  geom_text(size=2.5,hjust = 0.5, vjust = -1)+
  geom_point(size=1.5)+
  facet_wrap(~ts)+
  scale_x_continuous(breaks = 1:19)+
  scale_y_continuous(limits= c(0, 110),breaks = seq(0,100,20))+
  theme(panel.grid.minor.x=element_blank())+
  labs(x="Event ID", y="Peak rainfall [mm/h]")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10, face="bold"))+
  theme(strip.text.x = element_text(size = 8))


#4. plotting coefficient of variation
theme_set(theme_bw())
means= data.frame(ts=rep(c("Mean - 2 min","Mean - 30 min"),2),x1=rep(c(0,10),each=2),x2=rep(c(10,100),each=2),
                  y1=c(6.62,1.74,3.06,1.21),y2=c(6.62,1.74,3.06,1.21))
ggplot(event.peaks, aes(x = pred, y = cod)) +
  geom_point(aes(shape=ts),size=2.5)+
  scale_shape(solid = FALSE)+
  scale_x_log10() +
  geom_vline(xintercept = 10, linetype = "longdash")+
  geom_segment(data=means, aes(x=x1,xend=x2,y=y1,yend=y2,linetype = ts))+
  annotate("text", label = "< 10 mm/h", x = 7, y = 15,size=4.5)+
  annotate("text", label = "> 10 mm/h", x = 14.5, y = 15,size=4.5)+
  theme(panel.grid.minor.x=element_blank())+
  labs(x="Predicted peak rainfall [mm/h]", y="Coefficient of variation [%]")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10, face="bold"),
        legend.position = c(.85, .75),legend.title=element_blank(),
        legend.key = element_blank(),legend.text.align=0)+
        scale_color_manual(labels=c("min","max"),colour=C("blue","orange"))

#plotting: Additional work for referee Comments
theme_set(theme_bw())
means= data.frame(ts=rep(c("Mean - 2 min","Mean - 30 min"),2),x1=rep(c(0,10),each=2),x2=rep(c(10,100),each=2),
                  y1=c(6.62,1.74,3.06,1.21),y2=c(6.62,1.74,3.06,1.21))
ggplot(event.peaks, aes(x = pred, y = cod)) +
  geom_point(aes(shape=ts),size=2.5)+
  scale_x_log10() +
  geom_vline(xintercept = 10, linetype = "longdash")+
  geom_segment(data=means, aes(x=x1,xend=x2,y=y1,yend=y2,linetype = ts))+
  annotate("text", label = "< 10 mm/h", x = 7, y = 15,size=4.5)+
  annotate("text", label = "> 10 mm/h", x = 14.5, y = 15,size=4.5)+
  theme(panel.grid.minor.x=element_blank())+
  labs(x="Predicted peak rainfall [mm/h]", y="Coefficient of variation [%]")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10, face="bold"),
        legend.position = c(.85, .75),legend.title=element_blank(),
        legend.key = element_blank(),legend.text.align=0)+
  scale_color_manual(labels=c("min","max"),colour=C("blue","orange"))

#5. plotting prediction vs variance
#reading and formatting data for ggplot
allevents.5min=read.csv("./1_Data/allevents_5min.csv")
allevents.5min$ts="5 min"
allevents.5min$CV=allevents.5min$sd/allevents.5min$pred*100

allevents.30min=read.csv("./1_Data/allevents_30min.csv")
allevents.30min$ts="30 min"
allevents.30min$CV=allevents.30min$sd/allevents.30min$pred*100
#allevents.5min$Date = as.POSIXct(allevents.5min$Date,format='%Y-%m-%d %H:%M:%S')
allevents.com=rbind(allevents.5min,allevents.30min)

#plotting only one time scale
theme_set(theme_bw())
ggplot(allevents.5min, aes(x = pred, y = CV))+
  geom_point(size=1.5)+
  scale_x_log10(breaks=c(1,5,10,100)) +
  geom_vline(xintercept = c(5,10), linetype = "longdash")+
  labs(x="Predicted rainfall [mm/h]", y="Coefficient of variation [%]")+
  #scale_x_discrete(limits=c(1,5,10,100))+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10, face="bold"))
#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#plotting multiple timescale
theme_set(theme_bw())
allevents.com$CV=allevents.com$sd/allevents.com$pred*100
ggplot(allevents.com, aes(x = pred, y = CV, colour=ts))+
  geom_point(size=1.5)+
  scale_x_log10() +
  labs(x="Predicted rainfall [mm/h]", y="Coefficient of variation [%]")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10, face="bold"))
#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#5. combining events
aaa.all=NULL
for (i in 1:14){
  aaa=read.csv(paste0(i,"_5.csv"))
  aaa.all=rbind(aaa.all,aaa)
}
