library(calibrate)
library(quantreg)

setwd('..')
source("./2_Scripts/weightRain.R")

# © WaPUG 2002 WaPUG Code of Practice for the Hydraulic Modelling of Sewer Systems Version 3.001
# Page 38 of 69
# i)    cond1: The total depth should be greater than 5 mm (thres_cond1).
# ii)   cond2: The range of total durations of the rainfall should vary. Ideally one storm should
#       have a duration (thres_cond2 in min) equal to ½ of the time of concentration (Tc) of the system, one
#       storm equal to Tc, and the third more than twice the time Tc.
# iii)  cond3: The rainfall intensity should be greater than 6 mm/hour (thres_cond3) for more than 4 minutes.
#       (Have to check with timestep of the timeseries selected. The default is set to 6, assuming the timestep is 4 min)
# iv)   cond4: The period between events should be sufficient for the flow to return to dry weather
#       conditions. Not conidered here, check it manually from filtered events

avr.int=5 #averaging interval of rainfall intensity (min)

#reading rainfall time series and extracting the events using timeSerToEvents
rainfall=read.csv(paste0("./1_Data/int_brad_",avr.int,".csv"))
rainfall$avr=rowMeans(rainfall[,2:9])
rainfall.zoo=read.zoo(rainfall,sep = ",", header = TRUE, index = 1,tz = "",format = "%d/%m/%Y %H:%M")
allEvents=timeSerToEvents(timeSer=rainfall[,c(1,10)],thres_cond1=5,thres_cond2=30,
                          thres_cond3=5,plot=TRUE) #keep it same for all avr interval inclduing 2 min
tsInt=rainfall.zoo[,-9] # change this for 2 min avr interval


#1. culative to intensity
#start
eventData=read.csv("./1_Data/eventData.csv",header=TRUE)
k=c(17)  #selEvents
j=c("2 min") #timeScale

#raingauges
gauMean=TRUE #if the mean value of pair raingauge should be used or not
rgNo=seq(1,15,2) #raingauge numbers  (doesn't make any dofferent if TRUE is selected above) 

#reading from csv files and extracting data for the event date and time
fileName=paste(eventData[k,2])
eventStart=as.POSIXct(eventData[k,3],tz="",format = "%d/%m/%Y %H:%M")
eventEnd=as.POSIXct(eventData[k,4],tz="",format = "%d/%m/%Y %H:%M")
eventDur=as.numeric(difftime(eventEnd, eventStart, units="hours"))
tsInt=csvToInten("./1_Data/",fileName,eventStart,eventEnd,timeScale=j)

tsInt=as.data.frame(tsInt)
write.csv(tsInt, "./1_Data/int_1min_2013c.csv")

#averaging of two rainguages for a particular timescale
tsInt=read.zoo("./1_Data/int_1min_2013.csv",sep = ",", header = TRUE, 
               index = 1,tz = "", format = "%d/%m/%Y %H:%M")
rgNo=seq(1,15,2)
ts=1#time sclae in min
fileName="int_brad2013_"

avrIntAll=IntTS(tsInt,ts)
avrIntAll=get_avr(avrIntAll,gauMean=TRUE,rgNo)
avrIntAll=as.data.frame(avrIntAll)
colnames(avrIntAll)=1:length(rgNo)
write.csv(avrIntAll,paste0("./1_Data/",fileName,ts,".csv"))
#end


#2. extracting events at different timescales except for 2min
int.events=c(1:19)
for (i in int.events) {
  
  aaa=window(tsInt,start=allEvents[[2]][1,i],end=allEvents[[2]][2,i])
  
  aaa.new=aaa[c(-1,-nrow(aaa)),]
  
  bbb=rollapply(aaa.new,3,mean,by=3,partial=TRUE) # for 15 min
  if (nrow(aaa.new)>nrow(bbb)*3) {   #to add last row
    ab=colSums(tail(aaa,nrow(aaa[-1,])-nrow(bbb)*3))/3
    ab.last=zoo(rbind(ab),end(bbb) + 15*60)
    bbb=rbind(bbb,ab.last)
    }
  
  ab.end=zoo(rbind(rep(0,8)),end(bbb) + 15*60)
  ab.start=zoo(rbind(rep(0,8)),start(bbb) + -15*60)
  bbb=rbind(ab.start, bbb, ab.end)
  
  ccc=rollapply(aaa.new,6,mean,by=6) # for 30 min
  if (nrow(aaa.new)>nrow(ccc)*6) {   #to add last row
    ac=colSums(tail(aaa,nrow(aaa[-1,])-nrow(ccc)*6))/6
    ac.last=zoo(rbind(ac),end(ccc) + 30*60)
    ccc=rbind(ccc,ac.last)
  }
  
  ac.end=zoo(rbind(rep(0,8)),end(ccc) + 30*60)
  ac.start=zoo(rbind(rep(0,8)),start(ccc) + -30*60)
  ccc=rbind(ac.start, ccc, ac.end)
  
  write.csv(as.data.frame(aaa),paste0("./3_Results/",i,"_",5,".csv"))
  write.csv(as.data.frame(bbb),paste0("./3_Results/",i,"_",15,".csv"))
  write.csv(as.data.frame(ccc),paste0("./3_Results/",i,"_",30,".csv"))
  plot(aaa,plot.type="single",col=rainbow(8),main=i,
       xlab=paste0("Time (",as.Date(index(aaa[1,])),")"),
       ylab="Intensity [mm/hr]")
  grid(nx=NA,ny=NULL)
}

# 3. extracting events for 2 min
int.events=c(1:19)
for (i in int.events) {
  ddd=window(tsInt,start=allEvents[[2]][1,i],end=allEvents[[2]][2,i])
  write.csv(as.data.frame(ddd),paste0("./3_Results/",i,"_",2,".csv"))
  plot(ddd,plot.type="single",col=rainbow(8),main=i,
       xlab=paste0("Time (",as.Date(index(aaa[1,])),")"),
       ylab="Intensity [mm/hr]")
  grid(nx=NA,ny=NULL)
}


#4. just plotting measurments from all the stations 
int.events=c(1:14)
for (i in int.events) {
  aaa=window(tsInt,start=allEvents[[2]][1,i],end=allEvents[[2]][2,i])
  plot(aaa,plot.type="single",col=rainbow(8),main=i,
       xlab=paste0("Time (",as.Date(index(aaa[1,])),")"),
       ylab="Intensity [mm/hr]")
  grid(nx=NA,ny=NULL)
}

##just used for extra work 
##start
bbb=NULL
nnn=NULL
i=c(6)
start=allEvents[[2]][1,i]
end=allEvents[[2]][2,i]

aaa=window(tsInt,start=start,end=end)
write.csv(aaa,paste0(i,"_",avr.int,".csv"))
plot(aaa,plot.type="single",col=rainbow(8),xlab="Time",ylab="Intensity [mm/hr]")
grid(nx=NA,ny=NULL)
##end

#5. summary table of the events                                                                                                       

#dates and duration
int.events=c(1:14)
event.times=allEvents$`Time of the events`
event.dates=as.Date(sapply(event.times[1,],function (x) as.Date(x)))
event.durs=apply(event.times,2,function (x) difftime(x[2],x[1],tz="", "hours"))

#calculating avergae intensity
event.int.all=allEvents$`Event timeseries`
event.int=rapply(event.int.all, function (x) mean(x[,1]))

#calulating summary statistics of peaks of the events
event.peaks=rapply(event.int.all, function (x) index(x[which(x[,1]==max(x[,1]))]),
                   how="list")
peak.data=rapply(event.peaks,function (x) rainfall.zoo[x,])
peak.data=as.data.frame(matrix(peak.data,ncol=9,byrow=TRUE))
colnames(peak.data)=c(1:8,"mean")
peak.data$sd=apply(peak.data[,1:8],1,sd)
peak.data$max=apply(peak.data[,1:8],1,max)
peak.data$min=apply(peak.data[,1:8],1,min)
peak.data$range=c(diff(apply(peak.data[,1:8],1,range)))
peak.data$per.range=peak.data$range/peak.data$mean*100


#table
sum.tab=data.frame(int.events,event.dates,event.durs,event.int,peak.data$mean,peak.data$sd,
                     peak.data$max,peak.data$min,peak.data$range,peak.data$per.range)
sum.tab
write.csv(sum.tab,"./3_Results/summary_table.csv")

#6. daily rainfall
library(xts)

cum.rf=read.zoo("./1_Data/complete_min_2012.csv",sep = ",", header = TRUE, 
               index = 1,tz = "", format = "%d/%m/%Y %H:%M")
rgNo=seq(1,15,2)
cum.rf=get_avr(cum.rf,gauMean=TRUE,rgNo)
cum.rf=as.xts(cum.rf)
daily.cum=cum.rf[.indexhour(cum.rf) %in% c(23) & .indexmin(cum.rf) %in% c(59)]
daily.cum=diff(daily.cum)
avr.daily.cum=zoo(rowMeans(daily.cum), time(daily.cum))
# plot(avr.daily.cum,type="l",lwd="2", xlab="",ylab="Daily rainfall [mm]",
#      xlim=as.POSIXct(c("2012-04-01","2012-09-01")),ylim=c(0,45))
# text(as.POSIXct("2012-04-01"), 35,  labels=c("2012"),cex=2)
# abline(v=as.POSIXct(event.dates+1),lty="dashed")
avr.daily.cum.plot=fortify(avr.daily.cum)
ggplot(avr.daily.cum.plot, aes(x = Index, y = avr.daily.cum)) +
  geom_line(size=1)+
  scale_y_continuous(limits = c(0, 45))+
  scale_x_datetime(limits = as.POSIXct(c("2012-04-01","2012-09-01")))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_vline(xintercept=as.numeric(as.POSIXct(c(event.dates+1))),linetype=2)+
  labs(x="", y="Daily rainfall [mm]")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  annotate("text", label = "2012", x = as.POSIXct("2012-04-08"), y = 40, size=7)

cum.rf=read.zoo("./1_Data/complete_min_2013.csv",sep = ",", header = TRUE, 
                index = 1,tz = "", format = "%d/%m/%Y %H:%M")
rgNo=seq(1,15,2)
cum.rf=get_avr(cum.rf,gauMean=TRUE,rgNo)
cum.rf=as.xts(cum.rf)
daily.cum=cum.rf[.indexhour(cum.rf) %in% c(23) & .indexmin(cum.rf) %in% c(59)]
daily.cum=diff(daily.cum)
avr.daily.cum=zoo(rowMeans(daily.cum), time(daily.cum))
avr.daily.cum.plot=fortify(avr.daily.cum)
ggplot(avr.daily.cum.plot, aes(x = Index, y = avr.daily.cum)) +
  geom_line(size=1)+
  scale_y_continuous(limits = c(0, 45))+
  scale_x_datetime(limits = as.POSIXct(c("2013-04-01","2013-09-01")))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_vline(xintercept=as.numeric(as.POSIXct(c(event.dates+1))),linetype=2)+
  labs(x="", y="Daily rainfall [mm]")+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=9,face="bold"))+
  annotate("text", label = "2013", x = as.POSIXct("2013-04-08"), y = 40, size=4.5)

#7. extracting peak
events=c(1:19)
ts=c(2,5,15,30)

for (j in ts) {
  peaks_comb=NULL
  for (i in events) {
    xxx=read.csv(paste0("./3_Results/",i,"_",j,".csv"))
    xxx$mean=rowMeans(xxx[,-1])
    peaks=xxx[which(xxx$mean==max(xxx$mean)),]
    peaks_comb=rbind(peaks_comb,peaks)
  }
  write.csv(peaks_comb,paste0("./3_Results/peaks_",j,".csv"))
  read.csv(paste0("./3_Results/peaks_",j,".csv"))[,-c(1,11)]
  
}


#8.  1 min intensity to 't'min intensity
t=30
eve.ts=tsInt=read.zoo("./1_Data/int_brad_1.csv",sep = ",", header = TRUE, 
                      index = 1,tz = "", format = "%d/%m/%Y %H:%M")
int.t=IntTS(eve.ts,t)
write.csv(as.data.frame(int.t),paste0("./1_Data/int_brad_",t,".csv"))



#bin
#saving a specific event
# S <- split(seq_len(nrow(bbb)), rep.int(seq_along(nnn), nnn))
# S <-lapply(S, function(x) bbb[x, , drop = FALSE])
# write.csv(as.data.frame(S[[1]]),"temp.csv")


if (length(a)>length(b)*3){
  c=sum(tail(a,length(a)-length(b)*3))/3
  }




