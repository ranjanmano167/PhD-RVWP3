library(rgdal)
library(spatstat)
library(ggplot2)
library(grid)
library(gridExtra)


source("./2_Scripts/gen_util.R")

###thiessen weights###
#rgCoord: raingauge coordinates (x,y respectively in matrix form)
#rgId: which raingauges should be conidered (row numbers of rgCoord)
#bdryBas: coordinates of the catchment bounday 
#         (x,y respectively in matrix form with headers 'x' and 'y')
#bdryThi: coordinates of thiessen polygon bounday 
#         (x,y respectively in matrix form with headers 'x' and 'y')

thiWeights=function(rgCoord,rgId=c(1:8),bdryBas,bdryThi) {
  rgSel=rgCoord[rgId,]
  W1 <- owin(poly=bdryThi)
  w2<- owin(poly=bdryBas)
  X <- as.ppp(rgSel,W=W1)
  #plot(X) 
  dX <- dirichlet(X)
  plot(dX,main="Thiessen Polygons")
  intSec=lapply(tiles(dX),function (x) intersect.owin(w2,x,fatal=FALSE))
  #intSec[[which(sapply(intSec,is.null))]]=owin(c(0,0),c(0,0))
  clr=rainbow(length(intSec))
  for (i in 1:length(intSec)){
    if (is.null(intSec[[i]])==TRUE){
      intSec[[i]]=owin(c(0,0),c(0,0))
    }
    plot(intSec[[i]],add=TRUE, main=NULL, col=clr[i])
  }
  cents <- as.data.frame(t(sapply(tiles(dX),centroid.owin)))
  text(cents,labels=rgId)
  thiArea=sapply(intSec,area.owin)
  names(thiArea)=rgId
  thiWt=thiArea/sum(thiArea)
  names(thiWt)=rgId
  result=list(rgId,thiArea,thiWt)
  return(result)
}

#owin object to polygon (https://stat.ethz.ch/pipermail/r-sig-geo/2009-May/005781.html)
owin2Polygons <- function(x, id="1") {
  stopifnot(is.owin(x))
  x <- as.polygonal(x)
  closering <- function(df) { df[c(seq(nrow(df)), 1), ] }
  pieces <- lapply(x$bdry,
                   function(p) {
                     Polygon(coords=closering(cbind(p$x,p$y)),
                             hole=is.hole.xypolygon(p))  })
  z <- Polygons(pieces, id)
  return(z)
}

thiWeights=function(rgCoord,rgId=c(5,7),bdryBas,bdryThi) {
  rgSel=rgCoord[rgId,]
  W1 <- owin(poly=bdryThi)
  w2<- owin(poly=bdryBas)
  X <- as.ppp(rgSel,W=W1)
  #plot(X) 
  dX <- dirichlet(X)
  plot(dX,main="Thiessen Polygons")
  intSec=lapply(tiles(dX),function (x) intersect.owin(w2,x,fatal=FALSE))
  #intSec[[which(sapply(intSec,is.null))]]=owin(c(0,0),c(0,0))
  clr=rainbow(length(intSec))
  for (i in 1:length(intSec)){
    if (is.null(intSec[[i]])==TRUE){
      intSec[[i]]=owin(c(0,0),c(0,0))
    }
    plot(intSec[[i]],add=TRUE, main=NULL, col=clr[i])
  }
  cents <- as.data.frame(t(sapply(tiles(dX),centroid.owin)))
  text(cents,labels=rgId)
  points(rgSel)
  thiArea=sapply(intSec,area.owin)
  names(thiArea)=rgId
  thiWt=thiArea/sum(thiArea)
  names(thiWt)=rgId
  
  #calculating distribution factor
  thiPoly=lapply(intSec,owin2Polygons)
  thiCoor=lapply(thiPoly,fortify)
  
  farDistAll=NULL
  rgDistAll=NULL
  yyy=0
  for (i in 1:length(thiCoor)) {
    #calculating maxmimum distance from centroid
    longs=thiCoor[[i]][,1]
    lats=thiCoor[[i]][,2]
    dist=sqrt((longs-cents[i,1][[1]])^2+(lats-cents[i,2][[1]])^2)
    farDist=max(dist)
    farDistAll=cbind(farDistAll,farDist)
    
    #calculating distance between centroid and raingauge
    if (length(thiCoor)==1) {
      rgDist=sqrt((rgSel[1]-cents[i,1][[1]])^2+(rgSel[2]-cents[i,2][[1]])^2)
    } else {
      rgDist=sqrt((rgSel[i,1]-cents[i,1][[1]])^2+(rgSel[i,2]-cents[i,2][[1]])^2)
    }
    
    rgDistAll=cbind(rgDistAll,rgDist)
    
    xxx=rgDist/farDist*thiArea[i]
    yyy=sum(yyy,xxx)
  }
  
  DisFactor=1-yyy/sum(thiArea)
  
  result=list(rgId,thiArea,thiWt,DisFactor)
  return(result)
}

###weighted rainfall###
#rainData:rainfall timeseries wwhere each column represent different raingauges
#weights: weights as a numeric vector  (can be less than number of columns in 
#         rainData, see below)
#rgId:raingauge numbers for which the weights are assigned above 
#     (assumes the rest of the weights are zero)
weightRain=function(rainData, weights,rgId=c(1:8)){
  weightsAdj=rep(0,ncol(rainData))
  weightsAdj[rgId]=weights
  
  bbb=NULL
  for (i in rgId){
    aaa=rainData[,i]*weightsAdj[i]
    bbb=cbind(bbb,aaa)
  }
  ccc=rowSums(bbb)
  return(ccc)
}

###Selecting rainfall events, which satisfies WAPUG guidelines, from a timeseries ###

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

#this function can be used to extract events from any timeseries. For example it could be used for 
#discharge timeseries and in that case the input are
#thres_cond1: cumulative discharge threshold i.e. only the events which  result in 
#             total ruoff of thres_hold1 or above will be considered (units m3/h)
#thres_cond2: Duration of the runoff event i.e. only the events which  lasts for 
#             thres_holds or longer will be considered (units min)
#thres_cond3: minimum runoff threshold  i.e. only the events which exceeds this 
#             runoff threshold will be considered (units min)



#timeSer=rainfall timeseries (Date and time(%d/%m/%Y %H:%M), Intensties(mm/hr))
#plot=TRUE : if the filetered events should be plotted (Intensity vs Time)

timeSerToEvents=function(timeSer,thres_cond1=5,thres_cond2=30,thres_cond3=6,
                         plot=TRUE) {
  timeSer$Index=1:nrow(timeSer)
  timeSerZoo=read.zoo(timeSer,header = TRUE, index = 1, tz = "", 
                      format = "%d/%m/%Y %H:%M")
  timeStep=as.numeric(difftime(index(timeSerZoo)[2],index(timeSerZoo)[1], 
                               units="min"))
  timeSerZoo1=timeSerZoo[timeSerZoo[,1]>0,]
  timeSerZoo1$d=c(1,diff(timeSerZoo1[,2]))
  events=split(timeSerZoo1[timeSerZoo1[,3]==1,], 
               cumsum(timeSerZoo1[,3]!=1)[timeSerZoo1[,3]==1])
  names(events)=1:length(events)
  events=lapply(events,function (x) x[,c(1,2)])
  eventsCom=lapply(events, function(x) tryCatch(rbind(timeSerZoo[(as.numeric(x[,2])[1]-2):
                                                                   (as.numeric(x[,2])[1]-1),],x,
                                                      timeSerZoo[(as.numeric(tail(x,1)[,2])+1),]),
                                                error=function(e) x))
  
  #filtering events based on conditioned defined
  cond1=sapply(eventsCom,function (x) (sum(x[,1])*(timeStep/60))>thres_cond1)
  eventsF1=eventsCom[cond1]
  cond2=sapply(eventsF1,function (x) (nrow(x)-1)*timeStep>thres_cond2)
  eventsF12=eventsF1[cond2]
  cond3=sapply(eventsF12,function (x) (max(x[,1])>thres_cond3))                                           
  eventsF123=eventsF12[cond3]
  
  sel_events=eventsF123
  names(sel_events)=1:length(sel_events)
  
  
  if (plot==TRUE) {
    lapply(seq_along(sel_events), function(idx) plot(sel_events[[idx]][,1],
                                                     main=paste0("event-",names(sel_events[idx])),
                                                     xlab=paste0("Time (",as.Date(index(sel_events[[idx]][1,])),")"),
                                                     ylab="Intensity [mm/hr]"))
  }
  
  event_time=lapply(sel_events,function(x) index(c(x[1,],x[nrow(x),])))
  event_time=as.data.frame(event_time)
  rownames(event_time)=c("start","end")
   
  result=list(timeStep,event_time,sel_events)
  names(result)=c("Time step (min)","Time of the events", "Event timeseries" )
  return (result)
  
}

###producing runoffpeak plot and subplot of corrosponding rainfall event###
#runoffEvents=runoff events from number of different raingauge combinations 
#rainfallEvent=corrosponding rainfall event
#numComb=number of runs corrosponding to number of raingauges considered

peakPlot=function (runoffEvents,rainfallEvent,numComb=c(8,28,56,70,56,28,8,1)) {
  timeStep=as.numeric(difftime(index(rainfallEvent)[2],index(rainfallEvent)[1], 
                               units="min"))
  
  
  #plotting all runoff events
#   runoffPlot=autoplot.zoo(runoffEvents,facet=NULL,linetype =NULL)+ xlab(paste0("Time (",as.Date(index(runoffEvents[1,])),")"))+
#     ylab("Runoff [m3/s]")+ theme(legend.position = "none")+l
#   print(runoffPlot)
  clr=rainbow(ncol(runoffEvents))
  plot(runoffEvents[,1],ylim=c(0, max(runoffEvents)),xlab=paste0("Time (",as.Date(index(runoffEvents[1,])),")"),
       ylab="Runoff [cum]") #
  allPlots=sapply(seq_along(2:ncol(runoffEvents)),function(x) lines(runoffEvents[,x],col=clr[x]))
  
  #extracting the peak runoff
  maxAll=apply(runoffEvents,2,max)
  maxvsnum=matrix(c(rep(1:length(numComb),numComb),maxAll),ncol=2,byrow=FALSE)
  
  #calculating max percentage different 
  dddd=NULL
  ffff=NULL
  for (k in 1:length(numComb)){
    aaaa=range(maxvsnum[which(maxvsnum[,1]==k),2])
    bbbb=round(diff(aaaa)/mean(aaaa)*100,1)
    cccc=cbind(k,bbbb)
    dddd=rbind(dddd,cccc)
    eeee=max(aaaa)
    ffff=rbind(ffff,eeee)
  }
  
  #plotting peakrunoff vs num.of.raingauges (add theme_bw(base_size = 18) after plot1 for b&w theme)
  options(warn=-1)
  plot1=qplot(maxvsnum[,1],maxvsnum[,2]*1000,xlab="Number of raingauges",
              ylab="Peak runoff [L/s]")
  mainplot=plot1+
    geom_text(label=maxvsnum[,2]*1000)+#,x=maxvsnum[,1],y=maxvsnum[,2]*1000)+
    scale_x_continuous(breaks=c(1:8))+ 
    theme(panel.grid.minor = element_blank())+
    annotate(geom="text",x=dddd[,1],y=ffff*1000,label=paste0(dddd[,2],"%"),
             hjust = -0.5,col="red")+ 
    annotate(geom="text",x=7,y=max(maxvsnum[,2]*1000),
             label="x% - Maximum Pecentage Difference",col="red")
  
  #plotting rainfall event
  duration=seq(from=0,by=timeStep,length.out=length(rainfallEvent))
  plot2=qplot(duration,coredata(rainfallEvent),
                  geom=c("line"),xlab=NULL,ylab=NULL)
  subplot=plot2+theme_bw()+ggtitle("Rainfall Intensity [mm/h]\n Vs Duration [min]")+
    theme_classic()+theme(plot.title = element_text(face="bold",lineheight=.8,size=12))+
    annotate(geom="point",x=duration[which(rainfallEvent==max(rainfallEvent))],
                             y=max(rainfallEvent), colour ="red", size = 3) 
    
    
  vp= viewport(x = 0.82, y = 0.25, width = 0.28, height = 0.25)
  print(mainplot)
  print(subplot,vp=vp)
return(list(maxAll,max(coredata(rainfallEvent))))
 
}