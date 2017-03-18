library(zoo)
library(abind)
library(sp)

###reading from CSV files and extracting data for a particular event###
eventTS=function(directory,filename,eventStartDatetime,eventEndDatetime){
  tsComplete=read.zoo(paste0(directory,filename),sep = ",", header = TRUE, 
                      index = 1, tz = "", format = "%d/%m/%Y %H:%M")
  tsEvent=window(tsComplete, start=as.POSIXct(eventStartDatetime),
                 end=as.POSIXct(eventEndDatetime))
  return(tsEvent)
}
###The following function is appplicable when the data is rainfall intensity

###intensity time series for a given time scale###
IntTS=function(eventTS,timeScale=60){
  tsInt=rollapply(eventTS,width=timeScale,FUN=mean,by=timeScale,align="right")
  return(tsInt)
}


###The following functions are applicable when rainfall data is cumulative 

###extracting time series for a particular timescale###
extractTS=function(tsEvent,timeScale) {
  aaa=index(tsEvent)
  bbb=cut(aaa,breaks=timeScale)
  ccc=levels(bbb)
  ddd=lapply(ccc,function(x) tsEvent[x==aaa])
  eee=do.call("rbind",ddd)
  return(eee)
}

###Cumulative to step accumulated###
cumToStAcc=function(tsCum) {
  tsStAcc=diff(tsCum)
  return(tsStAcc)
}

###step accumulated to intenstity###
stAccToInt=function(tsStAcc, timeScale) {
  tsInt=tsStAcc*(60/as.numeric(strsplit(timeScale," ")[[1]][[1]]))
}

###CSV to timescaled intensities (All the above function are combined into one)###
csvToInten=function(directory,filename,eventStartDatetime,eventEndDatetime,timeScale) {
  tsEvent=eventTS(directory,filename,eventStartDatetime,eventEndDatetime)
  tsTimescale=extractTS(tsEvent,timeScale)
  tsTimescaleStAcc=cumToStAcc(tsTimescale)
  tsTimescaleInt=stAccToInt(tsTimescaleStAcc,timeScale)
  return(tsTimescaleInt)
}

###coordinates to distance matrix#####
co_to_dis=function(coord) {
  aaa=as.matrix(coord)
  bbb=spDists(aaa)
  return(bbb)
}

###averaging of two rainguages####
get_avr=function(tsInt,gauMean=TRUE,rgNo=seq(1,15,2)) {
  if(gauMean==TRUE) {
    #averaging of two rainguages
    avrIntAll=NULL
    for (i in rgNo) {
      avrInt=apply(tsInt[,i:(i+1)],1,mean)
      avrIntAll=cbind(avrIntAll,avrInt)
    } 
  } else {
    #using single rain gauges
    avrIntAll=tsInt[,rgNo]
  }
  return(avrIntAll)
}

#zoo to data frame
zoo.to.data.frame <- function(x, index.name="Date") {
  stopifnot(is.zoo(x))
  xn <- if(is.null(dim(x))) deparse(substitute(x)) else colnames(x)
  setNames(data.frame(index(x), x, row.names=NULL), c(index.name,xn))
}


#bin
#***to convert time series to 'xts'(zoo) object***
# library(xts)
# converttoZoo=function(timeseries,timeColumn=1) {
#   xtsTimeser<-xts(timeseries[,-1],order.by=timeseries[,timeColumn])
#   return(xtsTimeser)
# }