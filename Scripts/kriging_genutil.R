#ver.2 - samp.percentage added.

#susetting data according to thresholds and format it for variogram prediction
geostat_data=function(data,coord,lthres=4,uthres=12, samp.per=100) {
  n=ncol(data)-1
  rowMean=rowMeans(data[,-1])
  data.sub1=data[which(rowMean>lthres),]
  
  rowMean=rowMeans(data.sub1[,-1])
  data.sub2=data.sub1[which(rowMean<uthres),]
  
  #sampling
  sel.rows=sample(1:nrow(data.sub2),samp.per*nrow(data.sub2)/100)
  data.sub3=data.sub2[sort(sel.rows),]
  #preparing the data for kriging
  data.rain=data.sub3
  data.rainZ=data.rain[,-1]
  data.rainZ=as.vector(t(data.rainZ))
  data.rainT=data.rain[,1]
  data.rainT=rep(data.rainT,each=n)
  
  coord.x=rep(coord[,1],nrow(data.rain))
  coord.y=rep(coord[,2],nrow(data.rain))
  
  data.kriging=data.frame(matrix(NA,nrow=length(data.rainZ),
                                 ncol=4))            
  data.kriging[,1]=coord.x
  data.kriging[,2]=coord.y
  data.kriging[,3]=data.rainT
  data.kriging[,4]=data.rainZ
  colnames(data.kriging)=c("x","y","time","int")
  return(data.kriging)
 }



