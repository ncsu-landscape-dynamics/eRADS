
## Define Function Input Variable ##
# inf is the infested file in raster format, 
# buffer is the buffer size around each treatment pixel
# width is the buffer size within which the total number of uninfested host would be derived
# ply is the polygon converted from the infested raster, the main use of ply is to exclude the infested host
# pixelArea is the area of each pixel

########### Method 1 -- Randomly Select Infested Pixel for Treatment ################
## Randomly select infested pixels based on budget and cost per unit ##

random_pixel=function(inf,buffer,budget,cost_per_meter_sq){
  
  area=budget/cost_per_meter_sq
  
  inf2=inf
  k<- which(inf2[]==0)
  inf2[k]=NA
  
  inf2[inf2==0]=NA
  
  rd=sampleRandom(inf2,size=n,na.rm=T,asRaster=T)
  rd_ply=rasterToPolygons(rd,dissolve=F)
  ply_bf=buffer(rd_ply,width=buffer,dissolve=F)
  
  n=floor(area/pixelArea)
  
  ply_bf$Cumu_Area=0
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply_bf[ply_bf$Cumu_Area> area,]
  
  if (df>0){
    crds=gCentroid(nontr[1,])
    crds_bf=buffer(crds,width=sqrt(df/pi))
    treatment=gUnion(treatment,crds_bf)
  }
  
  treatmentRa=rasterize(treatment,inf,background=0,field=1,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
}
# excuation example
# treatRd_Pixel= random_pixel(inf,budget,cost_per_meter_sq,pixelArea)



########### Method 2 -- Select infested pixels for treatment based on total number 
        ## of uninfested host within the given width of buffer area of each pixel ##

Thost <- function(inf,width=1000,ply,host){ # this function will be called by the next function
  
  inf2=inf
  
  
  k<- which(inf2[]==0)
  inf2[k]=NA
  inf2[inf2==0]=NA
  
  raPly=rasterToPolygons(inf2,na.rm=T,dissolve=F)
  
  ply_buf=buffer(raPly,width,dissolve=F)
  ply2=gUnionCascaded(ply)
  ply2=gBuffer(ply2,width=0)
  
  library(parallel)
  library(foreach)
  library(doParallel)
  library(rgeos)
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  ply_buf2=foreach (i=1:dim(pts)[1],.combine=rbind,.errorhandling = "pass",  .export=c('gDifference'), .packages=c('raster',"rgdal","rgeos")) %dopar% {
    #gDifference(plybuf[i,],ply[i,])
    buff_only=gBuffer(ply_buf[i,],width=0)
    buff_only2=gDifference(buff_only, ply2)
    extract(host,buff_only2,na.rm=T,fun=sum)
  }
  
  stopCluster(cl)
  
  hst=as.data.frame(ply_buf2)
  
  raPly$Thst=hst[1:dim(raPly)[1],1]
  raPly2=raPly[order(raPly$Thst,decreasing = T),]
  
  return(raPly2)
}

treat_pixel=function(inf,width=1000,ply,host, budget, buffer, cost_per_meter_sq,pixelArea){
  
  area=budget/cost_per_meter_sq
  
  
  raPly2=Thost(inf,width,ply,host)
  ply_bf=buffer(raPly2,width=buffer,dissolve=F)
  
  n=floor(area/pixelArea)
  ply_bf$Cumu_Area=0
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply_bf[ply_bf$Cumu_Area> area,]
  
  if (df>0){
    crds=gCentroid(nontr[1,])
    crds_bf=buffer(crds,width=sqrt(df/pi))
    treatment=gUnion(treatment,crds_bf)
  }
  
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  return(treatmentLs)
}
## excuation example ##
## treat_Pixel= treat_pixel(inf,width=1000,ply,host,budget,cost_per_meter_sq,pixelArea)



########## Method 3 -- Select highest infested pixel  ###########

Hinfest_pixel = function(inf,buffer,budget,cost_per_meter_sq,pixelArea){
  
  area=budget/cost_per_meter_sq
  
  inf2=inf
  k<- which(inf2[]==0)
  inf2[k]=NA
  inf2[inf2==0]=NA
  
  ra_ply=rasterToPolygons(inf2,na.rm = T,dissolve = F)
  ra_ply$infest_level=extract(inf2,ra_ply,fun=mean)
  
  ra_ply2=ra_ply[order(ra_ply$infest_level,decreasing = T),]
  ply_bf=buffer(ra_ply2,width=buffer,dissolve=F)
  
  
  n=floor(area/pixelArea)
  ply_bf$Cumu_Area=0
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply_bf[ply_bf$Cumu_Area> area,]
  
  if (df>0){
    crds=gCentroid(nontr[1,])
    crds_bf=buffer(crds,width=sqrt(df/pi))
    treatment=gUnion(treatment,crds_bf)
  }
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  return(treatmentLs)
}
