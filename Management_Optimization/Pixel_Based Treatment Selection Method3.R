
## Define Function Input Variable ##
# inf: the infested file in raster format, 
# buffer: the buffer size around each treatment pixel
# width: the buffer size within which the total number of uninfested host would be derived
# ply: the polygon converted from the infested raster, the main use of ply is to exclude the infested host
# pixelArea: the area of each pixel
# distance_classes: only used for Method 5 to indicate how many classes your want to cluster the infested pixel based on 
  #  their distance to the wavefront.
# Wcoef: mean weather coefficient during the simulation period
# Nstep: possible number of reproduction times 

########### Method 1 -- Randomly Select Infested Pixel for Treatment ################
## Randomly select infested pixels based on budget and cost per unit ##

random_pixel=function(inf,buffer,budget,cost_per_meter_sq,pixelArea){
  
  area=budget/cost_per_meter_sq
  
  inf2=inf
  k<- which(inf2[]==0)
  inf2[k]=NA
  
  inf2[inf2==0]=NA
  n=floor(area/pixelArea)+2
  
  rd=sampleRandom(inf2,size=n,na.rm=T,asRaster=T)
  rd_ply=rasterToPolygons(rd,dissolve=F)
  ply_bf=buffer(rd_ply,width=buffer,dissolve=F)
  
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

# idealy, the width should be choosen based on the dispersal ability of targeted species
Thost <- function(ra_ply,width,ply,host){
  
  
  ply_bf=buffer(ra_ply,width,dissolve=F)
  ply2=gUnionCascaded(ply)
  ply2=gBuffer(ply2,width=0)
  
  library(parallel)
  library(foreach)
  library(doParallel)
  library(rgeos)
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  hostN=foreach (i=1:dim(ply_bf)[1],.combine=rbind,.errorhandling = "pass",  .export=c('gDifference'), .packages=c('raster',"rgdal","rgeos")) %dopar% {
    #gDifference(plybuf[i,],ply[i,])
    buff_only=gBuffer(ply_bf[i,],width=0)
    buff_only2=gDifference(buff_only, ply2)
    va=try(extract(host,buff_only2,na.rm=T,fun=sum))
    if (class(va)=="try-error"){va=0}
    else {va=extract(host,buff_only2,na.rm=T,fun=sum)}
    return(va) 
  }
  
  stopCluster(cl)
  
  hst=as.data.frame(unlist(hostN))
  
  ra_ply$Thst=hst[,1]
  return(ra_ply)
}

threat_pixel=function(inf,width=1000,ply,host, budget, buffer, cost_per_meter_sq,pixelArea){
  
  area=budget/cost_per_meter_sq
  
  ra_ply=rasterToPolygons(inf,fun=function(x){x>0},dissolve = F)
  
  ra_ply=Thost(ra_ply,width,ply,host)
  ra_ply2=ra_ply[order(ra_ply$Thst,decreasing = T),]

  ply_bf=buffer(ra_ply2,width=buffer,dissolve=F)
  n=floor(area/pixelArea)
  ply_bf$Cumu_Area=0
  
  n=floor(area/pixelArea)+2
  
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
## treatTrt_Pixel= threat_pixel(inf,width=1000,ply,host,budget,cost_per_meter_sq,pixelArea)


########## Method 3 -- Select highest infested pixel  ###########

Hinfest_pixel = function(inf,buffer,budget,cost_per_meter_sq,pixelArea){
  
  area=budget/cost_per_meter_sq
  
  ra_ply=rasterToPolygons(inf,fun=function(x){x>0},dissolve = F)
  
  # derive buffer size 
  library(parallel)
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  infest_level=foreach (i=1:length(ra_ply),.errorhandling = "pass", .combine=rbind, .export=c('extract'), .packages=c('raster',"rgdal")) %dopar% {
    extract(inf,ra_ply[i,],fun=mean)
  }
  
  stopCluster(cl)
  
  ra_ply$infest_level=as.data.frame(infest_level)[,1]
  
  ra_ply2=ra_ply[order(ra_ply$infest_level,decreasing = T),]
  ply_bf=buffer(ra_ply2,width=buffer,dissolve=F)
  
  
  n=floor(area/pixelArea)+2
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
## treatTrt_Pixel= Hinfest_pixel(inf,buffer,budget,cost_per_meter_sq,pixelArea)


########## Method 4 -- Select infested pixel at wave front ###########

# get wavefront
wvfrt = function(inf){
  ## Derive wave front 
  
  library(concaveman)
  
  pts=rasterToPoints(inf,fun=function(x){x>0},spatial = T)
  
  crds=pts@coords
  concave=concaveman(crds,concavity = 1)
  p=Polygon(concave)
  ps=Polygons(list(p),1)
  sps=SpatialPolygons(list(ps))
  crs(sps)=crs(pts)
  sps=SpatialPolygonsDataFrame(sps,data=data.frame(ID=1))
  
  # lns is used as the wave front
  lns=as(sps,"SpatialLinesDataFrame")
  return(lns)
}

wvfrt_pixel = function(inf,buffer,budget, cost_per_meter_sq, pixelArea){
  
  # get wavefront, lns is spatial Line
  lns=wvfrt(inf)
  
  ra_ply=rasterToPolygons(inf,fun=function(x){x>0},dissolve = F)
  distance=gDistance(lns, ra_ply, byid=T)
  
  ra_ply$dis_frtwv=distance[1:length(ra_ply)]
  ra_ply2=ra_ply[order(ra_ply$dis_frtwv,decreasing = F),]
  
  ply_bf= buffer(ra_ply2,width=buffer,dissolve=F)
  
  area=budget/cost_per_meter_sq
  n=floor(area/pixelArea)
  ply_bf$Cumu_Area=0
  
  n=floor(area/pixelArea)+2
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
  
  treatmentRa=rasterize(treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  return(treatmentLs)
  
}

######### Method 5  -- wavefront method weighted by hazard ########

wvfrt = function(inf){
  ## Derive wave front 
  
  library(concaveman)
  
  pts=rasterToPoints(inf,fun=function(x){x>0},spatial = T)
  
  crds=pts@coords
  concave=concaveman(crds,concavity = 1)
  p=Polygon(concave)
  ps=Polygons(list(p),1)
  sps=SpatialPolygons(list(ps))
  crs(sps)=crs(pts)
  sps=SpatialPolygonsDataFrame(sps,data=data.frame(ID=1))
  
  # lns is used as the wave front
  lns=as(sps,"SpatialLinesDataFrame")
  return(lns)
}

wvfrtHzd = function(inf,width=1000, distance_classes=6, host,buffer,budget,  cost_per_meter_sq, pixelArea){
  
  # get wavefront, lns is spatial Line
  lns=wvfrt(inf)
  
  ra_ply=rasterToPolygons(inf,fun=function(x){x>0},dissolve = F)
  distance=gDistance(lns, ra_ply, byid=T)
  
  ra_ply$dis_frtwv=distance[1:length(ra_ply)]
  ra_ply2=ra_ply[order(ra_ply$dis_frtwv,decreasing = F),]
  ra_ply2=Thost(ra_ply2,width,ply,host)
  
  # based on the distance to wave front, classify all infested pixels/polygons into 
  # 6 classes using kmeans
  library(dplyr)
  centers <- kmeans(ra_ply2$dis_frtwv, centers = distance_classes)$centers
  # order the centers
  centers <- sort(centers)
  group <- kmeans(ra_ply2$dis_frtwv, centers = centers)$cluster
  ra_ply2$group=group
  
  
  # resign group name, pixels/polygons on wavefront have small group number
  ply_bf=buffer(ra_ply2,width=buffer,dissolve=F)
  
  
  n=floor(area/pixelArea)+2
  ply_bf$Cumu_Area=0
  
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  
  ply_select=ply_bf[ply_bf$Cumu_Area< area & ply_bf$Cumu_Area!=0,]
  ply_Nonselect=ply_bf[ply_bf$Cumu_Area> area,]
  
  #get the group of the first pixel/polyon in the non selected pixel/polyon
  group_firstNsel=ply_Nonselect$group[1]
  
  
  ply_groups1=ply_bf[ply_bf$group<=group_firstNsel-1,]
  ply_groups2=ply_bf[ply_bf$group==group_firstNsel,]
  ply_groups2b= ply_groups2[order(ply_groups2$Thst,decreasing = T),]
  
  ply_groups=rbind(ply_groups1,ply_groups2b)
  
  for (i in 1:length(ply_groups)){
    trt=gUnionCascaded(ply_groups[1:i,])
    ply_groups$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_groups[ply_groups$Cumu_Area<=area,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply_groups[ply_groups$Cumu_Area>area,]
  
  if (df>0){
    crds=gCentroid(nontr[1,])
    crds_bf=buffer(crds,width=sqrt(df/pi))
    treatment=gUnion(treatment,crds_bf)
  }
  
  treatmentRa=rasterize(treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa)) 
 
  return(treatmentLs)
    
}



######## Method 6 --  based on number of uninfested host within different  ####
        ## size of buffer for each pixel, the buffer size is determined based on 
        ## dispersal kernel and climatic factors


## Calculate buffer size round each pixel
range_m6 = function(inf, dispersal_rate,reproductive_rate,Nstep,Wcoef){
  
  
  a=dispersal_rate/2
  rep=reproductive_rate
  n=Nstep-1
  
  propa=inf*((rep*Wcoef)^n)
  
  dis_buf=function(a,propa){
    #y=(exp(-x/a)/(2*pi*a))
    #pden=propa*y
    x=-a*log(2*pi*a/propa)
    return(x)
  }
  
  buf= dis_buf(a,propa)
  buf[buf<0]=0.01
  
  return(buf)
}

## derive number of uninfested host within the buffer size for each pixel
host_m6= function(inf,dispersal_rate,reproductive_rate,Nstep,Wcoef,ply,host){
  ra_ply=rasterToPolygons(inf,fun=function(x){x>0},dissolve = F)
  buf=range_m6(inf, dispersal_rate,reproductive_rate,Nstep,Wcoef)
  
  # derive buffer size 
  library(parallel)
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  buf_size=foreach (i=1:length(ra_ply),.errorhandling = "pass", .combine=rbind, .export=c('extract'), .packages=c('raster',"rgdal")) %dopar% {
    extract(buf,ra_ply[i,],fun=mean)
  }
  
  stopCluster(cl)
  
  ra_ply$buf_size=as.data.frame(buf_size)[,1]
  ply_bf=buffer(ra_ply,width=ra_ply$buf_size,dissolve=F)
  
  
  # calculate number of uninfested host within the buffer for each infested pixel/polygon
  ply_bf2= Thost(ply_bf,width=0,ply,host)
  ra_ply$Nhost=ply_bf2$Thst
  
  return(ra_ply)
}

## select treated pixel/poly based on the number of threatened host

treat_m6 = function(inf,dispersal_rate,reproductive_rate,Nstep,Wcoef,ply,host, budget, buffer,cost_per_meter_sq,pixelArea){
  ra_ply= host_m6(inf,dispersal_rate,reproductive_rate,Nstep,Wcoef,ply,host)
  ra_ply2=ra_ply[order(ra_ply$Nhost,decreasing = T),]
  
  area=budget/cost_per_meter_sq
  
  ply_bf=buffer(ra_ply2,width=buffer,dissolve=F)
  n=floor(area/pixelArea)+2
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
