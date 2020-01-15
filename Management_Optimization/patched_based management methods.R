#####################  Functions for patch/polygon based methods ####################

###################### Raster to polygons ##############################

ra2plyMerge <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',pypath="C:\\Program Files\\GDAL\\gdal_polygonize.py") {
  
  require(raster)
  require(rgdal)
  
  x=reclassify(x,c(1,Inf,1))
  
  
  dir=getwd()
  on.exit(setwd(dir))
  setwd(dirname(pypath))
  
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists))
      stop(sprintf('File already exists: %s',
                   toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  
  writeRaster(x, {f <- tempfile(fileext='.tif')})
  rastpath <- normalizePath(f)
  
  system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                  pypath, rastpath, gdalformat, outshape)))
  shp <- readOGR(dirname(outshape), layer = basename(outshape))
  crs(shp)=crs(x)
  
  shp=shp[shp$DN!=0,]
  
  polys=shp
  
  require(rgeos)
  require(maptools)
  
  disMatr=gDistance(polys,polys,byid=T)
  disMatr=as.data.frame(disMatr)
  n=dim(disMatr)[1]
  
  polys$id=1:n
  
  for (i in 1:n){
    polys$id[disMatr[,i]==0]=polys$id[i]
  }
  
  groups=polys$id
  polys2=unionSpatialPolygons(polys,groups)
  n2=length(polys2)
  
  polys3=SpatialPolygonsDataFrame(polys2,data=data.frame(ID=1:n2), match.ID = F)
  return(polys3)
  
  setwd(dir)
}



########################################### Method 1 -- Random selection of infested patches/polys ###############################
random_pch <- function(inf,ply,budget, cost_per_meter_sq, buffer){
  
  area=budget/cost_per_meter_sq
  pixelArea=xres(inf)*yres(inf)
  n=floor(area/pixelArea)+1
  
  ## Random selection
  set.seed(random_seed)
  ply2=ply[sample(1:length(ply),length(ply)),]
  
  ply_bf=buffer(ply2, width=buffer,dissolve=F)
  ply_bf$area=area(ply_bf)
  
  ply_bf$Cumu_Area=ply_bf$area*0
  ply_bf$index=ply_bf$area*0
  
  
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply2[ply_bf$Cumu_Area> area,]
  
  n2=floor(df/pixelArea)+1
  
  if (df>=9){
    
    ply_div= nontr[1,]
    
    div_ra= rasterize(ply_div,inf,background=NA,field=1)
    rd=sampleRandom(div_ra,size=n2,na.rm=T,asRaster=T)
    rd_ply=rasterToPolygons(rd,dissolve=F)
    
    div_bf=buffer(rd_ply,width=buffer,dissolve=F)
    
    if (length(rd_ply)==1){
      crds=gCentroid(rd_ply)
      crds_bf=buffer(crds,width=sqrt(df/pi))
      treatment=gUnion(treatment,crds_bf)
    } else {
      
      div_bf$Cumu_Area=0
      
      for (i in 1:length(div_bf)){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      if (df2>=9){
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=sqrt(df2/pi))
        treatment=gUnion(treatment,crds_bf)
      }
    }
  }
  
  treatmentRa=rasterize(treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
}



###################################### Method 2 -- Select infested patches/polys based on number of uninfested host within ##############
## a given width of constant buffer 
#### Method 2 -- Select infested patches/polys based on number of uninfested host within ####
## a given width of constant buffer.  2 methods are associated with this approach: one is based on the number of nearby uninfested host, 
## and one is based on the (number of nearby uninfested host)/(area of infested patch) 

threat_host <- function(ply,width,all_infested_poly,host, buffer){
  
  plybuf=buffer(ply,width=width,dissolve=F)
  ply2=gUnionCascaded(all_infested_poly)
  ply2=gBuffer(ply2,width=0)
  
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  host_number=foreach (i=1:dim(ply)[1],.errorhandling = "pass",  .export=c('gDifference'), .packages=c('raster',"rgdal","rgeos")) %dopar% {
    #gDifference(plybuf[i,],ply[i,])
    buff_only=gBuffer(plybuf[i,],width=0)
    buff_only=gDifference(plybuf[i,], ply2)
    
    va=try(extract(host,buff_only,na.rm=T,fun=sum))
    if (class(va)=="try-error"){va=0}
    else {va=extract(host,buff_only,na.rm=T,fun=sum)}
    return(va) 
    
  }
  
  stopCluster(cl)
  
  TtHost=as.data.frame(unlist(host_number))
  
  ply$threat_host_number = TtHost[,1]
  ply_bf=buffer(ply,width=buffer,dissolve=F)
  ply_bf$areaBF=area(ply_bf)
  ply$BCratio=ply$threat_host_number/ply_bf$areaBF
  
  return(ply)
}
#ply=threat_host(ply,width,all_infested_poly=ply,host, buffer ,n_cores=10)
# need to run the above code as input for ply, the two functions were seperated to improve efficiency

# based on number of uninfested host within 1000m
hzd1Nu_pch <- function(inf,ply=ply,host, budget, cost_per_meter_sq, buffer,width){
  
  area=budget/cost_per_meter_sq
  pixelArea=xres(inf)*yres(inf)
  n=floor(area/pixelArea)+1
  
  ply2 = ply[order(ply$threat_host_number,decreasing = T),]
  
  ply_bf=buffer(ply2,width=buffer,dissolve=F)
  ply_bf$areaBF=area(ply_bf)
  
  ply_bf$Cumu_Area=ply_bf$areaBF*0
  ply_bf$index=ply_bf$areaBF*0
  
  
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply2[ply_bf$Cumu_Area> area,]
  
  if (df>9){
    
    ply_div= nontr[1,]
    div_ra= rasterize(ply_div,inf,background=NA,field=1)
    div_ply=rasterToPolygons(div_ra,na.rm=T)
    
    if (length(div_ply)==1){
      
      df2=df
      
      r=0
      for (i in 1:999){
        
        r=r+sqrt(df2/pi)
        crds=gCentroid(div_ply)
        crds_bf=buffer(crds,width=r)
        treatment=gUnion(treatment,crds_bf)
        df2=area-area(treatment)
        if (df2<=0){break}
      }
      
    } else {
      
      ra_ply=threat_host(div_ply,width=width,all_infested_poly=ply,host,buffer)
      ra_ply2=ra_ply[order(ra_ply$threat_host_number,decreasing = T),]
      
      div_bf=buffer(ra_ply2,width=buffer,dissolve=F)
      
      
      n2=length(div_bf)
      div_bf$Cumu_Area=0
      
      for (i in 1:n2){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      r=0
      for (i in 1:999){
        
        r=r+sqrt(df2/pi)
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=r)
        treatment=gUnion(treatment,crds_bf)
        df2=area-area(treatment)
        if (df2<=0){break}
      }
      
    }
  }
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
  
}

## based on (number of uninfested host within 1000m)/(area of infested patch)
hzd1Rt_pch <- function(inf,ply=ply,host, budget, cost_per_meter_sq, buffer,width){
  
  area=budget/cost_per_meter_sq
  pixelArea=xres(inf)*yres(inf)
  n=floor(area/pixelArea)+1
  
  ply2 = ply[order(ply$BCratio,decreasing = T),]
  
  ply_bf=buffer(ply2,width=buffer,dissolve=F)
  ply_bf$areaBF=area(ply_bf)
  
  ply_bf$Cumu_Area=ply_bf$areaBF*0
  ply_bf$index=ply_bf$areaBF*0
  
  
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply2[ply_bf$Cumu_Area> area,]
  
  if (df>=9){
    
    ply_div= nontr[1,]
    div_ra= rasterize(ply_div,inf,background=NA,field=1)
    div_ply=rasterToPolygons(div_ra,na.rm=T)
    
    if (length(div_ply)==1){
      
      df2=df
      
      r=0
      for (i in 1:999){
        
        r=r+sqrt(df2/pi)
        crds=gCentroid(div_ply)
        crds_bf=buffer(crds,width=r)
        treatment=gUnion(treatment,crds_bf)
        df2=area-area(treatment)
        if (df2<=0){break}
      }
      
    } else {
      
      ra_ply=threat_host(div_ply,width,all_infested_poly=ply,host,buffer)
      ra_ply2=ra_ply[order(ra_ply$BCratio,decreasing = T),]
      
      div_bf=buffer(ra_ply2,width=buffer,dissolve=F)
      
      
      n2=length(div_bf)
      div_bf$Cumu_Area=0
      
      for (i in 1:n2){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      r=0
      for (i in 1:999){
        
        r=r+sqrt(df2/pi)
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=r)
        treatment=gUnion(treatment,crds_bf)
        df2=area-area(treatment)
        if (df2<=0){break}
      }
      
    }
  }
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
  
}


#################################### Method 3 -- Select highest infested patches/polys #################################
Hinfest_pch = function(inf,ply, buffer,budget,cost_per_meter_sq){
  
  pixelArea= xres(inf)*yres(inf)
  
  area=budget/cost_per_meter_sq
  
  # derive buffer size 
  library(parallel)
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  infest_level=foreach (i=1:length(ply),.errorhandling = "pass", .combine=rbind, .export=c('extract'), .packages=c('raster',"rgdal")) %dopar% {
    extract(inf,ply[i,],fun=mean)
  }
  
  stopCluster(cl)
  
  ply$infest_level=as.data.frame(unlist(infest_level))[,1]
  
  ply2=ply[order(ply$infest_level,decreasing = T),]
  ply_bf=buffer(ply2,width=buffer,dissolve=F)
  
  
  n=floor(area/pixelArea)+1
  ply_bf$Cumu_Area=0
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply2[ply_bf$Cumu_Area> area,]
  
  if (df>=9){
    
    ply_div= nontr[1,]
    div_ra= rasterize(ply_div,inf,background=NA,field=1)
    div_ply=rasterToPolygons(div_ra,na.rm=T)
    
    if (length(div_ply)==1){
      crds=gCentroid(div_ply)
      crds_bf=buffer(crds,width=sqrt(df/pi))
      treatment=gUnion(treatment,crds_bf)
    } else {
      
      div_ply$infest_level=extract(inf,div_ply,fun=mean,na.rm=T)[,1]
      div_ply2=div_ply[order(div_ply$infest_level,decreasing = T),]
      
      div_bf=buffer(div_ply2,width=buffer,dissolve=F)
      
      
      n2=length(div_bf)
      div_bf$Cumu_Area=0
      
      for (i in 1:n2){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      if (df2>=9){
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=sqrt(df2/pi))
        treatment=gUnion(treatment,crds_bf)
      }
    }  
  }
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
  
}



#################################################### Method 4 -- Select infested patches/polys at wave front #############################
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

wvfrt_pch = function(inf,ply, buffer,budget, cost_per_meter_sq){
  
  pixelArea=xres(inf)*yres(inf)
  # get wavefront, lns is spatial Line
  lns=wvfrt(inf)
  
  distance=gDistance(lns, ply, byid=T)
  ply$dis_frtwv=distance[1:length(ply)]
  ply2=ply[order(ply$dis_frtwv,decreasing = F),]
  
  ply_bf= buffer(ply2,width=buffer,dissolve=F)
  
  area=budget/cost_per_meter_sq
  n=floor(area/pixelArea) +1
  ply_bf$Cumu_Area=0
  
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply2[ply_bf$Cumu_Area> area,]
  
  if (df>=9){
    
    ply_div= nontr[1,]
    div_ra= rasterize(ply_div,inf,background=NA,field=1)
    div_ply=rasterToPolygons(div_ra,na.rm=T)
    
    if (length(div_ply)==1){
      crds=gCentroid(div_ply)
      crds_bf=buffer(crds,width=sqrt(df/pi))
      treatment=gUnion(treatment,crds_bf)
    } else {
      
      div_ply$dis_frtwv=gDistance(lns, div_ply, byid=T)[,1]
      div_ply2=div_ply[order(div_ply$dis_frtwv,decreasing = F),]
      
      div_bf=buffer(div_ply2,width=buffer,dissolve=F)
      
      
      n2=length(div_bf)
      div_bf$Cumu_Area=0
      
      for (i in 1:n2){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      if (df2>=9){
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=sqrt(df2/pi))
        treatment=gUnion(treatment,crds_bf)
      }
    }
  }
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
  
}



#################################################### Method 5 -- wavefront method weighted by hazard ########################
wvfrtHzd_pch = function(inf,ply, all_infested_poly=ply,width=1000, distance_classes=6, host,buffer,budget,  cost_per_meter_sq){
  
  pixelArea= xres(inf)*yres(inf)
  area=budget/cost_per_meter_sq
  
  # get wavefront, lns is spatial Line
  lns=wvfrt(inf)
  
  distance=gDistance(lns, ply, byid=T)
  
  ply$dis_frtwv=distance[1:length(ply)]
  ply2=ply[order(ply$dis_frtwv,decreasing = F),]
  ply2=threat_host(ply2,width,all_infested_poly=ply, host,  buffer)
  
  # based on the distance to wave front, classify all infested pixels/polygons into 
  # 6 classes using kmeans
  #library(dplyr)
  centers <- kmeans(ply2$dis_frtwv, centers = distance_classes)$centers
  # order the centers
  centers <- sort(centers)
  group <- kmeans(ply2$dis_frtwv, centers = centers)$cluster
  ply2$group=group
  
  
  ply_bf=buffer(ply2,width=buffer,dissolve=F)
  
  
  n=floor(area/pixelArea)+1
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
  ply_groups2b= ply_groups2[order(ply_groups2$threat_host_number,decreasing = T),]
  
  ply_groups=rbind(ply_groups1,ply_groups2b)
  
  for (i in 1:length(ply_groups)){
    trt=gUnionCascaded(ply_groups[1:i,])
    ply_groups$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_groups[ply_groups$Cumu_Area<=area,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply_groups[ply_groups$Cumu_Area>area,]
  
  if (df>=9){
    
    div_ply= buffer(nontr[1,],width=-buffer,dissolve=F)
    div_ra=rasterize(div_ply,inf,field=1,background=NA)
    div_ply=rasterToPolygons(div_ra,na.rm=T)
    
    if (length(div_ply)==1){
      crds=gCentroid(div_ply)
      crds_bf=buffer(crds,width=sqrt(df/pi))
      treatment=gUnion(treatment,crds_bf)
    } else {
      
      div_ply=threat_host(div_ply,width,all_infested_poly=ply,host,buffer)
      div_ply2=div_ply[order(div_ply$threat_host_number,decreasing = T),]
      
      div_bf=buffer(div_ply2,width=buffer,dissolve=F)
      
      
      n2=length(div_bf)
      div_bf$Cumu_Area=0
      
      for (i in 1:n2){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      if (df2>=9){
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=sqrt(df2/pi))
        treatment=gUnion(treatment,crds_bf)
      }
    }  
  }
  
  treatmentRa=rasterize(treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa)) 
  
  return(treatmentLs)
  
}



############################################## Method 6 -- based on number of uninfested host within different  #####################
## size of buffer for each patch/polygons, the buffer size is determined based on 
## dispersal kernel and clamatic factors. 2 methods are associated with this approach: one is based on the number of nearby uninfested host, 
## and one is based on the (number of nearby uninfested host)/(area of infested patch) 

## Calculate buffer size around each patch/polygon
buffSize <- function(ply,inf,Wcoef, dispersal_rate,reproductive_rate,Nstep){
  
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  
  ply_vars=foreach (i=1:dim(ply)[1], .multicombine = T,.combine=rbind, .packages=c('raster',"rgdal")) %dopar% {
    t=extract(Wcoef,ply[i,],fun=mean)
    d=extract(inf,ply[i,],fun=sum)
    list(t,d)
  }
  
  stopCluster(cl)
  
  Wcoef=as.data.frame(matrix(unlist(ply_vars[,1]),ncol=1))
  Tinf=as.data.frame(matrix(unlist(ply_vars[,2]),ncol=1))
  
  
  
  a=dispersal_rate/2
  rep=reproductive_rate
  n=Nstep-1
  
  propa=Tinf*((rep*Wcoef)^n)
  
  dis_buf=function(a,propa){
    #y=(exp(-x/a)/(2*pi*a))
    #pden=propa*y
    x=-a*log(2*pi*a/propa)
    return(x)
  }
  
  buf= dis_buf(a,propa)
  buf[buf<0]=0.0001
  
  n2=dim(buf)[1]
  ply$buffer=buf$V1[1:n2]
  
  return (ply)
}

## derive buffer area and extract number of uninfested host within the buffer ##
threat_host2 <- function(ply,all_infested_poly,inf,Wcoef, dispersal_rate,reproductive_rate,Nstep,host,buffer){
  
  ply=buffSize(ply,inf,Wcoef, dispersal_rate,reproductive_rate,Nstep)
  plybuf=buffer(ply,width=ply$buffer,dissolve=F)
  ply2=gUnionCascaded(all_infested_poly)
  ply2=gBuffer(ply2,width=0)
  
  
  cl <- makeCluster(detectCores()-4)
  registerDoParallel(cl)
  
  host_number=foreach (i=1:dim(ply)[1],.errorhandling = "pass",  .export=c('gDifference'), .packages=c('raster',"rgdal","rgeos")) %dopar% {
    #gDifference(plybuf[i,],ply[i,])
    buff_only=gBuffer(plybuf[i,],byid=T,width=0)
    buff_only=gDifference(buff_only, ply2)
    
    va=try(extract(host,buff_only,na.rm=T,fun=sum))
    if (class(va)=="try-error"){va=0}
    else {va=extract(host,buff_only,na.rm=T,fun=sum)}
    return(va) 
    
  }
  
  stopCluster(cl)
  
  ply$threat_host_number=as.data.frame(unlist(host_number))[,1]
  ply_bf2=gBuffer(ply,byid=T,width=buffer)
  ply$area=area(ply_bf2)
  ply$BCratio=ply$threat_host_number/ply$area
  
  ## Rank based on benefit cost ratio
  ply_RankRatio=ply[order(ply$BCratio,decreasing = T),]
  
  ## Rank based on number of threatened host
  ply_RankNhost=ply[order(ply$threat_host_number,decreasing = T),]
  
  ls=list(ply_RankRatio,ply_RankNhost)
  
  return(ls)
}

## derive treatment locations, run the below 3 lines before run the treat_m6Ply function

#rgbf=threat_host2(ply,all_infested_poly=ply,inf,Wcoef, dispersal_rate,reproductive_rate,Nstep,host,buffer)
#ply_ratio=rgbf[[1]]
#ply_nhost=rgbf[[2]]
# use ply_ratio or ply_nhost as the input for ply
hzd2Rt_pch <- function(ply,all_infested_poly,inf,budget,cost_per_meter_sq,buffer,Wcoef=Tcoef,dispersal_rate,reproductive_rate,Nstep){
  
  area=budget/cost_per_meter_sq
  n=floor(area/pixelArea)+1
  
  
  ply_bf=buffer(ply,width=buffer,dissolve=F)
  ply_bf$areaBF=area(ply_bf)
  ply_bf$Cumu_Area=ply_bf$areaBF*0
  ply_bf$index=ply_bf$areaBF*0
  
  
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply[ply_bf$Cumu_Area> area,]
  
  if (df>=9){
    
    ply_div= nontr[1,]
    div_ra= rasterize(ply_div,inf,background=NA,field=1)
    div_ply=rasterToPolygons(div_ra,na.rm=T)
    
    if (length(div_ply)==1){
      crds=gCentroid(div_ply)
      crds_bf=buffer(crds,width=sqrt(df/pi))
      treatment=gUnion(treatment,crds_bf)
    } else {
      
      # at this step, the return 1 equal return 2, as area for all pixels is the same
      ra_ply=threat_host2(div_ply,all_infested_poly,inf,Wcoef, dispersal_rate,reproductive_rate,Nstep,host,buffer)[[1]]
      div_bf=buffer(ra_ply,width=buffer,dissolve=F)
      
      
      n2=length(div_bf)
      div_bf$Cumu_Area=0
      
      for (i in 1:n2){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      if (df2>=9){
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=sqrt(df2/pi))
        treatment=gUnion(treatment,crds_bf)
      }
    }
    
  }
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
}
hzd2Nu_pch <- function(ply,all_infested_poly,inf,budget,cost_per_meter_sq,buffer,Wcoef=Tcoef,dispersal_rate,reproductive_rate,Nstep){
  
  area=budget/cost_per_meter_sq
  n=floor(area/pixelArea)+1
  
  
  ply_bf=buffer(ply,width=buffer,dissolve=F)
  ply_bf$areaBF=area(ply_bf)
  ply_bf$Cumu_Area=ply_bf$areaBF*0
  ply_bf$index=ply_bf$areaBF*0
  
  
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  df=area-area(treatment)
  nontr=ply[ply_bf$Cumu_Area> area,]
  
  if (df>=9){
    
    ply_div= nontr[1,]
    div_ra= rasterize(ply_div,inf,background=NA,field=1)
    div_ply=rasterToPolygons(div_ra,na.rm=T)
    
    if (length(div_ply)==1){
      crds=gCentroid(div_ply)
      crds_bf=buffer(crds,width=sqrt(df/pi))
      treatment=gUnion(treatment,crds_bf)
    } else {
      
      # at this step, the return 1 equal return 2, as area for all pixels is the same
      ra_ply=threat_host2(div_ply,all_infested_poly,inf,Wcoef, dispersal_rate,reproductive_rate,Nstep,host,buffer)[[2]]
      div_bf=buffer(ra_ply,width=buffer,dissolve=F)
      
      
      n2=length(div_bf)
      div_bf$Cumu_Area=0
      
      for (i in 1:n2){
        trt2=gUnionCascaded(div_bf[1:i,])
        div_bf$Cumu_Area[i]=area(trt2)
      }
      
      treatment1b=div_bf[div_bf$Cumu_Area<= df & div_bf$Cumu_Area!=0,]
      if (length(treatment1b)==0){
        treatment=treatment} else {
          treatment=gUnionCascaded(gUnion(treatment,treatment1b))}
      
      df2=area-area(treatment)
      nontr2=div_bf[div_bf$Cumu_Area> df,]
      
      if (df2>=9){
        crds=gCentroid(nontr2[1,])
        crds_bf=buffer(crds,width=sqrt(df2/pi))
        treatment=gUnion(treatment,crds_bf)
      }
    }
    
  }
  
  treatmentRa=rasterize( treatment,inf,field=1,background=0,getCover=T)
  treatmentLs=list(as.matrix(treatmentRa))
  
  return(treatmentLs)
}
