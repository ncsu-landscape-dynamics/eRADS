## Step 1 -- rank each pixel based on infestation potential  ##
# inf: infestation raster; 
# host: host map;  
# np: number of step, this parameter doesn't really matter, but need to be a postive number
# a : natural dispersal scale
# rep: reproduction rate
# Wcoef: mean weather coefficient raster of the simulated time period (mean for every year or the whole simulation period)
# range_buffer: a number for distance in meter, for a given infested pixel, only consider susceptible hosts within this range,
            # this number should be a reasonable large number impacted by the natural dispersal rate
# ncore: number of cores used for calculation
ip_rank = function(inf, host, np, a, rep, Wcoef, range_buffer,ncore){
  
  # only counts susceptible hosts for the infestation potential 
  host_uninf = host
  host_uninf[inf>0]=0
  
  # covert infestation and host raster to points to get raster values and 
  # to easier calculate distance between infested pixel and susceptible host
  inf_pts = rasterToPoints(inf,fun=function(x){x>0} ,spatial=T)
  names(inf_pts)="inf_level"
  wc_va = extract(Wcoef,inf_pts)
  inf_pts$inf_level = inf_pts$inf_level 
  
  host_pts = rasterToPoints(host_uninf, fun=function(x){x>0}, spatial=T)
  names(host_pts)="host"
  
  n = length(inf_pts)
  
  ## for each infested pixel, calculate the infestation potential for each susceptible host within the range_buffer
  ## and sum up
  cl = makeCluster(ncore)
  registerDoParallel(cl) 
  
  ip_va = foreach(i=1:n, .combine=rbind, .export=c("host_pts","inf_pts"), .packages = c('raster','rgdal','geosphere','rgeos'))%dopar%{
    dist= gDistance(inf_pts[i,], host_pts,byid=T)
    
    dist2 = 1/(dist^2 + a^2)
    dist2[dist>range_buffer] = 0
    
    # so the calculation of infestation potential is: 
    #(dispersal kernel)* (infestation level)* (reproduction rate) * (weather coefficient of the infested pixel) * (number of step) * (host abundence)
    infe_pot = inf_pts$inf_level[i]*rep*wc_va[i]*np*dist2*host_pts$host
   
    
    infe_potS = sum(infe_pot)
    return(infe_potS)
  }
  
  stopCluster(cl)
  
  inf_pts$ip = as.vector(ip_va)
  inf_ply = rasterToPolygons(inf, fun=function(x){x>0},na.rm = T)
  inf_ply$ip = inf_pts$ip
  inf_ply2 = inf_ply[order(inf_ply$ip,decreasing = T),]
  #inf_pts2 = spTransform(inf_pts2, crs(inf))
  return(inf_ply2)
}

# Step 2 -- select pixels for treatment based on rank of infestation potential #
# ip_rank_ply: return of the ip_rank function (polygon dataframe), each polygon in the polygon dataframe is one infested pixel
              # the polygons are already ranked based on the infestation potential (ip) in a way that the first polgyon has the highest ip
# inf : infestation raster
# budget, cost_per_meter_sq : the same variables in your functions
# buffer: treatment buffer (in meter)
ip_treat = function(ip_rank_ply, budget, buffer, cost_per_meter_sq, inf){
  
  area=budget/cost_per_meter_sq
  pixelArea=xres(inf)*yres(inf)
  
  
  ply_bf=buffer(ip_rank_ply,width=buffer,dissolve=F)
  n=floor(area/pixelArea)+2
  ply_bf$Cumu_Area=0
  
  # 
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=area(trt)
  }
  
  # select pixels together with the treatment buffer whoes total area is not larger than the budget allowed
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  # if the selected treatment area is smaller than budget allowed, select part of the next infested pixel
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
