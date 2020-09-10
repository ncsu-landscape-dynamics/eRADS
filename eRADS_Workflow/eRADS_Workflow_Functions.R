#### Step 0 -- Delineate Surveillance Range ####
#### Function 1 -- format points ####
# function to reproject points into simple mercator projection
# and to only keep presences found in the U.S. 

pts_format = function(presence_pts){
  load("us_County48.RData")
  load("usall.RData")
  
  pre = spTransform(presence_pts, crs(usall))   
  tfs=as.data.frame(matrix(gIntersects(pre, usall, byid=T), ncol=1))
  pre=pre[tfs$V1==T,]
  
  return(pre)
}

#### Function 2 -- determine delineate range  ####
# coarse_estimation_of_annual_spread_rate in meter 
# buffer_size: size of buffer within which surveillance is desired around a presence in meter
# road_only: Boolean variable TRUE/FALSE, use true if only want to survey locations near to road
# road_buffer: float, when road_only = T, size of road buffer outside which survey wouldn't be conducted
delineate_range = function(presence_pts, coarse_estimation_of_annual_spread_rate, buffer_size, road_only, road_buffer){
  
  ## get centroid of major invasion populations ##
  ctd=gCentroid(presence_pts)
  ## get distance between all presences and infestation centroid ##
  dis=gDistance(ctd,presence_pts,byid=T)
  ## get maximum distance between centroid and presences
  max_dis = max(dis)
  
  ## calculate minimum surveillance range
  min_surv_size = max_dis + buffer_size
  min_surv = gBuffer(ctd, width = min_surv_size)
  
  ## calculate the suggested surveillance range
  sug_surv_size = max_dis + max(coarse_estimation_of_annual_spread_rate, buffer_size)
  sug_surv = gBuffer(ctd, width = sug_surv_size)
  
  if (road_only==T){
    load("Roads.RData")
    rds = crop(rds, extent(sug_surv))
    rds_bf = buffer(rds, width = road_buffer)
    
    sug_surv = gIntersection(sug_surv, rds_bf)
    min_surv = gIntersection(min_surv, rds_bf)
  }
  
  ls = list(ctd, min_surv, sug_surv)
  names(ls) = c("Invasion_centroid", "Minimum_surveillance_range", "Suggested_surveillance_range")
  return(ls)
  
}


#### Function 3 -- suggest surveillance density for delineation ####

# assume the number of pest found in each presence point is instored in the attribute named number_pest
# if number_pest attribute is missing, assume 1 pest for each presence point
# confi - float (0-1), desired detection confidence
# sensi - float (o-1), sensitivity of detection, be conservative on this estimation 
delineate_density = function(confi,presence_pts,sensi){
  
  # get convex hull of the presence points
  coords=presence_pts@coords
  concave=concaveman(coords)
  
  p=Polygon(concave)
  ps=Polygons(list(p),1)
  sps=SpatialPolygons(list(ps))
  crs(sps)=crs(presence_pts)
  plys=SpatialPolygonsDataFrame(sps,data=as.data.frame("id"))
  
  # estimate observed infestation density 
  if (isFALSE("nmbr_ps"%in%names(presence_pts))){presence_pts$nmbr_ps = 1}
  presence_pts$nmbr_ps = as.numeric(presence_pts$nmbr_ps)
  
  pop_den = sum(presence_pts$nmbr_ps)/(gArea(plys)/1000000)
  #
  surv_den=log(1-confi)/log(1-pop_den*sensi)
  names(surv_den) = "Surveillance_density"
  return(surv_den)
}


####### Step 1 -- estimate dispersal parameters based on first 2 or 3 years' sampling data  ########
# Function 1 input data: 1) host map 
# 2) folder path where all the presence points shapefiles are located
# each shapefile includes all presence points for one year, and the shapefile should be named properly so
# the list.file() function in R list the file for earliest year first and the latest year last

# this function return a rasterstack, and each raster represent the infestation in one year
pts_raster = function(inf_pts_folder, host, inf_raster_folder){
  
  pts_files = list.files(inf_pts_folder, ".shp$")
  
  inf_st = stack()
  for (i in 1:length(pts_files)){
    pts = readOGR(inf_pts_folder, strsplit(pts_files[i],".shp")[[1]][1])
    pts = spTransform(pts, crs(host))
    
    if (isFALSE("nmbr_ps"%in%names(pts))){pts$nmbr_ps = 1}
    pts$nmbr_ps = as.numeric(pts$nmbr_ps)
    
    pts_inf = rasterize(pts,host, background=0, field="nmbr_ps" ,fun="count")
    raster_name1 = paste("InfestationYear", as.character(i),".tif",sep="")
    raster_name2 = paste(inf_raster_folder, raster_name1, sep="")
    writeRaster(pts_inf, raster_name2)
    
    inf_st = stack(inf_st, pts_inf)
  }
  names(inf_st)= unlist(lapply(1:length(pts_files), function(x){paste("year",as.character(x),sep="")}))
  
  #write_name = paste(inf_raster_stack_folder,"infestations.tif",sep="")
  #writeRaster(inf_st, write_name)
  return(inf_st)
}

#### Function 2 -- estimate dispersal parameters for running PoPS model Using Chris' ABC Method####
parameters <- abc_calibration(infected_years_file, 
                              number_of_observations,
                              prior_number_of_observations,
                              prior_means, prior_cov_matrix, 
                              params_to_estimate = c(T, T, T, T, F, F),
                              number_of_generations = 20,
                              generation_size = 5000, 
                              checks = c(1000,5000, 100, 200),
                              infected_file,
                              host_file, 
                              total_plants_file, 
                              temp = FALSE, 
                              temperature_coefficient_file = "", 
                              precip = FALSE, 
                              precipitation_coefficient_file = "", 
                              model_type = "SI", 
                              latency_period = 0,
                              time_step = "month", 
                              season_month_start = 3,
                              season_month_end = 10, 
                              start_date = '2017-03-01',
                              end_date = '2018-12-30', 
                              use_lethal_temperature = FALSE, 
                              temperature_file = "",
                              lethal_temperature = -15.87, 
                              lethal_temperature_month = 1,
                              mortality_on = FALSE, 
                              mortality_rate = 0, 
                              mortality_time_lag = 0, 
                              management = FALSE, 
                              treatment_dates = c("2017-03-01"), 
                              treatments_file = "",
                              treatment_method = "ratio",
                              natural_kernel_type = "cauchy", 
                              anthropogenic_kernel_type = "cauchy",
                              natural_dir = "NONE", 
                              natural_kappa = 0, 
                              anthropogenic_dir = "NONE", 
                              anthropogenic_kappa = 0,
                              pesticide_duration = c(0), 
                              pesticide_efficacy = 1.0,
                              mask = NULL, 
                              success_metric = "number of locations", 
                              output_frequency = "year",
                              movements_file = "", 
                              use_movements = FALSE)


parameters$posterior_means

####### Step 2 -- estimate potential infestation intensity and infestation risk using PoPS ######
# input data: all data needed to run the pops_model

infestation_estimation = function(use_lethal_temperature, lethal_temperature,
                                  lethal_temperature_month, infected, susceptible, total_plants,
                                  mortality_on, mortality_tracker, mortality, treatment_maps,
                                  treatment_years, weather, temperature, weather_coefficient, ew_res,
                                  ns_res, num_rows, num_cols, time_step, reproductive_rate,
                                  mortality_rate = 0, mortality_time_lag = 2L,
                                  season_month_start = 1L, season_month_end = 12L, start_time = 2018,
                                  end_time = 2018, treatment_month = 12L, treatment_method = "ratio",
                                  natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                                  use_anthropogenic_kernel = FALSE, percent_natural_dispersal = 0,
                                  natural_distance_scale = 21, anthropogenic_distance_scale = 0,
                                  natural_dir = "NONE", natural_kappa = 0,
                                  anthropogenic_dir = "NONE", anthropogenic_kappa = 0){
  
  intensity_st = stack()
  risk_st = stack()
  
  for (i in 1:100){
    
    random_seed= i*121+5003
    
    data <- pops_model(random_seed, use_lethal_temperature, lethal_temperature,
                       lethal_temperature_month, infected, susceptible, total_plants,
                       mortality_on, mortality_tracker, mortality, treatment_maps,
                       treatment_years, weather, temperature, weather_coefficient, ew_res,
                       ns_res, num_rows, num_cols, time_step, reproductive_rate,
                       mortality_rate = 0, mortality_time_lag = 2L,
                       season_month_start = 1L, season_month_end = 12L, start_time = 2018,
                       end_time = 2018, treatment_month = 12L, treatment_method = "ratio",
                       natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                       use_anthropogenic_kernel = FALSE, percent_natural_dispersal = 0,
                       natural_distance_scale = 21, anthropogenic_distance_scale = 0,
                       natural_dir = "NONE", natural_kappa = 0,
                       anthropogenic_dir = "NONE", anthropogenic_kappa = 0)
    
    
    inf=data$infected[[1]]
    inf=as.data.frame(inf)
    inf=as(inf,"matrix")
    inf=raster(inf,template=inf)
    
    intensity_st = stack(intensity_st, inf)
    
    inf[inf>=1]=1
    inf[inf<1]=0
    risk_st = stack(risk_st, inf)
    
    
  }
  
  pred_intensity = mean(intensity_st)
  pred_risk = sum(risk_st)
  
  ls = list(pred_intensity, pred_risk)
  names(ls) = c("Estimated_density_map", "Estimated_risk")
  return(ls)
}


####### Step 3 -- Surveillance optimization #######

#### define important functions ####

## Detection probability/confidence function given surveillance density, estimated infestation
## population density, and detection sensitivity ##
p_dtc = function(surv_den,p_den,sensi){
  p_dtc=1-(1-sensi)^(surv_den*p_den)
  names(p_dtc)="Detection_confidence"
  return(p_dtc)
}

## Surveillance density for each pixel function given desired detection confidence, estimated infestation
## population density, and detection sensitivity##
surv_den = function(confi,p_den,sensi){
  #surv_den=log(1-confi, base = exp(1))/log(1-p_den*sens, base = exp(1))
  p_den[p_den==0]=NA
  surv_den=log(1-confi)/(p_den*log(1-sensi))
  names(surv_den)="Surveillance_density"
  return(surv_den)
}

## Select given percentage, e.g. 95%, of locations with highest infestation intensity
percent_location = function(density_map, percentage){
  
  number_locations = sum(getValues(density_map/density_map),na.rm = T)
  number_select_locations = floor(number_locations*percentage)
  
  pts=rasterToPoints(density_map,fun=function(x){x>0},spatial = T)
  names(pts) = 'density'
  
  pts2 = pts[order(pts$density,decreasing = T),]
  pts2_select = pts2[1:number_select_locations,]
  
  surv_locations = rasterize(pts2_select, density_map, field=1, background=0)
  surv_locations = surv_locations*density_map
  density_map = surv_locations
  
  names(density_map) = "Inf_density_for_selected_pixels"
  return(density_map)
}

## Select locations with highest infestation intensity, which totally include given percentage 
## of estimated infestation population
percent_population = function(density_map, percentage){
  
  all_population = sum(getValues(density_map), na.rm=T)
  percent_population = all_population*percentage
  
  pts=rasterToPoints(density_map,fun=function(x){x>0},spatial = T)
  names(pts) = 'density'
  
  pts2 = pts[order(pts$density,decreasing = T),]
  pts2$cumulative_population = sapply(1:length(pts2),FUN=function(i)sum(pts2$density[1:i]))
  pts3 = pts2[pts2$cumulative_population <= percent_population,]
  
  surv_locations = rasterize(pts3, density_map, field=1, background=0)
  surv_locations = surv_locations*density_map
  density_map = surv_locations
  names(density_map) = "Inf_density_for_selected_pixels"
  
  return(density_map)
}

## optimal percent location and population and detection confidencedence ##
confidence_all = function(confidence,density_map,sensi,budget,cost_per_survey){
  
  number_survey = floor(budget/cost_per_survey)
  
  survDen = surv_den(confidence,density_map,sensi)
  surv_den_values = getValues(survDen)
  surv_den_values = as.data.frame(surv_den_values)
  surv_den_values$ID = 1:dim(surv_den_values)[1]
  surv_den_values$Pden = getValues(density_map)
  
  surv_den_values2 = surv_den_values[order(surv_den_values$surv_den_values,decreasing = F,na.last = NA),]
  surv_den_values2 = surv_den_values2[surv_den_values2$surv_den_values>0,]
  surv_den_values2$cumu_number = sapply(1:dim(surv_den_values2)[1],function(i) sum(surv_den_values2$surv_den_values[1:i]))
  surv_den_values3 = surv_den_values2[surv_den_values2$cumu_number<= number_survey & surv_den_values2$cumu_number>0,]
  
  number_select = dim(surv_den_values3)[1]
  number_all = sum(getValues(density_map/density_map),na.rm = T)
  
  all_conf_loca = (number_select/number_all)*confidence
  all_conf_popu = (sum(surv_den_values3$Pden)/sum(getValues(density_map),na.rm = T))*confidence
  percent_loca = (number_select/number_all)
  percent_popu = (sum(surv_den_values3$Pden)/sum(getValues(density_map),na.rm = T))
  confall = list(all_conf_loca, percent_loca, all_conf_popu, percent_popu)
  names(confall) = c("Location_mean_confidence","Percent_location","Population_mean_confidence","Percent_population")
  
  return(confall)
}
opticonfidence = function(density_map,sensi,budget,cost_per_survey,confidenceMin=0.50,confidenceMax=0.99, number_core){
  require(parallel)
  require(foreach)
  require(doParallel)
  
  cl = makeCluster(number_core)
  registerDoParallel(cl)
  
  confi = seq(confiMin,confiMax,0.01)
  n = length(confi)
  
  confAll = foreach(i=1:n, .combine = rbind,.export = c("confidence_all","sensi","surv_den"),.packages = c("raster","rgdal","rgeos")) %dopar% {
    confall = confidence_all(confi[i],density_map,sensi,budget,cost_per_survey)
    result = as.data.frame(matrix(c(confall[[1]],confall[[2]],confall[[3]],confall[[4]]),nrow=1,ncol=4))
    colnames(result) = c("All_Confidence_Location","Percent_Location","All_Confidence_Population","Percent_Population")
    
    return(result)
  }
  
  stopCluster(cl) 
  results = as.data.frame(confAll) 
  results$confi = seq(confiMin,confiMax,0.01)
  optiConfi_location = results[which.max(results$All_Confidence_Location),]$confi
  optiConfi_population = results[which.max(results$All_Confidence_Population),]$confi
  
  optiConfi_location_percent = results[which.max(results$All_Confidence_Location),]$Percent_Location
  optiConfi_population_percent = results[which.max(results$All_Confidence_Population),]$Percent_Population
  
  location_ls = list(optiConfi_location, optiConfi_location_percent)
  names(location_ls) = c("Optimal_location_confidence", "Optimal_location_percentage")
  population_ls = list(optiConfi_population, optiConfi_population_percent)
  names(population_ls) = c("Optimal_population_confidence", "Optimal_population_percentage")
  
  ls = list(location_ls,population_ls)
  names(ls) = c("Based_on_location","Based_on_population")
  return(ls)
}


## Estimate infestation potential of each pixel if get infested ##
infestation_potential = function(density_map, host_map,  dispersal_scale, Wcoef=NA, range_buffer=2000, percentage_location, number_core){
  
  # only counts susceptible hosts for the infestation potential 
  host_uninf = host_map
  host_uninf[density_map>0]=0
  
  # covert infestation and host raster to points to get raster values and 
  # to easier calculate distance between infested pixel and susceptible host
  inf_pts = rasterToPoints(density_map,fun=function(x){x>0} ,spatial=T)
  names(inf_pts)="inf_level"
  inf_pts$inf_level = inf_pts$inf_level 
  # extract weather_coefficient for each infested pixel/points     
  #if(is.na(Wcoef)){Wcoef = host_map*0+1}
  wc_va = extract(Wcoef,inf_pts)
  
  host_pts = rasterToPoints(host_uninf, fun=function(x){x>0}, spatial=T)
  names(host_pts)="host"
  
  n = length(inf_pts)
  
  ## for each infested pixel, calculate the infestation potential for each susceptible host within the range_buffer
  ## and sum up
  cl = makeCluster(number_core)
  registerDoParallel(cl) 
  
  ip_va = foreach(i=1:n, .combine=rbind, .packages = c('raster','rgdal','geosphere','rgeos'))%dopar%{
    dist= gDistance(inf_pts[i,], host_pts,byid=T)
    
    dist2 = 1/(dist^2 + dispersal_scale^2)
    dist2[dist>range_buffer] = 0
    
    # so the calculation of infestation potential is: 
    infe_pot = inf_pts$inf_level[i]*wc_va[i]*dist2*host_pts$host
    
    infe_potS = sum(infe_pot)
    return(infe_potS)
  }
  
  stopCluster(cl)
  
  inf_pts$ip = as.vector(ip_va)
  inf_ply = rasterToPolygons(density_map, fun=function(x){x>0},na.rm = T)
  inf_ply$ip = inf_pts$ip
  inf_ply2 = inf_ply[order(inf_ply$ip,decreasing = T),]
  
  select_location = floor(length(inf_ply2)*percentage_location)
  inf_ply2 = inf_ply2[1:select_location, ]
  ip_map = rasterize(inf_ply2, density_map, field="ip", background=NA)
  
  #inf_pts2 = spTransform(inf_pts2, crs(inf))
  names(ip_map)="Infestation_potential_map"
  return(ip_map)
}


#### Method 1 -- Risk Based Strategy 1 -- Regular Surveillance ####
# budget: float - budget for surveillance
# cost_per_survey: float 
# risk_map: raster file
# risk_threshold: float -- regions with risk lower than this threshold won't be surveyed

surv_strategy1 = function(budget, cost_per_survey, risk_map, risk_threshold, road_only, road_buffer){
  
  risk_map[risk_map<risk_threshold] = 0 # set pixel values with risk lower than threshold to 0
  number_survey = floor(budget / cost_per_survey)
  
  
  # limit surveillance region around road if road_only = T
  if (road_only == T){
    
    load("Roads.RData")
    rds = crop(rds, extent(risk_map))
    rds_bf = gBuffer(rds, byid=T , width = road_buffer)
    
    risk_map = crop(risk_map, rds_bf)
    risk_map = mask(risk_map, rds_bf)
    
  }
  
  
  number_pixel = sum(getValues(risk_map/risk_map), na.rm=T)
  pixel_area = xres(risk_map)*yres(risk_map)/1000000
  surv_area = number_pixel * pixel_area
  surv_density = number_survey / surv_area
  
  ls = list(surv_density, risk_map)
  names(ls)=c("Surveillance_density", "Risk_of_surveillance_pixels")
  return(ls)
}


#### Method 2 -- Risk Based Strategy 2 -- Risk-Weighted Surveillance ####
surv_strategy2 = function(budget, cost_per_survey, risk_map, risk_threshold ,road_only, road_buffer, risk_level=T, number_risk_level){
  
  risk_map[risk_map<risk_threshold] = 0 # set pixel values with risk lower than threshold to 0
  number_survey = floor(budget / cost_per_survey)
  
  
  # limit surveillance region around road if road_only = T
  if (road_only == T){
    
    load("Roads.RData")
    rds = crop(rds, extent(risk_map))
    rds_bf = gBuffer(rds, byid=T , width = road_buffer)
    
    risk_map = crop(risk_map, rds_bf)
    risk_map = mask(risk_map, rds_bf)
    
  }
  
  # approach 1 -- assign surveillance density exactly based on risk #
  surv_Density_map1 = risk_map * number_survey/ sum(getValues(risk_map),na.rm=T)
  
  # approach 2 -- classify risk into different level #
  if(risk_level==T){
    length_interval = 100/number_risk_level 
    breaks= c(seq(risk_threshold, 100, by = length_interval),100)
    
    n = length(breaks)-1
    for (i in 1:n){
      
      risk_map = reclassify(risk_map, c(breaks[i],breaks[i+1], i), include.lowest= T, right=F)  
      
    }
    
    surv_Density_map2 = risk_map * number_survey/ sum(getValues(risk_map),na.rm=T)
    
  }
  
  
  if (risk_level==T){
    result = surv_Density_map2
  } else ( result = surv_Density_map1)
  
  names(result) = "Surveillance_density"
  return(result)
}


#### Method 3 -- Infestation Density Based Strategy 1 -- Regular Surveillance ####
surv_strategy3 = function( budget, cost_per_survey, density_map, density_threshold, road_only, road_buffer){
  
  density_map[density_map<density_threshold] = 0 # set pixel values with risk lower than threshold to 0
  number_survey = floor(budget / cost_per_survey)
  
  
  # limit surveillance region around road if road_only = T
  if (road_only == T){
    
    load("Roads.RData")
    rds = crop(rds, extent(risk_map))
    rds_bf = gBuffer(rds, byid=T , width = road_buffer)
    
    density_map = crop(density_map, rds_bf)
    density_map = mask(density_map, rds_bf)
    
  }
  
  
  number_pixel = sum(getValues(density_map/density_map), na.rm=T)
  pixel_area = xres(density_map)*yres(density_map)/1000000
  surv_area = number_pixel * pixel_area
  surv_density = number_survey / surv_area
  
  ls = list(surv_density, density_map)
  names(ls) = c("Surveillance_denisty", "Infestation_density_for_surveillance_pixels")
  return(ls)
}

#### Method 4 -- Infestation Density Based Strategy 2 -- Same Detection Confidence Method####
surv_strategy4 = function(use_budget=T, detection_confidence=0.95 , budget, cost_per_survey, density_map, density_threshold, road_only, road_buffer, sensi){
  
  density_map[density_map<density_threshold] = 0 # set pixel values with risk lower than threshold to 0
  density_map[density_map==0]=NA
  
  # limit surveillance region around road if road_only = T
  if (road_only == T){
    
    load("Roads.RData")
    rds = crop(rds, extent(risk_map))
    rds_bf = gBuffer(rds, byid=T , width = road_buffer)
    
    density_map = crop(density_map, rds_bf)
    density_map = mask(density_map, rds_bf)
    
  }
  
  
  if (use_budget==F){
    survDensity = surv_den(detection_confidence, density_map, sensi)
    n_survey = sum(getValues(survDensity), na.rm = T)
    
    ls = list(n_survey, survDensity)
    names(ls) = c("Total_number_surveillance", "Surveillance_density_map") 
  } else {
    
    number_survey = floor(budget / cost_per_survey)
    
    survDensity = surv_den(0.5, density_map, sensi)
    survDensity2 = survDensity * number_survey / sum(getValues(survDensity), na.rm = T)
    
    confi_map = p_dtc(survDensity2, density_map, sensi)
    confi_value = getValues(confi_map)
    confi_value[is.na(confi_value)] = 0
    confi_value = max(confi_value)
    
    ls = list(confi_value, survDensity2)
    names(ls) = c("Confidence_of_surveillance_pixels", "Surveillance_density")
  }
  
  return(ls)
}

#### Method 5 -- Infestation Density Based Strategy 3 -- Prioritize percentage locations with higher infestation density ####
surv_strategy5 = function(use_budget=T, detection_confidence=0.95 , budget, cost_per_survey, density_map, density_threshold, road_only, road_buffer, percentage, location_percentage = T, sensi){
  
  all_number_pixel = sum(getValues(density_map/density_map), na.rm=T)  
  density_map[density_map<density_threshold] = 0 # set pixel values with risk lower than threshold to 0
  number_survey = floor(budget / cost_per_survey)
  
  # selection locations for surveillance
  if (location_percentage==T){
    density_map = percent_location(density_map, percentage)
  } else{
    density_map = percent_population(density_map, percentage )
  }
  
  
  # limit surveillance region around road if road_only = T
  if (road_only == T){
    
    load("Roads.RData")
    rds = crop(rds, extent(risk_map))
    rds_bf = gBuffer(rds, byid=T , width = road_buffer)
    
    density_map = crop(density_map, rds_bf)
    density_map = mask(density_map, rds_bf)
    
  }
  
  if (use_budget==F){
    survDensity = surv_den(detection_confidence, density_map, sensi)
    n_survey = sum(getValues(survDensity), na.rm = T)
    
    ls = list(n_survey, survDensity)
    names(ls) = c("Total_number_surveillance", "Surveillance_density_map") 
  }  else {
    # conduct regular survey for the selected locations
    number_survey_pixel = sum(getValues(density_map/density_map),na.rm=T)
    select_location_percentage = number_survey_pixel/all_number_pixel
    survey_area = number_survey_pixel*xres(density_map)*yres(density_map)/1000000
    regular_surv_den = number_survey/survey_area
    
    
    # conduct safe detection confidence for the selected locations
    same_confi_surv_results = surv_strategy4(use_budget=T, detection_confidence, budget, cost_per_survey, density_map, density_threshold=0, road_only=F, road_buffer=0, sensi)
    same_detection_confi = same_confi_surv_results[[1]]
    same_confi_surv_den = same_confi_surv_results[[2]]
    
    if (location_percentage==T){
      ls = list(density_map, regular_surv_den, same_detection_confi, same_confi_surv_den)
      names(ls) = c("Infestation_density_of_surveillance_pixels", "Regular_surveillance_density", "Sameconfi_detection_confi", "Surveillance_denisty_sameConfi")
    } else {ls = list(select_location_percentage, density_map, regular_surv_den, same_detection_confi, same_confi_surv_den)
    names(ls) = c("Percentage_location_selected", "Infestation_density_of_surveillance_pixels", "Regular_surveillance_density", "Sameconfi_detection_confi", "Surveillance_denisty_sameConfi")
    }
  }
  
  return(ls)
}


#### Method 6 -- Infestation Density Based Strategy 4 -- Survey location with estimated infestation density higher than a given threshold ####
surv_strategy6 = function(use_budget=T, detection_confidence=0.95 , budget, cost_per_survey, density_map, density_threshold, road_only, road_buffer, sensi){
  
  number_survey = floor(budget / cost_per_survey)
  all_number_pixel = sum(getValues(density_map/density_map),na.rm = T)
  density_map[density_map<density_threshold]=0
  
  # limit surveillance region around road if road_only = T
  if (road_only == T){
    
    load("Roads.RData")
    rds = crop(rds, extent(risk_map))
    rds_bf = gBuffer(rds, byid=T , width = road_buffer)
    
    density_map = crop(density_map, rds_bf)
    density_map = mask(density_map, rds_bf)
    
  }
  
  if (use_budget==F){
    survDensity = surv_den(detection_confidence, density_map, sensi)
    n_survey = sum(getValues(survDensity), na.rm = T)
    
    ls = list(n_survey, survDensity)
    names(ls) = c("Total_number_surveillance", "Surveillance_density_map") 
  } else {
    # conduct regular survey for the selected locations
    number_survey_pixel = sum(getValues(density_map/density_map),na.rm=T)
    select_location_percentage = number_survey_pixel/all_number_pixel
    survey_area = number_survey_pixel*xres(density_map)*yres(density_map)/1000000
    regular_surv_den = number_survey/survey_area
    
    
    # conduct safe detection confidence for the selected locations
    same_confi_surv_results = surv_strategy4(use_budget=T, detection_confidence, budget, cost_per_survey, density_map, density_threshold=0, road_only=F, road_buffer=0, sensi)
    same_detection_confi = same_confi_surv_results[[1]]
    same_confi_surv_den = same_confi_surv_results[[2]]
    
   
    # conduct regular survey for the selected locations
    number_survey_pixel = sum(getValues(density_map/density_map),na.rm=T)
    select_location_percentage = number_survey_pixel/all_number_pixel
    survey_area = number_survey_pixel*xres(density_map)*yres(density_map)/1000000
    regular_surv_den = number_survey/survey_area
    
    
    # conduct same detection confidence for the selected locations
    same_confi_surv_results = surv_strategy4(use_budget=T, detection_confidence, budget, cost_per_survey, density_map, density_threshold=0, road_only=F, road_buffer=0, sensi)
    same_detection_confi = same_confi_surv_results[[1]]
    same_confi_surv_den = same_confi_surv_results[[2]]
    
    ls = list(select_location_percentage, density_map, regular_surv_den, same_detection_confi, same_confi_surv_den)
    names(ls) = c("Percentage_location_selected", "Infestation_density_of_surveillance_pixels", "Regular_surveillance_density", "Sameconfi_detection_confi", "Surveillance_denisty_sameConfi")
    
  }
  
  return(ls)
}


#### Method 7 -- Risk and Infestation Density Based Strategy -- Survey location with higher infestation potential (i.e. potential of causing more infestation) ####
surv_strategy7 = function(use_budget=T, detection_confidence=0.95 , budget, cost_per_survey, risk_map, density_map, density_threshold, percentage_location, road_only, road_buffer, sensi, dispersal_scale,Wcoef=NA, range_buffer=2000){
  
  library(doParallel)
  library(parallel)
  
  density_map[density_map<density_threshold] = 0
  risk_values = getValues(risk_map)
  risk_values[is.na(risk_values)]=0
  density_map = density_map*risk_map
  number_survey = floor(budget / cost_per_survey)
  
  if(is.na(Wcoef)){
    Wcoef = density_map*0 + 1
  } 
  
  ip_map = infestation_potential(density_map, host_map, 100, Wcoef, range_buffer = 5000,  percentage_location ,number_core=5)
  density_map = ip_map
  
  # limit surveillance region around road if road_only = T
  if (road_only == T){
    
    load("Roads.RData")
    rds = crop(rds, extent(risk_map))
    rds_bf = gBuffer(rds, byid=T , width = road_buffer)
    
    density_map = crop(density_map, rds_bf)
    density_map = mask(density_map, rds_bf)
    
  }
  
  if (use_budget==F){
    survDensity = surv_den(detection_confidence, density_map, sensi)
    n_survey = sum(getValues(survDensity), na.rm = T)
    
    ls = list(n_survey, survDensity)
    names(ls) = c("Total_number_surveillance", "Surveillance_density_map") 
  } else {
    # conduct regular survey for the selected locations
    number_survey_pixel = sum(getValues(density_map/density_map),na.rm=T)
    
    survey_area = number_survey_pixel*xres(density_map)*yres(density_map)/1000000
    regular_surv_den = number_survey/survey_area
    regular_surv_confi = p_dtc(regular_surv_den, density_map, sensi)
    
    
    # conduct safe detection confidence for the selected locations
    same_confi_surv_results = surv_strategy4(use_budget=T, detection_confidence, budget, cost_per_survey, density_map, density_threshold=0, road_only=F, road_buffer=0, sensi)
    same_detection_confi = same_confi_surv_results[[1]]
    same_confi_surv_den = same_confi_surv_results[[2]]
    
    ls = list(density_map, regular_surv_den, same_detection_confi, same_confi_surv_den)
    names(ls) = c("Infestation_density_of_surveillance_pixels", "Regular_surveillance_density", "Sameconfi_detection_confi", "Surveillance_denisty_sameConfi")
    
  }
   
  return(ls)
}


####### Step 4 -- Evaluate efficacy  of each method ######

#### Evaluate on the estimated potential infestation intensity map (can be conducted before getting the sample data) ####
# use this function to get the detection confidence curve of all strategies with different bugget level #
strategy_evaluation = function(density_map, density_threshold=0, surv_den, sensi){
  density_map[density_map< density_threshold]=0
  
  confidence_map = p_dtc(surv_den, density_map, sensi)
  confidence_map[density_map==0]=NA
  confidence_map[density_map>0 & is.na(confidence_map)]=0
  mean_confidence =mean(getValues(confidence_map), na.rm=T)
  names(mean_confidence)="Mean_confidence"
  return(mean_confidence)
}

#### Evaluate on the kernel density map derived from actual sample data ####
# function 1 -- estimate kernel density map from presence points #

kernel_density = function(presence_pts){
  
  library(spatialEco)
  library(sparr)
  
  if (isFALSE("nmbr_ps"%in%names(presence_pts))){presence_pts$nmbr_ps = 1}
  presence_pts$nmbr_ps = as.numeric(presence_pts$nmbr_ps)
  
  ow = owin(xrange = range(presence_pts@coords[,1]), yrange= range(presence_pts@coords[,2]), xy=presence_pts@coords)
  pp = ppp(presence_pts@coords[,1], presence_pts@coords[,2], window = ow)
  bw = LIK.density(pp, hlim=c(100, 10000), resolution = 1000)
  inf = rasterize(presence_pts, host, field = "nmbr_ps", fun='count', background=0)
  inf_values = getValues(inf)
  
  kd = sp.kde(presence_pts,y=as.numeric(presence_pts$nmbr_ps), bw=bw, host, standardize = T)
  kd_values = getValues(kd)
  
  kd = kd*max(inf_values) / max(kd_values, na.rm=T)
  names(kd)="Kernel_density_map"
  return(kd)
}
# surv_den can be a single numeric number or a map #
presence_absence_confidence = function(kernel_density_map, presence_pts, absence_pts, surv_density, sensi){
  
  #kernel_density_map = kernel_density(presence_pts)
  confidence_map = p_dtc(surv_density, kernel_density_map, sensi)
  
  presence_confi = extract(confidence_map, presence_pts)
  absence_confi = extract(confidence_map, absence_pts)
  
  presence_confi[is.na(presence_confi)]=0
  absence_confi[is.na(absence_confi)]=0
  
  presence_confi_mean = mean(presence_confi)
  absence_confi_mean = mean(absence_confi)
  
  confidence_map[kernel_density_map==0] = NA
  confidence_map[kernel_density_map>= 1 & is.na(confidence_map)] = 0
  all_confi = getValues(confidence_map)
  mean_confi = mean(all_confi, na.rm=T)
  
  ls = c(presence_confi_mean, absence_confi_mean, mean_confi)
  names(ls)=c("Mean_confi_presences", "Mean_confi_absences", "Mean_confi_allPixels")
  return(ls)
}


###### Step 5 -- Management Optimization and Eradication Feasibility ######
## This step need finer resolution data (e.g. 100-m)##

infestation_potential_rank = function(inf_map, host_map,  dispersal_scale, Wcoef=NA, range_buffer, number_core){
  
  # only counts susceptible hosts for the infestation potential 
  host_uninf = host_map
  host_uninf[inf_map>0]=0
  
  # covert infestation and host raster to points to get raster values and 
  # to easier calculate distance between infested pixel and susceptible host
  inf_pts = rasterToPoints(inf_map,fun=function(x){x>0} ,spatial=T)
  names(inf_pts)="inf_level"
  inf_pts$inf_level = inf_pts$inf_level 
  # extract weather_coefficient for each infested pixel/points     
  if(is.na(Wcoef)){Wcoef = host_map*0+1}
  wc_va = extract(Wcoef,inf_pts)
  
  host_pts = rasterToPoints(host_uninf, fun=function(x){x>0}, spatial=T)
  names(host_pts)="host"
  
  n = length(inf_pts)
  
  ## for each infested pixel, calculate the infestation potential for each susceptible host within the range_buffer
  ## and sum up
  cl = makeCluster(number_core)
  registerDoParallel(cl) 
  
  ip_va = foreach(i=1:n, .combine=rbind, .packages = c('raster','rgdal','geosphere','rgeos'))%dopar%{
    dist= gDistance(inf_pts[i,], host_pts,byid=T)
    
    dist2 = 1/(dist^2 + dispersal_scale^2)
    dist2[dist>range_buffer] = 0
    
    # so the calculation of infestation potential is: 
    infe_pot = inf_pts$inf_level[i]*wc_va[i]*dist2*host_pts$host
    
    infe_potS = sum(infe_pot)
    return(infe_potS)
  }
  
  stopCluster(cl)
  
  inf_pts$ip = as.vector(ip_va)
  inf_ply = rasterToPolygons(inf_map, fun=function(x){x>0},na.rm = T)
  inf_ply$ip = inf_pts$ip
  inf_ply2 = inf_ply[order(inf_ply$ip,decreasing = T),]
  
  return(inf_ply2)
}
ip_treat = function(ip_rank_ply, budget, buffer, cost_per_meter_sq, inf_map){
  
  area=budget/cost_per_meter_sq
  pixelArea=xres(inf)*yres(inf)
  
  
  ply_bf=buffer(ip_rank_ply,width=buffer,dissolve=F)
  n=floor(area/pixelArea)+2
  ply_bf$Cumu_Area=0
  
  if (length(ply_bf)< n){n = length(ply_bf)}
  
  # 
  for (i in 1:n){
    trt=gUnionCascaded(ply_bf[1:i,])
    ply_bf$Cumu_Area[i]=gArea(trt)
  }
  
  # select pixels together with the treatment buffer whoes total area is not larger than the budget allowed
  treatment=ply_bf[ply_bf$Cumu_Area<= area & ply_bf$Cumu_Area!=0,]
  treatment=gUnionCascaded(treatment)
  
  # if the selected treatment area is smaller than budget allowed, select part of the next infested pixel
  df=area-gArea(treatment)
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

contain_feasibility = function(inf_map, host_map,estimate_contain_budget =T,  buffer, Wcoef=NA, use_lethal_temperature, lethal_temperature,
                               lethal_temperature_month, infected, susceptible, total_plants,
                               mortality_on, mortality_tracker, mortality, treatment_maps = treat_list,
                               treatment_years, weather, temperature, weather_coefficient, ew_res,
                               ns_res, num_rows, num_cols, time_step, reproductive_rate,
                               mortality_rate = 0, mortality_time_lag = 2L,
                               season_month_start = 1L, season_month_end = 12L, start_time = 2018,
                               end_time = 2018, treatment_month = 12L, treatment_method = "ratio",
                               natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                               use_anthropogenic_kernel = FALSE, percent_natural_dispersal = 0,
                               natural_distance_scale = 21, anthropogenic_distance_scale = 0,
                               natural_dir = "NONE", natural_kappa = 0,
                               anthropogenic_dir = "NONE", anthropogenic_kappa = 0){
  
  current_inf_number = sum(getValues(inf/inf), na.rm = T)    
  
  ip_rank_ply = ip_rank(inf,host, natural_distance_scale,Wcoef,5000 ,number_core = 20)
  ip2 = ip_treat(ip_rank_ply, budget, buffer, cost_per_meter_sq, inf)
  
  df = data.frame(matrix(0, nrow=20,ncol=1))
  colnames(df)="number_inf"
  
  for (j in 1:20){
    random_seed= j*1000+50
    
    
    dataIp <- pops_model(random_seed, use_lethal_temperature, lethal_temperature,
                         lethal_temperature_month, infected, susceptible, total_plants,
                         mortality_on, mortality_tracker, mortality, treatment_maps,
                         treatment_years, weather, temperature, weather_coefficient, ew_res,
                         ns_res, num_rows, num_cols, time_step, reproductive_rate,
                         mortality_rate = 0, mortality_time_lag = 2L,
                         season_month_start = 1L, season_month_end = 12L, start_time = 2018,
                         end_time = 2018, treatment_month = 12L, treatment_method = "ratio",
                         natural_kernel_type = "cauchy", anthropogenic_kernel_type = "cauchy",
                         use_anthropogenic_kernel = FALSE, percent_natural_dispersal = 0,
                         natural_distance_scale = 21, anthropogenic_distance_scale = 0,
                         natural_dir = "NONE", natural_kappa = 0,
                         anthropogenic_dir = "NONE", anthropogenic_kappa = 0)
    
    df[j,1]=dataIp$area_infected[1]/pixelArea
    #df1[j,11]=dataIp$number_infected[1]
    
    #rate[j,1:4] = dataIp$rates[[1]]
    print(j)
  }
  
  inf_number_afterTrt = median(df$number_inf)
  if (inf_number_afterTrt <= current_inf_number){
    print("Current annual budget can contain the spread")
  } else {print("Current annual budget can not contain the spread")}
  
  
  if (estimate_contain_budget == T){
    budgets = unlist(lapply(1:10, function(x){budget + 5000000*x}))
    
    for (i in 1:length(budgets)){
      
      df = data.frame(matrix(0, nrow=20,ncol=1))
      colnames(df)="number_inf"
      
      budget_now = budgets[i]
      
      ip_rank_ply = ip_rank(inf_Ip,host, natural_distance_scale[[1]],Wcoef,rg ,20)
      ip2 = ip_treat(ip_rank_ply, budget, buffer, cost_per_meter_sq, inf)
      
      
      for (j in 1:20){
        random_seed= j*1000+50
        
        
        dataIp <- pops_model(random_seed = random_seed, 
                             use_lethal_temperature = use_lethal_temperature, 
                             lethal_temperature = lethal_temperature, 
                             lethal_temperature_month = lethal_temperature_month,
                             
                             use_movements=FALSE,
                             movements=list(0,0,0,0,0),
                             movements_dates= start_time,
                             
                             exposed=temperature[1:2],
                             model_type_="SEI",
                             latency_period=1,
                             
                             
                             infected = infected_species_Ip[[1]],
                             susceptible = susceptible_species_Ip[[1]],
                             total_plants = total_pl[[1]],
                             mortality_on = mortality_on,
                             mortality_tracker = infected_species_Ip[[1]]*0,
                             mortality = infected_species_Ip[[1]]*0,
                             treatment_maps = ip2,
                             
                             treatment_dates = c("2019-03-01"),
                             pesticide_duration=c(0),
                             resistant = infected_species_Ip[[1]]*0,
                             weather = weather,
                             temperature = temperature,
                             weather_coefficient = wc,
                             ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                             time_step = time_step, reproductive_rate = reproductive_rate[[1]],
                             mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                             season_month_start = season_month_start, season_month_end = season_month_end,
                             start_date = start_time, end_date = end_time,
                             treatment_method = "all infected", 
                             
                             natural_kernel_type = natural_kernel_type[[1]], anthropogenic_kernel_type = anthropogenic_kernel_type[[1]], 
                             use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal[[1]],
                             natural_distance_scale = natural_distance_scale[[1]], anthropogenic_distance_scale = anthropogenic_distance_scale[[1]], 
                             natural_dir = natural_dir[[1]], natural_kappa = natural_kappa[[1]],
                             anthropogenic_dir = anthropogenic_dir[[1]], anthropogenic_kappa = anthropogenic_kappa[[1]],output_frequency = "year")
        
        
        df[j,1]=dataIp$area_infected[1]/pixelArea
        #df1[j,11]=dataIp$number_infected[1]
        
        #rate[j,1:4] = dataIp$rates[[1]]
        print(j)
      }
      
      if (median(df$number_inf) <= current_inf_number){
        print(paste("Annual budget can contain the spread ", as.character(budget_now/1000000), " million can contain the spread." ))
        break  
      }
      
      
    }
    
  }
}