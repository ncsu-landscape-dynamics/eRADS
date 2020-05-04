######################## 2. Estimate Possible Infestation Intensity #######

library(spatstat)
library(spatialEco)
setwd("Q:/Shared drives/APHIS  Projects/eRADS/Biweekly meetings/EGVM/")
ctd=readOGR(".","Y09Ctd2")
pts=readOGR(".","EGVM_proj")
ra=raster(extent(pts),crs=crs(pts),resolution=1000,vals=1)
pts_ra=rasterize(SpatialPoints(pts),ra,fun="count")

#### approach 1 -- kernel density (very first year of invasion)  ####
kd=sp.kde(pts,y=1,bw=5000,ra)
sum(getValues(kd),na.rm = T)
plot(kd)
plot(kd,ext=c(500000,600000,4200000,4300000))
plot(pts,add=T,col="blue")


#### approach 2 -- pops simulation (following years of invasion) ####
inf=kd
no=list(as.matrix(inf*0))
inf_Itensi=stack()
inf_Risk=stack()

for (i in 1:100){
  
  random_seed=i*2901+20
  dataNo <- pops_model(random_seed = random_seed, 
                       use_lethal_temperature = use_lethal_temperature, 
                       lethal_temperature = lethal_temperature, 
                       lethal_temperature_month = lethal_temperature_month,
                       
                       infected = infected_speci[[1]],
                       susceptible = susceptible_speci[[1]],
                       total_plants = total_pl[[1]],
                       mortality_on = mortality_on,
                       mortality_tracker = infected_speci[[1]]*0,
                       mortality = infected_speci[[1]]*0,
                       treatment_maps = no,
                       
                       treatment_dates = c("2019-03-01"),
                       pesticide_duration=c(0),
                       resistant = infected_speci[[1]]*0,
                       weather = weather,
                       temperature = temperature,
                       weather_coefficient = wc,
                       ew_res = ew_res, ns_res = ns_res, num_rows = num_rows, num_cols = num_cols,
                       time_step = time_step, reproductive_rate = reproductive_rate[[2]],
                       mortality_rate = mortality_rate, mortality_time_lag = mortality_time_lag,
                       season_month_start = season_month_start, season_month_end = season_month_end,
                       start_date = start_time, end_date = end_time,
                       treatment_method = "all infected", 
                       
                       natural_kernel_type = natural_kernel_type[[1]], anthropogenic_kernel_type = anthropogenic_kernel_type[[1]], 
                       use_anthropogenic_kernel = use_anthropogenic_kernel, percent_natural_dispersal = percent_natural_dispersal[[1]],
                       natural_distance_scale = natural_distance_scale[[2]], anthropogenic_distance_scale = anthropogenic_distance_scale[[1]], 
                       natural_dir = natural_dir[[1]], natural_kappa = natural_kappa[[1]],
                       anthropogenic_dir = anthropogenic_dir[[1]], anthropogenic_kappa = anthropogenic_kappa[[1]],output_frequency = "year")
  
  infNo=dataNo$infected
  infNo=as.data.frame(infNo)
  infNo=as(infNo,"matrix")
  infNo=raster(infNo,template=inf)
  inf_Itensi=stack(inf_Itensi,infNo)
  
  infNo=reclassify(infNo,c(0,1,0))
  infNo=reclassify(infNo,c(1,Inf,1))
  inf_Risk=stack(inf_Risk,infNo)
}
inf_more=mean(inf_Itensi)
inf_more=round(reclassify(inf_more,c(0,1,0)))

ply0=ra2plyMerge(inf)
ply1=ra2plyMerge(inf_more)

gdiff=gDifference(ply1,ply0)
infnew=crop(inf_more,gdiff)
infnew=mask(infnew,gdiff)





