
# use follow code to install R package "concaveman"
#library(devtools)
#install_github("joelgombin/concaveman")


# Load packages
library(concaveman)
library(rgdal)
library(raster)
library(rgeos)
library(outliers)
library(rgdal)
library(PoPS)
library(spatialEco)

# set to the folder where you store the files attached
setwd("Q:/Shared drives/APHIS  Projects/eRADS/eRADS_workflow/eRADS_workflow/caseStudy")

####### Step 0 -- decide surveillance range to delineate infestation  ########
# Function 1&2 input data: first year point presences #
# Function 3 input data: first year point presences and estimation of surveillance sensitivity #

# read presence points #
load("us_County48.RData")
load("Roads.RData")
pts_dir = "./PreAbe"
pts_nm17 = "2017_Positive_ECFF_Trap_Sites"
pts_nm19 = "2019_Positive_ECFF_Trap_Sites"

pre_17 = readOGR(pts_dir,pts_nm17)
pre_17 = pts_format(pre_17)

pre_19 = readOGR(pts_dir,pts_nm19)
pre_19 = pts_format(pre_19)

# visulaize presence locations#
plot(pre_19, col="red")
plot(pre_17, add=T ,col="blue")
plot(ctys,add=T,border="blue")

#### Function 2 -- determine delineation range at the very beginning of invasion  ####

# coarse_estimation_of_annual_spread_rate in meter 
# buffer_size: size of buffer within which surveillance is desired around a presence in meter
# road_only: Boolean variable TRUE/FALSE, use true if only want to survey locations near to road
# road_buffer: float, when road_only = T, size of road buffer outside which survey wouldn't be conducted

sug_deli = delineate_range(pre_17, coarse_estimation_of_annual_spread_rate = 10000, buffer_size = 5000, road_only = T, road_buffer = 300)
ctd = sug_deli[[1]]
min_surv = sug_deli[[2]]
sug_surv = sug_deli[[3]]

# visualization of suggested range
plot(sug_surv,border="blue")
plot(min_surv, add=T,border="orange")
plot(ctd, add=T, pch=4, cex=1.5, col="red")
plot(pre_17, add=T)
plot(ctys, add=T, border="black")

#### Function 3 -- suggest surveillance density for delineation ####

# confi - float (0-1), desired detection confidence
# sensi - float (o-1), sensitivity of detection, be conservative on this estimation 
# presence_pts, spatialpointsdataframe, presence points of the invasive species 
surv_density = delineate_density(confi=0.95, presence_pts = pre_17, sensi=0.15)
surv_density


####### Step 1 -- estimate dispersal parameters based on first 2 or 3 years' sampling data  ########
# Function 1 input data: 1) host map 
# 2) folder path where all the presence points shapefiles are located
# each shapefile includes all presence points for one year, and the shapefile should be named properly so
# the list.file() function in R list the file for earliest year first and the latest year last

host_file = "./host.tif"
host = raster(host_file)
host = projectRaster(host, crs=crs(rds), res=1000)
inf_pts_folder="./presences/." # don't forget the dot at the end
inf_raster_folder = "./inf_rasters/" # don't forget slash at the end


#### Function 1 -- convert point presence to infestation raster ####
# inf_pts_folder: folder where the presence points (shapefile) are located, host: host raster data, 
# each shapefile includes all presence points for one year, and the shapefile should be named properly so
# the list.file() function in R list the file for earliest year first and the latest year last
# inf_raster_folder: folder path where you want to save the raster files converted from the points data

inf_st = pts_raster(inf_pts_folder = inf_pts_folder, host = host, inf_raster_folder = inf_raster_folder)
par(mfrow=c(1,2))
plot(raster(inf_st,1))
plot(raster(inf_st,2))

#### Function 2 -- estimate dispersal parameters for running PoPS model Using Chris' ABC Method####

# run the abc_calibration function to estimate parameters -- this would take hours or days to complete#

host_file = "./host.tif"
inf_raster_folder = "./inf_rasters/" # don't forget slash at the end
inf_files = list.files(inf_raster_folder,".tif$")
inf_files = paste(inf_raster_folder, inf_files,sep="")
infected_years_file = inf_files
number_of_observations = 600
prior_number_of_observations= 0.1
prior_means = c(20, 70, 0.97, 10000, 0, 0)
prior_cov_matrix = matrix(data=0.1, nrow=6, ncol=6)
infected_file = inf_files[1]
host = raster(host_file)
#host = projectRaster(host, crs=crs(rds), res=1000)

total_plants_file = "./hostTotal.tif"

parameters <- abc_calibration(infected_years_file, 
                              number_of_observations,
                              prior_number_of_observations,
                              prior_means, prior_cov_matrix, 
                              params_to_estimate = c(T, T, T, T, F, F),
                              number_of_generations = 2,
                              generation_size = 500, 
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


parameters
parameters$posterior_means



####### Step 2 -- estimate potential infestation intensity and infestation risk using PoPS ######
# input data: all data needed to run the pops_model

# this step use pops, need to check your PoPS version to set parameters; Basically this step use parameters and maps estiamted/calculated from step 1 to estimate 
# the potential invasion next year. 


####### Step 3 -- Surveillance optimization #######

risk_map_file = "./estimations/Risk_2019.tif"
risk_map = raster(risk_map_file)
density_map_file = "./estimations/Intensity_2019.tif"
density_map = raster(density_map_file)

### Note:  set road_only = T is only want to survey regions near to road, and 
### road_buffer parameter determines the pixels to survey, and pixels whose distance to road is further than this value wouldn't be surveyed
### sensi: surveillance sensitivity parameter
# budget: float - budget for surveillance
# cost_per_survey: float 
# risk_map: raster file
# risk_threshold: float -- regions with risk lower than this threshold won't be surveyed
# density_threshold: float --  regions with estimated infestation density lower than this threshold won't be surveyed


#### Method 1 -- Risk Based Strategy 1 -- Regular Surveillance ####

# function  use case
# risk_threhold paramter: pixels with predicted risk lower than this value wouldn't be surveyed

surv_results1 = surv_strategy1(budget= 1000000, cost_per_survey = 70, risk_map = risk_map, risk_threshold = 5, road_only = T, road_buffer = 300)
surv_density1 = surv_results1[[1]] # output is a float
risk_surv_Region1 = surv_results1[[2]] # output is a map
surv_density1
risk_surv_Region1

# function to select optimal confidence level and selection percentage for strategy 5
opticonfidence(density_map = density_map, sensi = 0.15, budget = 1000000, cost_per_survey = 70, confidenceMin = 0.3, confidenceMax = 0.9, number_core = 5 )


#### Method 2 -- Risk Based Strategy 2 -- Risk-Weighted Surveillance ####
# risk_level parameter: whether want to classiy the risk map into number_risk_level risk levels before conducting surveillance optmization
surv_results2 = surv_strategy2(budget= 1000000, cost_per_survey = 70, risk_map, risk_threshold = 5, road_only = T, road_buffer = 300, risk_level=T, number_risk_level = 5)
surv_density2 = surv_results2 # output is a map
surv_density2

#### Method 3 -- Infestation Density Based Strategy 1 -- Regular Surveillance ####
surv_results3 = surv_strategy3(budget= 1000000, cost_per_survey = 70, density_map, density_threshold = 2, road_only = T, road_buffer= 300)
surv_density3 = surv_results3[[1]] # output is a float
risk_surv_Region3 = surv_results3[[2]] # output is a map
surv_density3


#### Method 4 -- Infestation Density Based Strategy 2 -- Same Detection Confidence Method####
surv_results4 = surv_strategy4(budget= 1000000, cost_per_survey = 70, density_map, density_threshold = 2, road_only = T, road_buffer= 300, sensi= 0.15)
detc_confidence = surv_results4[[1]] # output is a float
surv_density_map4 = surv_results4[[2]] # output is a map
detc_confidence
plot(surv_density_map4)

#### Method 5 -- Infestation Density Based Strategy 3 -- Prioritize percentage locations with higher infestation density ####
# location_percentage = T -- select  percentage of pixels with highest estimated density
# location_percentage = F -- select pixels with highest estimated density which in total have percentage of all estimated infestation population

surv_results5 = surv_strategy5(budget= 1000000, cost_per_survey = 70, density_map = density_map, density_threshold = 2, road_only = T, road_buffer= 300, percentage =  0.95, location_percentage = T ,sensi= 0.15)
surv_results5 
density_map5 = surv_results5[[1]] # output is a float
regular_surv_den = surv_results5[[2]]
same_surv_confi = surv_results5[[3]]
same_surv_confi_density = surv_results5[[4]]

#### Method 6 -- Infestation Density Based Strategy 4 -- Survey location with estimated infestation density higher than a given threshold ####
surv_results6=surv_strategy6(budget = 1000000, cost_per_survey = 70, density_map = density_map, density_threshold = 1, road_only=F, road_buffer=0, sensi=0.15)
surv_results6


#### Method 7 -- Risk and Infestation Density Based Strategy -- Survey location with higher infestation potential (i.e. potential of causing more infestation) ####
host_map_file = "./host.tif"
host_map = raster(host_map_file)
host = host_map
range_buffer = 5000 

surv_results7=surv_strategy7(budget = budget, cost_per_survey = cost_per_survey, risk_map = risk_map, density_map = density_map, percentage_location=0.9, road_only=F, road_buffer=0, sensi=0.15)
surv_results7


####### Step 4 -- Evaluate efficacy  of each method ######
#setwd("../eRADS_workflow/eRADS_workflow/caseStudy/")
budgets = seq(500000, 1500000, by=500000)
risk_map_file = "./estimations/Risk_2019.tif"
risk_map = raster(risk_map_file)
density_map_file = "./estimations/Intensity_2019.tif"
density_map = raster(density_map_file)
risk_map = projectRaster(risk_map, crs=crs(rds), res=1000)
density_map = projectRaster(density_map, crs=crs(rds), res=1000)
host_map = raster("./host.tif")
cost_per_survey = 70
sensi= 0.17

df1 = data.frame(matrix(0, ncol=11, nrow=length(budgets)))
colnames(df1) = c("Budget","Risk_S1", "Risk_S2", "Dens_S3_reg", "Dens_S4_sameConfi", "Dens_S5_Perct_reg","Dens_S5_Perct_sameConfi", 
                  "Dens_S6_thresh_reg","Dens_S6_thresh_sameConfi", "Dens_S7_IP_reg","Dens_S7_IP_sameConfi")

colnames(df1) = c("Sensi","Risk_S1", "Risk_S2", "Dens_S3_reg", "Dens_S4_sameConfi", "Dens_S5_Perct_reg","Dens_S5_Perct_sameConfi", 
                  "Dens_S6_thresh_reg","Dens_S6_thresh_sameConfi", "Dens_S7_IP_reg","Dens_S7_IP_sameConfi")
sensis = c(0.1, 0.3, 0.5, 0.7)
for (i in 1:length(sensis)){
  budget = budgets[2]
  sensi = sensis
  # strategy 1 
  surv_results1 = surv_strategy1(budget, cost_per_survey, risk_map, risk_threshold= 5, road_only = F, 300)
  surv_density1 = surv_results1[[1]] # output is a float
  confi1 = strategy_evaluation(density_map, density_threshold=0.2, surv_density1, sensi)
  
  # strategy 2
  surv_results2 = surv_strategy2(budget, cost_per_survey, risk_map, risk_threshold= 5, road_only = F, 300, risk_level = T ,number_risk_level=10)
  surv_density2 = surv_results2
  confi2 = strategy_evaluation(density_map, density_threshold=0.2,surv_density2, sensi)
  
  # strategy 3
  surv_results3 = surv_strategy3(budget, cost_per_survey, density_map, density_threshold=0.2, road_only = F, 500)
  surv_density3 = surv_results3[[1]] # output is a float
  confi3 = strategy_evaluation(density_map,density_threshold=0.2, surv_density3, sensi)
  
  # strategy 4
  surv_results4 = surv_strategy4(use_budget=T, detection_confidence=0.9 ,budget, cost_per_survey, density_map, density_threshold=0.2, road_only = F, road_buffer = 500, sensi)
  surv_density4 = surv_results4[[2]] # this method directly outputs detection confidence
  confi4 = strategy_evaluation(density_map,density_threshold=0.2, surv_density4, sensi)
  
  # strategy 5
  surv_results5 = surv_strategy5(use_budget=T, detection_confidence=0.9 ,budget, cost_per_survey, density_map, density_threshold=0.2, road_only = F, 300, percentage=0.9, location_percentage = T , sensi)
  regular_surv_den5 = surv_results5[[2]]
  same_surv_confi_den5 = surv_results5[[4]]
  confi5_reg = strategy_evaluation(density_map, density_threshold=0.2,regular_surv_den5, sensi)
  confi5_sameConfi = strategy_evaluation(density_map, density_threshold=0.2,same_surv_confi_den5, sensi)
  
  
  # strategy 6
  surv_results6 = surv_strategy6(use_budget=T, detection_confidence=0.9 ,budget, cost_per_survey, density_map, density_threshold = 0.2, road_only=F, road_buffer=0, sensi=sensi)
  regular_surv_den6 = surv_results6[[3]]
  same_surv_confi_den6 = surv_results6[[5]]
  confi6_reg = strategy_evaluation(density_map,density_threshold=0.5, regular_surv_den6, sensi)
  confi6_sameConfi = strategy_evaluation(density_map, density_threshold=0.5,same_surv_confi_den6, sensi)
  
  # strategy 7
  surv_results7 = surv_strategy7(use_budget=T, detection_confidence=0.9 ,budget, cost_per_survey, risk_map, density_map, density_threshold = 0.2, percentage_location=0.9, road_only=F, road_buffer=0, sensi=0.15)
  regular_surv_den7 = surv_results7[[2]]
  same_surv_confi_den7 = surv_results7[[4]]
  confi7_reg = strategy_evaluation(density_map, density_threshold=0.2,regular_surv_den7, sensi)
  confi7_sameConfi = strategy_evaluation(density_map,density_threshold=0.2, same_surv_confi_den7, sensi)
  
  df1[i, 2:11] = c(confi1, confi2, confi3, confi4, confi5_reg, confi5_sameConfi, confi6_reg, confi6_sameConfi,confi7_reg, confi7_sameConfi)
  df1[i, 1] = budget
}
#df1$Sensi=sensis
df1$Budget = budgets

## Visualization 
library(ggplot2)
library(reshape2)
par(mfrow=c(1,1))
df1_re = melt(df1, id.vars="Budget", variable.name="strategies")
head(df1_re)
ggplot(df1_re, aes(Budget, value)) + geom_line(aes(colour = strategies, linetype = strategies),size=1) +
 geom_point(aes(colour = strategies)) +theme_classic() 


#### Evaluate on the estimated potential infestation intensity map (can be conducted before getting the sample data) ####
# use this function to get the detection confidence curve of all strategies with different bugget level #
#### Evaluate on the kernel density map derived from actual sample data ####

# use case -- evaluate all strategies with different budget level #
presence_pts_folder = "./evaluations/." 
presence_pts_name = "ECFF_2019"
absence_pts_folder = "./evaluations/."
absence_pts_name = "y19absence"

risk_map_file = "./estimations/Risk_2019.tif"
risk_map = raster(risk_map_file)
density_map_file = "./estimations/Intensity_2019.tif"
density_map = raster(density_map_file)

presence_pts = readOGR(presence_pts_folder, presence_pts_name)
presence_pts = spTransform(presence_pts, crs(rds))
absence_pts = readOGR(absence_pts_folder, absence_pts_name)
absence_pts = spTransform(absence_pts, crs(rds))


kernel_density_map = kernel_density(presence_pts)

budgets = seq(500000, 1500000, by=500000)
cost_per_survey = 70
sensi= 0.17

df_meanAll = data.frame(matrix(0, ncol=11, nrow=length(budgets)))
colnames(df_meanAll) = c("Budget","Risk_S1", "Risk_S2", "Dens_S3_reg", "Dens_S4_sameConfi", "Dens_S5_Perct_reg","Dens_S5_Perct_sameConfi", 
                         "Dens_S6_thresh_reg","Dens_S6_thresh_sameConfi", "Dens_S6_IP_reg","Dens_S6_IP_sameConfi")
df_pre = data.frame(matrix(0, ncol=11, nrow=length(budgets)))
colnames(df_pre) = c("Budget","Risk_S1", "Risk_S2", "Dens_S3_reg", "Dens_S4_sameConfi", "Dens_S5_Perct_reg","Dens_S5_Perct_sameConfi", 
                     "Dens_S6_thresh_reg","Dens_S6_thresh_sameConfi", "Dens_S6_IP_reg","Dens_S6_IP_sameConfi")
df_abe = data.frame(matrix(0, ncol=11, nrow=length(budgets)))
colnames(df_abe) = c("Budget","Risk_S1", "Risk_S2", "Dens_S3_reg", "Dens_S4_sameConfi", "Dens_S5_Perct_reg","Dens_S5_Perct_sameConfi", 
                     "Dens_S6_thresh_reg","Dens_S6_thresh_sameConfi", "Dens_S6_IP_reg","Dens_S6_IP_sameConfi")

for (i in 1:length(budgets)){
  budget = budgets[i]
  
  # strategy 1 
  surv_results1 = surv_strategy1(budget, cost_per_survey, risk_map, risk_threshold= 5, road_only = F, 300)
  surv_density1 = surv_results1[[1]] # output is a float
  kernel_density_map[kernel_density_map<0.2]=0
  confi1 = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, surv_density1, sensi)
  
  # strategy 2
  surv_results2 = surv_strategy2(budget, cost_per_survey, risk_map, risk_threshold= 5, road_only = F, 300, risk_level = T ,number_risk_level=10)
  surv_density2 = surv_results2
  confi2 = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, surv_density2, sensi)
  
  # strategy 3
  surv_results3 = surv_strategy3(budget, cost_per_survey, density_map, density_threshold=0, road_only = F, 500)
  surv_density3 = surv_results3[[1]] # output is a float
  confi3 = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, surv_density3, sensi)
  
  # strategy 4
  surv_results4 = surv_strategy4(budget, cost_per_survey, density_map, density_threshold=0, road_only = F, 500, sensi)
  surv_density4 = surv_results4[[2]] # this method directly outputs detection confidence
  confi4 = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, surv_density4, sensi)
  
  # strategy 5
  surv_results5 = surv_strategy5(budget, cost_per_survey, density_map, density_threshold=0, road_only = F, 300, percentage=0.9, location_percentage = T , sensi)
  regular_surv_den5 = surv_results5[[2]]
  same_surv_confi_den5 = surv_results5[[4]]
  confi5_reg = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, regular_surv_den5, sensi)
  confi5_sameConfi = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, same_surv_confi_den5, sensi)
  
  
  # strategy 6
  surv_results6=surv_strategy6(budget, cost_per_survey, density_map, density_threshold = 0.5, road_only=F, road_buffer=0, sensi=sensi)
  regular_surv_den6 = surv_results6[[3]]
  same_surv_confi_den6 = surv_results6[[5]]
  confi6_reg = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, regular_surv_den6, sensi)
  confi6_sameConfi = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, same_surv_confi_den6, sensi)
  
  # strategy 7
  surv_results7=surv_strategy7(budget, cost_per_survey, risk_map, density_map, percentage_location=0.9, road_only=F, road_buffer=0, sensi=0.15)
  regular_surv_den7 = surv_results7[[2]]
  same_surv_confi_den7 = surv_results7[[4]]
  confi7_reg = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, regular_surv_den7, sensi)
  confi7_sameConfi = presence_absence_confidence(kernel_density_map, presence_pts, absence_pts, same_surv_confi_den7, sensi)
  
  df_meanAll[i, 2:11] = c(confi1[3], confi2[3], confi3[3], confi4[3], confi5_reg[3], confi5_sameConfi[3], confi6_reg[3], confi6_sameConfi[3], confi7_reg[3], confi7_sameConfi[3])
  df_meanAll[i, 1] = budget
  
  df_pre[i, 2:11] = c(confi1[1], confi2[1], confi3[1], confi4[1], confi5_reg[1], confi5_sameConfi[1], confi6_reg[1], confi6_sameConfi[1], confi7_reg[1], confi7_sameConfi[1])
  df_pre[i, 1] = budget
  
  df_abe[i, 2:11] = c(confi1[2], confi2[2], confi3[2], confi4[2], confi5_reg[2], confi5_sameConfi[2], confi6_reg[2], confi6_sameConfi[2], confi7_reg[2], confi7_sameConfi[2])
  df_abe[i, 1] = budget
}
df_meanAll
df_pre
df_abe

par(mfrow=c(1,1))
df_pre_re = melt(df_pre, id.vars="Budget", variable.name="strategies")
head(df_pre_re)
ggplot(df_pre_re, aes(Budget, value)) + geom_line(aes(colour = strategies, linetype= strategies), size=1) +
  geom_point(aes(colour = strategies, shape = strategies)) + theme_classic() 



###### Step 5 -- Management Optimization and Eradication Feasibility ######
# data are 100-m resolution #
## this step takes a long time to run the simulations -- this step also use PoPS_model, need to prepar data to run PoPS##
inf_file = "./feasibility/inf.tif"
host_file = "./feasibility/host.tif"
inf_map = raster(inf_file)
host_map = raster(host_file)
dispersal_scale = 55.55
cost_per_meter_sq = 1
buffer = 0
budget = 1000000
#### convert infestation raster into polygons (one polygon covers one infestation pixel), and rank all polygons based on the infestation potential ####
### the return is polygons and the first polygon has highest infestation potential ###
# Wcoef - mean weather coefficient raster 
# range_buffer - the size of buffer in meter, it's assumed that each infested pixel have little impact on locations outside of this buffer in the following dispersal year. 
# people can use the coarse estimation of annual local dispersal rate for this parameter 

ip_rank_ply = infestation_potential_rank(inf_map = inf_map, host_map = host_map, dispersal_scale = dispersal_scale, Wcoef = NA, range_buffer = 2000, number_core = 5)


# select locations for treatment based on the infestation potential of each pixels, the return is a treatment map list
treat_list = ip_treat(ip_rank_ply = ip_rank_ply, budget = budget, buffer = buffer, cost_per_meter_sq = cost_per_meter_sq, inf_map = inf_map)
sum(unlist(treat_list))

# assess feasibility to decrease the infestation, if  estimate_contain_budget = T, also estimate the budget needed to decrease the infestation, most of the
# parameters here are the parameters needed to run pops
contain_feasibility(inf_map, host_map, estimate_contain_budget =F,  buffer=0, Wcoef=NA, use_lethal_temperature, lethal_temperature,
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
                    anthropogenic_dir = "NONE", anthropogenic_kappa = 0)
