
#### **************Contact wliang5@ncsu.edu anytime for questions, concerns, or comments *******************####



#### Step 0 Functions -- Delineate Surveillance Range ####

#### pts_format ####
# this function will only keep presence points in the United States,
# and reprojects the points to North America Albers Equal Area Conic projection
#' @return a spatial points dataframe
#' @param presence_pts  spatial points where the pest were detected
# important notice: it's assumed that the abundance of pest for each point is saved in the column named "nmbr_ps"

#### raster_format ####
# this function reproject raster data into North America Albers Equal Area Conic projection and 1000m resolution
#' @return a raster map
#' @param raster_map  a raster map to be reprojected


#### delineate_range ####
# this function suggests a minimum surveillance range  #
# and a larger surveillance range to delineate infestation based on current record 
#' @return  a list with 3 elements: 1st element is the centroid, 2nd and 3rd elements are minimum surveillance range and suggested surveillance range, respectively.
#' # important notice: it's assumed that the abundance of pest for each point is saved in the column named "nmbr_ps"
#' @param coarse_estimation_of_annual_spread_rate in meter. If based on current infestation status, 
#' # it's likely that the invasive species was already present last year, it's suggested that
#' #   the radius of the surveillance range can be the annual spread rate to account for the potential dispersal
#' @param  buffer_size in meter size of buffer within which surveillance is desired around a presence
#' @param  road_only T (TRUE) or F (FALSE), only survey regions along the road if T
#' @param road_buffer in meter, when road_only = T, the surveillance will only be conducted in regions whose distance 
#'   to road is not further than road_buffer

#### delineate_density ####
# this function suggest surveillance density for delineating the infestation based on first year infestation record,
# however, the suggested surveillance density is likely higher than needed as the reported infestation density 
# is likely lower than the real density for the first year
#' @return a single number, which is the suggested surveillance density
#' @param confi: desired detection confidence, scale 0-0.999 
#' @param presence_pts  spatial points where the pest were detected
#' # important notice: it's assumed that the abundance of pest for each point is saved in the column named "nmbr_ps"
#' @param sensi: survey sensitivity


#### Step 1 Functions ####

#### pts_raster ####
# this function covert presence points into raster infestation data #
#' @return multiple raster file with each one represent one years infestation
#' @param inf_pts_folder folder path where all the presence points shapefiles are located
# each shapefile includes all presence points for one year, and the shapefile should be named properly so
# the list.file() function in R list the file for earliest year first and the latest year last
# again it's assumed that the abundance of pest for each point is saved in the column named "nmbr_ps"
#' @param host host map and the generated annual infestation raster will have the 
#' same projection, extent, and resolution as the host map
#' @param inf_raster_folder folder path where the converted infestation raster will be saved

#### abc_calibration ####
# refer to PoPS package abc_calibration function


#### Step 2 Functions ####
#### infestation_estimation ####
# this function estimate the potential infestation intensity and risk using pops_model from PoPS pacakge
# see pops_model from PoPS for details on parameters


#### Step 3 Functions ####

#### p_dtc  ####
# this function calculates the detection probability for a spatial unit given the surveillance density, 
# population density, and sensitivity per survey
#' @return a raster map or a single value (based on input), indicating the estimated detection confidence
#' @param surv_den surveillance density, can be a single  number (double) or a raster map
#' @param p_den population density per spatial unit, can be single number (double) or a raster map
#' @param sensi surveillance sensitivity, value between 0-1

#### surv_den ####
# this function calculate the surveillance density to achieve desired detection confidence for given population density and surveillance sensitivity
#' @return a raster map or a single value (based on input), indicating the suggested surveillance density
#' @param confi desired detection confidence, this can be single number (double) or a raster map 
#' @param p_den population density per spatial unit, can be single number (double) or a raster map
#' @param sensi surveillance sensitivity, value between 0-1

### percent_location, this function is designed to be used by surv_strategy5 ###
# this function select the given percentage of highest (predicted) infested locations 
#' @return a raster map, the same as input parameter density_map except it sets the value of non-selected pixels to 0
#' @param density_map estimated infestation/population density raster map
#' @param percentage double between 0-1, percentage of estimated invaded locations 

### percent_population, this function is designed to be used by surv_strategy5 ###
# this function select the highest (predicted) infested locations which in total include given percentage of all (predicted) invaded populations
#' @return a raster map, the same as input parameter density_map except it sets the value of non-selected pixels to 0
#' @param density_map estimated infestation/population density raster map
#' @param percentage double between 0-1, percentage of total invaded populations

### confidence_all & opticonfidence ###
### confidence_all, designed to  be used internally by the below opticonfidence function### 
### opticonfidence, designed to be used optionally by surv_strategy5 to calculate the optimal detection confidence 
# for selected surveillance pixels when using same detection confidence strategies ###
#' @return a list with 2 elements, the first element includes optimal detection confidence and location percentage for selection of percentage location,
#'        # and the 2nd element includes the same information for selection based on population
#' @param density_map estimated  infestation/population density raster map
#' @param budget, double, budget in dollar
#' @param cost_per_survey, double, estimation of cost per survey site per year, this cost includes all survey cost such as trap cost, labor cost, and facility cost
#' @param sensi, a single value between 0-1, surveillance sensitivity
#' @param confidenceMin & @param confidenceMax, float, 0.1 - 0.999, as the opticonfidence function explore all possibly desired detection confidence for each pixel, 
#' and then return the confidence which can give highest overall detection confidence over all infested area, these two parameters are the user
#' accepted minimum and maximum detection confidence for individual survyed pixels
#' @param numer_core, integer, number of core to be used to calculate the optimal confidence


#### infestation_potential, designed for surv_strategy7 ####
# this function calculate the infestation potential (IP) of each (predicted) infestation pixel and return an IP raster map for 
# for given percentage of (predicted) infested pixels with highest IP
#' @return a raster map with raster value indicates the IP of the infested pixel
#' @param density_map estimated  infestation/population density raster map
#' @param host_map raster map for host
#' @param dispersal_scale the natural dispersal scale, the same as used in pops
#' @param Wcoef mean weather coefficient raster map, if NA, the weather coefficient will be set to 1 for all locations
#' @param range_buffer  in meter, a single number; for a given (estimated) invaded pixel,
#'  susceptible host within this buffer is considered as potentially can be invaded because of the foci,
#' local dispersal rate can be used as the buffer size; experiment showed 2000 or 5000 can be a good value
#' @param percentage_location a single value between 0-1, the output IP raster only includes IP values for the given percentage of pixels,  please note this method prioritize locations with higher IP
#' @param number_core number of cores used to calculate IP


#### Surv_strategy 1-7  ####
    ##*************please note the default map resolution is 1km **************##
    ##************ essentially, surv_strategy3  + surv_strategy4 is surv_strategy6 if use the same density_threshold 
#' @param budget, double, budget in dollar
#' @param cost_per_survey, double, estimation of cost per survey site per year, this cost includes all survey cost such as trap cost, labor cost, and facility cost
#' @param road_only, value T or F, if T only survey pixels near to road
#' @param road_buffer, a single numeric value in meter, if road_only=T, only pixels whose minimum distance to roads is no larger than the value of road_buffer will be surveyed

## Strategy 1-2 ##
#' @param risk_map, raster map,value between 0-100, for strategy 1 and strategy 2, estimated invasion risk map
#' @param risk_threshold, a single value between 0-100, for strategy 1 and strategy 2, pixels with risk lower than this value will be ignored
#' @param risk_level, T or F, parameter for strategy 2, if T, the risk map will be reclassified into a risk_level map 
#' @param number_risk_level, an integer, parameter for strategy 2, if risk_level=T, the risk map will be reclassified into a risk_level map with number_risk_level levels of risk

#' @param density_map, raster map, estimated infestation density map, for strategy 3-7
#' @param density_threshold, a single value between 0-100, pixels with estimated infestation lower than this value will be ignored
#' @param sensi, a value between 0-1, per surveillance sensitivity, which is defined as the probability of detecting an invasion
#'        # when the surveillance density is 1 and invasion/infestation density is 1   

## Strategy 4-7 ##
#' @param use_budget, boolean, T or F, default T
#' @param detection_confidence, a single value between 0-1, if use_budget=F, the desired detection confidence for surveillance pixels 

## Strategy 5 ##
#' @param percentage, a value between 0-1, used for surv_strategy5,percentage of total infested locations or total infestation populations selected for surveillance
#' @param location_percentage, Boolean variable, T or F, used for surv_strategy5, if T, the strategy will select the given percentage of estiamted infested locations for surveillance,
#'          if F, the strategy will select locations which in total have given percentage of estimated infestation populations
#'          either way, this strategy prioritize locations with higher estimated infestation density/intensity      
## comments: based on earlier experiences, given the same percentage value for percent_location and percent_population, the percent_location function usually select more locations 

## Strategy 7 ##
#' @param dispersal_scale the natural dispersal scale, the same as used in pops
#' @param Wcoef mean weather coefficient raster map, if NA, the weather coefficient will be set to 1 for all locations
#' @param range_buffer  in meter, a single number; for a given (estimated) invaded pixel,
#'  susceptible host within this buffer is considered as potentially can be invaded because of the foci,
#' local dispersal rate can be used as the buffer size; experiment showed 2000 or 5000 can be a good value
#' @param percentage_location a single value between 0-1, the output IP raster only includes IP values for the given percentage of pixels,  please note this method prioritize locations with higher IP
#' @param number_core number of cores used to calculate IP

## returns of strategy1-7 ##
# Strategy 1 #
#' @return a list with 2 elements: 1st element is a single number indicating surveillance density,
#' 2nd element return the risk raster map (same as input risk_map) but only includes non-zero values for selected surveillance locations
# Strategy 2 #
#' @return a surveillance density map
# Strategy 3 #
#' @return a list with 2 elements: 1st element is a raster map indicating surveillance density,
#' 2nd element return the infestation density  map (same as input density_map) but only includes non-zero values for selected surveillance locations

# Return of Strategy 4-7 if use_budget = F #
#' @return a list with 2 elements: 1st element is a single number indicating the total number of surveillance (Trap) needed,
#'         and the 2nd element is a surveillance density raster map   

# Return of Strategy 4-7 if use_budget = T #      
# Strategy 4 #
#' @return a list with 2 elements: 1st element is a single value indicating the detection confidence for surveillance pixels,
#' 2nd element return the surveillance density map 

# Strategy 5 #
# if location_percentage = T
#' @return a list with 4 elements: 1st element is the infestation density  map (same as input density_map) but only includes non-zero values for selected surveillance locations
#'                                 2nd element is a single number - the surveillance density for regular surveillance strategy
#'                                 3rd element is a single number -detection confidence for same detection confidence strategy
#'                                 4th element is a raster map - surveillance density map for same detection confidence strategy  
# if location_percentage = F
#' @return a list with 5 elements: 1st element is a single number - percentage of pixels selected for surveillance
#'                                 2nd element is the infestation density  map (same as input density_map) but only includes non-zero values for selected surveillance locations
#'                                 3rd element is a single number - the surveillance density for regular surveillance strategy
#'                                 4th element is a single number -detection confidence for same detection confidence strategy
#'                                 5th element is a raster map - surveillance density map for same detection confidence strategy  

# Strategy 6 #
#' @return a list with 5 elements: 1st element is a single number - percentage of pixels selected for surveillance
#'                                 2nd element is the infestation density  map (same as input density_map) but only includes non-zero values for selected surveillance locations
#'                                 3rd element is a single number - the surveillance density for regular surveillance strategy
#'                                 4th element is a single number -detection confidence for same detection confidence strategy
#'                                 5th element is a raster map - surveillance density map for same detection confidence strategy  

# Strategy 7 #
#' @return a list with 4 elements: 1st element is the infestation density  map (same as input density_map) but only includes non-zero values for selected surveillance locations
#'                                 2nd element is a single number - the surveillance density for regular surveillance strategy
#'                                 3rd element is a single number -detection confidence for same detection confidence strategy
#'                                 4th element is a raster map - surveillance density map for same detection confidence strategy  


#### Step 4 -- Evaluate efficacy  of each method Functions ####

### strategy_evaluation ###
# this function conducts evaluation on the estimated potential infestation intensity map (can be conducted before getting the next year sample data) #
# use this function to get the detection confidence curve of all strategies with different budget levels or sensitivity values #
#' @return a single value - mean detection confidence over the whole infested area
#' @param density_map, raster map, estimated infestation density map, for strategy 3-7
#' @param density_threshold, a single value between 0-100, pixels with estimated infestation lower than this value will be ignored
            # please note that this value should keep the same as the value used for estimating the surveillance density value or map
#' @param surv_den, a single value or a raster map, surveillance density map, this should be an output from surv_strategy functions
#' @param sensi, a single value between 0-1, surveillance sensitivity

###  kernel_density ###
# this function convert the input spatial presence point to a kernel density raster map, which is an estimation of invasion density over the geospatial extent of the input point data
#' @return a raster map estimating the possible infestation density 
#' @param presence_pts  spatial points where the pest were detected
# important notice: it's assumed that the abundance of pest for each point is saved in the column named "nmbr_ps"

### presence_absence_confidence ###
# this function calculate the mean detection confidence for the field sampling points, could be presence points or absence points, or both    
#' @returna list with 3 elements: mean detection confidence for presence locations, mean detection confidence for absence locations, and mean detection confidence for the whole invaded areas based on the kernel density map
#' @param kernel_density_map output from the kernel_density function
#' @param presence_pts  spatial points where the pest were detected
# important notice: it's assumed that the abundance of pest for each point is saved in the column named "nmbr_ps"
#' @param absence_pts  spatial points where the pest were not detected
#' @param surv_den, a single value or a raster map, surveillance density map, this should be an output from surv_strategy functions
#' @param sensi, a single value between 0-1, surveillance sensitivity


#### Step 5 FUnctions ####
## ** Functions in this step take a long time to run **##

### infestation_potential_rank ###
# this function calculate the infestation potential (ip) of each infested pixel, and 
#' @return a spatial polygon dataframe, and the polygons in the dataframe are ranked decreasingly based on their ip value
# each polygon covers one infested pixel
#' @param inf_map, an infestation raster map
#' @param host_map, a raster host map
#' @param dispersal_scale the natural dispersal scale, the same as used in pops
#' @param Wcoef mean weather coefficient raster map, if NA, the weather coefficient will be set to 1 for all locations
#' @param range_buffer  in meter, a single number; for a given (estimated) invaded pixel,
#'  susceptible host within this buffer is considered as potentially can be invaded because of the foci,
#' local dispersal rate can be used as the buffer size; experiment showed 2000 or 5000 can be a good value
#' @param number_core number of cores used to calculate IP

### ip_treat ###
# this function selects locations for treatment and prioritizes pixels with higher ip values
#' @return  a treatment list which can be used as the treatment map for pops_model
#' @param ip_rank_ply a spatial polygon dataframe, the return from the infestation_potential_rank function
#' @param budget Double, budget in dollar for treatment
#' @param buffer Double, buffer size in meter, this is the buffer around selected pixels for treatment which will also be tretment
#' @param cost_per_meter_sq Double, cost for treating per square meter
#' @param inf_map infestation raster map

### contain_feasibility ###
# this function assess whether a given budget can decrease the infestation or not, 
# if not, it can also evaluate how much budget can decrease the infestation if set estimate_contain_budget = T
#' @return a sentence stating if the containment is possible or not; 
#'   if estimate_contain_budget = T, it also returns an estimation of how much budget to decrease infestation if the containment is not possible
#' @param inf_map infestation raster map
#' @param host_map, a raster host map
#' @param estimate_contain_budget, Boolean T or F, if T, return an estimation on how much budget on treatment can decrease current invasion
#' @param buffer Double, buffer size in meter, this is the buffer around selected pixels for treatment which will also be treated
# for all other parameters, refer to pops_model in PoPS package



