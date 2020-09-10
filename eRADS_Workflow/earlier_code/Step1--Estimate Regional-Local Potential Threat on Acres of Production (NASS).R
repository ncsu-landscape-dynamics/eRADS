######################## 0. Estimate Regional/Local Potential Production Impact ###############

setwd("Q:/Shared drives/APHIS  Projects/shared resources/data/USDA_NASS_Data")
#NASS=read.csv("2017USDA_NASS_Census_Data.csv")
load("NASS.RData")
dim(NASS)
head(NASS)

states="CALIFORNIA"
croptypes=c("GRAPES")
production_types=c("GRAPES - ACRES BEARING")

subsetNASS= function(NASS,croptypes,production_types){
  
  dt1= subset(NASS,is.element(NASS$COMMODITY_DESC,croptypes))
  dt2= subset(dt1,is.element(dt1$SHORT_DESC,production_types))
  
  return(dt2)
}

culti= subsetNASS(NASS,croptypes,production_types)
national=culti[culti$AGG_LEVEL_DESC=="NATIONAL",]
state=culti[culti$STATE_NAME=="CALIFORNIA" & culti$AGG_LEVEL_DESC=="STATE",]
county=culti[culti$STATE_NAME=="CALIFORNIA" & culti$COUNTY_NAME=="NAPA" & culti$AGG_LEVEL_DESC=="COUNTY",]
   
sum(state$VALUE,na.rm = T)/sum(national$VALUE,na.rm = T)
sum(county$VALUE,na.rm = T)/sum(state$VALUE,na.rm = T)


