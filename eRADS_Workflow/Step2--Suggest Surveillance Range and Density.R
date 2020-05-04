################ 1  -- decide surveillance range to estimate infestation range ###################

setwd("Q:/Shared drives/APHIS  Projects/eRADS/Biweekly meetings/EGVM/")
library(rgdal)
library(raster)
library(rgeos)
library(outliers)
library(rgdal)

#### Step 1 -- decide range size -- ####
napa=readOGR("Q:/Shared drives/APHIS  Projects/eRADS/caseStudy","Napa")
pre=readOGR("Q:/Shared drives/APHIS  Projects/eRADS/caseStudy/EGVM/.","Y09All")

napa=spTransform(napa,crs(pre))
plot(napa)
plot(pre,add=T,col="red")

dt=as.data.frame(pre@coords)
cl=kmeans(dt,centers=2,nstart=2)
pre2=pre[cl$cluster==modal(cl$cluster),]

plot(pre2,add=T,col="blue")

ctd=gCentroid(pre2)
plot(pre2,pch=4,cex=1.5)
plot(ctd,pch=18,add=T,cex=2,col="red")
dis=gDistance(ctd,pre,byid=T)
summary(dis)

pre3=gDifference(pre,pre2)
plot(pre3,add=T,col="orange")
ctd2=gCentroid(pre3)

dis2=gDistance(ctd,ctd2)
dis2

spreadRate=20000

if (length(pre3)>2){
  buffer=dis2+spreadRate
} else {buffer=dis2}

setwd("Q:/Shared drives/APHIS  Projects/eRADS/caseStudy/")
survMin=gBuffer(ctd,byid=F,width=dis2)
surv=gBuffer(ctd,byid=F,width=buffer)
cty=readOGR("EGVM//CA_counties/.","CA_counties")

#### Step 2 -- load host and suggest locations for surveillance ####
host=raster("Q:/Shared drives/APHIS  Projects/eRADS/caseStudy/host_grape.tif")
surv2=spTransform(surv,crs(host))
survMin2=spTransform(survMin,crs(host))
cty2=spTransform(cty,crs(host))
pre_tr=spTransform(pre,crs(host))

host2=crop(host,surv2)
host2=mask(host2,surv2)
host_ply=rasterToPolygons(host2)


plot(surv2,border="orange")
plot(survMin2,add=T,border="red")
plot(cty2,add=T)  
plot(host_ply,add=T,col="purple",border="purple")
plot(pre_tr,add=T,col="blue")



#### Step 3 -- suggest surveillance density ####
# Surveillance density #
surv_den = function(confi,p_den,sensi){
  #surv_den=log(1-confi, base = exp(1))/log(1-p_den*sens, base = exp(1))
  p_den2=reclassify(p_den,c(0,Inf,1))
  area=sum(getValues(p_den2),na.rm = T)
  
  surv_den=log(1-confi)/log(1-p_den*sensi)
  surv_den=reclassify(surv_den,c(-Inf,0,0),include.lowest=T)
  return(surv_den)
}

inf=host2*2/host2
den=surv_den(0.99,inf,0.1)
plot(den)


