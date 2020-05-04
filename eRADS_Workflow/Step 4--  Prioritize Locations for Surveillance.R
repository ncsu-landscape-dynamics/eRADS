########################## 4. Prioritize Locations for Surveillance #############################
## Detection probability ##
p_dtc = function(surv_den,p_den,sensi){
  p_dtc=1-(1-sensi)^(surv_den*p_den)
  return(p_dtc)
}

## Surveillance density ##
surv_den = function(confi,p_den,sensi){
  #surv_den=log(1-confi, base = exp(1))/log(1-p_den*sens, base = exp(1))
  p_den[p_den==0]=NA
  surv_den=log(1-confi)/(p_den*log(1-sensi))
  
  return(surv_den)
}

host2=ra2plyMerge(host)
host2bf=buffer(host2,width=2000)
host2bf=gUnionCascaded(host2bf)

#### Method 1 -- Select all locations with same detection confidence ####
den=surv_den(0.90,pden,0.1)
sum(getValues(den),na.rm = T)/sum(getValues(pden/pden),na.rm=T)
plot(den,main="Trap Density")
plot(cty,add=T,border="blue")
plot(host2bf,add=T,border="red")
den2=crop(den,host2bf)
den2=mask(den2,host2bf)
sum(getValues(den2),na.rm = T)*77

#### Method 2 -- Select desired % of locations with highest infestation ####
percent_locat = function(p_den,percentage,confi,sensi){
  
  n_all=sum(getValues(p_den/p_den),na.rm=T)
  
  pts=rasterToPoints(p_den,fun=function(x){x>0},spatial = T)
  
  pts$avl=extract(p_den,pts)
 
  seleN=round(n_all*percentage)
  
  pts2=pts[order(pts$avl,decreasing = T),]
  
  pts3=pts2[1:seleN,]
  ra2=rasterize(pts3,p_den,field=1,background=0)
  
  p_selec=ra2*p_den
  survDen=surv_den(confi,p_selec,sensi)
  
  survTT=sum(getValues(survDen),na.rm = T)
  n=sum(getValues(p_den/p_den),na.rm=T)
  survMden=survTT/n
  
  results=list(survDen,survTT,survMden)
  
  return(results)
}

results=percent_locat(pden,0.95,0.95,0.1)
plot(results[[1]])
plot(cty,add=T,border="blue")
results[[2]]
results[[3]]
den=results[[1]]
den2=crop(den,host2bf)
den2=mask(den2,host2bf)
sum(getValues(den2),na.rm = T)*77


#### Method 3 -- Select desired % of populations with highest infested locations prioritized ####
percent_pop = function(p_den,percentage,confi,sensi){
  
  nt=sum(getValues(p_den/p_den),na.rm=T)
  
  pts=rasterToPoints(p_den,fun=function(x){x>0},spatial = T)
  
  pts$avl=extract(p_den,pts)
  ap=sum(getValues(p_den), na.rm=T)
  ctp1=ap*percentage
  
  pts2=pts[order(pts$avl,decreasing = T),]
  
  pts2$V3=sapply(1:length(pts2),FUN=function(i)sum(pts2$avl[1:i]))
  pts3=pts2[pts2$V3<=ctp1,]
  
  ra2=rasterize(pts3,p_den,field=1,background=0)
  
  p_selec=ra2*p_den
  survDen=surv_den(confi,p_selec,sensi)
  pixel_percent=length(pts3)/nt
  
  survTT=sum(getValues(survDen),na.rm = T)
  survMden=survTT/nt
  
  results=list(p_selec,survDen,pixel_percent,survTT,survMden)
  return(results)
}

results=percent_pop(pden,0.95,0.95,0.1)
survD=results[[2]]
results[[3]]
results[[4]]
results[[5]]
plot(survD)
plot(cty,add=T,border="blue")

den=survD
den2=crop(den,host2bf)
den2=mask(den2,host2bf)
sum(getValues(den2),na.rm = T)*77


#### Method 4 -- Select locations with infestation higher than given threshold ####
pop_thresh = function(p_den=pden,threshold,confi,sensi){
  p_den2=reclassify(p_den,c(0,threshold,NA))
  survden= surv_den(confi,p_den2,sensi)
  
  n1=sum(getValues(p_den/p_den),na.rm = T)
  n2=sum(getValues(p_den2/p_den2),na.rm = T)
  ratio=n2/n1
  
  survTT=sum(getValues(survden),na.rm = T)
  survMden=survTT/n1
  
  results=list(ratio,survden,survTT,survMden)
  return(results)
}

results2=pop_thresh(pden,0.5,0.95,0.1)
results2[[1]]
survD=results2[[2]]
results2[[3]]
results2[[4]]
plot(survD)
plot(cty,add=T,border="blue")

den=survD
den2=crop(den,host2bf)
den2=mask(den2,host2bf)
sum(getValues(den2),na.rm = T)*77

#### Method 5 -- Prioritize locations based on (probability), (epidemic size), and (hazard),#### 
##and select % locations ##
no=list(as.matrix(inf*0))
inf_Itensi=stack()
inf_Risk=stack()

for (i in 1:10){
  
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
risk=mean(inf_Risk)

rank_selectLoca = function(inf,inf_more,Wcoef=Tcoef,dispersal_rate=360,reprodution_rate=20,Nstep=8,host,risk,percentage=p){
  
  ply0=ra2plyMerge(inf)
  ply1=ra2plyMerge(inf_more)
  
  gdiff=gDifference(ply1,ply0)
  infnew=crop(inf_more,gdiff)
  infnew=mask(infnew,gdiff)
  
  ts=infnew/infnew
  nt=sum(getValues(ts),na.rm=T)
  
  
  new_ply=ra2plyMerge(infnew)
  
  rt=host_m6(infnew,dispersal_rate,reproduction_rate,Nstep,Wcoef,ply1,host)
  
  
  pixelArea=xres(inf)*yres(inf)
  prob1=extract(risk,rt,buffer=0,fun=sum,na.rm=T)
  ply_Infsize=extract(host,rt,buffer=0,fun=sum,na.rm=T)
  
  
  surv_rank1=prob1*ply_Infsize
  surv_rank2=prob1*rt$Nhost
  surv_rank3=prob1*(ply_Infsize+rt$Nhost)
  
  
  rt$rank_Epidsize=surv_rank1
  rt$rank_Hzd=surv_rank2
  rt$rank_EpidsizeHzd=surv_rank3
  
  rank_Epidsize=rt[order(rt$rank_Epidsize,decreasing = T),]
  rank_Hzd=rt[order(rt$rank_Hzd,decreasing = T),]
  rank_EpidsizeHzd=rt[order(rt$rank_EpidsizeHzd,decreasing = T),]
  
  nselec=round(nt*p)
  
  rank_Epidsize=rank_Epidsize[1:nselec,]
  rank_Hzd=rank_Hzd[1:nselec,]
  rank_EpidsizeHzd=rank_EpidsizeHzd[1:nselec,]
  
  ls=list(rank_Epidsize,rank_Hzd,rank_EpidsizeHzd)
  
  return(ls)
}

ranks= rank_selectLoca(inf,inf_more,Wcoef = Tcoef,dispersal_rate = 360,reprodution_rate = 20,Nstep = 8,host,risk=risk)
epidsize=ranks[[1]]
hzd=ranks[[2]]
epidhzd=ranks[[3]]

ply=ra2plyMerge(inf)
ra_ply=rasterToPolygons(inf,fun=function(x) x>0)
infestation_potential=Thost(ra_ply,1000,ra_ply[1,],host)
infestation_potential$epidsize=infestation_potential$kde*infestation_potential$kde*(1+infestation_potential$Thst)
infestation_potential=infestation_potential[order(infestation_potential$epidsize,decreasing = T),]
plot(inf)




rank_selectPopu = function(inf,inf_more,Wcoef=Tcoef,dispersal_rate=360,reprodution_rate=20,Nstep=8,host,risk,percentage=p){
  
  ply0=ra2plyMerge(inf)
  ply1=ra2plyMerge(inf_more)
  
  gdiff=gDifference(ply1,ply0)
  infnew=crop(inf_more,gdiff)
  infnew=mask(infnew,gdiff)
  
  ts=infnew/infnew
  nt=sum(getValues(ts),na.rm=T)
  
  
  new_ply=ra2plyMerge(infnew)
  
  rt=host_m6(infnew,dispersal_rate,reproduction_rate,Nstep,Wcoef,ply1,host)
  
  
  pixelArea=xres(inf)*yres(inf)
  prob1=extract(risk,rt,buffer=0,fun=sum,na.rm=T)
  ply_Infsize=extract(host,rt,buffer=0,fun=sum,na.rm=T)
  
  
  surv_rank1=prob1*ply_Infsize
  surv_rank2=prob1*rt$Nhost
  surv_rank3=prob1*(ply_Infsize+rt$Nhost)
  
  
  rt$rank_Epidsize=surv_rank1
  rt$rank_Hzd=surv_rank2
  rt$rank_EpidsizeHzd=surv_rank3
  
  rt$infIntensity=unlist(extract(infnew,rt,na.rm=T))
  
  rank_Epidsize=rt[order(rt$rank_Epidsize,decreasing = T),]
  rank_Hzd=rt[order(rt$rank_Hzd,decreasing = T),]
  rank_EpidsizeHzd=rt[order(rt$rank_EpidsizeHzd,decreasing = T),]
  
  Popselec=sum(rt$infIntensity)*p
  
  avl2$V3=sapply(1:dim(avl2)[1],FUN=function(i)sum(avl2[1:i,1]))
  
  rank_Epidsize$Cumu_inf=sapply(1:length(rt),FUN= function(i)sum(rank_Epidsize$infIntensity[1:i]))
  rank_Hzd$Cumu_inf=sapply(1:length(rt),FUN= function(i)sum(rank_Hzd$infIntensity[1:i]))
  rank_EpidsizeHzd$Cumu_inf=sapply(1:length(rt),FUN= function(i)sum(rank_EpidsizeHzd$infIntensity[1:i]))
  
  rank_Epidsize=rank_Epidsize[rank_Epidsize$Cumu_inf<=Popselec,]
  rank_Hzd=rank_Hzd[rank_Hzd$Cumu_inf<=Popselec,]
  rank_EpidsizeHzd=rank_EpidsizeHzd[rank_EpidsizeHzd$Cumu_inf<=Popselec,]
  
  ls=list(rank_Epidsize,rank_Hzd,rank_EpidsizeHzd)
  
  return(ls)
}

ranks= rank_selectPopu(inf,inf_more,Wcoef = Tcoef,dispersal_rate = 360,reprodution_rate = 20,Nstep = 8,host,risk=risk)
epidsize=ranks[[1]]
hzd=ranks[[2]]
epidhzd=ranks[[3]]



 ######################################## 5. Budget Limit & Location Selection #######################################
#### Method 1 -- Select all locations with same detection confidence ####


den=surv_den(0.80,pden,0.1)
sum(getValues(den),na.rm = T)/sum(getValues(pden/pden),na.rm=T)
plot(den,main="Trap Density")
plot(cty,add=T,border="blue")
plot(host2bf,add=T,border="red")
den2=crop(den,host2bf)
den2=mask(den2,host2bf)
sum(getValues(den2),na.rm = T)*77


#### Method 2 & 3  Optimal Detection Confidence based on Locations or Populations ####
confi_all = function(confi,p_den,sensi,budget,cost_per_survey){
  
  n=budget/cost_per_survey
  survDen=surv_den(confi,p_den,sensi)
  avl=getValues(survDen)
  avl=as.data.frame(avl)
  avl$ID=1:dim(avl)[1]
  avl$Pden=getValues(p_den)
  avl2=avl[order(avl$avl,decreasing = F,na.last = NA),]
  avl2=avl2[avl2$avl>0,]
  avl2$cumu_number=sapply(1:dim(avl2)[1],function(i) sum(avl2$avl[1:i]))
  
  avl3=avl2[avl2$cumu_number<= n & avl2$cumu_number>0,]
  
  n_sel=dim(avl3)[1]
  n_all=sum(getValues(p_den/p_den),na.rm = T)
  
  conf_loca= (n_sel/n_all)*confi
  conf_popu= (sum(avl3$Pden)/sum(getValues(p_den),na.rm = T))*confi
  perc_loca= (n_sel/n_all)
  
  confall=list(conf_loca,conf_popu,perc_loca)
  return(confall)
}

optiConfi = function(p_den,sensi,budget,cost_per_survey,confiMin=0.10,confiMax=0.99){
  require(parallel)
  require(foreach)
  require(doParallel)
  
  cl=makeCluster(detectCores()-6)
  registerDoParallel(cl)
  
  confi=seq(confiMin,confiMax,0.01)
  n=length(confi)
  
  confAll= foreach(i=1:n, .combine = rbind,.export = c("confi_all","sensi","surv_den"),.packages = c("raster","rgdal","rgeos")) %dopar% {
    confall=confi_all(confi[i],p_den,sensi,budget,cost_per_survey)
    result=as.data.frame(matrix(c(confall[[1]],confall[[2]]),nrow=1,ncol=2))
    return(result)
  }
  
  stopCluster(cl) 
  results=as.data.frame(confAll) 
  colnames(results)=c("Location","Population")
  results$confi=seq(confiMin,confiMax,0.01)
  optiConfi_location=results[which.max(results$Location),]$confi
  optiConfi_population=results[which.max(results$Population),]$confi
  
  ls=list(results,optiConfi_location,optiConfi_population)
  return(ls)
}

opti=optiConfi(p_den,sensi=0.1,budget=2000000,cost_per_survey=77,confiMin=0.30,confiMax=0.99)
optiV=opti[[1]]
opti[[2]]
opti[[3]]
confi=seq(0.30,0.99,0.01)
opti$confi=confi
par(mfrow=c(1,1))

plot(optiV$Location~optiV$confi,font.lab=2,type="b",pch=2,ylab="Detection Confidence of All Locations",col="blue",xlab="Detection Confidence for Each Location",main="")
plot(optiV$Population~confi,xlab="Detection Confidence for Each Location",ylab="Detection Confidence of All Populations",font.lab=2,type="b",pch=2,col="blue")


results=percent_locat(pden,0.7166287,0.97,0.1)
plot(results[[1]])
plot(cty2,add=T,border="blue")
results[[2]]
results[[3]]
den=results[[1]]
den2=crop(den,host2bf)
den2=mask(den2,host2bf)
sum(getValues(den2),na.rm = T)*77

#### Method 4 selection locations with intensity higher than a threshold ###
results2=pop_thresh(pden,0.5,0.95,0.1)
results2[[1]]
surv_den=results2[[2]]
surv_den2=surv_den/(sum(getValues(surv_den),na.rm = T)/25974)
pdetc=p_dtc(surv_den2,pden,0.1)
plot(surv_den2)
plot(cty2,add=T,border="blue")

#### Regular surveillance strategy ####
reg_surv=function(p_den,sensi,budget,cost_per_survey){
  
  n=sum(getValues(p_den/p_den),na.rm = T)
  nSurv=budget/cost_per_survey
  
  survDen=nSurv/n
  conf=p_dtc(survDen,p_den,sensi)
  meanconf=sum(getValues(conf),na.rm=T)/n
  results=list(survDen,conf,meanconf)
}

regSurv=reg_surv(p_den,sensi,budget,cost_per_survey)
survDen=regSurv[[1]]
confi=regSurv[[2]]
meanconfi=regSurv[[3]]
survDen
plot(confi)
plot(cty2,add=T,border="blue")
summary(confi)





#### Method 5 ####
infestation_potential

confi_all = function(confi,p_den,sensi,budget,cost_per_survey,infestation_potential){
  
  n=budget/cost_per_survey
  survDen=surv_den(confi,p_den,sensi)
  avl=extract(survDen,infestation_potential)
  avl=as.data.frame(unlist(avl))
  avl[is.na(avl)]=0
  avl2=avl
  colnames(avl2)="avl"
  avl2$cumu_number=sapply(1:dim(avl2)[1],function(i) sum(avl2$avl[1:i]))
  
  avl3=avl2[avl2$cumu_number<= n & avl2$cumu_number>0,]
  
  n_sel=dim(avl3)[1]
  n_all=sum(getValues(p_den/p_den),na.rm = T)
  
  conf_loca= (n_sel/n_all)*confi
  perc_loca= (n_sel/n_all)
  
  confall=list(conf_loca,perc_loca)
  return(confall)
}

optiConfi = function(p_den,sensi,budget,cost_per_survey,confiMin=0.10,confiMax=0.99,infestation_potential){
  require(parallel)
  require(foreach)
  require(doParallel)
  
  cl=makeCluster(detectCores()-6)
  registerDoParallel(cl)
  
  confi=seq(confiMin,confiMax,0.01)
  n=length(confi)
  
  confAll= foreach(i=1:n, .combine = rbind,.export = c("confi_all","sensi","surv_den","infestation_potential"),.packages = c("raster","rgdal","rgeos")) %dopar% {
    confall=confi_all(confi[i],p_den,sensi,budget,cost_per_survey,infestation_potential)
    result=as.data.frame(matrix(c(confall[[1]],confall[[2]]),nrow=1,ncol=2))
    return(result)
  }
  
  stopCluster(cl) 
  results=as.data.frame(confAll) 
  colnames(results)=c("Location_Confi","Location_Percent")
  results$confi=seq(confiMin,confiMax,0.01)
  optiConfi_location=results[which.max(results$Location_Confi),]$confi
  
  ls=list(results,optiConfi_location)
  return(ls)
}

opti=optiConfi(p_den,sensi=0.1,budget=2000000,cost_per_survey=77,confiMin=0.30,confiMax=0.99,infestation_potential)
optiV=opti[[1]]
opti[[2]]
opti[[3]]
confi=seq(0.30,0.99,0.01)
opti$confi=confi
par(mfrow=c(1,1))
plot(optiV$Location_Confi~optiV$confi,xlab="Detection Confidence for Each Location",ylab="Detection Confidence of All Populations",font.lab=2,type="b",pch=2,col="blue")

survD=surv_den(0.78,pden,0.1)
n2=round(length(infestation_potential)*0.9157175)
infestation_potential2=infestation_potential[1:n2,]
survD2=crop(survD,infestation_potential2)
survD2=mask(survD2,infestation_potential2)
plot(survD2)
plot(cty2,add=T,border="blue")
sum(getValues(survD2),na.rm=T)
