######################### 6. Eradication Feasibility ######################
multi_year_pops=function(number_of_year,inf,annual_budget,cost_per_meter_sq){
  
  ny=number_of_year
  inf_number=sum(getValues(inf),na.rm=T)
  
  
  area=data.frame(matrix(0,nrow=ny,ncol=1))
  number=data.frame(matrix(0,nrow=ny,ncol=1))
  
  for (j in 1:ny){
    
    inf=inf
    
    infected_speci[[1]]=as.matrix(inf)
    susceptible_speci[[1]]=as.matrix(host - inf)
    
    susceptible_species <- list(as.matrix(host - inf))
    
    treat_ment=treat_m6(inf,dispersal_rate=360,reproductive_rate=20,Nstep=9,Wcoef,ply,host, budget, buffer,cost_per_meter_sq)
    #treat_ment=list(as.matrix(inf)*0)
    ls=list()
    
    for (i in 1:20){
      random_seed=i*10808+i
      
      #treatment_date=treatmens[i]
      
      data <- pops_model(random_seed = random_seed, 
                         use_lethal_temperature = use_lethal_temperature, 
                         lethal_temperature = lethal_temperature, 
                         lethal_temperature_month = lethal_temperature_month,
                         
                         infected = infected_speci[[1]],
                         susceptible = susceptible_speci[[1]],
                         total_plants = total_pl[[1]],
                         mortality_on = mortality_on,
                         mortality_tracker = infected_speci[[1]]*0,
                         mortality = infected_speci[[1]]*0,
                         treatment_maps = treat_ment,
                         
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
      
      
      infect=as.matrix(as.data.frame(data$infected[1]))
      ls[[i]]=infect
    }
    
    inf_matrix=apply(simplify2array(ls),1:2,median)
    inf=raster(inf_matrix,template=inf)
    
    area[j,1]=sum(getValues(inf/inf),na.rm=T)
    number[j,1]=sum(getValues(inf),na.rm=T)
    
    
  }
  
  
  list=list(area,number,inf)
  return(list)
  
}



## Total budget ##
erad_possiT= function(inf,total_budget, cost_per_meter_sq){
  
  pixelArea=xres(inf)*yres(inf)
  ttInf=sum(getValues(inf/inf),na.rm=T)
  
  erad_Area=pixelArea*ttInf
  Ttreat_Area= total_budget/cost_per_meter_sq
  
  minBu=erad_Area*cost_per_meter_sq/1000
  
  percent=Ttreat_Area*100/erad_Area
  percent=round(percent,digits=2)
  if (percent>=100){
    cat(paste("You have budget to treat ",as.character(percent),"% of the total infested area.",sep=""), "\n")
    cat("There is a possibity for eradication!!!")
  } else {
    cat(paste("You have budget to treat ",as.character(percent),"% of the total infested area.",sep=""), "\n")
    cat(paste("More budget is needed. The minimum budget for eradication is ",as.character(minBu),"K dollar.",sep=""))}
}
erad_possiT(inf,5000000,1.37)

surd=surv_den(0.999,inf,0.1)
surv_cost=sum(getValues(surd),na.rm=T)*77
treat_cost=sum(getValues(surd/surd),na.rm=T)*300/0.004
surv_cost+treat_cost

#### Based on annual budget ####

## simulate management result for 3 years, if keep increasing, eradication not feasibile 
## simulate management result for 3 years, if keep decreasing, keep simulating till eradication
erad_possiA = function(inf, annual_budget){
  results=multi_year_pops(3,inf,150000,1.37)
  
  inf_area=results[[1]]
  diff=inf_area[-1,1]-inf_area[-dim(inf_area)[1],1]
  
  if (names(which.max(table(diff>0)))==T){cat("The annual budget is not enough for eradication. The total infected area is increasing with current budget.")}
  else {cat("Great! The total infected area is decreasing with current annual budget.")}
}
results=erad_possiA(inf,annual_budget=500000)

results=multi_year_pops(3,inf,150000,1.37)
number2=results[[1]]
number=results[[1]]
plot(number[1:3,1]~as.integer(c(1:3)),font.lab=2,col="blue",xlab="Invasion Year",ylab="Number of Infested Cell",main="Number of Infested Cell after Treatment")

trt=c(636,698,888,1395)
notrt=c(636,948,2098,5838)
year=c(0,1,2,3)
plot(notrt~year,font.lab=2,col="red",ylim=c(600,6200),xlab="Invasion Year",ylab="Number of Infested Cell",main="Number of Infested Cell")
points(year,trt,col="blue")

l1=lm(trt~year+I(year^2))
lno=lm(notrt~year+I(year^2))

x=as.data.frame(seq(0,3,0.01))
colnames(x)="year"
y1=predict(l1,x)
y2=predict(lno,x)

lines(x$year,y1,col="blue")
lines(x$year,y2,col="red")
legend(x=0.01,y=5500,col=c("red","blue"),pch=1,lty=1,legend=c("No Treatment","Optimal Treatment"))

# minimum first year budget to make the infeastion have a decreasing trend
pops= function(trt,inf){
  
  trt=gUnionCascaded(trt)
  toRa=rasterize(trt,inf,field=1,background=0)
  
  trt=list(as.matrix(toRa))
  df=data.frame(matrix(0,nrow=20,ncol=1))
  
  for (j in 1:20){
    
    random_seed= j*1000+50
    
    
    data <- pops_model(random_seed = random_seed, 
                       use_lethal_temperature = use_lethal_temperature, 
                       lethal_temperature = lethal_temperature, 
                       lethal_temperature_month = lethal_temperature_month,
                       
                       infected = infected_speci[[1]],
                       susceptible = susceptible_speci[[1]],
                       total_plants = total_pl[[1]],
                       mortality_on = mortality_on,
                       mortality_tracker = infected_speci[[1]]*0,
                       mortality = infected_speci[[1]]*0,
                       treatment_maps = trt,
                       
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
    
    area=data$area_infected
    area=area/10000
    
    df[j,1]=area
  }
  
  mean=mean(df[1:20,1])
  return(mean)
}

treatment_ranks= host_m6(inf,dispersal_rate=360,reproductive_rate=20,Nstep=9,Wcoef,ply,host)
treatment_ranks=treatment_ranks[order(treatment_ranks$Nhost,decreasing = T),]

reproductive_rate[[2]]=20
minBudget = function(rank,buffer,cost_per_meter_sq){
  
  treatment_ranks=rank
  treatment_buffer=gBuffer(treatment_ranks,byid=T,width = buffer)
  
  tn=length(rank)
  #perc=round(c(tn*0.1,tn*0.2,tn*0.3,tn*0.4,tn*0.5,tn*0.6,tn*0.7,tn*0.8,tn*0.9))
  perc=round(c(tn*0.2,tn*0.4,tn*0.6,tn*0.8))
  
  areas=sapply(perc, function(x) pops(treatment_buffer[1:x,],inf))
  results=areas<(tn*0.95)
  #results=areas<(tn)
  
  decrease=results[which(results==T)]
  ntrue=length(decrease)
  seqs=seq(0,0.8,0.2)
  
  if (ntrue==0){newtrys=round(c(tn*0.85,tn*0.9,tn*0.95))}
  else {
    newseqs=seq(seqs[5-ntrue],seqs[6-ntrue],0.05)
    newtrys=round(tn*newseqs)
    newtrys=newtrys[-1]
  }
  
  areas2=sapply(newtrys, function(x) pops(treatment_buffer[1:x,],inf))
  results2=areas2<(tn*0.95)
  #results2=areas2<(tn)
  decrease2=newtrys[results2==T]
  
  minTreat=treatment_buffer[1:decrease2[1],]
  minTreat=gUnionCascaded(minTreat)
  minBudget=area(minTreat)*cost_per_meter_sq
  minBudget2=round(minBudget/1000)
  #cat(paste("The minimum first year budget to decrease infestation is ",as.character(minBudget),"k dollar.",sep=""))
  print.noquote(paste("The minimum budget to decrease infestation is ",as.character(minBudget2),"k dollar.",sep=""))
  return(minBudget)
}

minBudget(treatment_ranks,buffer=0,cost_per_meter_sq = 1.37)
minBudget(treatment_ranks,buffer=50,cost_per_meter_sq = 1.37)
minBudget(treatment_ranks,buffer=150,cost_per_meter_sq = 1.37)

