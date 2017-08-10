#----------------------------------------------------------------------------------------------------------------
#- This is a collection of many functions that do the actual work of data analysis and plotting. These functions
#    are called by just a few lines of code in "main script.R" to recreate the analyses and figures.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------------
standard.error <- function(dat,na.rm=F,...){
  if(na.rm==T){
    dat <- subset(dat,is.na(dat)==F)
  }
  std <- sd(dat)
  n <- length(dat)
  se <- std/sqrt(n)
  return(se)
}
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
#- function to partition the hourly net flux observations into GPP and Ra using an Arrhenious function
partitionHourlyFluxCUE_arr <- function(data.hr.gf=data.hr.gf,Ea=57.69,lagdates,leafRtoTotal = 1,leafRreduction=0){
  rvalue = 8.134
  require(data.table)
  
  #-- convert mmolCO2 s-1 to gC hr-1
  data.hr.gf$FluxCO2_g <- with(data.hr.gf,FluxCO2*60*60/1000*12.0107)
  data.hr.gf$period <- ifelse(data.hr.gf$PAR>2,"Day","Night")
  
  #-- partition day-time net C exchange into GPP and Ra, similar to how it is done in eddy-covariance.
  #-- create a series of dates
  date.vec <- seq.Date(from=min(data.hr.gf$Date),to=max(data.hr.gf$Date),by="day")
  
  #-- estimate R-Tref and Tref for each date for each chamber
  #lagDates <- 3 # establish how many prior days to include
  RTdat <- expand.grid(Date=date.vec,chamber=levels(data.hr.gf$chamber))
  RTdat$Tref <- RTdat$R_Tref <- NA
  
  
  print("Partitioning Net CO2 fluxes into GPP and Ra")
  #- set up progress bar to track that this is working
  pb <- txtProgressBar(min = 0, max = nrow(RTdat), style = 3)
  
  for (i in 1:nrow(RTdat)){
    #- trial a data.table alternative to speed this up. The filter thing was actually slower.
    
    #dat <- dplyr::filter(data.hr.gf,chamber==RTdat$chamber[i],Date <= RTdat$Date[i], Date >= (RTdat$Date[i]-lagDates),period =="Night")
    #RTdat$Tref[i] <- mean(dat$Tair_al,na.rm=T)
    #RTdat$R_Tref[i] <- mean(dat$FluxCO2_g,na.rm=T)
    
    inds <- which(data.hr.gf$chamber==RTdat$chamber[i] & data.hr.gf$Date <= RTdat$Date[i] & data.hr.gf$Date >= (RTdat$Date[i]-lagdates) & data.hr.gf$period =="Night" )
    RTdat$Tref[i] <- mean(data.hr.gf$Tair_al[inds],na.rm=T)
    RTdat$R_Tref[i] <- mean(data.hr.gf$FluxCO2_g[inds],na.rm=T)
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  
  RTdat$Tref_K <- with(RTdat,Tref+273.15)
  
  #-- merge these reference data into the gap-filled flux dataframe to estimate Ra during the daytime, and hence GPP
  data.hr.gf3 <- merge(data.hr.gf,RTdat,by=c("Date","chamber"))
  #data.hr.gf3$Ra_est <- with(data.hr.gf3,R_Tref*Q10^((Tair_al-Tref)/10)) # estimate respiration rate. This is a negative number.
  data.hr.gf3$Ra_est <- with(data.hr.gf3,R_Tref*exp((Ea*1000/(rvalue*Tref_K))*(1-Tref_K/(Tair_al+273.15)))) # estimate respiration rate. This is a negative number.
  data.hr.gf3$Ra_est <- ifelse(data.hr.gf3$period=="Day",
                              data.hr.gf3$Ra_est-leafRreduction*leafRtoTotal*data.hr.gf3$Ra_est,
                              data.hr.gf3$Ra_est) # estimate respiration rate. This is a negative number. If it's day, subtract 30% from the leaf R fraction
  

  data.hr.gf3$GPP <- ifelse(data.hr.gf3$period=="Night",0,data.hr.gf3$FluxCO2_g-data.hr.gf3$Ra_est)
  data.hr.gf3$Ra <- ifelse(data.hr.gf3$period=="Night",data.hr.gf3$FluxCO2_g,data.hr.gf3$Ra_est)
  
  return(data.hr.gf3)
}
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#- function to return daily sums of GPP and Ra from hourly data
returnCUE.day <- function(dat=data.hr.p){
  
  # #- calculate leaf-area specific rates of GPP, convert to umol CO2 m-2 s-1
  dat$GPP_la <- with(dat,GPP/leafArea)
  dat$GPP_la_umol <- with(dat,GPP_la/12*1*10^6/60/60)
  dat$PAR_mol <- dat$PAR*60*60*1*10^-6

  #- create a date variable that moves the window for "night" observations, to sum observations for the night period following a day period
  dat$Date2 <- as.Date(dat$DateTime)
  dat$hour <- hour(dat$DateTime)
  earlynights <- which(dat$hour < 12 & dat$period =="Night")
  dat$Date2[earlynights] <- dat$Date2[earlynights]-1
  
  #- Identify and name the drought/watered treatments
  drought.chamb = unique(dat$chamber[ dat$Water_treatment %in% as.factor("drydown")])
  dat$chamber_type = as.factor( ifelse(dat$chamber %in% drought.chamb, "drought", "watered") )
  
  #- create daily sums 
  data.day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment,chamber_type,period),
                              GPP = sum(GPP,na.rm=T),
                              Ra = sum(Ra,na.rm=T),
                              FluxCO2_g=sum(FluxCO2_g,na.rm=T),
                              Tair_al=mean(Tair_al,na.rm=T),
                              VPDair=max(VPDair,na.rm=T),
                              PAR=sum(PAR_mol,na.rm=T),
                              leafArea=mean(leafArea,na.rm=T))
  data.day <- as.data.frame(data.day)
  data.day <- data.day[with(data.day,order(Date2,chamber)),]
  names(data.day)[1] <- "Date"
  
  
  Tair_day <- dplyr::summarize(group_by(dat,Date2,chamber,T_treatment,chamber_type),
                               Tair_al=mean(Tair_al,na.rm=T))
  names(Tair_day)[5] <- "Tair_24hrs"                   
  
  data.day_tair <- as.data.frame(Tair_day)
  data.day <- as.data.frame(data.day)
  
  #- merge night and day data
  data.day2d <- subset(data.day,period=="Day")[,c("Date","chamber","T_treatment","GPP","Ra","FluxCO2_g","PAR","Tair_al","VPDair","leafArea","chamber_type")]
  names(data.day2d)[4:8] <- c("GPP","Raday","Cgain","PAR","T_day")
  data.day2n <- subset(data.day,period=='Night')[,c("Date","chamber","T_treatment","Ra","FluxCO2_g","Tair_al","VPDair","chamber_type")]
  names(data.day2n)[4:7] <- c("Ranight","Closs","T_night","VPDair_night")
  names(data.day2n)[1] <- c("Date")
  
  #- merge data such that the nightly data are one day ahead of the daily data (i.e., does tonight's respiration depend on today's photosynthesis?)
  cue.day1 <- subset(merge(data.day2d,data.day2n,by=c("Date","chamber","T_treatment","chamber_type")))
  cue.day <- merge(cue.day1,data.day_tair,by.x=c("Date","chamber","T_treatment","chamber_type"),by.y=c("Date2","chamber","T_treatment","chamber_type"))
  cue.day$Ra <- with(cue.day,-1*(Raday+Ranight))
  cue.day$RtoA <- with(cue.day,Ra/GPP)
  cue.day$CUE <- with(cue.day,(1-(Ra/GPP)))
  cue.day$GPP_la <- with(cue.day,GPP/leafArea)
  cue.day$Ra_la <- with(cue.day,Ra/(leafArea)) #g C m-2 day-1
  
  # remove a few days from the beginning of the dataset in C09. Problem with flux data.
  toremove <- which(cue.day$Ra_la < 1.5 & cue.day$chamber=="C07" & cue.day$Date < as.Date("2013-10-3"))
  cue.day <- cue.day[-toremove,] 
  
  # #- merge in the key
  # key <- read.csv("data/WTC_TEMP_CM_TREATKEY_20121212-20140528_L1_v1.csv")
  # key$chamber_type = as.factor( ifelse(key$chamber %in% drought.chamb, "drought", "watered") )
  # 
  # cue.day$T_treatment <- ifelse(cue.day$T_treatment=="ambient","ambient","elevated")
  # cue.day2 <- merge(cue.day,key,by=c("chamber","Date","T_treatment","chamber_type"))
  
  # cue.day.trt <- summaryBy(GPP+GPP_la+Ra+Ra_la+RtoA+PAR+T_day+T_night+Tair_24hrs+VPDair+VPDair_night~Date+T_treatment+Water_treatment,data=cue.day2,FUN=c(mean,standard.error))
  # cue.day.trt <- summaryBy(GPP+Ra+leafArea~Date+T_treatment+chamber_type,data=cue.day2,FUN=c(mean,standard.error,length))
  # return(list(cue.day2,cue.day.trt))
  cue.day.trt <- summaryBy(GPP+Ra+leafArea~Date+T_treatment+chamber_type,data=cue.day,FUN=c(mean,standard.error))
  return(list(cue.day,cue.day.trt))
}
#----------------------------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
#-- Code to return tree mass for every day of the experiment.
#- this code from WTC_Rallometry_script.R in the WTC3 folder. I'm trying to estimate tree mass at the beginning of the experiment
returnTreeMass <- function(plotson=F){
  
  #- get an estimate of d2h for every day of the experiment by interpolation
  size <- returnd2h()
  
  #- get an estimate of volume and mass for wood and bark for each measurement day
  print("Calculating stem volume")
  vol <- getvol()
  
  #- interpolate volume for every day of the experiment
  print("Interpolating stem volume")
  vol.all <- gapfillvol(vol)
  
  
  #- get branch mass
  branchMass <- interpolateBranchDW(vol.all=vol.all,plotson=0) 
  stem_branch_mass <- merge(vol.all,branchMass,by=c("chamber","Date"))
  
  
  
  #- now get an estimate of leaf mass for every day of the experiment
  # get the interpolated leaf areas (laDaily)
  
  dat.hr <- data.frame(data.table::fread("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv"))
  dat.hr$DateTime <- as.POSIXct(dat.hr$DateTime,format="%Y-%m-%d %H:%M:%S",tz="GMT")
  dat.hr$Date <- as.Date(dat.hr$DateTime)
  dat.hr$chamber <- as.factor(dat.hr$chamber)
  laDaily <- summaryBy(leafArea~Date+chamber,data=dat.hr,keep.names=T,na.rm=T)
  
  # get a canopy-weighted SLA estimate
  #setwd("//ad.uws.edu.au/dfshare/HomesHWK$/30035219/My Documents/Work/HFE/WTC3/gx_wtc3/")
  harvest <- read.csv("data/WTC_TEMP_CM_HARVEST-CANOPY_20140526-20140528_L1_v1.csv")
  leafmass <- summaryBy(TotalLeafDM~chamber,data=harvest,FUN=sum,keep.names=F)
  harvest2 <- merge(harvest,leafmass,by="chamber")
  harvest2$weights <- with(harvest2,TotalLeafDM/TotalLeafDM.sum)
  SLA <- plyr::ddply(harvest2,~chamber,summarise, SLA=weighted.mean(SLA,weights))
  
  # merge SLA and leaf area to estimate leaf mass for each day
  leafMass <- merge(laDaily,SLA,by=c("chamber"))
  leafMass$leafMass <- with(leafMass,leafArea/SLA*10000)
  
  treeMass1 <- merge(stem_branch_mass,leafMass,by=c("chamber","Date"))
  treeMass1$boleMass <- with(treeMass1,mass_wood+mass_bark)
  treeMass1$totMass <- rowSums(treeMass1[,c("boleMass","branchMass","leafMass")])
  treeMass <-subset(treeMass1,select=c("chamber","T_treatment","Water_treatment","Date","Measurement","Days_since_transplanting","branchMass","boleMass","leafMass","totMass","SLA"))
  
  print("Done. Returned treeMass")
  return(treeMass)
}
#----------------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#- Read and process the diameter and height dataset, return a dataframe of gapfilled size data.
returnd2h <- function(){
  
  #------------------------------------------------------------------------------------------------------------
  # get the tree size data
  size <- read.csv("data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
  size$chamber_n <- as.numeric(substr(size$chamber,start=2,stop=3))
  size$DateTime <- as.Date(size$DateTime)
  
  #do some processing to get d2h
  size2 <- base::subset(size,select=c("DateTime","Days_since_transplanting","chamber_n","chamber","T_treatment",
                                      "Water_treatment","Plant_height","Stem_number"))
  size2$diam <- size[,18]
  size2$diam_15 <- size[,16]
  size3 <- subset(size2,chamber_n<=12 & Stem_number==1)
  
  
  
  # tree 11 is different, because the height of stem 1 is not the actual true height. put the max height of the two stems in as the Plant_height
  tree11 <- subset(size2,chamber_n==11)
  tree11.max <- summaryBy(Plant_height~chamber_n+DateTime+Days_since_transplanting,data=tree11,FUN=max,na.rm=T,keep.names=T)
  size3[which(size3$chamber_n==11),"Plant_height"] <- tree11.max$Plant_height
  
  
  #- convert diameters to cm, heights to m
  size3$diam <- size3$diam/10
  size3$diam_15 <- size3$diam_15/10
  size3$Plant_height <- size3$Plant_height/100
  
  #get rid of NA's
  size4 <- subset(size3,diam>0)
  size4 <- size4[with(size4,order(DateTime,chamber_n)),]
  
  #get just the data with diameter at 15
  size_small <- subset(size3,diam_15>0 & year(DateTime)<2014)
  #------------------------------------------------------------------------------------------------------------
  
  
  
  #------------------------------------------------------------------------------------------------------------
  # gap fill diameter and height data
  
  # create dataframe for all days
  alldates <- rep(seq.Date(from=as.Date(range(size4$DateTime)[1]),to=as.Date(range(size4$DateTime)[2]),by="day"),12)
  chamber_n <- c(rep(1,length(alldates)/12),rep(2,length(alldates)/12),rep(3,length(alldates)/12),rep(4,length(alldates)/12),
                 rep(5,length(alldates)/12),rep(6,length(alldates)/12),rep(7,length(alldates)/12),rep(8,length(alldates)/12),
                 rep(9,length(alldates)/12),rep(10,length(alldates)/12),rep(11,length(alldates)/12),rep(12,length(alldates)/12))
  datedf <- data.frame(chamber_n=chamber_n,DateTime=alldates)                                                                                             
  
  #merge data in with dataframe of all days
  size5 <-merge(size4,datedf,all=T,by=c("chamber_n","DateTime"))
  
  #break across list, gapfill list
  size6 <- zoo(size5)
  size6$Days_since_transplanting <- na.approx(size6$Days_since_transplanting)
  size6$Plant_height <- na.approx(size6$Plant_height)
  size6$diam <- na.approx(size6$diam)
  
  
  # get it back to a normal dataframe
  size7 <- numericdfr(fortify.zoo(size6))
  size7$Date <- as.Date(size7$DateTime)
  size7$d2h <- with(size7,(diam/10)^2*Plant_height)
  
  # put some of the other bits back together
  size7$chamber <- as.factor(paste0("C",sprintf("%02.0f",size7$chamber_n)))
  size7$T_treatment <- as.factor(ifelse(size7$chamber_n %% 2 ==1,"ambient","elevated"))
  size7$Water_treatment <- "control"
  size7$Index <- NULL
  size7$Stem_number <- NULL
  
  #- establish a treatment key for the Water_treatment variable
  key <- data.frame(chamber=levels(size7$chamber),
                    Water_treatment=c("drydown","control","drydown","drydown","control","drydown",
                                      "control","drydown","control","control","drydown","control"))
  
  size_before <- subset(size7,Date<as.Date("2014-02-4"))
  size_after <- subset(size7,Date>=as.Date("2014-02-4"))
  size_after$Water_treatment <- NULL
  size_after2 <- merge(size_after,key,by="chamber")
  
  #- combined dataframes from before and after the drought began
  size8 <- rbind(size_before,size_after2)
  
  #- clean up dataframe for output
  size_out <- size8[,c("Date","chamber","T_treatment","Water_treatment","diam","Plant_height","d2h")]
  
  
  return(size_out)
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Function for converting all columns in a dataframe to numeric
numericdfr <- function(dfr){
  for(i in 1:ncol(dfr)){
    options(warn=-1)
    oldna <- sum(is.na(dfr[,i]))
    num <- as.numeric(as.character(dfr[,i]))
    if(sum(is.na(num)) > oldna)
      next
    else
      dfr[,i] <- num
  }
  options(warn=0)
  return(dfr)
}

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# get an estimate of bole volume for every measurement day
getvol <- function(){
  
  
  # get the data. the "interpolated 38" csv was manual editied in excel to linearly interpolate the UPPER stem estimates
  #   on measurement date 38 as the mean of measurements #37 and #39.
  #downloadCSV(filename="Wdata/from HIEv/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV,topath="data/from HIEv/")
  #size <- read.csv("data/from HIEv/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
  size <- read.csv("WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2_interpolated38.csv")
  names(size)[15:47] <- paste("d",substr(names(size)[15:47],start=2,stop=5),sep="")
  size <- size[,1:49]
  size$DateTime <- as.Date(size$DateTime)
  
  #ignore the reference plots
  size2 <- subset(size,T_treatment !="reference")
  size2$T_treatment <- factor(size2$T_treatment)
  size2$Water_treatment <- factor(size2$Water_treatment)
  size2$chamber <- factor(size2$chamber)
  
  #melt into "long" format
  sLong <- reshape2::melt(size2,measure.vars=names(size)[15:47],value.name="diam",variable.name="height")
  sLong$height_n <- as.numeric(substr(sLong$height,start=2,stop=4))
  
  sLong <- subset(sLong,height_n>=65)
  
  
  
  # estimate volume
  # split long dataframe into a list for each chamber and each measurement date
  sLong$chamber_date <- paste(sLong$chamber,sLong$Measurement,sep="_")
  sLong.l <- split(sLong,sLong$chamber_date)
  
  #-- add a number for the ground (assume no taper) and a value for the maximum tree height (assume diameter of 0.1mm)
  for (i in 1:length(sLong.l)){
    # add a line to the dataframe for the diameter at floor height
    firstline <- sLong.l[[i]][1,]  
    firstline$height_n[1] <- 0 #. Edited to give Mike an estimate of total tree volume to the ground.
    
    # add a line to the dataframe for the tiny diameter at total plant height
    lastline <- sLong.l[[i]][nrow(sLong.l[[i]]),]
    lastline$height_n[1] <- lastline$Plant_height[1]
    lastline$diam[1] <- 0.1 
    sLong.l[[i]] <- rbind(firstline,sLong.l[[i]],lastline)
    sLong.l[[i]] <- sLong.l[[i]][which(is.na(sLong.l[[i]]$diam)==F),] # get rid of NA's
  }
  treevol(sLong.l[[100]]) # example of volume calculation for a single observation of a single tree
  vols <- do.call(rbind,lapply(sLong.l,FUN=treevol))
  return(vols)
}

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
# estimate tree wood and bark volume and wood and bark mass from a series of measurements of diameter, and the height of that diameter measurement
treevol <- function(dat){
  dat2 <- subset(dat,diam>0)
  
  #-- get the bark depth estimates (along with wood and bark densities)
  density <- getBarkWoodDensity()
  density$barkDepth <- with(density,(diamoverbark-diamunderbark)/2)
  density$log_barkDepth <- log10(density$barkDepth)
  density$log_diamoverbark <- log10(density$diamoverbark)
  lm_bd <- lm(log_barkDepth~log_diamoverbark,data=density) # fit the harvest data
  dens.coefs <- unname(coef(lm_bd))
  dat2$bark_depth <- 10^(dens.coefs[1]+dens.coefs[2]*log10(dat2$diam)) # apply bark depths to estimate inner diameters
  dat2$diam_inner <- with(dat2,diam-bark_depth)  
  
  # estimate density from diameter
  bd.coef <- unname(coef(lm(bark_density~diamoverbark,data=density)))
  wd.coef <- unname(coef(lm(wooddensity~diamoverbark,data=density)))
  
  
  #loop over each stem segment and calculate volume as the fustrum of a cone. Assume the (unmeasured) bottom bit is a cylinder.
  vol <- 0
  vol_inner <- 0
  vol.cum <- 0
  vol.cum.inner <- 0
  mass_wood <- vol_wood <- 0
  mass_bark <- vol_bark <- 0
  for (i in 1:nrow(dat2)){
    # get the vertical segment length
    if(i==1){h1 <- 0}
    if(i > 1){h1 <- dat2[i-1,"height_n"]}
    h2 <- dat2[i,"height_n"]
    h <- h2-h1
    
    # get the radii in cm
    if(i==1){r1 <- dat2[i,"diam"]/20}
    if(i > 1){r1 <- dat2[i-1,"diam"]/20}
    r2 <- dat2[i,"diam"]/20
    
    #outer volume
    vol[i] <- pi*h/3* (r1^2+r1*r2+r2^2) # volume in cm3, fustrum of a cone. See http://jwilson.coe.uga.edu/emt725/Frustum/Frustum.cone.html
    
    
    # get the inner radii in cm
    
    if(i==1){r1.inner <- dat2[i,"diam_inner"]/20}
    if(i > 1){r1.inner <- dat2[i-1,"diam_inner"]/20}
    r2.inner <- dat2[i,"diam_inner"]/20
    
    vol_inner[i] <- pi*h/3* (r1.inner^2+r1.inner*r2.inner+r2.inner^2) # volume in cm3, fustrum of a cone. See http://jwilson.coe.uga.edu/emt725/Frustum/Frustum.cone.html
    
    # get wood volume and mass
    vol_wood[i] <- vol_inner[i]
    mass_wood[i] <- vol_wood[i]*(wd.coef[1]+wd.coef[2]*mean(r1.inner,r2.inner)) #wood mass is volume time density
    
    # get bark volume and mass
    vol_bark[i] <- vol[i] - vol_inner[i]
    mass_bark[i] <- vol_bark[i]*(bd.coef[1]+bd.coef[2]*mean(r1.inner,r2.inner)) #wood mass is volume time density
  }
  
  
  df.out <- dat2[1,c(1,2,3,4,5,6)]
  df.out$diam <- max(dat2$diam)
  df.out$height <- max(dat2$Plant_height)
  df.out$vol <- sum(vol,na.rm=T)
  df.out$vol_wood <- sum(vol_wood,na.rm=T)
  df.out$vol_bark <- sum(vol_bark,na.rm=T)
  df.out$mass_wood <- sum(mass_wood,na.rm=T)
  df.out$mass_bark <- sum(mass_bark,na.rm=T)
  
  return(df.out)
}

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------
#- function to return the measured stemwood and bark densities
getBarkWoodDensity <- function(){
  
  #--- Bark density data
  bark <- read.csv("data/WTC_TEMP_CM_BARKDENSITY_20140528_L1.csv")
  
  
  #--- wood density data
  wood <- read.csv("data/WTC_TEMP_CM_WOODDENSITY_20140528_L1.csv")
  
  wood$chamber_n <- as.numeric(substr(wood$chamber,start=2,stop=3))
  wood$T_treatment <- as.factor(ifelse(wood$chamber_n %% 2 == 1, "ambient","elevated"))
  wood2 <- subset(wood,Layer != "G Tip")
  wood2$Layer <- factor(wood2$Layer)
  wood2$position <- as.factor(ifelse(wood2$Layer == "Base","low",
                                     ifelse(wood2$Layer == "Middle","mid","top")))
  
  #-- merge wood and bark data
  dense <- merge(wood2,bark,by=c("chamber","position","T_treatment"))
  
  dense_out <- dense[,c("chamber","position","T_treatment","wooddensity","bark_density","diamoverbark","diamunderbark")]
  return(dense_out)
}

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# gapfill stem volume (for daily data)
gapfillvol <- function(dat){
  
  
  
  dat$chamber_n <- as.numeric(substr(dat$chamber,start=2,stop=3))
  
  # create dataframe for all days
  alldates <- rep(seq.Date(from=as.Date(range(dat$DateTime)[1]),to=as.Date(range(dat$DateTime)[2]),by="day"),12)
  chamber_n <- c(rep(1,length(alldates)/12),rep(2,length(alldates)/12),rep(3,length(alldates)/12),rep(4,length(alldates)/12),
                 rep(5,length(alldates)/12),rep(6,length(alldates)/12),rep(7,length(alldates)/12),rep(8,length(alldates)/12),
                 rep(9,length(alldates)/12),rep(10,length(alldates)/12),rep(11,length(alldates)/12),rep(12,length(alldates)/12))
  datedf <- data.frame(chamber_n=chamber_n,DateTime=alldates)                                                                                             
  
  #merge data in with dataframe of all days, sort it
  dat2 <- merge(dat,datedf,all=T,by=c("chamber_n","DateTime"))
  dat2 <- dat2[with(dat2,order(chamber_n,DateTime)),]
  
  #break across list, gapfill list
  dat3 <- zoo(dat2)
  dat3$Days_since_transplanting <- na.approx(dat3$Days_since_transplanting)
  dat3$vol <- na.approx(dat3$vol)
  dat3$vol_wood <- na.approx(dat3$vol_wood)
  dat3$vol_bark <- na.approx(dat3$vol_bark)
  dat3$mass_wood <- na.approx(dat3$mass_wood)
  dat3$mass_bark <- na.approx(dat3$mass_bark)
  
  
  
  # get it back to a normal dataframe
  dat4 <- numericdfr(fortify.zoo(dat3))
  dat4$DateTime <- as.Date(dat4$DateTime)
  dat4$Date <- dat4$DateTime
  
  # put some of the other bits back together
  dat4$chamber <- as.factor(paste0("C",sprintf("%02.0f",dat4$chamber_n)))
  dat4$T_treatment <- as.factor(ifelse(dat4$chamber_n %% 2 ==1,"ambient","elevated"))
  dat4$Index <- NULL
  dat4$Stem_number <- NULL
  
  
  #- establish a treatment key for the Water_treatment variable
  key <- data.frame(chamber=levels(dat4$chamber),
                    Water_treatment=c("drydown","control","drydown","drydown","control","drydown",
                                      "control","drydown","control","control","drydown","control"))
  
  size_before <- subset(dat4,Date<as.Date("2014-02-4"))
  size_before$Water_treatment <- "control"
  
  size_after <- subset(dat4,Date>=as.Date("2014-02-4"))
  size_after$Water_treatment <- NULL
  size_after2 <- merge(size_after,key,by="chamber")
  
  #- combined dataframes from before and after the drought began
  dat5 <- rbind(size_before,size_after2)
  
  #- clean up for exporting
  dat5_out <- dat5[,c("chamber","Date","T_treatment","Water_treatment","Measurement","Days_since_transplanting","vol","vol_wood","vol_bark","mass_wood","mass_bark")]
  
  return(dat5_out)
  
}

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Estimate branch biomass through time in WTC3.
interpolateBranchDW <- function(vol.all=vol.all,plotson=0){
  
  bc <- read.csv("data/WTC_TEMP_CM_BRANCHCENSUS_20130910-20140516_L0_v1.csv")
  bc$Date <- as.POSIXct(bc$Date,format="%d/%m/%Y",tz="GMT")
  bc$Date <- as.Date(bc$Date)
  bc$branchid <- as.factor(paste(bc$chamber,bc$branchnumber,sep="-"))
  
  #--------------------------------------------------------------------------------------------------
  # okay, so now can we estimate branch mass allometrically at the harvest? 
  branches <- read.csv("data/WTC_TEMP_CM_GX-RBRANCH_20140513-20140522_L1_v1.csv")
  
  br.allom <- lm(log10(branch_dw)~log10(diam.m),data=branches) # fit log-log branch allometry
  if(plotson==1){
    # windows(12,12)
    x11(12,12)
    plotBy(branch_dw~diam.m|T_treatment,data=branches,log="xy",pch=15,col=c("black","red"),cex.lab=1.5,cex.axis=1.3,
           xlab= "Branch diameter (mm)",ylab="Branch mass (g)");abline(br.allom)
  }
  br.coef <- coef(br.allom) # get the coefficients
  
  #apply allometry to the branch census dataset
  bc$branch_dw <- 10^(br.coef[1]+br.coef[2]*log10(bc$branchdiameter))
  bc <- subset(bc,branch_dw>0)
  
  #get the layer heights
  layers <- read.csv("data/WTC_TEMP_CM_LAYERHEIGHT_20140519_L1.csv")
  cuts <- c(0,mean(layers$toplayer1fromfloor),mean(layers$toplayer2fromfloor))
  layer.df <- data.frame(layer=c("low","mid","top"),layerlimit=cuts)
  
  bc$layer <- ifelse(bc$heightinsertion<=cuts[2],"low",ifelse(bc$heightinsertion<=cuts[3],"mid","top"))
  #sum up estimate of branch mass for each date for each chamber
  bc.sum <- summaryBy(branch_dw~chamber+Date+layer,data=bc,FUN=sum,keep.names=T)
  palette(c("black","blue","green","forestgreen","darkmagenta","chocolate1","cyan","darkgrey","darkgoldenrod1",
            "brown1","darkred","deeppink"))
  #compare branch allometry number to the harvest number
  harvest <- read.csv("data/WTC_TEMP_CM_HARVEST-CANOPY_20140526-20140528_L1_v1.csv")
  harvest.sum <- summaryBy(BranchDM~chamber+layer,data=harvest,FUN=sum,keep.names=T)
  harvest.sum$branch_allom <- subset(bc.sum,Date==as.Date("2014-05-16"))$branch_dw
  
  if(plotson==1){
    X11(12,12)
    plotBy(branch_allom~BranchDM|layer,data=harvest.sum,xlab="Measured branch mass (g)",ylab="Allometric estimate (g)",cex.lab=1.3,xlim=c(0,1000),ylim=c(0,1000),pch=15);abline(0,1)
    textxy(X=harvest.sum$BranchDM,Y=harvest.sum$branch_allom,labs=harvest.sum$chamber,cex=0.9)
  }
  #--------------------------------------------------------------------------------------------------
  
  #--------------------------------------------------------------------------------------------------
  
  
  bc.sum <- summaryBy(branch_dw~chamber+Date,data=bc,FUN=sum,keep.names=T)
  
  size2 <- merge(bc.sum,vol.all,by=c("Date","chamber"),all=T)
  
  #there is a nice relationship between tree (stem) volume and total branch mass
  if (plotson==1){
    X11(12,12)
    plotBy(branch_dw~vol|chamber,pch=15,type="b",cex=1.2,data=size2,log="xy",xlim=c(300,21000),
           xlab="Stem volume (cm3)",ylab="Branch mass (g)",cex.lab=1.3)
  }
  
  #- fit log-linear regression for each chamber, but all at once
  lm.all <- lm(log10(branch_dw)~log10(vol)*chamber,data=size2)
  
  # fit log- linear regression for each chamber
  dat.l <- split(size2,size2$chamber)
  lms1 <- list()
  for (i in 1:length(dat.l)){
    tofit <- dat.l[[i]]
    lms1[[i]] <- lm(log10(branch_dw)~log10(vol),data=subset(tofit,branch_dw>0))
    if (plotson==1) abline(lms1[[i]],col=i)
  }
  lms1 <- as.data.frame(do.call(rbind,lapply(lms1,coef)))
  lms1$chamber <- levels(size2$chamber)
  names(lms1) <- c("int","slope","chamber")
  
  
  # merge allometric estimates into the large size dataframe
  size3 <- merge(size2,lms1,by="chamber",all=T)
  size3$branch_dw_est <- with(size3,10^(int+slope*log10(vol)))
  #--------------------------------------------------------------------------------------------------
  
  if (plotson==1){
    X11(12,12)
    plotBy(branch_dw_est~Date|chamber,data=size3,ylab="Branch mass (g)",cex.lab=1.3)
  }
  
  
  out <- size3[,c("chamber","Date","branch_dw_est")]
  out <- out[with(out,order(chamber,Date)),]
  names(out)[3] <- "branchMass"
  return(out)
}
#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#- plot metrics of size over time.

plotsize <- function(output=T){
  
  
  
  #------------------------------------------------------------------------------------------------------------
  # #- get initial tree size 
  size <- read.csv("data/WTC_TEMP_CM_TREE-HEIGHT-DIAMETER_20121120-20140527_L1_V2.CSV")
  size$chamber_n <- as.numeric(substr(size$chamber,start=2,stop=3))
  size$DateTime <- as.Date(size$DateTime)
  
  summaryBy(Plant_height~T_treatment,data=subset(size,DateTime==as.Date("2012-12-12")),FUN=c(mean,standard.error))
  #------------------------------------------------------------------------------------------------------------
  
  
  
  
  #------------------------------------------------------------------------------------------------------------
  #- Get the three direct observations of leaf area. They happened on 9 Sept 2013, 10 Feb 2014, and
  #    during the harvest at ~25 May 2014.
  treeMass <- read.csv("data/WTC_TEMP_CM_WTCFLUX_20130914-20140526_L2_V2.csv")
  treeMass.sub <- subset(treeMass,as.Date(DateTime) %in% as.Date(c("2013-09-14","2014-02-10","2014-05-27")))
  treeMass.sub$Date <- as.Date(treeMass.sub$DateTime)
  leafArea <- summaryBy(leafArea~Date+chamber+T_treatment+Water_treatment,data=treeMass.sub,FUN=c(mean),keep.names=T)
  leafArea1 <- summaryBy(leafArea~Date+T_treatment+Water_treatment,data=leafArea,FUN=c(mean,standard.error))
  
  #-------------------
  #- get an estimate of volume and mass for wood and bark for each measurement day
  vol <- getvol()
  vol$diam <- vol$diam/10      # convert to cm
  vol$height <- vol$height/100 # convert to m
  vol$vol <- vol$vol/10000     # convert to m3
  
  size.m <- summaryBy(diam+height+vol~DateTime+T_treatment+Water_treatment,
                      data=subset(vol,Days_since_transplanting>0),FUN=c(mean,standard.error))
  
  
  
  # windows(40,40);
  par(mfrow=c(2,2),mar=c(0,5.5,0,3),oma=c(5,2,2,0),las=1,cex.axis=1.5)
  palette(c("#377EB8","#E41A1C"))
  
  #- plot diameter
  xlims <-as.Date(c("2013-3-1","2014-6-1"))
  yearlims <-as.Date(c("2013-1-1","2014-1-1"))
  plotBy(diam.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="control"),pch=16,type="o",ylim=c(0,8),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="control")$DateTime,
               y=subset(size.m,Water_treatment=="control")$diam.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="control")$diam.se,direction="updown",
               col=c("#377EB8","#E41A1C"))       
  plotBy(diam.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="drydown"),pch=1,type="o",add=T,legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="drydown")$DateTime,
               y=subset(size.m,Water_treatment=="drydown")$diam.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="drydown")$diam.se,direction="updown",
               col=c("#377EB8","#E41A1C"))    
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=F)
  title(ylab=expression(Diameter~(cm)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  #text(x=as.Date("2013-9-13"),y=8.8,labels="Flux measurements begin",xpd=NA,cex=1.3)
  legend("topleft",pch=c(16,16,1,1),lty=c(1),col=c("#377EB8","#E41A1C"),seg.len=1.5,
         legend=c("A-Wet","W-Wet","A-Dry","W-Dry"),bty="n",cex=1.2)
  legend("bottomright","a",bty="n",inset=0.002,cex=1.2)
  
  
  #-- plot leaf area over time
  plotBy(leafArea.mean~Date|T_treatment,data=subset(leafArea1,Water_treatment=="control"),pch=16,type="p",ylim=c(0,25),cex=1.5,
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(leafArea1,Water_treatment=="control")$Date,
               y=subset(leafArea1,Water_treatment=="control")$leafArea.mean,barlen=0.0,
               SE=subset(leafArea1,Water_treatment=="control")$leafArea.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  plotBy(leafArea.mean~Date|T_treatment,data=subset(leafArea1,Water_treatment=="drydown"),pch=1,type="p",add=T,legend=F,cex=1.5)
  adderrorbars(x=subset(leafArea1,Water_treatment=="drydown")$Date,
               y=subset(leafArea1,Water_treatment=="drydown")$leafArea.mean,barlen=0.0,
               SE=subset(leafArea1,Water_treatment=="drydown")$leafArea.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=F)
  title(ylab=expression(Total~leaf~area~(m^2)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","c",bty="n",inset=0.002,cex=1.2)
  
  
  #- plot height
  plotBy(height.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="control"),pch=16,type="o",ylim=c(0,11),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="control")$DateTime,
               y=subset(size.m,Water_treatment=="control")$height.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="control")$height.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  plotBy(height.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="drydown"),pch=1,type="o",add=T,legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="drydown")$DateTime,
               y=subset(size.m,Water_treatment=="drydown")$height.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="drydown")$height.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=T,las=2,
            format="%b")
  title(ylab=expression(Height~(m)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","b",bty="n",inset=0.002,cex=1.2)
  
  
  #- plot stem volume
  plotBy(vol.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="control"),pch=16,type="o",ylim=c(0,2),
         xlim=xlims,
         xaxt="n",yaxt="n",xlab="",ylab="",legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="control")$DateTime,
               y=subset(size.m,Water_treatment=="control")$vol.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="control")$vol.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  plotBy(vol.mean~DateTime|T_treatment,data=subset(size.m,Water_treatment=="drydown"),pch=1,type="o",add=T,legend=F)
  adderrorbars(x=subset(size.m,Water_treatment=="drydown")$DateTime,
               y=subset(size.m,Water_treatment=="drydown")$vol.mean,barlen=0.0,
               SE=subset(size.m,Water_treatment=="drydown")$vol.se,direction="updown",
               col=c("#377EB8","#E41A1C"))
  
  magaxis(side=c(2,4),labels=c(1,1),frame.plot=T,majorn=3,las=1)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="month"),tcl=0.25,labels=F)
  axis.Date(side=1,at=seq.Date(from=xlims[1],to=xlims[2],by="quarter"),tcl=0.75,labels=T,las=2,
            format="%b")
  title(ylab=expression(Stem~volume~(m^3)),cex.lab=2)
  abline(v=as.Date("2013-9-13"),lty=2)
  legend("bottomright","d",bty="n",inset=0.002,cex=1.2)
  
  
  
  if(output==T) dev.copy2pdf(file="output/treeSize.pdf")
}
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
