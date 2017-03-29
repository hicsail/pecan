############################
#
# script to generate driving data for SDGVM 
# from an ASCII .csv
#
# Awalker
# June 2013
#
############################

rm(list=ls())
library(ncdf4)

root <- '/mnt/disk2/'
dir  <- 'Research_Projects/FACE_modelling/Phase_3/'
wdd  <- paste(root,dir,'data/processed/',sep='')
wdo  <- paste(root,dir,'SDGVM/data/',sep='')

#input files
site <- c('NDFF','KSCO','PHAC','RHIN','DUKE','ORNL','PHAC_TAtrmt')
sia  <- 1:4
sia  <- 5:6
sia  <- c(3,7)

#site lat and long
lat  <- c(  36.77 ,  28.63 ,   41.18 ,  45.68 , 35.97 ,  35.9,   41.18 )
long <- c(-115.96 , -80.70 , -104.19 , -89.63 , 79.08 , 84.33, -104.19 )

#start year of the industrial period spin-up
ind_year <- 1850

#output directories files
out_dirs   <- c('standard','par','avgpar')
out_cols   <- list(1:6,1:7,c(1:6,8))
clim_ofile <- 'site.dat'

#days in month arrays
d_in_m    <- c(31,28,31,30,31,30,31,31,30,31,30,31)
d_in_m_ly <- c(31,29,31,30,31,30,31,31,30,31,30,31)

#column subscripts to translate .txt file column order to SDGVM file column order
# YEAR DOY Temp Rain RH SWdown
cia    <- c(1,2,4,3,5,9)
cnames <- c('YEAR','MOY','DOM','T','Rain','RH','SWR','meanSWR')

#function to write the readme file
write.readme<-function(lat=lat,lon=long,years){
  file <- 'readme.dat'
  loc  <- data.frame(lat=lat,lon=lon)
  write.table('SITED',file,row.names=F,col.names=F,quote=F)
  write.table('lat and long',file,append=T,row.names=F,col.names=F,quote=F)
  write.table(loc,file,append=T,row.names=F,col.names=F,quote=F)
  write.table('initial and final years of the dataset',file,append=T,row.names=F,col.names=F,quote=F)
  write.table(years,file,append=T,row.names=F,col.names=F,quote=F)
}

#function to determine leap years
f_leap_yr<-function(y){
  ifelse(y%%400==0,leap_yr<-T,
         ifelse(y%%100==0,leap_yr<-F,
                ifelse(y%%4==0,leap_yr<-T,leap_yr<-F)))
  leap_yr
}

# #function
# f_Q_to_RH<-function(i,df,Qair,Tair,psurf){
#   #convert specific humidity to RH - use with lapply
#   #Tair in K, psurf in Pa, Qair i
#   
#   Qair <- df$Qair[i]
#   Tair <- df$Tair[i]
#   psurf <- df$PSurf[i]
#   
#   #calculate vapour pressure 
#   vp   <- psurf / (0.622/Qair + 1)
#   #use buck 1981 eq to calculate sat vp
#   svp  <- 6.1121 * exp(17.502*(Tair-273.15) / (240.9+Tair-273.15)) * 100
#   
#   vp/svp*100
# }

#function to read and process data
read_process <- function(ifile,cia){
  header <- read.table(ifile,header=T,nrows=1)
  mydata <- read.table(ifile,skip=2)
  names(mydata) <- names(header)
  head(mydata)
  
  #convert to SDGVM
  daily_data     <- mydata[cia]
  # convert T to oC
  daily_data[,3] <- daily_data[,3] - 273.15
  #convert SW
  head(daily_data)
  daily_data$SWdown <- daily_data$SWdown * 1e6 / (3600 * 24) 
  
  #create month and DOM array including leap years - for 2012 to 2023 files
  days <- length(daily_data[,1])
  time <- matrix(nrow=days,ncol=3)
  meanSW <- 1:days 
  
  i<-1
  for(y in daily_data$YEAR[1]:daily_data$YEAR[days]){
    dy <- 1
    for(m in 1:12){
      for(d in 1:d_in_m[m]){
        time[i,]  <- c(y,m,d)
#        meanSW[i] <- xSW[dy,2]
        i  <-  i+1
        dy <- dy+1
      }
#       if(f_leap_yr(y)&m==2) {time[i,]<-c(y,2,29) ; meanSW[i]<-xSW[dy,2] ; i<-i+1 }
      if(f_leap_yr(y)&m==2) {time[i,]<-c(y,2,29) ; i<-i+1 }
    }
  }
#   data.frame(time,daily_data[,3:6],meanSW)
  data.frame(time,daily_data[,3:6],daily_data$SWdown)
}

#function to create industrial period metdata
ind_metdata <- function(sind=1850,real_data){
  sexp <- min(real_data$YEAR)
  eexp <- max(real_data$YEAR)
  eind <- sexp - 1
          
  #create start of ind period to experiment start yr-1 by repeating experimental period data in sequence
  ry <- sexp
  for(y in sind:eind){
    df      <- subset(real_data,YEAR==ry)
    df$YEAR <- y
    
    if(f_leap_yr(ry)&!f_leap_yr(y)){
      #remove feb 29th on non-leap years
      df <- df[-which(df$MOY==2&df$DOM==29),]
    } else if(!f_leap_yr(ry)&f_leap_yr(y)){
      #add feb 29th for leap years by repeating feb 28th
      if(f_leap_yr(y)) {df<-df[c(1:59,59:365),];df[60,3]<-29}            
    }    
    
    ifelse(y==sind , ind_data<-df , ind_data<-rbind(ind_data,df) )
    ry <- ry+1
    if(ry==eexp) ry <- sexp
  }
  
  ind_data
}


#################################
#process met data
# s <- 7
for(s in sia){
  setwd(paste(wdd,site[s],'/',sep=''))
  if(site[s]=='KSCO') setwd(paste(wdd,site[s],'_1996','/',sep=''))
  if(substr(site[s],1,4)=='PHAC') setwd(paste(wdd,site[s],'_160630','/',sep=''))
  
  print('')
  ifile <- paste(site[s],'_forcing_d.txt',sep='')
  print(ifile)
  
  #read and process data
  daily_data <- read_process(ifile,cia)
  names(daily_data)[1:3] <- cnames[1:3]
      
  #add industrial period data
  ind_data   <- ind_metdata(real_data=daily_data)
  
  #add chamber data to KSCO site 
  if(site[s]=='KSCO') {
    #replace data with inside chamber data @ KSCO
    setwd(paste(wdd,site[s],'_chamber_1996/',sep=''))
    print('')
    ifile <- paste(site[s],'_chamber_forcing_d.txt',sep='')
    print(ifile)
    daily_data_c     <- read_process(ifile,cia)
    names(daily_data_c)[1:3] <- cnames[1:3]
    daily_data_c[,5] <- daily_data[,5]
    daily_data       <- daily_data_c
  }
  
  #combine industrial period data with experimental period 
  daily_data <- rbind(ind_data,daily_data)
  names(daily_data)  <- cnames
  
  #write climate data
  setwd(wdo)
  dir.create('clim')
  setwd(paste(wdo,'clim/',sep=''))
  dir.create(site[s])
 
  # for(d in 1:length(out_dirs)){
  for(d in 2) {
    setwd(paste(wdo,'clim/',site[s],sep=''))
    dir.create(out_dirs[d])
    setwd(paste(wdo,'clim',site[s],out_dirs[d],sep='/'))

    #standard datset
    # if(d==1)     {   write.table(format(daily_data[,out_cols[[d]]],width=11,scientific=F),clim_ofile,row.names=F,col.names=F,quote=F)
    #                  years<-data.frame(start=ind_year,end=daily_data$YEAR[length(daily_data[,1])])                   
    #                  write.readme(years=years,lat=lat[s],lon=long[s]) 
    # }
    #dataset with SW radiation
    else if(d==2){   write.table(format(daily_data[,out_cols[[d]]],width=11,scientific=F),clim_ofile,row.names=F,col.names=F,quote=F)
                     years<-data.frame(start=ind_year,end=daily_data$YEAR[length(daily_data[,1])])
                     write.readme(years=years,lat=lat[s],lon=long[s])
    }
    #dataset with mean SW radiation
    # else         {   write.table(format(daily_data[,out_cols[[d]]],width=11,scientific=F),clim_ofile,row.names=F,col.names=F,quote=F)
    #                  years<-data.frame(start=ind_year,end=daily_data$YEAR[length(daily_data[,1])])
    #                  write.readme(years=years,lat=lat[s],lon=long[s])
    # }
  }
}
