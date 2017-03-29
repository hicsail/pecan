###############################
#
# Create CO2 and Ndep files for SDGVM and CLM4
# eCO2 modelling
#
# AWalker 
# Jan 2014
#
###############################

rm(list=ls())
library(ncdf4)

do_co2  <- T
do_ndep <- F

CLM     <- F
SDGVM   <- T

wd    <- '/mnt/disk2/Research_Projects/FACE_modelling/Phase_3/data/processed/'
site  <- c('NDFF','KSCO','PHAC','RHIN','DUKE','ORNL')
siten <- c('NDF' ,'KSC' ,'PHA' ,'RHI' ,'DUK' ,'ORN')
slat  <- c(126.77,118.63,131.18,135.68,125.97,125.9) # degrees N
slon  <- c(64.04,99.3,75.1,90.47,100.92,95.67)       # degrees E - from date line?
slong <- slon + 180                                  # degrees E - from gm
sia   <- 5:6
# sia   <- 5:6
sia   <- 3
sowd  <- '/mnt/disk2/Research_Projects/FACE_modelling/Phase_3/SDGVM/data/co2/'
nwd   <- '/mnt/disk2/Research_Projects/FACE_modelling/Phase_3/data/Ndeposition/'


#functions
f_leap_yr<-function(y){
  ifelse(y%%400==0,leap_yr<-T,
         ifelse(y%%100==0,leap_yr<-F,
                ifelse(y%%4==0,leap_yr<-T,leap_yr<-F)))
  leap_yr
}


#create month and DOM array for normal and leap years
#days in month arrays
d_in_m    <- c(31,28,31,30,31,30,31,31,30,31,30,31)
d_in_m_ly <- c(31,29,31,30,31,30,31,31,30,31,30,31)
for(m in 1:12) { 
  ifelse(m==1,ms<-rep(m,each=d_in_m[m]),ms<-c(ms,rep(m,each=d_in_m[m])))
  ifelse(m==1,ds<-1:d_in_m[m],ds<-c(ds,1:d_in_m[m]))
  ifelse(m==1,msl<-rep(m,each=d_in_m_ly[m]),msl<-c(msl,rep(m,each=d_in_m_ly[m])))
  ifelse(m==1,dsl<-1:d_in_m_ly[m],dsl<-c(dsl,1:d_in_m_ly[m]))
}


#open CO2 data files
s <- 3
for (s in sia){
  print(site[s])
  if(do_co2){
    
    #open CO2 template file
    setwd(wd)
    nctemplate <- nc_open('fco2_datm_1765-2007_c100614.nc')
    print(nctemplate)
    co2_temp        <- ncvar_get(nctemplate,'CO2')

    
    #open mean treatment CO2 file
    mean_co2   <- read.table('site_mean_CO2.txt',header=T)
    
    
    ### Process CO2 data ###
    #####
    wdd     <- paste(wd,site[s],sep='')
    if(site[s]=='KSCO') wdd <- paste(wd,site[s],'_1996',sep='')
    if(site[s]=='PHAC') wdd <- paste(wd,site[s],'_150327',sep='')
    setwd(wdd)
    
    
    #open historical co2 file
    hco2file <- paste(site[s],'_forcing_y.txt',sep='')
    top      <- read.table(hco2file,header=T)
    hco2     <- read.table(hco2file,skip=2)
    names(hco2) <- names(top)
    
    
    #open treatment co2 file
    co2file  <- paste(site[s],'_forcing_d.nc',sep='')
    mynetcdf <- nc_open(co2file)
    print(mynetcdf)
    atm_press<- mean(ncvar_get(mynetcdf,'Psurf')) * 1e-6
    aco2     <- ncvar_get(mynetcdf,'aCO2')
    eco2     <- ncvar_get(mynetcdf,'eCO2')
    years    <- ncvar_get(mynetcdf,'YEAR')
    doy      <- ncvar_get(mynetcdf,'DOY')
    aco2_nl  <- aco2[-which(doy==366)]
    eco2_nl  <- eco2[-which(doy==366)]
    lat      <- ncvar_get(mynetcdf,'nav_lat')
    lon      <- ncvar_get(mynetcdf,'nav_lon')
    nc_close(mynetcdf)
    
    
    #process co2 data
    sr          <- which(hco2$YEAR==1765)
    er          <- which(hco2$YEAR==min(years)-1)
    co2data_nl  <- rep(hco2$Global_Mean_CO2[sr:er],each=365)
    co2data     <- rep(hco2$Global_Mean_CO2[which(hco2$YEAR==1850)],365)
    for(y in 1851:(min(years)-1)){
      if(f_leap_yr(y)) co2data <- c(co2data,rep(hco2$Global_Mean_CO2[which(hco2$YEAR==y)],366))
      else             co2data <- c(co2data,rep(hco2$Global_Mean_CO2[which(hco2$YEAR==y)],365))
    }
    
    
    #SDGVM data
    #create time arrays for SDGVM full year sequence
    if(SDGVM){
      f_leap_yr(1850)
      time <- cbind(1850,ms,ds)
      for(y in 1851:max(years)){
        if(f_leap_yr(y)) time <- rbind(time,cbind(y,msl,dsl))
        else             time <- rbind(time,cbind(y,ms,ds))
      }
      d <- length(c(co2data,eco2))
      aco2data    <- data.frame(time,c(co2data,aco2)*atm_press)  
      eco2data    <- data.frame(time,c(co2data,eco2)*atm_press)  
      tindco2data <- data.frame(time,rep(285,d)*atm_press)  
      tambco2data <- data.frame(time,rep(mean_co2[mean_co2$site==site[s]&mean_co2$treatment=='a',3],d)*atm_press)  
      telvco2data <- data.frame(time,rep(mean_co2[mean_co2$site==site[s]&mean_co2$treatment=='e',3],d)*atm_press)  
      setwd(sowd)
      write.table(format(aco2data,width=4),paste(site[s],'_co2_amb_daily.dat',sep=''),row.names=F,col.names=F,quote=F)
      write.table(format(eco2data,width=4),paste(site[s],'_co2_ele_daily.dat',sep=''),row.names=F,col.names=F,quote=F)
      write.table(format(tindco2data,width=4),paste(site[s],'_co2_tind_daily.dat',sep=''),row.names=F,col.names=F,quote=F)
      write.table(format(tambco2data,width=4),paste(site[s],'_co2_tamb_daily.dat',sep=''),row.names=F,col.names=F,quote=F)
      write.table(format(telvco2data,width=4),paste(site[s],'_co2_tele_daily.dat',sep=''),row.names=F,col.names=F,quote=F)
    }
    
    
    #CLM data
    if(CLM){
      setwd(wdd)
      days <- 365 * (max(years)-1764)
      
      dnv   <- ncdim_def('nv',units='',vals=1:4,longname='')
      dlon  <- ncdim_def('lon',units='',vals=1,longname='')
      dlat  <- ncdim_def('lat',units='',vals=1,longname='')
      dtime <- ncdim_def('Time',units='days since 1765-01-01',vals=1:days,calendar='noleap',longname='Time')
      
      lonc <- ncvar_def('lonc',units='degrees E',dim=list(dlon,dlat),prec='double') 
      latc <- ncvar_def('latc',units='degrees N',dim=list(dlon,dlat),prec='double')
      lonv <- ncvar_def('lonv',units='',dim=list(dlon,dlat,dnv),prec='double')
      latv <- ncvar_def('latv',units='',dim=list(dlon,dlat,dnv),prec='double')
      mask <- ncvar_def('mask',units='',dim=list(dlon,dlat),prec='double')
      frac <- ncvar_def('frac',units='',dim=list(dlon,dlat),prec='double')
      area <- ncvar_def('area',units='',dim=list(dlon,dlat),prec='double')
      CO2  <- ncvar_def('CO2',units='ppmv',dim=list(dlon,dlat,dtime),prec='float',longname='CO2 concentration')
      time <- ncvar_def('time',units='days since 1765-01-01',dim=list(dtime),prec='integer',longname='time')
      
      vars    <- list(lonc,latc,lonv,latv,mask,frac,area,CO2,time)
      co2file <- paste('fco2_datm_1765-',max(years),'_US-',siten[s],'_control.nc',sep='')
      newnc   <- nc_create(co2file,vars)
      
      ncatt_put(newnc,varid='CO2',attname='coordinate',attval='latc lonc')
      ncatt_put(newnc,varid='time',attname='calendar',attval='noleap')
      
      #  RCODE = NF90_DEF_VAR(NCID, 'date', NF90_INT, dimid(4), varidout(10))
      #  RCODE = NF90_PUT_ATT(NCID, varidout(10), 'long_name', 'current date as yyyymmdd')
      
      print(newnc)
      ncvar_put(newnc,varid='lonc',vals=180)
      ncvar_put(newnc,varid='latc',vals=0)
      ncvar_put(newnc,varid='lonv',vals=c(0,360,360,0))
      ncvar_put(newnc,varid='latv',vals=c(90,90,-90,-90))
      ncvar_put(newnc,varid='mask',vals=1)
      ncvar_put(newnc,varid='frac',vals=1)
      ncvar_put(newnc,varid='area',vals=12.5663706143592)
      ncvar_put(newnc,varid='CO2',vals=c(co2data_nl,aco2_nl))
      ncvar_put(newnc,varid='time',vals=1:days)
      nc_close(newnc)
      
      
      #copy new ncdf file
      system(paste('cp ',co2file,paste('fco2_datm_1765-',max(years),'_US-',siten[s],'.nc',sep='')))
      system(paste('cp ',co2file,paste('fco2_datm_1765-',max(years),'_US-',siten[s],'_tind.nc',sep='')))
      system(paste('cp ',co2file,paste('fco2_datm_1765-',max(years),'_US-',siten[s],'_tamb.nc',sep='')))
      system(paste('cp ',co2file,paste('fco2_datm_1765-',max(years),'_US-',siten[s],'_telv.nc',sep='')))
      
      
      #open and re-write co2 in eCO2 netcdf file
      eco2nc <- nc_open(paste('fco2_datm_1765-',max(years),'_US-',siten[s],'.nc',sep=''),write=T)
      ncvar_put(eco2nc,varid='CO2',vals=c(co2data_nl,eco2_nl))
      nc_close(eco2nc)
      
      d <- length(c(co2data_nl,eco2_nl))
      
      eco2nc <- nc_open(paste('fco2_datm_1765-',max(years),'_US-',siten[s],'_tind.nc',sep=''),write=T)
      ncvar_put(eco2nc,varid='CO2',vals=rep(285,d))
      nc_close(eco2nc)
      
      eco2nc <- nc_open(paste('fco2_datm_1765-',max(years),'_US-',siten[s],'_tamb.nc',sep=''),write=T)
      ncvar_put(eco2nc,varid='CO2',vals=rep(mean_co2[mean_co2$site==site[s]&mean_co2$treatment=='a',3],d))
      nc_close(eco2nc)
      
      eco2nc <- nc_open(paste('fco2_datm_1765-',max(years),'_US-',siten[s],'_telv.nc',sep=''),write=T)
      ncvar_put(eco2nc,varid='CO2',vals=rep(mean_co2[mean_co2$site==site[s]&mean_co2$treatment=='e',3],d))
      nc_close(eco2nc)
      
    }
  }

  if(do_ndep){
    ### Process N deposition data ###
    #####
    
    # it appears that althopugh in the CLM ndep .nc files ndep units are g/m2/yr
    # the units are really kg/ha/yr
    
    #for CLM
    setwd(nwd)

    #open mean treatment Ndep file
    mean_ndep  <- read.table('site_mean_Ndep.txt',header=T)
    
    #read N dep
    ndep <- read.csv('ndep_face_sites_exp.csv',header=T,sep=';')
    c    <- which(names(ndep)==site[s])
    r    <- which(ndep[,1]==2023)
#     r    <- 157
    
    #     ndepfile  <- 'fndep_clm_simyr1849-2023_1x1pt_AU-EuF_fromXYK.nc'  
    ndepfile  <- 'fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc' # lat res is 1.875
    ndep_root <- paste('fndep_clm_simyr1849-2023_1x1pt_US-',siten[s],sep='')
        
    #copy new ncdf file
    system(paste('cp ',ndepfile,paste(ndep_root,'.nc',sep='')))
    system(paste('cp ',ndepfile,paste(ndep_root,'_tind.nc',sep='')))
    system(paste('cp ',ndepfile,paste(ndep_root,'_ttre.nc',sep='')))
    
    ndepnc <- nc_open(paste(ndep_root,'.nc',sep=''),write=T)
    #calculate grid reference integer
    nlat <- ceiling(slat[s]/(180/length(ncvar_get(ndepnc,'lat'))))
    nlon <- ceiling(slong[s]/(360/length(ncvar_get(ndepnc,'lon'))))
    #add new nedep data to correct grid point    
    prior <- ncvar_get(ndepnc,varid='NDEP_year',start=c(nlon,nlat,1),count=c(1,1,157+1))
    for(nlon in 1:144){
      for(nlat in 1:96){
        ncvar_put(ndepnc,varid='NDEP_year',vals=ndep[c(1,1:r),c]/10,start=c(nlon,nlat,1),count=c(1,1,r+1))        
      }
    }
    #ncdays <- ncvar_get(ndepnc,varid='time')
    ncdays  <- seq(674885,738395,365)
    ncyears <- c(1849:2023)
    ncvar_put(ndepnc,varid='time',vals=ncdays)            
    ncvar_put(ndepnc,varid='YEAR',vals=ncyears)            
    nc_close(ndepnc)
    
    #     a      <- length(ncvar_get(ndepnc,'NDEP_year')) 
    ndepnc <- nc_open(paste(ndep_root,'_tind.nc',sep=''),write=T)
    ncvar_put(ndepnc,varid='NDEP_year',vals=rep(ndep[1,c]/10,r+1),start=c(nlon,nlat,1),count=c(1,1,r+1))
    ncvar_put(ndepnc,varid='time',vals=ncdays)            
    ncvar_put(ndepnc,varid='YEAR',vals=ncyears)            
    nc_close(ndepnc)
    
    ndepnc <- nc_open(paste(ndep_root,'_ttre.nc',sep=''),write=T)
    ncvar_put(ndepnc,varid='NDEP_year',vals=rep(mean_ndep[mean_ndep$site==site[s],2]/10,r+1),start=c(nlon,nlat,1),count=c(1,1,r+1))
    ncvar_put(ndepnc,varid='time',vals=ncdays)            
    ncvar_put(ndepnc,varid='YEAR',vals=ncyears)            
    nc_close(ndepnc)
  }  
  
#site loop
}


#QC
#################
s <- 6
if(do_co2){
  setwd(paste(wd,site[s],sep='/'))
  type <- 2
  lab <- c('_control','','_tind','_tamb','_telv')
  eyear <- c(2008,2006,2013,2008,2006,2008)
  eco2nc <- nc_open(paste('fco2_datm_1765-',eyear[s],'_US-',siten[s],lab[type],'.nc',sep=''))
  #eco2nc
  plot(ncvar_get(eco2nc,'CO2')~ncvar_get(eco2nc,'time'))
  nc_close(eco2nc)  
}

if(do_ndep){
  library(lattice)

  l <- 175
  xyplot(ndep[1:l,2]+ndep[1:l,3]+ndep[1:l,4]+ndep[1:l,5]+ndep[1:l,6]+ndep[1:l,7]+ndep[1:l,8]~1:l,auto.key=T)
    
  setwd(nwd)
  type <- 1
  nlab <- c('','_tind','_ttre')
  ndepfile  <- 'fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc' 
  ndepncdefault <- nc_open(ndepfile,write=T)
  
  s <- 6
  for(s in sia){
    ndep_root <- paste('fndep_clm_simyr1849-2023_1x1pt_US-',siten[s],sep='')
    print(ndep_root)
    ndepnc <- nc_open(paste(ndep_root,nlab[type],'.nc',sep=''),write=T)
    #calculate grid reference integer
    nlat <- ceiling(slat[s]/(180/length(ncvar_get(ndepnc,'lat'))))
    nlon <- ceiling(slong[s]/(360/length(ncvar_get(ndepnc,'lon'))))
    print(
      xyplot(c(ncvar_get(ndepncdefault,'NDEP_year',start=c(nlon,nlat,1),count=c(1,1,158)),
             ncvar_get(ndepnc,'NDEP_year',start=c(nlon,nlat,1),count=c(1,1,175)))~
             c(ncvar_get(ndepncdefault,'YEAR'),ncvar_get(ndepnc,'YEAR')),groups=c(rep(1:2,each=158),rep(2,17)),
             main=site[s],auto.key=T)
    )
    print(
      xyplot(ncvar_get(ndepnc,varid='NDEP_year',start=c(nlon,nlat,1),count=c(1,1,175))~
               ncvar_get(ndepnc,varid='YEAR'),main=site[s],auto.key=T)
    )
    nc_close(ndepnc)      
  }

  ncvar_get(ndepnc,varid='NDEP_year',start=c(101,12,1),count=c(1,1,175))
  ncvar_get(ndepnc,varid='NDEP_year',start=c(145,nlat,1),count=c(1,1,175))
  levelplot(ncvar_get(ndepnc,varid='NDEP_year',start=c(1,1,1),count=c(144,96,1)))
  
}

# length(ncvar_get(ndepnc,'NDEP_year',start=c(nlon,nlat,1),count=c(1,1,175)))
# length(ncvar_get(ndepnc,'YEAR'))
# length(c(ncvar_get(ndepncdefault,'NDEP_year',start=c(nlon,nlat,1),count=c(1,1,158)),rep(NA,17)))

# ndepfile  <- 'fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc' 
# ndepnc    <- nc_open(ndepfile,write=T)
# nd        <- ncvar_get(ndepnc,varid='NDEP_year')
