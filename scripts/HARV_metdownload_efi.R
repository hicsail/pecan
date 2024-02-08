#!/usr/bin/env Rscript
#load libraries
# pecan_home example: /pecan
pecan_home <- Sys.getenv("PECAN_HOME")
# nc_dir example: /projectnb/dietzelab/ahelgeso/NOAA_met_data/
nc_dir <- Sys.getenv("NC_DIR")
# clim_dir example: /projectnb/dietzelab/ahelgeso/NOAA_met_data_CH1/
clim_dir <- Sys.getenv("CLIM_DIR")
source(paste0(pecan_home, "/modules/data.atmosphere/R/download_noaa_gefs_efi.R"))
source(paste0(pecan_home, "/modules/data.atmosphere/R/noaa_gefs_efi_helper.R"))

library(dplyr)
option_list = list(optparse::make_option("--start.date",
        default = toString(Sys.Date()-1),
        type="character"),
    optparse::make_option("--jumpback",
        default = 10,
        type="integer"),
    optparse::make_option("--jumpback.date",
        type="character")
)
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
startdate <- as.Date(args$start.date)
if(!is.null(args$jumpback.date))
  args$jumpback = as.numeric(startdate - as.Date(args$jumpback.date))
if(is.null(args$jumpback))
  args$jumpback = 10
jumpback = args$jumpback

#set fcn inputs
site.lat = 42.5
site.lon = -72.15
sitename <- "HARV"
siteid <- 646
runDays <- format(seq(from=as.POSIXlt(startdate) - lubridate::days(jumpback),
        to=as.POSIXlt(as.character(startdate)),
        by="days"),"%Y-%m-%d")
print("args set")

download_data <- function(runDays, t, sitename, site, nc_dir, clim_dir) {
  print(runDays[t])
  
  #download met using EFI fcns
  tmp = download.NOAA_GEFS_EFI(sitename = sitename,
      outfolder = nc_dir,
      start_date = runDays[t],
      site.lat = site.lat,
      site.lon = site.lon)
  #set up path for met2model
  output_path <- file.path(nc_dir, "noaa", "NOAAGEFS_1hr", sitename, runDays[t], "00")
  ########## Met2Model For SIPNET ##############
  outfolder = file.path(clim_dir, sitename, runDays[t])
  if(!dir.exists(outfolder)){dir.create(outfolder, recursive = TRUE)}
  
  in.path = output_path
  in.prefix = list.files(output_path)
  
  end_date = as.Date(runDays[t]) + lubridate::days(35)
  
  for(l in 1:length(in.prefix)){
    
    print(paste0("VARS in.path: ", in.path))
    print(paste0("VARS : in.prefix[l]", in.prefix[l]))
    print(paste0("VARS : outfolder", outfolder))
    print(paste0("VARS : runDays[t]", runDays[t]))
    print(paste0("VARS : end_date", end_date))
    PEcAn.SIPNET::met2model.SIPNET(in.path = in.path, 
        in.prefix = in.prefix[l], 
        outfolder = outfolder, 
        start_date = runDays[t], 
        end_date = end_date,
        overwrite = FALSE,
        verbose = FALSE, 
        year.fragment = TRUE) 
    
  } 
  print(paste0("Look for data in: ", outfolder))
}

for(t in seq_along(runDays)){
  tryCatch({
      tryCatch({
          download_data(runDays, t, sitename, site, nc_dir, clim_dir)
        }
        , error = function(ex) {
          if(!grepl("Path does not exist", toString(ex))) {
            num_seconds <- 10
            print(paste0("An error occured, will try again in ", num_seconds, " seconds: ", ex))
            Sys.sleep(num_seconds)
            download_data(runDays, t, sitename, site, nc_dir, clim_dir)
          }
        }
      )
    }
    , error = function(ex) {
      num_seconds <- 20
      print(paste0("An error occured, will try again in ", num_seconds, " seconds: ", ex))
      Sys.sleep(num_seconds)
      download_data(runDays, t, sitename, site, nc_dir, clim_dir)
    }
  )
  
}
print("done")

