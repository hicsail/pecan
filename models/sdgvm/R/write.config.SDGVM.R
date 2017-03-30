#-------------------------------------------------------------------------------
# Copyright (c) 2017 PEcAn Project
# All rights reserved. This program and the accompanying materials
# are made available under the terms of the 
# University of Illinois/NCSA Open Source License
# which accompanies this distribution, and is available at
# http://opensource.ncsa.illinois.edu/license.html
#-------------------------------------------------------------------------------

##-------------------------------------------------------------------------------------------------#
##' Writes a SDGVM config file.
##'
##' Requires a pft xml object, a list of trait values for a single model run,
##' and the name of the file to create
##'
##' @name write.config.SDGVM
##' @title Write SDGVM configuration files
##' @param defaults list of defaults to process
##' @param trait.values vector of samples for a given trait
##' @param settings list of settings from pecan settings file
##' @param run.id id of run
##' @return configuration file for SDGVM for given run
##' @author Mike Dietze, Rob Kooper
##' 
##' @export
##-------------------------------------------------------------------------------------------------#
write.config.SDGVM <- function(defaults, trait.values, settings, run.id) {

  # find out where to write run/ouput
  rundir <- file.path(settings$host$rundir, run.id)
  outdir <- file.path(settings$host$outdir, run.id)
  
  #-----------------------------------------------------------------------
  # create launch script (which will create symlink)
  if (!is.null(settings$model$jobtemplate) && file.exists(settings$model$jobtemplate)) {
    jobsh <- readLines(con = settings$model$jobtemplate, n = -1)
  } else {
    jobsh <- readLines(con = system.file("template.job", package = "PEcAn.SDGVM"), n = -1)
  }
  
  # create host specific settings
  hostsetup <- ""
  if (!is.null(settings$model$prerun)) {
    hostsetup <- paste(hostsetup, sep = "\n", paste(settings$model$prerun, collapse = "\n"))
  }
  if (!is.null(settings$host$prerun)) {
    hostsetup <- paste(hostsetup, sep = "\n", paste(settings$host$prerun, collapse = "\n"))
  }
  
  hostteardown <- ""
  if (!is.null(settings$model$postrun)) {
    hostteardown <- paste(hostteardown, sep = "\n", paste(settings$model$postrun, collapse = "\n"))
  }
  if (!is.null(settings$host$postrun)) {
    hostteardown <- paste(hostteardown, sep = "\n", paste(settings$host$postrun, collapse = "\n"))
  }
  
  # create job.sh
  jobsh <- gsub("@HOST_SETUP@", hostsetup, jobsh)
  jobsh <- gsub("@HOST_TEARDOWN@", hostteardown, jobsh)
  
  jobsh <- gsub("@OUTDIR@", outdir, jobsh)
  jobsh <- gsub("@RUNDIR@", rundir, jobsh)
  jobsh <- gsub("@RUNID@", run.id, jobsh)
  jobsh <- gsub("@BINARY@", settings$model$binary, jobsh)
  writeLines(jobsh, con = file.path(settings$rundir, run.id, "job.sh"))
  Sys.chmod(file.path(settings$rundir, run.id, "job.sh"))
  
  #-----------------------------------------------------------------------
  ### Copy default INPUT files to rundir
  if (!is.null(settings$model$config) && dir.exists(settings$model$config)) {
    input.dir <- settings$model$config
  } else {
    input.dir <- file.path(system.file(package = "PEcAn.SDGVM"),"inputs")
  }
  system2("cp", args = paste0(input.dir, "/* ", rundir))
  
  ## read general inputs
  input.file <- file.path(rundir,"input.dat")
  input.text <- readLines(con = input.file, n = -1)
  
  ## (2) MET PATH
  input.text[2] <- dirname(settings$run$inputs$met$path)
  
  ## (4) GLOBAL SOIL
  input.text[4] <- dirname(settings$run$inputs$soil$path)
  
  ## (5) GLOBAL LAND USE
  input.text[5] <- dirname(settings$run$inputs$lulc$path)
  
  ## (6 & 7) CO2 PATH & CO2 SWITCH
  if(is.null(settings$run$inputs$co2)){ ## no CO2 provided
    input.text[7] <- as.character(abs(as.numeric(input.text[7]))) ## default conc
  } else { ## CO2 file provided
    input.text[6] <- settings$run$inputs$co2$path  ## includes filename
    input.text[7] <- as.character(-abs(as.numeric(input.text[7]))) ## negative flag
  }

  ## (8) DEFAULT DATA DIR
  input.text[8] <- dirname(settings$run$inputs$defaults$path)
  
  ## (9) model version
  
  ## (10, 11) switches
  
  ## (13) RESTART READ DIR
  input.text[13] <- outdir
  
  ## (14) OUTPUT WRITE DIR
  input.text[14] <- outdir
  
  ## (16) RUN TYPE
  ## Spin-up:		1 1 0
  ## Run:       0 0 0
  ## Spin + run: 1 0 0
  input.text[16] <- "1 0 0"
  
  ## (17) START/END YEAR
  input.text[17] <- paste(lubridate::year(as.Date(settings$run$start.date)),
                          lubridate::year(as.Date(settings$run$end.date)))
  
  ## parse spaces to find blocks, fill in from bottom-up
  breaks <- which(trimws(input.text) == "")
  
  ## LAND COVER TRANSITIONS
  nb <- dplyr::last(breaks)
  # for now, leaving the defaults
  
  ## SOIL
  # for now, leaving the defaults
  nb <- breaks[length(breaks)-1]
  
  ## LAT/LON
  nb <- breaks[length(breaks)-2]
  input.text[nb+1] <- paste(settings$run$site$lat,settings$run$site$lon)
  
  ## FIRE
  # for now, leaving the defaults
  nb <- breaks[length(breaks)-3]
  
  ## LAND COVER CLASSES
  # needs to be more sophisticated, but start with equal area
  nb <- breaks[length(breaks)-4]
  npft <- length(trait.values) - 1
  pft.names <- names(trait.values)
  input.text[nb+1] <- paste("1",paste(rep(100/npft,npft),pft.names,collapse = " "))
  
  ## PFTS
  parm.names <- c("name","mixing_parameter","photosynthetic_pathway","FOO","FOO","FOO",
                  "AGEMX","wood_density","FOO","FOO","SLA","leaf_turnover_rate",
                  "wood_turnover_rate","root_turnover_rate","FOO","FOO","FOO","FOO",
                  "FOO","FOO","FOO","FOO","FOO","FOO","FOO","FOO",
                  "FOO","clumping_factor","Vcmax","FOO","Jmax","FOO",
                  "cuticular_cond","stomatal_slope")
  nb <- breaks[length(breaks)-5:4] + c(1,-1)
  default.pfts <- read.table(textConnection(input.text[nb[1]:nb[2]]),header=FALSE)
  colnames(default.pfts) <- parm.names
  updated.pfts <- default.pfts
  pft.lut <- read.table(system.file("PFTs.txt",package = "PEcAn.SDGVM"),header=TRUE,stringsAsFactors = FALSE) ## PFT Look up table to map BETY pfts to defaults
  for(i in 1:npft){
    if(pft.names[i] %in% default.pfts$name){
      ## Update PFTs that are among default list
      j = which(default.pfts$name == pft.names[i])
    } else { ## Add new PFTs
      
      ## which PFT is the 'template'
      if(pft.names[i] %in% pft.lut$PFT){
        j = which(pft.names[i] == pft.lut$PFT)
        j = which(default.pfts$name == pft.lut$Default[j])
      } else {
        ## if match not found, set to a default
        j = which(default.pfts$name == "Dc_Bl")
      }
      
      ## clone defaults to a new row
      updated.pfts <- rbind(updated.pfts,default.pfts[j,])
      j = nrow(updated.pfts)
      updated.pfts$name[j] = pft.names[i]
      
    }
    
    ## update individual parameter values
    pft <- trait.values[[i]]
    for(k in seq_along(pft)){
      
      ## photosynthetic_pathway
      # translate 3/4/5 for C3/C4/CAM to 0=C4, 1=C3
    
      ## SLA m2/kg -> m2/g DW
      updated.pfts$SLA[j] = udunits2::ud.convert(pft$SLA,"m2/kg","m2/g") 
    
      ## leaf_turnover_rate (yr-1) -> leaf life span (days)
    
      ## wood_turnover_rate (yr-1) -> stem life span (days)
    
      ## root_turnover_rate (yr-1) -> root life span (days)
    
      ## Convert Jmax to Jmax/Vcmax
    
      ## stomatal slope and intercept units? (Bety g0 = umol H2O m-2 s-1)
    
    } ## END loop over traits
    
  }
  
  ## write updated PFTs
  
  ## write new PFTs
  #input.text <- append(input.text,PFT.string,after = 28)
  
  ## write general input file
  writeLines(input.text, con = input.file)  
  
  ###########################################
  ## write additional parameters in param.dat
  param.file <- file.path(rundir,"param.dat")
  param.text <- readLines(con = param.file, n = -1)
  
  writeLines(param.text, con = param.file)  
  
  ## write additional parameters in misc_params.dat
  
} # write.config.SDGVM
