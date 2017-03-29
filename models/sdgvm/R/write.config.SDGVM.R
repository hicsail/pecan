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
##' @param trait.samples vector of samples for a given trait
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
  
  ## (4) GLOBAL SOIL
  
  ## (5) GLOBAL LAND USE
  
  ## (6) CO2 PATH
  
  ## (7) CO2 SWITCH
  
  ## (8) DEFAULT DATA DIR
  
  ## LAND COVER TRANSITIONS
  
  ## SOIL
  
  ## FIRE
  
  ## LAND COVER CLASSES
  
  ## PFTS
  #input.text <- append(input.text,PFT.string,after = 28)
  
  
  ## write general input file
  writeLines(input.text, con = input.file)  
  
} # write.config.SDGVM
