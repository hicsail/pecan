#!/usr/bin/Rscript --vanilla
# Modified from Alexis Helgeson's assim.sequential/inst/restart_SDAworkflow_scripts/SDA_Workflow_NA.R
#
# ----------------------------------------------------------------------
#------------------------------------------ Load required libraries-----
# ----------------------------------------------------------------------
.libPaths(c("/projectnb/dietzelab/dietze/test-pecan/R/library",.libPaths()))
library("PEcAn.all")
library("PEcAn.settings")
library("PEcAn.utils")
library("PEcAn.data.remote")
library("PEcAnAssimSequential")
library("REddyProc")
library("tidyverse")
library("furrr")
library("R.utils")
library("dynutils")
library('nimble')
library("sp")
library("sf")
library("lubridate")
#plan(multisession)
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

# ----------------------------------------------------------------------------------------------
#------------------------------------------Prepared SDA Settings -----
# ----------------------------------------------------------------------------------------------
## parse start date
option_list = list(optparse::make_option("--start.date",
                                         default = Sys.Date(),
                                         type="character"),
                   optparse::make_option("--prev",
                                         type="character"),
                   optparse::make_option("--settings",
                                         type="character")
                   )
args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
#args$start.date = "2022-05-18 00:00:00"
#args$prev = "/projectnb/dietzelab/dietze/hf_landscape_SDA/test02/FOF2022-05-17/"

## set dates
start.date = lubridate::as_date(args$start.date)
prev.date  = basename(args$prev) |> substr(start=4,stop=13) |> lubridate::as_date()
end.date <- start.date + lubridate::days(35)

#------------------------------------------------------------------------------------------------
#------------------------------------------ Preparing the pecan xml -----------------------------
#------------------------------------------------------------------------------------------------
restart <- list()
restart$filepath <- args$prev
set = readRDS(args$settings)
projectdir = set$outdir

# --------------------------------------------------------------------------------------------------
#---------------------------------------------- NA DATA -------------------------------------
# --------------------------------------------------------------------------------------------------

#initialize obs.mean/cov NAs & dates
# note:
#   for t=1: SDA doesn’t use metSplit, uses met.start, met.end
#            for iterative forecast restart, this is PREV
#   for t>1: metSplit uses obs.times t-1 -> t
site.ids <- PEcAn.settings::papply(set,function(x)(x$run$site$id)) %>% unlist() %>% as.character()
nsite = length(site.ids)

NAdata = data.frame(date = c(rep(format(start.date, "%Y-%m-%d %H:%M:%S", tz = "GMT"),nsite),
                             rep(format(end.date, "%Y-%m-%d %H:%M:%S", tz = "GMT"),nsite)), ##coerce to datetimes to avoid pumpkin rule
                    site_id = rep(site.ids,times=2),
                    data = rep(NA,nsite*2))
obs.mean <- obs.cov <- split(NAdata, NAdata$date)
date.obs <- names(obs.mean)
obs.mean <- purrr::map(
  names(obs.mean),
  function(namesl){
    split(
      obs.mean[[namesl]],
      obs.mean[[namesl]]$site_id) %>%
      purrr::map(
        ~.x[3] %>%
          stats::setNames(c("LAI")) %>%
          `row.names<-`(NULL))
  }
) %>% stats::setNames(date.obs)

obs.cov <- purrr::map(
  names(obs.cov),
  function(namesl){
    purrr::map(
      split(
        obs.cov[[namesl]],
        obs.cov[[namesl]]$site_id),
      ~.x[3]^2 %>%
        unlist %>%
        diag(nrow = 1, ncol = 1))
  }
) %>% stats::setNames(date.obs)

#add start.cut to restart list
restart$start.cut <- as_datetime(start.date)
restart$start.cut <- format(restart$start.cut, "%Y-%m-%d %H:%M:%S", tz = "GMT")


#-----------------------------------------------------------------------------------------------
#------------------------------------------ Fixing the settings --------------------------------
#-----------------------------------------------------------------------------------------------
#Using the found dates to run - this will help to download mets
for(s in seq_along(set)){
  set[[s]]$run$start.date = start.date
  set[[s]]$run$end.date   = end.date
  set[[s]]$run$site$met.start = prev.date  ## time period for t == 1, prev forecast
  set[[s]]$run$site$met.end   = start.date   
}

# Setting dates in assimilation tags - This will help with preprocess split in SDA code
set$state.data.assimilation$start.date <-as.character(start.date)
set$state.data.assimilation$end.date <-as.character(max(names(obs.mean)))

# --------------------------------------------------------------------------------------------------
#---------------------------------------------- PEcAn Workflow -------------------------------------
# --------------------------------------------------------------------------------------------------
#info
set$info$date <- paste0(format(Sys.time(), "%Y/%m/%d %H:%M:%S"))
next.oldir <- paste0(format(Sys.time(), "%Y-%m-%d-%H-%M"))
#Update/fix/check settings. Will only run the first time it's called, unless force=TRUE
#set <- PEcAn.settings::prepare.settings(set, force = TRUE)
## TODO: make sure settings are prepared; for remote, make sure to set host directories

## outdirs
set$outdir = file.path(set$outdir,paste0("FNA",start.date,"/"))
set$rundir = file.path(set$outdir,"run")
set$modeloutdir = file.path(set$outdir,"out")
set$pfts$pft$outdir = file.path(set$outdir,"pft")
set$host$rundir <- set$rundir
set$host$outdir <- set$modeloutdir
set$host$folder <- set$modeloutdir
dir.create(set$outdir,showWarnings = FALSE)
dir.create(set$rundir,showWarnings = FALSE)
dir.create(set$modeloutdir,showWarnings = FALSE)
dir.create(set$pfts$pft$outdir,showWarnings = FALSE)

#manually add in clim files 
path = "/projectnb/dietzelab/ahelgeso/NOAA_met_data_CH1/noaa_clim/HARV/" ## hack
met_paths <- list.files(path = file.path(path, start.date), full.names = TRUE, pattern = ".clim")
#met_paths <- list.files(path = file.path(settings$run$settings.1$inputs$met$path, start.date), full.names = TRUE, pattern = ".clim")
if(purrr::is_empty(met_paths)){
  print(paste("SKIPPING: NO MET FOR",start.date))
##  cat(as.character(start.date),sep="\n",file=file.path(dirname(set$outdir),"NO_MET"),append=TRUE) ## add to list of dates missing met
## TODO: for forecast, don't update NO_MET as files may still arrive; only update for reanalysis
  stop_quietly()
}
met_paths = as.list(met_paths)
names(met_paths) = rep("path",nsite)
for(s in seq_along(set)){
  set[[s]]$run$inputs$met$source = "GEFS" 
  set[[s]]$run$inputs$met$path = met_paths
}

#add run ids from previous sda to settings object to be passed to build X
prev_run_ids = list.files(file.path(restart$filepath, "out"))
run_id = as.data.frame(strsplit(prev_run_ids,"-")) %>% t()
colnames(run_id) <- c("pre","ens","site")
rownames(run_id) <- NULL
run_id = as.data.frame(run_id) %>% 
  dplyr::mutate(folder=prev_run_ids,id = paste0("id",.data$ens)) %>% 
  dplyr::group_by(site)
###settings$runs$id = run_id
for(s in seq_along(set)){
  site_run_id = run_id |> dplyr::filter(site == as.list(set[[s]]$run$site$id)[[1]])
  set[[s]]$run$id =  as.list(site_run_id$folder)
  names(set[[s]]$run$id) = site_run_id$id
}

## job.sh
if(is.null(set$model$jobtemplate)) set$model$jobtemplate = file.path(projectdir,"template.job")

#save restart object
save(restart, next.oldir, args, file = file.path(set$outdir, "restart.Rdata"))
#run sda function
sda.enkf.multisite(settings = set, 
                   obs.mean = obs.mean, 
                   obs.cov = obs.cov, 
                   Q = NULL, 
                   restart = restart, 
                   forceRun = TRUE, 
                   keepNC = TRUE, 
                   run_parallel = FALSE,
                   parallel_qsub = FALSE,
                   control = list(trace = TRUE,
                                  FF = FALSE,
                                  interactivePlot = FALSE,
                                  TimeseriesPlot = FALSE,
                                  BiasPlot = FALSE,
                                  plot.title = NULL,
                                  facet.plots = FALSE,
                                  debug = FALSE,
                                  pause = FALSE,
                                  Profiling = FALSE,
                                  OutlierDetection=FALSE,
                                  free_run = FALSE))  ## seems to be defined twice




