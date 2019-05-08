# #############################################################################
# ExtremeEvents.r
#
# Authors: Björn Brötz (DLR, Germany)
#          Marit Sandstad (CICERO, Norway)
#          Christian W. Mohr (CICERO, Norway)
# #############################################################################
# Description
#    Calculate extreme events with plotting functionality
#
# Modification history
#    2019 0506-hard_jo  : conversion to ESMValTool2
#    2018 1006-A_cwmohr : observation read and sorting fixes
#    2018 1003-A_cwmohr : correcting r.interface output for observation data.
#    2018 0725-A_cwmohr : modification of timeseries_main() and climdex selection
#    2018 0615-A_cwmohr : more clean up of code
#    2018 0131-A_laue_ax: clean-up of code, adaptation to ESMValTool standards
#                        added tagging support
#    2017 0920-A_sand_ma: modification to include plotting
#    2016 0414-A_broe_bj: written
# ############################################################################

library(climdex.pcic.ncdf)  # nolint
library(tools)
library(yaml)
library(ncdf4)
library(ncdf4.helpers)

cdo <- function(command, args = "", input = "", options = "", output = "",
                stdout = "", noout = F) {
  if (args != "") args <- paste0(",", args)
  if (stdout != "") {
    stdout <- paste0(" > '", stdout, "'")
    noout <- T
  }
  if (input[1] != "") {
    for (i in 1:length(input)) {
      input[i] <- paste0("'", input[i], "'")
    }
    input <- paste(input, collapse = " ")
  }
  output0 <- output
  if (output != "") {
    output <- paste0("'", output, "'")
  } else if (!noout) {
    output <- tempfile()
    output0 <- output
  }
  argstr <- paste0(
    options, " ", command, args, " ", input, " ", output,
    " ", stdout
  )
  print(paste("cdo", argstr))
  ret <- system2("cdo", args = argstr)
  if (ret != 0) {
    stop(paste("Failed (", ret, "): cdo", argstr))
  }
  return(output0)
}

nco <- function(cmd, argstr) {
  ret <- system2(cmd, args = argstr)
  if (ret != 0) {
    stop(paste("Failed (", ret, "): ", cmd, " ", argstr))
  }
}

provenance_record <- function(infile) {
  xprov <- list(
    ancestors = infile,
    authors = list("broe_bj", "sand_ma", "mohr_cw", "hard_jo"),
    references = list("zhang-2011"),
    projects = list("crescendo", "c3s-magic"),
    caption = "Extreme events indices",
    statistics = list("other"),
    realms = list("atmos"),
    themes = list("phys"),
    domains = list("global")
  )
  return(xprov)
}

diag_scripts_dir <- Sys.getenv("diag_scripts")
source(paste0(diag_scripts_dir, "/extreme_events/cfg_climdex.r"))  # nolint
source(paste0(diag_scripts_dir, "/extreme_events/cfg_extreme.r"))  # nolint
source(paste0(diag_scripts_dir,
         "/extreme_events/common_climdex_preprocessing_for_plots.r"))  # nolint
source(paste0(diag_scripts_dir,
         "/extreme_events/make_timeseries_plot.r"))  # nolint
source(paste0(diag_scripts_dir,
         "/extreme_events/make_Glecker_plot2.r"))  # nolint

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
for (myname in names(settings)) {
  temp <- get(myname, settings)
  assign(myname, temp)
}

list0 <- yaml::read_yaml(settings$input_files[1])
# extract metadata
models_name <- unname(sapply(list0, "[[", "dataset"))
models_ensemble <- unname(sapply(list0, "[[", "ensemble"))
models_start_year <- unname(sapply(list0, "[[", "start_year"))
models_end_year <- unname(sapply(list0, "[[", "end_year"))
models_experiment <- unname(sapply(list0, "[[", "exp"))
models_project <- unname(sapply(list0, "[[", "project"))
diag_base <- unname(sapply(list0, "[[", "diagnostic"))[1]
#### Correction r.interface output correction ####
models_experiment[models_experiment == "No_value"] <- "No-values"

variables <- c()
climofiles <- c()
models <- c()
metadata <- c()

# loop over variables
for (i in 1:length(settings$input_files)) {
  metadata <- yaml::read_yaml(settings$input_files[i])
  models_name <- unname(sapply(metadata, "[[", "dataset"))
  short_name <- unname(sapply(metadata, "[[", "short_name"))
  variables <- c(variables, short_name)
  models <- c(models, models_name)
  climofiles <- c(climofiles, names(metadata))
}

field_type0 <- "T2Ds"
# associated to first climofile
print(paste(diag_base, ": starting routine"))

# create working dirs if they do not exist
work_dir <- settings$work_dir
regridding_dir <- settings$run_dir
plot_dir <- settings$plot_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)
dir.create(plot_dir, recursive = T, showWarnings = F)

# setup provenance file and list
provenance_file <- paste0(regridding_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

if (anyNA(base_range)) {
  stop("Please choose a base_range!")
}
model_range <- c(
  max(strtoi(models_start_year)),
  min(strtoi(models_end_year))
)
if ( (base_range[1] < max(strtoi(models_start_year))) |
  (base_range[2] > min(strtoi(models_end_year)))) {
  stop(paste(
    "Base range", base_range[1], "-", base_range[2],
    "outside available model data period",
    model_range[1], "-", model_range[2]
  ))
}
print(paste("Base range:", base_range[1], "-", base_range[2]))

if (anyNA(regrid_dataset)) {
  regrid_dataset <- reference_datasets[1]
  print(paste(
    "Regrid dataset not set, choosing first reference dataset:",
    regrid_dataset
  ))
}

## Find earlier climdex indices in work folder
climdex_files <- list.files(path = work_dir, pattern = "ETCCDI")

# Fix input files removing bounds
print("Removing bounds from preprocessed files")
for (i in 1:length(climofiles)) {
  tmp <- tempfile()
  nco("ncks", paste(
    "-C -O -x -v lat_bnds,lon_bnds,time_bnds",
    climofiles[i], tmp
  ))
  nco("ncatted", paste("-O -a bounds,time,d,,", tmp))
  nco("ncatted", paste("-O -a bounds,lat,d,,", tmp))
  nco("ncatted", paste("-O -a bounds,lon,d,,", tmp))
  nco("ncatted", paste0("-O -a coordinates,", variables[i], ",d,, ", tmp))
  file.copy(tmp, climofiles[i], overwrite = TRUE)
  unlink(tmp)
}

##
## At this stage climdex indices are calculated. This process is extremely tedious and check points are in place to check whether the indicese are already produced.
## If the climdex files are there, then this process is skipped. Delete the climdex files from the work folder if you wish to have the climdex indices recalculated.
##
for (model_idx in c(1:length(models_name))) {
  author.data <- list(institution = "None", institution_id = "None")
  template <- paste("var_timeres_", models_name[model_idx], "_",
    models_experiment[model_idx], "_",
    models_ensemble[model_idx], "_",
    models_start_year[model_idx],
    "01-", models_end_year[model_idx],
    "12.nc",
    sep = "", collapse = ""
  )
  print("")
  print(paste0(">>>>>>>> Template name: ", template))
  print("")

  idx_select <- unique(c(timeseries_idx, gleckler_idx))

  ## Check point for existing files
  climdex_file_check <- paste0(
    idx_select, "_",
    models_name[model_idx], "_",
    models_experiment[model_idx], "_",
    models_ensemble[model_idx], "_",
    models_start_year[model_idx], "-",
    models_end_year[model_idx]
  )
  check_control <- vector("logical", length(climdex_file_check))
  n <- 0
  for (chck in climdex_file_check) {
    n <- n + 1
    tmp <- length(grep(chck, climdex_files))
    check_control[n] <- (tmp > 0)
  }
  print(check_control)

  if (!all(check_control)) {
    print("")
    print(paste0(">>>>>>> Producing Indices for ", models_name[model_idx]))
    print(climofiles[models == models_name[model_idx]])
    print("")
    infiles <- climofiles[models == models_name[model_idx]]
    indices <- sub("ETCCDI.*", "", idx_select)
    create.indices.from.files(infiles,  # nolint
      work_dir, template, author.data,
      base.range = base_range,
      parallel = climdex_parallel,
      verbose = TRUE, max.vals.millions = 20,
      climdex.vars.subset = indices
    )

    # Set provenance for output files
    # Get new list of files after computation
    infiles <- climofiles[models == models_name[model_idx]]
    xprov <- provenance_record(list(infiles))
    climdex_files <- list.files(
      path = work_dir, pattern = "ETCCDI",
      full.names = TRUE
    )
    for (chck in climdex_file_check) {
      fname <- grep(chck, climdex_files, value = TRUE)
      provenance[[fname]] <- xprov
    }
  }
}

if (write_plots) {
  #############################
  # A climdex processing section is needed here for observation data.
  # CMORized observation data found in the obs directory,
  # has it's climdex indices calculated,
  # which are then placed in the work/ExtremeEvents directory
  #############################

  ## Splitting models from observations

  ###################################
  #### Produce time series plots ####
  ###################################

  if (anyNA(analysis_range)) {
    analysis_range[1] <- max(strtoi(models_start_year))
    analysis_range[2] <- min(strtoi(models_end_year))
    print(paste(
      "Analysis range not defined, assigning model range:",
      analysis_range[1], "-", analysis_range[2]
    ))
  }
  if ( (analysis_range[1] < max(strtoi(models_start_year))) |
    (analysis_range[2] > min(strtoi(models_end_year)))) {
    stop(paste(
      "Analysis range", analysis_range[1], "-", analysis_range[2],
      "outside available model data period",
      model_range[1], "-", model_range[2]
    ))
  }
  print(paste("Analysis range:", analysis_range[1], "-", analysis_range[2]))

  # These are forced here for testing

  print("------ Model datasets ------")
  print(setdiff(models_name, reference_datasets))
  print("---- Reference datasets ----")
  print(reference_datasets)
  print("----------------------------")
  if (ts_plt) {
    print("")
    print(paste0(">>>>>>>> TIME SERIE PROCESSING INITIATION"))
    plotfiles <- timeseries_main(
      path = work_dir, idx_list = timeseries_idx,
      model_list = setdiff(models_name, reference_datasets),
      obs_list = reference_datasets, plot_dir = plot_dir,
      normalize = normalize,
      start_yr = analysis_range[1], end_yr = analysis_range[2]
    )
    xprov <- provenance_record(climofiles)
    for (fname in plotfiles) {
      provenance[[fname]] <- xprov
    }
  }

  ###############################
  #### Produce Gleckler plot ####
  ###############################
  if (glc_plt) {
    print("")
    print(paste0(">>>>>>>> GLECKLER PROCESSING INITIATION"))

    ## Check if Gleckler Array already exists
    nidx <- length(gleckler_idx) # number of indices
    nmodel <- length(models_name) # number of models
    nobs <- length(reference_datasets) # number of observations
    arrayname <- paste0(
      "Gleckler-Array_", nidx, "-idx_",
      nmodel, "-models_", nobs, "-obs", ".RDS"
    )
    arraydirname <- paste0(plot_dir, "/", diag_base, "/", arrayname)
    if (glc_arr) {
      if (file.exists(arraydirname)) {
        file.remove(arraydirname)
      }
      promptinput <- "y"
    }

    if (file.exists(arraydirname)) {
      promptinput <- "n"
    } else {
      promptinput <- "y"
    }

    #### Running gleckler_main ####
    plotfiles <- gleckler_main(
      path = work_dir, idx_list = gleckler_idx,
      model_list = setdiff(models_name, reference_datasets),
      obs_list = reference_datasets,
      plot_dir = plot_dir, promptinput = promptinput,
      start_yr = analysis_range[1], end_yr = analysis_range[2]
    )

    xprov <- provenance_record(list(climofiles))
    for (fname in plotfiles) {
      provenance[[fname]] <- xprov
    }
  }
}

# Write provenance to file
write_yaml(provenance, provenance_file)