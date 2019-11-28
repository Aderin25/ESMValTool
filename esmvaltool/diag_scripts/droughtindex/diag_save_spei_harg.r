library(yaml)
library(ncdf4)
library(SPEI)
library(RColorBrewer) # nolint

leap_year <- function(year) {
  return(ifelse( (year %% 4 == 0 & year %% 100 != 0) |
                  year %% 400 == 0, TRUE, FALSE))
}

getnc <- function(yml, m, lat = FALSE) {
  print("yml[m][[1]]$filename")
  print(yml[m][[1]]$filename)
  id <- nc_open(yml[m][[1]]$filename, readunlim = FALSE)
  if (lat){
     v <- ncvar_get(id, "lat")
  }else{
   v <- ncvar_get(id, yml[m][[1]]$short_name)
   if (yml[m][[1]]$short_name == "tas") v <- v - 273.15
   if (yml[m][[1]]$short_name == "tas") print(v)
   if (yml[m][[1]]$short_name == "pr"){
     time <- ncvar_get(id, "time")
     tcal <- ncatt_get(id, "time", attname = "calendar")
     tunits <- ncatt_get(id, "time", attname = "units")
     tustr <- strsplit(tunits$value, " ")
     stdate <- as.Date(time[1], origin = unlist(tustr)[3])
     nddate <- as.Date(time[length(time)], origin = unlist(tustr)[3])
     if (tcal$value == "365_day"){
       # Correct for missing leap years in nddate
       diff <- as.numeric(nddate - stdate, units = "days")
       dcorr <- floor( (diff / 365 - diff / 365.25) * 365.25)
       nddate <- nddate + dcorr
     }
     if (tcal$value == "360_day"){
       v <- v * 30 * 24 * 3600.
     }else{
       cnt <- 1
       monarr <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
       date <- stdate
       while (date <= nddate){
         year <- as.numeric(substr(date, 1, 4))
         lpyear <- leap_year(year)
         month <- as.numeric(substr(date, 6, 7))
         mdays <- monarr[month]
         pdays <- mdays
         if (month == 2 & lpyear == TRUE){
           pdays <- 29
           if (tcal$value != "365_day"){
             mdays <- 29
           }else{
             mdays <- 28
           }
         }
         v[,,cnt] <- v[,,cnt] * mdays * 24 * 3600.
         date <- date + pdays
         cnt <- cnt + 1
       }
     }
   }
  }
  nc_close(id)
  return(v)
}

ncwritespei <- function(yml, m, data, wdir){
  fnam <- strsplit(yml[m][[1]]$filename, "/")[[1]]
  pcs <- strsplit(fnam[length(fnam)], "_")[[1]]
  pcs[which(pcs == yml[m][[1]]$short_name)] <- "spei"
  onam <- paste0(wdir, "/", paste(pcs, collapse = "_"))
  ncid_in <- nc_open(yml[m][[1]]$filename)
  xdim <- ncid_in$dim[["lon"]]
  ydim <- ncid_in$dim[["lat"]]
  tdim <- ncid_in$dim[["time"]]
  allatt <- ncatt_get(ncid_in, "pr")
  fillvalue <- ncatt_get(ncid_in,"pr","_FillValue")
  globat <- ncatt_get(ncid_in, 0)
  fillfloat <- 1.e+20
  as.single(fillfloat)
  var_spei <- ncvar_def("spei", "1", list(xdim, ydim, tdim), fillfloat)
  idw <- nc_create(onam, var_spei)
  ncvar_put(idw, "spei", data)
  cntatt <- 1
  for (thisattname in names(globat)){
    ncatt_put(idw, 0, thisattname, globat[[cntatt]])
    cntatt <- cntatt + 1
  }
  nc_close(idw)
  nc_close(ncid_in)
  return(onam)
}

whfcn <- function(x, ilow, ihigh){
  return(length(which(x >= ilow & x < ihigh)))
}


dohargreaves <- function(vmin, vmax, v1, lat){
 print("Estimating PET with Penman method.")
 dpet <- vmin * NA
 d <- dim(dpet)
 for (i in 1:d[2]){
  tmp1 <- vmin[,i,] - 273.15
  tmp2 <- vmax[,i,] - 273.15
  tmp3 <- v1[,i,]
  #tmp4 <- vr[,i,] * (86400.0 / 1e6) # W/(m2) to MJ/(m2 d)
  tmp6 <- hargreaves(t(tmp1), t(tmp2), lat = rep(lat[i], d[1]), Pre = t(tmp3), na.rm = TRUE)
  d2 <- dim(tmp6)
  tmp6 <- as.numeric(tmp6)
  dim(tmp6) <- d2
  dpet[,i,] <- t(tmp6)
 }
 return(dpet)
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
wdir <- params$work_dir
rundir <- params$run_dir
dir.create(wdir, recursive = TRUE)
pdir <- params$plot_dir
dir.create(pdir, recursive = TRUE)
nrvars <- length(params$input_files)
print("nrvars")
print(nrvars)


for (vnr in 1:nrvars){
  metadata <- read_yaml(params$input_files[vnr])
  modfile <- names(metadata)
  var_input <- read_yaml(params$input_files[vnr])
  if (var_input[1][[1]]$short_name == "pr") {
    var1_input <- var_input
    modfile1 <- names(metadata)
  }
  if (var_input[1][[1]]$short_name == "tasmin") {
    tasmin_input <- var_input
    modfiletasmin <- names(metadata)
  }
  if (var_input[1][[1]]$short_name == "tasmax") {
    tasmax_input <- var_input
    modfiletasmax <- names(metadata)
  }
#  if (var_input[1][[1]]$short_name == "rsds") {
#    rs_input <- var_input
#    rad_flag = TRUE
#    modfilers <- names(metadata)
#  }#
}

nmods <- length(names(var1_input))

fillfloat <- 1.e+20
as.single(fillfloat)

# setup provenance file and list
provenance_file <- paste0(rundir, "/", "diagnostic_provenance.yml")
provenance <- list()

# histbrks <- c(-99999, -2, -1.5, -1, 1, 1.5, 2, 99999)
# histnams <- c("Extremely dry", "Moderately dry", "Dry",
#               "Neutral",
#               "Wet", "Moderately wet", "Extremely wet")
refnam <- var1_input[1][[1]]$reference_dataset

n <- 1
while (n <= nmods){
  if (var1_input[n][[1]]$dataset == refnam) break
  n <- n + 1
}
nref <- n
lat <- getnc(var1_input, nref, lat = TRUE)
if (max(lat) > 90){
  print(paste0("Latitude must be [-90,90]: min=",
  min(lat), " max=", max(lat)))
  stop("Aborting!")
}

ref <- getnc(var1_input, nref, lat = FALSE)
refmsk <- apply(ref, c(1, 2), FUN = mean, na.rm = TRUE)

ref2 <- getnc(tasmin_input, nref, lat = FALSE)
refmsk2 <- apply(ref2, c(1, 2), FUN = mean, na.rm = TRUE)
refmsk2[refmsk2 > 1e29] <- NA

ref3 <- getnc(tasmax_input, nref, lat = FALSE)
refmsk3 <- apply(ref3, c(1, 2), FUN = mean, na.rm = TRUE)
refmsk3[refmsk3 > 1e29] <- NA

refmsk[is.na(refmsk2)] <- NA
refmsk[is.na(refmsk3)] <- NA
refmsk[refmsk > 1e29] <- NA
refmsk[!is.na(refmsk)] <- 1

xprov <- list(ancestors = list(""),
              authors = list("berg_pe"),
              references = list("vicente10jclim"),
              projects = list("c3s-magic"),
              caption = "",
              statistics = list("other"),
              realms = list("atmos"),
              themes = list("phys"),
              domains = list("global"))

# histarr <- array(NA, c(nmods, length(histnams)))
for (mod in 1:nmods){
   lat <- getnc(var1_input, mod, TRUE)
   v1 <- getnc(var1_input, mod, FALSE)
   v2 <- getnc(tasmin_input, mod, FALSE)
   v3 <- getnc(tasmax_input, mod, FALSE)
   pet <- dohargreaves(v2, v3, v1, lat)
   pme <- v1 - pet
   print(var1_input[mod][[1]]$cmor_table)
   d <- dim(pme)
   pme_spei <- pme * NA
   for (i in 1:d[1]){
     wh <- which(!is.na(refmsk[i,]))
     if (length(wh) > 1){
       tmp <- pme[i,wh,]
       pme_spei[i,wh,] <- t(spei(t(tmp), 6, na.rm = TRUE)$fitted)
     }
   }
   pme_spei[is.infinite(pme_spei)] <- fillfloat
   pme_spei[is.na(pme_spei)] <- fillfloat
   pme_spei[pme_spei > 10000] <- fillfloat
   filename <- ncwritespei(var1_input, mod, pme_spei, wdir)

   # pme_spei[is.infinite(pme_spei)] <- NA
   # pme_spei[pme_spei > 10000] <- NA
   # hist_spei <- array(NA, c(d[1], d[2], length(histbrks) - 1))
   # for (nnh in 1:(length(histbrks) - 1)){
   #   hist_spei[,,nnh] <- apply(pme_spei, c(1, 2), FUN = whfcn,
   #                            ilow = histbrks[nnh],
   #                            ihigh = histbrks[nnh + 1])
   # }
   # filename <- ncwritenew(var1_input, mod, hist_spei, wdir, histbrks)
   # Set provenance for output files
   xprov$caption <- "SPEI index per grid point."
   xprov$ancestors <- list(modfile1[mod], modfiletasmin[mod], modfiletasmax[mod])
   # provenance[[filename]] <- xprov
   for (t in 1:d[3]){
      tmp <- pme_spei[,,t]
      tmp[is.na(refmsk)] <- NA
      pme_spei[,,t] <- tmp
   }
   # pme_spei[is.infinite(pme_spei)] <- NA
   # pme_spei[pme_spei > 10000] <- NA
   # Weight against latitude
   # h <- c(1:length(histnams)) * 0
   # for (j in 1:d[2]){
   #   h <- h + hist(pme_spei[j,,], breaks = histbrks,
   #                 plot = FALSE)$counts * cos(lat[j] * pi / 180.)
   # }
   # histarr[mod,] <- h / sum(h, na.rm = TRUE)
}
# filehist <- paste0(params$work_dir, "/", "histarr.rsav")
# save(histarr, file = filehist)
# plot_file <- paste0(params$plot_dir, "/", "histplot.png")
# xprov$caption <- "Global latitude-weighted histogram of SPEI index."
# xprov$ancestors <- list(modfile1, modfile2)
# xprov[["plot_file"]] <- plot_file
provenance[[filename]] <- xprov
write_yaml(provenance, provenance_file)

# bhistarr <- array(NA, c(nmods - 1, 7))
# marr <- c(1:nmods)[c(1:nmods) != nref]
# cnt <- 1
# for (m in marr){
#   bhistarr[cnt,] <- histarr[m,] - histarr[nref,]
#   cnt <- cnt + 1
# }
# parr <- c(nref, marr)

# mnam <- c(1:nmods) * NA
# for (m in 1:nmods) mnam[m] <- var1_input[m][[1]]$dataset

# qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual",] # nolint
# col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, # nolint
#                            rownames(qual_col_pals)))
# cols <- c("black", sample(col_vector, nmods - 1))

# png(plot_file, width = 1000, height = 500)
#  par(mfrow = c(2, 1), oma = c(3, 3, 3, 13), mar = c(2, 1, 1, 1))
#  barplot(histarr[parr,], beside = 1, names.arg = histnams,
#          col = cols, xaxs = "i")
#  box()
#  mtext("Probability", side = 2, line = 2.1)
#  barplot(bhistarr, beside = 1, names.arg = histnams,
#          col = cols[2:nmods], xaxs = "i")
#  box()
#  mtext("Absolute difference", side = 2, line = 2.1)
#  mtext("Standardized precipitation-evapotranspiration index",
#         outer = TRUE, cex = 2, font = 2)
#  par(fig = c(0.8, .95, 0.1, 0.9), new = T, oma = c(1, 1, 1, 1) * 0,
#      mar = c(0, 0, 0, 0))
#  legend("topright", mnam[parr], fill = cols)
# dev.off()
