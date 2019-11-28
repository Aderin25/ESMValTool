library(yaml)
library(ncdf4)
library(SPEI)
library(RColorBrewer) # nolint

getnc <- function(yml, m, lat = FALSE) {
  id <- nc_open(yml[m][[1]]$filename, readunlim = FALSE)
  if (lat){
    v <- ncvar_get(id, "lat")
  }else{
    v <- ncvar_get(id, yml[m][[1]]$short_name)
  }
  nc_close(id)
  return(v)
}

ncwritespi <- function(yml, m, data, wdir){
  fnam <- strsplit(yml[m][[1]]$filename, "/")[[1]]
  pcs <- strsplit(fnam[length(fnam)], "_")[[1]]
  pcs[which(pcs == yml[m][[1]]$short_name)] <- "spi"
  onam <- paste0(wdir, "/", paste(pcs, collapse = "_"))
  ncid_in <- nc_open(yml[m][[1]]$filename)
  # var <- ncid_in$var[[yml[m][[1]]$short_name]]
  xdim <- ncid_in$dim[["lon"]]
  ydim <- ncid_in$dim[["lat"]]
  tdim <- ncid_in$dim[["time"]]
  allatt <- ncatt_get(ncid_in, "pr")
  print("allatt")
  print(allatt)
  fillvalue <- ncatt_get(ncid_in,"pr","_FillValue")
  print("fillvalue")
  print(fillvalue)
  globat <- ncatt_get(ncid_in, 0)
  print("globat")
  print(globat)
  # hdim <- ncdim_def("bins", "level", bins[1:(length(bins) - 1)])
  # hdim2 <- ncdim_def("binsup", "level", bins[2:length(bins)])
  fillfloat <- 1.e+20
  as.single(fillfloat)
  var_spi <- ncvar_def("spi", "1", list(xdim, ydim, tdim), fillfloat)
  idw <- nc_create(onam, var_spi)
  ncvar_put(idw, "spi", data)
  # ncatt_put(idw, "spi", "_FillValue", fillfloat)
  cntatt <- 1
  for (thisattname in names(globat)){
    print("thisattname")
    print(thisattname)
    print("globat[[cntatt]]")
    print(globat[[cntatt]])
    ncatt_put(idw, 0, thisattname, globat[[cntatt]])
    cntatt <- cntatt + 1
  }
  nc_close(idw)
  nc_close(ncid_in)
}

whfcn <- function(x, ilow, ihigh){
  return(length(which(x >= ilow & x < ihigh)))
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
wdir <- params$work_dir
dir.create(wdir, recursive = TRUE)
pdir <- params$plot_dir
dir.create(pdir, recursive = TRUE)
var1_input <- read_yaml(params$input_files[1])
print("var1_input")
print(var1_input)
nmods <- length(names(var1_input))

fillfloat <- 1.e+20
as.single(fillfloat)

histbrks <- c(-99999, -2, -1.5, -1, 1, 1.5, 2, 99999)
histnams <- c("Extremely dry", "Moderately dry", "Dry",
              "Neutral",
              "Wet", "Moderately wet", "Extremely wet")
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
refmsk[refmsk > 10000] <- fillfloat
refmsk[!is.na(refmsk)] <- 1

# histarr <- array(1.e+20, c(nmods, length(histnams)))
for (mod in 1:nmods){
   v1 <- getnc(var1_input, mod)
   print("var1_input[mod][[1]]$cmor_table")
   print(var1_input[mod][[1]]$cmor_table)
   d <- dim(v1)
   v1_spi <- array(fillfloat, dim=d)
   for (i in 1:d[1]){
     wh <- which(!is.na(refmsk[i,]))
     if (length(wh) > 0){
       tmp <- v1[i,wh,]
       v1_spi[i,wh,] <- t(spi(t(tmp), 6, na.rm = TRUE,
                        distribution = "Gamma")$fitted)
     }
   }
   v1_spi[is.infinite(v1_spi)] <- fillfloat
   v1_spi[is.na(v1_spi)] <- fillfloat
   v1_spi[v1_spi > 10000] <- fillfloat
   ncwritespi(var1_input, mod, v1_spi, wdir)
   # hist_spi <- array(1.e+20, c(d[1], d[2], length(histbrks) - 1))
   # for (nnh in 1:(length(histbrks) - 1)){
     # hist_spi[,,nnh] <- apply(v1_spi, c(1, 2), FUN = whfcn,
     #                          ilow = histbrks[nnh],
     #                          ihigh = histbrks[nnh + 1])
   # }
   # ncwritenew(var1_input, mod, hist_spi, wdir, histbrks)
   # Weight against latitude
   # h <- c(1:length(histnams)) * 0
   # for (j in 1:d[2]){
     # h <- h + hist(v1_spi[j,,], breaks = histbrks,
     #               plot = FALSE)$counts * cos(lat[j] * pi / 180.)
   # }
   # histarr[mod, ] <- h / sum(h, na.rm = TRUE)
}
# save(histarr, file = paste0(params$work_dir,
#                             "/", "histarr.rsav"))

# bhistarr <- array(1.e+20, c(nmods - 1, 7))
# marr <- c(1:nmods)[c(1:nmods) != nref]
# cnt <- 1
# for (m in marr){
#   bhistarr[cnt, ] <- histarr[m, ] - histarr[nref, ]
#   cnt <- cnt + 1
# }
# parr <- c(nref, marr)
# 
# mnam <- c(1:nmods) * 1.e+20
# for (m in 1:nmods) mnam[m] <- var1_input[m][[1]]$dataset
# 
# qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ] # nolint
# col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, # nolint
#                            rownames(qual_col_pals)))
# cols <- c("black", sample(col_vector, nmods - 1))
# 
# png(paste0(params$plot_dir, "/", "histplot.png"),
#     width = 1000, height = 500)
#  par(mfrow = c(2, 1), oma = c(3, 3, 3, 13), mar = c(2, 1, 1, 1))
#  barplot(histarr[parr, ], beside = 1, names.arg = histnams,
#          col = cols, xaxs = "i")
#  box()
#  mtext("Probability", side = 2, line = 2.1)
#  barplot(bhistarr, beside = 1, names.arg = histnams,
#          col = cols[2:nmods], xaxs = "i")
#  box()
#  mtext("Absolute difference", side = 2, line = 2.1)
#  mtext("Standardized precipitation index", outer = TRUE,
#        cex = 2, font = 2)
#  par(fig = c(0.8, .95, 0.1, 0.9), new = T, oma = c(0, 0, 0, 0),
#      mar = c(0, 0, 0, 0))
#  legend("topright", mnam[parr], fill = cols)
# dev.off()