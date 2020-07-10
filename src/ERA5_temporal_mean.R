#!/usr/bin/env Rscript
################################################################################
# ERA5_temporal_mean.R
################################################################################
# George Allen, June 2020

# data originally downloaded by Alberto Borges:
# "Wind speed and precipitation
# I've extracted the ERA-5 Land monthly product 25°N-25°S for 2010-2019.
# The idea is to compute wind speed (square root of (v-Wind ^2 + u-Wind^2)) 
# and then average the whole data-set per station (to have a climatology) 
# and then extract wind and precipitation for each lake in Hydrolakes."


library(ncdf4)
library(raster)
library(zoo)
library(RColorBrewer)


################################################################################ 
# Define file paths: 
wd = "E:/research/2020_06_11_Borges_etal_tripical_lakes_CH4/git/Borges_etal_202x_tropic_lakes_CH4"
inDir = paste0(wd, "/in")
outDir = paste0(wd, "/out")
ncpath = paste0(inDir, "/ERA5/ERA5_u-wind_v-wind_precipitation_air_temp.nc")


################################################################################
# read in ERA5 netCDF (originally downloaded by Alberto Borges):
ncIn = nc_open(ncpath)
ncDim = ncIn$var$u10$varsize
lon = ncvar_get(ncIn,"longitude")
lat = ncvar_get(ncIn,"latitude")
lonRange = range(lon)
latRange = range(lat)

# generate blank output:
ncOut_u10 = ncOut_v10 = ncOut_t2m = ncOut_tp = matrix(NA, ncDim[2], ncDim[1])


#############################################################################
start_time = Sys.time()

for (i in 1:(ncDim[2]-1)){
  nc = ncvar_get(ncIn, "u10", 
                 start=c(1, i, 1), # starting index of netCDF to read in 
                 count=c(ncDim[1], i, ncDim[3])) # end
  
  ncOut_u10[i,] = apply(nc, 1, mean, na.omit=T)
  print(i)
}

Sys.time() - start_time 
image(t(ncOut_u10), col=rev(brewer.pal(10,"RdBu")))



for (i in 1:tr$n) {
  v <- getValuesBlock(s, row=tr$row[i], nrows=tr$nrows[i])
  mmean = apply(v, 1, mean, na.rm=T)
  bmean <- writeValues(bmean, mmean, tr$row[i])
}


#############################################################################
start_time = Sys.time()

for (i in 1:(ncDim[2]-1)){
  nc = ncvar_get(ncIn, "v10", 
                 start=c(1, i, 1), # starting index of netCDF to read in 
                 count=c(ncDim[1], i, ncDim[3])) # end
  
  ncOut_v10[i,] = apply(nc, 1, mean, na.omit=T)
  print(i)
}

Sys.time() - start_time 

image(t(ncOut_v10), col=rev(brewer.pal(10,"RdBu")))


#############################################################################
start_time = Sys.time()

for (i in 1:(ncDim[2]-1)){
  nc = ncvar_get(ncIn, "t2m", 
                 start=c(1, i, 1), # starting index of netCDF to read in 
                 count=c(ncDim[1], i, ncDim[3])) # end
  
  ncOut_t2m[i,] = apply(nc, 1, mean, na.omit=T)
  print(i)
}

Sys.time() - start_time 

image(t(ncOut_t2m), col=rev(brewer.pal(10,"RdBu")))




#############################################################################
start_time = Sys.time()

for (i in 1:(ncDim[2]-1)){
  nc = ncvar_get(ncIn, "tp", 
                 start=c(1, i, 1), # starting index of netCDF to read in 
                 count=c(ncDim[1], i, ncDim[3])) # end
  
  ncOut_tp[i,] = apply(nc, 1, mean, na.omit=T)
  print(i)
}

Sys.time() - start_time 

image(t(ncOut_tp), col=rev(brewer.pal(10,"RdBu")))


#############################################################################

ncExtent = ncIn

r <- raster(ncOut_t2m)
# replace with correct coordinates
extent(r) <- c(lonRange, latRange)
writeRaster(r, paste0(outDir, "/", 'temp01.tif'), overwrite=T)

r <- raster(
  ncOut_t2m,
  xmn = lonRange[1], xmx=lonRange[2],
  ymn = latRange[1], ymx=latRange[2], 
  crs = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
)




I = ncvar_get(ncIn, "u10")
I_slice <- I[1,1,]
lon <- ncvar_get(ncIn,"longitude")
lat <- ncvar_get(ncIn,"latitude")

# prep for display:
lat <- rev(lat)
I_slice <- I_slice[, ncol(I_slice):1]
image(lon, lat, I_slice, col=rev(brewer.pal(10,"RdBu")))





# get temporal mean of wind (n-s), wind (e-2), temperature, and precipitation. 
stack_u10 = stack(ncpath, varname="u10")



### i renamed your "asdatadates" object for simplicity
dates = as.Date(ncin$dim$time$vals/24, origin="1900-01-01")


pixel.ts = zoo(I, dates) 
out = as.numeric(aggregate(x=pixel.ts, , mean, na.rm=T))
out[is.nan(out)] = NA     

temporal_mean_stack <- function(x) {
  pixel.ts = zoo(x, dates) 
  out = as.numeric(aggregate(x=pixel.ts, rep(1, 120), fun=mean, na.rm=T))
  out[is.nan(out)] = NA     
  return(out)
}



means_stack = calc(stack_u10, temporal_mean_stack)


print(ncin)

# get longitude and latitude
lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"latitude")
nlat <- dim(lat)
head(lat)

print(c(nlon,nlat))

time <- ncvar_get(ncin,"time")
time

tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

####
# get temperature
tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)



################################################################################
# apply thresholds to calculate flow occurence for each GRADES segment:
require(ncdf4)
require(foreign)
require(sf)
require(ggplot2)
library(zoo)

# read in zero flow threshold table: 
# bin is range in aridity, LQ=0.25 quantile, MD=0.5, UQ=0.75:
thresholds = read.csv(paste0(outDir, "/tab/thresholds.csv"), header=T)

# calculate the following for each year in GRADES:
# % flow occurence 
# standard deviation of flow occurence
# coefficient of variability
# temporal trends <- too slow for now

# Join aridity from layers CSV to GRADES shapefiles based on COMID:
start_time = Sys.time() 
#layers = read.csv(paste0(inDir, "/geospatial_layers/data_global.csv"), header=T)
Sys.time() - start_time
# 25 secs


  start = 1
  end = ncIn$var$Q$varsize[2]
  
  j = 1
  interval  = 50 
  year = interval # assuming no leap years
  
  
  # read in river shapefile dbf file:
  dbf = foreign::read.dbf(mhPaths[i])
  
  # join aridity index to dbf based on COMID:
  m = match(dbf$COMID, layers$COMID)
  dbf$AI = layers$AI[m]
  # interpolate nas:
  dbf$AI = zoo::na.approx(dbf$AI)
  dbf$AI[dbf$AI<0] = 0
  
  # generate a zero flow threshold based on each AI value:
  #breaks = c(seq(0, 1.2, 0.2), 10) # this should be defined at top of code
  LQ = MQ = UQ = cut(dbf$AI, breaks)
  levels(LQ) = thresholds$LQ
  levels(MQ) = thresholds$MQ
  levels(UQ) = thresholds$UQ
  LQ = as.numeric(as.character(LQ))
  MQ = as.numeric(as.character(MQ))
  UQ = as.numeric(as.character(UQ))
  
  # create three tables, one for each flow threshold:
  dbf_LQ = data.frame(dbf)
  dbf_MQ = data.frame(dbf)
  dbf_UQ = data.frame(dbf)
  
  start_time = Sys.time()
  # read in GRADES netCDF in yearly chunks (to prevent memory overflow):
  while (start < end){
    
    print(paste("chunk", j))
    
    #print("reading in nc file...")
    
    # for last year: 
    if ((start+interval) > ncIn$var$Q$varsize[2]){
      interval = ncIn$var$Q$varsize[2] - start
    }
    
    # start_time = Sys.time()
    nc = ncvar_get(ncIn, "Q", 
                   start=c(1, start), # starting index of netCDF to read in 
                   count=c(ncIn$var$Q$varsize[1], interval)) # ending index of netCDF to read in 
    # Sys.time() - start_time
    #1.7 sec
    
    colName = paste0("ch_", substr(paste0("00", j), nchar(j), nchar(j)+2))
    
    print("Calculating flow permanence...")
    # a column giving flow persistence for each quartile:
    # start_time = Sys.time()
    
    
    # subtract no-flow threshold and count number of no flows:
    nc_LQ = sweep(nc, 1, LQ)
    nc_MQ = sweep(nc, 1, MQ)
    nc_UQ = sweep(nc, 1, UQ)
    
    # sum up the number of flows greater than the threshold: 
    # this could probably be sped up by applying a single threshold function over the entire 
    # table and then using row.sums() function. 
    dbf_LQ[, colName] = apply(nc_LQ, 1, function(x){sum(x > 0)}) / interval
    dbf_MQ[, colName] = apply(nc_MQ, 1, function(x){sum(x > 0)}) / interval
    dbf_UQ[, colName] = apply(nc_UQ, 1, function(x){sum(x > 0)}) / interval
    # Sys.time() - start_time
    # 9.2 sec
    
    
    
    start = start+interval
    j = j+1
    
  }
  
  
  Sys.time() - start_time
  
  # calculate mean, sd, and cov of flow occurence:
  chunkCols = grep("^ch", names(dbf_LQ))
  dbf_LQ$fPermLQ = rowMeans(dbf_LQ[ , chunkCols])
  dbf_MQ$fPermMQ = rowMeans(dbf_MQ[ , chunkCols])
  dbf_UQ$fPermUQ = rowMeans(dbf_UQ[ , chunkCols])
  
  # standard devations (of the chunks):
  dbf_LQ$fPermLQ_sd = apply(dbf_LQ[ , chunkCols], 1, sd)
  dbf_MQ$fPermMQ_sd = apply(dbf_MQ[ , chunkCols], 1, sd)
  dbf_UQ$fPermUQ_sd = apply(dbf_UQ[ , chunkCols], 1, sd)
  
  # coefficient of variation:
  dbf_LQ$fPermLQ_cov = dbf_LQ$fPermLQ_sd/dbf_LQ$fPermLQ
  dbf_MQ$fPermMQ_cov = dbf_MQ$fPermMQ_sd/dbf_MQ$fPermMQ
  dbf_UQ$fPermUQ_cov = dbf_UQ$fPermUQ_sd/dbf_UQ$fPermUQ
  
  
  
  # take least square trends:
  # chunkTab = dbf_LQ[, grep("^ch", names(dbf_LQ))]
  # regTab = apply(chunkTab, 1, 
  #       function(x){
  #         reg = lm(ind ~ x)
  #         return(data.frame(intercept = reg$coefficients[[1]], slope = reg$coefficients[[2]]))
  #       }
  #     ) 
  
  
  # take median estimated trends:
  # library(zyp)
  # theil_sen = function(x){
  #   sen = zyp.sen(x ~ ind)
  #   return(data.frame(intercept = sen$coefficients[[1]], slope = sen$coefficients[[2]]))
  # }
  # 
  # chunkTab = dbf_LQ[, grep("^ch", names(dbf_LQ))]
  # ind = 1:length(chunkCols)
  # senDF = apply(chunkTab, 1, theil_sen) 
  
  # clean up table:
  # dbf_LQ = dbf_LQ[, -grep("^fPerm", names(dbf_LQ))]
  # dbf_MQ = dbf_MQ[, -grep("^fPerm", names(dbf_MQ))]
  # dbf_UQ = dbf_UQ[, -grep("^fPerm", names(dbf_UQ))]
  
  
  # match COMIDs to sorted shapefile tables: 
  dbf_sorted = foreign::read.dbf(mhSortpaths[i])
  
  sortInd = match(dbf_sorted$COMID, dbf$COMID)
  
  outDbf_LQ = dbf_LQ[sortInd, grep("^fPerm", names(dbf_LQ))]
  outDbf_MQ = dbf_MQ[sortInd, grep("^fPerm", names(dbf_MQ))]
  outDbf_UQ = dbf_UQ[sortInd, grep("^fPerm", names(dbf_UQ))]
  
  # Calculate IQR:
  IQR = outDbf_LQ$fPermLQ - outDbf_UQ$fPermUQ
  
  dbfOut = cbind(dbf_sorted, 
                 AI=dbf$AI[sortInd], 
                 LQ_thresh=LQ[sortInd],  
                 MQ_thresh=MQ[sortInd],  
                 UQ_thresh=UQ[sortInd], 
                 outDbf_LQ, outDbf_MQ, outDbf_UQ, IQR)
  
  # write out shapefile attribute table:
  foreign::write.dbf(dbfOut, mhOutpaths[i])

Sys.time() - start_time


