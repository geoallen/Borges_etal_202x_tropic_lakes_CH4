#!/usr/bin/env Rscript
################################################################################
# Borges_etal_202x_tropicalLakesCH4emission_finalTableJoin.R
################################################################################
# George Allen, June-July 2020


# join data from Tomin to lake polygons 
# global lakeshed polygons to lake polygons (geographic area) 
# join coverage data from dbfs to lake polygons 
# join lake polygons for additional attributes $hylak_id 



################################################################################ 
library("foreign")


wd = "E:/research/2020_06_11_Borges_etal_tripical_lakes_CH4/git/Borges_etal_202x_tropic_lakes_CH4"
inDir = paste0(wd, "/in")
outDir = paste0(wd, "/out")


# Toming DOC path:
tomingPath = paste0(inDir, "/Toming_DOC/Toming_et_al_global_DOC_DATA.txt")

# lake polygon shapefile path:
lakePolypath = paste0(inDir, "/hydroLAKES/hydroLAKES_polys/HydroLAKES_polys_v10_tropicsClip.dbf")
lakePolyOutpath = paste0(outDir, "/HydroLAKES_polys/HydroLAKES_polys_v10_tropicsClip.dbf")

# lakeShed path:
lakeShedsDBFpath = paste0(inDir, "/lakeSheds/hydroSHEDS_polygon_tropics.dbf")
lakeshedOutpath = paste0(outDir, "/lakeSheds/hydroSHEDS_polygon_tropics.dbf")

# lakeshed GIEMS paths:
shedCoverageDir = paste0(inDir, "/lakeShed_landCover")
allShedCoveragePaths = list.files(shedCoverageDir, ".dbf", recursive=T, full.names=T)
shedCoverageOutpath = paste0(outDir, "/lakeShed_landCover/lakeShed_GIEMS.csv")

# worldclim climate data:
precipPath = paste0(inDir, "/hydroLAKES_worldclim_mean_DBFs/HydroLAKES_precip.dbf")
tempPath = paste0(inDir, "/hydroLAKES_worldclim_mean_DBFs/HydroLAKES_temp.dbf")
windPath = paste0(inDir, "/hydroLAKES_worldclim_mean_DBFs/HydroLAKES_wind.dbf")

# lake pourPoint shapefile path:
lakePPpath = paste0(inDir, "/hydroLAKES/hydroLAKES_points/HydroLAKES_points_v10_tropicsClip.dbf")










################################################################################  
# Join Toming DOC data to lake polygon:

# read in Toming DOC data:
toming = read.table(tomingPath, header=F, sep="\t")

# update column names:
names(toming) = c("hylak_id", "name", "country", "continent", 
                  "Toming_5", "Toming_6", "Toming_7", "Toming_8")
joinTab = toming[, grep("Toming", names(toming))]


# read in lake polygon data: 
lakePoly = read.dbf(lakePolypath)

# join toming DOC data to lake polygon table:
j = match(lakePoly$Hylak_id, toming$hylak_id)

lakePoly = data.frame(lakePoly, joinTab[j,])

# "Toming_5", "Toming_6", "Toming_7", "Toming_8" --> refer to columns in Toming table



################################################################################ 
# join climatology data to lake polygons:

precip = read.dbf(precipPath)
j = match(lakePoly$Hylak_id, precip$Hylak_id)
prec_mon_mm = precip$MEAN[j]
lakePoly = data.frame(lakePoly, prec_mon_mm)


temp = read.dbf(tempPath)
j = match(lakePoly$Hylak_id, temp$Hylak_id)
temp_C = temp$MEAN[j]
lakePoly = data.frame(lakePoly, temp_C)

wind = read.dbf(windPath)
j = match(lakePoly$Hylak_id, wind$Hylak_id)
wind_mps  = wind$MEAN[j]
lakePoly = data.frame(lakePoly, wind_mps)



################################################################################  
# Join lakeShed areas to lake polygon:

# read in lakeshed data:
lakeshed = read.dbf(lakeShedsDBFpath)

# rename columns: 
names(lakeshed) = c("Id", "hylak_id", "lakeShedArea_km2")

# joing lakeshed area to lake polygon table:
j = match(lakePoly$Hylak_id, lakeshed$hylak_id)

lakeShedArea_km2 = lakeshed$lakeShedArea_km2[j]

lakePoly = data.frame(lakePoly, lakeShedArea_km2)



################################################################################  
# calculate lake shed wetland coverage 
# (~1 min run time --> export written to disk):

# concatenate all contents into a global tropical table:
for (i in 1:length(allShedCoveragePaths)){ 
  
  covDF = read.dbf(allShedCoveragePaths[i])
  
  if (i == 1){
    covDF_glob = covDF
  }else{ 
    covDF_glob = rbind(covDF_glob, covDF)
  }
  
}

# rename column headers:
names(covDF_glob) = sub("VALUE", "GIEMS", names(covDF_glob))
names(covDF_glob)[1] = "hylak_id"

# average land cover of duplicated hylak_id:
dupIDs = covDF_glob$hylak_id[duplicated(covDF_glob$hylak_id)]

for (i in 1:length(dupIDs)){
  j = which(dupIDs[i] == covDF_glob$hylak_id)
  
  if (length(j)==1){ next }
  
  # change first instance of duplicated ID to mean, delete others:
  covDF_glob[j[1],] = colMeans(covDF_glob[j,])
  
  covDF_glob[j[-1],] = NA
  
}

dim(covDF_glob)
covDF_glob = na.omit(covDF_glob)
dim(covDF_glob)


# normalize land cover data so they add to unity:
covDF_vals = covDF_glob[,-1]
rowSums = rowSums(covDF_vals)
covDF_vals_norm = sweep(covDF_vals, 1, "/", STATS=rowSums(covDF_vals))
covDF_glob[, grep("GIEMS", names(covDF_glob))] = covDF_vals_norm


write.csv(covDF_glob, shedCoverageOutpath, row.names=F)



################################################################################
# join lake shed GIEMS land cover data to lake polygons

# lake shed land cover data:
covDF = read.csv(shedCoverageOutpath, header=T)

joinTab = covDF[, grep("GIEMS", names(covDF))]

j = match(lakePoly$Hylak_id, covDF$hylak_id)

lakePoly = data.frame(lakePoly, joinTab[j,])

# write out joined lake poly table:
write.dbf(lakePoly, lakePolyOutpath)



write.csv(lakePoly, sub(".dbf", ".csv", lakePolyOutpath), row.names=F)




################################################################################
# Join attributes to the lake shed shapefile (for display purposes):

# read in lakeshed data:
lakeshed = read.dbf(lakeShedsDBFpath)
names(lakeshed) = c("Id", "hylak_id", "lakeShedArea_km2")


# read in lake poly out table:
lakePoly = read.dbf(lakePolyOutpath)

joinTab = lakePoly[, grep("GIEMS", names(lakePoly), invert=T)]

j = match(lakeshed$hylak_id, lakePoly$Hylak_id)

lakeshed = data.frame(lakeshed, joinTab[j,])



# lake shed land cover data intentionally seperately than from lake poly table 
# b/c GIEM basins might have been dropped in lake poly table join above:
covDF = read.csv(shedCoverageOutpath, header=T)

joinTab = covDF[, grep("GIEMS", names(covDF))]


# join GIEMS data to lake sheds:
j = match(lakeshed$hylak_id, covDF$hylak_id)

lakeshed = data.frame(lakeshed, joinTab[j,])


# write out joined lake shed table:
write.dbf(lakeshed, lakeshedOutpath)







