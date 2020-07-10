#!/usr/bin/env Rscript
################################################################################
# Borges_etal_202x_tropicalLakesCH4emission_finalTableJoin.R
################################################################################
# George Allen, June-July 2020



################################################################################ 
library("foreign")


wd = "E:/research/2020_06_11_Borges_etal_tripical_lakes_CH4/git/Borges_etal_202x_tropic_lakes_CH4"
inDir = paste0(wd, "/in")
outDir = paste0(wd, "/out")



# Toming DOC path:
tomingPath = paste0(inDir, "/Toming_DOC/Toming_et_al_global_DOC_DATA.txt")

# lake polygon shapefile path:
lakePolypath = paste0(inDir, "/hydroLAKES/hydroLAKES_polys/HydroLAKES_polys_v10_tropicsClip.dbf")


# lakeShed path:
lakeShedsDBFpath = paste0(inDir, "/lakeSheds/hydroSHEDS_polygon_tropics.dbf")


# lakeshed GIEMS paths:
shedCoverageDir = paste0(inDir, "/lakeShed_landCover")

contList = c("af", "as", "au", "ca", "eu", "na", "sa")
covList = c("GIEMS", "TCI", "WTD", "RFW")
outShedCoverageCSVpaths = paste0(outDir, "/CSV/", contList, ".csv")
allShedCoveragePaths = list.files(shedCoverageDir, ".dbf", recursive=T, full.names=T)


# lake pourPoint shapefile path:
lakePPpath = paste0(inDir, "/hydroLAKES/hydroLAKES_points/HydroLAKES_points_v10_tropicsClip.dbf")


# want to make a map of lakesheds and their land cover type --> need to join lake attributes to lake sheds
# want to make a map of lake polygons and their climatology --> need to join shed attributes to lakes

# join data from Tomin to lake polygons --> $hylak_id 
# global lakeshed polygons to lake polygons (get geographic area) --> $gridcode 
# join coverage data from dbfs to lake polygons  --> $HydroLAKES_id


# join lake polygons for additional attributes $hylak_id 


# remove reservoirs from lakes table and send to Alberto 

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
# calculate lake shed wetland coverage: 

for (i in 1:length(contList)){ #i = 1
  
  contPaths = grep(contList[i], allShedCoveragePaths, value=T)
  
  for (j in 1:length(covList)){ #j = 1
  
    covPath = grep(covList[j], contPaths, value=T)
    
    covDF = foreign::read.dbf(covPath)
    
    # first column of data frame is the index of watershed. 
    # remove first column and normalize other columns:
    covDF_vals = covDF[,-1]
    rowSums = rowSums(covDF_vals)
    covDF_vals_norm = sweep(covDF_vals, 1, "/", STATS=rowSums(covDF_vals))
    
    names(covDF_vals_norm) = paste0(covList[j], "_", 0:(ncol(covDF_vals_norm)-1))
    
    if (j==1){ 
      covDF_out = data.frame(covDF[,1], covDF_vals_norm)
      names(covDF_out)[1] = "hylak_id"
    }else{
      
      # HydroLAKES IsD match between coverage layers:
      covDF_out = data.frame(covDF_out, covDF_vals_norm)
    }
  }
  
  write.csv(covDF_out, outShedCoverageCSVpaths[i], row.names=F)
}

# rbind all land cover tables together:
for (i in 1:length(outShedCoverageCSVpaths)){
  covDF = read.csv(outShedCoverageCSVpaths[i], header=T)
  
  if (i == 1){
    covDF_glob = covDF
  }else{ 
    covDF_glob = rbind(covDF_glob, covDF)
  }
  
}

# average land cover of duplicated hybas:
dupI = which(duplicated(covDF_glob$hylak_id))

for (i in 1:length(dupI)){
  j = which(covDF_glob$hylak_id[dupI[i]] == covDF_glob$hylak_id)
  colMeans(covDF_glob[j,])
}



################################################################################
# join lake shed GIEMS land cover data to lake polygons




length()
  
  # only need to use GIEMS data:
  giems_cols = c(1, grep("GIEMS", names(covDF)))
  
  joinTab = covDF[, giems_cols]
  
  # join toming DOC data to lake polygon table:
  j = match(lakePoly$Hylak_id, covDF$hylak_id)
  
  if (i == 1){
    lakePoly = data.frame(lakePoly, joinTab[j,])
  }else{
    names(lakePoly)
  }
  
}



################################################################################ 
# join coverage data to lakeShed shapefiles:
lSheds = read.dbf(lakeShedsDBFpath)

names(lSheds)[names(lSheds) == 'gridcode'] = "hylak_id"
bindTab = as.data.frame(array(NA, c(nrow(lSheds, 5))))

for (i in 1:length(outShedCoverageCSVpaths)){
  covDF = read.csv(outShedCoverageCSVpaths[i], header=T)
  
  # only need to use GIEMS data:
  giems_cols = c(1, grep("GIEMS", names(covDF)))
  
  #joinInd = match(lSheds$hylak_id, covDF$HydroLAKES_ID)
  
  joinInd = match(covDF$HydroLAKES_ID, lSheds$hylak_id)
  


  bindTab[joinInd] =  covDF[,giems_cols]
  
  lSheds$gridcode[62077]
  covDF[joinInd] 
  
}



lSheds$gridcode

################################################################################ 
# join lake area data: 







