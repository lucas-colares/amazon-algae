# ##########################################################################
# ############## Landscape analysis for Palheta et al. (2023) ##############
# ################### By: Lucas Colares ####################################
# ################### Date: 20-08-2023 #####################################
# ##########################################################################

# Load required packages for spatial and raster data analysis
library(raster)  # Handles raster datasets
library(sf)      # Handles spatial vector data (shapefiles)

# Load raster data (land cover for Brazil in 2018 and 2019)
raster("datasets/rasters/brasil_coverage_2019.tif") -> brazil2019
raster("datasets/rasters/brasil_coverage_2018.tif") -> brazil2018

# Load shapefile containing municipality boundaries in Pará (PA)
st_read("datasets/shapefiles/PA_Municipios_2022.shp") -> PAshp

# Load coordinates data from CSV file
read.csv("datasets/coordinates.csv", sep=";", row.names = 1) -> coords

# Match coordinate rows with the HII (Human Influence Index) data
coords[match(rownames(HII), rownames(coords)), ] -> coords
coords$IIH <- HII$IIH_total  # Add HII total values to the coordinates

# Save updated coordinates file
write.csv(coords, "datasets/coordinates.csv", fileEncoding = "latin1", row.names = F)

# Crop raster data to the extent of the Pará (PA) shapefile
raster::crop(brazil2019, raster::extent(PAshp)) -> PAcov2019
raster::crop(brazil2018, raster::extent(PAshp)) -> PAcov2018

# Remove original full Brazil rasters to free memory
rm(brazil2019, brazil2018)

# Create a spatial points dataframe using the coordinate data
spdf <- SpatialPointsDataFrame(coords = coords, data = coords, 
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Create buffers of various sizes around spatial points
Buffers <- {}  # Initialize an empty list
for(x in c(seq(30,90,10), seq(100,1000,100))) {
  st_buffer(st_as_sf(spdf), dist = x) -> Buffers[[length(Buffers) + 1]]
}
names(Buffers) = c(seq(30,90,10), seq(100,1000,100))  # Assign buffer names

# Crop the raster data to the extent of the largest buffer
raster::crop(PAcov2019, extent(Buffers[[length(Buffers)]])) -> FinalCov2019
raster::crop(PAcov2018, extent(Buffers[[length(Buffers)]])) -> FinalCov2018

# Remove intermediate cropped rasters to free memory
rm(PAcov2019, PAcov2018)
gc()  # Run garbage collection

# Extract land cover values within buffers
AllValues2019 <- {}  # Initialize empty list
for(x in 1:length(Buffers)) {
  for(y in 1:nrow(Buffers[[x]])) {
    mask2019 = raster::mask(FinalCov2019, Buffers[[x]][y,])  # Mask raster using buffer
    Values2019 <- raster::getValues(mask2019)  # Extract values within buffer
    
    # Reclassify certain land cover categories
    gsub(paste0("^", c(19,39,20,40,62,41,36,46,47,48), "$", collapse = "|"), 18, Values2019) -> Values2019
    
    # Convert values to frequency table
    data.frame(table(na.omit(Values2019))) -> Freq
    colnames(Freq)[1] = "Cat"  # Rename category column
    
    # Store results with buffer and site information
    data.frame(Freq, Site = rownames(Buffers[[x]][y,]), Buffer = names(Buffers)[x]) -> AllValues2019[[length(AllValues2019) + 1]]
  }
}

# Combine all extracted values into a single dataframe
do.call("rbind", AllValues2019) -> AllValues2019

# Remove category 33 from the dataset
AllValues2019[!AllValues2019$Cat == 33, ] -> AllValues2019

# Create a base dataframe to store proportions of land cover categories
data.frame(matrix(NA, nrow = length(unique(AllValues2018$Site)), 
                  ncol = length(unique(c(unique(AllValues2019$Cat), unique(AllValues2018$Cat)))))) -> DataFrameBase
colnames(DataFrameBase) = unique(c(unique(AllValues2019$Cat), unique(AllValues2018$Cat)))
rownames(DataFrameBase) = unique(AllValues2019$Site)

# Normalize land cover data by buffer
Final2019 <- {}  # Initialize empty list
for(x in names(Buffers)) {
  AllValues2019[AllValues2019$Buffer == x, ] -> SelBuf  # Select data for buffer size x
  NewData = DataFrameBase  # Copy base dataframe
  
  for(u in rownames(NewData)) {
    SelBuf[SelBuf$Site == u, ] -> SelSite  # Select data for specific site
    SelSite$Freq / sum(SelSite$Freq) -> SelSite$Freq  # Normalize frequency
    
    # Fill the dataframe with normalized frequencies
    NewData[rownames(NewData) == u, ][, match(SelSite$Cat, colnames(NewData))] <- SelSite$Freq
  }
  Final2019[[length(Final2019) + 1]] <- NewData  # Store processed data
}
names(Final2019) = names(Buffers)  # Assign buffer names

# Replace NA values with 0
for(x in 1:length(Final2019)) {
  Final2019[[x]][is.na(Final2019[[x]])] <- 0
}

### Select scale of effect using stepwise regression
Adjs <- {}  # Initialize empty list
for(x in 1:length(Final2019)) {
  lm(HII$IIH_total ~ ., data = Final2019[[x]]) -> LMIHH  # Linear model
  step(LMIHH, direction = "both") -> LMIHH  # Stepwise model selection
  summary(LMIHH)  # Display model summary
  
  # Store adjusted R-squared values for model selection
  data.frame(Buffer = names(Final2019)[x], RAdj = RsquareAdj(LMIHH)[[2]]) -> Adjs[[length(Adjs) + 1]]
}
do.call("rbind", Adjs) -> Adjs
Adjs[order(Adjs$RAdj, decreasing = T), ]  # Sort results by highest adjusted R-squared

# Calculate species richness (S) and density (Dens)
rowSums(decostand(comm[, 1:ncol(comm) - 1], "pa")) -> S  # Species richness
rowSums(comm[, 1:ncol(comm) - 1]) -> Dens  # Species density

# Select scale of effect for species richness using Poisson regression
Adjs_Rich <- {}
for(x in 1:length(Final2019)) {
  glm(S ~ ., data = Final2019[[x]], family = "poisson") -> LMIHH  # Poisson regression
  step(LMIHH, direction = "both") -> LMIHH  # Stepwise selection
  summary(LMIHH)  # Display model summary
  
  # Store adjusted R-squared values
  data.frame(Buffer = names(Final2019)[x], 
             RAdj = with(summary(LMIHH), 1 - deviance/null.deviance)) -> Adjs_Rich[[length(Adjs_Rich) + 1]]
}
do.call("rbind", Adjs_Rich) -> Adjs_Rich
Adjs_Rich[order(Adjs_Rich$RAdj, decreasing = T), ]  # Sort results

# Select scale of effect for species density using log-transformed linear regression
Adjs_Dens <- {}
for(x in 1:length(Final2019)) {
  lm(log10(Dens) ~ ., data = Final2019[[x]]) -> LMIHH  # Log-transformed regression
  step(LMIHH, direction = "forward") -> LMIHH  # Forward stepwise selection
  summary(LMIHH)  # Display model summary
  
  # Store adjusted R-squared values
  data.frame(Buffer = names(Final2019)[x], RAdj = RsquareAdj(LMIHH)[[2]]) -> Adjs_Dens[[length(Adjs_Dens) + 1]]
}
do.call("rbind", Adjs_Dens) -> Adjs_Dens
Adjs_Dens[order(Adjs_Dens$RAdj, decreasing = T), ]  # Sort results
