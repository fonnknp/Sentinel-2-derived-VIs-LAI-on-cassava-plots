# ==========================================================================
# Section 2.3.3. Vegetation indices (VIs) computation
# Zonal median extraction of Sentinel-2 vegetation indices per cassava plot and MAP.
# =========================================================================

# Load libraries
library(tidyterra)   # Tidyverse-compatible tools for terra spatial objects
library(tidyverse)   # Data manipulation and visualization
library(ggspatial)   # Add map elements (scale bar, north arrow) to ggplot
library(terra)       # Raster and vector spatial data processing

# List all GeoTIFF raster files from the Google Earth Engine (GEE) export folder
listfiles <- list.files(
  "Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2022-01-01_2022-12-31",
  pattern = "*\\T47PRT.tif$",
  full = T
)

# Select only 5 images with acquisition dates closest to ground-LAI measurement dates (Table 1)
# Index [4]  = 16 Jan 2022 → closest to LAI measurement on 17-18 Jan (2 MAP)
# Index [7]  = 25 Feb 2022 → closest to LAI measurement on 23-24 Feb (3 MAP)
# Index [11] = 27 Mar 2022 → closest to LAI measurement on 19-20 Mar (4 MAP)
# Index [15] = 26 Apr 2022 → closest to LAI measurement on 21-22 Apr (5 MAP)
# Index [16] = 20 Jun 2022 → closest to LAI measurement on 25-26 Jun (7 MAP)
# Note: No image met quality criteria in May 2022 (6 MAP), so this month is excluded
listfiles <- listfiles[c(4, 7, 11, 15, 16)]

# Stack all 5 raster files into a multi-layer SpatRaster object
# Each file contains multiple bands/VIs, so the total number of layers = bands × dates
rast_stack <- rast(listfiles)

# Define a function to compute zonal statistics (median) for each plot polygon
zonal_median <- function(rast_stack) {
  
  # Read the shapefile containing boundaries of all 47 cassava plots
  plots47_shp <- vect("Data/Nongbua-shapefile/Nongbua_47plots/Nongbua_47plots.shp")
  
  # Get the total number of layers in the raster stack
  n <- nlyr(rast_stack)
  
  # Initialize an empty list to store results from each layer
  result_list <- list()
  
  # Loop through each layer (i.e., each band/VI for each date)
  for (i in 1:n) {
    
    # Extract the i-th layer (one band/VI from one date)
    stack <- rast_stack[[i]]
    
    # Compute the median of all pixel values within each plot boundary
    # Median is used instead of mean to reduce the influence of outlier pixels
    # (e.g., from mixed pixels at plot edges or residual cloud artifacts)
    zonal.median <- terra::extract(stack, plots47_shp, fun = 'median', ID = T)
    
    # Count the number of non-NA pixels within each plot boundary
    # This serves as a quality check to ensure sufficient spatial coverage per plot
    zonal.count <- terra::extract(stack, plots47_shp, fun = function(x) sum(!is.na(x)))
    
    # Replace numeric IDs with plot names from the shapefile attribute table
    zonal.median$ID <- plots47_shp$Name
    zonal.count$ID <- plots47_shp$Name
    
    # Rename columns for clarity
    colnames(zonal.median) <- c("Plot", "Median.index")
    colnames(zonal.count) <- c("Plot", "Pixel.count")
    
    # Merge median values and pixel counts by plot name
    zonal_data <- merge(zonal.median, zonal.count, by = "Plot")
    
    # Extract the 8-digit date string (e.g., "20220116") from the raster source filename
    zonal_data$date <- str_extract(sources(stack), "\\d{8}")
    
    # Convert the date string to Date format and calculate Months After Planting (MAP)
    # MAP calculation: extract the month number (characters 6-7 of the date) and add 1
    # This works because cassava was planted around December 2021, so:
    # Jan (01) + 1 = 2 MAP, Feb (02) + 1 = 3 MAP, ..., Jun (06) + 1 = 7 MAP
    zonal_data <- zonal_data %>%
      mutate(
        date = ymd(date),
        MAP = as.double(substr(date, 6, 7)) + 1
      )
    
    # Add the band/VI name from the raster layer name as a new column
    zonal_data$Index <- names(stack)
    
    # Append this layer's results to the list
    result_list[[i]] <- zonal_data
  }
  
  # Combine all individual dataframes into one final dataframe
  final_result <- do.call(rbind, result_list)
  return(final_result)
}

# Run the zonal median function on the stacked rasters
# Then select and reorder columns for a clean output format
ALL_INDEX_zonal <- zonal_median(rast_stack) %>%
  select(MAP, date, Plot, Index, Pixel.count, Median.index)

# Export the results as a CSV file for subsequent statistical analysis
# (e.g., linear mixed-effects models, stage-specific regressions)
write.csv(ALL_INDEX_zonal, "Data/all_index_zonal_median.csv", row.names = F)
###########################################################################
