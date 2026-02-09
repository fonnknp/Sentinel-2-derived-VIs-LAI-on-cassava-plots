library(tidyterra)
library(tidyverse)
library(ggspatial)
library(terra)

#create list files to collect raster file names
listfiles <- list.files("Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2022-01-01_2022-12-31", pattern ="*\\T47PRT.tif$", full = T)

#list files selecting the dates that are closest to the ground-LAI measurement dates 
#[4] "Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2022-01-01_2022-12-31/20220116T034059_20220116T035009_T47PRT.tif"
#[7] "Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2022-01-01_2022-12-31/20220225T033719_20220225T035209_T47PRT.tif"
#[11] "Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2022-01-01_2022-12-31/20220327T033529_20220327T034838_T47PRT.tif"
#[15] "Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2022-01-01_2022-12-31/20220426T033529_20220426T034806_T47PRT.tif"
#[16] "Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2022-01-01_2022-12-31/20220620T033551_20220620T034626_T47PRT.tif"

#import raster files of each MAP, selecting the dates that are closest to the ground-LAI measurement dates
listfiles <- listfiles[c(4,7,11,15,16)]

#stack all rasters 
rast_stack <- rast(listfiles)

#create a function to compute zonal statistics based on median for each image in the stack raster
zonal_median <- function(rast_stack){
  plots47_shp <- vect("Data/Nongbua-shapefile/Nongbua_47plots/Nongbua_47plots.shp") #cassava plots - vector
  n <- nlyr(rast_stack) #the number of layers
  result_list <- list() #create a list to contain results
  for (i in 1:n){
    stack <- rast_stack[[i]] #stack images in each day
    
    zonal.median <- terra::extract(stack,plots47_shp,fun='median',ID=T) #zonal statistic by median
    zonal.count <- terra::extract(stack, plots47_shp, fun = function(x) sum(!is.na(x)))
    
    zonal.median$ID <- plots47_shp$Name #change index to plot names
    zonal.count$ID <- plots47_shp$Name #change index to plot names
    
    colnames(zonal.median) <- c("Plot","Median.index") # change column names
    colnames(zonal.count) <- c("Plot","Pixel.count") # change column names
    
    zonal_data <- merge(zonal.median, zonal.count, by = "Plot")
    
    zonal_data$date <- str_extract(sources(stack), "\\d{8}")
    zonal_data <- zonal_data %>%
      mutate(date = ymd(date),
             MAP = as.double(substr(date, 6, 7)) + 1) 
    
    zonal_data$Index <- names(stack) #create index column
    result_list[[i]] <- zonal_data #collect results into the list
  }
  # Combine all results into a single dataframe
  final_result <- do.call(rbind, result_list)
  return(final_result)
}

#compute zonal analysis by median for all bands----
ALL_INDEX_zonal <- zonal_median(rast_stack) %>% 
  select(MAP,date,Plot,Index,Pixel.count,Median.index)

write.csv(ALL_INDEX_zonal, "Data/all_index_zonal_median.csv", row.names = F)

##############################################################
#M1_list <- list.files("Data/GEE-image/Nongbua_S2-SR-Level2A-Harmonized_allBands_VIs_Res10m_2021-10-01_2023-01-01/", pattern ="\\d{4}(01)\\d{2}T.*\\T47PRT.tif$", full = T)
#M1_stack <- rast(M1_list)
#M1_SeLI <- M1_stack[[grep("^SeLI$", names(M1_stack))]]

#rast_stack[[21]]
#nlyr(rast_stack)
#names(rast_stack)
#zonal.med.test <- extract(rast_stack[[21]],plots47_shp,fun='median',ID=T)
#zonal.med.test$ID <- plots47_shp$Name
#colnames(zonal.med.test) <- c("Plot","Value")
#zonal.med.test$date <- str_extract(sources(rast_stack[[21]]), "\\d{8}")
#zonal.med.test$Index <- names(rast_stack[[21]])
#


