## ---------------------------------------------------------------------------------------------------------------------------
rm(list = ls())

# Load the necessary libraries
library(raster)
library(sf)
library(ggplot2)
library(exactextractr)
library(dplyr)


## ---------------------------------------------------------------------------------------------------------------------------
spain_shpfile <- st_read("../data/raw/georef-spain-comunidad-autonoma/georef-spain-comunidad-autonoma-millesime.shp")

# Subset the regions of interest: Catalunya and Aragon
CT <- spain_shpfile %>% 
  filter(acom_iso316 == "CT")

AR <- spain_shpfile %>% 
  filter(acom_iso316 == "AR")

# Calculate shared border between both regions (computing the intersection between both boundaries)
shared_border <- st_intersection(st_boundary(CT), st_boundary(AR))

# Plot the shared border
ggplot() + 
  geom_sf(data = CT, fill = "white", color = "black", alpha = 0.5) +
  geom_sf(data = AR, fill = "white", color = "black", alpha = 0.5) +
  geom_sf(data = shared_border, color = "red", size = 2)

## ---------------------------------------------------------------------------------------------------------------------------
#test <- raster("../data/VIIRS-Land_v001_JP113C1_NOAA-20_20231002_c20240124060036.nc")
#plot(test)


## ---------------------------------------------------------------------------------------------------------------------------
process_geojson_files <- function(region) {
  # Directory containing the exported GeoJSON files
  output_dir <- paste0("../data/test_polygons/", region, "/")

  # List all GeoJSON files in the directory
  geojson_files <- list.files(output_dir, pattern = "\\.geojson$", full.names = TRUE)

  # Load all GeoJSON files into a list
  multipolygons_list <- lapply(geojson_files, st_read)

  # Combine all multipolygons into a single sf object
  multipolygons <- do.call(rbind, multipolygons_list)

  # Check for and repair invalid geometries
  multipolygons <- st_make_valid(multipolygons)

  # Remove empty geometries
  multipolygons <- multipolygons[!st_is_empty(multipolygons), ]

  return(multipolygons)
}

plot_multipolygons <- function(multipolygons) {
  ggplot() +
    geom_sf(data = CT) +
    geom_sf(data = AR) +
    geom_sf(data = multipolygons, fill = "blue", color = "black", alpha = 0.5) +
    ggtitle("Multipolygons from Masks") +
    theme_minimal()
}



## ---------------------------------------------------------------------------------------------------------------------------
# Load and process the GeoJSON files for the specified region
multipolygons_ar <- process_geojson_files("ar")
multipolygons_cat <- process_geojson_files("cat")


## ---------------------------------------------------------------------------------------------------------------------------
# Plot the multipolygons
#plot_multipolygons(multipolygons_ar)
#plot_multipolygons(multipolygons_cat)


## ---------------------------------------------------------------------------------------------------------------------------
process_nc_files <- function(nc_files, region, multipolygons, output_directory, dates) {
  # Add a polygon column to the multipolygons dataframe
  multipolygons$polygon <- 1:nrow(multipolygons)
  
  # Initialize a list to store NDVI values for each date
  ndvi_values_list <- list()
  
  for (i in seq_along(nc_files)) {
    file_path <- nc_files[i]
    
    # Load the raster data from the NetCDF file
    raster_data <- raster(file_path)
    
    # Transform the Spain polygon to the CRS of the raster
    spain <- st_transform(spain_shpfile, crs(raster_data))
    
    # Crop the raster to the extent of Spain
    ndvi_spain <- crop(raster_data, extent(spain_shpfile))
    
    # Mask the raster with the Spain polygon
    ndvi_spain <- mask(ndvi_spain, spain_shpfile)
    
    # Check and transform the CRS of the polygons to match the raster CRS if necessary
    if (st_crs(multipolygons) != crs(raster_data)) {
      multipolygons <- st_transform(multipolygons, crs(raster_data))
    }
    
    # Extract the average NDVI values for each multipolygon using exactextractr
    average_ndvi_values <- exact_extract(ndvi_spain, multipolygons, fun = "mean")
    
    # Store the NDVI values in the list
    ndvi_values_list[[i]] <- average_ndvi_values
  }
  
  # Combine NDVI values with the multipolygons dataframe
  for (i in seq_along(dates)) {
    ndvi_col_name <- paste0("ndvi", format(dates[i], "%Y%m"))
    multipolygons[[ndvi_col_name]] <- ndvi_values_list[[i]]
  }
  
  # Remove multipolygons with NA values for all NDVI dates
  multipolygons <- multipolygons %>% filter(rowSums(is.na(select(., starts_with("ndvi")))) != length(dates))
  
  # Calculate centroids of the multipolygons
  centroids <- st_centroid(multipolygons)
  
  # Ensure centroids retain the polygon and NDVI columns
  centroids <- centroids %>% 
    select(polygon, starts_with("ndvi"))
  
  # Save the centroids to a shapefile
  output_file <- paste0(output_directory, "/", region, "/ndvi.shp")
  st_write(centroids, output_file, delete_layer = TRUE)
  
  cat("Processed and saved:", output_file, "\n")
}


## ---------------------------------------------------------------------------------------------------------------------------
# Directory paths and file listing
nc_directory <- "../data/ndvi/"
output_directory <- "../data/processed_ndvi/"

# List all NetCDF files in the directory, sorted by name
nc_files <- list.files(nc_directory, pattern = "\\.nc$", full.names = TRUE)
nc_files <- sort(nc_files)

# Define the dates corresponding to the files
dates <- seq(as.Date("2023-09-01"), as.Date("2024-06-01"), by = "month")

# Print the file names to verify order matches dates
cat("NetCDF files:\n")
print(nc_files)
cat("Corresponding dates:\n")
print(dates)

# Process the NetCDF files
process_nc_files(nc_files, region = "ar", multipolygons_ar, output_directory, dates)
process_nc_files(nc_files, region = "cat", multipolygons_cat, output_directory, dates)


## ---------------------------------------------------------------------------------------------------------------------------
ndvi_cat <- st_read("../data/processed_ndvi/cat/ndvi.shp")
ndvi_ar <- st_read("../data/processed_ndvi/ar/ndvi.shp")



