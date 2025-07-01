resample_dem <- function(dem,dem_ref){

  values(dem)[values(dem)==0] <- NA
  dem <- raster::projectRaster(dem,crs = "EPSG:31983")
  dem <- raster::resample(dem,dem_ref)
  return(dem)

}

extract_river <- function(dem,lon,lat,target_area,radius,n_points, snap_dist){

  # Define the initial outlet and radius (make sure your coordinates are in a projected system for distance accuracy)
  initial_outlet <- c(x = lon, y = lat)  # example coordinates in meters

  #-------------------------------------------------------------------------------
  # Create temporary directory for auxiliary files in the river network analysis
  wd <- tempdir()

  #-------------------------------------------------------------------------------eu p
  # River network analysis

  values(dem)[values(dem)==0] <- NA

  writeRaster(dem,paste0(wd,"/","dem",".tif"), overwrite = T )

  # Prepare DEM for Hydrology Analyses

  wbt_breach_depressions_least_cost(
    dem = "dem.tif",
    output = "dem_breached.tif",
    dist = 5,
    fill = TRUE,
    wd = wd,
    verbose_mode = FALSE)

  wbt_d8_flow_accumulation(input = "dem_breached.tif",
                           output = "D8FA.tif",
                           wd = wd,
                           verbose_mode = FALSE)


  wbt_d8_pointer(dem = "dem_breached.tif",
                 output = "D8pointer.tif",
                 wd = wd,
                 verbose_mode = FALSE)

  threshold <- mean(values(rast(paste0(wd,"//D8FA.tif"))),na.rm=T)

  wbt_extract_streams(flow_accum = "D8FA.tif",
                      output = "raster_streams.tif",
                      threshold = threshold,
                      wd = wd,
                      verbose_mode = FALSE)

  # Create flow accumulation and pointer grids

  # Define the function to generate a random point
  f <- function() {
    # Generate two random numbers and sum them
    u <- runif(1) + runif(1)
    # Generate a random angle between 0 and 2*pi
    t <- runif(1, min = 0, max = 2 * pi)
    # Apply the condition to compute r
    r <- if (u > 1) 2 - u else u
    # Return the (x, y) coordinate
    c(x = r * cos(t), y = r * sin(t))
  }

  # Generate 10,000 random points
  points <- t(replicate(n_points, f()))


  # Compute the x and y coordinates of the random sample points
  random_points <- data.frame(
    x = lon + points[,1]*radius,
    y = lat + points[,2]*radius
  )


  st_as_sf(random_points,
           coords = c("x","y"), crs = "EPSG:31983") -> pour_point

  plot(dem)
  plot(pour_point,add=T)

  st_write(pour_point,paste0(wd,"\\pourpoints.shp"), delete_layer  = TRUE)

  wbt_jenson_snap_pour_points(pour_pts = "pourpoints.shp",
                              streams = "raster_streams.tif",
                              output =  "snappedpp.shp",
                              snap_dist = 50,
                              wd = wd,
                              verbose_mode = FALSE) #careful with this! Know the units of your data


  outlets_snapped <- st_read(file.path(wd, "snappedpp.shp"), quiet = TRUE)

  wbt_watershed(d8_pntr = "D8pointer.tif",
                pour_pts = "snappedpp.shp",
                output = "catchment.tif",
                wd = wd,
                verbose_mode = FALSE)

  # Vectorize catchments

  catchments <- st_sfc(crs = "EPSG:31983")
  for(id in outlets_snapped$FID){
    # Write outlet shapefile for the current id
    st_write(outlets_snapped[outlets_snapped$FID == id, ],
             file.path(wd, paste0(id, "_outlet.shp")),
             delete_layer = TRUE, quiet = TRUE)

    # Run watershed delineation using WhiteboxTools
    whitebox::wbt_watershed(d8_pntr = "D8pointer.tif",
                            pour_pts = paste0(id, "_outlet.shp"),
                            output = paste0(id, "_catchment.tif"),
                            wd = wd)

    # Read the delineated catchment raster
    drainage <- read_stars(file.path(wd, paste0(id, "_catchment.tif")))

    # Convert the raster catchment to vector polygons
    contours <- tryCatch({
      st_contour(drainage, breaks = 1) |>
        st_geometry() |>
        st_cast("POLYGON")
    }, error = function(e){
      message(paste("Error in st_contour for FID", id, ":", e$message))
      NA
    })

    # Choose the largest polygon by area
    contours <- tryCatch({
      contours[which.max(st_area(contours))]
    }, error = function(e){
      message(paste("Error in selecting max area for FID", id, ":", e$message))
      NA
    })

    # If contours extraction was successful, add to catchments
    if(!is.na(contours)[1]){
      new_sf <- st_sf(id = id,
                      geometry = st_sfc(contours, crs = st_crs(catchments)))
      catchments <- rbind(catchments, new_sf)
    }
  }

  compare_areas <- data.frame(
    name = catchments$id,
    expected_area =  target_area,
    computed_area = st_area(catchments) |> units::set_units("km^2") |> round(1)
  )
  compare_areas %>%
    mutate(
      error = abs(expected_area-as.numeric(computed_area))
    ) %>%
    subset(
      error < 10
    ) %>% arrange(error) -> compare_areas

  catchments <- catchments %>% subset(id == compare_areas$name[1])

  print(tm_shape(dem) +
          tm_raster(palette = terrain.colors(10))+
          tm_shape(catchments)+
          tm_borders("black", lwd = .5))

  streams <- rast(paste0(wd,paste0("\\raster_streams.tif")))

  wbt_raster_streams_to_vector(
    streams = "raster_streams.tif",
    d8_pntr = "D8pointer.tif",
    output = "streams.shp",
    wd = wd,
    verbose_mode = FALSE
  )

  catchments <- st_transform(catchments,crs = "EPSG:31983")

  stream_newtork <- st_read(paste0(wd,"\\streams.shp"),crs="EPSG:31983")
  stream_newtork <- st_intersection(stream_newtork,catchments)
  plot(stream_newtork)

  stream_newtork <- stream_newtork %>% filter(st_is(. , "LINESTRING"))
  sf::st_write(stream_newtork,paste0(wd,"\\stream.shp"), delete_layer = T)

  wbt_vector_stream_network_analysis(
    streams =  "stream.shp",
    output = "stream_analysis.shp",
    snap = 0.1,
    wd = wd,
    verbose_mode = F
  )

  # Stream network analysis

  stream_newtork <- st_read(paste0(wd,"\\stream.shp"))
  streams_analysis <- st_read(paste0(wd,"\\stream_analysis.shp"))


  streams_new <- cbind(stream_newtork,
                       st_drop_geometry(streams_analysis))

  dem2 <- crop(dem, extent(catchments))
  dem3 <- mask(dem2, catchments)

  return(list(centerline_sf = streams_new,dem = dem3))

}


carve_dem <- function(dem,centerline_sf){

  elv <- terra::extract(rast(dem),centerline_sf,method = "bilinear",cells=TRUE)
  colnames(elv) <- c("ID","Elv","cell")

  strahler_max <- max(centerline_sf$STRAHLER)

  centerline_sf <- centerline_sf %>% mutate(carve = ifelse(STRAHLER == strahler_max,3,
                                                           ifelse(STRAHLER == strahler_max-1,2,0.5)))

  centerline_sf <- merge(centerline_sf, elv, by.x="FID",by.y="ID")
  centerline_sf <- centerline_sf %>%
    mutate(Elv = Elv-carve)

  dem[centerline_sf$cell] <- centerline_sf$Elv

  dem <- projectRaster(dem,crs = "EPSG:31983")

  return(dem)

}

dem_flow_correction <- function(dem_carved){
  wd <- tempdir()
  writeRaster(dem_carved,file.path(wd,"dem.tif"),overwrite=T)

  # 1. Generate Flow Direction and Flow Accumulation with whitebox
  wbt_fill_depressions(dem = 'dem.tif', output = "filled_dem.tif",wd = wd)
  wbt_d8_pointer(dem = "filled_dem.tif", output = "flow_dir.tif",wd = wd)
  wbt_d8_flow_accumulation(input = "filled_dem.tif", output = "flow_acc.tif", out_type = "cells",wd = wd)
  wbt_extract_streams(flow_accum = 'flow_acc.tif',threshold = 1000, output = 'streams.tif',wd=wd)

  wbt_stream_link_identifier(d8_pntr = "flow_dir.tif",streams = 'streams.tif',output="stream_id.tif",wd=wd)
  wbt_strahler_stream_order(d8_pntr = "flow_dir.tif",streams = 'streams.tif',output = "strahler.tif",wd =wd)

  # 2. Define stream network (threshold: e.g., >1000 upstream cells)
  acc <- rast(file.path(wd,"flow_acc.tif"))
  streams <- acc > 3000
  plot(streams)
  writeRaster(streams, file.path(wd,"stream_raster.tif"), overwrite = TRUE)

  # 3. Load everything into terra
  flowdir <- rast(file.path(wd,"flow_dir.tif"))
  flowacc <- rast(file.path(wd,"flow_acc.tif"))
  streams <- rast(file.path(wd,"stream_raster.tif"))
  streams_id <- rast(file.path(wd,"stream_id.tif"))
  strahler <- rast(file.path(wd,"strahler.tif"))

  # 4. Get cell indices of stream pixels
  stream_cells <- which(values(streams) == 1)


  # 5. Extract values to iterate

  flowdir_vals <- values(flowdir)[stream_cells]
  flowacc_vals <- values(flowacc)[stream_cells]
  dem_vals <- values(dem_carved)[stream_cells]
  stream_id_vals <- values(streams_id)[stream_cells]
  strahler_vals <- values(strahler)[stream_cells]

  # Create a data frame of stream pixels sorted by flow accumulation (upstream first)
  stream_df <- data.frame(
    stream_cells = stream_cells,
    id = stream_id_vals,
    strahler = strahler_vals,
    flowdir = flowdir_vals,
    flowacc = flowacc_vals,
    elevation = dem_vals
  ) %>%
    arrange(id)


  ids <- sort(unique(stream_id_vals))

  new_elv <- list()


  # Check and identify streams that for a sequence


  min_slope <- 1e-2  # you can tweak this to be stronger or weaker

  remove_bumps <- function(x) {
    slope <- (max(x)-min(x))/length(x)
    for (i in 2:length(x)) {
      if (x[i] > x[i - 1]) {
        x[i] <- x[i - 1]-min_slope
      }
    }
    return(x)
  }

  for(i in 1:length(ids)){
    str <- stream_df[stream_df$id==ids[i],]
    if(nrow(str) < 5){
      next
    }
    else{
      if(min(str$strahler)==1){
        str %>%
          arrange(flowacc) %>%
          mutate(new_elv = remove_bumps(elevation)) -> str

        new_elv[[i]] <- str


      }

      # For Strahler higher than 1, the junctions must be observed

      else if(min(str$strahler)==2){
        str %>%
          arrange(flowacc) %>%
          mutate(new_elv = remove_bumps(elevation)) -> str

        new_elv[[i]] <- str

      }
      else if(min(str$strahler)==3){
        str %>%
          arrange(flowacc) %>%
          mutate(new_elv = remove_bumps(elevation)) -> str

        new_elv[[i]] <- str
      }
      else if(min(str$strahler)==4){
        str %>%
          arrange(flowacc) %>%
          mutate(new_elv = remove_bumps(elevation)) -> str

        new_elv[[i]] <- str

      }
      else{
        str %>%
          arrange(flowacc) %>%
          mutate(new_elv = remove_bumps(elevation)) -> str

        new_elv[[i]] <- str

      }
    }
  }

  plot(new_elv[[48]]$elevation,type='l',col='blue',lwd=2)
  lines(new_elv[[48]]$new_elv,col='red',lwd=2)


  do.call("rbind",new_elv) -> new_elv

  dem_corrected <- dem_carved

  values(dem_corrected)[new_elv$stream_cells] <- new_elv$new_elv

  return(dem_corrected)

}

#' Run ForwardMole 1D Python Model from R
#'
#' Configures Python environment using `reticulate`, sources the ForwardMole 1D Python script,
#' and runs the `FM1D` function to perform 1D hydrological processing on the DEM and hydrograph data.
#'
#' @param workspace Character. Base working directory path where output folders reside. Default is `"."`.
#' @param name Character. Base name prefix used to build DEM and hydrograph file names.
#' @param python_path Character. Path to the Python executable to use with `reticulate`. Default is `"C:/Users/marco/anaconda3/python.exe"`.
#' @param python_script_path Character. Path to the Python script to source. Default is `"Data/forward_mole_1D_2.py"`.
#'
#' @return Invisible `NULL`. Side effects include executing the Python FM1D function, which processes and saves outputs in the specified workspace.
#'
#' @details
#' This function is a wrapper to run a Python-based hydrological model from within R. It sets the Python environment, sources the Python script,
#' and calls the `FM1D` function with constructed paths based on the input parameters.
#'
#' Make sure the Python script `forward_mole_1D_2.py` is accessible at the specified path and the Python environment has all necessary dependencies installed.
#'
#' @references
#' Rápalo, L. M. C., Gomes Jr, M. N., & Mendiondo, E. M. (2024).
#' Developing an open-source flood forecasting system adapted to data-scarce regions:
#' A digital twin coupled with hydrologic-hydrodynamic simulations. *Journal of Hydrology*.
#' https://doi.org/10.1016/j.jhydrol.2024.131929
#'
#' Source code: \url{https://github.com/luiscastillo1993/Forward-Mole}
#'
#' @examples
#' \dontrun{
#' ForwardMole1D(
#'   workspace = "C:/Projects/HydroModel",
#'   name = "basin01",
#'   python_path = "C:/Users/marco/anaconda3/python.exe",
#'   python_script_path = "Data/forward_mole_1D_2.py"
#' )
#' }
#'
#' @export
ForwardMole1D <- function(workspace = ".",
                          name,
                          python_path = "C:/Users/marco/anaconda3/python.exe",
                          python_script_path = "Data/forward_mole_1D_2.py") {
  # Validate inputs
  if (missing(name) || !nzchar(name)) {
    stop("Argument 'name' must be provided and non-empty.")
  }
  if (!file.exists(python_script_path)) {
    stop("Python script not found at: ", python_script_path)
  }
  if (!file.exists(python_path)) {
    warning("Python executable not found at ", python_path, ". Reticulate may use default Python.")
  }

  # Configure Python environment
  reticulate::use_python(python_path, required = TRUE)
  py_config <- reticulate::py_config()
  message("Using Python at: ", py_config$python)

  # Source the Python script
  reticulate::source_python(python_script_path)

  # Construct paths for FM1D function
  output_path <- file.path(workspace, "output2")
  dem_file <- paste0(name, "_carved")
  hydro_file <- paste0(name, "_centerline")

  # Call the Python FM1D function
  if (!exists("FM1D")) {
    stop("FM1D function not found after sourcing ", python_script_path)
  }

  FM1D(path2 = output_path,
       dem_name = dem_file,
       hydro_name = hydro_file)

  invisible(NULL)
}
