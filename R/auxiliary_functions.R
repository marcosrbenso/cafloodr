#' Create a Controller File for CAFlood Simulations
#'
#' Generates a text-based controller file required by the CAFlood hydrodynamic model. This file defines simulation parameters, model input/output files, and numerical solver settings.
#'
#' @param Input Character. Path to the folder where the controller file will be saved, typically the model input directory.
#' @param run_name Character. Optional. Full descriptive name of the simulation. Default is `"This is a CAFLOOD simulation"`.
#' @param run_short_name Character. Required short name used to label outputs and log files.
#' @param time_start Numeric. Simulation start time in seconds. Default is `0`.
#' @param time_end Numeric. Simulation end time in seconds. Default is `86400` (1 day).
#' @param max_DT Numeric. Maximum dynamic time step in seconds. Controls adaptive time stepping.
#' @param min_DT Numeric. Minimum dynamic time step in seconds. Prevents instability by limiting smallest time step.
#' @param alpha_fraction Numeric. Weighting factor for implicit/explicit solution of the continuity equation. Typically between 0 and 1.
#' @param max_iterations Integer. Maximum number of iterations for the solver. Default is `1000000`.
#' @param roughness_global Numeric. Manning's roughness coefficient (dimensionless), applied globally across the domain.
#' @param ignore_WD Numeric. Water depth threshold (in meters) below which cells are considered dry.
#' @param tolerance Numeric. Convergence tolerance for the numerical solver.
#' @param slope_tolerance Numeric. Minimum slope value used to avoid flat surface numerical instability.
#' @param boundary_ele Numeric. Elevation used for boundary conditions.
#' @param dem_path Character. Filename of the Digital Elevation Model (DEM) file used for topography input.
#' @param dem_name Character. Internal name or reference for the DEM (e.g., `"dem_main"`).
#' @param rain_path Character. Filename of the rainfall input file.
#' @param flow_path Character. Filename of the inflow (discharge or water level) input file, if applicable.
#' @param output_period Numeric. Time step in seconds for writing output files. Default is `300`.
#' @param Raster_WD_Tolerance Numeric. Water depth threshold (in meters) used for raster outputs to identify wet cells.
#' @param Upstream_Reduction Numeric. Coefficient (between 0 and 1) to reduce flow influence from upstream cells in CAFlood’s routing scheme.
#'
#' @return Writes a controller text file to the specified `Input` folder. The file contains all necessary configuration settings for the CAFlood model to run a simulation.
#'
#' @details
#' This function is part of a wrapper workflow for running CAFlood simulations in R. It builds a properly formatted controller file (`.txt`) containing references to DEMs, rainfall inputs, flow boundaries, and numerical solver options. The file is read by the CAFlood executable during simulation.
#'
#' Ensure all input files (`dem_path`, `rain_path`, `flow_path`) exist in the specified directory and are formatted according to CAFlood requirements.
#'
#' @examples
#' \dontrun{
#' caflood_controller(
#'   Input = "C:/Projects/Input",
#'   run_short_name = "flood_run01",
#'   roughness_global = 0.035,
#'   boundary_ele = 580,
#'   dem_path = "dem_sp.tif",
#'   dem_name = "DEM_SP",
#'   rain_path = "rain_event_2021.txt",
#'   flow_path = "inflow_event_2021.txt"
#' )
#' }
#'
#' @export
caflood_controller <- function(
    Input,
    run_name = "This is a CAFLOOD simulation",
    run_short_name,
    time_start = 0,
    time_end = 86400,
    max_DT = 60,
    min_DT = 0.01,
    alpha_fraction = 0.1,
    max_iterations = 1000000,
    roughness_global,
    ignore_WD = 0.0001,
    tolerance = 0.0001,
    slope_tolerance = 0.528,
    boundary_ele,
    dem_path,
    dem_name,
    rain_path,
    flow_path,
    output_period = 300,
    Raster_WD_Tolerance = 0.01,
    Upstream_Reduction = 0.5
){

  writeLines(
    c(paste0("Simulation Name			,",run_name ),
      paste0("Short Name (for outputs)			,",run_short_name),
      "Version	   			, 1,0,0",
      "Model Type			, WCA2Dv2",
      paste0("Time Start (seconds)		,", time_start),
      paste0("Time End   (seconds)		,",  time_end),
      paste0("Max DT (seconds)		,",max_DT),
      paste0("Min DT (seconds)		,", min_DT),
      paste0("Update DT (seconds)		,", 60),
      paste0("Alpha (Fraction DT 0.0-1.0)	,",alpha_fraction),
      paste0("Max Iterations			, ",max_iterations),
      paste0("Roughness Global		, ",roughness_global),
      paste0("Ignore WD (meter)		, ", ignore_WD),
      paste0("Tolerance (meter)		, ", tolerance),
      paste0("Slope Tolerance (%)		, ", slope_tolerance),
      paste0("Boundary Ele (Hi/Closed-Lo/Open), ", boundary_ele),
      paste0("Elevation ASCII			, ", dem_path),
      paste0("Rain Event CSV			, ", rain_path),
      paste0("Water Level Event CSV		, "),
      paste0("Inflow Event CSV		, ", flow_path),
      paste0("Time Plot CSV			, ", paste0(dem_name,"_WLpoints.csv,"), paste0(dem_name,"_VELpoints.csv")),
      paste0("Raster Grid CSV			, ", "WDraster.csv,", "VELraster.csv"),
      paste0("Output Console			, ", 'true'),
      paste0("Output Period (s)		, ", output_period),
      paste0("Output Computation Time		, ", 'true'),
      paste0("Check Volumes	   		, ", 'true'),
      paste0("Remove Proc Data (No Pre-Proc)	, ", 'true'),
      paste0("Remove Pre-Proc Data		, ", 'true'),
      paste0("Raster VEL Vector Field		, ", 'true'),
      paste0("Raster WD Tolerance (meter)	, ", Raster_WD_Tolerance),
      paste0("Update Peak Every DT		, ", "false"),
      paste0("Expand Domain                   , ", "false"),
      paste0("Ignore Upstream			, ", 'true'),
      paste0("Upstream Reduction (meter)	, ",Upstream_Reduction)),
    paste0(Input,run_short_name,"_","control.csv")
  )

}


#' Create Rainfall Time Series File for CAFlood
#'
#' Generates a rainfall input file formatted for use with the CAFlood hydrodynamic model. The file includes rainfall intensity values over specified time steps and is saved to the input directory.
#'
#' @param Input Character. Directory where the rainfall file will be saved.
#' @param Event_name Character. Name or identifier of the rainfall event, used to name the output file.
#' @param timesteps Numeric vector. Time steps in seconds, representing the duration of each rainfall measurement interval.
#' @param rainfall Numeric vector. Rainfall intensity values (in mm/hr) corresponding to each timestep.
#'
#' @return Saves a formatted rainfall text file in the specified `Input` directory. This file can be used directly by the CAFlood model.
#'
#' @details
#' The rainfall file produced by this function is formatted according to CAFlood’s expected input structure. Each row typically includes a time and corresponding rainfall intensity value. The rainfall data must align with the model's timestep configuration.
#' Ensure that `timesteps` and `rainfall` are of equal length, and that rainfall values are consistent with the desired temporal resolution.
#'
#' @examples
#' \dontrun{
#' make_it_rain(
#'   Input = "C:/Projects/Input",
#'   Event_name = "rain_event_2021",
#'   timesteps = seq(0, 3600, by = 300),
#'   rainfall = c(0, 10, 15, 5, 0)
#' )
#' }
#'
#' @export
make_it_rain <- function(Input,Event_name,timesteps,rainfall){
  writeLines(
    c(paste0("Event Name,", Event_name),
      paste(c("Rain Intensity (mm/hr)",rainfall),sep=',',collapse = ","),
      paste(c("Time Stop (seconds)",timesteps),sep=',',collapse = ","),
      "Area (tlx tly brx bry),"),
    paste0(Input,Event_name,"_","rain",".csv")
  )
}


#' Create Inflow Time Series for CAFlood
#'
#' Generates the inflow configuration for the CAFlood model based on water level time series and location coordinates. The output is typically saved in the model's input folder.
#'
#' @param Input Character. Path to the input directory where the output inflow file will be saved.
#' @param Event_name Character. Name or identifier of the hydrological event, used for naming the output file.
#' @param timesteps Numeric vector. Sequence of time steps (e.g., in seconds) for which water level data is available.
#' @param water_level Numeric vector or matrix. Water level values (in meters) corresponding to each timestep. Should be the same length as `timesteps`, or a matrix with rows equal to `length(timesteps)` if multiple points.
#' @param tlx Numeric. X (longitude or easting) coordinate of the inflow point.
#' @param tly Numeric. Y (latitude or northing) coordinate of the inflow point.
#'
#' @return Writes an inflow configuration file (usually `.txt`) into the specified `Input` directory. This file is used by the CAFlood model to simulate inflows at specific locations.
#'
#' @details
#' This function formats and saves a structured time series file required by CAFlood to simulate inflows at specific grid cells or coordinates. The inflow location is determined by the `tlx` and `tly` coordinates, and the flow magnitude is derived from the `water_level` time series.
#' The file is typically named based on the `Event_name` and includes headers required by CAFlood.
#'
#' @examples
#' \dontrun{
#' make_it_flow(
#'   Input = "C:/Projects/Input",
#'   Event_name = "storm_event_01",
#'   timesteps = seq(0, 3600, by = 300),
#'   water_level = c(0, 0.1, 0.2, 0.15, 0),
#'   tlx = 345678.9,
#'   tly = 7654321.0
#' )
#' }
#'
#' @export
make_it_flow <- function(Input,Event_name,timesteps,water_level,tlx,tly){

  writeLines(
    c(paste0("Event Name,", Event_name),
      paste(c("Inflow (cumecs)",water_level),sep=',',collapse = ","),
      paste(c("Time Stop (seconds)",timesteps),sep=',',collapse = ","),
      paste(c("Zone (tlx tly w h)", c(tlx,tly,1,1)),sep=',',collapse = ",")),
    paste0(Input,Event_name,"_","WaterLevelBC",".csv")
  )

}



#' Run the CAFlood Model
#'
#' Executes the CAFlood hydrodynamic flood model using a specified controller file, input directory, and output directory.
#'
#' @param caflood_path Character. Full file path to the CAFlood executable (e.g., `"C:/Models/CAFlood/caflood.exe"`).
#' @param Input Character. Path to the input folder containing model setup files (e.g., DEM, rainfall, control points, etc.).
#' @param Output Character. Path to the output folder where simulation results will be stored.
#' @param controller Character. Path to the control file (e.g., `.txt`) that contains simulation settings and file references.
#'
#' @return This function runs for its side effects and does not return a value. It executes the CAFlood model, generating outputs in the specified `Output` directory. If implemented, it may log progress or return execution status.
#'
#' @details
#' The `make_it_run()` function is a simple wrapper to call the CAFlood executable via system command. The controller file should include all necessary parameters and file paths for the simulation.
#' Ensure that all required input files (DEM, rainfall, control points, etc.) are correctly referenced in the controller file and exist in the `Input` directory.
#'
#' @examples
#' \dontrun{
#' make_it_run(
#'   caflood_path = "C:/Models/CAFlood/caflood.exe",
#'   Input = "C:/Projects/Input",
#'   Output = "C:/Projects/Output",
#'   controller = "C:/Projects/Input/controller.txt"
#' )
#' }
#'
#' @export
make_it_run <- function(caflood_path, Input, Output, controller) {
  # Ensure full paths are used
  caflood_path <- normalizePath(caflood_path)
  exe_path <- file.path(caflood_path, "caflood.exe")

  # Build the arguments
  args <- c("/WCA2D", shQuote(Input), shQuote(controller), shQuote(Output))

  # Run the command using system2
  system2(command = exe_path, args = args, stdout = TRUE, stderr = TRUE)
}




#' Create and Snap Control Points to Stream Network
#'
#' Generates a shapefile of control (outlet) points based on longitude and latitude inputs, optionally snapping them to the nearest location on a stream network within a defined distance.
#'
#' @param Input Character. Directory path where the output control point shapefile will be saved.
#' @param Lon Numeric vector. Longitudes of the control points.
#' @param Lat Numeric vector. Latitudes of the control points.
#' @param id Character or numeric vector. Unique identifiers for each control point.
#' @param dem SpatRaster or RasterLayer. Digital Elevation Model used for elevation sampling and snapping logic.
#' @param dem_name Character. Name to be used in the output shapefile, typically associated with the DEM used.
#' @param Event_name Character. Name of the hydrological event or simulation, used for naming or metadata.
#' @param period Character or numeric. The simulation period or time interval associated with the event.
#' @param stream_network Spatial object (e.g., `sf` or `SpatVector`). Stream network used for snapping the control points.
#' @param snap_dist Numeric. Maximum distance (in meters) for snapping control points to the stream network.
#'
#' @return Saves a shapefile of the snapped control points to the specified `Input` directory. May also return a spatial object with snapped points if implemented.
#'
#' @details
#' This function is commonly used to define outlet or monitoring points for hydrological modeling. It ensures that the provided coordinates are aligned with the stream network by snapping them within a buffer distance defined by `snap_dist`. The resulting control points are compatible with CAFlood and other similar flood routing models.
#'
#' @examples
#' \dontrun{
#' make_control_points(
#'   Input = "C:/Projects/Input",
#'   Lon = c(-47.885, -47.891),
#'   Lat = c(-22.004, -22.008),
#'   id = c("P1", "P2"),
#'   dem = rast("C:/Projects/Input/DEM.tif"),
#'   dem_name = "dem_sp",
#'   Event_name = "storm_event_01",
#'   period = "2021-12-01",
#'   stream_network = st_read("C:/Projects/Hydro/streams.shp"),
#'   snap_dist = 30
#' )
#' }
#'
#' @export
make_control_points <- function(Input,
                                Lon,Lat,id,
                                dem,
                                dem_name,
                                Event_name,
                                period,
                                stream_network,
                                snap_dist){


  inputfile <- tempdir()

  dem <- raster::projectRaster(dem,crs="EPSG:31983")
  raster::writeRaster(dem,paste0(inputfile,"\\dem.tif"),overwrite=TRUE)
  raster::writeRaster(dem, paste0(Input,dem_name,".asc"),format = "ascii",overwrite=TRUE)
  streams <- st_read(stream_network)

  st_write(streams,paste0(inputfile,"\\streams.shp"),
           delete_dsn = T,
           delete_layer = T)

  wbt_fill_burn(
    dem = "dem.tif",
    streams = "streams.shp",
    output = "dem_filled.tif",
    wd = inputfile
  )

  wbt_d8_flow_accumulation(input = paste0(inputfile,"\\dem_filled.tif"),
                           output = paste0(inputfile,"\\D8FA.tif"))

  wbt_d8_pointer(dem = paste0(inputfile,"\\dem_filled.tif"),
                 output = paste0(inputfile,"\\D8pointer.tif"))

  ppoints <- data.frame(
    Lon = Lon,
    Lat = Lat
  )

  ppointsSP <- SpatialPoints(ppoints, proj4string = CRS("+init=epsg:31982"))
  ppointsSP <- st_as_sf(ppointsSP)

  st_write(ppointsSP,
           paste0(inputfile,"\\pourpoints.shp"),
           delete_dsn = T,
           delete_layer = T)


  threshold <- 2000

  wbt_extract_streams(flow_accum = paste0(inputfile,"\\D8FA.tif"),
                      output = paste0(inputfile,"\\raster_streams.tif"),
                      threshold = threshold)

  wbt_jenson_snap_pour_points(pour_pts = paste0(inputfile,"\\pourpoints.shp"),
                              streams = paste0(inputfile,"\\raster_streams.tif"),
                              output = paste0(inputfile,"\\snappedpp.shp"),
                              snap_dist = 100) #careful with this! Know the units of your data

  pp <- shapefile(paste0(inputfile,"\\snappedpp.shp"))

  coords <- coordinates(pp)

  ### Write VEL points
  writeLines(
    c(paste0("Time Plot Name		,", Event_name),
      "Physical Variable	, VEL",
      paste(c("Points Name",id),sep = ",",collapse = ','),
      paste(c("Points X Coo", coords[,1]) ,sep = ",",collapse = ','),
      paste(c("Points Y Coo", coords[,2])  ,sep = ",",collapse = ','),
      paste("Period (seconds)", period,sep = ",")),
    paste0(Input,dem_name,"_","VELpoints",".csv")
  )

  ### Write WL points
  writeLines(
    c(paste0("Time Plot Name		,", Event_name),
      "Physical Variable	, WL",
      paste(c("Points Name",id),sep = ",",collapse = ','),
      paste(c("Points X Coo", coords[,1]) ,sep = ",",collapse = ','),
      paste(c("Points Y Coo", coords[,2])  ,sep = ",",collapse = ','),
      paste("Period (seconds)", period,sep = ",")),
    paste0(Input,dem_name,"_","WLpoints",".csv")
  )


  writeLines(
    c(
      "Raster Grid Name	, Teste LiDAR 0.5 m Water Depth Raster Grid",
      paste("Physical Variable","WD",sep = ","),
      paste("Peak","true",sep = ","),
      paste("Period (seconds)","0",sep = ",")),
    paste0(Input,"WDraster",".csv")
  )

  writeLines(
    c(
      "Raster Grid Name	, Teste LiDAR 0.5 m Water Depth Raster Grid",
      paste("Physical Variable","VEL",sep = ","),
      paste("Peak","true",sep = ","),
      paste("Period (seconds)","0",sep = ",")),
    paste0(Input,"VELraster",".csv")

  )

}



#' Calculate Effective Precipitation Using the SCS Curve Number Method
#'
#' Computes effective precipitation (also known as runoff) from total precipitation using
#' the USDA Soil Conservation Service (SCS) Curve Number method.
#'
#' @param P Numeric vector or scalar. Total precipitation (mm).
#' @param CN Numeric vector or scalar. Curve Number, a parameter representing land use, soil type, and antecedent moisture condition (typically between 30 and 100).
#' @param lambda Numeric scalar. Initial abstraction coefficient (typically 0.2 as per standard SCS method).
#' @param SS Numeric scalar or vector. Optional. Pre-existing soil storage (mm), default is 0.
#'
#' @return A numeric vector or scalar representing effective precipitation (mm).
#'
#' @details The SCS Curve Number method estimates the amount of runoff generated from a rainfall event.
#' The potential maximum retention (S) is calculated from the Curve Number:
#' \deqn{S = (25400 / CN) - 254}
#' The initial abstraction (Ia), representing losses before runoff begins (e.g., infiltration, interception), is:
#' \deqn{Ia = \lambda \cdot S}
#' If \eqn{P > Ia}, the effective precipitation (Q) is calculated as:
#' \deqn{Q = \frac{(P - Ia)^2}{(P - Ia + S)}}
#' Otherwise, runoff is 0.
#'
#' @examples
#' effective_precipitation(P = 50, CN = 75, lambda = 0.2)
#' effective_precipitation(P = c(10, 30, 60), CN = 85, lambda = 0.05)
#'
#' @export
effective_precipitation <- function(P,CN,lambda,SS=0){
  #CN<-CN+(imperviousness)*(98-CN)
  S = 25400/CN-254
  Ia = lambda*S
  Prec.acc <- cumsum(P)-SS
  Prec.eff <- ifelse(Prec.acc <= Ia,0,(Prec.acc-Ia)^2/(Prec.acc-Ia+S))
  Prec.eff.diff <- c(0,diff(Prec.eff,1))
  return(Prec.eff.diff)
}

#' Run CAFlood Hydrological Simulation
#'
#' Executes the CAFlood model with specified parameters and input datasets for hydrodynamic flood modeling.
#'
#' @param param List. A list of model parameters, typically including calibration and simulation settings (e.g., `time_start`, `time_end`, `max_DT`, `roughness_global`, etc.).
#' @param caflood_path Character. Full path to the CAFlood executable or script to be run.
#' @param Input Character. Directory path containing input data files required for the simulation (e.g., DEM, rainfall, etc.).
#' @param Output Character. Directory where simulation output will be saved.
#' @param event_name Character. Name or ID of the rainfall/flood event being simulated; used for labeling outputs.
#' @param dem_name Character. Filename (without extension) of the DEM to be used, assumed to be in the `Input` directory.
#' @param path2 Character. Working directory path where temporary or auxiliary files may be stored or processed.
#' @param outlet_path Character. File path to a shapefile containing outlet or control points. Defaults to a specific São Paulo shapefile.
#' @param rain_path Character. File path to the rainfall input data (e.g., CSV or TXT format).
#' @param stream_network Character. File path to the stream network shapefile or raster used to delineate flow directions and connectivity.
#' @param snap_dist Numeric. Maximum distance (in meters) for snapping outlet points to the stream network.
#'
#' @return This function is typically run for its side effects. It writes simulation results to the `Output` directory and may return logs, paths, or a summary object if implemented.
#'
#' @details
#' This function prepares and runs the CAFlood model by integrating elevation data (DEM), rainfall time series, and river network structures. It also manages outlet snapping and may perform pre-processing for flow direction and accumulation.
#' The simulation results typically include water depth maps, hydrographs at outlet points, and diagnostics files.
#'
#' @examples
#' \dontrun{
#' param <- list(time_start = 0, time_end = 3600, max_DT = 5, roughness_global = 0.035)
#' cafloodr(
#'   param = param,
#'   caflood_path = "C:/Models/CAFlood/caflood.exe",
#'   Input = "C:/Projects/Input",
#'   Output = "C:/Projects/Output",
#'   event_name = "storm_event_01",
#'   dem_name = "dem_sp",
#'   path2 = "C:/Projects/Working",
#'   outlet_path = "C:/Projects/Shapefiles/outlets.shp",
#'   rain_path = "C:/Projects/Rain/rain_data.csv",
#'   stream_network = "C:/Projects/Hydro/stream.shp",
#'   snap_dist = 30
#' )
#' }
#'
#' @export
cafloodr <- function(param,
                     caflood_path,
                     Input,
                     Output,
                     event_name,
                     dem_name,
                     path2,
                     outlet_path = "C:/Projetos/01_DTI-A/02_Input/01_SP/01_Hidro_hora/estacoes.shp",
                     rain_path,
                     stream_network,
                     snap_dist
){

  alpha_fraction = param$alpha_fraction
  roughness_global = param$roughness_global
  slope_tolerance = param$slope_tolerance
  boundary_ele = param$boundary_ele
  Raster_WD_Tolerance = param$Raster_WD_Tolerance
  Upstream_Reduction = param$Upstream_Reduction

  Event_name <- paste0("E_",event_name)
  run_name <- paste(Event_name,dem_name,sep='_')
  run_short_name <- paste(Event_name,dem_name,sep="_")

  dem_path <- paste0(Input,path2,dem_name,".tif")

  stream_network_path <- paste0(stream_network,".shp")


  rain <- read.csv(rain_path)
  colnames(rain)[1] <- "timesteps"

  #rain$prec <- na_kalman(rain$prec)

  make_it_rain(Input = Input,
               Event_name = Event_name,
               timesteps = rain$timesteps,
               rainfall  = rain$prec*6)


  outlets <- st_read(outlet_path)
  outlets <- st_transform(outlets,crs = 31983)

  outlets <- outlets |>
    subset(ID %in% c(11,157,1000858,1000886))


  make_control_points(Input = Input,
                      Lon = st_coordinates(outlets)[,1],
                      Lat = st_coordinates(outlets)[,2],
                      id = outlets$ID,
                      dem = raster(dem_path,crs = "EPSG:31983"),
                      dem_name = dem_name,
                      Event_name = Event_name,
                      period = 600,
                      stream_network = stream_network_path,
                      snap_dist = snap_dist
  )

  options(scipen = 999)

  caflood_controller(
    Input,
    run_name = run_name,
    run_short_name,
    time_start = 0,
    time_end = max(rain$timesteps)*2,
    max_DT = 60,
    min_DT = 0.01,
    alpha_fraction = alpha_fraction,
    max_iterations = 100000000,
    roughness_global = roughness_global,
    ignore_WD = 0.0001,
    tolerance = 0.0001,
    slope_tolerance = slope_tolerance,
    boundary_ele = boundary_ele,
    dem_path = paste0(dem_name,".asc"),
    dem_name = dem_name,
    rain_path = paste0(Event_name,"_rain.csv"),
    flow_path = NULL,
    output_period = 300,
    Raster_WD_Tolerance = Raster_WD_Tolerance,
    Upstream_Reduction = Upstream_Reduction
  )

  make_it_run(caflood_path,
              Input,
              Output,
              controller = paste0(run_short_name,"_control.csv"))


}

#' Run CAFlood Hydrological Simulation Version 2
#'
#' Executes the CAFlood model with specified parameters and input datasets for hydrodynamic flood modeling.
#'
#'
#' @param param List. Model parameters including: `alpha_fraction`, `roughness_global`, `slope_tolerance`, `boundary_ele`, `Raster_WD_Tolerance`, `Upstream_Reduction`.
#' @param caflood_path Character. Path to the CAFlood executable.
#' @param Input Character. Directory containing model input files.
#' @param Output Character. Directory to store model output.
#' @param event_name Character. Name or ID of the event; used for file naming.
#' @param dem_name Character. Base name of the DEM file (without extension).
#' @param path2 Character. Subdirectory inside `Input` where the DEM is located.
#' @param outlet_path Character. Path to outlet shapefile (default: predefined path).
#' @param rain_path Character. Path to rainfall CSV file.
#' @param stream_network Character. Path to stream network shapefile (without `.shp` extension).
#' @param snap_dist Numeric. Max distance to snap outlets to stream network (meters).
#'
#' @return Logical. `TRUE` if the simulation was triggered successfully, otherwise `FALSE`.
#'
#' @export
cafloodr_A01 <- function(param,
                         caflood_path,
                         Input,
                         Output,
                         event_name,
                         dem_name,
                         path2,
                         outlet_path = "C:/Projetos/01_DTI-A/02_Input/01_SP/01_Hidro_hora/estacoes.shp",
                         rain_path,
                         stream_network,
                         snap_dist) {

  # --- Validate inputs ---
  stopifnot(is.list(param), file.exists(caflood_path), dir.exists(Input), dir.exists(Output))
  if (!file.exists(rain_path)) stop("Rainfall file not found.")
  if (!file.exists(outlet_path)) stop("Outlet shapefile not found.")
  if (!file.exists(paste0(stream_network, ".shp"))) stop("Stream network file not found.")

  dem_path <- file.path(Input, path2, paste0(dem_name, ".tif"))
  if (!file.exists(dem_path)) stop("DEM file not found at: ", dem_path)

  # --- Extract required parameters with defaults ---
  get_param <- function(name, default = NULL) {
    if (!name %in% names(param)) {
      warning(paste("Parameter", name, "not found in `param`. Using default:", default))
      return(default)
    }
    param[[name]]
  }

  alpha_fraction <- get_param("alpha_fraction", 0.1)
  roughness_global <- get_param("roughness_global", 0.035)
  slope_tolerance <- get_param("slope_tolerance", 0.5)
  boundary_ele <- get_param("boundary_ele", 500)
  Raster_WD_Tolerance <- get_param("Raster_WD_Tolerance", 0.01)
  Upstream_Reduction <- get_param("Upstream_Reduction", 0.5)

  # --- Setup names ---
  Event_name <- paste0("E_", event_name)
  run_name <- paste(Event_name, dem_name, sep = "_")
  run_short_name <- paste(Event_name, dem_name, sep = "_")

  # --- Process rainfall ---
  rain <- read.csv(rain_path)
  if (!all(c("prec") %in% names(rain))) stop("Rainfall data must contain a `prec` column.")
  rain$prec <- forecast::na_kalman(rain$prec)
  rain$timesteps <- seq(10, length.out = nrow(rain), by = 10) * 60

  make_it_rain(
    Input = Input,
    Event_name = Event_name,
    timesteps = rain$timesteps,
    rainfall = rain$prec * 6 # Convert mm/10min to mm/h
  )

  # --- Prepare outlet points ---
  outlets <- sf::st_read(outlet_path, quiet = TRUE) |>
    sf::st_transform(crs = 31983)

  make_control_points(
    Input = Input,
    Lon = sf::st_coordinates(outlets)[, 1],
    Lat = sf::st_coordinates(outlets)[, 2],
    id = outlets$ID,
    dem = terra::rast(dem_path),
    dem_name = dem_name,
    Event_name = Event_name,
    period = 600,
    stream_network = paste0(stream_network, ".shp"),
    snap_dist = snap_dist
  )

  # --- Create controller file ---
  caflood_controller(
    Input = Input,
    run_name = run_name,
    run_short_name = run_short_name,
    time_start = 0,
    time_end = max(rain$timesteps),
    max_DT = 60,
    min_DT = 0.01,
    alpha_fraction = alpha_fraction,
    max_iterations = 1e6,
    roughness_global = roughness_global,
    ignore_WD = 0.0001,
    tolerance = 0.0001,
    slope_tolerance = slope_tolerance,
    boundary_ele = boundary_ele,
    dem_path = paste0(dem_name, ".asc"),
    dem_name = dem_name,
    rain_path = paste0(Event_name, "_rain.csv"),
    flow_path = NULL,
    output_period = 300,
    Raster_WD_Tolerance = Raster_WD_Tolerance,
    Upstream_Reduction = Upstream_Reduction
  )

  # --- Run model ---
  controller_file <- file.path(Input, paste0(run_short_name, "_control.csv"))
  if (!file.exists(controller_file)) {
    warning("Controller file not created: ", controller_file)
    return(FALSE)
  }

  make_it_run(
    caflood_path = caflood_path,
    Input = Input,
    Output = Output,
    controller = controller_file
  )

  message("CAFlood simulation launched successfully.")
  return(TRUE)
}


