#' Countryside-SAR
#' @description
#' A function to perform a complete base SAR analysis or a countryside SAR (cSAR) analysis using binary species presence/absence data. It contains two analysis pathways: a nested "circles" approach, or a hierarchical "clusters" approach. By choosing the "circles" pathway, the function samples species data in stepwise increasing circles, while "clusters" assigns standardized squares to each sampling point and then groups them in clusters of increasing size based on proximity. The sampled data is then aggregated and used for the SAR or cSAR analyis, the latter includes habitat affinity values of the species groups to different habitat types. 
#'
#' @param data Binary species data table of the following structure: the first column is the location ID, second and third columns are longitude and latitude values of the sampling locations, columns 4 and onward contain binary presence/absence data of the species.
#' @param crs Coordinate reference system (CRS) of the sampling location.
#' @param method Parameter to select the sampling method: "circles" or "clusters". Pathway "circles" samples species data stepwise within a defined hull around all sampling points based on the parameter 'radius', beginning with a random starting point. Pathway "clusters" assigns each sampling point a standardized square with the extent of 'square_size'. The squares are then grouped according to 'cluster_sizes' based on their proximity to one another. The species count of all circles or clusters is then aggregated into one result table.
#' @param radius Defines the radius size as well as the total amount of circles for method: "circles". E.g. c(2000 * 1:10), sample for 10 circles with an extent of 2000 units each.
#' @param break_threshold Initiates break-protocol for method "circles" based on the proportion of the circular vector that lies inside the convex hull. E.g. break_threshold = 0.9 -> if less than 90 % of the circle lies inside of the hull, stop.
#' @param custom_hull Import a polygon hull for method "circles". If method = "clusters" the function will ignore the imported hull and auto-generate a hull instead. If custom_hull = NULL, the function auto-generates a hull for method "circles".
#' @param square_size The size of the initial square buffer for each sampling point for method "clusters".
#' @param cluster_sizes Numerical vector that defines the amount of levels as well as the amount of sampling points within each cluster of each level for the hierarchical "clusters" approach. E.g. c(1, 4, 16, 64, 256) -> 5 levels for 256 sampling points; first level 256 clusters, second level 64 clusters, etc.
#' @param habitat Land-use raster of the sampling location, a .tif file.
#' @param habitat_names Character vector with the names of the land-use types 'habitat' land-use raster.
#' @param classification Species classification file. A table with a first column for species names and the following columns as binary (0/1) group indicators.
#' @param groups Character vector to define the group columns used for analyis. If NULL = use all columns except first.
#' @param seed Optional seed for reproducibility.
#' @param transform_to_utm If TRUE, transforms geographic coordinates (long/lat) to UTM projection. It automatically detects the appropriate UTM zone based on mean longitude. If data was sampled in polar, equatorial or very large regions, use an appropriate projection instead of UTM.
#' @param target_crs Optional numeric EPSG code for coordinate transformation. If provided it overrides the automatic UTM zone detection. Use this when you need a specific national or local projection instead of UTM.
#' @param warn_projection Defaults to TRUE, will warn the user about possible projection issues.
#' @param n_runs Option to run the analyis multiple times, default value = 1. Only available for method "circles" as "clusters" is deterministic, multiple runs produce the same result.
#'
#' @return A list containing the method used, the number of runs n_runs, the sampling data of each run, the sf-transformed input data with an added geometry column as well as information about the convex hull. The sampling data of each run contains a results_table of aggregated habitat area data and species richness data. It further contains the SAR analysis result and linear model summary as well as geometry data of each circle or cluster and species data within each circle or cluster level.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- countryside_sar(
#'   data = mydata,
#'   crs = 3763,
#'   method = "circles",
#'   radius = 2000 * 1:10,
#'   habitat = myraster,
#'   habitat_names = c("Forest", "Agriculture", "Shrubland"),
#'   habitat_codes = 1:3,
#'   classification = myclassif,
#'   groups = "Forest",
#'   seed = 123
#' )
#' }
countryside_sar <- function(
    data,
    crs,
    method = c("circles", "clusters"),
    # Circles parameters
    radius = NULL,
    break_threshold = 0.5,
    custom_hull = NULL,
    # Clusters parameters
    square_size = NULL,
    cluster_sizes = NULL,
    # Habitat raster
    habitat = NULL,
    habitat_names = NULL,
    habitat_codes = NULL,
    # Species classification
    classification = NULL,
    groups = NULL,
    seed = NULL,
    # Coordinate transformation options
    transform_to_utm = FALSE,
    target_crs = NULL,
    warn_projection = TRUE,
    n_runs = 1
) {

  # create vector for habitat codes of the raster
  habitat_codes <- seq_along(habitat_names)
    
  # Optional: Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)

  # If needed: Coordinate transformation
  if (transform_to_utm) {

    # Store original data
    original_data <- data

    # Check if coordinates are already projected
    if (any(abs(data$long) > 180) || any(abs(data$lat) > 90)) {
      stop("Coordinates appear to already be in a projected CRS. ",
           "Set transform_to_utm = FALSE or check your data.")
    }

    # Create sf object with user-provided CRS
    points_sf_original <- sf::st_as_sf(data, coords = c("long", "lat"), crs = crs)

    # Determine target CRS
    if (is.null(target_crs)) {
      # Auto-detect UTM zone
      mean_lon <- mean(data$long, na.rm = TRUE)
      mean_lat <- mean(data$lat, na.rm = TRUE)
      utm_zone <- floor((mean_lon + 180) / 6) + 1
      utm_zone <- max(1, min(60, utm_zone))
      hemisphere <- ifelse(mean_lat >= 0, "N", "S")

      if (hemisphere == "N") {
        target_crs <- 32600 + utm_zone
      } else {
        target_crs <- 32700 + utm_zone
      }

      # Projection limitation warnings:
      if (warn_projection) {
        # Check if data spans UTM zones
        lon_range <- range(data$long, na.rm = TRUE)
        utm_zone_min <- floor((min(lon_range) + 180) / 6) + 1
        utm_zone_max <- floor((max(lon_range) + 180) / 6) + 1

        if (utm_zone_max - utm_zone_min >= 1) {
          warning("Data spans multiple UTM zones (zones ", utm_zone_min, " to ", utm_zone_max, ").\n",
                  "  Using UTM zone ", utm_zone, " (EPSG:", target_crs, ") for the entire dataset.\n",
                  "  This approach may cause distortion at the edges. For better accuracy, please consider:\n",
                  "  - Using a national projection\n",
                  "  - Specifying target_crs manually with an appropriate projection")
        }

        # Check for high latitudes
        if (max(abs(data$lat)) > 70) {
          warning("Data includes high latitudes (>70°). UTM projection may have significant distortion.\n",
                  "  Consider using a polar projection for better accuracy.")
        }

        # Check for very large extent
        lon_span <- diff(range(data$long))
        lat_span <- diff(range(data$lat))
        if (lon_span > 20 || lat_span > 20) {
          warning("Data covers a large area (>20°). UTM projection may have significant distortion.\n",
                  "  Consider using an equal-area projection for better accuracy.")
        }
      }

    } else {
      # User specified target_crs - validate it's numeric
      if (!is.numeric(target_crs)) {
        stop("target_crs must be a numeric EPSG code (e.g., 3763 for Portugal TM06)")
      }

      # Optional warning if using non-UTM projection
      if (warn_projection && !(target_crs %in% c(32601:32660, 32701:32760))) {
        warning("You specified a non-UTM projection (EPSG:", target_crs, ").\n",
                "  Ensure this projection is appropriate for your study area.")
      }
    }

    # Transform coordinates
    points_sf_transformed <- sf::st_transform(points_sf_original, crs = target_crs)

    # Extract transformed coordinates
    coords <- sf::st_coordinates(points_sf_transformed)

    # Rebuild data frame with transformed coordinates
    data <- data.frame(
      locationID = data[[1]],
      long = coords[, "X"],
      lat = coords[, "Y"],
      data[, 4:ncol(data)]
    )

    # Update CRS to the target
    crs <- target_crs
  }

  #---------------------------- 1. Input validation ----------------------------
  method <- match.arg(method)

  # dataframe
  if (!"long" %in% names(data) || !"lat" %in% names(data))
    stop("Data must contain 'long' and 'lat' columns.")

  if (ncol(data) < 4)
    stop("Data must have at least 4 columns: locationID, long, lat, and binary species data.")

  # crs
  if (is.null(crs)) stop("Coordinate reference system (crs) must be provided.")

  # n_runs
  if (!is.numeric(n_runs) || n_runs < 1 || n_runs != round(n_runs)) {
    stop("n_runs must be a positive integer")
  }

  # Method specific validation
  if (method == "circles") {
    if (is.null(radius)) stop("'radius' is required for method = 'circles'.")
    if (!is.numeric(radius) || any(radius <= 0))
      stop("'radius' must be a positive numeric vector.")

    # Validate custom hull if provided
    if (!is.null(custom_hull)) {
      if (!inherits(custom_hull, c("sf", "sfc"))) {
        stop("custom_hull must be an sf or sfc object")
      }

      # Is hull CRS = points CRS
      hull_crs <- sf::st_crs(custom_hull)
      if (hull_crs != sf::st_crs(crs)) {
        warning("Custom hull has different CRS than data. Attempting to transform...")
        custom_hull <- sf::st_transform(custom_hull, crs = crs)
      }

      # make sure custom hull is a polygon
      if (!all(sf::st_geometry_type(custom_hull) %in% c("POLYGON"))) {
        stop("custom_hull must be a polygon geometry")
      }
    }
  } else {
    if (is.null(square_size) || is.null(cluster_sizes))
      stop("'square_size' and 'cluster_sizes' are required for method = 'clusters'.")

    if (!is.numeric(square_size) || square_size <= 0)
      stop("'square_size' must be a positive number.")

    if (!is.numeric(cluster_sizes) || any(cluster_sizes <= 0))
      stop("'cluster_sizes' must be a positive numeric vector.")

    # Warning if custom hull is provided for clusters method
    if (!is.null(custom_hull)) {
      warning("custom_hull is ignored for method = 'clusters'")
    }

    # For clusters n_runs > 1 doesn't make sense
    if (n_runs > 1) {
      warning("n_runs > 1 is ignored for method = 'clusters' (clustering is deterministic)")
      n_runs <- 1
    }
  }

  # Habitat specific validation
  if (is.null(habitat) || is.null(habitat_names) || is.null(habitat_codes))
    stop("'habitat', 'habitat_names', and 'habitat_codes' are required.")

  if (length(habitat_names) != length(habitat_codes))
    stop("Length of 'habitat_names' must equal length of 'habitat_codes'.")

  if (is.null(classification))
    stop("'classification' of species groups must be provided.")

  if (ncol(classification) < 2)
    stop("'classification' must have at least two columns: species names and one group column.")

  species_name_col <- names(classification)[1]
  group_cols <- if (is.null(groups)) names(classification)[-1] else groups

  missing_cols <- setdiff(group_cols, names(classification))
  if (length(missing_cols) > 0)
    stop("The following group columns are not in classification: ",
         paste(missing_cols, collapse = ", "))

  #---------------------------- 2. Species group construction ------------------
  # Species data
  sp_cols <- 4:ncol(data) # species data starts at column 4
  sp_names <- colnames(data)[sp_cols]
  n_sp <- length(sp_cols) # number of species
  clean_sp <- tolower(trimws(sp_names)) # clean species names for uniform formatting (lowercase, no space at the end or beginning)

  # Classification data
  classif_sp <- classification[[species_name_col]]
  clean_classif <- tolower(trimws(classif_sp))

  # Lists for species groups and group names
  sp_groups <- list()
  grp_names <- character()

  for (col in group_cols) {
    vals <- classification[[col]]

    if (is.logical(vals)) {
      idx <- which(vals)
    } else if (is.numeric(vals)) {
      idx <- which(vals == 1)
    } else {
      stop("Group column '", col, "' must be logical or numeric (0/1).")
    }

    sp_clean <- clean_classif[idx]
    positions <- na.omit(match(sp_clean, clean_sp))

    if (length(positions) > 0) {
      sp_groups <- c(sp_groups, list(positions))
      grp_names <- c(grp_names, paste0("Sp_", col))
    }
  }

  for (i in seq_along(sp_groups)) {
    if (any(sp_groups[[i]] < 1 | sp_groups[[i]] > n_sp))
      stop("Error: group indices out of range.")
  }

  #---------------------------- 3. Helper functions ----------------------------

  filter_points_in_expanding_circles <- function(points_sf,
                                                 radius_vector,
                                                 convex_hull,
                                                 break_threshold) {
    # Randomly selected point of sampling start
    selected_point <- points_sf[sample(1:nrow(points_sf), 1), ]
    points_within_circles <- list()

    # Loop through the radius_vector
    for (radius in radius_vector) {
      circle <- sf::st_geometry(sf::st_buffer(selected_point, dist = radius))
      intersection <- sf::st_intersection(circle, convex_hull)
      circle_area <- as.numeric(sf::st_area(circle))
      intersection_area <- as.numeric(sf::st_area(intersection))

      if (intersection_area / circle_area < break_threshold) {
        break
      }

      points_in_circle <- points_sf[sf::st_intersects(points_sf, circle, sparse = FALSE), ]
      points_within_circles[[paste0("radius_", radius)]] <-
        list(points = points_in_circle, circle = circle)
    }
    return(points_within_circles)
  }


  summarize_samples <- function(samples,
                                polygons,
                                habitat_raster,
                                habitat_names,
                                habitat_values,
                                species_groups,
                                species_group_names)
  {
    # Calculate the "Other" code as max habitat value + 1 for reclass
    other_code <- max(habitat_values) + 1

    # Add "Other" to habitat_names and habitat_values (if not already present)
    if (!"Other" %in% habitat_names) {
      habitat_names <- c(habitat_names, "Other")
      habitat_values <- c(habitat_values, other_code)
    }

    # Initialize an empty data frame for the results (added Polygon_Area for comparison with calculated Area_Total)
    results_df <- data.frame(matrix(ncol = length(habitat_names)+
                                      length(species_group_names)+3,
                                    nrow = 0))
    colnames(results_df) <- c(habitat_names, "Area_Total", species_group_names, "Sp_Total", "Polygon_Area")

    # Iterate over each area sample (e.g. group of sites)
    for (i in seq_along(samples))
    {
      sample <- samples[[i]]
      polygon <- polygons[[i]]

      # Calculate polygon area
      polygon_area <- sf::st_area(polygon)

      # Crop raster to polygon extent
      habitat_cropped <- terra::crop(habitat_raster, terra::vect(polygon))

      # Reclassify NA cells to "Other" code
      habitat_cropped[is.na(habitat_cropped)] <- other_code

      # mask raster to polygon extent
      habitat_masked <- terra::mask(habitat_cropped, terra::vect(polygon))

      # Calculate the area of each habitat type
      habitat_df <- terra::freq(habitat_masked, bylayer = FALSE)

      # Ensure all possible values are included (now includes other_code)
      all_values_df <- data.frame(value = habitat_values, count = 0)

      # Merge with actual frequency data, replacing 0 where missing
      habitat_df <- merge(all_values_df, habitat_df, by = "value", all.x = TRUE)

      # Fill NA counts with 0
      habitat_df$count <- ifelse(is.na(habitat_df$count.y), 0, habitat_df$count.y)

      # Drop unnecessary column
      habitat_df <- habitat_df[, c("value", "count")]

      habitat_df$area <- habitat_df$count * terra::res(habitat_raster)[1] * terra::res(habitat_raster)[2]

      # Store the results
      results_df[i, seq_along(habitat_names)] <- habitat_df$area
      results_df[i, length(habitat_names)+1] <- sum(results_df[i, seq_along(habitat_names)])

      # Add the polygon area
      results_df[i, "Polygon_Area"] <- as.numeric(polygon_area)

      # Subset species occurrences for each group
      total_species <- 0
      for (k in seq_along(species_groups))
      {
        group_species <- sf::st_drop_geometry(sample)[, species_groups[[k]] + 1]
        species_present <- colSums(group_species > 0)
        num_species <- sum(species_present > 0)
        results_df[i, length(habitat_names) + 1 + k] <- num_species
        total_species <- total_species + num_species
      }

      # Store the total number of species
      results_df[i, length(habitat_names) + length(species_groups) + 2] <- total_species
    }

    return(results_df)
  }


  create_squares <- function(points_sf,
                             width)
  {
    if (!inherits(points_sf, "sf") || !inherits(sf::st_geometry(points_sf), "sfc_POINT"))
      stop("Input must be an sf object with point geometries.")
    # Calculate half-width (to shift the square corners)
    half_width <- width / 2
    squares_sf <- sf::st_as_sf(sf::st_buffer(points_sf, dist = half_width, endCapStyle = "SQUARE"))
    return(squares_sf)
  }


  filter_points_in_clusters <- function(points_sf,
                                        squares_sf,
                                        cluster_size_vector)
  {
    npoints <- nrow(points_sf)
    n_clusters_vector <- npoints %/% cluster_size_vector
    n_clusters_vector[n_clusters_vector == 0] <- 1  # when whole landscape, npoints < cluster_size
    points_within_clusters <- list()

    for (i in seq_along(cluster_size_vector)) #here
    {
      n_clusters <- n_clusters_vector[i] #here
      size_val <- cluster_size_vector[i] #here

      if (n_clusters == npoints)
      {
        points_in_clusters <- split(points_sf, 1:npoints)
        clusters_convex_hulls <- split(sf::st_geometry(squares_sf), 1:npoints)

      } else
      {
        coords <- sf::st_coordinates(points_sf)
        kmeans_result <- kmeans(coords, centers = n_clusters, iter.max = 100, nstart = 25)
        cluster_assignments <- kmeans_result$cluster

        points_in_clusters <- list()
        squares_in_clusters <- list()

        for (c in 1:n_clusters)
        {
          cluster_idx <- which(cluster_assignments == c)
          points_in_clusters[[c]] <- points_sf[cluster_idx, ]
          squares_in_clusters[[c]] <- squares_sf[cluster_idx, ]
        }

        clusters_convex_hulls <- list()
        for (j in seq_along(squares_in_clusters))
        {
          n_sq <- nrow(squares_in_clusters[[j]])

          if (n_sq == 0)
          {
            warning("Cluster ", j, " has zero squares – skipping")
            next
          }

          merged <- sf::st_union(squares_in_clusters[[j]])
          hull <- sf::st_convex_hull(merged)
          clusters_convex_hulls[[j]] <- hull
        }
      }

      level_name <- paste0("size_", size_val)
      points_within_clusters[[level_name]] <-
        list(points = points_in_clusters, chulls = clusters_convex_hulls)
    }

    return(points_within_clusters)
  }


  # ---- SAR analysis function ----
  analyze_sar <- function(results_table,
                          method_used)
  {
    # Initialize results and check validity for further analysis
    sar_results <- list()
    sar_results$valid <- FALSE
    sar_results$message <- ""

    # Check if there is any data for analysis
    if (nrow(results_table) == 0) {
      sar_results$message <- "No data for SAR analysis"
      return(sar_results)
    }

    # Clean species data: filter out rows with zero/NA species
    valid_rows <- results_table$Sp_Total > 0 & !is.na(results_table$Sp_Total) &
      is.finite(results_table$Sp_Total)

    if (sum(valid_rows) < 2) {
      sar_results$message <- paste("Insufficient non-zero samples for SAR analysis (need at least 2, have",
                                   sum(valid_rows), ")")
      return(sar_results)
    }

    # Use only valid rows
    valid_table <- results_table[valid_rows, ]

    # Make sure area data is valid for analysis
    area_valid <- valid_table$Area_Total > 0 & !is.na(valid_table$Area_Total) &
      is.finite(valid_table$Area_Total)

    if (sum(area_valid) < 2) {
      sar_results$message <- "Insufficient valid area values for SAR analysis"
      return(sar_results)
    }

    valid_table <- valid_table[area_valid, ]

    # Log-log transformation
    log_area <- log(valid_table$Area_Total)
    log_sp <- log(valid_table$Sp_Total)

    # Check for infinite values
    if (any(!is.finite(log_area)) || any(!is.finite(log_sp))) {
      sar_results$message <- "Infinite values in log-transformed data"
      return(sar_results)
    }

    # Store transformed data
    sar_results$log_area <- log_area
    sar_results$log_sp <- log_sp
    sar_results$n_samples_used <- nrow(valid_table)
    sar_results$valid <- TRUE

    # Fit linear model with error handling
    tryCatch( # continue calculation despite error
    {
      sar_results$lm_model <- lm(log_sp ~ log_area)
      sar_results$lm_summary <- summary(sar_results$lm_model)
      sar_results$message <- "SAR analysis completed successfully"
    }, error = function(e)
      {
      sar_results$message <- paste("Error fitting linear model:", e$message)
      sar_results$valid <- FALSE
      }
    )

    # Store original data used
    sar_results$sar_data <- valid_table

    return(sar_results)
  }


  #---------------------------- 4. Main processing -----------------------------
  points_sf <- sf::st_as_sf(data,
                            coords = c("long", "lat"),
                            crs = crs)

  if (method == "circles")
  {
    # Use custom hull if provided, otherwise create from points
    if (!is.null(custom_hull))
    {
      convex_hull <- custom_hull
    } else
      {
      convex_hull <- sf::st_convex_hull(sf::st_union(points_sf))
      }

    # Run multiple iterations for circles method
    runs <- list()

    for (run in 1:n_runs)
    {
      samples <- filter_points_in_expanding_circles(points_sf,
                                                    radius,
                                                    convex_hull,
                                                    break_threshold)
      polygons <- lapply(samples, `[[`, "circle")
      samples_points <- lapply(samples, `[[`, "points")

      results_table <- summarize_samples(
        samples = samples_points,
        polygons = polygons,
        habitat_raster = habitat,
        habitat_names = habitat_names,
        habitat_values = habitat_codes,
        species_groups = sp_groups,
        species_group_names = grp_names
      )

      # Perform SAR analysis
      sar_analysis <- analyze_sar(results_table,
                                  method)

      # Store results
      runs[[run]] <- list(
        run = run,
        results_table = results_table,
        sar_analysis = sar_analysis,
        samples = samples,
        polygons = polygons
      )
    }

    # results
    res <- list(
      method = method,
      n_runs = n_runs,
      runs = runs,
      points_sf = points_sf,
      convex_hull = convex_hull,
      hull_source = ifelse(is.null(custom_hull), "derived", "custom")
    )

  } else {
    # Clusters method (deterministic, always single run)
    squares_sf <- create_squares(points_sf,
                                 square_size)

    samples <- filter_points_in_clusters(points_sf,
                                         squares_sf,
                                         cluster_sizes)

    samples_points <- lapply(samples, `[[`, "points")
    clusters_chulls <- lapply(samples, `[[`, "chulls")

    # Flatten the nested lists for summarize_samples
    samples_flat <- unlist(samples_points, recursive = FALSE)
    polygons_flat <- unlist(clusters_chulls, recursive = FALSE)

    results_table <- summarize_samples(
      samples = samples_flat,
      polygons = polygons_flat,
      habitat_raster = habitat,
      habitat_names = habitat_names,
      habitat_values = habitat_codes,
      species_groups = sp_groups,
      species_group_names = grp_names
    )

    # Perform SAR analysis
    # if or both for normal sar and csar analysis
    sar_analysis <- analyze_sar(results_table,
                                method)

    # result table
    res <- list(
      method = method,
        results_table = results_table,
        sar_analysis = sar_analysis,
        samples = samples,
        squares_sf = squares_sf,
        clusters_chulls = clusters_chulls
      )
      points_sf = points_sf
  }

  return(res)
}
