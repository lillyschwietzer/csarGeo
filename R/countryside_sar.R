#' Countryside-SAR
#' @description
#' A function to perform a complete classic SAR analysis or a countryside SAR (cSAR) analysis. It contains two analysis pathways: a nested "circles" approach, or a hierarchical "clusters" approach. Using "circles", the function samples species data in step wise increasing circles, while "clusters" groups it in clusters of increasing size based on their proximity to one another. The sampled data is then aggregated and used for the SAR or cSAR analysis, the latter includes habitat affinity values of the species groups to different habitat types.
#'
#' @details
#' The function implements two complementary sampling approaches:
#'
#' \strong{Circles method:}
#' Starting from a randomly selected point, circles expand outward at radii
#' defined by \code{radius}. Sampling stops when the proportion of circle area
#' falling within the convex hull drops below \code{break_threshold}. This
#' approach is stochastic; use \code{seed} and \code{n_runs} for reproducibility.
#'
#' \strong{Clusters method:}
#' Each sampling point is buffered with a square of size \code{square_size}.
#' These squares are then grouped hierarchically using k-means clustering
#' based on spatial proximity. The \code{cluster_sizes} vector defines the
#' target number of points per cluster at each level.
#'
#' For a detailed explanation, please colsult the vignette of the csarGeo package.
#'
#' @param data A data frame containing binary species data with the following
#'   required columns:
#'   \itemize{
#'     \item Column 1: Location ID (character or numeric)
#'     \item Column 2: Longitude or X coordinate (numeric). Can be decimal
#'       degrees (e.g., -8.5) or projected coordinates (e.g., 500000) depending
#'       on the CRS.
#'     \item Column 3: Latitude or Y coordinate (numeric). Can be decimal
#'       degrees (e.g., 40.2) or projected coordinates (e.g., 4500000) depending
#'       on the CRS.
#'     \item Columns 4+: Species presence/absence (0 = absent, 1 = present)
#'   }
#' @param method Character string specifying the sampling method. Options:
#'   \itemize{
#'     \item \code{"circles"}: Samples species data within expanding circles
#'       from a random starting point. Sampling is limited by a convex hull defined by
#'       \code{custom_hull}.
#'     \item \code{"clusters"}: Assigns each sampling point a square of size
#'       \code{square_size}.
#'   }
#' @param radius Numeric vector defining the radii when \code{method = "circles"}. E.g. \code{radius = 2000 * 1:10}, sample in 10 circles with an extent of 2000 units each.
#' @param break_threshold Numeric value ranging from 0 - 1, defines break-protocol-sensitivity for \code{method = "circles"} based on the proportion of the circular vector that lies inside the convex hull. E.g. \code{break_threshold = 0.9} -> if less than 90 % of the circle lies inside of the hull, stop sampling. Defaults to 0.5.
#' @param custom_hull Optional polygon hull import (sf or sfc object) for \code{method "circles"} to define the study area boundary. If \code{custom_hull = NULL}, the function auto-generates a hull for \code{method = "circles"}. For \code{method = "clusters"} the function will ignore the imported hull and auto-generate a hull instead.
#' @param square_size Numeric value defining the size of the square buffer created around each sampling point for \code{method = "clusters"}. Required if \code{method = "clusters"}.
#' @param cluster_sizes Numeric vector defining the hierarchical clustering levels for \code{method = "clusters"}. Each value specifies the approximate number of points per cluster at a given level and should therefore be divisors of the total number of points. The length of \code{cluster_sizes} defines the number of clustering levels. E.g. \code{cluster_sizes = c(1, 4, 16, 64, 256)}, creates 5 levels of 256, 64, 16, 4 and 1 clusters each. Required if \code{method = "clusters"}.
#' @param habitat Land-use raster input, SpatRaster or file path to a raster file (e.g. a .tif) containing land-use/land-cover classification of the study area. Categorical values in the raster should correspond to the habitat types defined in \code{habitat_names}.
#' @param habitat_names Character vector specifying the names of habitat types corresponding to the values in the \code{habitat} raster, acts as a legend to the land-use raster. E.g. , \code{habitat_names = c("Forest", "Agriculture", "Shrubland")} if raster values 1, 2, 3 represent these classes.
#' @param classification A data frame defining species classifications with the following structure:
#'   \itemize{
#'     \item Column 1: Species name (must match the name used in \code{data})
#'     \item Column 2+: Species group (binary, 0 = doesn't belong to group, 1 = belongs to group, multiple groups per species possible)
#'   }
#' @param groups Character vector specifying which columns from \code{classification} to use as species groups. If \code{groups = NULL} (default), all group columns are used for analysis.
#' @param seed Optional integer for reproducibility, defaults to \code{NULL}.
#' @param transform_to_utm transforms geographic coordinates (longitude/latitude) to UTM projection. If data was sampled in polar, equatorial or across very large regions, use an appropriate projection with \code{target_crs} instead of UTM to avoid distortion. Defaults to \code{FALSE}.
#' @param n_runs Integer value, number of iterations for \code{method = "circles"}. Each run uses a different random starting point for the expanding circles. For \code{method = "clusters"}, this parameter is ignored as clustering is deterministic. Default value \code{n_runs = 1}.
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
#'   classification = myclassif,
#'   groups = "Forest",
#'   seed = 123
#' )
#' }
countryside_sar <- function(
    data,
    method = c("circles", "clusters"),
    radius = NULL,
    break_threshold = 0.5,
    custom_hull = NULL,
    square_size = NULL,
    cluster_sizes = NULL,
    habitat = NULL,
    habitat_names = NULL,
    classification = NULL,
    groups = NULL,
    seed = NULL,
    n_runs = 1
) {

  habitat_codes <- seq_along(habitat_names)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(habitat)) stop("Habitat raster with CRS information must be provided.")
  crs <- terra::crs(habitat)
  if (is.na(crs) || crs == "") stop("Raster has no CRS information.")

  #---------------------------- 1. Input validation ----------------------------
  method <- match.arg(method)

  if (!"long" %in% names(data) || !"lat" %in% names(data))
    stop("Data must contain 'long' and 'lat' columns.")
  if (ncol(data) < 4)
    stop("Data must have at least 4 columns: locationID, long, lat, and binary species data.")
  if (is.null(crs)) stop("No Coordinate reference system (crs) provided.")
  if (!is.numeric(n_runs) || n_runs < 1 || n_runs != round(n_runs))
    stop("n_runs must be a positive integer")

  if (method == "circles") {
    if (is.null(radius)) stop("'radius' is required for method = 'circles'.")
    if (!is.numeric(radius) || any(radius <= 0))
      stop("'radius' must be a positive numeric vector.")
    if (!is.null(custom_hull)) {
      if (!inherits(custom_hull, c("sf", "sfc")))
        stop("custom_hull must be an sf or sfc object")
      hull_crs <- sf::st_crs(custom_hull)
      if (hull_crs != sf::st_crs(crs)) {
        warning("Custom hull has different CRS than data. Attempting to transform...")
        custom_hull <- sf::st_transform(custom_hull, crs = crs)
      }
      if (!all(sf::st_geometry_type(custom_hull) %in% c("POLYGON")))
        stop("custom_hull must be a polygon geometry")
    }
  } else {
    if (is.null(square_size) || is.null(cluster_sizes))
      stop("'square_size' and 'cluster_sizes' are required for method = 'clusters'.")
    if (!is.numeric(square_size) || square_size <= 0)
      stop("'square_size' must be a positive number.")
    if (!is.numeric(cluster_sizes) || any(cluster_sizes <= 0))
      stop("'cluster_sizes' must be a positive numeric vector.")
    if (!is.null(custom_hull))
      warning("custom_hull is ignored for method = 'clusters'")
    if (n_runs > 1) {
      warning("n_runs > 1 is ignored for method = 'clusters' (clustering is deterministic)")
      n_runs <- 1
    }
  }

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

  #---------------------------- 3. Helper functions ----------------------------

  if (method == "circles") {

    filter_points_in_expanding_circles <- function(points_sf,
                                                   radius_vector,
                                                   convex_hull,
                                                   break_threshold) {
      selected_point <- points_sf[sample(1:nrow(points_sf), 1), ]
      points_within_circles <- list()

      for (radius in radius_vector) {
        circle <- sf::st_geometry(sf::st_buffer(selected_point, dist = radius))
        intersection <- sf::st_intersection(circle, convex_hull)
        circle_area <- as.numeric(sf::st_area(circle))
        intersection_area <- as.numeric(sf::st_area(intersection))

        if (intersection_area / circle_area < break_threshold) break

        points_in_circle <- points_sf[sf::st_intersects(points_sf, circle, sparse = FALSE), ]
        points_within_circles[[paste0("radius_", radius)]] <-
          list(points = points_in_circle, circle = circle)
      }
      return(points_within_circles)
    }

  } else {

    create_squares <- function(points_sf, width) {
      if (!inherits(points_sf, "sf") || !inherits(sf::st_geometry(points_sf), "sfc_POINT"))
        stop("Input must be an sf object with point geometries.")
      half_width <- width / 2
      squares_sf <- sf::st_as_sf(sf::st_buffer(points_sf, dist = half_width, endCapStyle = "SQUARE"))
      return(squares_sf)
    }

    filter_points_in_clusters <- function(points_sf, squares_sf, cluster_size_vector) {
      npoints <- nrow(points_sf)
      n_clusters_vector <- npoints %/% cluster_size_vector
      n_clusters_vector[n_clusters_vector == 0] <- 1
      points_within_clusters <- list()

      for (i in seq_along(cluster_size_vector)) {
        n_clusters <- n_clusters_vector[i]
        size_val   <- cluster_size_vector[i]

        if (n_clusters == npoints) {
          points_in_clusters    <- split(points_sf, 1:npoints)
          clusters_convex_hulls <- split(sf::st_geometry(squares_sf), 1:npoints)
        } else {
          coords             <- sf::st_coordinates(points_sf)
          kmeans_result      <- kmeans(coords, centers = n_clusters, iter.max = 100, nstart = 25)
          cluster_assignments <- kmeans_result$cluster
          points_in_clusters  <- list()
          squares_in_clusters <- list()

          for (c in 1:n_clusters) {
            cluster_idx            <- which(cluster_assignments == c)
            points_in_clusters[[c]]  <- points_sf[cluster_idx, ]
            squares_in_clusters[[c]] <- squares_sf[cluster_idx, ]
          }

          clusters_convex_hulls <- list()
          for (j in seq_along(squares_in_clusters)) {
            n_sq <- nrow(squares_in_clusters[[j]])
            if (n_sq == 0) { warning("Cluster ", j, " has zero squares – skipping"); next }
            merged <- sf::st_union(squares_in_clusters[[j]])
            hull   <- sf::st_convex_hull(merged)
            clusters_convex_hulls[[j]] <- hull
          }
        }

        level_name <- paste0("size_", size_val)
        points_within_clusters[[level_name]] <-
          list(points = points_in_clusters, chulls = clusters_convex_hulls)
      }
      return(points_within_clusters)
    }

    extract_species_positions <- function(species_habitat_matrix, species_site_matrix) {

      # Get habitat group names (all columns except species name column)
      habitat_names <- colnames(species_habitat_matrix[, -1])

      # Get species column names from site matrix excluding locationID (col 1)
      # species_site_matrix still has locationID as col 1
      site_species_cols <- colnames(sf::st_drop_geometry(species_site_matrix))[-1]

      habitat_positions <- list()

      for (habitat in habitat_names) {

        # Use [[1]] and [[habitat]] to extract as plain vectors (works for tibbles too)
        species_in_habitat <-
          species_habitat_matrix[[1]][species_habitat_matrix[[habitat]] == 1]

        # Find column positions relative to species-only columns
        species_positions <- which(site_species_cols %in% species_in_habitat)

        habitat_positions[[habitat]] <- species_positions
      }

      return(habitat_positions)
    }
  }

  # ── Shared helpers ───────────────────────────────────────────────────────────

  summarize_samples <- function(samples,
                                polygons,
                                habitat_raster,
                                habitat_names,
                                habitat_values,
                                species_groups,
                                species_group_names) {
    other_code <- max(habitat_values) + 1

    if (!any(tolower(habitat_names) == "other")) {
      habitat_names  <- c(habitat_names, "Other")
      habitat_values <- c(habitat_values, other_code)
    }

    results_df <- data.frame(matrix(ncol = length(habitat_names) +
                                      length(species_group_names) + 3,
                                    nrow = 0))
    colnames(results_df) <- c(habitat_names, "Area_Total", species_group_names, "Sp_Total", "Polygon_Area")

    for (i in seq_along(samples)) {
      sample       <- samples[[i]]
      polygon      <- polygons[[i]]
      polygon_area <- sf::st_area(polygon)

      habitat_cropped                      <- terra::crop(habitat_raster, terra::vect(polygon))
      habitat_cropped[is.na(habitat_cropped)] <- other_code
      habitat_masked                       <- terra::mask(habitat_cropped, terra::vect(polygon))
      habitat_df                           <- terra::freq(habitat_masked, bylayer = FALSE)

      all_values_df  <- data.frame(value = habitat_values, count = 0)
      habitat_df     <- merge(all_values_df, habitat_df, by = "value", all.x = TRUE)
      habitat_df$count <- ifelse(is.na(habitat_df$count.y), 0, habitat_df$count.y)
      habitat_df     <- habitat_df[, c("value", "count")]
      habitat_df$area <- habitat_df$count *
        terra::res(habitat_raster)[1] * terra::res(habitat_raster)[2]

      results_df[i, seq_along(habitat_names)] <- habitat_df$area
      results_df[i, length(habitat_names) + 1] <- sum(results_df[i, seq_along(habitat_names)])
      results_df[i, "Polygon_Area"] <- as.numeric(polygon_area)

      total_species <- 0
      for (k in seq_along(species_groups)) {
        # + 1 skips locationID column (col 1) in the sf object
        group_species   <- sf::st_drop_geometry(sample)[, species_groups[[k]] + 1]
        species_present <- colSums(group_species > 0)
        num_species     <- sum(species_present > 0)
        results_df[i, length(habitat_names) + 1 + k] <- num_species
        total_species   <- total_species + num_species
      }
      results_df[i, length(habitat_names) + length(species_groups) + 2] <- total_species
    }
    return(results_df)
  }

  analyze_sar <- function(results_table, method_used) {
    sar_results <- list(valid = FALSE, message = "")

    if (nrow(results_table) == 0) {
      sar_results$message <- "No data for SAR analysis"
      return(sar_results)
    }

    valid_rows <- results_table$Sp_Total > 0 & !is.na(results_table$Sp_Total) &
      is.finite(results_table$Sp_Total)

    if (sum(valid_rows) < 2) {
      sar_results$message <- paste("Insufficient non-zero samples for SAR analysis (need at least 2, have",
                                   sum(valid_rows), ")")
      return(sar_results)
    }

    valid_table <- results_table[valid_rows, ]
    area_valid  <- valid_table$Area_Total > 0 & !is.na(valid_table$Area_Total) &
      is.finite(valid_table$Area_Total)

    if (sum(area_valid) < 2) {
      sar_results$message <- "Insufficient valid area values for SAR analysis"
      return(sar_results)
    }

    valid_table <- valid_table[area_valid, ]
    log_area    <- log(valid_table$Area_Total)
    log_sp      <- log(valid_table$Sp_Total)

    if (any(!is.finite(log_area)) || any(!is.finite(log_sp))) {
      sar_results$message <- "Infinite values in log-transformed data"
      return(sar_results)
    }

    sar_results$log_area       <- log_area
    sar_results$log_sp         <- log_sp
    sar_results$n_samples_used <- nrow(valid_table)
    sar_results$valid          <- TRUE

    tryCatch({
      sar_results$lm_model   <- lm(log_sp ~ log_area)
      sar_results$lm_summary <- summary(sar_results$lm_model)
      sar_results$message    <- "SAR analysis completed successfully"
    }, error = function(e) {
      sar_results$message <- paste("Error fitting linear model:", e$message)
      sar_results$valid   <- FALSE
    })

    sar_results$sar_data <- valid_table
    return(sar_results)
  }

  analyze_countryside_sar <- function(csar_data, habitat_names, species_group_names) {
    csar_results <- list(valid = FALSE, message = "")

    if (nrow(csar_data) < 2) {
      csar_results$message <- "Insufficient data for countryside SAR (need at least 2 samples)"
      return(csar_results)
    }

    tryCatch({
      if (!requireNamespace("sars", quietly = TRUE)) {
        csar_results$message <- "Package 'sars' is required"
        return(csar_results)
      }
      csar_results$model <- sars::sar_countryside(
        data      = csar_data,
        modType   = "power",
        gridStart = "partial",
        habNam    = habitat_names,
        spNam     = species_group_names
      )
      csar_results$valid   <- TRUE
      csar_results$message <- "Countryside SAR completed successfully"
    }, error = function(e) {
      csar_results$message <- paste("Error:", e$message)
      csar_results$valid   <- FALSE
    })

    return(csar_results)
  }

  #---------------------------- 4. Main processing -----------------------------
  points_sf <- sf::st_as_sf(data, coords = c("long", "lat"), crs = crs)

  if (method == "circles") {

    # ----- Step 1: Build species groups -----
    sp_groups <- list()
    grp_names <- character()

    for (col in group_cols) {
      species_in_group <- classification[classification[[col]] == 1, ][[species_name_col]]
      positions        <- match(species_in_group, colnames(data)[4:ncol(data)])
      positions        <- positions[!is.na(positions)]
      if (length(positions) > 0) {
        sp_groups <- c(sp_groups, list(positions))
        grp_names <- c(grp_names, paste0("Sp_", col))
      }
    }
    species_group_names <- grp_names

    # ----- Step 2: Define boundary -----
    convex_hull <- if (!is.null(custom_hull)) {
      custom_hull
    } else {
      sf::st_convex_hull(sf::st_union(points_sf))
    }

    # ----- Step 3: Run iterations -----
    runs <- list()

    for (run in 1:n_runs) {
      samples        <- filter_points_in_expanding_circles(points_sf, radius,
                                                           convex_hull, break_threshold)
      polygons       <- lapply(samples, `[[`, "circle")
      samples_points <- lapply(samples, `[[`, "points")

      results_table <- summarize_samples(
        samples             = samples_points,
        polygons            = polygons,
        habitat_raster      = habitat,
        habitat_names       = habitat_names,
        habitat_values      = habitat_codes,
        species_groups      = sp_groups,
        species_group_names = grp_names
      )

      sar_analysis <- analyze_sar(results_table, method)

      runs[[run]] <- list(
        run           = run,
        results_table = results_table,
        sar_analysis  = sar_analysis,
        samples       = samples,
        polygons      = polygons
      )
    }

    # ----- Step 4: Average SAR results across runs -----
    if (n_runs > 1) {

      valid_runs <- which(sapply(runs, function(x) x$sar_analysis$valid))

      if (length(valid_runs) > 0) {
        slopes     <- sapply(valid_runs, function(i) coef(runs[[i]]$sar_analysis$lm_model)[2])
        intercepts <- sapply(valid_runs, function(i) coef(runs[[i]]$sar_analysis$lm_model)[1])
        r_squared  <- sapply(valid_runs, function(i) runs[[i]]$sar_analysis$lm_summary$r.squared)

        avg_sar_results <- list(
          valid            = TRUE,
          message          = paste("Average SAR results from", length(valid_runs), "valid runs"),
          n_valid_runs     = length(valid_runs),
          avg_slope        = mean(slopes),
          avg_intercept    = mean(intercepts),
          avg_r_squared    = mean(r_squared),
          sd_slope         = sd(slopes),
          sd_intercept     = sd(intercepts),
          sd_r_squared     = sd(r_squared),
          slopes           = slopes,
          intercepts       = intercepts,
          r_squared_values = r_squared
        )
      } else {
        avg_sar_results <- list(valid = FALSE, message = "No valid SAR runs to average")
      }

      # ----- Step 5: Return results (multiple runs) -----
      res <- list(
        method          = method,
        n_runs          = n_runs,
        runs            = runs,
        avg_sar_results = avg_sar_results,
        points_sf       = points_sf,
        convex_hull     = convex_hull,
        hull_source     = ifelse(is.null(custom_hull), "derived", "custom")
      )

    } else {

      # ----- Step 5: Return results (single run) -----
      res <- list(
        method      = method,
        n_runs      = n_runs,
        runs        = runs,
        sar_results = runs[[1]]$sar_analysis,
        points_sf   = points_sf,
        convex_hull = convex_hull,
        hull_source = ifelse(is.null(custom_hull), "derived", "custom")
      )
    }

  } else {

    # ----- Clusters path -----
    squares_sf <- create_squares(points_sf, square_size)
    samples    <- filter_points_in_clusters(points_sf, squares_sf, cluster_sizes)

    classif_subset <- if (!is.null(groups)) {
      classification[, c(species_name_col, group_cols)]
    } else {
      classification
    }

    # Pass full points_sf (locationID stays as col 1 so + 1 offset in summarize_samples is correct)
    species_groups      <- extract_species_positions(species_habitat_matrix = classif_subset,
                                                     species_site_matrix    = points_sf)
    species_group_names <- names(species_groups)

    samples_flat  <- unlist(lapply(samples, `[[`, "points"), recursive = FALSE)
    polygons_flat <- unlist(lapply(samples, `[[`, "chulls"), recursive = FALSE)

    results_table <- summarize_samples(
      samples             = samples_flat,
      polygons            = polygons_flat,
      habitat_raster      = habitat,
      habitat_names       = habitat_names,
      habitat_values      = habitat_codes,
      species_groups      = species_groups,
      species_group_names = species_group_names
    )

    sar_analysis  <- analyze_sar(results_table, method)
    csar_data     <- results_table[, c(habitat_names, species_group_names)]
    csar_analysis <- analyze_countryside_sar(csar_data, habitat_names, species_group_names)

    res <- list(
      method          = method,
      results_table   = results_table,
      sar_analysis    = sar_analysis,
      csar_analysis   = csar_analysis,
      samples         = samples,
      squares_sf      = squares_sf,
      clusters_chulls = lapply(samples, `[[`, "chulls"),
      points_sf       = points_sf
    )
  }

  return(res)
}
