#' Visualize countryside_sar results
#'
#' @description
#' A function to visualize the results of the function countryside_sar. It can create a "spatial" plot for the sampling process of either method and a "sar" plot for the results of the SAR-analysis and linear model results.
#'
#'
#' @param result Output object from countryside_sar function.
#' @param plot_type Type of plot: "spatial" for sampling design, "sar" for species-area relationship, or "both" in a grid.
#'
#' @return A plot or a grid of plots of the given parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_countryside_sar(res, plot_type = "both")
#' }
visual_sar <- function(result,
                      plot_type = NULL, # ("map","sar","csar")
                      # Circles
                      plot_all_runs = TRUE
                      plot_run_n = NULL,
                      # Clusters
                      plot_all_levels = TRUE,
                      plot_level = NULL) 
{

  main_title

#----------------- Input infos ------------------
 method <- result[["method"]]

#------------------------ 1) Helper Functions ---------------------------

  if (method == circles) {
  
prepare_circles_for_plot <- function(res_circles,  # data prep
                                     run_number = 1) {
  
  # Extract samples + geometries for specified run
  samples <- res_circles[["runs"]][[run_number]][["samples"]]
  circle_polygons <- lapply(samples, \(x) x$circle)
  
  # Return list for plot function
  return(list(
    points = res_circles[["points_sf"]],
    circles = circle_polygons,
    hull = res_circles[["convex_hull"]]
  ))
}
  
    #-------------- Circles Map ---------------

    plot_spatial_circles <- function(points_sf, 
                                     circles, 
                                     convex_hull) {
  plot(st_geometry(points_sf), main = "Sampling Circles")
  plot(convex_hull, border = "red", lwd = 2, add = TRUE)
  for (circle in circles) { 
    plot(circle, border = "blue", add = TRUE)
  }
}

    plot(st_geometry(points)) # st_geometry extracts only geometry components from points -> just coordinates
plot(convex_hull, border = "red", lwd = 2, add=TRUE)

for (i in seq_along(pt_in_circles)) # loop through the number of circles created
{
  print(pt_in_circles[[i]]$points) 
  plot(st_geometry(pt_in_circles[[i]]$circle),  border = "blue", add=TRUE) 
}
    
  #---------------- Circles SAR ---------------
plot_sar_circles <- function(sar_results) {
  # Contains: log_area, log_sp, lm_model, etc.
  
  # Extract log-transformed data
  log_area <- sar_results[["log_area"]]
  log_sp <- sar_results[["log_sp"]]
  
  # Create plot
  plot(log_area, log_sp, 
       xlab = "log(Area)", 
       ylab = "log(Species Richness)",
       main = "Species-Area Relationship (SAR)",
       pch = 16)
  
  # Add regression line
  abline(sar_results[["lm_model"]], col = "red", lwd = 2)
  
  # Add R-squared value
  r2 <- sar_results[["lm_summary"]][["r.squared"]]
  legend("bottomright", 
         legend = c(paste("R² =", round(r2, 3))),
         bty = "n")
}
    
  } else { #----------- clusters method functions -----------------
    
    #---------------- Clusters Map -----------------
 plot_spatial_clusters <- function(clusters_chulls, points, main_title = NULL) {
  
  for (level in seq_along(clusters_chulls)) {
    # Use provided title or default
    if (is.null(main_title)) {
      main_title <- paste("Clustering Level", level)
    }
    
    plot(sf::st_geometry(points), 
         main = main_title,
         cex = 0.5, pch = 16)
    
    hulls <- clusters_chulls[[level]]
    if (length(hulls) > 0) {
      colors <- rainbow(length(hulls))
      for (i in seq_along(hulls)) {
        plot(hulls[[i]], 
             border = "black", 
             col = adjustcolor(colors[i], alpha.f = 0.5),
             add = TRUE)
      }
    }
  }
}

    #---------------- Clusters SAR ---------------
      plot_sar_clusters <- function(sar_results) {
      # Contains: log_area, log_sp, lm_model, etc.
      
      # Extract log-transformed data
      log_area <- sar_results[["log_area"]]
      log_sp <- sar_results[["log_sp"]]
      
      # Create plot
      plot(log_area, log_sp, 
           xlab = "log(Area)", 
           ylab = "log(Species Richness)",
           main = "Species-Area Relationship (SAR) - Clusters",
           pch = 16)
      
      # Add regression line
      abline(sar_results[["lm_model"]], col = "red", lwd = 2)
      
      # Add R-squared value
      r2 <- sar_results[["lm_summary"]][["r.squared"]]
      legend("bottomright", 
             legend = c(paste("R² =", round(r2, 3))),
             bty = "n")
    }
           
    #----------- Affinity Heatmap # library(ggplot2) library(reshape2)
    plot_heat_clusters <- function(affinity_data) {

  # Convert list to matrix if needed
  affinity_mat <- as.matrix(affinity_data)
  
  # Convert to long format for ggplot
  affinity_long <- melt(affinity_mat, varnames = c("Species Group", "Habitat"))
  
  # Create heatmap with legend on left
  ggplot(affinity_long, aes(x = Habitat, y = Species_Group, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(heat.colors(10)), 
                         name = "Affinity",
                         limits = c(0, 1)) +
    labs(title = "Species Habitat Affinity") +
    theme_minimal() +
    theme(legend.position = "right",  
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = 0.5)  # Adjust tile shape
}
  }
           
  #---------------------- 2) Main Processing ---------------------

  #----------------------------- Circles -------------------------------
if (method == "circles"){
  
  n_runs <- res_circles[["n_runs"]]  # number of runs
  
  if (plot_all_runs == TRUE) { # Plot all runs that were made
    
    # Calculate grid dimensions
    n_cols <- ceiling(sqrt(n_runs))
    n_rows <- ceiling(n_runs / n_cols)
    
    # Set up plotting window
    par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))
    
    # Loop through each run
    for (run in 1:n_runs) {
      # Extract circle data 
      circles_for_plot <- res_circles[["runs"]][[run]][["samples"]]
      circle_geometries <- lapply(circles_for_plot, \(x) x$circle)
      
      # Plot one map per run
      plot_spatial_circles(
        points_sf = res_circles[["points_sf"]],
        circles_list = circle_geometries,
        convex_hull = res_circles[["convex_hull"]],
        main_title = paste("Run", run) 
      )
    }
    
    # Reset plotting parameters
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
    
  } if (plot_run_n == run_number) {
  
    circles_for_plot <- res_circles[["runs"]][[1]][["samples"]]
    circle_geometries <- lapply(circles_for_plot, \(x) x$circle)
    
    plot_spatial_circles(
      points_sf = res_circles[["points_sf"]],
      circles_list = circle_geometries,
      convex_hull = res_circles[["convex_hull"]],
      main_title = "Sampling Circles"
    )
  }

  #------------------------- sar -----------------------------

 n_runs <- res_circles[["n_runs"]]  # number of runs
  
  if (n_runs > 1) {

  
  plot_sar_circles(
} else { # only 1 run

  }
                                
} else { #--------------------- clusters ----------------------------

# map plot
points_combined <- do.call(rbind, res_clusters[["samples"]][["size_1"]][["points"]])

  # determine level count for plotting parameters
  n_levels <- length(res_clusters[["samples"]]) 
  
  if (n_levels > 1) { # Plot all levels
    n_cols <- ceiling(sqrt(n_levels))
    n_rows <- ceiling(n_levels / n_cols)

    # set plotting grid
    par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))

    level_names <- names(res_clusters[["clusters_chulls"]])

for (level in 1:n_levels) {
  plot_spatial_clusters(
    clusters_chulls = res_clusters[["clusters_chulls"]][level],
    points = points_combined
  )
  title(level_names[level]) 
}

    # reset plotting parameters
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
    
  } else { # Plot only 1 level
    plot_spatial_clusters(
      clusters_chulls = res_clusters[["clusters_chulls"]],
      points = points_combined
    )
  }

  #--------------------------- sar plot -----------------------------
 # Extract SAR results from clusters object
  sar_results <- res_clusters[["sar_analysis"]]
  
  # Plot SAR
  plot_sar_clusters(sar_results)
}

  #----------------------- Heatmap affinities -----------------------
  # extract affinity data and create matrix
  affinity_matrix <- csar_res$model$affinity
  affinity_df <- do.call(rbind, affinity_matrix)

  plot_heat_clusters <- function(affinity_df)
  
}
