#' Visualize countryside_sar results
#'
#' @description
#' A function to visualize the results of the function countryside_sar. It can create a "spatial" plot for the sampling process of either method and a "sar" plot for the results of the SAR-analysis and linear model results.
#'
#'
#' @param result Output object from countryside_sar function.
#' @param plot_type Type of plot: "spatial" for sampling design, "sar" for species-area relationship, or "both" in a grid.
#' @param habitat_raster Optional parameter to add a land-use raster background.
#' @param cluster_level For clusters method, defines which level to plot (e.g., "size_256")
#' @param plot_all_levels For "clusters" method, plot all levels at once in an arranged grid. If plot_all_levels = TRUE it is not possible to define a value for cluster_level.
#' @param ncol description
#' @param point_size Optional parameter to change the point size of the plot.
#' @param point_col Optional parameter to change the colour of the points.
#' @param hull_col Optional parameter to change the hull colour of the plot.
#'
#' @return No return value, function creates plots.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_countryside_sar(res, plot_type = "both")
#' }
visual_sar <- function(result) {


  method <- result(method)
  #--------------------- Data Manipulation: Circles -------------------------

prepare_circles_for_plot <- function(res_circles, 
                                     run_number = 1,
                                     circle_color = "blue",
                                     highlight_radii = NULL) {
  
  # Extract samples for specified run
  samples <- res_circles[["runs"]][[run_number]][["samples"]]
  
  # Extract circle geometries
  circle_polygons <- lapply(samples, function(x) x$circle)
  
  # Optional: Name them by radius for reference
  names(circle_polygons) <- names(samples)
  
  # Return a clean list with everything the plot function needs
  return(list(
    points = res_circles[["points_sf"]],
    circles = circle_polygons,
    hull = res_circles[["convex_hull"]]
  ))
}
  
  #------------------------ 1) Helper Functions ---------------------------
  
  if (method == "circles") {
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
    
  } else { # clusters method functions
    
    # Clusters Map
    plot_spatial_clusters <- function (clusters_chulls,
                                      points){
       for (size in seq_along(clusters_chulls)) # loops through each level of hulls
{
  # Set up colors
  colors <- rainbow(length(clusters_chulls[[size]]))
  plot(st_geometry(points)) # new plot for each level of clustering, Add remaining polygons with different colors to existing plot
  for (i in 1:length(clusters_chulls[[size]])) {
    plot(clusters_chulls[[size]][[i]], col = colors[i], border = "black", add = TRUE)
  }
 }
}


    # Clusters SAR
    plot_clusters_sar <- function(Area_Total,
                                  Sp_Total){
    plot(res_clusters_sar$log_area,res_clusters_sar$log_sp) # logscale plot of area vs. species richness
    abline(sar=lm(res_clusters_sar$log_sp~(res_clusters_sar$log_area))
}


    # Clusters cSAR
    plot_clusters_csar

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
if (method == "circles"){
  
  n_runs <- res_circles[["n_runs"]]  # number of runs
  
  if (n_runs > 1) { # if circles ran more than once 
    
    # Calculate grid dimensions
    n_cols <- ceiling(sqrt(n_runs))
    n_rows <- ceiling(n_runs / n_cols)
    
    # Set up plotting window
    par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))
    
    # Loop through each run
    for (run in 1:n_runs) {
      # Extract circle data 
      circles_for_plot <- res_circles[["runs"]][[run]][["samples"]]
      circle_geometries <- lapply(circles_for_plot, function(x) x$circle)
      
      # Plot with run-specific title
      plot_spatial_circles(
        points_sf = res_circles[["points_sf"]],
        circles_list = circle_geometries,
        convex_hull = res_circles[["convex_hull"]],
        main_title = paste("Run", run) 
      )
    }
    
    # Reset plotting parameters
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
    
  } else { # Single run 
    circles_for_plot <- res_circles[["runs"]][[1]][["samples"]]
    circle_geometries <- lapply(circles_for_plot, function(x) x$circle)
    
    plot_spatial_circles(
      points_sf = res_circles[["points_sf"]],
      circles_list = circle_geometries,
      convex_hull = res_circles[["convex_hull"]],
      main_title = "Sampling Circles"
    )
  }

                                
} else {

# map plot


  # sar plot


  # csar plot


  # Heatmap affinities
  # extract affinity data and create matrix
  affinity_matrix <- csar_res$model$affinity
  affinity_df <- do.call(rbind, affinity_matrix)

  plot_heat_clusters <- function(affinity_df)
  
}
