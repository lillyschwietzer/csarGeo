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
#' @param area_unit 
#' @param plot_runs Number of runs done using method "Circles". Not available for "Clusters".
#' @param show_average description
#' @param average_color description
#' @param average_lwd description
#' @param show_formula description
#' @param show_slope description
#' @param show_intercept description
#' @param show_p_value description
#' @param show_r_squared description
#' @param show_adj_r_squared description
#' @param show_run_legend description
#'
#' @return No return value, function creates plots.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_countryside_sar(res, plot_type = "both")
#' }
visual_sar <- function(result,
                       method = NULL,

                       
                       plot_type = "both",
                       habitat_raster = NULL,
                       cluster_level = "size_256",
                       plot_all_levels = FALSE,
                       ncol = NULL,
                       # Spatial plot options
                       point_size = 16,
                       point_col = "black",
                       hull_col = "red",
                       area_unit = "m",
                       # Multiple runs options
                       plot_runs = FALSE,
                       show_average = TRUE,
                       average_color = "red",
                       average_lwd = 3,
                       # SAR plot options
                       show_formula = TRUE,
                       show_slope = TRUE,
                       show_intercept = FALSE,
                       show_p_value = FALSE,
                       show_r_squared = TRUE,
                       show_adj_r_squared = FALSE,
                       show_run_legend = TRUE) {

  # ------------ Data Access------------  # or better in main?
  if (...) { # plot map
    access following data from res: xyz
    then assign parameters
  }


  #--------------------- Cluster Data

  
  res_clusters_sar <- res_clusters[["sar_analysis"]][["lm_model"]][["model"]]
  
  #------------------------ 1) Helper Functions ---------------------------
  
  if (method == "circles") {
    # Circles Map
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
  print(pt_in_circles[[i]]$points) # access points for circle i
  plot(st_geometry(pt_in_circles[[i]]$circle),  border = "blue", add=TRUE) # access circle geometry for circle
}
    
  # Circles SAR
  plot_sar_circles <- function(Area_Total,
                               Sp_Total){
    plot(log(res$Area_Total),log(res$Sp_Total)) # logscale plot of area vs. species richness
    abline(sar=lm(log(res$Sp_Total)~log(res$Area_Total)))
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
 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        if (plot_type == "spatial" || plot_type == "both")
  {

    if (result$method == "circles")
    {
      # Circles method default single plot
      if (plot_all_levels)
      {
        warning("plot_all_levels = TRUE is not applicable for circles method. Showing single plot.")
      }

      # Plot points
      plot(sf::st_geometry(result$points_sf),
           main = "Expanding Circles Sampling",
           pch = point_size, cex = 0.7, col = point_col)

      # Add convex hull
      plot(result$convex_hull, border = hull_col, lwd = 2, add = TRUE)

      # Add circles
      if (!is.null(result$runs))
      {
        # For structure with runs (even single run)
        for (run_idx in seq_along(result$runs))
          {
          if (!is.null(result$runs[[run_idx]]$samples))
            {
            for (sample_idx in seq_along(result$runs[[run_idx]]$samples))
              {
              if (!is.null(result$runs[[run_idx]]$samples[[sample_idx]]$circle))
                {
                plot(sf::st_geometry(result$runs[[run_idx]]$samples[[sample_idx]]$circle),
                     border = "blue", add = TRUE)
                }
              }
             }
          }
       } else
        {
        # Fallback for older structure
        for (i in seq_along(result$samples))
         {
          if (!is.null(result$samples[[i]]$circle))
          {
            plot(sf::st_geometry(result$samples[[i]]$circle),
                 border = "blue", add = TRUE)
          }
         }
        }

    } else # clusters
      {

      if (plot_all_levels) { # option 1: plot all levels in one grid layout
        # Access samples through runs[[1]] for clusters method
        all_levels <- names(result$runs[[1]]$samples)
        n_levels <- length(all_levels)

        # Check if there are levels to plot
        if (n_levels == 0)
        {
          warning("No cluster levels found to plot")
          return(invisible(result))
        }

        # Determine grid layout
        if (is.null(ncol))
        {
          ncol <- ceiling(sqrt(n_levels))
        }
        nrow <- ceiling(n_levels / ncol)

        # Ensure nrow and ncol are valid positive integers
        nrow <- max(1, nrow)
        ncol <- max(1, ncol)

        # Set up plotting area
        par(mfrow = c(nrow, ncol), mar = c(2, 2, 3, 1))

        # Plot each level
        for (level_name in all_levels)
        {
          # Access through runs[[1]]
          hulls_at_level <- result$runs[[1]]$samples[[level_name]]$chulls
          n_clusters <- length(hulls_at_level)

          # Get points at that level
          points_list <- result$runs[[1]]$samples[[level_name]]$points

          if (is.null(points_list))
          {
            warning("No points found for level ", level_name)
            next
          }

          # Try to combine points
          points_at_level <- tryCatch({
            if (inherits(points_list, "sf"))
            {
              points_list
            } else if (is.list(points_list))
              {
              do.call(rbind, points_list)
              } else
                {
              warning("Unexpected structure for points at level ", level_name)
              next
                }
          }, error = function(e)
            {
            warning("Error combining points for level ", level_name, ": ", e$message)
            next
            })

          # Create empty plot first
          plot(sf::st_geometry(points_at_level),
               main = paste(level_name, "- Clusters:", n_clusters),
               pch = point_size, cex = 0.5, col = point_col,
               cex.main = 0.9, type = "n")

          # Add filled cluster polygons
          if (n_clusters > 0)
          {
            colors <- rainbow(n_clusters)
            for (i in 1:n_clusters)
            {
              if (!is.null(hulls_at_level[[i]]))
              {
                plot(hulls_at_level[[i]], col = colors[i], border = "black", add = TRUE)
              }
            }
          }

          # Add points
          plot(sf::st_geometry(points_at_level),
               pch = point_size, cex = 0.5, col = point_col, add = TRUE)

          # Add point count info
          mtext(paste(nrow(points_at_level), "points"), side = 3, line = 0, cex = 0.7)
        }

        # Reset plotting parameters
        par(mfrow = c(1, 1))

      } else { # option 2: plot only one level
        # Access through runs[[1]]
        level_idx <- which(names(result$runs[[1]]$samples) == cluster_level)
        if (length(level_idx) == 0)
        {
          warning("Cluster level '", cluster_level, "' not found. Using first level.")
          level_idx <- 1
          cluster_level <- names(result$runs[[1]]$samples)[1]
        }

        hulls_at_level <- result$runs[[1]]$samples[[level_idx]]$chulls
        n_clusters <- length(hulls_at_level)

        # Get points at that level
        points_list <- result$runs[[1]]$samples[[level_idx]]$points

        if (is.null(points_list))
        {
          stop("No points found for level ", cluster_level)
        }

        # Try to combine points
        points_at_level <- tryCatch(
        {
          if (inherits(points_list, "sf"))
          {
            points_list
          } else if (is.list(points_list))
            {
            do.call(rbind, points_list)
            } else
              {
            stop("Unexpected structure for points at level ", cluster_level)
              }
        }, error = function(e)
          {
          stop("Error combining points for level ", cluster_level, ": ", e$message)
          })

        # Create empty plot first
        plot(sf::st_geometry(points_at_level),
             main = paste(cluster_level, "-", n_clusters, "clusters"),
             pch = point_size, cex = 0.7, col = point_col,
             type = "n")

        # Add filled cluster polygons
        if (n_clusters > 0)
        {
          colors <- rainbow(n_clusters)
          for (i in 1:n_clusters)
          {
            if (!is.null(hulls_at_level[[i]]))
            {
              plot(hulls_at_level[[i]], col = colors[i], border = "black", add = TRUE)
            }
          }
        }

        # Add points on top
        plot(sf::st_geometry(points_at_level),
             pch = point_size, cex = 0.7, col = point_col, add = TRUE)
      }
    }
  }

  #-------------------------- 2.) Plot SAR regression line(s) --------------------

  if (plot_type == "sar" || plot_type == "both")
  {
    # Create a new plot or add to existing
    if (plot_type == "both")
    {
      if (plot_all_levels && result$method == "clusters")
      {
        # If we already did multi-panel spatial plots, create new device for SAR
        dev.new()
      } else
        {
        dev.new()
        }
    }

    # If result has more than 1 run -> plot multiple runs with grey lines + red average
    has_multiple_runs <- !is.null(result$n_runs) && result$n_runs > 1 && plot_runs

    if (has_multiple_runs)
    {
      all_x <- c()
      all_y <- c()

      # Determine plot limits from data
      for (run in 1:result$n_runs)
      {
        if (result$runs[[run]]$sar_analysis$valid)
        {
          all_x <- c(all_x, result$runs[[run]]$sar_analysis$log_area)
          all_y <- c(all_y, result$runs[[run]]$sar_analysis$log_sp)
        }
      }

      # Create empty plot with limits
      x_label <- paste0("log(Area) (", area_unit, "\u00B2)")
      y_label <- "log(Species Richness)"
      main_title <- paste("Species-Area Relationship (", result$n_runs, " runs)")

      plot(1, type = "n",
           xlim = range(all_x), ylim = range(all_y),
           xlab = x_label, ylab = y_label, main = main_title)

      # Storage for regression statistics
      slopes <- c()
      intercepts <- c()
      r_squareds <- c()

      # First pass: plot all regression lines as grey dotted lines
      for (run in 1:result$n_runs)
      {
        if (result$runs[[run]]$sar_analysis$valid)
        {
          model <- result$runs[[run]]$sar_analysis$lm_model

          # Regression line as grey dotted
          abline(model, col = "grey40", lwd = 1, lty = 3)

          # Store statistics for average
          slopes <- c(slopes, coef(model)[2])
          intercepts <- c(intercepts, coef(model)[1])
          r_squareds <- c(r_squareds, summary(model)$r.squared)
        }
      }

      # Second pass: plot all points in black
      for (run in 1:result$n_runs)
      {
        if (result$runs[[run]]$sar_analysis$valid)
        {
          x <- result$runs[[run]]$sar_analysis$log_area
          y <- result$runs[[run]]$sar_analysis$log_sp

          # Add data points
          points(x, y, pch = 16, col = "black", cex = 0.7)
        }
      }

      # Calculate and plot average as red line if requested
      if (show_average && length(slopes) > 0)
      {
        avg_slope <- mean(slopes)
        avg_intercept <- mean(intercepts)

        # Create average line across x range
        x_range <- range(all_x)
        avg_line_x <- seq(x_range[1], x_range[2], length.out = 100)
        avg_line_y <- avg_intercept + avg_slope * avg_line_x

        # Plot average as red line
        lines(avg_line_x, avg_line_y, col = average_color, lwd = average_lwd, lty = 1)

        # Add average statistics
        avg_r2 <- mean(r_squareds)
        avg_stats <- paste("Avg slope =", round(avg_slope, 3),
                           "| Avg R² =", round(avg_r2, 3))

        # Add average info to plot
        mtext(avg_stats, side = 3, line = 0.5, cex = 0.8)
      }

      # Optional: Add legend
      if (show_run_legend)
      {
        legend("bottomright",
               legend = c("Individual runs",
                          "Data points",
                          "Average of individual runs"),
               col = c("grey40", "black", average_color),
               lty = c(3, NA, 1),
               pch = c(NA, 16, NA),
               lwd = c(1, NA, average_lwd),
               bg = "white", cex = 0.8)
      }

      grid()

    } else
      {
      # Single run plot (original behavior)
      if (is.null(result$sar_analysis$valid) || !result$sar_analysis$valid)
      {
        # If no valid SAR, try to get from first run if multiple runs exist
        if (!is.null(result$n_runs) && result$n_runs > 0 &&
            !is.null(result$runs[[1]]$sar_analysis$valid) &&
            result$runs[[1]]$sar_analysis$valid)
        {
          # Use first run's SAR analysis
          x_vals <- result$runs[[1]]$sar_analysis$log_area
          y_vals <- result$runs[[1]]$sar_analysis$log_sp
          lm_model <- result$runs[[1]]$sar_analysis$lm_model
          lm_summary <- result$runs[[1]]$sar_analysis$lm_summary
          valid_sar <- TRUE
        } else
          {
          x_vals <- result$sar_analysis$log_area
          y_vals <- result$sar_analysis$log_sp
          lm_model <- result$sar_analysis$lm_model
          lm_summary <- result$sar_analysis$lm_summary
          valid_sar <- result$sar_analysis$valid
          }
      } else
        {
        x_vals <- result$sar_analysis$log_area
        y_vals <- result$sar_analysis$log_sp
        lm_model <- result$sar_analysis$lm_model
        lm_summary <- result$sar_analysis$lm_summary
        valid_sar <- result$sar_analysis$valid
        }

      if (valid_sar)
        {
        # create labels
        x_label <- paste0("log(Area) (", area_unit, "\u00B2)")
        y_label <- "log(Species Richness)"
        main_title <- paste("Species-Area Relationship (", result$method, ")")

        # Plot with arguments
        graphics::plot(x = x_vals, y = y_vals,
                       xlab = x_label,
                       ylab = y_label,
                       main = main_title,
                       pch = 16, col = "black")

        # Add regression line
        graphics::abline(lm_model, col = "red", lwd = 2)

        # Extract model statistics
        coefs <- coef(lm_model)
        intercept <- round(coefs[1], 3)
        slope <- round(coefs[2], 3)
        r_squared <- round(lm_summary$r.squared, 3)
        adj_r_squared <- round(lm_summary$adj.r.squared, 3)
        p_value <- round(lm_summary$coefficients[2, 4], 4)

        # Build legend text based on user choices
        legend_text <- c()

        if (show_formula)
        {
          formula_text <- paste("log(S) =", intercept, "+", slope, "× log(A)")
          legend_text <- c(legend_text, formula_text)
        }

        if (show_slope)
        {
          legend_text <- c(legend_text, paste("slope =", slope))
        }

        if (show_intercept)
        {
          legend_text <- c(legend_text, paste("intercept =", intercept))
        }

        if (show_r_squared)
        {
          legend_text <- c(legend_text, paste("R\u00B2 =", r_squared))
        }

        if (show_adj_r_squared)
        {
          legend_text <- c(legend_text, paste("adj R\u00B2 =", adj_r_squared))
        }

        if (show_p_value)
        {
          # Format p-value nicely
          if (p_value < 0.0001)
          {
            p_text <- "p < 0.0001"
          } else
            {
            p_text <- paste("p =", p_value)
            }
          legend_text <- c(legend_text, p_text)
        }

        # Add legend if there's anything to show
        if (length(legend_text) > 0)
        {
          graphics::legend("topleft",
                           legend = legend_text,
                           bty = "n", cex = 0.9)
        }

        graphics::grid()
      } else
        {
        # Create empty plot with message
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = "Species-Area Relationship")
        text(1, 1, "SAR analysis not available\nor invalid", cex = 1.2)
        }
    }
  }

  # If both plot types and all levels -> create a summary text
  if (plot_type == "both" && plot_all_levels && result$method == "clusters")
  {
    cat("\n=== Summary of Cluster Levels ===\n")
    for (level_name in names(result$runs[[1]]$samples))
    {
      n_clusters <- length(result$runs[[1]]$samples[[level_name]]$chulls)
      n_points <- sum(sapply(result$runs[[1]]$samples[[level_name]]$points, nrow))
      cat(sprintf("%s: %d clusters, %d points\n", level_name, n_clusters, n_points))
    }

    # Show SAR parameters by level
    cat("\n=== SAR Parameters by Level ===\n")
    for (level in unique(result$runs[[1]]$results_table$Cluster_Level))
    {
      level_data <- result$runs[[1]]$results_table[result$runs[[1]]$results_table$Cluster_Level == level, ]
      if (nrow(level_data) > 1)
      {
        lm_level <- lm(log(Sp_Total) ~ log(Area_Total), data = level_data)
        cat(sprintf("%s: slope = %.3f, R\u00B2 = %.3f, n = %d\n",
                    level, coef(lm_level)[2], summary(lm_level)$r.squared, nrow(level_data)))
      }
    }
  }

  # Return invisible for potential further use
  invisible(result)
}
