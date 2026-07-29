#' Visualize countryside_sar() function results
#'
#' @description
#' A function to visualize the results of \code{countryside_sar()} from the csarGeo-Package.
#' It can generate three plot types: "map", a plot of the polygon hull with the sampling locations and the chosen sampling method, "sar" for the base SAR analysis and "csar", only available for a result of method "clusters" to create an affinity heat map of the species groups to different habitats.
#' The SAR results can be plotted as a single run result, a grid of multiple runs or as an average SAR.
#' The plots of the sampling maps can either be done per level or run or in a grid.
#'
#' @param result Output object from countryside_sar function() from the csarGeo-Package.
#' @param plot_type Defines the type of plot. Three possible options
#'   \itemize{
#'     \item \code{map} (Map of the polygon hull + sampling locations with the chosen sampling method. Adds sampling circles for method = 'circles' or plots the cluster )
#'     \item \code{sar} (log-log scaled linear regression)
#'     \item \code{csar} (heatmap of species group affinities to habitat types)
#'     }
#' \code{map} and \code{sar} are available for results of both \code{method = circles} and \code{method = clusters}. \code{plot_type = csar} is only available for results of \code{method = clusters}.
#' @param plot_all_runs Plot all runs of "circles" into one grid if n_runs > 1. Defaults to TRUE.
#' @param plot_run_n Plot only the result of the n-th run of "circles". Defaults to NULL.
#' @param plot_average If TRUE, plots the average SAR result if n_runs > 1. Defaults to FALSE.
#' @param plot_all_levels Plot all clustering levels of method "clusters" into one grid. Defaults to TRUE.
#' @param plot_level_n Plot only the n-th clustering level as a solo plot. Defaults to NULL.
#'
#' @return Creates either a single plot or a grid of multiple plots. The function returns the input result object invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- countryside_sar(
#'   data = species_data,
#'   method = "circles",
#'   radius = 2000 * 1:10,
#'   habitat = land_use,
#'   habitat_names = c("Forest", "Agriculture", "Shrubland"),
#'   classification = classes_clusters,
#'   groups = c("Forest_Sp", "Grassland_Sp", "generalists_Sp")
#' )
#'
#' # Map plot
#' visuals_sar(res, plot_type = "map")
#'
#' # SAR plot
#' visuals_sar(res, plot_type = "sar")
#'
#' # Average SAR (circles only)
#' visuals_sar(res, plot_type = "sar", plot_average = TRUE)
#'
#' res <- countryside_sar(
#'   data = data,
#'   method = "clusters",
#'   square_size = 2000,
#'   cluster_sizes = c(1, 4, 16, 64, 256),
#'   habitat = lu1995,
#'   habitat_names = c("Forest", "Agriculture", "Shrubland"),
#'   classification = classes_clusters,
#'   groups = c("Forest_Sp", "Grassland_Sp", "generalists_Sp")
#' )
#'
#' # CSAR heatmap (clusters only)
#' visuals_sar(test_clusters, plot_type = "csar")
#' }
visuals_sar <- function(result,
                        plot_type = NULL, # "map", "sar", "csar",
                        # Circles
                        plot_all_runs = TRUE,
                        plot_run_n = NULL,
                        plot_average = FALSE,
                        # Clusters
                        plot_all_levels = TRUE,
                        plot_level_n = NULL) {

  #----------------- Input validation ------------------
  if (!plot_type %in% c("map", "sar", "csar"))
  {
    stop("plot_type must be 'map', 'sar', or 'csar'")
  }

  method <- result[["method"]]

  # Check if SAR exists (for both circles and clusters)
  has_sar <- (!is.null(result[["sar_analysis"]]) &&
                !is.null(result[["sar_analysis"]]$valid) &&
                result[["sar_analysis"]]$valid) ||
    (!is.null(result[["avg_sar_results"]]) &&
       !is.null(result[["avg_sar_results"]]$valid) &&
       result[["avg_sar_results"]]$valid)

  # Check if cSAR exists
  has_csar <- !is.null(result[["csar_analysis"]]) &&
    !is.null(result[["csar_analysis"]]$valid) &&
    result[["csar_analysis"]]$valid

  if (plot_type == "sar" && !has_sar) {
    stop("SAR analysis not available in this result object.")
  }

  #------------------------ 1) Helper Functions ---------------------------

  # Circles data preparation
  prepare_circles_for_plot <- function(res_circles, run_number = 1)
  {
    samples <- res_circles[["runs"]][[run_number]][["samples"]]
    circle_polygons <- lapply(samples, `[[`, "circle")

    return(list(
      points = res_circles[["points_sf"]],
      circles = circle_polygons,
      hull = res_circles[["convex_hull"]]
    ))
  }

  if (method == "circles")
  {

    # Circles Map
    plot_spatial_circles <- function(points_sf, circles, convex_hull, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Sampling Circles"
      plot(sf::st_geometry(points_sf), main = main_title)
      plot(convex_hull, border = "red", lwd = 2, add = TRUE)
      # add sampling circles
      for (circle in circles)
      {
        plot(circle, border = "blue", add = TRUE)
      }
    }

    # Circles SAR - Single run
    plot_sar_circles <- function(sar_results, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Circles"

      log_area <- sar_results[["log_area"]]
      log_sp <- sar_results[["log_sp"]]

      intercept <- coef(sar_results[["lm_model"]])[1]
      slope <- coef(sar_results[["lm_model"]])[2]
      r_squared <- sar_results[["lm_summary"]]$r.squared

      plot(log_area, log_sp,
           xlab = "log(Area)",
           ylab = "log(Species Richness)",
           main = main_title,
           pch = 16)

      grid()

      abline(sar_results[["lm_model"]], col = "red", lwd = 2)

      legend_text <- paste0(
        "log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
        "Slope (z) = ", round(slope, 3), "\n",
        "R² = ", round(r_squared, 3)
      )
      legend("topleft", legend = legend_text, bty = "n", cex = 0.9)
    }

    # Circles SAR - Multiple runs in grid
    plot_sar_circles_grid <- function(sar_results, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "SAR"

      log_area <- sar_results[["log_area"]]
      log_sp <- sar_results[["log_sp"]]

      intercept <- coef(sar_results[["lm_model"]])[1]
      slope <- coef(sar_results[["lm_model"]])[2]
      r_squared <- sar_results[["lm_summary"]]$r.squared

      plot(log_area, log_sp,
           xlab = "",
           ylab = "",
           main = main_title,
           pch = 16,
           axes = FALSE)

      # Add custom axes
      box()
      axis(1, labels = FALSE)
      axis(2, labels = FALSE)

      grid()
      abline(sar_results[["lm_model"]], col = "red", lwd = 2)

      legend_text <- paste0(
        "log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
        "z = ", round(slope, 3), "\n",
        "R² = ", round(r_squared, 3)
      )
      legend("topleft",
             legend = legend_text,
             bty = "n",
             cex = 0.8,
             inset = c(-0.02, -0.02))
    }

    # Average Circles SAR
    plot_average_sar_circles <- function(avg_sar_results, runs, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Average of Runs"

      intercept <- avg_sar_results[["avg_intercept"]]
      slope <- avg_sar_results[["avg_slope"]]
      r_squared <- avg_sar_results[["avg_r_squared"]]
      n_valid <- avg_sar_results[["n_valid_runs"]]

      # Determine overall range for axes
      all_log_area <- unlist(lapply(runs, function(x) x$sar_analysis$log_area))
      all_log_sp <- unlist(lapply(runs, function(x) x$sar_analysis$log_sp))

      plot(all_log_area, all_log_sp,
           xlab = "log(Area)",
           ylab = "log(Species Richness)",
           main = main_title,
           pch = 16,
           col = adjustcolor("black", alpha.f = 0.4),
           cex = 0.7)

      grid()

      # Add average regression line
      x_vals <- range(all_log_area, na.rm = TRUE)
      y_vals <- intercept + slope * x_vals
      lines(x_vals, y_vals, col = "red", lwd = 3)

      legend_text <- paste0(
        "log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
        "Slope (z) = ", round(slope, 3), "\n",
        "R² = ", round(r_squared, 3), "\n",
        "Valid runs = ", n_valid
      )
      legend("topleft",
             legend = legend_text,
             bty = "n",
             cex = 0.8)
    }


  } else
  { # method == "clusters"

    # Clusters Map
    plot_spatial_clusters <- function(points, hulls, level_name = NULL, main_title = NULL)
    {
      if (is.null(main_title))
      {
        main_title <- if (!is.null(level_name)) level_name
        else "Clustering Pattern"
      }

      plot(sf::st_geometry(points),
           main = main_title,
           cex = 0.5, pch = 16)

      if (length(hulls) > 0)
      {
        colors <- rainbow(length(hulls))
        for (i in seq_along(hulls))
        {
          plot(hulls[[i]],
               border = "black",
               col = adjustcolor(colors[i], alpha.f = 0.5),
               add = TRUE)
        }
      }
    }

    # Clusters SAR
    plot_sar_clusters <- function(sar_results, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Clusters"

      log_area <- sar_results[["log_area"]]
      log_sp <- sar_results[["log_sp"]]

      intercept <- coef(sar_results[["lm_model"]])[1]
      slope <- coef(sar_results[["lm_model"]])[2]
      r_squared <- sar_results[["lm_summary"]]$r.squared

      plot(log_area, log_sp,
           xlab = "log(Area)",
           ylab = "log(Species Richness)",
           main = main_title,
           pch = 16)

      grid()

      abline(sar_results[["lm_model"]], col = "red", lwd = 2)

      legend_text <- paste0(
        "log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
        "Slope (z) = ", round(slope, 3), "\n",
        "R² = ", round(r_squared, 3)
      )
      legend("topleft",
             legend = legend_text,
             bty = "n",
             cex = 0.8)
    }

    # Affinity Heatmap
    plot_heatmap <- function(csar_results)
    {

      # Extract and convert affinity data
      affinity_list <- csar_results$model$affinity
      affinity_mat <- do.call(rbind, affinity_list)
      affinity_long <- reshape2::melt(affinity_mat,
                                      varnames = c("Species_Group", "Habitat"),
                                      value.name = "Affinity")

      # Add text color based on background brightness
      affinity_long$text_color <- ifelse(affinity_long$Affinity > 0.5, "white", "black")

      # Create heatmap
      ggplot2::ggplot(affinity_long,
                      ggplot2::aes(x = Habitat, y = Species_Group, fill = Affinity)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::scale_fill_viridis_c(
          name = NULL,
          limits = c(0, 1),
          option = "viridis",
          direction = -1,
          guide = "none"
        ) +
        ggplot2::geom_text(ggplot2::aes(label = round(Affinity, 3),
                                        color = text_color),
                           size = 4.5) +
        ggplot2::scale_color_identity() +
        ggplot2::labs(title = "Species Habitat Affinity",
                      x = "Habitat Type",
                      y = "Species Group") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          panel.grid = ggplot2::element_blank()
        )
    }
  }

  #---------------------- 2) Main Processing ---------------------

  if (method == "circles")
  {
    n_runs <- result[["n_runs"]]

    # ----------- Map Plot -----------
    if (plot_type == "map") {

      if (plot_all_runs == TRUE && n_runs > 1)
      {
        # Plot all runs in grid
        n_cols <- ceiling(sqrt(n_runs))
        n_rows <- ceiling(n_runs / n_cols)
        par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))

        for (run in 1:n_runs)
        {
          circles_data <- prepare_circles_for_plot(result, run)
          plot_spatial_circles(
            points_sf = circles_data$points,
            circles = circles_data$circles,
            convex_hull = circles_data$hull,
            main_title = paste("Run", run)
          )
        }
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

      } else
      {
        # Plot single run
        run_to_plot <- if (!is.null(plot_run_n)) plot_run_n else 1
        circles_data <- prepare_circles_for_plot(result, run_to_plot)
        plot_spatial_circles(
          points_sf = circles_data$points,
          circles = circles_data$circles,
          convex_hull = circles_data$hull,
          main_title = paste("Run", run_to_plot)
        )
      }
    }

    # ----------- SAR PLOT -----------
    if (plot_type == "sar")
    {

      if (plot_average == TRUE)
      {
        # Plot the averaged SAR result across all runs
        avg_sar_results <- result[["avg_sar_results"]]
        runs <- result[["runs"]]

        if (avg_sar_results$valid)
        {
          plot_average_sar_circles(avg_sar_results, runs)
        } else
        {
          cat("Average SAR analysis not valid\n")
          cat("Message:", avg_sar_results$message, "\n")
        }

      } else if (!is.null(plot_run_n))
      {
        # Plot a specific run
        sar_results <- result[["runs"]][[plot_run_n]][["sar_analysis"]]

        if (sar_results$valid)
        {
          plot_sar_circles(sar_results, main_title = paste("SAR - Run", plot_run_n))
        } else
        {
          cat("SAR analysis not valid for run", plot_run_n, "\n")
          cat("Message:", sar_results$message, "\n")
        }

      } else if (n_runs > 1)
      {
        # Grid for multiple runs
        n_cols <- ceiling(sqrt(n_runs))
        n_rows <- ceiling(n_runs / n_cols)

        # Set up grid with margins
        par(mfrow = c(n_rows, n_cols),
            mar = c(0.5, 0.5, 2, 0.5),
            oma = c(4, 4, 2, 2))

        for (run in 1:n_runs)
        {
          sar_results <- result[["runs"]][[run]][["sar_analysis"]]
          if (sar_results$valid)
          {
            plot_sar_circles_grid(sar_results, main_title = paste("Run", run))
          } else
          {
            plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(1, 1, paste("Run", run, "\nInvalid"), col = "red", cex = 1.2)
          }
        }

        # Add shared axis labels
        mtext("log(Area)", side = 1, outer = TRUE, line = 2.5, cex = 1.2)
        mtext("log(Species Richness)", side = 2, outer = TRUE, line = 2.5, cex = 1.2)

        # Reset plotting parameters
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

      } else
      {
        # Single run - with full axis labels
        run_to_plot <- if (!is.null(plot_run_n)) plot_run_n else 1
        sar_results <- result[["runs"]][[run_to_plot]][["sar_analysis"]]

        if (sar_results$valid)
        {
          plot_sar_circles(sar_results, main_title = paste("SAR - Run", run_to_plot))
        } else
        {
          cat("SAR analysis not valid for run", run_to_plot, "\n")
          cat("Message:", sar_results$message, "\n")
        }
      }
    }

  } else
  { # method == "clusters"

    # Extract data
    points_combined <- do.call(rbind, result[["samples"]][["size_1"]][["points"]])
    clusters_chulls <- result[["clusters_chulls"]]
    n_levels <- length(clusters_chulls)
    level_names <- names(clusters_chulls)

    # ----------- MAP PLOT -----------
    if (plot_type == "map")
    {

      if (plot_all_levels == TRUE && n_levels > 1)
      {
        # Plot all levels in grid
        n_cols <- ceiling(sqrt(n_levels))
        n_rows <- ceiling(n_levels / n_cols)
        par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))

        for (level in 1:n_levels)
        {
          display_name <- gsub("size_", "Size ", level_names[level])

          plot_spatial_clusters(
            points = points_combined,
            hulls = clusters_chulls[[level]],
            level_name = display_name
          )
        }
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

      } else
      {
        # Plot single level
        if (!is.null(plot_level_n))
        {
          target_name <- paste0("size_", plot_level_n)
          level_idx <- which(level_names == target_name)

          if (length(level_idx) == 0)
          {
            stop("Cluster size '", plot_level_n, "' not found. Available sizes: ",
                 paste(gsub("size_", "", level_names), collapse = ", "))
          }
          level_to_plot <- level_idx
        } else
        {
          level_to_plot <- 1
        }

        display_name <- gsub("size_", "Size ", level_names[level_to_plot])

        plot_spatial_clusters(
          points = points_combined,
          hulls = clusters_chulls[[level_to_plot]],
          level_name = display_name
        )
      }
    }

    # ----------- SAR PLOT -----------
    if (plot_type == "sar")
    {
      sar_results <- result[["sar_analysis"]]

      if (sar_results$valid)
      {
        plot_sar_clusters(sar_results)
      } else
      {
        cat("SAR analysis not valid.\n")
        cat("Message:", sar_results$message, "\n")
      }
    }

    # ----------- CSAR HEATMAP -----------
    if (plot_type == "csar")
    {
      if (has_csar)
      {
        print(plot_heatmap(result[["csar_analysis"]]))
      } else
      {
        cat("CSAR analysis not available in this result object.\n")
      }
    }
  }

  # Return invisible
  invisible(result)
}
