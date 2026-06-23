#' Visualize countryside_sar results
#'
#' @description
#' A function to visualize the results of \code{countryside_sar()} from the csarGeo-Package. It can generate three plot types: "map", a plot of the polygon hull with the sampling locations and the chosen sampling method, "sar" for the base SAR analysis and "csar", only available for a result of method "clusters" to create an affinity heat map of the species groups to different habitats.
#'
#'
#' @param result Output object from countryside_sar function.
#' @param plot_type Defines the type of plot. Three possible options
#'   \itemize{
#'     \item \code{map} (map of the polygon hull + sampling locations with the chosen sampling method)
#'     \item \code{sar} (log-log scaled linear regression)
#'     \item \code{csar} (heatmap of species group affinities to habitat types)
#'     }
#' \code{map} and \code{sar} are available for results of both \code{method = circles} and \code{method = clusters}. \code{plot_type = csar} is only available for results of \code{method = clusters}.
#' @param plot_all_runs Plot all runs of "circles" into one grid. Defaults to TRUE.
#' @param plot_run_n Plot only the result of the n-th run of "circles". Defaults to NULL.
#' @param plot_all_levels Plot all clustering levels of method "clusters" into one grid. Defaults to TRUE.
#' @param plot_level_n Plot only the n-th clustering level as a solo plot. Defaults to NULL.
#'
#' @return Creates either a single plot or a grid of multiple plots.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_countryside_sar(res, plot_type = "map")
#' }
visuals_sar <- function(result,
                       plot_type = NULL, # "map", "sar", "csar",
                       # Circles
                       plot_all_runs = TRUE,
                       plot_run_n = NULL,
                       # Clusters
                       plot_all_levels = TRUE,
                       plot_level_n = NULL) {

  #----------------- Input validation ------------------
  if (!plot_type %in% c("map", "sar", "csar")) {
    stop("plot_type must be 'map', 'sar', or 'csar'")
  }

  method <- result[["method"]]

  # Check if csar exists for heatmap
  has_csar <- !is.null(result[["csar_analysis"]]) &&
    !is.null(result[["csar_analysis"]]$valid) &&
    result[["csar_analysis"]]$valid

  if (plot_type == "csar" && !has_csar) {
    stop("CSAR analysis not available in this result object.")
  }

  #------------------------ 1) Helper Functions ---------------------------

  # Circles data preparation
  prepare_circles_for_plot <- function(res_circles, run_number = 1) {
    samples <- res_circles[["runs"]][[run_number]][["samples"]]
    circle_polygons <- lapply(samples, `[[`, "circle")

    return(list(
      points = res_circles[["points_sf"]],
      circles = circle_polygons,
      hull = res_circles[["convex_hull"]]
    ))
  }

  if (method == "circles") {

    # Circles Map
    plot_spatial_circles <- function(points_sf, circles, convex_hull, main_title = NULL) {
      if (is.null(main_title)) main_title <- "Sampling Circles"
      plot(sf::st_geometry(points_sf), main = main_title)
      plot(convex_hull, border = "red", lwd = 2, add = TRUE)
      for (circle in circles) {
        plot(circle, border = "blue", add = TRUE)
      }
    }

    # Circles SAR
    plot_sar_circles <- function(sar_results, main_title = NULL) {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Circles"

      log_area <- sar_results[["log_area"]]
      log_sp <- sar_results[["log_sp"]]

      plot(log_area, log_sp,
           xlab = "log(Area)",
           ylab = "log(Species Richness)",
           main = main_title,
           pch = 16)

      abline(sar_results[["lm_model"]], col = "red", lwd = 2)

      r2 <- sar_results[["lm_summary"]][["r.squared"]]
      legend("bottomright", legend = paste("R² =", round(r2, 3)), bty = "n")
    }

  } else { # method == "clusters"

    # Clusters Map
    plot_spatial_clusters <- function(points, hulls, level_name = NULL, main_title = NULL) {
      if (is.null(main_title)) {
        main_title <- if (!is.null(level_name)) level_name else "Clustering Pattern"
      }

      plot(sf::st_geometry(points),
           main = main_title,
           cex = 0.5, pch = 16)

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

    # Clusters SAR
    plot_sar_clusters <- function(sar_results, main_title = NULL) {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Clusters"

      log_area <- sar_results[["log_area"]]
      log_sp <- sar_results[["log_sp"]]

      plot(log_area, log_sp,
           xlab = "log(Area)",
           ylab = "log(Species Richness)",
           main = main_title,
           pch = 16)

      abline(sar_results[["lm_model"]], col = "red", lwd = 2)

      r2 <- sar_results[["lm_summary"]][["r.squared"]]
      legend("bottomright", legend = paste("R² =", round(r2, 3)), bty = "n")
    }

    # Affinity Heatmap
    plot_heatmap <- function(csar_results) {

      # Extract and convert affinity data
      affinity_list <- csar_results$model$affinity
      affinity_mat <- do.call(rbind, affinity_list)

      # Convert to long format for ggplot
      affinity_long <- reshape2::melt(affinity_mat,
                                      varnames = c("Species_Group", "Habitat"),
                                      value.name = "Affinity")

      # Create heatmap
      ggplot2::ggplot(affinity_long,
                      ggplot2::aes(x = Habitat, y = Species_Group, fill = Affinity)) +
        ggplot2::geom_tile(color = "white", size = 0.5) +
        ggplot2::scale_fill_viridis_c(
          name = "Affinity",
          limits = c(0, 1),
          option = "viridis",
          direction = -1
        ) +
        ggplot2::geom_text(ggplot2::aes(label = round(Affinity, 3)),
                           color = "white", size = 3.5) +
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

  if (method == "circles") {

    n_runs <- result[["n_runs"]]

    # ----------- Map Plot -----------
    if (plot_type == "map") {

      if (plot_all_runs == TRUE && n_runs > 1) {
        # Plot all runs in grid
        n_cols <- ceiling(sqrt(n_runs))
        n_rows <- ceiling(n_runs / n_cols)
        par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))

        for (run in 1:n_runs) {
          circles_data <- prepare_circles_for_plot(result, run)
          plot_spatial_circles(
            points_sf = circles_data$points,
            circles = circles_data$circles,
            convex_hull = circles_data$hull,
            main_title = paste("Run", run)
          )
        }
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

      } else {
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
    if (plot_type == "sar") {
      run_to_plot <- if (!is.null(plot_run_n)) plot_run_n else 1
      sar_results <- result[["runs"]][[run_to_plot]][["sar_analysis"]]

      if (sar_results$valid) {
        plot_sar_circles(sar_results, main_title = paste("SAR - Run", run_to_plot))
      } else {
        cat("SAR analysis not valid for run", run_to_plot, "\n")
        cat("Message:", sar_results$message, "\n")
      }
    }

  } else { # method == "clusters"

    # Extract data
    points_combined <- do.call(rbind, result[["samples"]][["size_1"]][["points"]])
    clusters_chulls <- result[["clusters_chulls"]]
    n_levels <- length(clusters_chulls)

    # ----------- MAP PLOT -----------
    if (plot_type == "map") {

      if (plot_all_levels == TRUE && n_levels > 1) {
        # Plot all levels in grid
        n_cols <- ceiling(sqrt(n_levels))
        n_rows <- ceiling(n_levels / n_cols)
        par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))

        level_names <- names(clusters_chulls)
        for (level in 1:n_levels) {
          plot_spatial_clusters(
            points = points_combined,
            hulls = clusters_chulls[[level]],
            level_name = level_names[level]
          )
        }
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

      } else {
        # Plot single level
        level_to_plot <- if (!is.null(plot_level_n)) plot_level_n else 1
        level_names <- names(clusters_chulls)
        plot_spatial_clusters(
          points = points_combined,
          hulls = clusters_chulls[[level_to_plot]],
          level_name = level_names[level_to_plot]
        )
      }
    }

    # ----------- SAR PLOT -----------
    if (plot_type == "sar") {
      sar_results <- result[["sar_analysis"]]

      if (sar_results$valid) {
        plot_sar_clusters(sar_results)
      } else {
        cat("SAR analysis not valid.\n")
        cat("Message:", sar_results$message, "\n")
      }
    }

    # ----------- CSAR HEATMAP -----------
    if (plot_type == "csar") {
      if (has_csar) {
        print(plot_heatmap(result[["csar_analysis"]]))
      } else {
        cat("CSAR analysis not available in this result object.\n")
      }
    }
  }

  # Return invisible
  invisible(result)
}
