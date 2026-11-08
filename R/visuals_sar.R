#' Visualize countryside_sar() function results
#'
#' @description
#' A function to visualize the results of \code{countryside_sar()} from the csarGeo-Package.
#' It can generate six plot types: "map", a plot of the polygon hull with the
#' sampling locations and the chosen sampling method; "sar" for the base SAR analysis;
#' and, only for \code{method = "clusters"}, "csar" for effective area scaling,
#' "fit" for a predicted vs. observed scatterplot, "affinities" for a heatmap of
#' species-group affinities to habitats, and "habitat_proportion" for richness
#' along a single-habitat proportion gradient.
#' The SAR results can be plotted as a single run result, a grid of multiple runs or
#' as an average SAR. The plots of the sampling maps can either be done per level or run
#' or in a grid.
#'
#' @param result Output object from countryside_sar function() from the csarGeo-Package.
#' @param plot_type Defines the type of plot. Six possible options:
#'   \itemize{
#'     \item \code{map} (Map of the polygon hull + sampling locations with the chosen sampling method. Adds sampling circles for method = 'circles' or plots the clusters)
#'     \item \code{sar} (log-log scaled linear regression of total area vs total richness)
#'     \item \code{csar} (log-log effective area scaling per species group; clusters only)
#'     \item \code{fit} (predicted vs. observed total richness scatterplot; clusters only)
#'     \item \code{affinities} (heatmap of species group affinities to habitat types; clusters only)
#'     \item \code{habitat_proportion} (predicted richness along a focal habitat proportion gradient; if \code{focal_habitat = NULL}, a grid across all habitats; clusters only)
#'     }
#' \code{map} and \code{sar} are available for results of both \code{method = "circles"} and
#' \code{method = "clusters"}. All other types are only available for \code{method = "clusters"}.
#' @param plot_all_runs Plot all runs of "circles" into one grid if n_runs > 1. Defaults to TRUE.
#' @param plot_run_n Plot only the result of the n-th run of "circles". Defaults to NULL.
#' @param plot_average If TRUE, plots the average SAR result if n_runs > 1. Defaults to FALSE.
#' @param plot_all_levels Plot all clustering levels of method "clusters" into one grid. Defaults to TRUE.
#' @param plot_level_n Plot only the n-th clustering level as a solo plot. Defaults to NULL.
#' @param focal_habitat For \code{plot_type = "habitat_proportion"}, the habitat whose
#'   proportion is varied along the x-axis. If \code{NULL} (default), a grid of all
#'   habitats is produced.
#' @param total_area For \code{plot_type = "habitat_proportion"}, the fixed total area at
#'   which habitat proportions are varied. Defaults to the mean Area_Total of the data.
#'
#' @return Creates either a single plot or a grid of multiple plots. The function returns
#'   the input result object invisibly.
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
#' # cSAR effective area scaling
#' visuals_sar(res, plot_type = "csar")
#'
#' # cSAR model fit (predicted vs observed)
#' visuals_sar(res, plot_type = "fit")
#'
#' # Affinity heatmap
#' visuals_sar(res, plot_type = "affinities")
#'
#' # Habitat proportion gradient (grid of all habitats)
#' visuals_sar(res, plot_type = "habitat_proportion")
#'
#' # Habitat proportion gradient (single focal habitat)
#' visuals_sar(res, plot_type = "habitat_proportion", focal_habitat = "Forest")
#' }
visuals_sar <- function(result,
                        plot_type = NULL, # "map", "sar", "csar", "fit", "affinities", "habitat_proportion"
                        # Circles
                        plot_all_runs = TRUE,
                        plot_run_n = NULL,
                        plot_average = FALSE,
                        # Clusters
                        plot_all_levels = TRUE,
                        plot_level_n = NULL,
                        # cSAR
                        focal_habitat = NULL,
                        total_area = NULL) {

  #--------------------------- Input validation --------------------------------
  valid_types <- c("map", "sar", "csar", "fit", "affinities", "habitat_proportion")
  if (!plot_type %in% valid_types)
  {
    stop("plot_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  method <- result[["method"]]

  # Check if SAR exists (for circles and clusters)
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

  if (plot_type %in% c("csar", "fit", "affinities", "habitat_proportion")) {
    if (method != "clusters") {
      stop("plot_type '", plot_type, "' is only available for method = 'clusters'.")
    }
    if (!has_csar) {
      stop("cSAR analysis not available in this result object.")
    }
  }

  #-------------------------- 1) Helper Functions ------------------------------

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

  #------- Extract valid h params and habitat columns for csar plots -----------
  extract_csar_layout <- function(res)
  {
    fits <- res[["csar_analysis"]][["model"]][["fits"]]
    cf0  <- stats::coef(fits[[1]])
    valid_h <- cf0[grep("^h[0-9]+$", names(cf0))]
    n_hab   <- sum(!is.na(valid_h))
    habitats <- names(res[["results_table"]])[1:n_hab]
    list(fits = fits, n_hab = n_hab, habitats = habitats)
  }

  # Per-group prediction helper for fit/csar plots
  predict_all_groups <- function(res, area_mat, fits, n_hab)
  {
    sapply(names(fits), function(grp) {
      cf <- stats::coef(fits[[grp]])
      h  <- cf[grep("^h[0-9]+$", names(cf))]; h <- h[!is.na(h)]
      z  <- cf[["z"]]
      if (length(h) != n_hab || is.na(z)) return(rep(0, nrow(area_mat)))
      eff <- as.numeric(area_mat %*% h)
      ifelse(eff > 0, eff^z, 0)
    })
  }

  # ---------------- cSAR fit: predicted vs. observed scatterplot ---------------
  plot_csar_observed_vs_predicted <- function(res)
  {
    rt   <- res[["results_table"]]
    lay  <- extract_csar_layout(res)
    keep <- rt$Sp_Total > 0 & is.finite(rt$Sp_Total)
    obs  <- rt$Sp_Total[keep]
    area_mat <- as.matrix(rt[keep, lay$habitats, drop = FALSE])

    pred_mat   <- predict_all_groups(res, area_mat, lay$fits, lay$n_hab)
    pred_total <- rowSums(pred_mat, na.rm = TRUE)
    r2     <- summary(stats::lm(obs ~ pred_total))$r.squared
    rmse   <- sqrt(mean((obs - pred_total)^2))
    slope  <- coef(stats::lm(obs ~ pred_total))[2]

    rng <- range(c(obs, pred_total), na.rm = TRUE)
    plot(pred_total, obs,
         xlab = "Predicted Total Richness",
         ylab = "Observed Total Richness",
         main = "cSAR - Predicted vs. Observed",
         pch = 16, xlim = rng, ylim = rng,
         col = adjustcolor("black", alpha.f = 0.5))
    grid()
    abline(0, 1, col = "red", lwd = 2, lty = 2)
    legend("topleft",
           legend = paste0(
             "R\u00B2 = ", round(r2, 3), "\n",
             "RMSE = ", round(rmse, 2), "\n",
             "Slope = ", round(slope, 2)
           ),
           bty = "n", cex = 0.9)
    invisible(list(observed = obs, predicted = pred_total, r_squared = r2))
  }

  # ---------------- cSAR: affinity heatmap ----------------
  plot_heatmap <- function(csar_results)
  {
    affinity_list <- csar_results$model$affinity
    affinity_mat  <- do.call(rbind, affinity_list)
    affinity_long <- reshape2::melt(affinity_mat,
                                    varnames   = c("Species_Group", "Habitat"),
                                    value.name = "Affinity")
    affinity_long$text_color <- ifelse(affinity_long$Affinity > 0.5, "white", "black")

    ggplot2::ggplot(affinity_long,
                    ggplot2::aes(x = Habitat, y = Species_Group, fill = Affinity)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_viridis_c(name = NULL, limits = c(0, 1),
                                    option = "viridis", direction = -1,
                                    guide = "none") +
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
        plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
        panel.grid  = ggplot2::element_blank()
      )
  }

  # ----------------------- cSAR: effective area -------------------------------
  plot_csar_effective_area <- function(res)
  {
    rt  <- res[["results_table"]]
    lay <- extract_csar_layout(res)
    area_mat <- as.matrix(rt[, lay$habitats, drop = FALSE])
    eff_mat  <- predict_all_groups(res, area_mat, lay$fits, lay$n_hab)
    colnames(eff_mat) <- names(lay$fits)

    df_list <- list()
    for (grp in colnames(eff_mat)) {
      sp_col <- grp
      if (!sp_col %in% names(rt)) sp_col <- gsub("^Sp_", "", grp)

      obs  <- rt[[sp_col]]
      eff  <- eff_mat[, grp]
      keep <- obs > 0 & eff > 0 & is.finite(obs) & is.finite(eff)

      df_list[[grp]] <- data.frame(
        Group          = grp,
        LogEffArea     = log(eff[keep]),
        LogObsRichness = log(obs[keep])
      )
    }
    df <- do.call(rbind, df_list)
    df$Group <- factor(df$Group, levels = colnames(eff_mat))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = LogEffArea, y = LogObsRichness)) +
      ggplot2::geom_point(alpha = 0.5, size = 1.5, color = "black") +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                           color = "red", se = FALSE, linewidth = 1) +
      ggplot2::facet_wrap(~ Group, scales = "free") +
      ggplot2::labs(
        title    = "countryside SAR - Effective Area",
        subtitle = "Effective Area = \u03A3 (affinity \u00D7 habitat area) for each species group",
        x        = "log(Effective Area)",
        y        = "log(Observed Species Richness)"
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30", size = 10),
        strip.text    = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing.x  = ggplot2::unit(2.5, "lines"),
        panel.spacing.y  = ggplot2::unit(1.5, "lines"),
        plot.margin      = ggplot2::margin(15, 15, 15, 15)
      )

    stat_df <- do.call(rbind, lapply(levels(df$Group), function(g) {
      sub <- df[df$Group == g, ]
      fit <- stats::lm(LogObsRichness ~ LogEffArea, data = sub)
      sm  <- summary(fit)
      r2  <- sm$r.squared
      slope <- sm$coefficients[2, 1]
      se    <- sm$coefficients[2, 2]

      data.frame(
        Group = g,
        Label = paste0("R\u00B2 = ", round(r2, 2), "\n",
                       "z = ", round(slope, 2), " \u00B1 ", round(se, 2)),
        x     = min(sub$LogEffArea),
        y     = max(sub$LogObsRichness)
      )
    }))

    p <- p + ggplot2::geom_text(
      data = stat_df,
      ggplot2::aes(x = x, y = y, label = Label),
      hjust = 0, vjust = 1,
      inherit.aes = FALSE,
      size = 4, fontface = "italic"
    )

    print(p)
    invisible(eff_mat)
  }

  # ------------- cSAR: habitat proportion gradient ----------------------------
  plot_csar_habitat_proportion <- function(res, focal_habitat = NULL, total_area = NULL)
  {
    rt  <- res[["results_table"]]
    lay <- extract_csar_layout(res)

    if (!is.null(focal_habitat) && !focal_habitat %in% lay$habitats)
      stop("focal_habitat must be one of: ", paste(lay$habitats, collapse = ", "))
    if (is.null(total_area)) total_area <- max(rt$Area_Total, na.rm = TRUE)

    props  <- seq(0, 1, length.out = 101)
    groups <- names(lay$fits)
    df_list <- list()

    for (hab in lay$habitats) {
      others <- setdiff(lay$habitats, hab)

      area_grid <- matrix(0, nrow = length(props), ncol = lay$n_hab,
                          dimnames = list(NULL, lay$habitats))
      area_grid[, hab] <- props * total_area
      for (o in others) area_grid[, o] <- (1 - props) * total_area / length(others)

      pred_grid <- predict_all_groups(res, area_grid, lay$fits, lay$n_hab)

      for (grp in groups) {
        df_list[[paste(hab, grp)]] <- data.frame(
          Habitat    = hab,
          Group      = grp,
          Proportion = props,
          Richness   = pred_grid[, grp]
        )
      }
    }

    df <- do.call(rbind, df_list)

    if (!is.null(focal_habitat)) {
      df <- df[df$Habitat == focal_habitat, ]
      title_text <- paste0("cSAR - Richness vs. ", focal_habitat, " Proportion")
    } else {
      title_text <- "cSAR - Richness vs. Habitat Proportion"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Proportion, y = Richness, color = Group)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::facet_wrap(~ Habitat, scales = "fixed") +
      ggplot2::scale_color_viridis_d(option = "viridis") +
      ggplot2::labs(
        title    = title_text,
        subtitle = "Effect of habitat proportion on predicted species richness",
        x        = "Habitat Proportion",
        y        = "Predicted Species Richness",
        color    = "Species Group"
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30", size = 10),
        strip.text    = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing.x  = ggplot2::unit(2.5, "lines"),
        panel.spacing.y  = ggplot2::unit(1.5, "lines"),
        plot.margin      = ggplot2::margin(15, 15, 15, 15)
      )

    print(p)
    invisible(df)
  }

  if (method == "circles")
  {
    # Circles Map
    plot_spatial_circles <- function(points_sf, circles, convex_hull, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Sampling Circles"
      plot(sf::st_geometry(points_sf), main = main_title)
      plot(convex_hull, border = "red", lwd = 2, add = TRUE)
      for (circle in circles) plot(circle, border = "blue", add = TRUE)
    }

    # Circles SAR - Single run (solid black points)
    plot_sar_circles <- function(sar_results, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Circles"
      log_area <- sar_results[["log_area"]]
      log_sp   <- sar_results[["log_sp"]]
      intercept <- coef(sar_results[["lm_model"]])[1]
      slope     <- coef(sar_results[["lm_model"]])[2]
      r_squared <- sar_results[["lm_summary"]]$r.squared

      plot(log_area, log_sp, xlab = "log(Area)", ylab = "log(Species Richness)",
           main = main_title, pch = 16, col = "black")   # <-- solid black
      grid()
      abline(sar_results[["lm_model"]], col = "red", lwd = 2)
      legend("topleft",
             legend = paste0("log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
                             "Slope (z) = ", round(slope, 3), "\n",
                             "R² = ", round(r_squared, 3)),
             bty = "n", cex = 0.9)
    }

    # Circles SAR - Multiple runs in grid (solid black points)
    plot_sar_circles_grid <- function(sar_results, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "SAR"
      log_area <- sar_results[["log_area"]]
      log_sp   <- sar_results[["log_sp"]]
      intercept <- coef(sar_results[["lm_model"]])[1]
      slope     <- coef(sar_results[["lm_model"]])[2]
      r_squared <- sar_results[["lm_summary"]]$r.squared

      plot(log_area, log_sp, xlab = "", ylab = "", main = main_title,
           pch = 16, axes = FALSE, col = "black")   # <-- solid black
      box(); axis(1, labels = FALSE); axis(2, labels = FALSE)
      grid(); abline(sar_results[["lm_model"]], col = "red", lwd = 2)
      legend("topleft",
             legend = paste0("log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
                             "z = ", round(slope, 3), "\n",
                             "R² = ", round(r_squared, 3)),
             bty = "n", cex = 0.8, inset = c(-0.02, -0.02))
    }

    # Average Circles SAR (with opacity, no CI band)
    plot_average_sar_circles <- function(avg_sar_results, runs, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Average of Runs"
      intercept <- avg_sar_results[["avg_intercept"]]
      slope     <- avg_sar_results[["avg_slope"]]
      r_squared <- avg_sar_results[["avg_r_squared"]]
      n_valid   <- avg_sar_results[["n_valid_runs"]]
      sd_slope  <- avg_sar_results[["sd_slope"]]

      all_log_area <- unlist(lapply(runs, function(x) x$sar_analysis$log_area))
      all_log_sp   <- unlist(lapply(runs, function(x) x$sar_analysis$log_sp))

      plot(all_log_area, all_log_sp,
           xlab = "log(Area)", ylab = "log(Species Richness)", main = main_title,
           pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 0.8)  # <-- opacity kept
      grid()

      x_vals <- range(all_log_area, na.rm = TRUE)
      y_vals <- intercept + slope * x_vals
      lines(x_vals, y_vals, col = "red", lwd = 3)

      legend("topleft",
             legend = paste0("log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
                             "Avg slope (z) = ", round(slope, 3), " \u00B1 ", round(sd_slope, 3), "\n",
                             "Avg R\u00B2 = ", round(r_squared, 3), "\n",
                             "Valid runs = ", n_valid),
             bty = "n", cex = 0.8)
    }

  } else
  { # method == "clusters"

    # Clusters Map
    plot_spatial_clusters <- function(points, hulls, level_name = NULL, main_title = NULL)
    {
      if (is.null(main_title))
        main_title <- if (!is.null(level_name)) level_name else "Clustering Pattern"
      plot(sf::st_geometry(points), main = main_title, cex = 0.5, pch = 16)
      if (length(hulls) > 0) {
        colors <- rainbow(length(hulls))
        for (i in seq_along(hulls)) {
          plot(hulls[[i]], border = "black",
               col = adjustcolor(colors[i], alpha.f = 0.5), add = TRUE)
        }
      }
    }

    # Clusters SAR
    plot_sar_clusters <- function(sar_results, main_title = NULL)
    {
      if (is.null(main_title)) main_title <- "Species-Area Relationship (SAR) - Clusters"
      log_area  <- sar_results[["log_area"]]
      log_sp    <- sar_results[["log_sp"]]
      intercept <- coef(sar_results[["lm_model"]])[1]
      slope     <- coef(sar_results[["lm_model"]])[2]
      r_squared <- sar_results[["lm_summary"]]$r.squared

      plot(log_area, log_sp, xlab = "log(Area)", ylab = "log(Species Richness)",
           main = main_title, pch = 16,
           col = adjustcolor("black", alpha.f = 0.5))
      grid()
      abline(sar_results[["lm_model"]], col = "red", lwd = 2)
      legend("topleft",
             legend = paste0("log(S) = ", round(intercept, 3), " + ", round(slope, 3), "x\n",
                             "Slope (z) = ", round(slope, 3), "\n",
                             "R² = ", round(r_squared, 3)),
             bty = "n", cex = 0.8)
    }
  }

  #----------------------------- 2) Main Processing ----------------------------

  if (method == "circles")
  {
    n_runs <- result[["n_runs"]]

    # ----------- Map Plot -----------
    if (plot_type == "map") {
      if (plot_all_runs == TRUE && is.null(plot_run_n) && n_runs > 1)
        {
        n_cols <- ceiling(sqrt(n_runs)); n_rows <- ceiling(n_runs / n_cols)
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

    # ----------- SAR Plot -----------
    if (plot_type == "sar") {
      if (plot_average == TRUE) {
        avg_sar_results <- result[["avg_sar_results"]]
        runs <- result[["runs"]]
        if (avg_sar_results$valid) {
          plot_average_sar_circles(avg_sar_results, runs)
        } else {
          cat("Average SAR analysis not valid\n")
          cat("Message:", avg_sar_results$message, "\n")
        }
      } else if (!is.null(plot_run_n)) {
        sar_results <- result[["runs"]][[plot_run_n]][["sar_analysis"]]
        if (sar_results$valid) {
          plot_sar_circles(sar_results, main_title = paste("SAR - Run", plot_run_n))
        } else {
          cat("SAR analysis not valid for run", plot_run_n, "\n")
          cat("Message:", sar_results$message, "\n")
        }
      } else if (n_runs > 1 && is.null(plot_run_n)) {
        n_cols <- ceiling(sqrt(n_runs)); n_rows <- ceiling(n_runs / n_cols)
        par(mfrow = c(n_rows, n_cols), mar = c(0.5, 0.5, 2, 0.5), oma = c(4, 4, 2, 2))
        for (run in 1:n_runs) {
          sar_results <- result[["runs"]][[run]][["sar_analysis"]]
          if (sar_results$valid) {
            plot_sar_circles_grid(sar_results, main_title = paste("Run", run))
          } else {
            plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
            text(1, 1, paste("Run", run, "\nInvalid"), col = "red", cex = 1.2)
          }
        }
        mtext("log(Area)", side = 1, outer = TRUE, line = 2.5, cex = 1.2)
        mtext("log(Species Richness)", side = 2, outer = TRUE, line = 2.5, cex = 1.2)
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))
      } else {
        run_to_plot <- if (!is.null(plot_run_n)) plot_run_n else 1
        sar_results <- result[["runs"]][[run_to_plot]][["sar_analysis"]]
        if (sar_results$valid) {
          plot_sar_circles(sar_results, main_title = paste("SAR - Run", run_to_plot))
        } else {
          cat("SAR analysis not valid for run", run_to_plot, "\n")
          cat("Message:", sar_results$message, "\n")
        }
      }
    }

  } else
  { # method == "clusters"

    # Extract data
    clusters_chulls <- result[["clusters_chulls"]]
    n_levels        <- length(clusters_chulls)
    level_names     <- names(clusters_chulls)

    # ----------- Map Plot -----------
    if (plot_type == "map") {
      if (plot_all_levels == TRUE && is.null(plot_level_n) && n_levels > 1) {
        n_cols <- ceiling(sqrt(n_levels)); n_rows <- ceiling(n_levels / n_cols)
        par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))
        for (level in 1:n_levels) {
          display_name <- gsub("size_", "Size ", level_names[level])
          pts <- do.call(rbind, result[["samples"]][[level]][["points"]])
          plot_spatial_clusters(
            points = pts,
            hulls = clusters_chulls[[level]],
            level_name = display_name
          )
        }
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
      } else {
        if (!is.null(plot_level_n)) {
          target_name_a <- paste0("size_", plot_level_n)
          target_name_b <- as.character(plot_level_n)
          level_idx <- which(level_names %in% c(target_name_a, target_name_b))
          if (length(level_idx) == 0) {
            stop("Cluster size '", plot_level_n, "' not found. Available sizes: ",
                 paste(gsub("size_", "", level_names), collapse = ", "))
          }
          level_to_plot <- level_idx
        } else {
          level_to_plot <- 1
        }
        display_name <- gsub("size_", "Size ", level_names[level_to_plot])
        if (display_name == level_names[level_to_plot]) {
          display_name <- paste0("Size ", level_names[level_to_plot])
        }
        pts <- do.call(rbind, result[["samples"]][[level_to_plot]][["points"]])
        plot_spatial_clusters(
          points = pts,
          hulls = clusters_chulls[[level_to_plot]],
          level_name = display_name
        )
      }
    }

    # ----------- SAR Plot -----------
    if (plot_type == "sar") {
      sar_results <- result[["sar_analysis"]]
      if (sar_results$valid) {
        plot_sar_clusters(sar_results)
      } else {
        cat("SAR analysis not valid.\n")
        cat("Message:", sar_results$message, "\n")
      }
    }

    # ----------- cSAR Effective Area -----------
    if (plot_type == "csar") {
      plot_csar_effective_area(result)
    }

    # ----------- cSAR Observed vs Predicted -----------
    if (plot_type == "fit") {
      plot_csar_observed_vs_predicted(result)
    }

    # ----------- Affinities Heatmap -----------
    if (plot_type == "affinities") {
      print(plot_heatmap(result[["csar_analysis"]]))
    }

    # ----------- Habitat Proportion -----------
    if (plot_type == "habitat_proportion") {
      plot_csar_habitat_proportion(result, focal_habitat, total_area)
    }
  }

  # Return invisible
  invisible(result)
}
