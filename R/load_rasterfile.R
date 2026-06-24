#' Load Example Data SpatRaster
#'
#' Downloads land use raster "coa95_raster_20.tif", which contains data from 1995, from the GitHub release of the csarGeo package. It caches the raster locally to prevent re-downloading.
#'
#' @param force Logical value. If \code{force = TRUE}, force re-download even if the file already exists in the cache. \code{FALSE} by default.
#'
#' @returns A SpatRaster object (from the \code{terra} package) containing the 1995 land use/land cover classification with 20m resolution and EPSG:3763.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the raster (downloads if not cached)
#' land_use <- load_rasterfile()
#'
#' # Force re-download
#' land_use <- load_rasterfile(force = TRUE)
#' }
load_rasterfile <- function(force = FALSE)
{
  cache_dir <- tools::R_user_dir("csarGeo", which = "cache")
  lu <- file.path(cache_dir, "coa95_raster_20.tif")

  # Only download if file doesn't exist OR force = TRUE
  if (!file.exists(lu) || force)
  {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    download.file(
      url = "https://github.com/lillyschwietzer/csarGeo/releases/download/v1.0.0_raster/coa95_raster_20.tif",
      destfile = lu,
      mode = "wb"
    )
  }

  terra::rast(lu)
}
