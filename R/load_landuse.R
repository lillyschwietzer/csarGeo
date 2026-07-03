#' Load the 1995 Land Use Raster
#'
#' Loads the built-in 1995 land use raster (coa95_raster_20.tif) default data of
#' the csarGeo package. The raster contains land use/land cover classification
#' for the study area with 20m resolution and EPSG:3763 (Portugal TM06) CRS.
#'
#' @return A SpatRaster object (from the \code{terra} package) containing the
#'   1995 land use classification with the following habitat codes:
#'   \itemize{
#'     \item 1: Forest
#'     \item 2: Agriculture
#'     \item 3: Shrubland
#'     \item NA: Other / Unclassified
#'   }
#'
#' @examples
#' \dontrun{
#' # Load the raster
#' lu1995 <- load_lu1995()
#'
#' # Plot the raster
#' terra::plot(lu1995, main = "Land Use 1995")
#'
#' # Use in analysis
#' res <- countryside_sar(
#'   data = species_data,
#'   habitat = lu1995,
#'   habitat_names = c("Forest", "Agriculture", "Shrubland"),
#'   method = "clusters"
#' )
#' }
#'
#' @export
load_lu1995 <- function() {
  file_path <- system.file("extdata", "coa95_raster_20.tif", package = "csarGeo")
  if (file_path == "") {
    stop("Raster file not found. Please re-install the package.")
  }
  terra::rast(file_path)
}
