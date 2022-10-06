# Copyright (C) 2021-2022 Koen Derks

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Draw a Fractal Flame
#'
#' @description This function implements the fractal flame algorithm.
#'
#' @usage canvas_flame(colors, background = "#fafafa",
#'              iterations = 1000000, zoom = 1, resolution = 1000,
#'              variations = NULL, blend = TRUE,
#'              post = FALSE, final = FALSE, extra = FALSE,
#'              verbose = FALSE)
#'
#' @param colors      a string or character vector specifying the color(s) used for the artwork.
#' @param background  a character specifying the color used for the background.
#' @param iterations  a positive integer specifying the number of iterations of the algorithm.
#' @param zoom        a positive value specifying the amount of zooming.
#' @param resolution  resolution of the artwork in pixels per row/column. Increasing the resolution increases the quality of the artwork but also increases the computation time exponentially.
#' @param variations  an integer (vector) specifying the variations to be included in the flame. The default \code{NULL} includes a random number of variations. See the details section for more information about possible variations.
#' @param blend       logical. Whether to blend the variations (\code{TRUE}) or pick a unique variation in each iteration (\code{FALSE}).
#' @param post        logical. Whether to apply a post transformation in each iteration.
#' @param final       logical. Whether to apply a final transformation in each iteration.
#' @param extra       logical. Whether to apply an additional post transformation after the final transformation. Only has an effect when \code{final = TRUE}.
#' @param verbose     logical. Whether to print information.
#'
#' @details           The \code{variation} argument can be used to include specific variations into the flame. See the appendix in the references for examples of all variations. Possible variations are:
#'
#' \itemize{
#'  \item{\code{0}: Linear}
#'  \item{\code{1}: Sine}
#'  \item{\code{2}: Spherical}
#'  \item{\code{3}: Swirl}
#'  \item{\code{4}: Horsehoe}
#'  \item{\code{5}: Polar}
#'  \item{\code{6}: Handkerchief}
#'  \item{\code{7}: Heart}
#'  \item{\code{8}: Disc}
#'  \item{\code{9}: Spiral}
#'  \item{\code{10}: Hyperbolic}
#'  \item{\code{11}: Diamond}
#'  \item{\code{12}: Ex}
#'  \item{\code{13}: Julia}
#'  \item{\code{14}: Bent}
#'  \item{\code{15}: Waves}
#'  \item{\code{16}: Fisheye}
#'  \item{\code{17}: Popcorn}
#'  \item{\code{18}: Exponential}
#'  \item{\code{19}: Power}
#'  \item{\code{20}: Cosine}
#'  \item{\code{21}: Rings}
#'  \item{\code{22}: Fan}
#'  \item{\code{23}: Blob}
#'  \item{\code{24}: PDJ}
#'  \item{\code{25}: Fan2}
#'  \item{\code{26}: Rings2}
#'  \item{\code{27}: Eyefish}
#'  \item{\code{28}: Bubble}
#'  \item{\code{29}: Cylinder}
#'  \item{\code{30}: Perspective}
#'  \item{\code{31}: Noise}
#'  \item{\code{32}: JuliaN}
#'  \item{\code{33}: JuliaScope}
#'  \item{\code{34}: Blur}
#'  \item{\code{35}: Gaussian}
#'  \item{\code{36}: RadialBlur}
#'  \item{\code{37}: Pie}
#'  \item{\code{38}: Ngon}
#'  \item{\code{39}: Curl}
#'  \item{\code{40}: Rectangles}
#'  \item{\code{41}: Arch}
#'  \item{\code{42}: Tangent}
#'  \item{\code{43}: Square}
#'  \item{\code{44}: Rays}
#' }
#'
#' @return A \code{ggplot} object containing the artwork.
#'
#' @references \url{https://flam3.com/flame_draves.pdf}
#'
#' @author Koen Derks, \email{koen-derks@hotmail.com}
#'
#' @keywords artwork canvas
#'
#' @seealso \code{colorPalette}
#'
#' @examples
#' \donttest{
#' set.seed(2)
#'
#' # Simple example
#' canvas_flame(colors = colorPalette("origami"))
#'
#' # Advanced example (no-blend sine and spherical variations)
#' canvas_flame(colors = colorPalette("origami"), variations = c(1, 2), blend = FALSE)
#' }
#'
#' @export

canvas_flame <- function(colors, background = "#fafafa",
                         iterations = 1000000, zoom = 1, resolution = 1000,
                         variations = NULL, blend = TRUE,
                         post = FALSE, final = FALSE, extra = FALSE,
                         verbose = FALSE) {
  .checkUserInput(
    resolution = resolution, background = background
  )
  iterations <- iterations + 20
  varNames <- .getVariationNames()
  noVariations <- length(varNames)
  user <- FALSE
  if (is.null(variations)) {
    user <- TRUE
    v <- 0:(noVariations - 1)
    variations <- sample(x = v, size = sample(1:5, size = 1), replace = FALSE)
  } else if (min(variations) < 0 || max(variations) > (noVariations - 1)) {
    stop("'variations' must be between 0 and ", (noVariations - 1))
  }
  if (verbose) {
    cat("\nVariation:", paste(varNames[variations + 1], collapse = " + "), "\n")
    catp <- if (post) "Post transformation" else NULL
    catc <- if (final) "Final transformation" else NULL
    cate <- if (extra) "Additional post transformation" else NULL
    cat("Effect:", paste(c("Affine transformation", catp, catc, cate), collapse = " + "), "\n")
  }
  nvariations <- length(variations)
  nfunc <- sample(1:10, size = 1)
  w_i <- stats::runif(nfunc, 0, 1)
  if (user) {
    v_ij <- matrix(1, nrow = nfunc, ncol = nvariations)
  } else {
    v_ij <- matrix(stats::runif(nfunc * nvariations, min = 0, max = 1), nrow = nfunc, ncol = nvariations)
  }
  for (i in 1:nrow(v_ij)) {
    v_ij[i, ] <- v_ij[i, ] / sum(v_ij[i, ])
  } 
  v_params <- c(stats::runif(1, 0, 1), stats::runif(1, -1, 0), stats::runif(1, 1, 10), # blob.high, blob.low, blob.waves
                stats::runif(1, 0, 1), stats::runif(1, 0, 1), stats::runif(1, 0, 1), stats::runif(1, 0, 1), # padj.a, pdj.b, pdj.c, pdj.d
                stats::runif(1, 0, 1), # rings2.val
                stats::runif(1, 1, pi), stats::runif(1, 0, 1), # perspective.angle, perspective.dist
                stats::runif(1, 1, 5), stats::runif(1, 0, 10), # juliaN.power, juliaN.dist
                stats::runif(1, 1, 5), stats::runif(1, 0, 10), # juliaScope.power, juliaScope.dist
                stats::runif(1, 1, pi), stats::runif(1, 1, 5), # radialBlur.angle, v_36
                sample(1:10, size = 1), stats::runif(1, 1, pi), stats::runif(1, 1, 5), # pie.slices, pie.rotation, pie.thickness
                stats::runif(1, 1, 4), 2 * pi / sample(3:10, size = 1), sample(2:10, size = 1), stats::runif(1, 0, 1), # ngon.power, ngon.sides, ngon.corners, ngon.circle
                stats::runif(1, 0, 1), stats::runif(1, 0, 1), # curl.c1, curl.c2
                stats::runif(1, 2, 50), stats::runif(1, 2, 50), # rectangles.x, rectangles.y
                stats::runif(1, 0, 100), # v_41
                stats::runif(1, 0, 10)) # v_44
  df <- iterate_flame(
    iterations = iterations,
    variations = variations,
    point = stats::runif(2, -1, 1),
    w_i = w_i / sum(w_i),
    mat_coef = matrix(stats::runif(nfunc * 6, min = -1, max = 1), nrow = nfunc, ncol = 6),
    blend_variations = blend,
    v_ij = v_ij,
    v_params = v_params,
    transform_p = post,
    p_coef = matrix(stats::runif(nfunc * 6, min = -1, max = 1), nrow = nfunc, ncol = 6),
    transform_f = final,
    f_coef = stats::runif(6, min = -1, max = 1),
    transform_e = extra,
    e_coef = stats::runif(6, min = -1, max = 1)
  )
  df <- df[!is.infinite(df[["x"]]) & !is.infinite(df[["y"]]), ]
  df <- df[!is.na(df[["x"]]) & !is.na(df[["y"]]), ]
  if (nrow(df) == 0) {
    stop("The algorithm did not converge")
  }
  center <- c(stats::median(df[["x"]]), stats::median(df[["y"]]))
  spanx <- diff(quantile(df[["x"]], probs = c(0.1, 0.9))) * (1 / zoom)
  xbins <- seq(center[1] - spanx, center[1] + spanx, length.out = resolution + 1)
  spany <- diff(quantile(df[["y"]], probs = c(0.1, 0.9))) * (1 / zoom)
  ybins <- seq(center[2] - spany, center[2] + spany, length.out = resolution + 1)
  canvas <- color_flame(
    canvas = matrix(0, nrow = resolution + 1, ncol = resolution + 1),
    x = df[["x"]], y = df[["y"]], binsx = xbins, binsy = ybins
  )
  # TODO: stop if canvas is (almost) empty
  full_canvas <- .unraster(canvas, c("x", "y", "z"))
  full_canvas$z[full_canvas$z != 0] <- log(full_canvas$z[full_canvas$z != 0], base = 1.2589)
  full_canvas$z[full_canvas$z == 0] <- NA
  artwork <- ggplot2::ggplot(data = full_canvas, mapping = ggplot2::aes(x = x, y = y, fill = z)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_gradientn(colors = colors, na.value = background)
  artwork <- aRtsy:::theme_canvas(artwork, background = background)
  return(artwork)
}

.getVariationNames <- function() {
  x <- c(
    "Linear",
    "Sine",
    "Spherical",
    "Swirl",
    "Horsehoe",
    "Polar",
    "Handkerchief",
    "Heart",
    "Disc",
    "Spiral",
    "Hyperbolic",
    "Diamond",
    "Ex",
    "Julia",
    "Bent",
    "Waves",
    "Fisheye",
    "Popcorn",
    "Exponential",
    "Power",
    "Cosine",
    "Rings",
    "Fan",
    "Blob",
    "PDJ",
    "Fan2",
    "Rings2",
    "Eyefish",
    "Bubble",
    "Cylinder",
    "Perspective",
    "Noise",
    "JuliaN",
    "JuliaScope",
    "Blur",
    "Gaussian",
    "RadialBlur",
    "Pie",
    "Ngon",
    "Curl",
    "Rectangles",
    "Arch",
    "Tangent",
    "Square",
    "Rays"
  )
  return(x)
}
