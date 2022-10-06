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
#' @usage canvas_flame(colors, background = "#000000",
#'                       iterations = 1000000, resolution = 1000,
#'                       variations = NULL, blend = TRUE, 
#'                       post = FALSE, final = FALSE, extra = FALSE,
#'                       verbose = FALSE)
#'
#' @param colors      a string or character vector specifying the color(s) used for the artwork.
#' @param background  a character specifying the color used for the background.
#' @param iterations  a positive integer specifying the number of iterations of the algorithm.
#' @param resolution  resolution of the artwork in pixels per row/column. Increasing the resolution increases the quality of the artwork but also increases the computation time exponentially.
#' @param variations  an integer (vector) specifying the variations to be included in the flame. The default \code{NULL} includes all variations. See the details section for more information about possible variations.
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
#' }
#'
#' @export

canvas_flame <- function(colors, background = "#000000",
                         iterations = 1000000, resolution = 1000,
                         variations = NULL, blend = TRUE,
                         post = FALSE, final = FALSE, extra = FALSE,
                         verbose = FALSE) {
  .checkUserInput(
    resolution = resolution, background = background
  )
  iterations <- iterations + 20
  varNames <- .getVariationNames()
  noVariations <- length(varNames)
  if (is.null(variations)) {
    v <- 0:(noVariations - 1)
    variations <- sample(x = v, size = sample(1:length(v), size = 1))
  } else if (min(variations) < 0 || max(variations) > (noVariations - 1)) {
    stop("'variations' must be between 0 and ", (noVariations - 1)) 
  }
  if (verbose) {
    cat("\nVariation:", paste(varNames[variations + 1], collapse = " + "), "\n")
	catp <- if (post) "Post transformation" else NULL
	catc <- if (final) "Final transformation" else NULL
	cate <- if (extra) "Final transformation" else NULL
	cat("Effect:", paste(c("Affine transformation", catp, catc, cate), collapse = " + "), "\n")
  }
  nvariations <- length(variations)
  nfunc <- sample(1:10, size = 1)
  w_i <- stats::runif(nfunc, 0, 1)
  v_ij <- matrix(stats::runif(nfunc * nvariations, min = 0, max = 1), nrow = nfunc, ncol = nvariations)
  for (i in 1:nrow(v_ij)) {
    v_ij[i, ] <- v_ij[i, ] / sum(v_ij[i, ])
  }
  df <- iterate_flame(
    iterations = iterations,
    variations = variations,
    point = stats::runif(2, -1, 1),
    w_i = w_i / sum(w_i),
    mat_coef = matrix(stats::runif(nfunc * 6, min = -1, max = 1), nrow = nfunc, ncol = 6),
    blend_variations = blend,
    v_ij = v_ij,
    transform_p = post,
    p_coef = matrix(stats::runif(nfunc * 6, min = -1, max = 1), nrow = nfunc, ncol = 6),
    transform_f = final,
	f_coef = stats::runif(6, min = -1, max = 1),
    transform_e = extra,
    e_coef = stats::runif(6, min = -1, max = 1)
  )
  df <- df[!is.infinite(df[["x"]]) & !is.infinite(df[["y"]]), ]
  df <- df[!is.na(df[["x"]]) & !is.na(df[["y"]]), ]
  center <- c(stats::median(df[["x"]]), stats::median(df[["y"]]))
  spanx <- diff(quantile(df[["x"]], probs = c(0.25, 0.75)))
  xbins <- seq(center[1] - spanx, center[1] + spanx, length.out = resolution + 1)
  spany <- diff(quantile(df[["y"]], probs = c(0.25, 0.75)))
  ybins <- seq(center[2] - spany, center[2] + spany, length.out = resolution + 1)
  canvas <- color_flame(
    canvas = matrix(0, nrow = resolution + 1, ncol = resolution + 1),
    x = df[["x"]], y = df[["y"]], binsx = xbins, binsy = ybins
  )
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
  x <- c("Linear", 
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
         "Julia")
  return(x)
}
