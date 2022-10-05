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
#' @description This function draws fractal flames.
#'
#' @usage canvas_flame(colors, background = "#fafafa", iterations = 500000,
#'                       resolution = 1000, span = 2, weights = NULL)
#'
#' @param colors      a string or character vector specifying the color(s) used for the artwork.
#' @param background  a character specifying the color used for the background.
#' @param iterations  a positive integer specifying the number of iterations of the algorithm.
#' @param resolution  resolution of the artwork in pixels per row/column. Increasing the resolution increases the quality of the artwork but also increases the computation time exponentially.
#' @param span        a value indicating the width.
#' @param weights     weights for the different variations.
#'
#' @return A \code{ggplot} object containing the artwork.
#'
#' @references \url{https://en.wikipedia.org/wiki/Fractal_flame}
#'
#' @author Koen Derks, \email{koen-derks@hotmail.com}
#'
#' @keywords artwork canvas
#'
#' @seealso \code{colorPalette}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#'
#' # Simple example
#' canvas_flame(colors = colorPalette("jungle"))
#' }
#'
#' @export

canvas_flame <- function(colors, background = "#fafafa", iterations = 500000,
                         resolution = 1000, span = 2, weights = NULL) {
  .checkUserInput(
    iterations = iterations, resolution = resolution, background = background
  )
  nfunctions <- 9
  if (is.null(weights)) {
    weights <- runif(nfunctions)
  } else {
    if (length(weights) != nfunctions) {
      stop(paste0("'weights' must be a vector of length ", nfunctions))
    }
  }
  df <- iterate_flame(
    iterations = iterations,
    weights = weights + 0.01,
    point = stats::runif(2, -1, 1),
    coef = stats::runif(6, min = -10, max = 10)
  )
  df <- df[!is.infinite(df$x) & !is.infinite(df$y), ]
  df <- df[!is.na(df$x) & !is.na(df$y), ]
  canvas <- matrix(0, nrow = resolution + 1, ncol = resolution + 1)
  center <- c(stats::median(df$x), stats::median(df$y))
  row.bins <- seq(center[2] - span, center[2] + span, length.out = resolution + 1)
  col.bins <- seq(center[1] - span, center[1] + span, length.out = resolution + 1)
  for (i in 1:nrow(df)) {
    indx <- findInterval(df[i, 1], col.bins)
    if (indx == 0 | indx == length(col.bins)) {
      next
    }
    indy <- findInterval(df[i, 2], row.bins)
    if (indy == 0 | indy == length(row.bins)) {
      next
    }
    canvas[indx, indy] <- canvas[indx, indy] + 1
  }
  full_canvas <- .unraster(canvas, c("x", "y", "z"))
  full_canvas$z[full_canvas$z != 0] <- log(full_canvas$z[full_canvas$z != 0])
  full_canvas$z[full_canvas$z == 0] <- NA
  artwork <- ggplot2::ggplot(data = full_canvas, mapping = ggplot2::aes(x = x, y = y, fill = z)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_gradientn(colors = colors, na.value = background)
  artwork <- aRtsy:::theme_canvas(artwork, background = background)
  return(artwork)
}
