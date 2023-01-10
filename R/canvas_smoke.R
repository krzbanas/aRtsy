# Copyright (C) 2021-2023 Koen Derks

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

#' Draw Rainbow Smoke
#'
#' @description This function creates rainbow smoke.
#'
#' @usage canvas_smoke(colors, resolution = 100, distance = c("minimum", "average"))
#'
#' @param colors      a string or character vector specifying the color(s) used for the artwork.
#' @param resolution  resolution of the artwork in pixels per row/column. Increasing the resolution increases the quality of the artwork but also increases the computation time exponentially.
#' @param distance    an character specifying whether to take the neighbor with the minimum distance or the smallest average distance.
#'
#' @return A \code{ggplot} object containing the artwork.
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
#' canvas_smoke(colors = colorPalette("random", 1024), resolution = 250)
#' }
#'
#' @export
canvas_smoke <- function(colors, resolution = 100, distance = c("minimum", "average")) {
  .checkUserInput(resolution = resolution)
  distance <- match.arg(distance)
  all_colors <- length(colors) == 1 && colors[1] == "all"
  color_mat <- .getColorMat(colors, all_colors)
  algorithm <- switch(distance,
    "minimum" = 0,
    "average" = 1
  )
  canvas <- array(-1, c(resolution, resolution, 3))
  canvas <- iterate_smoke(canvas, algorithm, all_colors, color_mat)
  full_canvas <- as.data.frame(expand.grid(x = 1:resolution, y = 1:resolution))
  full_canvas[["col"]] <- grDevices::rgb(red = canvas[, , 1], green = canvas[, , 2], blue = canvas[, , 3], maxColorValue = 255)
  artwork <- ggplot2::ggplot(data = full_canvas, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_raster(interpolate = TRUE, fill = full_canvas[["col"]])
  artwork <- theme_canvas(artwork)
  return(artwork)
}

.getColorMat <- function(colors, all_colors) {
  if (all_colors) {
    return(matrix(NA, 1, 1))
  } else {
    return(t(as.matrix(col2rgb(colors))))
  }
}
