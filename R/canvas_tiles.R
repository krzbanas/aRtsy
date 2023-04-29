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

#' Draw Portuguese Tiles
#'
#' @description This function uses a reaction diffusion algorithm in an attempt to draw a Portuguese-styled tiling pattern.
#'
#' @usage canvas_tiles(colors, background = "#ffffff", iterations = 1000,
#'              duplicates = 5, col.line = "#000000", resolution = 100)
#'
#' @param colors      a string or character vector specifying the color(s) used for the artwork.
#' @param background  a character specifying the color of the background.
#' @param duplicates  a non-negative integer specifying how many times the tile should be duplicated.
#' @param col.line    a character specifying the color of the tile borders.
#' @param iterations  a positive integer specifying the number of iterations of the algorithm.
#' @param resolution  resolution of the artwork in pixels per row/column. Increasing the resolution increases the quality of the artwork but also increases the computation time exponentially.
#'
#' @return A \code{ggplot} object containing the artwork.
#'
#' @references \url{https://en.wikipedia.org/wiki/Reaction–diffusion_system} \url{https://www.karlsims.com/rd.html}
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
#' canvas_tiles(colors = colorPalette("bell"), iterations = 5000)
#' }
#'
#' @export

canvas_tiles <- function(colors, background = "#ffffff", iterations = 1000,
                         duplicates = 5, col.line = "#000000", resolution = 100) {
  .checkUserInput(
    resolution = resolution, background = background, iterations = iterations
  )
  rate_a <- 1
  rate_b <- 0.5
  type <- sample.int(4, size = 1)
  feed_rate <- switch(type,
    "1" = 0.0545,
    "2" = 0.055,
    "3" = 0.029,
    "4" = 0.03
  )
  kill_rate <- switch(type,
    "1" = 0.062,
    "2" = 0.062,
    "3" = 0.057,
    "4" = 0.06
  )
  cex.line <- 7 / duplicates
  ones <- matrix(1, resolution, resolution)
  conv_mat <- matrix(c(
    0.05, 0.2, 0.05,
    0.2, -1, 0.2,
    0.05, 0.2, 0.05
  ), nrow = 3)
  binary <- matrix(sample(c(0, 1), size = resolution * resolution, replace = TRUE, prob = c(0.99, 0.01)), ncol = resolution)
  tile <- array(c(ones, binary), dim = c(resolution, resolution, 2))
  tile <- draw_tile(tile, conv_mat, rate_a, rate_b, feed_rate, kill_rate, iterations)[, , 2]
  tile <- cbind(tile, tile[, rev(seq_len(ncol(tile)))])
  tile <- rbind(tile, tile[rev(seq_len(nrow(tile))), ])
  column <- tile
  for (i in seq_len(duplicates)) {
    column <- rbind(column, tile)
  }
  canvas <- column
  for (i in seq_len(duplicates)) {
    canvas <- cbind(canvas, column)
  }
  rownames(canvas) <- colnames(canvas) <- seq_len(nrow(canvas))
  full_canvas <- .unraster(canvas, names = c("x", "y", "z"))
  full_canvas$z <- as.numeric(as.factor(cut(full_canvas$z, breaks = length(colors) + 1)))
  artwork <- ggplot2::ggplot(data = full_canvas, ggplot2::aes(x = x, y = y, fill = z)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_gradientn(colours = c(background, colors)) +
    ggplot2::coord_cartesian(xlim = c(0, (resolution * 2) * (duplicates + 1)), ylim = c(0, (resolution * 2) * (duplicates + 1)))
  if (duplicates > 0 && !is.null(col.line)) {
    lineDataX <- data.frame(x = (resolution * 2) * (0:(duplicates + 1)), xend = (resolution * 2) * (0:(duplicates + 1)), y = rep(0, duplicates + 2), yend = rep((resolution * 2) * (duplicates + 1), duplicates + 2))
    lineDataY <- data.frame(y = (resolution * 2) * (0:(duplicates + 1)), yend = (resolution * 2) * (0:(duplicates + 1)), x = rep(0, duplicates + 2), xend = rep((resolution * 2) * (duplicates + 1), duplicates + 2))
    lineDataX2 <- data.frame(x = resolution * (0:((duplicates * 2) + 1)), xend = resolution * (0:((duplicates * 2) + 1)), y = rep(0, (duplicates * 2) + 2), yend = rep((resolution * 2) * ((duplicates * 2) + 1), (duplicates * 2) + 2))
    lineDataY2 <- data.frame(y = resolution * (0:((duplicates * 2) + 1)), yend = resolution * (0:((duplicates * 2) + 1)), x = rep(0, (duplicates * 2) + 2), xend = rep((resolution * 2) * ((duplicates * 2) + 1), (duplicates * 2) + 2))
    artwork <- artwork + ggplot2::geom_segment(data = lineDataX2, mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend), col = col.line, inherit.aes = FALSE, size = cex.line * 0.1) +
      ggplot2::geom_segment(data = lineDataY2, mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend), col = col.line, inherit.aes = FALSE, size = cex.line * 0.1) +
      ggplot2::geom_segment(data = lineDataX, mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend), col = col.line, inherit.aes = FALSE, size = cex.line * 0.75) +
      ggplot2::geom_segment(data = lineDataY, mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend), col = col.line, inherit.aes = FALSE, size = cex.line * 0.75) +
      ggplot2::geom_segment(data = lineDataX, mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend), col = background, inherit.aes = FALSE, size = cex.line * 0.3) +
      ggplot2::geom_segment(data = lineDataY, mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend), col = background, inherit.aes = FALSE, size = cex.line * 0.3)
  }
  artwork <- aRtsy::theme_canvas(artwork)
  return(artwork)
}
