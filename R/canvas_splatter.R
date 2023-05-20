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

#' Draw Paint Splatters
#'
#' @description This function draws paint splatters on a canvas using a particle
#'   system.
#'
#' @usage canvas_splatter(
#'   colors,
#'   background = "#fafafa",
#'   iterations = 100,
#'   n = 100,
#'   lwd = 0.1,
#'   resolution = 500
#' )
#'
#' @param colors      a character (vector) specifying the color(s) used for the
#'   artwork.
#' @param background  a character specifying the color used for the background.
#' @param iterations  a positive integer specifying the number of iterations of
#'   the algorithm.
#' @param n           a positive integer specifying the number of particles.
#' @param lwd         expansion factor for the line width.
#' @param resolution  resolution of the artwork in pixels per row/column.
#'   Increasing the resolution increases the quality of the artwork but also
#'   increases the computation time exponentially.
#'
#' @return A \code{ggplot} object containing the artwork.
#'
#' @references \url{https://mattdesl.svbtle.com/generative-art-with-nodejs-and-canvas}
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
#' canvas_splatter(colors = colorPalette("mixer1"))
#' }
#'
#' @export

canvas_splatter <- function(colors, background = "#fafafa", iterations = 100,
                            n = 100, lwd = 0.1, resolution = 500) {
  time <- 0
  seed <- sample.int(.Machine$integer.max, 1)
  heightMap <- .noise(dims = c(resolution, resolution), type = "perlin", limits = c(0, 255))
  particles <- vector("list", n)
  particles <- lapply(particles, .splatter_reset, resolution, c(background, colors))
  canvas <- data.frame(
    x = numeric(iterations * n),
    y = numeric(iterations * n),
    xend = numeric(iterations * n),
    yend = numeric(iterations * n),
    color = character(iterations * n),
    z = numeric(iterations * n),
    width = numeric(iterations * n)
  )
  for (i in seq_len(iterations)) {
    time <- time + 1
    for (j in seq_len(n)) {
      particle <- particles[[j]]
      x <- particle[["xpos"]]
      y <- particle[["ypos"]]
      fx <- min(max(round(x), 0), resolution - 1)
      fy <- min(max(round(y), 0), resolution - 1)
      heightValue <- heightMap[fy, fx] / 255
      pS <- min(max(heightValue, 0.00001), 0.0001)
      angle <- ambient::gen_simplex(x = fx * pS, y = fy * pS, z = particle[["duration"]] + time, seed = seed) * pi * 2
      # Calculate particle speed
      speed <- particle[["speed"]] + min(max(1 - heightValue, 0), 0.1)
      # Update particle velocity
      velocity <- c(particle[["xvel"]], particle[["yvel"]]) + c(cos(angle), sin(angle))
      velocity <- velocity / sqrt(sum(c(particle[["xvel"]]^2, particle[["yvel"]]^2)))
      particle[["xvel"]] <- velocity[1]
      particle[["yvel"]] <- velocity[2]
      # Update particle position
      move <- velocity * speed
      particle[["xpos"]] <- particle[["xpos"]] + move[1]
      particle[["ypos"]] <- particle[["ypos"]] + move[2]
      # Calculate line width
      s2 <- stats::runif(1, 0.0001, 0.05)
      r <- particle[["radius"]] * abs(ambient::gen_simplex(x * s2, y * s2, particle[["duration"]] + time))
      width <- r * min(max(heightValue, 0.01), lwd) * (particle[["time"]] / particle[["duration"]])
      # Record position
      index <- (i - 1) * n + 1 + (j - 1)
      canvas[index, "x"] <- x
      canvas[index, "y"] <- y
      canvas[index, "xend"] <- particle[["xpos"]]
      canvas[index, "yend"] <- particle[["ypos"]]
      canvas[index, "color"] <- particle[["color"]]
      canvas[index, "z"] <- particle[["z"]]
      canvas[index, "width"] <- width
      # Update time
      particle[["time"]] <- particle[["time"]] + 1
      if (particle[["time"]] > particle[["duration"]] || (particle[["xpos"]] <= 0 || particle[["xpos"]] >= resolution) || (particle[["ypos"]] <= 0 || particle[["ypos"]] >= resolution)) {
        particle <- .splatter_reset(particles, resolution, c(background, colors))
      }
      # Update particle
      particles[[j]] <- particle
    }
  }
  canvas <- canvas[((canvas$xend > 0 & canvas$xend < resolution) & (canvas$yend > 0 & canvas$yend < resolution)), ]
  artwork <- ggplot2::ggplot(canvas) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend, group = z), linewidth = canvas[["width"]], color = canvas[["color"]], lineend = "round") +
    ggplot2::scale_x_continuous(limits = c(0, resolution)) +
    ggplot2::scale_y_continuous(limits = c(0, resolution))
  artwork <- aRtsy::theme_canvas(artwork, background = background)
  return(artwork)
}

.splatter_reset <- function(particles, resolution, cols) {
  startArea <- 0.5
  maxRadius <- 10
  particle <- list()
  scale <- resolution / 2
  pos <- .splatter_sphere(c(0, 0), stats::runif(1, 0, scale * startArea))
  if (all(sapply(particles, is.null))) {
    particle[["z"]] <- 1
  } else {
    nums <- unlist(sapply(particles, "[[", "z"))
    particle[["z"]] <- max(nums) + 1
  }
  particle[["xpos"]] <- pos[1] + resolution / 2
  particle[["ypos"]] <- pos[2] + scale
  particle[["radius"]] <- stats::runif(1, 0.01, maxRadius)
  particle[["duration"]] <- stats::runif(1, 1, 500)
  particle[["time"]] <- stats::runif(1, 0, particle[["duration"]])
  particle[["xvel"]] <- stats::runif(1, -1, 1)
  particle[["yvel"]] <- stats::runif(1, -1, 1)
  particle[["speed"]] <- stats::runif(1, 0.5, 1)
  particle[["color"]] <- cols[sample.int(length(cols), 1)]
  return(particle)
}

.splatter_sphere <- function(out, scale) {
  r <- stats::runif(1) * 2.0 * pi
  out <- numeric(2)
  out[1] <- cos(r) * scale
  out[2] <- sin(r) * scale
  return(out)
}
