// Copyright (C) 2021-2022 Koen Derks

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>

arma::umat shuffle(arma::umat& A) {
  for (int i=0; i < A.n_rows; ++i) {
    A.row(i) = arma::shuffle( A.row(i), 1 );
  }
  return A;
}

double color_difference(const Rcpp::IntegerVector& c1,
                        const Rcpp::IntegerVector& c2) {
  int r = c1[0] - c2[0];
  int g = c1[1] - c2[1];
  int b = c1[2] - c2[2];
  double diff = sqrt(pow(r, 2) + pow(g, 2) + pow(b, 2));
  return diff;
}

void mark_neighbors(const Rcpp::IntegerVector& point,
                    arma::umat& available,
                    const arma::umat& colored,
                    const int& resolution) {
  for (int dy = -1; dy <= 1; ++dy) {
    int ny = point[1] + dy;
    if (ny == -1 || ny == resolution)
      continue;
    for (int dx = -1; dx <= 1; ++dx) {
      if (dx == 0 && dy == 0)
        continue;
      int nx = point[0] + dx;
      if (nx == -1 || nx == resolution)
        continue;
      if (colored.at(ny, nx) != 1) {
        available.at(ny, nx) = 1;
      }
    }
  }
}

Rcpp::IntegerVector min_diff(const arma::cube& canvas,
                             const arma::umat& available,
                             const arma::umat& colored,
                             const Rcpp::IntegerVector& color,
                             const int& resolution) {
  Rcpp::IntegerVector point(2);
  Rcpp::IntegerVector neighborcolor(3);
  // find the pixel with the smallest difference between the current color and all of it's neighbors' colors
  int smallestDifference = 99999999;
  for (int y = 0; y < resolution; ++y) {
    for (int x = 0; x < resolution; ++x) {
      // skip any that arent' available or that are already colored
      if (available.at(y, x) != 1 || colored.at(y, x) == 1)
        continue;
      // loop through its neighbors
      int smallestDifferenceAmongNeighbors = 99999999;
      for (int dy = -1; dy <= 1; ++dy) {
        if (y + dy == -1 || y + dy == resolution)
          continue;
        for (int dx = -1; dx <= 1; dx++) {
          if (x == 0 && y == 0)
            continue;
          if (x + dx == -1 || x + dx == resolution)
            continue;
          int nx = x + dx;
          int ny = y + dy;
          // skip any neighbors that don't have a color
          if (colored.at(ny, nx) != 1)
            continue;
          neighborcolor[0] = canvas.at(ny, nx, 0);
          neighborcolor[1] = canvas.at(ny, nx, 1);
          neighborcolor[2] = canvas.at(ny, nx, 2);
          int difference = color_difference(neighborcolor, color);
          if (difference < smallestDifferenceAmongNeighbors) {
            smallestDifferenceAmongNeighbors = difference;
          }
        }
      }
      if (smallestDifferenceAmongNeighbors < smallestDifference || (smallestDifferenceAmongNeighbors == smallestDifference && R::runif(0, 1) < 0.5)) {
        smallestDifference = smallestDifferenceAmongNeighbors;
        point[0] = x;
        point[1] = y;
      }
    }
  }
  return point;
}

Rcpp::IntegerVector avg_diff(const arma::cube& canvas,
                             const arma::umat& available,
                             const arma::umat& colored,
                             const Rcpp::IntegerVector& color,
                             const int& resolution)  {
  Rcpp::IntegerVector point(2);
  Rcpp::IntegerVector neighborcolor(3);
  int neighborCount, neighborColorDifferenceTotal, difference;
  int smallestAverageDifference = 99999999;
  for (int y = 0; y < resolution; y++) {
    for (int x = 0; x < resolution; x++) {
      // skip any that arent' available or that are already colored
      if (available.at(y, x) != 1 || colored.at(y, x) == 1)
        continue;
      neighborCount = 0;
      neighborColorDifferenceTotal = 0;
      // loop through its neighbors
      for (int dy = -1; dy <= 1; dy++) {
        if (y + dy == -1 || y + dy == resolution)
          continue;
        for (int dx = -1; dx <= 1; dx++) {
          if (x == 0 && y == 0)
            continue;
          if (x + dx == -1 || x + dx == resolution)
            continue;
          int nx = x + dx;
          int ny = y + dy;
          // skip any neighbors that don't already have a color
          if (colored.at(ny, nx) != 1)
            continue;
          ++neighborCount;
          neighborcolor[0] = canvas.at(ny, nx, 0);
          neighborcolor[1] = canvas.at(ny, nx, 1);
          neighborcolor[2] = canvas.at(ny, nx, 2);
          difference = color_difference(neighborcolor, color);
          neighborColorDifferenceTotal += difference;
        }
      }
      int averageDifferenceAmongNeighbors = neighborColorDifferenceTotal / neighborCount;
      if (averageDifferenceAmongNeighbors < smallestAverageDifference || (averageDifferenceAmongNeighbors == smallestAverageDifference && R::runif(0, 1) < 0.5)) {
        smallestAverageDifference = averageDifferenceAmongNeighbors;
        point[0] = x;
        point[1] = y;
      }
    }
  }
  return point;
}

Rcpp::IntegerVector get_point(arma::cube& canvas,
                              arma::umat& available,
                              arma::umat& colored,
                              Rcpp::IntegerVector& color,
                              int resolution,
                              int algorithm) {
  switch (algorithm) {
    case 0:
      return min_diff(canvas, available, colored, color, resolution);
    case 1:
      return avg_diff(canvas, available, colored, color, resolution);
  }
}

arma::umat createPaletteRGB(const int& resolution) {
  const int color_count = pow(resolution, 2);
  const int numcolors = ceil(pow(color_count, 1.0/3));
  arma::umat colors(color_count, 3);
  int i = 0;
  for (int b = 0; b < numcolors; b++) {
    for (int g = 0; g < numcolors; g++) {
      for (int r = 0; r < numcolors; r++) {
        colors.at(i, 0) = r * 255 / (numcolors - 1);
        colors.at(i, 1) = g * 255 / (numcolors - 1);
        colors.at(i, 2) = b * 255 / (numcolors - 1);
        ++i;
        if (i == color_count)
          return colors;
      }
    }
  }
}

arma::umat createPaletteGBR(const int& resolution) {
  const int color_count = pow(resolution, 2);
  const int numcolors = ceil(pow(color_count, 1.0/3));
  arma::umat colors(color_count, 3);
  int i = 0;
  for (int r = 0; r < numcolors; r++) {
    for (int b = 0; b < numcolors; b++) {
      for (int g = 0; g < numcolors; g++) {
        colors.at(i, 0) = r * 255 / (numcolors - 1);
        colors.at(i, 1) = g * 255 / (numcolors - 1);
        colors.at(i, 2) = b * 255 / (numcolors - 1);
        ++i;
        if (i == color_count)
          return colors;
      }
    }
  }
}

arma::umat createPaletteBRG(const int& resolution) {
  const int color_count = pow(resolution, 2);
  const int numcolors = ceil(pow(color_count, 1.0/3));
  arma::umat colors(color_count, 3);
  int i = 0;
  for (int g = 0; g < numcolors; g++) {
    for (int r = 0; r < numcolors; r++) {
      for (int b = 0; b < numcolors; b++) {
        colors.at(i, 0) = r * 255 / (numcolors - 1);
        colors.at(i, 1) = g * 255 / (numcolors - 1);
        colors.at(i, 2) = b * 255 / (numcolors - 1);
        ++i;
        if (i == color_count)
          return colors;
      }
    }
  }
}

arma::umat createPalette(const int& resolution) {
  arma::umat palette;
  int pick = floor(R::runif(0, 3));
  switch (pick) {
  case 0:
    palette = createPaletteRGB(resolution);
  case 1:
    palette = createPaletteGBR(resolution);
  case 2:
    palette = createPaletteBRG(resolution);
  }
  palette = shuffle(palette);
  return palette;
}

arma::umat samplePalette(const int resolution,
                         arma::umat color_mat) {
  const int color_count = pow(resolution, 2);
  const int numcolors = ceil(pow(color_count, 1.0/3));
  arma::umat palette(color_count, 3);
  int pick;
  for (int i = 0; i < color_count; ++i) {
    pick = floor(R::runif(0, color_mat.n_rows));
    palette.at(i, 0) = color_mat.at(pick, 0);
	palette.at(i, 1) = color_mat.at(pick, 1);
	palette.at(i, 2) = color_mat.at(pick, 2);
  }
  return palette;
}

// [[Rcpp::export]]
arma::cube iterate_smoke(arma::cube& canvas,
                         const int& algorithm,
						 bool all_colors,
						 arma::umat color_mat) {
  bool first = true;
  const int resolution = canvas.n_rows;
  arma::umat colors;
  if (all_colors) {
    colors = createPalette(resolution);
  } else {
	colors = samplePalette(resolution, color_mat);
  }
  arma::umat available(resolution, resolution, arma::fill::zeros);
  arma::umat colored(resolution, resolution, arma::fill::zeros);
  Rcpp::IntegerVector color(3);
  Rcpp::IntegerVector point(2);
  for (int i = 0; i < colors.n_rows; ++i) {
    Rcpp::checkUserInterrupt();
    // Retreive the current color
    color[0] = colors.at(i, 0); // Red
    color[1] = colors.at(i, 1); // Green
    color[2] = colors.at(i, 2); // Blue
    if (first) {
      // Use the starting point
      point[0] = floor(R::runif(0, resolution));
      point[1] = floor(R::runif(0, resolution));
      first = false;
    } else {
      // Use the point with the lowest neighbor dissimilarity
      point = get_point(canvas, available, colored, color, resolution, algorithm);
    }
    // Set available and colored
    available.at(point[1], point[0]) = 0;
    colored.at(point[1], point[0]) = 1;
    // Color the pixel
    canvas.at(point[1], point[0], 0) = color[0];
    canvas.at(point[1], point[0], 1) = color[1];
    canvas.at(point[1], point[0], 2) = color[2];
    // Mark neighbors as available
    mark_neighbors(point, available, colored, resolution);
  }
  return canvas;
}
