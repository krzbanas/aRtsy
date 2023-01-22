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
// [[Rcpp::depends(RcppArmadillo)]]

arma::umat shuffle(arma::umat& A) {
  for (int i=0; i < A.n_rows; ++i) {
    A.row(i) = arma::shuffle(A.row(i), 1);
  }
  return A;
}

arma::umat create_palette_rgb(const int& resolution) {
  const int color_count = resolution * resolution;
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
        if (i == color_count) {
          return colors;
        }
      }
    }
  }
  return colors;
}

arma::umat create_palette_gbr(const int& resolution) {
  const int color_count = resolution * resolution;
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
        if (i == color_count) {
          return colors;
        }
      }
    }
  }
  return colors;
}

arma::umat create_palette_brg(const int& resolution) {
  const int color_count = resolution * resolution;
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
        if (i == color_count) {
          return colors;
        }
      }
    }
  }
  return colors;
}

arma::umat create_palette(const int& resolution) {
  arma::umat palette;
  const int pick = floor(R::runif(0, 3));
  switch (pick) {
    case 0:
      palette = create_palette_rgb(resolution);
    case 1:
      palette = create_palette_gbr(resolution);
    case 2:
      palette = create_palette_brg(resolution);
  }
  palette = shuffle(palette);
  return palette;
}

arma::umat sample_palette(const int& resolution,
                          arma::umat& color_mat) {
  const int color_count = resolution * resolution;
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

arma::umat get_palette(const int& resolution,
                       const bool& all_colors,
                       arma::umat& color_mat) {
  arma::umat palette;
  if (all_colors) {
    palette = create_palette(resolution);
  } else {
    palette = sample_palette(resolution, color_mat);
  }
  return palette;
}

double color_difference(const Rcpp::IntegerVector& c1,
                        const Rcpp::IntegerVector& c2) {
  const int r = c1[0] - c2[0];
  const int g = c1[1] - c2[1];
  const int b = c1[2] - c2[2];
  return (r*r + g*g + b*b) >> 1; // Bit-shifting
}

void mark_neighbors(arma::cube& canvas,
                    const Rcpp::IntegerVector& point) {
  const int resolution = canvas.n_rows;
  for (int dy = -1; dy <= 1; ++dy) {
    int ny = point[1] + dy;
    if (ny == -1 || ny == resolution) {
      continue;
    }
    for (int dx = -1; dx <= 1; ++dx) {
      if (dx == 0 && dy == 0) {
        continue;
      }
      int nx = point[0] + dx;
      if (nx == -1 || nx == resolution) {
        continue;
      }
      if (canvas.at(ny, nx, 4) != 1) {
        canvas.at(ny, nx, 3) = 1;
      }
    }
  }
}

Rcpp::IntegerVector min_diff(const arma::cube& canvas,
                             const Rcpp::IntegerVector& color) {
  Rcpp::IntegerVector point(2);
  Rcpp::IntegerVector neighborcolor(3);
  const int resolution = canvas.n_rows;
  // find the pixel with the smallest difference between the current color and all of it's neighbors' colors
  int smallestDifference = 99999999;
  for (int y = 0; y < resolution; ++y) {
    for (int x = 0; x < resolution; ++x) {
      // skip any that arent' available or that are already colored
      if (canvas.at(y, x, 3) != 1 || canvas.at(y, x, 4) == 1) {
        continue;
      }
      // loop through its neighbors
      int smallestDifferenceAmongNeighbors = 99999999;
      for (int dy = -1; dy <= 1; ++dy) {
        if (y + dy == -1 || y + dy == resolution) {
          continue;
        }
        for (int dx = -1; dx <= 1; dx++) {
          if (x == 0 && y == 0) {
            continue;
          }
          if (x + dx == -1 || x + dx == resolution) {
            continue;
          }
          int nx = x + dx;
          int ny = y + dy;
          // skip any neighbors that don't have a color
          if (canvas.at(ny, nx, 4) != 1) {
            continue;
          }
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

Rcpp::IntegerVector min_avg_diff(const arma::cube& canvas,
                                 const Rcpp::IntegerVector& color)  {
  Rcpp::IntegerVector point(2);
  Rcpp::IntegerVector neighborcolor(3);
  int neighborCount, neighborColorDifferenceTotal, averageDifferenceAmongNeighbors, difference;
  int smallestAverageDifference = 99999999;
  const int resolution = canvas.n_rows;
  for (int y = 0; y < resolution; y++) {
    for (int x = 0; x < resolution; x++) {
      // skip any that arent' available or that are already colored
      if (canvas.at(y, x, 3) != 1 || canvas.at(y, x, 4) == 1) {
        continue;
      }
      neighborCount = 0;
      neighborColorDifferenceTotal = 0;
      // loop through its neighbors
      for (int dy = -1; dy <= 1; dy++) {
        if (y + dy == -1 || y + dy == resolution) {
          continue;
        }
        for (int dx = -1; dx <= 1; dx++) {
          if (x == 0 && y == 0) {
            continue;
          }
          if (x + dx == -1 || x + dx == resolution) {
            continue;
          }
          int nx = x + dx;
          int ny = y + dy;
          // skip any neighbors that don't already have a color
          if (canvas.at(ny, nx, 4) != 1) {
            continue;
          }
          ++neighborCount;
          neighborcolor[0] = canvas.at(ny, nx, 0);
          neighborcolor[1] = canvas.at(ny, nx, 1);
          neighborcolor[2] = canvas.at(ny, nx, 2);
          difference = color_difference(neighborcolor, color);
          neighborColorDifferenceTotal += difference;
        }
      }
      averageDifferenceAmongNeighbors = 99999999;
      if (neighborCount > 0) {
        averageDifferenceAmongNeighbors = neighborColorDifferenceTotal / neighborCount;
      }
      if (averageDifferenceAmongNeighbors < smallestAverageDifference || (averageDifferenceAmongNeighbors == smallestAverageDifference && R::runif(0, 1) < 0.5)) {
        smallestAverageDifference = averageDifferenceAmongNeighbors;
        point[0] = x;
        point[1] = y;
      }
    }
  }
  return point;
}

Rcpp::IntegerVector get_point(const arma::cube& canvas,
                              const Rcpp::IntegerVector& color,
                              const int& algorithm) {
  Rcpp::IntegerVector point(2);
  if (algorithm == 0) {
    point = min_diff(canvas, color);
  } else {
    point = min_avg_diff(canvas, color);
  }
  return point;
}

// [[Rcpp::export]]
arma::cube iterate_smoke(arma::cube& canvas,
                         const int& algorithm,
                         const bool& all_colors,
                         arma::umat& color_mat) {
  const int resolution = canvas.n_rows;
  Rcpp::IntegerVector color(3), point(2);
  arma::umat colors = get_palette(resolution, all_colors, color_mat);
  for (int i = 0; i < colors.n_rows; ++i) {
    Rcpp::checkUserInterrupt();
    color = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(colors.row(i)));
    if (i == 0) {
      point[0] = floor(R::runif(0, resolution));
      point[1] = floor(R::runif(0, resolution));
    } else {
      point = get_point(canvas, color, algorithm);
    }
    canvas.at(point[1], point[0], 0) = color[0]; // Red
    canvas.at(point[1], point[0], 1) = color[1]; // Green
    canvas.at(point[1], point[0], 2) = color[2]; // Blue
    canvas.at(point[1], point[0], 3) = 0;        // Available
    canvas.at(point[1], point[0], 4) = 1;        // Filled
    mark_neighbors(canvas, point);
  }
  return canvas;
}
