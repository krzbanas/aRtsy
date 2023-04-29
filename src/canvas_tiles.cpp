// Copyright (C) 2021-2023 Koen Derks

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

int next_index(const int& resolution,
               const int& index) {
  if (index < 0)
    return index + resolution;
  if(index >= resolution)
    return index - resolution;
  return index;
}

arma::cube iterate_tile(const arma::cube& canvas,
                        const arma::mat& conv,
                        const double& diffa,
                        const double& diffb,
                        const double& feedrate,
                        const double& killrate){
  const int resolution = canvas.n_rows;
  arma::cube new_canvas(resolution, resolution, 2);
  for(int x = 0; x < resolution; ++x) {
    Rcpp::checkUserInterrupt();
    for(int y = 0; y < resolution; ++y){
      double asum = 0.0;
      double bsum = 0.0;
      for(int i = -1; i <= 1; ++i){
        for(int j = -1; j <= 1; ++j){
          const int x1 = next_index(resolution, x - i);
          const int y1 = next_index(resolution, y - j);
          asum += conv.at(i + 1, j + 1) * canvas.at(x1, y1, 0);
          bsum += conv.at(i + 1, j + 1) * canvas.at(x1, y1, 1);
        }
      }
      const double xpow = pow(canvas.at(x, y, 1), 2);
      const double kplusf = killrate + feedrate;
      new_canvas.at(x, y, 0) = canvas.at(x, y, 0) + diffa * asum - canvas.at(x, y, 0) * xpow + feedrate * (1 - canvas.at(x, y, 0));
      new_canvas.at(x, y, 1) = canvas.at(x, y, 1) + diffb * bsum + canvas.at(x, y, 0) * xpow - kplusf * canvas.at(x, y, 1);
    }
  }
  return new_canvas;
}

// [[Rcpp::export]]
arma::cube draw_tile(arma::cube& canvas,
                     const arma::mat& conv,
                     const double& diffa,
                     const double& diffb,
                     const double& feedrate,
                     const double& killrate,
                     const int& iterations){
  for(int i = 0; i < iterations; ++i){
    canvas = iterate_tile(canvas, conv, diffa, diffb, feedrate, killrate); 
  }
  return canvas;
}
