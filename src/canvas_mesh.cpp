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

Rcpp::IntegerVector int_seq(const int& first,
                            const int& last) {
  Rcpp::IntegerVector y(abs(last - first) + 1);
  if (first < last) {
   std::iota(y.begin(), y.end(), first);
  } else {
   std::iota(y.begin(), y.end(), last);
   std::reverse(y.begin(), y.end());
  }
  return y;
}

void shift_right(Rcpp::IntegerVector& x) {
  int x1 = x[0];
  x.erase(0);
  x.push_back(x1);
}

void shift_right(Rcpp::DoubleVector& x) {
  double x1 = x[0];
  x.erase(0);
  x.push_back(x1);
}

// [[Rcpp::export]]
Rcpp::DataFrame iterate_mesh(arma::mat& canvas,
                             const Rcpp::DoubleVector& points,
                             const Rcpp::DoubleVector& centers,
                             const int& iterations,
                             const int& start,
                             Rcpp::IntegerVector& order,
                             Rcpp::DoubleVector& radii,
                             Rcpp::DoubleVector& increase) {
  const int l = order.length();
  for (int i = 0; i < (iterations + 1); ++i) {
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    Rcpp::DoubleVector newy = start + centers[i] + radii * sin(points);
    const Rcpp::IntegerVector index = int_seq(i * l, i * l + (l- 1));
    for (int j = 0; j < l; ++j) {
      canvas.at(index[j], 0) = newy[j];
      canvas.at(index[j], 1) = order[j];
    }
    shift_right(order);
    radii += increase;
    shift_right(radii);
  }
  return canvas;
}