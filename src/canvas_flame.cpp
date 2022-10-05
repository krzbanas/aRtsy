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
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <iterator>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp ;

// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::DoubleVector double_seq(int first, int last) {
  Rcpp::DoubleVector y(abs(last - first) + 1);
  if (first < last) {
   std::iota(y.begin(), y.end(), first);
  } else {
   std::iota(y.begin(), y.end(), last);
   std::reverse(y.begin(), y.end());
  }
  return y;
}

Rcpp::DoubleVector affine(Rcpp::DoubleVector p, 
                         double a, 
                         double b, 
                         double c, 
                         double d, 
                         double e, 
                         double f) {
  Rcpp::DoubleVector x(2);
  x[0] = a * p[0] + b * p[1] + c;
  x[1] = d * p[0] + e * p[1] + f;
  return x;
}

Rcpp::DoubleVector variation(Rcpp::DoubleVector p, int i) {;
  Rcpp::DoubleVector x(2);
  double r = sqrt(pow(p[0], 2) + pow(p[1], 2));
  double theta = atan(p[0] / p[1]);
  if (i == 0) { // Linear
    x[0] = p[0];
    x[1] = p[1];
  } else if (i == 1) { // Sine
    x[0] = sin(p[0]);
    x[1] = sin(p[1]);
  } else if (i == 2) { // Sperical
    x[0] = p[0] / pow(r, 2);
    x[1] = p[1] / pow(r, 2);
  } else if (i == 3) { // Swirl
    x[0] = p[0] * sin(pow(r, 2)) - p[1] * cos(pow(r, 2));
    x[1] = p[0] * cos(pow(r, 2)) + p[1] * sin(pow(r, 2));
  } else if (i == 4) { // Horsehoe
    x[0] = 1 / r * (p[0] - p[1]) * (p[0] + p[1]);
    x[1] = 1 / r * 2 * p[0] * p[1];
  } else if (i == 5) { // Polar
    x[0] = theta / M_PI;
    x[1] = r - 1;
  } else if (i == 6) { // Handkerchief
    x[0] = r * sin(theta + r);
    x[1] = r * cos(theta - r);
  } else if (i == 7) { // Heart
    x[0] = r * sin(theta * r);
    x[1] = r * (-cos(theta * r));
  } else if (i == 8) { // Disc
    x[0] = theta / M_PI * sin(M_PI * r);
    x[1] = theta / M_PI * cos(M_PI * r);
  }
  return x;
}

Rcpp::DoubleVector posttransform(Rcpp::DoubleVector p, 
                         double alpha, 
                         double beta, 
                         double gamma, 
                         double delta, 
                         double epsilon, 
                         double zeta) {
  Rcpp::DoubleVector x(2);
  x[0] = alpha * p[0] + beta * p[1] + gamma;
  x[1] = delta * p[0] + epsilon * p[1] + zeta;
  return x;
}

// [[Rcpp::export]]
Rcpp::DataFrame iterate_flame(int iterations,
                              Rcpp::DoubleVector point,
                              Rcpp::DoubleVector w_i,
                              arma::mat v_ij,
                              arma::mat mat_coef,
                              arma::mat p_coef,
                              Rcpp::DoubleVector f_coef,
                              Rcpp::DoubleVector p2_coef) {
  int nvariations = v_ij.n_cols;
  int npoints = (iterations - 20);
  Rcpp::DoubleVector x(npoints);
  Rcpp::DoubleVector y(npoints);
  Rcpp::DoubleVector p(2);
  // TODO: Add color mixing
  for (int iter = 0; iter < iterations; iter++) {
    Rcpp::checkUserInterrupt();
    Rcpp::NumericVector i = RcppArmadillo::sample(double_seq(0, w_i.length() - 1), 1, false, w_i);
    Rcpp::DoubleVector newpoint(2);
    p = affine(point, mat_coef(i[0], 0), mat_coef(i[0], 1), mat_coef(i[0], 2), mat_coef(i[0], 3), mat_coef(i[0], 4), mat_coef(i[0], 5));
    for (int j = 0; j < nvariations; j++) {
      newpoint += v_ij(i[0], j) * variation(p, j);
    }
    point = newpoint;
    point = posttransform(point, p_coef(i[0], 0), p_coef(i[0], 1), p_coef(i[0], 2), p_coef(i[0], 3), p_coef(i[0], 4), p_coef(i[0], 5));
    point = affine(point, f_coef[0], f_coef[1], f_coef[2], f_coef[3], f_coef[4], f_coef[5]);
    point = posttransform(point, p2_coef[0], p2_coef[1], p2_coef[2], p2_coef[3], p2_coef[4], p2_coef[5]);
    if (iter > 19) {
      x[iter - 20] = point[0];
      y[iter - 20] = point[1];
    }
  }
  Rcpp::DataFrame flame = Rcpp::DataFrame::create(Rcpp::Named("x") = x,
                                                  Rcpp::Named("y") = y);
  return flame;
}
