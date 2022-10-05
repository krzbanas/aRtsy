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

Rcpp::IntegerVector int_seq(int first, int last) {
  Rcpp::IntegerVector y(abs(last - first) + 1);
  if (first < last) {
   std::iota(y.begin(), y.end(), first);
  } else {
   std::iota(y.begin(), y.end(), last);
   std::reverse(y.begin(), y.end());
  }
  return y;
}

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

// Linear function is implied as variation0

// Sine
Rcpp::DoubleVector variation1(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  x[0] = sin(p[0]);
  x[1] = sin(p[1]);
  return x;
}

// Spherical
Rcpp::DoubleVector variation2(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  double r = pow(p[0], 2) + pow(p[1], 2);
  x[0] = p[0] / r;
  x[1] = p[1] / r;
  return x;
}

// Swirl
Rcpp::DoubleVector variation3(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  double r = pow(p[0], 2) + pow(p[1], 2);
  x[0] = p[0] * sin(r) - p[1] * cos(r);
  x[1] = p[0] * cos(r) + p[1] * sin(r);
  return x;
}

// Horsehoe
Rcpp::DoubleVector variation4(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  double r = sqrt(pow(p[0], 2) + pow(p[1], 2));
  x[0] = 1 / r * (p[0] - p[1]) * (p[0] + p[1]);
  x[1] = 1 / r * 2 * p[0] * p[1];
  return x;
}

// Polar
Rcpp::DoubleVector variation5(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  double theta = atan(p[0] / p[1]);
  double r = sqrt(pow(p[0], 2) + pow(p[1], 2));
  x[0] = theta / M_PI;
  x[1] = r - 1;
  return x;
}

// Handkerchief
Rcpp::DoubleVector variation6(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  double theta = atan(p[0] / p[1]);
  double r = sqrt(pow(p[0], 2) + pow(p[1], 2));
  x[0] = r * sin(theta + r);
  x[1] = r * cos(theta - r);
  return x;
}

// Heart
Rcpp::DoubleVector variation7(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  double theta = atan(p[0] / p[1]);
  double r = sqrt(pow(p[0], 2) + pow(p[1], 2));
  x[0] = r * sin(theta * r);
  x[1] = r * (-cos(theta * r));
  return x;
}

// Disc
Rcpp::DoubleVector variation8(Rcpp::DoubleVector p) {
  Rcpp::DoubleVector x(2);
  double theta = atan(p[0] / p[1]);
  double r = sqrt(pow(p[0], 2) + pow(p[1], 2));
  x[0] = theta / M_PI * sin(M_PI * r);
  x[1] = theta / M_PI * cos(M_PI * r);
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
                              Rcpp::DoubleVector weights,
                              Rcpp::DoubleVector point,
                              Rcpp::DoubleVector coef) {
  int n = (iterations - 20);
  Rcpp::DoubleVector x(n);
  Rcpp::DoubleVector y(n);
  for (int i = 0; i < iterations; i++) {
    Rcpp::checkUserInterrupt();
	point = affine(point, coef[0], coef[1], coef[2], coef[3], coef[4], coef[5]);
	Rcpp::NumericVector choose = RcppArmadillo::sample(double_seq(0, weights.length() - 1), 1, false, weights);
	// if choose = 0, linear variation
	if (choose[0] == 1) {
	  point = variation1(point);
	} else if (choose[0] == 2) {
		point = variation2(point);
	} else if (choose[0] == 3) {
		point = variation3(point);
	} else if (choose[0] == 4) {
		point = variation4(point);
	} else if (choose[0] == 5) {
		point = variation5(point);
	} else if (choose[0] == 6) {
		point = variation6(point);
	} else if (choose[0] == 7) {
		point = variation7(point);
	} else if (choose[0] == 8) {
		point = variation8(point);
	}
    //point = posttransform(point, p_vec[0], p_vec[1], p_vec[2], p_vec[3], p_vec[4], p_vec[5]);
	//   point = affine(point, f_vec[0], f_vec[1], f_vec[2], f_vec[3], f_vec[4], f_vec[5]);
	// }
    if (i > 19) {
      x[i - 20] = point[0];
      y[i - 20] = point[1];
    }
  }
  Rcpp::DataFrame flame = Rcpp::DataFrame::create(Rcpp::Named("x") = x,
                                                  Rcpp::Named("y") = y);
  return flame;
}
