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

// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::DoubleVector variation(double x,
                             double y,
                             int i,
                             double a,
                             double b,
                             double c,
                             double d,
                             double e,
                             double f,
                             Rcpp::DoubleVector pparams) {;
  Rcpp::DoubleVector p(2);
  double r = sqrt(pow(x, 2) + pow(y, 2));
  double theta = atan(x / y);
  double phi = atan(y / x);
  if (i == 0) { // Linear
    p[0] = x;
    p[1] = y;
  } else if (i == 1) { // Sine
    p[0] = sin(x);
    p[1] = sin(y);
  } else if (i == 2) { // Sperical
    p[0] = x / pow(r, 2);
    p[1] = y / pow(r, 2);
  } else if (i == 3) { // Swirl
    p[0] = x * sin(pow(r, 2)) - y * cos(pow(r, 2));
    p[1] = x * cos(pow(r, 2)) + y * sin(pow(r, 2));
  } else if (i == 4) { // Horsehoe
    p[0] = (1 / r) * ((x - y) * (x + y));
    p[1] = (1 / r) * (2 * x * y);
  } else if (i == 5) { // Polar
    p[0] = theta / M_PI;
    p[1] = r - 1;
  } else if (i == 6) { // Handkerchief
    p[0] = r * sin(theta + r);
    p[1] = r * cos(theta - r);
  } else if (i == 7) { // Heart
    p[0] = r * sin(theta * r);
    p[1] = r * (-cos(theta * r));
  } else if (i == 8) { // Disc
    p[0] = theta / M_PI * sin(M_PI * r);
    p[1] = theta / M_PI * cos(M_PI * r);
  } else if (i == 9) { // Spiral
    p[0] = 1 / r * (cos(theta) + sin(r));
    p[1] = 1 / r * (sin(theta) + cos(r));
  } else if (i == 10) { // Hyperbolic
    p[0] = sin(theta) / r;
    p[1] = r * cos(theta);
  } else if (i == 11) { // Diamond
    p[0] = sin(theta) * cos(r);
    p[1] = cos(theta) * sin(r);
  } else if (i == 12) { // Ex
    double p0 = sin(theta + r);
    double p1 = cos(theta - r);
    p[0] = r * (pow(p0, 3) + pow(p1, 3));
    p[1] = r * (pow(p0, 3) - pow(p1, 3));
  } else if (i == 13) { // Julia
    double s = R::runif(0, 1);
    double Omega;
    if (s < 0.5) {
      Omega = 0;
    } else {
      Omega = M_PI;
    }
    p[0] = sqrt(r) * cos(theta / 2 + Omega);
    p[1] = sqrt(r) * sin(theta / 2 + Omega);
  } else if (i == 14) { // Bent
    if ((x >= 0) && (y >= 0)) {
      p[0] = x;
      p[1] = y;
    } else if ((x < 0) && (y >= 0)) {
      p[0] = 2 * x;
      p[1] = y;
    } else if ((x >= 0) && (y < 0)) {
      p[0] = x;
      p[1] = y / 2;
    } else if ((x < 0) && (y < 0)) {
      p[0] = 2 * x;
      p[1] = y / 2;
    }
  } else if (i == 15) { // Waves
    p[0] = x + b * sin(y / pow(c, 2));
    p[1] = y + e * sin(x / pow(f, 2));
  } else if (i == 16) { // Fisheye
    p[0] = (2 / (r + 1)) * y;
    p[1] = (2 / (r + 1)) * x;
  } else if (i == 17) { // Popcorn
    p[0] = x + c * sin(tan(3 * y));
    p[1] = y + f * sin(tan(3 * x));
  } else if (i == 18) { // Exponential
    p[0] = exp(x - 1) * cos(M_PI * y);
    p[1] = exp(x - 1) * sin(M_PI * y);
  } else if (i == 19) { // Power
    p[0] = pow(r, sin(theta)) * cos(theta);
    p[1] = pow(r, sin(theta)) * sin(theta);
  } else if (i == 20) { // Cosine
    p[0] = cos(M_PI * x) * cosh(y);
    p[1] = -(sin(M_PI * x) * sinh(y));
  } else if (i == 21) { // Rings
    double first = fmod(r + pow(c, 2), 2 * pow(c, 2)) - pow(c, 2) + r * (1 - pow(c, 2));
    p[0] = first * cos(theta);
    p[1] = first * sin(theta);
  } else if (i == 22) { // Fan
    double t = M_PI * pow(c, 2);
    if (fmod(theta + f, t) > (t / 2)) {
      p[0] = r * cos(theta - (t/2));
      p[1] = r * sin(theta - (t/2));
    } else {
      p[0] = r * cos(theta + (t/2));;
      p[1] = r * sin(theta + (t/2));
    }
  } else if (i == 23) { // Blob
    double first = r * (pparams[1] + ((pparams[0] - pparams[1])/ 2) * (sin(pparams[2] * theta) + 1));
    p[0] = first * cos(theta);
    p[1] = first * sin(theta);
  } else if (i == 24) { // PDJ
    p[0] = sin(pparams[3] * y) - cos(pparams[4] * x);
    p[1] = sin(pparams[5] * x) - cos(pparams[6] * y);
  } else if (i == 25) { // Fan2
    Rcpp::DoubleVector fanresult = variation(x, y, 22, a, b, c, d, e, f, pparams);
    double p1 = M_PI * pow(fanresult[0], 2);
    double p2 = fanresult[1];
    double t = theta + p2 - p1 * trunc((2 * theta * p2) / p1);
    if (t > (p1 / 2)) {
      p[0] = r * sin(theta - (p1 / 2));
      p[1] = r * cos(theta - (p1 / 2));
    } else {
      p[0] = r * sin(theta + (p1 / 2));
      p[1] = r * cos(theta + (p1 / 2));
    }
  } else if (i == 26) { // Rings2
    double pp = pow(pparams[7], 2);
    double t = r - 2 * pp * trunc((r + pp) / (2 * pp)) + r * (1 - pp);
    p[0] = t * sin(theta);
    p[1] = t * cos(theta);
  } else if (i == 27) { // Eyefish
    p[0] = (2 / (r + 1)) * x;
    p[1] = (2 / (r + 1)) * y;
  } else if (i == 28) { // Bubble
    p[0] = (4 / (pow(r, 2) + 4)) * x;
    p[1] = (4 / (pow(r, 2) + 4)) * y;	
  } else if (i == 29) { // Cylinder
    p[0] = sin(x);
    p[1] = y;
  } else if (i == 30) { // Perspective
    double first = (pparams[9] / (pparams[9] - y * sin(pparams[8])));
    p[0] = first * x;
    p[1] = first * (y * cos(pparams[8]));
  } else if (i == 31) { // Noise
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    p[0] = Psi1 * (x * cos(2 * M_PI * Psi2));
    p[1] = Psi1 * (y * sin(2 * M_PI * Psi2));
  } else if (i == 32) { // JuliaN
    double Psi = R::runif(0, 1);
    double p3 = trunc(fabs(pparams[11]) * Psi);
    double t = (phi + 2 * M_PI * p3) / pparams[11];
    p[0] = pow(r, pparams[11] / pparams[10]) * cos(t);
    p[1] = pow(r, pparams[11] / pparams[10]) * sin(t);
  } else if (i == 33) { // JuliaScope
    double Psi = R::runif(0, 1);
    int Lambda = floor(R::runif(0, 2));
    double p3 = trunc(fabs(pparams[12] * Psi));
    double t = (Lambda * phi + 2 * M_PI * p3) / pparams[12];
    p[0] = pow(r, pparams[13] / pparams[12]) * cos(t);
    p[1] = pow(r, pparams[13] / pparams[12]) * sin(t);
  } else if (i == 34) { // Blur
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    p[0] = Psi1 * cos(2 * M_PI * Psi2);
    p[1] = Psi1 * sin(2 * M_PI * Psi2);
  } else if (i == 35) { // Gaussian
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double Psi3 = R::runif(0, 1);
    double Psi4 = R::runif(0, 1);
    double Psi5 = R::runif(0, 1);
    p[0] = (Psi1 + Psi2 + Psi3 + Psi4 - 2) * cos(2 * M_PI * Psi5);
    p[1] = (Psi1 + Psi2 + Psi3 + Psi4 - 2) * sin(2 * M_PI * Psi5);
  } else if (i == 36) { // RadialBlur
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double Psi3 = R::runif(0, 1);
    double Psi4 = R::runif(0, 1);
    double p1 = pparams[14] * (M_PI / 2);
    double t1 = pparams[15] * (Psi1 + Psi2 + Psi3 + Psi4 - 2);
    double t2 = phi + t1 * sin(p1);
    double t3 = t1 * cos(p1) - 1;
    p[0] = (1 / pparams[15]) * (r * cos(t2) + t3 * x);
    p[1] = (1 / pparams[15]) * (r * sin(t2) + t3 * y);
  } else if (i == 37) { // Pie
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double Psi3 = R::runif(0, 1);
    double t1 = trunc(Psi1 * pparams[16] + 0.5);
    double t2 = pparams[17] + ((2 * M_PI) / pparams[16]) * (t1 + Psi2 * pparams[17]);
    p[0] = Psi3 * cos(t2);
    p[1] = Psi3 * sin(t2);
  } else if (i == 38) { // Ngon
    double t3 = phi - pparams[19] * floor(phi / pparams[19]);
    double t4;
    if (t3 > (pparams[19] / 2)) {
      t4 = t3;
    } else {
      t4 = t3 - pparams[19];
    }
    double k = (pparams[20] * ((1 / cos(t4) - 1) + pparams[21]) + pparams[22]) / pow(r, pparams[18]);
    p[0] = k * x;
    p[1] = k * y;
  } else if (i == 39) { // Curl
    double t1 = 1 + pparams[23] * x + pparams[24] * (pow(x, 2) - pow(y, 2));
    double t2 = pparams[23] * y + 2 * pparams[24] * x * y;
    double first = 1 / (pow(t1, 2) + pow(t2, 2));
    p[0] = first * (x * t1 + y * t2);
    p[1] = first * (y * t1 - x * t2);
  } else if (i == 40) { // Rectangles
    p[0] = (2 * floor(x / pparams[25]) + 1) * pparams[25] - x;
    p[1] = (2 * floor(y / pparams[26]) + 1) * pparams[26] - y;
  } else if (i == 41) { // Arch
    double Psi = R::runif(0, 1);
    p[0] = sin(Psi * M_PI * pparams[27]);
    p[1] = pow(sin(Psi * M_PI * pparams[27]), 2) / cos(Psi * M_PI * pparams[27]);
  } else if (i == 42) { // Tangent
    p[0] = sin(x) / cos(y);
    p[1] = tan(y);
  } else if (i == 43) { // Square
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    p[0] = Psi1 - 0.5;
    p[1] = Psi2 - 0.5;
  } else if (i == 44) { // Rays
    double Psi = R::runif(0, 1);
    double first = (pparams[28] * tan(Psi * M_PI * pparams[28])) / pow(r, 2);
    p[0] = first * cos(x);
    p[1] = first * sin(x);
  } else if (i == 45) { // Rays
    double Psi = R::runif(0, 1);
    p[0] = x * (cos(Psi * r * pparams[29]) + sin(Psi * r * pparams[29]));
    p[1] = x * (cos(Psi * r * pparams[29]) - sin(Psi * r * pparams[29]));
  } else if (i == 46) { // Secant
    p[0] = x;
    p[1] = 1 / (pparams[30] * cos(pparams[30] * r));
  } else if (i == 47) { // Twintrian
    double Psi = R::runif(0, 1);
    double t = log10(pow(sin(Psi * r * pparams[31]), 2)) + cos(Psi * r * pparams[31]);
    p[0] = x * t;
    p[1] = x * (t - M_PI * sin(Psi * r * pparams[31]));
  } else if (i == 48) { // Cross
    double first = sqrt(1 / pow(pow(x, 2) - pow(y, 2), 2));
    p[0] = first * x;
    p[1] = first * y;
  }
  return p;
}

// [[Rcpp::export]]
arma::cube iterate_flame(arma::cube canvas,
                         int iterations,
                         int resolution,
                         int edge,
                         bool blend,
                         bool weighted,
                         bool post,
                         bool final,
                         bool extra,
                         arma::mat colors,
                         Rcpp::DoubleVector functions,
                         Rcpp::DoubleVector funcWeights,
                         arma::mat funcPars,
                         Rcpp::DoubleVector variations,
                         arma::mat varWeights,
                         Rcpp::DoubleVector varParams,
                         arma::mat postPars,
                         Rcpp::DoubleVector finalPars,
                         Rcpp::DoubleVector extraPars) {
  int i, j, indx, indy;
  double x, y, xprev = R::runif(-1, 1), yprev = R::runif(-1, 1), c1 = R::runif(0, 1), c2 = R::runif(0, 1), c3 = R::runif(0, 1);
  Rcpp::DoubleVector tmp(2);
  for (int iter = 1; iter < iterations; iter++) {
    if ((iter % 100) == 0) {
      Rcpp::checkUserInterrupt();
    }
    // Pick an affine function to use according to their weights
    if (weighted) {
      i = Rcpp::sample(functions, 1, false, funcWeights)[0];
    } else {
      i = floor(R::runif(0, functions.length()));
    }
    // Apply the affine function to the current point
    x = funcPars(i, 0) * xprev + funcPars(i, 1) * yprev + funcPars(i, 2);
    y = funcPars(i, 3) * xprev + funcPars(i, 4) * yprev + funcPars(i, 5);
    // Apply the variation(s) to the point
    if (blend) {
      tmp.fill(0);
      for (int j = 0; j < variations.length(); j++) {
        tmp += varWeights(i, j) * variation(x, y, variations[j], funcPars(i, 0), funcPars(i, 1), funcPars(i, 2), funcPars(i, 3), funcPars(i, 4), funcPars(i, 5), varParams);
      }
    } else {
      if (weighted) {
        j = Rcpp::sample(variations, 1, false, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(varWeights.row(i))))[0];
      } else {
        j = floor(R::runif(0, varWeights.n_cols));
      }
      tmp = variation(x, y, variations[j], funcPars(i, 0), funcPars(i, 1), funcPars(i, 2), funcPars(i, 3), funcPars(i, 4), funcPars(i, 5), varParams);
    }
    x = tmp[0];
    y = tmp[1];
    // Apply a post transformation
    if (post) {
      tmp[0] = postPars(i, 0) * x + postPars(i, 1) * y + postPars(i, 2);
      tmp[1] = postPars(i, 3) * x + postPars(i, 4) * y + postPars(i, 5);
      x = tmp[0];
      y = tmp[1];
    }
    // Apply a final transformation
    if (final) {
      tmp[0] = finalPars[0] * x + finalPars[1] * y + finalPars[2];
      tmp[1] = finalPars[3] * x + finalPars[4] * y + finalPars[5];
      x = tmp[0];
      y = tmp[1];
      // Apply an additional post transformation
      if (extra) {
        tmp[0] = extraPars[0] * x + extraPars[1] * y + extraPars[2];
        tmp[1] = extraPars[3] * x + extraPars[4] * y + extraPars[5];
        x = tmp[0];
        y = tmp[1];
      }
    }
    // Update color channels for the current iteration
    c1 = (c1 + colors(i, 0)) / 2;
    c2 = (c2 + colors(i, 1)) / 2;
    c3 = (c3 + colors(i, 2)) / 2;
    // Color the four channels
    if (iter > 20) {
      indx = (x * resolution / (2 * edge)) + resolution / 2;
      if ((indx > 0) && (indx < resolution)) {
        indy = (y * resolution / (2 * edge)) + resolution / 2;
        if ((indy > 0) && (indy < resolution)) {
          canvas(indx, indy, 0) = canvas(indx, indy, 0) + 1;
          canvas(indx, indy, 1) = canvas(indx, indy, 1) + c1;
          canvas(indx, indy, 2) = canvas(indx, indy, 2) + c2;
          canvas(indx, indy, 3) = canvas(indx, indy, 3) + c3;
        }
      }
    }
    xprev = x;
    yprev = y;
  }
  return canvas;
}
