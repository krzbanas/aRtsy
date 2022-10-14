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

void transform(double& x,
               double& y,
               const double& a,
               const double& b,
               const double& c,
               const double& d,
               const double& e,
               const double& f) {
  const double newx = a * x + b * y + c;
  const double newy = d * x + e * y + f;
  x = newx;
  y = newy;
}

void variation(double& x,
               double& y,
               const int& i,
               const double& a,
               const double& b,
               const double& c,
               const double& d,
               const double& e,
               const double& f,
               const Rcpp::DoubleVector& pparams) {
  if (i == 0) { // Linear
    double newx = x;
    double newy = y;
    x = newx;
    y = newy;
  } else if (i == 1) { // Sine
    double newx = sin(x);
    double newy = sin(y);
    x = newx;
    y = newy;
  } else if (i == 2) { // Sperical
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = x / pow(r, 2);
    double newy = y / pow(r, 2);
    x = newx;
    y = newy;
  } else if (i == 3) { // Swirl
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = x * sin(pow(r, 2)) - y * cos(pow(r, 2));
    double newy = x * cos(pow(r, 2)) + y * sin(pow(r, 2));
    x = newx;
    y = newy;
  } else if (i == 4) { // Horsehoe
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = (1 / r) * ((x - y) * (x + y));
    double newy = (1 / r) * (2 * x * y);
    x = newx;
    y = newy;
  } else if (i == 5) { // Polar
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double newx = theta / M_PI;
    double newy = r - 1;
    x = newx;
    y = newy;
  } else if (i == 6) { // Handkerchief
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double newx = r * sin(theta + r);
    double newy = r * cos(theta - r);
    x = newx;
    y = newy;
  } else if (i == 7) { // Heart
    double theta = atan(x / y);
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = r * sin(theta * r);
    double newy = r * (-cos(theta * r));
    x = newx;
    y = newy;
  } else if (i == 8) { // Disc
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double newx = theta / M_PI * sin(M_PI * r);
    double newy = theta / M_PI * cos(M_PI * r);
    x = newx;
    y = newy;
  } else if (i == 9) { // Spiral
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double newx = 1 / r * (cos(theta) + sin(r));
    double newy = 1 / r * (sin(theta) + cos(r));
    x = newx;
    y = newy;
  } else if (i == 10) { // Hyperbolic
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double newx = sin(theta) / r;
    double newy = r * cos(theta);
    x = newx;
    y = newy;
  } else if (i == 11) { // Diamond
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double newx = sin(theta) * cos(r);
    double newy = cos(theta) * sin(r);
    x = newx;
    y = newy;
  } else if (i == 12) { // Ex
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double p0 = sin(theta + r);
    double p1 = cos(theta - r);
    double newx = r * (pow(p0, 3) + pow(p1, 3));
    double newy = r * (pow(p0, 3) - pow(p1, 3));
    x = newx;
    y = newy;
  } else if (i == 13) { // Julia
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double s = R::runif(0, 1);
    double Omega = (s < .05) ? 0 : M_PI;
    double newx = sqrt(r) * cos(theta / 2 + Omega);
    double newy = sqrt(r) * sin(theta / 2 + Omega);
    x = newx;
    y = newy;
  } else if (i == 14) { // Bent
    double newx, newy;
    if ((x >= 0) && (y >= 0)) {
      newx = x;
      newy = y;
    } else if ((x < 0) && (y >= 0)) {
      newx = 2 * x;
      newy = y;
    } else if ((x >= 0) && (y < 0)) {
      newx = x;
      newy = y / 2;
    } else if ((x < 0) && (y < 0)) {
      newx = 2 * x;
      newy = y / 2;
    }
    x = newx;
    y = newy;
  } else if (i == 15) { // Waves
    double newx = x + b * sin(y / pow(c, 2));
    double newy = y + e * sin(x / pow(f, 2));
    x = newx;
    y = newy;
  } else if (i == 16) { // Fisheye
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = (2 / (r + 1)) * y;
    double newy = (2 / (r + 1)) * x;
    x = newx;
    y = newy;
  } else if (i == 17) { // Popcorn
    double newx = x + c * sin(tan(3 * y));
    double newy = y + f * sin(tan(3 * x));
    x = newx;
    y = newy;
  } else if (i == 18) { // Exponential
    double newx = exp(x - 1) * cos(M_PI * y);
    double newy = exp(x - 1) * sin(M_PI * y);
    x = newx;
    y = newy;
  } else if (i == 19) { // Power
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double newx = pow(r, sin(theta)) * cos(theta);
    double newy = pow(r, sin(theta)) * sin(theta);
    x = newx;
    y = newy;
  } else if (i == 20) { // Cosine
    double newx = cos(M_PI * x) * cosh(y);
    double newy = -(sin(M_PI * x) * sinh(y));
    x = newx;
    y = newy;
  } else if (i == 21) { // Rings
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double first = fmod(r + pow(c, 2), 2 * pow(c, 2)) - pow(c, 2) + r * (1 - pow(c, 2));
    double newx = first * cos(theta);
    double newy = first * sin(theta);
    x = newx;
    y = newy;
  } else if (i == 22) { // Fan
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double t = M_PI * pow(c, 2);
    double newx = fmod(theta + f, t) > (t / 2) ? r * cos(theta - (t/2)) : r * cos(theta + (t/2));
    double newy = fmod(theta + f, t) > (t / 2) ? r * sin(theta - (t/2)) : r * sin(theta + (t/2));
    x = newx;
    y = newy;
  } else if (i == 23) { // Blob
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double first = r * (pparams[1] + ((pparams[0] - pparams[1])/ 2) * (sin(pparams[2] * theta) + 1));
    double newx = first * cos(theta);
    double newy = first * sin(theta);
    x = newx;
    y = newy;
  } else if (i == 24) { // PDJ
    double newx = sin(pparams[3] * y) - cos(pparams[4] * x);
    double newy = sin(pparams[5] * x) - cos(pparams[6] * y);
    x = newx;
    y = newy;
  } else if (i == 25) { // Fan2
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    variation(x, y, 22, a, b, c, d, e, f, pparams);
    double p1 = M_PI * pow(x, 2);
    double p2 = y;
    double t = theta + p2 - p1 * trunc((2 * theta * p2) / p1);
    double newx = t > (p1 / 2) ? r * sin(theta - (p1 / 2)) : r * sin(theta + (p1 / 2));
    double newy = t > (p1 / 2) ? r * cos(theta - (p1 / 2)) : r * cos(theta + (p1 / 2));
    x = newx;
    y = newy;
  } else if (i == 26) { // Rings2
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double theta = atan(x / y);
    double pp = pow(pparams[7], 2);
    double t = r - 2 * pp * trunc((r + pp) / (2 * pp)) + r * (1 - pp);
    double newx = t * sin(theta);
    double newy = t * cos(theta);
    x = newx;
    y = newy;
  } else if (i == 27) { // Eyefish
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = (2 / (r + 1)) * x;
    double newy = (2 / (r + 1)) * y;
    x = newx;
    y = newy;
  } else if (i == 28) { // Bubble
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = (4 / (pow(r, 2) + 4)) * x;
    double newy = (4 / (pow(r, 2) + 4)) * y;
    x = newx;
    y = newy;
  } else if (i == 29) { // Cylinder
    double newx = sin(x);
    double newy = y;
    x = newx;
    y = newy;
  } else if (i == 30) { // Perspective
    double first = (pparams[9] / (pparams[9] - y * sin(pparams[8])));
    double newx = first * x;
    double newy = first * (y * cos(pparams[8]));
    x = newx;
    y = newy;
  } else if (i == 31) { // Noise
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double newx = Psi1 * (x * cos(2 * M_PI * Psi2));
    double newy = Psi1 * (y * sin(2 * M_PI * Psi2));
    x = newx;
    y = newy;
  } else if (i == 32) { // JuliaN
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double phi = atan(y / x);
    double Psi = R::runif(0, 1);
    double p3 = trunc(fabs(pparams[11]) * Psi);
    double t = (phi + 2 * M_PI * p3) / pparams[11];
    double newx = pow(r, pparams[11] / pparams[10]) * cos(t);
    double newy = pow(r, pparams[11] / pparams[10]) * sin(t);
    x = newx;
    y = newy;
  } else if (i == 33) { // JuliaScope
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double phi = atan(y / x);
    double Psi = R::runif(0, 1);
    int Lambda = floor(R::runif(0, 2));
    double p3 = trunc(fabs(pparams[12] * Psi));
    double t = (Lambda * phi + 2 * M_PI * p3) / pparams[12];
    double newx = pow(r, pparams[13] / pparams[12]) * cos(t);
    double newy = pow(r, pparams[13] / pparams[12]) * sin(t);
    x = newx;
    y = newy;
  } else if (i == 34) { // Blur
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double newx = Psi1 * cos(2 * M_PI * Psi2);
    double newy = Psi1 * sin(2 * M_PI * Psi2);
    x = newx;
    y = newy;
  } else if (i == 35) { // Gaussian
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double Psi3 = R::runif(0, 1);
    double Psi4 = R::runif(0, 1);
    double Psi5 = R::runif(0, 1);
    double newx = (Psi1 + Psi2 + Psi3 + Psi4 - 2) * cos(2 * M_PI * Psi5);
    double newy = (Psi1 + Psi2 + Psi3 + Psi4 - 2) * sin(2 * M_PI * Psi5);
    x = newx;
    y = newy;
  } else if (i == 36) { // RadialBlur
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double phi = atan(y / x);
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double Psi3 = R::runif(0, 1);
    double Psi4 = R::runif(0, 1);
    double p1 = pparams[14] * (M_PI / 2);
    double t1 = pparams[15] * (Psi1 + Psi2 + Psi3 + Psi4 - 2);
    double t2 = phi + t1 * sin(p1);
    double t3 = t1 * cos(p1) - 1;
    double newx = (1 / pparams[15]) * (r * cos(t2) + t3 * x);
    double newy = (1 / pparams[15]) * (r * sin(t2) + t3 * y);
    x = newx;
    y = newy;
  } else if (i == 37) { // Pie
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double Psi3 = R::runif(0, 1);
    double t1 = trunc(Psi1 * pparams[16] + 0.5);
    double t2 = pparams[17] + ((2 * M_PI) / pparams[16]) * (t1 + Psi2 * pparams[18]);
    double newx = Psi3 * cos(t2);
    double newy = Psi3 * sin(t2);
    x = newx;
    y = newy;
  } else if (i == 38) { // Ngon
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double phi = atan(y / x);
    double p2 = 2 * M_PI / pparams[20];
    double t3 = phi - p2 * floor(phi / p2);
    double t4;
    if (t3 > (p2 / 2)) {
      t4 = t3;
    } else {
      t4 = t3 - p2;
    }
    double k = (pparams[21] * (1 / cos(t4) - 1) + pparams[22]) / pow(r, pparams[19]);
    double newx = k * x;
    double newy = k * y;
    x = newx;
    y = newy;
  } else if (i == 39) { // Curl
    double t1 = 1 + pparams[23] * x + pparams[24] * (pow(x, 2) - pow(y, 2));
    double t2 = pparams[23] * y + 2 * pparams[24] * x * y;
    double first = 1 / (pow(t1, 2) + pow(t2, 2));
    double newx = first * (x * t1 + y * t2);
    double newy = first * (y * t1 - x * t2);
    x = newx;
    y = newy;
  } else if (i == 40) { // Rectangles
    double newx = (2 * floor(x / pparams[25]) + 1) * pparams[25] - x;
    double newy = (2 * floor(y / pparams[26]) + 1) * pparams[26] - y;
    x = newx;
    y = newy;
  } else if (i == 41) { // Arch
    double Psi = R::runif(0, 1);
    double newx = sin(Psi * M_PI * pparams[27]);
    double newy = pow(sin(Psi * M_PI * pparams[27]), 2) / cos(Psi * M_PI * pparams[27]);
    x = newx;
    y = newy;
  } else if (i == 42) { // Tangent
    double newx = sin(x) / cos(y);
    double newy = tan(y);
    x = newx;
    y = newy;
  } else if (i == 43) { // Square
    double Psi1 = R::runif(0, 1);
    double Psi2 = R::runif(0, 1);
    double newx = Psi1 - 0.5;
    double newy = Psi2 - 0.5;
    x = newx;
    y = newy;
  } else if (i == 44) { // Rays
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double Psi = R::runif(0, 1);
    double first = (pparams[28] * tan(Psi * M_PI * pparams[28])) / pow(r, 2);
    double newx = first * cos(x);
    double newy = first * sin(x);
    x = newx;
    y = newy;
  } else if (i == 45) { // Blade
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double Psi = R::runif(0, 1);
    double newx = x * (cos(Psi * r * pparams[29]) + sin(Psi * r * pparams[29]));
    double newy = x * (cos(Psi * r * pparams[29]) - sin(Psi * r * pparams[29]));
    x = newx;
    y = newy;
  } else if (i == 46) { // Secant
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double newx = x;
    double newy = 1 / (pparams[30] * cos(pparams[30] * r));
    x = newx;
    y = newy;
  } else if (i == 47) { // Twintrian
    double r = sqrt(pow(x, 2) + pow(y, 2));
    double Psi = R::runif(0, 1);
    double t = log10(pow(sin(Psi * r * pparams[31]), 2)) + cos(Psi * r * pparams[31]);
    double newx = x * t;
    double newy = x * (t - M_PI * sin(Psi * r * pparams[31]));
    x = newx;
    y = newy;
  } else if (i == 48) { // Cross
    double first = sqrt(1 / pow(pow(x, 2) - pow(y, 2), 2));
    double newx = first * x;
    double newy = first * y;
    x = newx;
    y = newy;
  }
}

// [[Rcpp::export]]
arma::cube iterate_flame(arma::cube& canvas,
                         const int& iterations,
                         const int& resolution,
                         const int& edge,
                         const bool& blend,
                         const bool& weighted,
                         const bool& post,
                         const bool& final,
                         const bool& extra,
                         const arma::mat& colors,
                         const Rcpp::DoubleVector& functions,
                         const Rcpp::DoubleVector& funcWeights,
                         const arma::mat& funcPars,
                         const Rcpp::DoubleVector& variations,
                         const arma::mat& varWeights,
                         const Rcpp::DoubleVector& varParams,
                         const arma::mat& postPars,
                         const Rcpp::DoubleVector& finalPars,
                         const Rcpp::DoubleVector& extraPars,
                         const int& bsym) {
  const int nvar = variations.length(), nfunc = functions.length();
  double x = R::runif(-1, 1), y = R::runif(-1, 1), c1 = R::runif(0, 1), c2 = R::runif(0, 1), c3 = R::runif(0, 1);
  bool vary = !((nvar == 1) && (variations[0] == 0));
  for (int iter = 1; iter < iterations; ++iter) {
    if (iter % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
    // Pick an affine function to use and apply to the current point
    const int i = weighted ? Rcpp::sample(functions, 1, false, funcWeights)[0] : floor(R::runif(0, nfunc));
    transform(x, y, funcPars.at(i, 0), funcPars.at(i, 1), funcPars.at(i, 2), funcPars.at(i, 3), funcPars.at(i, 4), funcPars.at(i, 5));
    // Apply variations
    if (i < bsym) { // Functions with i < bsym are affine functions, the rest is symmtry functions so we skip
      if (vary) { // Do not vary if the only affine is linear
        if (blend) { // Blend variations
          double xc = 0, yc = 0;
          for (int j = 0; j < nvar; j++) {
            double& xp = x, yp = y;
            variation(xp, yp, variations[j], funcPars.at(i, 0), funcPars.at(i, 1), funcPars.at(i, 2), funcPars.at(i, 3), funcPars.at(i, 4), funcPars.at(i, 5), varParams);
            xc += varWeights.at(i, j) * xp;
            yc += varWeights.at(i, j) * yp;
          }
          x = xc, y = yc;
        } else { // Do not blend variations
          const int j = weighted ? Rcpp::sample(variations, 1, false, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(varWeights.row(i))))[0] : floor(R::runif(0, nvar));
          variation(x, y, variations[j], funcPars.at(i, 0), funcPars.at(i, 1), funcPars.at(i, 2), funcPars.at(i, 3), funcPars.at(i, 4), funcPars.at(i, 5), varParams);
        }
      }
      // Apply a post transformation
      if (post) {
        transform(x, y, postPars.at(i, 0), postPars.at(i, 1), postPars.at(i, 2), postPars.at(i, 3), postPars.at(i, 4), postPars.at(i, 5));
      }
      // Apply a final transformation
      if (final) {
        transform(x, y, finalPars[0], finalPars[1], finalPars[2], finalPars[3], finalPars[4], finalPars[5]);
        // Apply an additional post transformation
        if (extra) {
          transform(x, y, extraPars[0], extraPars[1], extraPars[2], extraPars[3], extraPars[4], extraPars[5]);
        }
      }
      // Update color channels for the current iteration
      c1 = (c1 + colors.at(i, 0)) / 2;
      c2 = (c2 + colors.at(i, 1)) / 2;
      c3 = (c3 + colors.at(i, 2)) / 2;
    }
    // Update the cube data structure
    if (iter > 20) {
      const int indx = (x * resolution / (2 * edge)) + resolution / 2;
      if ((indx >= 0) && (indx < resolution)) {
        const int indy = (y * resolution / (2 * edge)) + resolution / 2;
        if ((indy >= 0) && (indy < resolution)) {
          ++canvas.at(indx, indy, 0);
          canvas.at(indx, indy, 1) = canvas.at(indx, indy, 1) + c1;
          canvas.at(indx, indy, 2) = canvas.at(indx, indy, 2) + c2;
          canvas.at(indx, indy, 3) = canvas.at(indx, indy, 3) + c3;
        }
      }
    }
  }
  return canvas;
}
