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

double gen_simplex(double x, double y, double z, int seed) {
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("ambient");
    Rcpp::Function f = pkg["gen_simplex"];
    Rcpp::NumericVector result = f( Rcpp::Named("x", x), Rcpp::Named("y", y), Rcpp::Named("z", z), Rcpp::Named("seed", seed));
	double out = result[0];
	return out;
}

Rcpp::DoubleVector splatter_sphere(const double scale) {
  const double r = R::runif(0, 1) * 2 * M_PI;
  Rcpp::DoubleVector out(2);
  out[0] = cos(r) * scale;
  out[1] = sin(r) * scale;
  return out;
}

arma::mat init_particles(arma::mat particles,
                         const int resolution,
                         const int ncols) {
  double startArea = 0.5;
  double maxRadius = 10;
  double scale = resolution / 2;
  int nrows = particles.n_rows;
  for (int i = 0; i < nrows; ++i) {
    particles.at(i, 0) = i + 1; // z
    Rcpp::DoubleVector pos = splatter_sphere(R::runif(0, scale * startArea));
    particles.at(i, 1) = pos[0] + resolution / 2; // x-position
    particles.at(i, 2) = pos[1] + scale; // y-position
    particles.at(i, 3) = R::runif(0.01, maxRadius); // radius
    particles.at(i, 4) = R::runif(1, 500); // duration
    particles.at(i, 5) = R::runif(0, particles.at(i, 4)); // time
    particles.at(i, 6) = R::runif(-1, 1); // x-velocity
    particles.at(i, 7) = R::runif(-1, 1); // y-velocity
    particles.at(i, 8) = R::runif(0.5, 1); // speed
    particles.at(i, 9) = ceil(R::runif(0, ncols)); // color
  }
  return particles;
}

arma::mat reset_particle(arma::mat particles,
                         const int i,
                         const int resolution,
                         const int ncols) {
  double startArea = 0.5;
  double maxRadius = 10;
  double scale = resolution / 2;
  particles.at(i, 0) = arma::max(particles.col(0)) + 1; // z
  Rcpp::DoubleVector pos = splatter_sphere(R::runif(0, scale * startArea));
  particles.at(i, 1) = pos[0] + resolution / 2; // x-position
  particles.at(i, 2) = pos[1] + scale; // y-position
  particles.at(i, 3) = R::runif(0.01, maxRadius); // radius
  particles.at(i, 4) = R::runif(1, 500); // duration
  particles.at(i, 5) = R::runif(0, particles.at(i, 4)); // time
  particles.at(i, 6) = R::runif(-1, 1); // x-velocity
  particles.at(i, 7) = R::runif(-1, 1); // y-velocity
  particles.at(i, 8) = R::runif(0.5, 1); // speed
  particles.at(i, 9) = ceil(R::runif(0, ncols)); // color
  return particles;
}

// [[Rcpp::export]]
arma::mat cpp_splatter(const arma::mat& heightMap,
                       const int& iterations,
                       const int& n,
                       const int& resolution,
                       const int& ncols,
                       double& lwd) {
  arma::mat particles(n, 10);
  arma::mat canvas(iterations * n, 7);
  particles = init_particles(particles, resolution, ncols);
  const double seed = ceil(R::runif(0, 50000));
  int time = 0;
  for (int i = 0; i < iterations; ++i) {
    Rcpp::checkUserInterrupt();
    ++time;
    for (int j = 0; j < n; j++) {
      const double x = particles.at(j, 1);
      const double y = particles.at(j, 2);
      const double fx = fmin(fmax(round(x), 0), resolution - 1);
      const double fy = fmin(fmax(round(y), 0), resolution - 1);
      const double heightValue = heightMap.at(fy, fx) / 255;
      // Calculate line width
      const double s2 = R::runif(0.0001, 0.05);
      const double noise = gen_simplex(x * s2, y * s2, particles.at(j, 4) + time, seed);
      const double r = particles.at(j, 3) * fabs(noise);
      const double width = r * fmin(fmax(heightValue, 0.01), lwd) * (particles.at(j, 5) / particles.at(j, 4));
      // Calculate angle
      double pS = fmin(fmax(heightValue, 0.00001), 0.0001);
      const double angle = gen_simplex(fx * pS, fy * pS, particles.at(j, 4) + time, seed) * M_PI * 2;
      // Calculate speed
      const double speed = particles.at(j, 8) + fmin(fmax(1 - heightValue, 0), 0.1);
      // Update particle velocity
      Rcpp::DoubleVector velocity = {particles.at(j, 6) + cos(angle), particles.at(j, 7) + sin(angle)};
      velocity = velocity / sqrt(sum(pow(velocity, 2)));
      particles.at(j, 6) = velocity[0];
      particles.at(j, 7) = velocity[1];
      // Update particle position
      Rcpp::DoubleVector move = {velocity[0] * speed, velocity[1] * speed};
      particles.at(j, 1) = particles.at(j, 1) + move[0];
      particles.at(j, 2) = particles.at(j, 2) + move[1];
      // Record position
      const int index = i * n + j;
      canvas.at(index, 0) = x; // x position
      canvas.at(index, 1) = y; // y position
      canvas.at(index, 2) = particles.at(j, 1); // end x position
      canvas.at(index, 3) = particles.at(j, 2); // end y position
      canvas.at(index, 4) = particles.at(j, 9); // color
      canvas.at(index, 5) = particles.at(j, 0); // z
      canvas.at(index, 6) = width; // width
      // Update time
      ++particles.at(j, 5);
      if (particles.at(j, 5) > particles.at(j, 4) || (particles.at(j, 1) < 0 || particles.at(j, 1) > resolution) || (particles.at(j, 2) < 0 || particles.at(j, 2) > resolution)) {
        particles = reset_particle(particles, j, resolution, ncols);
      }
    }
  }
  return canvas;
}