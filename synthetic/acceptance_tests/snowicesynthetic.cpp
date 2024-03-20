
#include <math.h>
#include <stdio.h>

#include "snowicesynthetic.hpp"

static double cosineN(int N, double nx);

double tidesynthetic_horizontal_cosine1(double nx, double ny)
{
  return cosineN(1, nx);
}

double tidesynthetic_vertical_cosine1(double nx, double ny)
{
  return cosineN(1, ny);
}

double tidesynthetic_horizontal_cosine2(double nx, double ny)
{
  return cosineN(3, nx);
}
  
double tidesynthetic_vertical_cosine2(double nx, double ny)
{
  return cosineN(3, ny);
}

double tidesynthetic_horizontal_cosine3(double nx, double ny)
{
  return cosineN(5, nx);
}
  
double tidesynthetic_vertical_cosine3(double nx, double ny)
{
  return cosineN(5, ny);
}
  
static double cosineN(int N, double nx)
{
  return (cos(nx * M_PI * (double)N) + 1.0)/2.0;
}

double synthetic_gaussian_snow(double nx, double ny)
{
  double snow_blob1 = 0.5 * exp(-(pow(nx - 0.2, 2) + pow(ny - 0.2, 2)) / 0.1);
  double snow_blob2 = 0.3 * exp(-(pow(nx + 0.3, 2) + pow(ny + 0.3, 2)) / 0.15);
  return snow_blob1 + snow_blob2;
}

double synthetic_gaussian_ice(double nx, double ny)
{
  // Ice model with a larger central blob
  double ice_center = 6 * exp(-(pow(nx, 2) + pow(ny, 2)) / 0.4);

  // Ice model with dispersed pattern
  double ice_dispersed = 0.8 * exp(-((pow(nx - 0.5, 2) + pow(ny - 0.5, 2)) / 0.3));

  return ice_center + ice_dispersed;
}


double tidesynthetic_tas_sea(double nx, double ny)
{
  //
  // Smooth cubic increase with latitude
  //
  double eta = (ny + 1.0)/2.0;
  double eta2 = eta*eta;
  double eta3 = eta2*eta;

  return 12.0*eta2 - 8.0*eta3;
}

double tidesynthetic_tas_land(double nx, double ny)
{
  //
  // Fault line down centre in thirds
  //
  const double LEFT = -2.0;
  const double RIGHT = 5.0;

  double xi = (nx + 1.0)/2.0;
  double eta = (ny + 1.0)/2.0;
  
  if (eta < 0.33) {
    if (xi < 0.40) {
      return LEFT;
    } else {
      return RIGHT;
    }
  } else if (eta < 0.66) {

    double threshold = 0.40 + 0.20*(eta - 0.33)/0.33;
    
    if (xi < threshold) {
      return LEFT;
    } else {
      return RIGHT;
    }
  } else {
    if (xi < 0.60) {
      return LEFT;
    } else {
      return RIGHT;
    }
  }
}

