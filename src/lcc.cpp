/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/

#include <libgeo/lcc.h>
#include <libgeo/errors.h>

///  \file LCC.CPP Implementation of \ref Lambert "Lambert Conformal Conical"

using mlib::erc;

/*!
  \class Lambert
  Formulas from Snyder page 107-109
*/

#define MAX_ITER  10    //max number of iterations in reverse formulas

/*!
  Ctor for a %Lambert Conformal Conical projection
*/
Lambert::Lambert (const Params& pp) :
  Projection (pp)
{
  //check required parameters

  //at least one of north or south parallel must be specified
  if (!pp.contains (Params::Code::NORTH_PARALLEL)
   && !pp.contains (Params::Code::SOUTH_PARALLEL))
    throw erc (GeoError::PARAM, GeoErrors ());

  //if one is missing use the same value as the other one
  if (!pp.contains (Params::Code::NORTH_PARALLEL))
    par_.north_latitude (south_parallel ());
  else if (!pp.contains (Params::Code::SOUTH_PARALLEL))
    par_.south_latitude (north_parallel ());

  double m1 = cos (north_parallel ()) / sqrt (1 - ellipsoid ().e2 () 
    * sin (north_parallel ()) * sin (north_parallel ()));
  double t1 = tfunc (north_parallel());
  if (north_parallel() == south_parallel())
    n = sin (north_parallel());
  else
  {
    double m2 = cos (south_parallel()) / sqrt (1 - ellipsoid().e2 ()*sin (south_parallel())*sin (south_parallel()));
    n = log (m1 / m2) / log (t1 / tfunc (south_parallel()));
  }
  af_big = ellipsoid().a () * k0() * m1 / (n*pow (t1, n));
  rho0 = af_big * pow (tfunc (ref_latitude()), n);
  if (n == 0)
    erc (GeoError::PARAM, GeoErrors()).raise();
}

/*!
  Forward conversion from latitude/longitude to XY
*/
erc Lambert::GeoXY(double *x, double *y, double lat, double lon) const
{
  double theta = (lon - ref_longitude()) * n;
  double rho;
  if (fabs(lat) == M_PI_2)
  {
    if (n*lat > 0.)
      rho = 0;
    else
      return erc (GeoError::RANGE, GeoErrors());
  }
  else
    rho = af_big * pow(tfunc(lat), n);
  *x = rho * sin(theta);
  *y = rho0 - rho * cos(theta);
  fwd_adj (*x, *y);
  return erc::success;
}

/*!
  Inverse conversion from XY to geographical coordinates.
*/
erc Lambert::XYGeo(double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  double yprime = rho0-y;
  double rho = hypot (x,yprime);
  double theta;
  if (n < 0.)
  {
    rho *=-1;
    theta = atan(-x/(y-rho0));
  }
  else
    theta = atan(x/yprime);
  *lon = theta/n + ref_longitude();
  double t = pow( rho/af_big, 1/n);

  double e_ = ellipsoid().e();
  double phi_right = M_PI_2 - 2*atan(t);
  int iter =0;
  double delta = 10;
  while (delta > 1e-10 && iter++ < MAX_ITER)
  {
    double sinphi = sin(phi_right);
    double phi_left = M_PI_2 - 2*atan(t*pow((1-e_*sinphi)/(1+e_*sinphi),e_/2));
    delta = fabs(phi_left-phi_right);
    phi_right = phi_left;
  }
  *lat = phi_right;
  if (iter >= MAX_ITER)
    return erc (GeoError::NONCONV, GeoErrors());

  return erc::success;
}

double Lambert::k (double lat, double lon) const
{
  double rho = af_big * pow(tfunc(lat), n);
  return rho*n/(ellipsoid().a()*ellipsoid().m(lat));
}

/// helper function. Not to be confused with the Ellipsoid::t function
double Lambert::tfunc (double phi) const
{
  double sval = sin(phi);
  double t1 = (1. - sval) / (1. + sval);
  double e_ = ellipsoid().e();
  double t2 = pow((1. + e_ * sval) / (1. - e_ * sval), e_);
  return (sqrt(t1 * t2));
}

double Lambert::Convergence (double lat, double lon) const
{
  return -n*(lon-ref_longitude());
}
