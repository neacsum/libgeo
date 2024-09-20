/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <libgeo/albers.h>
#include <libgeo/errors.h>

///  \file ALBERS.CPP Implementation of \ref Albers "Albers Equal Area" projection

using namespace mlib;

/*!
  Constructor for a Albers Equal Area projection
*/
Albers::Albers (const Params& pp) :
  Projection (pp)
{
  //check required parameters
  if (!pp.contains (Params::Code::NORTH_PARALLEL)
   || !pp.contains (Params::Code::SOUTH_PARALLEL))
    throw erc(GeoError::PARAM, GeoErrors());

  npar_ = north_parallel ();
  spar_ = south_parallel ();
  reflat_ = ref_latitude ();
  reflon_ = ref_longitude ();

  double m1 = ellipsoid().m (npar_);
  double q1 = ellipsoid().q (spar_);
  if (npar_ == spar_)
    n = sin (npar_);
  else
  {
    double m2 = ellipsoid().m (spar_);
    double q2 = ellipsoid().q (spar_);
    n = (m1 * m1 - m2 * m2) / (q2 - q1);
  }
  c_big = m1 * m1 + n * q1;
  rho0 = ellipsoid().a () * sqrt (c_big - n * ellipsoid().q (reflat_)) / n;
  if (n == 0)
    erc (GeoError::PARAM, GeoErrors()).raise();
}

/*!
  Forward conversion from latitude/longitude to XY
*/
erc Albers::GeoXY (double *x, double *y, double lat, double lon) const
{
  double theta = (lon - reflon_) * n;
  double rho = ellipsoid().a () * sqrt (c_big - n * ellipsoid().q (lat)) / n;
  *x = rho * sin (theta);
  *y = rho0 - rho * cos (theta);
  fwd_adj (*x, *y);

  return erc::success;
}

/*!
  Inverse conversion from XY to geographical coordinates.
*/
erc Albers::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);

  double yprime = rho0 - y;
  double rho = hypot (x, yprime);
  double theta = atan (x / yprime);
  *lon = theta / n + reflon_;
  double t = rho*n/(ellipsoid().a ());
  double q = (c_big - t*t) / n;
  if (fabs (q) > 2)
    return erc (GeoError::RANGE, GeoErrors());

  double phi = asin (q / 2);
  int iter = 0;
  double delta = 10;
  while (fabs (delta) > 1e-7 && iter++ < MAX_ITER)
  {
    double sinphi = sin (phi);
    double v1 = 1 - ellipsoid().e2 () * sinphi * sinphi;
    delta = v1 * v1 / (2 * cos (phi))*(q / (1 - ellipsoid().e2 ()) - sinphi / v1 
      + ellipsoid().t (phi));
    phi += delta;
  }
  *lat = phi;
  if (iter >= MAX_ITER)
    return erc (GeoError::NONCONV, GeoErrors());

  return erc::success;
}

double Albers::h (double lat, double lon) const
{
  return sqrt(c_big - n * ellipsoid().q (lat)) / ellipsoid().m (lat);
}

double Albers::k (double lat, double lon) const
{
  if (fabs(lat) == M_PI_2)
    erc (GeoError::RANGE, GeoErrors()).raise();

  return 1/h(lat, lon);
}

