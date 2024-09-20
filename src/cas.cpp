/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <mlib/defs.h>
#include <libgeo/cas.h>
#include <mlib/ipow.h>

/*!
  \file CAS.CPP  Implementation of \ref Cassini "Cassini-Soldner" projection
*/


using mlib::erc;

/*!
  \class Cassini
  Formulas from Snyder page 95
*/


Cassini::Cassini (const Params& par)
  : Projection{ par }
  , s0{ ellipsoid ().lm (ref_latitude ()) }
{
}

erc Cassini::GeoXY (double *x, double *y, double lat, double lon) const
{
  double a, c, t, tmp;
  double r = ellipsoid().rn (lat);
  a = (lon - ref_longitude())*cos (lat);
  c = ellipsoid().ep2 ()*(tmp = cos (lat))*tmp;
  t = (tmp = tan (lat))*tmp;

  *x = r * (a - t * a*a*a / 6. - (8 - t + 8 * c)*t*(tmp = a * a)*tmp*a / 120.);
  *y = ellipsoid().lm (lat) - s0 + r * tan (lat)*(a*a / 2. + (5 - t + 6 * c)*(tmp = a * a)*tmp / 24.);

  fwd_adj (*x, *y);
  return erc::success;
}

erc Cassini::XYGeo (double x, double y, double *lat, double *lon) const
{
  double s1, mu, eps, phi1, d, t, t2, tmp;
  double e2 = ellipsoid().e2 ();

  inv_adj (x, y);
  s1 = (s0 + y);
  mu = s1 / (ellipsoid().a ()*(1 - e2 / 4. - 3 * e2*e2 / 64. - 5.*e2*e2*e2 / 256.));
  eps = (ellipsoid().a () - ellipsoid().b ()) / (ellipsoid().a () + ellipsoid().b ());

  tmp = eps * eps;
  phi1 = mu
    + (3. / 2. - 27 * eps*eps / 32.)*eps*sin (2 * mu)
    + tmp * (21. / 16. - 55 * tmp / 32.)*sin (4 * mu)
    + 1097.*tmp*tmp / 512.*sin (8 * mu);

  d = x / ellipsoid().rn (phi1);
  t = tan (phi1);
  t2 = t * t;

  tmp = d * d;
  *lat = phi1 - ellipsoid().rn (phi1)*t / ellipsoid().rm (phi1)*(tmp / 2. - (1. - 3.*t2)*tmp*tmp / 24.);
  *lon = ref_longitude() + 1 / cos (phi1)*(d - t2 * (tmp = d * d*d) / 3 + (1 + 3 * t2)*tmp*d*d / 15);
  return erc::success;
}

double Cassini::h (double lat, double lon) const
{
  return 1/sqrt(1-mlib::ipow(cos(lat)*sin(lon-ref_longitude()), 2));
}

