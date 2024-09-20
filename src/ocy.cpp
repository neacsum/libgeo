/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <mlib/defs.h>
#include <libgeo/ocy.h>

/*!
  \file OCY.CPP Implementation of \ref ObliqueCylindrical "Oblique Cylindrical"
  projection
*/

using mlib::erc;

ObliqueCylindrical::ObliqueCylindrical (const Params& par) :
  Projection (par)
{
  c1 = sqrt (1 + ellipsoid().e2 () / (1 - ellipsoid().e2 ()) * pow (cos (ref_latitude()), 4));
  chi0 = asin (sin (ref_latitude()) / c1);
  double b0prim = asin (ellipsoid().e () * sin (ref_latitude()));
  c2 = log (tan (M_PI_4 + chi0 / 2.)) - c1 * log (tan (M_PI_4 + ref_latitude() / 2.)) +
    c1 * ellipsoid().e () * log (tan (M_PI_4 + b0prim / 2));
  double n0 = ellipsoid().a () / sqrt (1 - ellipsoid().e2 () * sin (ref_latitude()) * sin (ref_latitude()));
  double m0 = ellipsoid().a () * (1 - ellipsoid().e2 ()) / pow (1 - ellipsoid().e2 () * sin (ref_latitude()) * sin (ref_latitude()), 1.5);
  map_radius = k0 () * sqrt (n0 * m0);
}

erc ObliqueCylindrical::GeoXY (double* x, double* y, double lat, double lon) const
{
  double l = c1 * (lon - ref_longitude());
  double bprim = asin (ellipsoid().e () * sin (lat));
  double tg = exp (c1 * log (tan (M_PI_4 + lat / 2)) - c1 * ellipsoid().e () * log (tan (M_PI_4 + bprim / 2)) + c2);
  double b = 2 * atan (tg) - M_PI_2;

  double bbar = asin (cos (chi0) * sin (b) - sin (chi0) * cos (b) * cos (l));
  double lbar = asin (cos (b) * sin (l) / cos (bbar));
  *x = map_radius * lbar;
  *y = map_radius * log (tan (M_PI_4 + bbar / 2));
  fwd_adj (*x, *y);
  return erc::success;
}

erc ObliqueCylindrical::XYGeo (double x, double y, double* lat, double* lon) const
{
  inv_adj (x, y);
  double xbar = exp (y / map_radius);
  double lbar = x / map_radius;
  double bbar = 2 * atan (xbar) - M_PI_2;
  double b = asin (cos (chi0) * sin (bbar) + sin (chi0) * cos (bbar) * cos (lbar));
  double l = asin (cos (bbar) * sin (lbar) / cos (b));

  double phi = ref_latitude();
  double delta = 100.;
  double temp1 = (log (tan (M_PI_4 + b / 2)) - c2) / c1;
  while (delta > 1e-7)
  {
    double bprim = asin (ellipsoid().e () * sin (phi));
    double tg = exp (temp1 + ellipsoid().e () * log (tan (M_PI_4 + bprim / 2)));
    double phi_new = 2 * atan (tg) - M_PI_2;
    delta = fabs (phi - phi_new);
    phi = phi_new;
  }
  *lat = phi;
  *lon = ref_longitude() + l / c1;
  return erc::success;
}

double ObliqueCylindrical::k (double lat, double lon) const
{
  double l = c1 * (lon - ref_longitude());
  double bprim = asin (ellipsoid().e () * sin (lat));
  double tg = exp (c1 * log (tan (M_PI_4 + lat / 2)) - c1 * ellipsoid().e () * log (tan (M_PI_4 + bprim / 2)) + c2);
  double b = 2 * atan (tg) - M_PI_2;

  double bbar = asin (cos (chi0) * sin (b) - sin (chi0) * cos (b) * cos (l));

  double kval = c1 * ellipsoid().rn (ref_latitude()) / ellipsoid().rn (lat) * cos (b) / (cos (lat) * cos (bbar));
  return kval;
}

