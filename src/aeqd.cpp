/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <mlib/defs.h>
#include <libgeo/aeqd.h>

/// \file aeqd.cpp Implementation of \ref AzEqDist "Azimuthal Equidistant"

using mlib::erc;

AzEqDist::AzEqDist (const Params& par) :
  Projection (par)
{
  sphi1 = sin (ref_latitude());
  cphi1 = cos (ref_latitude());
  n1 = ellipsoid().rn (ref_latitude());
  g = ellipsoid().ep ()*sphi1;
  north_polar = (fabs (ref_latitude() - M_PI_2) < 1.e-6);
  south_polar = (fabs (ref_latitude() + M_PI_2) < 1.e-6);
}

erc AzEqDist::GeoXY (double *x, double *y, double lat, double lon) const
{
  double cphi = cos (lat);
  double dl = lon - ref_longitude();
  double psi = atan ((1.0 - ellipsoid().e2 ())*tan (lat) 
    + ellipsoid().e2 ()*n1*sphi1 / (ellipsoid().rn (lat) * cphi));
  double az = atan2 (sin (dl), cphi1 * tan (psi) - sphi1 * cos (dl));
  double caz = cos (az);
  double saz = sin (az);

  double s;
  if (saz == 0.0)
    s = asin (cphi1*sin (psi) - sphi1 * cos (psi)) * ((caz < 0) ? -1. : 1.);
  else
    s = asin (sin (dl)*cos (psi) / saz);

  double h = ellipsoid().ep ()*cphi1* caz;
  double c = n1 * s*(1.0 - s * s*h*h*(1.0 - h * h) / 6.0 + s * s*s / 8.*g*h*(1. - 2.*h*h) +
    s * s*s*s / 120.0*(h*h*(4. - 7.*h*h) - 3 * g*g*(1. - 7.*h*h)) -
    s * s*s*s*s / 48.*g*h);

  *x = c * saz;
  *y = c * caz;
  fwd_adj (*x, *y);
  return erc::success;
}

erc AzEqDist::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  if (x == 0. && y == 0.)
  {
    //converting projection center - we already know the answer
    *lat = ref_latitude();
    *lon = ref_longitude();
  }
  else
  {
    double c = hypot (x, y);
    double az = atan2 (x, y);
    double caz = cos (az);
    double saz = sin (az);
    double A = -ellipsoid().ep2 ()*cphi1*cphi1*caz*caz;
    double B = 3.*ellipsoid().ep2 ()*(1. - A)*sphi1*cphi1*caz;
    double D = c / n1;
    double E = D - A * (1. + A)*D*D*D / 6. - B * (1. + 3.*A)*D*D*D*D / 24.;
    double F = 1. - A * E*E / 2. - B * E*E*E / 6.;
    double psi = asin (sphi1*cos (E) + cphi1 * sin (E)*caz);

    if (north_polar)
      *lon = ref_longitude() + atan2 (x, -y);
    else if (south_polar)
      *lon = ref_longitude() + atan2 (x, y);
    else
      *lon = lon_adjust (ref_longitude() + asin (saz*sin (E) / cos (psi)));
    *lat = atan ((1. - ellipsoid().e2 ()*F*sphi1 / sin (psi))*tan (psi) / (1 - ellipsoid().e2 ()));
  }
  return erc::success;
}

double AzEqDist::h (double lat, double lon) const
{
  return 1.;
}

double AzEqDist::k (double lat, double lon) const
{
  if (north_polar || south_polar)
    return (ellipsoid().lm (M_PI_2) - ellipsoid().lm (lat)) / (ellipsoid().a () * ellipsoid().m (lat));

  return 1.;
}

