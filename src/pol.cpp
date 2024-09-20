/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <libgeo/pol.h>
#include <libgeo/errors.h>

/// \file POL.CPP Implementation of Polyconic projection

using mlib::erc;
/*!
  \class Polyconic 

  Formulas from Snyder pag. 129-131
*/

Polyconic::Polyconic (const Params& par)
  : Projection (par)
{
  m0 = ellipsoid().lm (ref_latitude());
  c[0] = ((-0.01953125*ellipsoid().e2 () - 0.046875)*ellipsoid().e2 () - 0.250)* ellipsoid().e2 () + 1;
  c[1] = ((0.0439453125* ellipsoid().e2 () + 0.093750)* ellipsoid().e2 () + 0.375)* ellipsoid().e2 () * 2;
  c[2] = (0.0439453125* ellipsoid().e2 () + 0.05859375)* ellipsoid().e2 ()* ellipsoid().e2 () * 4;
  c[3] = -0.01139322916667* ellipsoid().e2 ()* ellipsoid().e2 ()* ellipsoid().e2 () * 6;
}

erc Polyconic::GeoXY (double *x, double *y, double lat, double lon) const
{
  double l = lon - ref_longitude();
  if (lat == 0.)
  {
    *x = ellipsoid().a ()*l;
    *y = -m0;
  }
  else
  {
    double ee = l * sin (lat);
    double n = ellipsoid().rn (lat) / tan (lat);
    *x = n * sin (ee);
    *y = ellipsoid().lm (lat) - m0 + n * (1 - cos (ee));
  }
  fwd_adj (*x, *y);
  return erc::success;
}

erc Polyconic::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  double x_a = x / ellipsoid().a ();
  if (y == -m0)
  {
    *lat = 0.;
    *lon = x_a + ref_longitude();
  }
  else
  {
    double aa = (m0 + y) / ellipsoid().a ();
    double bb = x_a * x_a + aa * aa;
    double cc;
    double phin;
    double phin1 = aa;
    double slat;
    int iter = 0;
    do {
      phin = phin1;
      cc = sqrt (1 - ellipsoid().e2 ()*(slat = sin (phin))*slat)*tan (phin);
      double ma = ellipsoid().lm (phin) / ellipsoid().a ();
      double mn1 = m1 (phin);
      phin1 = phin - (aa*(cc*ma + 1) - ma - 0.5*(ma*ma + bb)*cc) / (ellipsoid().e2 ()*sin (2 * phin)*(ma*ma + bb - 2 * aa*ma) / (4 * cc) + (aa - ma)*(cc*mn1 - 2 / sin (2 * phin)) - mn1);
    } while (fabs (phin1 - phin) > 1e-10 && iter++ < 10);

    if (iter >= MAX_ITER)
      return erc (GeoError::NONCONV, GeoErrors());

    *lat = phin1;
    *lon = asin (x*cc / ellipsoid().a ()) / sin (phin1) + ref_longitude ();
  }
  return erc::success;
}

double Polyconic::h (double lat, double lon) const
{
  double e2_ = ellipsoid().e2 ();
  double scale;
  double l = lon - ref_longitude();
  if (lat != 0.)
  {
    double slat = sin (lat);
    double clat = cos (lat);
    double ee = l * slat;
    double d = atan ((ee - sin (ee)) / (1 / (clat*clat) - cos (ee) - e2_ * slat*slat / (1 - e2_ * slat*slat)));
    scale = (1 - e2_ + 2 * (1 - e2_ * slat*slat)*sin (ee / 2)*sin (ee / 2) / (tan (lat)*tan (lat))) / ((1 - e2_)*cos (d));
  }
  else
    scale = (m1 (lat) + l * l / 2.) / (1 - e2_);

  return scale;
}

double Polyconic::m1 (double lat) const
{
  return c[0] - c[1] * cos (2 * lat) + c[2] * cos (4 * lat) + c[3] * cos (6 * lat);
}