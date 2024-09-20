/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <libgeo/ome.h>
#include <libgeo/errors.h>

///  \file OME.CPP Implementation of \ref ObliqueMercator "Oblique Mercator"

using mlib::erc;

ObliqueMercator::ObliqueMercator (const Params& par) 
  : Projection (par)
{
  B = sqrt (1 + ellipsoid().e2 ()*pow (cos (ref_latitude()), 4) / (1 - ellipsoid().e2 ()));
  double v1 = sqrt (1 - ellipsoid().e2 ());
  double v2 = 1 - ellipsoid().e2 ()*sin (ref_latitude())*sin (ref_latitude());
  A = ellipsoid().a ()* B * k0 () * v1 / v2;
  double D = B * v1 / (cos (ref_latitude())*sqrt (v2));
  double F = D + ((ref_latitude() >= 0) ? sqrt (D*D - 1) : -sqrt (D*D - 1));
  double t0 = exptau (ref_latitude());
  E = F * pow (t0, B);
  double G = (F - 1. / F) / 2;
  gamma0 = asin (sin (skew_azimuth ()) / D);
  lam1 = ref_longitude() - asin (G*tan (gamma0)) / B;
}

erc ObliqueMercator::GeoXY (double *x, double *y, double lat, double lon) const
{
  double Q = E / pow (exptau (lat), B);
  double S = (Q - 1. / Q) / 2.;
  double T = (Q + 1. / Q) / 2.;
  double V = sin (B*(lon - lam1));
  double U = (-V * cos (gamma0) + S * sin (gamma0)) / T;
  double v = A * log ((1. - U) / (1. + U)) / (2.*B);
  double u = A / B * atan2 (S*cos (gamma0) + V * sin (gamma0), cos (B*(lon - lam1)));
  deskew (u, v, x, y);
  return erc::success;
}

double ObliqueMercator::k (double lat, double lon) const
{
  double scale;
  double Q = E / pow (exptau (lat), B);
  double S = (Q - 1. / Q) / 2.;
  double T = (Q + 1. / Q) / 2.;
  double V = sin (B*(lon - lam1));
  double U = (-V * cos (gamma0) + S * sin (gamma0)) / T;
  double v = A * log ((1. - U) / (1. + U)) / (2.*B);
  double u = A / B * atan2 (S*cos (gamma0) + V * sin (gamma0), cos (B*(lon - lam1)));
  scale = A * cos (B*u / A)*sqrt (1 - ellipsoid().e2 ()*sin (lat)*sin (lat)) / ellipsoid().a ()*cos (lat)*cos (B*(lon - ref_longitude()));
  return scale;
}

erc ObliqueMercator::XYGeo (double x, double y, double *lat, double *lon) const
{
  double u, v;
  skew (&u, &v, x, y);
  return uvgeo (u, v, lat, lon);
}

/// helper function
double ObliqueMercator::exptau (double value) const
{
  value = sin (value);
  double t1 = (1. - value) / (1. + value);
  double e_val = ellipsoid().e ();						// don't calculate it every time
  double t2 = pow ((1. + e_val * value) / (1. - e_val * value), e_val);

  double result = sqrt (t1*t2);
  return result;
}


/// Conversion from uv (skewed coords) to geographic
erc ObliqueMercator::uvgeo (double u, double v, double *lat, double *lon) const
{
  double Q = exp (-B / A * v);
  double S = (Q - 1. / Q) / 2.;
  double T = (Q + 1. / Q) / 2.;
  double V = sin (B*u / A);
  double U = (V*cos (gamma0) + S * sin (gamma0)) / T;
  double t = pow (E / sqrt ((1. + U) / (1. - U)), 1. / B);
  double ll = M_PI / 2 - 2 * atan (t);
  int i = 0;
  double err = 1.;
  double e_ = ellipsoid().e ();
  while (err > 1e-9 && i < 30)
  {
    double ll1 = M_PI / 2 - 2 * atan (t*pow ((1 - e_ * sin (ll)) / (1 + e_ * sin (ll)), e_ / 2.));
    err = fabs (ll - ll1);
    ll = ll1;
    i++;
  }
  if (err > 1e-9)
    return erc(GeoError::NONCONV, GeoErrors());

  *lat = ll;
  *lon = lam1 - atan ((S*cos (gamma0) - V * sin (gamma0)) / cos (B*u / A)) / B;
  return erc::success;
}
