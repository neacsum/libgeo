/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/

/// \file MER.CPP Implementation of Mercator projection

#include <libgeo/mer.h>
#include <libgeo/errors.h>
using mlib::erc;

static double sqrarg__;
#define SQR(a) ((sqrarg__=(a)) == 0.?0.:sqrarg__*sqrarg__)

Mercator::Mercator (const Params& par) :
  Projection (par)
{
  if (ref_latitude() == M_PI / 2)
    erc (GeoError::PARAM, GeoErrors()).raise();
  sfeq = cos (ref_latitude())*ellipsoid().rn (ref_latitude());
}

erc Mercator::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  double t = exp (-y / sfeq);
  double phi = M_PI_2 - 2 * atan (t);
  double phidif = 1.;
  double e = ellipsoid().e ();
  while (phidif > .5e-10)
  {
    double sphi;
    double p5 = M_PI_2 - 2. * atan (t * pow ((1. - e * (sphi = sin (phi))) / (1. + e * sphi), e / 2.));
    phidif = fabs (phi - p5);
    phi = p5;
  }
  *lon = lon_adjust (ref_longitude() + x / sfeq);
  *lat = phi;
  return erc::success;
}

erc Mercator::GeoXY (double *x, double *y, double lat, double lon) const
{
  double tau, slat;
  if (fabs (lat) == M_PI_2)
    return erc (GeoError::RANGE, GeoErrors());
  lon = lon_adjust (lon - ref_longitude());
  *x = lon * sfeq;
  slat = sin (lat);
  double e = ellipsoid().e ();
  tau = sfeq / 2.*log ((1 + slat) / (1 - slat) * pow ((1 - e * slat) / (1 + e * slat), e));
  *y = tau;
  fwd_adj (*x, *y);
  return erc::success;
}

double Mercator::k (double lat, double lon) const
{
  return sqrt ((1 - ellipsoid().e2 ()*SQR (sin (lat))) / (1 - ellipsoid().e2 ()*SQR (sin (ref_latitude()))))*(cos (ref_latitude()) / cos (lat));
}

