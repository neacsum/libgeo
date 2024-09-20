/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/

///  \file TME.CPP Implementation of TransverseMercator projection



#include <libgeo/tme.h>
#include <libgeo/errors.h>
#include <libgeo/units.h>

using namespace mlib;

double pd[12] = { -3./2.,  9./16., 15./16., -15./32., -35./48., 315./512.,
                   3./2.,-27./32., 21./16., -55./32., 151./96., 1097./512.};

TransverseMercator::TransverseMercator (const Params& par) 
  : Projection (par)
{
  double a2, a4, a6, a8, c2, c4, c6, c8, omegaf;

  n = ellipsoid().f () / (2. - ellipsoid().f ());
  a2 = pd[0] * n + pd[1] * pow (n, 3);
  a4 = pd[2] * n * n + pd[3] * pow (n, 4);
  a6 = pd[4] * pow (n, 3);
  a8 = pd[5] * pow (n, 4);

  b0 = 2. * (a2 - 2 * a4 + 3 * a6 - 4. * a8);
  b2 = 8. * (a4 - 4. * a6 + 10. * a8);
  b4 = 32. * (a6 - 6. * a8);
  b6 = 128. * a8;


  r = ellipsoid().a () * (1. - n) * (1. - n * n) * (1. + 2.25 * n * n + 225. * pow (n, 4) / 64.);
  if (k0 () * r == 0.)
    erc (GeoError::PARAM, GeoErrors()).raise();
  omegaf = ref_latitude() + sin (ref_latitude()) * cos (ref_latitude()) *
    (b0 + b2 * pow (cos (ref_latitude()), 2) + b4 * pow (cos (ref_latitude()), 4) +
      b6 * pow (cos (ref_latitude()), 6));
  yvalue = -k0 () * omegaf * r;
  e2 = ellipsoid().f () * (2. - ellipsoid().f ());
  ep2 = e2 / (1. - e2);

  c2 = pd[6] * n + pd[7] * pow (n, 3);
  c4 = pd[8] * n * n + pd[9] * pow (n, 4);
  c6 = pd[10] * pow (n, 3);
  c8 = pd[11] * pow (n, 4);
  d0 = 2 * (c2 - 2 * c4 + 3 * c6 - 4 * c8);
  d2 = 8 * (c4 - 4 * c6 + 10 * c8);
  d4 = 32 * (c6 - 6 * c8);
  d6 = 128 * c8;
}

erc TransverseMercator::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  double phifp, rnfp, nufp2, tfp, tfp2, g2, g3, g4, g5, g6, g7, q, l;
  double omega = (y - yvalue) / (k0 () * r);
  if (fabs (omega) >= M_PI_2)
    return erc (GeoError::RANGE, GeoErrors());
  phifp = omega + sin (omega) * cos (omega) * (d0 + d2 * cos (omega) * cos (omega) +
    d4 * pow (cos (omega), 4) + d6 * pow (cos (omega), 6));
  if (fabs (phifp) >= M_PI_2)
    return erc (GeoError::RANGE, GeoErrors());
  rnfp = ellipsoid().a () / sqrt (1 - e2 * sin (phifp) * sin (phifp));
  nufp2 = ep2 * cos (phifp) * cos (phifp);
  tfp = tan (phifp);
  tfp2 = tfp * tfp;
  g2 = -tfp * (1 + nufp2) / 2;
  g3 = -(1 + 2 * tfp2 + nufp2) / 6;
  g4 = -(5 + 3 * tfp2 + nufp2 * (1 - 9 * tfp2) - 4 * nufp2 * nufp2) / 12;
  g5 = (5 + 28 * tfp2 + 24 * tfp2 * tfp2 + nufp2 * (6 + 8 * tfp2)) / 120;
  g6 = (61 + 90 * tfp2 + 45 * tfp2 * tfp2 + nufp2 *
    (46 - 252 * tfp2 - 90 * tfp2 * tfp2)) / 360;
  g7 = (61 + 662 * tfp2 + 1320 * tfp2 * tfp2 + 720 * pow (tfp, 6)) / 5040;
  q = x / (k0 () * rnfp);
  *lat = phifp + g2 * q * q * (1 + q * q * (g4 + g6 * q * q));
  if (*lat < -M_PI_2 || *lat > M_PI_2)
    return erc (GeoError::RANGE, GeoErrors());
  l = q * (1 + q * q * (g3 + q * q * (g5 + g7 * q * q)));
  *lon = lon_adjust (ref_longitude() + l / cos (phifp));
  return erc::success;
}

erc TransverseMercator::GeoXY (double *x, double *y, double lat, double lon) const
{
  double omega;
  double clat, nu2, l, tlat2, e3, e4, e5, e6, e7, xtrue, ytrue, rn;
  clat = cos (lat);
  nu2 = ep2 * clat * clat;
  if (-M_PI <= (lon - ref_longitude()) && (lon - ref_longitude()) <= M_PI)
    l = (lon - ref_longitude()) * clat;
  else if ((lon - ref_longitude()) > M_PI)
    l = (lon - ref_longitude() - 2 * M_PI) * clat;
  else
    l = (lon - ref_longitude() + 2 * M_PI) * clat;
  omega = lat + sin (lat) * clat * (b0 + b2 * clat * clat + b4 * pow (clat, 4) + b6 * pow (clat, 6));
  if (lat == M_PI / 2)
    return erc (GeoError::RANGE, GeoErrors());
  tlat2 = tan (lat);
  tlat2 *= tlat2;
  e3 = (1. - tlat2 + nu2) / 6.;
  e4 = (5. - tlat2 + nu2 * (9. + 4. * nu2)) / 12.;
  e5 = (5. - 18. * tlat2 + tlat2 * tlat2 + nu2 * (14. - 58. * tlat2)) / 120.;
  e6 = (61. - 58. * tlat2 + tlat2 * tlat2 + nu2 * (270. - 330. * tlat2)) / 360.;
  e7 = (61. - 479. * tlat2 + 179. * tlat2 * tlat2 - pow (tlat2, 3)) / 5040.;
  rn = ellipsoid().a () / sqrt (1 - e2 * pow (sin (lat), 2));
  xtrue = k0 () * rn * l * (1 + l * l * (e3 + l * l * (e5 + e7 * l * l)));
  ytrue = k0 () * (omega * r + rn * tan (lat) * l * l * (1. + l * l * (e4 + e6 * l * l)) / 2.);
  *x = xtrue;
  *y = ytrue + yvalue;

  fwd_adj (*x, *y);
  return erc::success;
}

double TransverseMercator::k (double lat, double lon) const
{
  double scale;
  double clat = cos (lat);
  double a_big2 = (lon - ref_longitude())*clat; a_big2 *= a_big2;
  double c_big = ep2 * clat*clat;
  double t_big = tan (lat); t_big *= t_big;  //tan(lat)^2

  double c[3] = {
    (1 + c_big) / 2,
    (5 - 4 * t_big + 42 * c_big + 13 * c_big*c_big - 28 * ep2) / 24,
    (61 - 148 * t_big + 16 * t_big*t_big) / 720
  };

  scale = 1 + (((c[2] * a_big2) + c[1])*a_big2 + c[0])*a_big2;
  scale *= k0 ();

  return scale;
}

double TransverseMercator::Convergence (double lat, double lon) const
{
  return atan (-sin (lat)*tan (lon - ref_longitude()));
}

// UTM Projection ------------------------------------------------------------
UTM::UTM (const Ellipsoid& ell, int zone, bool south)
  : TransverseMercator (
    Params{ ell }.k0 (0.9996).false_east (500000)
    .false_north (south ? 10000000. : 0.).ref_longitude ((zone * 6 - 183) * D2R))
{
  if (zone < 1 || zone >= 60)
    throw erc (GeoError::PARAM, GeoErrors ());
}
