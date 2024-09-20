/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/

/*!
  \file STEREO.CPP Implementation of Stereographic 
  and \ref PolarStereo "Polar Stereographic" projections  
*/

#include <mlib/defs.h>
#include <libgeo/stereo.h>

using namespace mlib;
/*!
  \class Stereographic

  Formulas from:<blockquote>
    <b>A manual for geodetic coordinate transformations in the maritime provinces</b>,
    D.B.Thomson, E.J.Krakiwsky, R.R.Steeves, Technical Report No 48, July
    1977, Department of Geodesy and Geomatics Engineering, University of New
    Brunswick, Canada.</blockquote>
*/

/*!
  Constructor for %Sterographic projection. 

  Calculate parameters which are constant for this projection:
  - c1, c2 - constants used in transformation from ellipsoid to conformal sphere
  - chi0, lam0s - spherical lat/lon of origin of projection.
  - map_radius - radius of the projection sphere
  - r0 - map diameter (adjusted by scale factor)
***********************************************************************/
Stereographic::Stereographic (const Params& par) :
  Projection (par)
{
  c1 = sqrt (1 + ellipsoid().e2 () / (1 - ellipsoid().e2 ())*pow (cos (ref_latitude()), 4));
  chi0 = asin (sin (ref_latitude()) / c1);
  lam0s = ref_longitude() * c1;
  c2 = tan (M_PI_4 + chi0 / 2.) / pow (tan (M_PI_4 + ref_latitude() / 2.)*
    pow ((1 - ellipsoid().e ()*sin (ref_latitude())) / (1 + ellipsoid().e ()*sin (ref_latitude())), ellipsoid().e () / 2.), c1);
  double n0 = ellipsoid().a () / sqrt (1 - ellipsoid().e2 ()*sin (ref_latitude())*sin (ref_latitude()));
  double m0 = ellipsoid().a ()*(1 - ellipsoid().e2 ()) / pow (1 - ellipsoid().e2 ()*sin (ref_latitude())*sin (ref_latitude()), 1.5);

  // TODO - Do I need both map_radius and r0?
  map_radius = sqrt (n0*m0);
  r0 = 2 * k0() * map_radius;
}

erc Stereographic::GeoXY (double *x, double *y, double lat, double lon) const
{
  double e = ellipsoid().e ();
  double das = c2 * pow (tan (M_PI_4 + lat / 2.) * pow ((1 - e*sin (lat)) / (1 + e*sin (lat)), e / 2), c1);
  double chi = 2 * atan (das) - M_PI_2;
  double dlam = lon * c1 - lam0s;
  double denom = 1 + sin (chi0) * sin (chi) + cos (chi0) * cos (chi) * cos (dlam);
  *x = r0 * cos (chi) * sin (dlam) / denom;
  *y = r0 * (sin (chi)*cos (chi0) - cos (chi)*sin (chi0)*cos (dlam)) / denom;
  fwd_adj (*x, *y);
  return erc::success;
}

erc Stereographic::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  x /= k0();
  y /= k0();

  double ss = hypot (x, y);
  if (ss > 1e-20)
  {
    x /= ss;
    y /= ss;
  }
  else
  {
    x = 1.;
    y = 0.;
  }

  //calculate spherical coordinates
  double del = 2 * atan (ss / (2 * map_radius));
  double chi = asin (sin (chi0)*cos (del) + cos (chi0)*sin (del)*y);
  double lams = lam0s + asin (x*sin (del) / cos (chi));

  //ellipsoidal coordinates
  double phi = chi;
  double dif = 1;
  while (fabs (dif) > 1e-12)
  {
    double fact = pow ((1 - ellipsoid().e ()*sin (phi)) / (1 + ellipsoid().e ()*sin (phi)), ellipsoid().e () / 2);
    double tanf = tan (M_PI_4 + phi / 2);
    dif = c2 * pow (tanf*fact, c1) - tan (M_PI_4 + chi / 2);
    dif /= c1 * c2*pow (tanf*fact, c1 - 1)*(fact / 2 / pow (cos (M_PI_4 + phi / 2), 2) -
      ellipsoid().e2 ()*sin (phi)*cos (phi)*tanf / (1 - ellipsoid().e2 ()*sin (phi)));
    phi -= dif;
  }
  *lat = phi;
  *lon = lams / c1;
  return erc::success;
}

double Stereographic::k (double lat, double lon) const
{
  double scale;

  double dlam = lon * c1 - lam0s;
  double e = ellipsoid().e ();
  double das = c2 * pow (tan (M_PI_4 + lat / 2.) * pow ((1 - e*sin (lat)) / (1 + e*sin (lat)), e / 2), c1);
  double chi = 2 * atan (das) - M_PI_2;
  double denom = 1 + sin (chi0) * sin (chi) + cos (chi0) * cos (chi) * cos (dlam);
  double a_big = r0 / denom;

  scale = a_big * cos (chi) / (ellipsoid().a ()* ellipsoid().m (lat));
  return scale;
}

/*==============================================================================
                      Polar Stereographic projection 
==============================================================================*/

PolarStereo::PolarStereo (const Params& par) 
  : Projection (par)
{
  double k = k0 ();
  if (fabs (fabs (ref_latitude()) - M_PI / 2) > 1e-5)
  {
    //projection specified by latitude of standard parallel
    double tf;
    double sr = sin (ref_latitude());
    double x = pow ((1 + ellipsoid().e ()*sr) / (1 - ellipsoid().e ()*sr), ellipsoid().e () / 2);
    if (ref_latitude() < 0)
      tf = tan (M_PI / 4 + ref_latitude() / 2) / x;      //South polar case
    else
      tf = tan (M_PI / 4 - ref_latitude() / 2)*x;         //North polar case
    double mf = cos (ref_latitude()) / sqrt (1 - ellipsoid().e2 ()*sr*sr);
    k = mf * sqrt (pow (1 + ellipsoid().e (), 1 + ellipsoid().e ()) * pow (1 - ellipsoid().e (), 1 - ellipsoid().e ())) / (2 * tf);
  }
  rho1 = 2 * ellipsoid().a () * k / sqrt (pow (1 + ellipsoid().e (), 1 + ellipsoid().e ()) * pow (1 - ellipsoid().e (), 1 - ellipsoid().e ()));
  sc[0] = ellipsoid().e2 ()*(1 / 2. + ellipsoid().e2 ()*(5. / 24. + ellipsoid().e2 ()*(1. / 12. + 13.* ellipsoid().e2 () / 360.)));
  sc[1] = ellipsoid().e2 ()* ellipsoid().e2 ()*(7. / 48. + ellipsoid().e2 ()*(29. / 240. + ellipsoid().e2 ()*811. / 11520.));
  sc[2] = ellipsoid().e2 ()* ellipsoid().e2 ()* ellipsoid().e2 ()*(7. / 120. + ellipsoid().e2 ()*81. / 1120.);
  sc[3] = ellipsoid().e2 ()* ellipsoid().e2 ()* ellipsoid().e2 ()* ellipsoid().e2 ()*4279. / 161280.;
}

erc PolarStereo::GeoXY (double *x, double *y, double lat, double lon) const
{
  *y = (signbit(ref_latitude())?1:-1) * rho (lat) * cos (lon - ref_longitude());
  *x = rho (lat) * sin (lon - ref_longitude());
  fwd_adj (*x, *y);
  return erc::success;
}

erc PolarStereo::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  double rhoprime = hypot (x, y);
  double tprime = rhoprime / rho1;
  double chi;
  if (ref_latitude() < 0)
  {
    chi = 2 * atan (tprime) - M_PI / 2;
    if (x == 0)
      *lon = ref_longitude();
    else
      *lon = lon_adjust (ref_longitude() + atan2 (x, y));
  }
  else
  {
    chi = M_PI / 2 - 2 * atan (tprime);
    if (x == 0)
      *lon = ref_longitude();
    else
      *lon = lon_adjust (ref_longitude() + atan2 (x, -y));
  }
  *lat = chi + sc[0] * sin (2 * chi) + sc[1] * sin (4 * chi) + sc[2] * sin (6 * chi) + sc[3] * sin (8 * chi);
  return erc::success;
}

double PolarStereo::k (double lat, double lon) const
{
  return rho (lat) / (ellipsoid().a ()* ellipsoid().m (lat));
}

double PolarStereo::rho (double lat) const
{
  double sphi = sin (lat);
  double t;
  double e = ellipsoid().e ();
  if (ref_latitude() < 0)
    t = tan (M_PI / 4 + lat / 2) / pow (((1 + e*sphi) / (1 - e*sphi)), e / 2); //south polar
  else
    t = tan (M_PI / 4 - lat / 2) * pow (((1 + e*sphi) / (1 - e*sphi)), e / 2);

  return t * rho1;
}