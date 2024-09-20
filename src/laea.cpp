/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <mlib/defs.h>
#include <libgeo/laea.h>
#include <mlib/ipow.h>
#include <libgeo/errors.h>

/*!
  \file LAEA.CPP Implementation of
  \ref AzEqArea "Lambert Azimuthal Equal Area" projection
*/

using mlib::erc;

/*! 
  \class AzEqArea
  
  Formulas from Snyder pag. 187-190
*/

/*!
  Ctor for a Azimuthal Equal Area projection
*/
AzEqArea::AzEqArea (const Params& par) :
  Projection (par)
{
  q_p = ellipsoid ().q (M_PI / 2);
  r_q = ellipsoid ().a () * sqrt (q_p / 2);
  q_1 = ellipsoid ().q (ref_latitude ());
  m_1 = ellipsoid ().m (ref_latitude ());
  polar = (M_PI_2 - fabs (ref_latitude ()) < 1.E-6);
  if (!polar)
  {
    beta_1 = ellipsoid ().beta (ref_latitude ());
    d = m_1 / (sqrt (q_p / 2) * cos (beta_1));
  }
  ncval = 1 - (1 - ellipsoid ().e2 ()) * ellipsoid ().t (M_PI_2);
}

/*!
  Forward conversion from latitude/longitude to XY
*/
erc AzEqArea::GeoXY (double* x, double* y, double lat, double lon) const
{
  auto dlon = lon - ref_longitude ();
  if (!polar)
  {
    double blat = ellipsoid ().beta (lat);
    double b = r_q * M_SQRT2 / sqrt (1 + sin (beta_1) * sin (blat)
                                  + cos (beta_1) * cos (blat) * cos (dlon));
    *x = b * d * cos (blat) * sin (dlon);
    *y = b / d * (cos (beta_1) * sin (blat) - sin (beta_1) * cos (blat) * cos (dlon));
  }
  else
  {
    double rho = ellipsoid ().a () * sqrt (q_p - ellipsoid ().q (lat));
    *x = rho * sin (dlon);
    *y = ((ref_latitude () < 0) ? rho : -rho) * cos (dlon);
  }
  fwd_adj (*x, *y);
  return erc::success;
}

/*!
  Inverse conversion from XY to geographical coordinates.
*/
erc AzEqArea::XYGeo (double x, double y, double* lat, double* lon) const
{
  double rho, q;

  inv_adj (x, y);
  if (!polar)
  {
    rho = sqrt (x * x / (d * d) + d * d * y * y);
    double c_e = 2 * asin (rho / (2 * r_q));
    q = q_p * (cos (c_e) * sin (beta_1) + d * y * sin (c_e) * cos (beta_1) / rho);
    *lon = ref_longitude () +
      atan2 (x * sin (c_e), d * rho * cos (beta_1) * cos (c_e) - d * d * y * sin (beta_1) * sin (c_e));
  }
  else
  {
    rho = hypot (x, y);
    q = q_p - rho * rho / (ellipsoid ().a () * ellipsoid ().a ());
    if (ref_latitude () < 0)
    {
      q = -q;
      *lon = ref_longitude () + atan2 (x, y);
    }
    else
      *lon = ref_longitude () + atan2 (x, -y);
  }

  double phi;
  if (fabs ((fabs (q) - ncval)) < 1e-7)
    phi = (q < 0) ? -M_PI_2 : M_PI_2;
  else
  {
    double phi_new = asin (q / 2);
    int iter = 0;
    do
    {
      phi = phi_new;
      double sphi = sin (phi);
      phi_new = phi + mlib::ipow ((1 - ellipsoid ().e2 () * mlib::ipow (sphi, 3)), 2) / (2 * cos (phi)) *
        (q / (1 - ellipsoid ().e2 ()) - sphi / (1 - ellipsoid ().e2 () * sphi * sphi) + ellipsoid ().t (phi));
    } while (fabs (phi_new - phi) > 1e-7 && iter < MAX_ITER);
    *lat = phi;
    if (iter >= MAX_ITER)
      return erc (GeoError::NONCONV, GeoErrors ());
  }
  *lat = phi;
  return erc::success;
}

double AzEqArea::h (double lat, double lon) const
{
  return 1 / k (lat, lon);
}

double AzEqArea::k (double lat, double lon) const
{
  return sqrt (q_p-ellipsoid().q(lat))/ellipsoid().m(lat);
}


