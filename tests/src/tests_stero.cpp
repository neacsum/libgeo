/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of tests project and licensed under MIT license.
*/
#include <utpp/utpp.h>
#include <mlib/defs.h>
#include <libgeo/libgeo.h>

/*
  Source:
    Coördinaattransformaties en kaartprojecties
    Formules en parameters
    3e herziene uitgave December 2000

    http://kadaster.nl
*/

struct Amersfoort : public Stereographic {
  Amersfoort () : Stereographic (Projection::Params (Ellipsoid::BESSEL_1841)
    .k0 (0.9999079)
    .ref_longitude (DMS (5, 23, 15.5))
    .ref_latitude (DMS (52, 9, 22.178))
    .false_east (155000)
    .false_north (463000)
  ) {}
};

TEST_FIXTURE (Amersfoort, Netherlands1)
{
  constexpr double 
    lat_in = DMS(52, 12, 34.567),
    lon_in = DMS(4, 23, 45.678),
    x_in = 87232.211,
    y_in = 469408.512;
  double lat_out, lon_out, x_out, y_out;

  CHECK_EQUAL (0, GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.001);
  CHECK_CLOSE (y_in, y_out, 0.001);

  CHECK_EQUAL (0, XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1E-9);
  CHECK_CLOSE (lon_in, lon_out, 1E-9);
}

/*
  EPSG Guidance note 7, page 77
*/
TEST_FIXTURE (Amersfoort, Nertherlands2)
{
  constexpr double
    lat_in = 53_deg,
    lon_in = 6_deg,
    x_in = 196105.283,
    y_in = 557057.739;
  double lat_out, lon_out, x_out, y_out;

  CHECK_EQUAL (0, GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.001);
  CHECK_CLOSE (y_in, y_out, 0.001);

  CHECK_EQUAL (0, XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1E-9);
  CHECK_CLOSE (lon_in, lon_out, 1E-9);
}
/*
  Source:
    A manual for geodetic coordinate transformations in the maritime provinces, 
    D.B. Thomson, E.J. Krakiwsky, R.R. Steeves, Technical Report No. 48, July 
    1977, Department of Geodesy and Geomatics Engineering, University of New 
    Brunswick, Canada.
    http://gge.unb.ca/Pubs/TR48.pdf
*/
TEST (NewBrunswick)
{
  Stereographic ost (
    Projection::Params (Ellipsoid::CLARKE_1866)
    .k0 (0.999912)
    .ref_longitude (-DM(66, 30))
    .ref_latitude (DM(46,30))
    .false_east (300000)
    .false_north (800000));

  const double 
    lat_ref = 47'03'24.644_dms,
    lon_ref = -65'29'03.453_dms,
    x_ref = 377164.887,
    y_ref = 862395.774;

  double lat, lon, x, y;

  //Forward transformation
  CHECK_EQUAL (0, ost.GeoXY (&x, &y, lat_ref, lon_ref));
  CHECK_CLOSE (x_ref, x, 0.001);
  CHECK_CLOSE (y_ref, y, 0.001);

  //Inverse transformation
  CHECK_EQUAL (0, ost.XYGeo (x_ref, y_ref, &lat, &lon));
  CHECK_CLOSE (lat_ref, lat, 1E-9);
  CHECK_CLOSE (lon_ref, lon, 1E-9);
}

/*
  Polar aspect with known k0. Source: Snyder p. 314, 317

  Synder gives x = -1573645.4. Manual calculation gives -1573645.25. The difference
  appears in the value of rho1. It is probably just a rounding error.
*/
TEST(PolarStereographic_k0)
{
  PolarStereo pst (
    Projection::Params (Ellipsoid::INTERNATIONAL)
    .k0 (0.994)
    .ref_longitude (-DEG (100))
    .ref_latitude (-DEG (90))
  );

  const double 
    lat_ref = -DEG(75),
    lon_ref = DEG(150),
    x_ref = -1573645.25,
    y_ref = -572760.1;

  double lat, lon, x, y;

  //Forward transformation
  CHECK_EQUAL (0, pst.GeoXY (&x, &y, lat_ref, lon_ref));
  CHECK_CLOSE (x_ref, x, 0.1);
  CHECK_CLOSE (y_ref, y, 0.1);

  //Inverse transformation
  CHECK_EQUAL (0, pst.XYGeo (x_ref, y_ref, &lat, &lon));
  CHECK_CLOSE (lat_ref, lat, 1E-7);
  CHECK_CLOSE (lon_ref, lon, 1E-7);
}

/* 
  Polar aspect with known phi_c not at the pole. Source Snyder p. 315
*/
TEST(PolarStereographic_phic)
{
  PolarStereo pst (
    Projection::Params (Ellipsoid::INTERNATIONAL)
    .ref_longitude (-DEG (100))
    .ref_latitude (-DEG (71))
  );

  const double 
    lat_ref = -DEG(75),
    lon_ref = DEG(150),
    x_ref = -1540033.6,
    y_ref = -560526.4,
    k_ref = 0.9896256;

  double lat, lon, x, y;

  //Forward transformation
  CHECK_EQUAL (0, pst.GeoXY (&x, &y, lat_ref, lon_ref));
  CHECK_CLOSE (x_ref, x, 0.1);
  CHECK_CLOSE (y_ref, y, 0.1);

  //Inverse transformation
  CHECK_EQUAL (0, pst.XYGeo (x_ref, y_ref, &lat, &lon));
  CHECK_CLOSE (lat_ref, lat, 1E-7);
  CHECK_CLOSE (lon_ref, lon, 1E-7);

  //Scale factor
  CHECK_CLOSE (k_ref, pst.k(lat_ref, lon_ref), 1e-7);
}

