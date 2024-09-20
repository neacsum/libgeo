/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of tests project and licensed under MIT license.
*/
#include <utpp/utpp.h>
#include <Windows.h>
#include <fstream>

#include <libgeo/libgeo.h>
#include <libgeo/units.h>

using mlib::erc;

TEST (Lambert_azimuthal_equal_area_oblique)
{
  AzEqArea proj (
    Projection::Params (Ellipsoid::CLARKE_1866)
      .ref_latitude (40_deg)
      .ref_longitude (-100_deg));
  constexpr double lat_in = 30_deg,
    lon_in = -110_deg,
    x_in = -965'932.1,
    y_in = -1'056'814.9;

  double lat_out, lon_out, x_out, y_out;
  CHECK_EQUAL (erc::success, proj.GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.1);
  CHECK_CLOSE (y_out, y_in, 0.1);

  CHECK_EQUAL (mlib::erc::success, proj.XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1e-7);
  CHECK_CLOSE (lon_in, lon_out, 1e-7);
}

TEST (Lambert_azimuthal_equal_area_polar)
{
  AzEqArea proj (
    Projection::Params (Ellipsoid::INTERNATIONAL)
      .ref_latitude (90_deg)
      .ref_longitude (-100_deg)
  );

  constexpr double lat_in = 80_deg,
    lon_in = 5_deg,
    x_in = 1'077'459.7,
    y_in = 288'704.5,
    k_in = 1.0038193,
    h_in = 0.9961952;

  double lat_out, lon_out, x_out, y_out, k_out, h_out;
  CHECK_EQUAL (mlib::erc::success, proj.GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.1);
  CHECK_CLOSE (y_out, y_in, 0.1);

  k_out = proj.k (lat_in, lon_in);
  h_out = proj.h (lat_in, lon_in);
  CHECK_CLOSE (k_in, k_out, 1e-6);
  CHECK_CLOSE (h_in, h_out, 1e-6);

  CHECK_EQUAL (erc::success, proj.XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1e-7);
  CHECK_CLOSE (lon_in, lon_out, 1e-7);
}

TEST (Azimuthal_equidistant_polar)
{
  AzEqDist proj (
    Projection::Params (Ellipsoid::INTERNATIONAL)
      .ref_latitude (90_deg)
      .ref_longitude (-100_deg));

  constexpr double lat_in = 80_deg,
    lon_in = 5_deg,
    x_in = 1'078'828.3,
    y_in = 289'071.2,
    k_in = 1.0050946;

  double lat_out, lon_out, x_out, y_out, k_out;
  CHECK_EQUAL (mlib::erc::success, proj.GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.1);
  CHECK_CLOSE (y_out, y_in, 0.1);

  k_out = proj.k (lat_in, lon_in);
  CHECK_CLOSE (k_in, k_out, 1e-6);

  CHECK_EQUAL (erc::success, proj.XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1e-7);
  CHECK_CLOSE (lon_in, lon_out, 1e-7);
}

TEST (Cassini_Trinidad)
{
  Cassini proj (Projection::Params (Ellipsoid::CLARKE_1858)
    .ref_latitude (DMS (10, 26, 30))
    .ref_longitude (-DM (61, 20))
    .unit (1_lkCla)
    .false_east (430000)
    .false_north (325000));

  constexpr double lat_in = 10_deg,
    lon_in = -62_deg,
    x_in = 66644.94,
    y_in = 82536.22;

  double lat_out, lon_out, x_out, y_out;
  CHECK_EQUAL (mlib::erc::success, proj.GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.1);
  CHECK_CLOSE (y_out, y_in, 0.1);

  CHECK_EQUAL (erc::success, proj.XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1e-7);
  CHECK_CLOSE (lon_in, lon_out, 1e-7);
}

TEST (unit_operators)
{
  double val = 123456.78;
  double v1 = DMS2deg (val);

  CHECK_CLOSE (12 + 34 / 60. + 56.78 / 3600, v1, 1e-14);
}

TEST_MAIN (int argc, char **argv)
{
  SetConsoleOutputCP (65001);
  if (argc == 1)
    return UnitTest::RunAllTests ();
  else
    return UnitTest::RunSuite (argv[1]);
}

