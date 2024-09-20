/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of tests project and licensed under MIT license.
*/
#include <utpp/utpp.h>
#include <mlib/defs.h>
#include <libgeo/libgeo.h>

using mlib::erc;

/* Paramters for tests in Snyder p. 296, 297
  Clarke-1866 with 2 standard parallels*/

struct SnyderFixture
{
  SnyderFixture ();
  Lambert lcc;
  double lat_check, lon_check;
  double x_check, y_check, k_check;
};

SnyderFixture::SnyderFixture ()
  : lcc (Projection::Params(Ellipsoid::CLARKE_1866)
    .ref_longitude (-96_deg)
    .ref_latitude (23_deg)
    .south_latitude(33_deg)
    .north_latitude (45_deg)
    )
  , lat_check (35_deg)
  , lon_check (-75_deg)
  , x_check (1894410.9)
  , y_check (1564649.5)
  , k_check (0.9970171)
{
}

TEST_FIXTURE (SnyderFixture, LCCSynder_GeoXY)
{
  double x_result, y_result;
  int ret;
  ret = lcc.GeoXY (&x_result, &y_result, lat_check, lon_check);
  CHECK_EQUAL (0, ret);
  CHECK_CLOSE (x_check, x_result, 0.1);
  CHECK_CLOSE (y_check, y_result, 0.1);
}

TEST_FIXTURE (SnyderFixture, LCCSynder_Scale)
{
  double k_result, h_result;
  k_result = lcc.k (lat_check, lon_check);
  CHECK_CLOSE (k_check, k_result, 1e-7);
  h_result = lcc.h (lat_check, lon_check);
  CHECK_CLOSE (k_check, h_result, 1e-7);
}


TEST_FIXTURE (SnyderFixture, LCCSynder_XYGeo)
{
  double lat_result, lon_result;
  int ret;
  ret = lcc.XYGeo (x_check, y_check, &lat_result, &lon_result);
  CHECK_EQUAL (0, ret);
  CHECK_CLOSE (lat_check, lat_result, 1e-7);
  CHECK_CLOSE (lon_check, lon_result, 1e-7);
}

//from EPSG Guidance Note 7 page 27
TEST (Michigan)
{
  Lambert proj (Projection::Params (Ellipsoid::CLARKE_1866)
    .ref_latitude (DM (43, 19))
    .ref_longitude (-DM (84, 20))
    .north_latitude (DM (44, 11))
    .south_latitude (DM (45, 42))
    .unit (1_ftUS)
    .k0 (1.0000382)
    .false_east (2'000'000));

  constexpr double lat_in = DM (43, 45), lon_in = -DM (83, 10),
    x_in = 2308335.75, y_in = 160210.48;
  double lat_out, lon_out, x_out, y_out;

  CHECK_EQUAL (erc::success, proj.GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.005);
  CHECK_CLOSE (y_out, y_in, 0.005);

  CHECK_EQUAL (mlib::erc::success, proj.XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1e-7);
  CHECK_CLOSE (lon_in, lon_out, 1e-7);
}

/* 
  Example from EPSG Guidance Note 7 page 26 but with reference longitude modified
  to account for the 29.2985" arc-seconds.
*/
TEST (Belge_Lambert_72)
{
  Lambert proj (Projection::Params (Ellipsoid::INTERNATIONAL)
    .ref_latitude (90_deg)
    .ref_longitude (DMS (4, 22, 02.952))
    .north_latitude (DM (51, 10))
    .south_latitude (DM (49, 50))
    .false_east (150000.01)
    .false_north (5'400'088.44));

  constexpr double lat_in = DMS (50, 40, 46.461), lon_in = DMS (5, 48, 26.533),
    x_in = 251763.20, y_in = 153034.13;
  double lat_out, lon_out, x_out, y_out;

  CHECK_EQUAL (erc::success, proj.GeoXY (&x_out, &y_out, lat_in, lon_in));
  CHECK_CLOSE (x_in, x_out, 0.005);
  CHECK_CLOSE (y_out, y_in, 0.005);

  CHECK_EQUAL (mlib::erc::success, proj.XYGeo (x_in, y_in, &lat_out, &lon_out));
  CHECK_CLOSE (lat_in, lat_out, 1e-7);
  CHECK_CLOSE (lon_in, lon_out, 1e-7);

}
