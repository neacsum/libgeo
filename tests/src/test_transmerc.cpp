/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of tests project and licensed under MIT license.
*/
#include <utpp/utpp.h>
#include <mlib/defs.h>
#include <libgeo/libgeo.h>

//Snyder pag. 269-270
TEST (TransverseMercator_Forward)
{
  double lat = 40.5 * D2R;
  double lon = -73.5 * D2R;
  double x, y, k;
  int ret;

  TransverseMercator tm(Projection::Params(Ellipsoid::CLARKE_1866)
    .k0(0.9996)
    .ref_longitude(-75_deg)
  );
  ret = tm.GeoXY (&x, &y, lat, lon);
  CHECK_EQUAL (0, ret);
  CHECK_CLOSE (127106.5, x, 0.1);
  CHECK_CLOSE (4484124.4, y, 0.1);

  k = tm.k(lat, lon);
  CHECK_CLOSE (0.9997989, k, 1e-7);

}

TEST (TransverseMercator_Inverse)
{
  double lat, lon;
  double x=127106.5, y=4484124.4;
  int ret;

  TransverseMercator tm (Projection::Params (Ellipsoid::CLARKE_1866)
    .k0 (0.9996)
    .ref_longitude (-75_deg)
  );
  ret = tm.XYGeo (x, y, &lat, &lon);
  CHECK_EQUAL (0, ret);
  CHECK_CLOSE (40.5, lat/D2R, 1e-6);
  CHECK_CLOSE (-73.5,lon/D2R, 1e-6);
}

//from old TEST.INI
TEST (TransverseMercator_own)
{
  double lat=45*D2R, 
    lon = -72*D2R,
    x = 423058.45, 
    y = 4985540.61;

  double latr, lonr, xr, yr;

  TransverseMercator tm(Projection::Params()
    .k0(0.9999)
    .ref_longitude (-73.5_deg)
    .false_east(304800)
   );
  tm.GeoXY (&xr, &yr, lat, lon);
  CHECK_CLOSE (x, xr, 0.01);
  CHECK_CLOSE (y, yr, 0.01);

  tm.XYGeo (x, y, &latr, &lonr);
  CHECK_CLOSE (latr, lat, 1e-6);
  CHECK_CLOSE (lonr, lon, 1e-6);

}

TEST (Transverse_Mercator_Convergence)
{
  TransverseMercator tm (
    Projection::Params(Ellipsoid::WGS_84)
    .k0(0.999933333)
    .ref_longitude(-DM(115, 45))
    .ref_latitude(DM(41, 40))
    .false_east(800000)
   );

  constexpr double lat=DMS(46, 20, 19.8681), lon = -DMS(117, 1, 50.8077),
    x = 701400.887, y = 519900.099, gamma = 0.926635181_deg;
  double latr, lonr, xr, yr, gammar;

  tm.GeoXY (&xr, &yr, lat, lon);
  CHECK_CLOSE (x, xr, 0.01);
  CHECK_CLOSE (y, yr, 0.01);

  tm.XYGeo (x, y, &latr, &lonr);
  CHECK_CLOSE (latr, lat, 1e-6);
  CHECK_CLOSE (lonr, lon, 1e-6);

  gammar = tm.Convergence (lat, lon);
  CHECK_CLOSE (gamma, gammar, 1e-10);
}

//From https://ngs.nooa.gov/NCAT
TEST (NOAA_NCAT)
{
  UTM utm{ Ellipsoid::WGS_84, 14 }; //UTM zone 14N

  double lat = DMS(39, 13, 26.7120),
    lon = -DMS(98, 32, 31.7454),
    x = 539520.632,
    y = 4341744.059,
    gamma = -DMS (0, 17, 22.30);
  double latr, lonr, xr, yr, gammar;

  utm.GeoXY (&xr, &yr, lat, lon);
  CHECK_CLOSE (x, xr, 0.01);
  CHECK_CLOSE (y, yr, 0.01);

  utm.XYGeo (x, y, &latr, &lonr);
  CHECK_CLOSE (latr, lat, 1e-6);
  CHECK_CLOSE (lonr, lon, 1e-6);

  gammar = utm.Convergence (lat, lon);
  CHECK_CLOSE (gamma, gammar, 1e-7);

}