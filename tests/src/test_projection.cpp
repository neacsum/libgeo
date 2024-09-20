/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of tests project and licensed under MIT license.
*/
#include <utpp/utpp.h>
#include <mlib/defs.h>
#include <libgeo/libgeo.h>


TEST (ProjectionParameters)
{
  Projection::Params par;
  par.k0 (0.9876)
    .ref_latitude (45_deg)
    .ref_longitude (-75_deg)
    .north_latitude (50_deg)
    .south_latitude (40_deg)
    .skew_azimuth (123_deg)
    .false_east (1234567)
    .false_north (7654321)
    .unit(.3048);

  PlateCarree proj(par);
  CHECK_EQUAL ((int)Projection::Method::EQC, proj.Id ());
  CHECK_EQUAL (A_WGS84, proj.ellipsoid().a ());
  CHECK_CLOSE (F_WGS84, proj.ellipsoid ().f (), 1e-8);
  CHECK_EQUAL (0.3048, proj.unit ());
  CHECK_EQUAL (0.9876, proj.k0 ());
  CHECK_EQUAL (-75_deg, proj.ref_longitude ());
  CHECK_EQUAL (45_deg, proj.ref_latitude ());
  CHECK_EQUAL (50_deg, proj.north_parallel ());
  CHECK_EQUAL (40_deg, proj.south_parallel ());
  CHECK_EQUAL (123_deg, proj.skew_azimuth ());
  CHECK_EQUAL (1234567, proj.false_east ());
  CHECK_EQUAL (7654321, proj.false_north ());
}

/* This function was used for many years. */
double OldLonAdjust( double lon )
{
	int sign = (lon>=0)?1:-1;

	while ( fabs(lon) > M_PI )
		lon -= 2*M_PI * sign;
	return lon;
}

TEST (Projection_LongitudeAdjustment)
{
  struct TestProjection : PlateCarree
  {
  public:
    TestProjection () : PlateCarree (Ellipsoid::WGS_84) {};
    Projection::lon_adjust;
  };

  TestProjection t;

  CHECK_CLOSE (0.5, t.lon_adjust (0.5+4*M_PI), 1E-15);

  double fval = 0.5+4*M_PI;

  CHECK_CLOSE (OldLonAdjust(fval), t.lon_adjust (fval), 1e-15);
  CHECK_CLOSE (OldLonAdjust(-fval), t.lon_adjust (-fval), 1e-15);

}