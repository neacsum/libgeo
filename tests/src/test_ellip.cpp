/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of tests project and licensed under MIT license.
*/
#include <utpp/utpp.h>
#include <libgeo/libgeo.h>

Ellipsoid wgs84;
Ellipsoid clarke66 (6378206.4, 1/294.978698200);

SUITE (Ellipsoid)
{
  TEST (SemiminorAxis)
  {
    CHECK_CLOSE (6356752.3142, wgs84.b (), 0.001);
  }

  TEST (FirstEccentricity)
  {
    CHECK_CLOSE (0.08181919084262, wgs84.e (), 1E-14);
  }

  TEST (FirstEccentricitySquared)
  {
    CHECK_CLOSE (0.00669437999014, wgs84.e2 (), 1E-14);
  }

  TEST (qAux)
  {
    CHECK_CLOSE (1.2792602, clarke66.q (40 * D2R), 1e-7);
  }

  TEST (AuthaliticLatitude)
  {
    double q = clarke66.q (40 * D2R);
    double qp = clarke66.q (M_PI / 2);

    CHECK_CLOSE (1.9954814, qp, 1e-7);
    CHECK_CLOSE (39.8722878, (clarke66.beta (40 * D2R)) / D2R, 1e-7);
  }
}


