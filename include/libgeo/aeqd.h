/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file AEQD.H Definition of \ref AzEqDist "Azimuthal Equidistant" projection

#include "projection.h"

/// Azimuthal Equidistant
class AzEqDist : public Projection
{
public:
  AzEqDist (const Params& par);
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;

  double h (double lat, double lon) const;
  double k (double lat, double lon) const;

  int Id () const { return (int)Projection::Method::AEQD; };

private:
  double sphi1, cphi1, n1, g;
  bool north_polar, south_polar;
};


