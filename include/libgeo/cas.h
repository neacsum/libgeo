/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file CAS.H Cassini projection definition

#include "projection.h"

/// Cassini-Soldner projection
class Cassini : public Projection
{
public:
  Cassini (const Params& pp);
  int Id () const { return (int)Projection::Method::CASS; };

  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;
  double h (double lat, double lon) const;
  double k (double lat, double lon) const { return 1.; };

private:
  double s0;
};

