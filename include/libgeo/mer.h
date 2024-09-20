/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file MER.H  Definition of Mercator projection

#include "projection.h"

/// %Mercator projection
class Mercator : public Projection
{
public:
  Mercator (const Params& pp);
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;

  double h (double lat, double lon) const { return k (lat, lon); };
  double k (double lat, double lon) const;

  int Id () const { return (int)Projection::Method::MER; };

private:
  double sfeq;    //radius of circle of parallel
};


