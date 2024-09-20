/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file POL.H Definition of Polyconic projection

#include "projection.h"

/// %Polyconic
class Polyconic : public Projection
{
public:
  Polyconic (const Params& pp);
  int Id () const { return (int)Projection::Method::POL; };

  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;
  double h (double lat, double lon) const;

  inline static int MAX_ITER = 10;     ///< max number of iterations in reverse formulas

private:
  double m1 (double lat) const;
  double m0;
  double c[4];
};

