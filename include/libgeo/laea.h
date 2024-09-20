/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file LAEA.H Definition of Lambert Azimuthal Equal Area projection

#include "projection.h"

/// Azimuthal Equal Area
class AzEqArea : public Projection
{
public:
  AzEqArea (const Params& pp);
  int Id () const { return (int)Projection::Method::LAEA; };

  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;

  double h (double lat, double lon) const;
  double k (double lat, double lon) const;

  inline static int MAX_ITER = 10;     ///< max number of iterations in reverse formulas

private:
  double r_q, q_p, q_1, m_1, beta_1, d;
  double ncval;
  bool polar;
};

