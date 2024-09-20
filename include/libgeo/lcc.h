/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file LCC.H Definition of \ref Lambert "Lambert Conformal Conical" projection

#include "projection.h"

/// %Lambert Conformal Conical
class Lambert : public Projection
{
public:
  Lambert (const Params& pp);
  int Id () const { return (int)Projection::Method::LCC; };

  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;
  double h (double lat, double lon) const { return k (lat, lon); };
  double k (double lat, double lon) const;

  double Convergence (double lat, double lon) const;

  inline static int MAX_ITER = 10;     ///< max number of iterations in reverse formulas

private:
  double n, af_big, rho0;
  double tfunc (double phi) const;
};

