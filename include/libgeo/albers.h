/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file ALBERS.H Definition of \ref Albers "Albers Equal Area" projection

#include "projection.h"

/// %Albers Equal Area
class Albers : public Projection
{
public:
  Albers (const Params& pp);
  int Id () const { return (int)Projection::Method::AEA; };

  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;

  double h (double lat, double lon) const;
  double k (double lat, double lon) const;

  inline static int MAX_ITER = 10;     ///< max number of iterations in reverse formulas

private:
  double npar_, spar_, reflat_, reflon_, n, c_big, rho0;
  using Projection::par_;
};

