/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

/*! 
  \file OCY.H Definition of \ref ObliqueCylindrical "Oblique Cylindrical"
  projection 
*/

#include "projection.h"

/// Oblique Cylindrical
class ObliqueCylindrical : public Projection
{
public:
  ObliqueCylindrical (const Params& pp);
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;

  double h (double lat, double lon) const { return k (lat, lon); };
  double k (double lat, double lon) const;

  int Id () const { return (int)Projection::Method::OCY; };

private:
  double map_radius, c1, c2, chi0;
};

