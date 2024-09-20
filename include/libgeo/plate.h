/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file PLATE.H  Definition of \ref PlateCarree "Plate Carée" projection

#include "projection.h"

///Plate Carree projection
class PlateCarree : public Projection
{
public:
  PlateCarree (const Params& pp);
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;

  double h (double lat, double lon) const { return 1.; };
  double k (double lat, double lon) const;

  int Id () const { return (int)Projection::Method::EQC; };
};

