/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

/// \file OME.H Definition of \ref ObliqueMercator "Oblique Mercator" projection

#include "projection.h"

///Oblique %Mercator
class ObliqueMercator : public Projection
{
public:
  ObliqueMercator (const Params& par);
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;

  double h (double lat, double lon) const { return k (lat, lon); };
  double k (double lat, double lon) const;

protected:
  virtual void deskew (double u, double v, double *x, double *y) const = 0;
  virtual void skew (double *u, double *v, double x, double y) const = 0;
  double exptau (double val) const;
  double gamma0;

private:
  mlib::erc uvgeo (double u, double v, double *lat, double *lon) const;
  double A, B, E, lam1;
};

