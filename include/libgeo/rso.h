/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file RSO.H Rectified Skew Orthomorphic projection definition

#include "ome.h"

/// Rectified Skew Orthomorphic
class RSO : public ObliqueMercator
{
public:
  RSO (const Params& par) : ObliqueMercator (par) {};
  int Id() const { return (int)Projection::Method::RSO; };

protected:
  virtual void deskew( double u, double v, double *x, double *y ) const;
  virtual void skew( double *u, double *v, double x, double y ) const;
};
