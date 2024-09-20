/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file HOM.H  Definition of \ref Hotine "Hotine Oblique Mercator" projection

#include "ome.h"

/// %Hotine Oblique %Mercator
class Hotine : public ObliqueMercator
{
public:
  Hotine (const Params& par) : ObliqueMercator (par) {}
  int Id () const {return (int)Projection::Method::HOM;}

protected:
  virtual void deskew (double u, double v, double *x, double *y) const;
  virtual void skew (double *u, double *v, double x, double y) const;
};
