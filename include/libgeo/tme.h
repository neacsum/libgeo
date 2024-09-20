/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

/*!  
    \file TME.H  Definition of \ref TransverseMercator "Transverse Mercator"
    and UTM projections
*/
#include "projection.h"

/// Transverse %Mercator projection
class TransverseMercator : public Projection
{
public:
  TransverseMercator (const Params& pp);
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;
  double k (double lat, double lon) const;
  const char *Name () const { return "TME"; };
  int Id () const { return (int)Projection::Method::TME; };
  double Convergence (double lat, double lon) const;
private:
  double n, e2, ep2;
  double d0, d2, d4, d6;
  double b0, b2, b4, b6;
  double r;
  double yvalue;
};



class UTM : public TransverseMercator
{
public:
  UTM (const Ellipsoid& ell, int zone, bool south = false);

};
