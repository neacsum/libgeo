/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

/*!  
  \file STEREO.H Definition of Sterographic 
  and \ref PolarStereo "Polar Stereographic" projections  
*/
#include "projection.h"

/// %Stereographic projection
class Stereographic : public Projection
{
public:
  Stereographic (const Params& pp);
  int Id () const { return (int)Projection::Method::STEREA; };

  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  double k (double lat, double lon) const;

private:
  double c1, c2, chi0, lam0s, r0;
  double map_radius;
};

/// Polar %Stereographic
class PolarStereo : public Projection
{
public:
  PolarStereo (const Params& pp);
  int Id () const { return (int)Projection::Method::PST; };

  mlib::erc GeoXY (double *x, double *y, double lat, double lon) const;
  mlib::erc XYGeo (double x, double y, double *lat, double *lon) const;
  double k (double lat, double lon) const;

private:
  double rho (double lat) const;
  double rho1;
  double sc[4];
};
