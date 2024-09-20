/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/

#include <mlib/defs.h>
#include <libgeo/hom.h>

///  \file HOM.CPP Implementation of \ref Hotine "Hotine Oblique Mercator"

/*!
  Rotation from uv coordinate system to XY. Also performs
  translation to false_east/false_north origin and unit change.
  (the unit for uv system is meters).
*/
void Hotine::deskew (double u, double v, double *x, double *y) const
{
  auto [sskew, cskew] = sincos(skew_azimuth());
  *x = cskew * v + sskew * u;
  *y = -sskew * v + cskew * u;
  fwd_adj (*x, *y);
}

/*!
  Rotation from xy coordinate system to uv
*/
void Hotine::skew (double *u, double *v, double x, double y) const
{
  auto [sskew, cskew] = sincos (skew_azimuth ());
  inv_adj (x, y);
  *v = cskew * x - sskew * y;
  *u = sskew * x + cskew * y;
}