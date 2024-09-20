/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <mlib/defs.h>
#include <libgeo/rso.h>

///  \file RSO.CPP Implementation of \ref RSO "Rectified Skew Orthomorphic" projection

/*!
  Rotation from uv coordinate system to XY. Also performs
  translation to false_east/false_north origin and unit change.
  (the unit for uv system is meters).
*/
void RSO::deskew (double u, double v, double *x, double *y) const
{
  auto [sskew, cskew] = sincos (gamma0);
  *x = cskew * v + sskew * u;
  *y = -sskew * v + cskew * u;
  fwd_adj (*x, *y);
}

/*!
  Rotation form xy coordinate system to uv
*/
void RSO::skew (double *u, double *v, double x, double y) const
{
  auto [sskew, cskew] = sincos (gamma0);
  inv_adj (x, y);
  *v = cskew * x - sskew * y;
  *u = sskew * x + cskew * y;
}
