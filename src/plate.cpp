/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <mlib/defs.h>
#include <libgeo/plate.h>

/// \file PLATE.CPP Implementation of \ref PlateCarree "Plate Carrée" projection

using mlib::erc;

PlateCarree::PlateCarree (const Params& par) :
  Projection (par)
{
}

erc PlateCarree::XYGeo (double x, double y, double *lat, double *lon) const
{
  inv_adj (x, y);
  *lat = y / ellipsoid().a ();
  *lon = ref_longitude() + x / (ellipsoid().a ()*cos (ref_latitude()));
  return erc::success;
}

erc PlateCarree::GeoXY (double *x, double *y, double lat, double lon) const
{
  lon = lon_adjust (lon - ref_longitude());
  *x = lon * ellipsoid().a () * cos (ref_latitude());
  *y = lat * ellipsoid().a ();
  fwd_adj (*x, *y);
  return erc::success;
}

double PlateCarree::k (double lat, double lon) const
{
  return cos (ref_latitude()) / cos (lat);
}

