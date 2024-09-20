/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <libgeo/libgeo.h>

///  \file PROJECTION.CPP Implementation of Projection class

using mlib::erc;

Projection::Params::Params (const Ellipsoid& ell)
  : ellip_ (ell)
{
}

Projection::Params::Params (Ellipsoid::well_known wk)
  : ellip_ (wk)
{
}

Projection::Params& Projection::Params::ellipsoid (const Ellipsoid& ell)
{
  ellip_ = ell;
  return *this;
}

Projection::Params& Projection::Params::ellipsoid (Ellipsoid::well_known wk)
{
  ellip_ = Ellipsoid (wk);
  return *this;
}




/*!
  \class Projection
  This class adds functions needed go between
  geographical coordinates and XY coordinates. It is an abstract class
  which serves as a base to different projection systems.
*/

/*!
  Protected constructor invoked by derived classes to initialize all projection parameters
*/
Projection::Projection (const Params& params)
  : par_ (params)
{
}

/*!
  Default scale factor function. Converts X/Y coordinates to lat/lon and calls
  h and k functions to calculate scale factor along latitude and meridian.
  Returned value is the geometric mean of the h and k values.
*/
erc Projection::Scale (double x, double y, double *scale)
{
  double lat, lon;
  try {
    XYGeo (x, y, &lat, &lon);
    *scale = sqrt (h (lat, lon)*k (lat, lon));
    return erc::success;
  }
  catch (erc& ec) {
    return ec;
  }
}


/*!
  Projection factory function creates a projection based on given parameters.
*/
Projection* Projection::create (Method id, const Params& par)
{
  switch (id)
  {
  case Method::AEA:   return new Albers (par);
  case Method::LAEA:   return new AzEqArea (par);
  case Method::AEQD:   return new AzEqDist (par);
  case Method::CASS:   return new Cassini (par);
  case Method::HOM:   return new Hotine (par);
  case Method::LCC:   return new Lambert (par);
  case Method::MER:   return new Mercator (par);
  case Method::OCY:   return new ObliqueCylindrical (par);
  case Method::STEREA:   return new Stereographic (par);
  case Method::PST:   return new PolarStereo (par);
  case Method::EQC:   return new PlateCarree (par);
  case Method::POL:   return new Polyconic (par);
  case Method::RSO:   return new RSO (par);
  case Method::TME:   return new TransverseMercator (par);
  }
  return NULL;
}
