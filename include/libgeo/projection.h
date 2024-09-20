/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

///  \file PROJECTION.H	Projection class definition

#include "ellip.h"

#include <mlib/errorcode.h>
#include <assert.h>
#include <map>

class Projection
{
public:
  ///%Projection parameters
  class Params
  {
  public:
    /// Parameter codes (matching EPSG codes)
    enum class Code {
      SCALE_FACTOR = 8805,
      REF_LATITUDE = 8801,
      REF_LONGITUDE = 8802,
      NORTH_PARALLEL = 8823,
      SOUTH_PARALLEL = 8824,
      SKEW_AZIMUTH = 8814,
      FALSE_EASTING = 8806,
      FALSE_NORTHING = 8807,
      UNIT_FACTOR = 1051
    };
    Params (const Ellipsoid& ell);
    Params (Ellipsoid::well_known wk = Ellipsoid::WGS_84);

    Params& ellipsoid (const Ellipsoid& ell);
    Params& ellipsoid (Ellipsoid::well_known wk);
    Params& k0 (double k);
    Params& ref_latitude (double phi);
    Params& ref_longitude (double lambda);
    Params& skew_azimuth (double alpha);
    Params& north_latitude (double phin);
    Params& south_latitude (double phis);
    Params& false_north (double y);
    Params& false_east (double x);
    Params& unit (double factor);

    bool contains (Code p) const;
  private:

    Ellipsoid ellip_;
    std::map<Code, double>val;

    friend class Projection;
  };

  /// Projection method
  enum class Method {
    LCC = 9802, ///< \ref Lambert "Lambert conformal conical"
    MER = 9804, ///< Mercator
    TME = 9807, ///< \ref TransverseMercator "Transverse Mercator"
    STEREA = 9809, ///< Oblique stereographic
    OCY = 9814, ///< Oblique cylindrical (Swiss and EOV systems)
    HOM = 9815, ///< \ref ObliqueMercator "Hotine Oblique Mercator" (Alaska)
    RSO = 9812, ///< \ref RSO "Rectified Skew Orthomorphic"
    AEQD = 9820, ///< Azimuthal Equidistant
    EQC = 9825, ///< Equidistant cylindrical (Plate-Carrée)
    AEA = 9822, ///< \ref Albers "Albers Equal Area"
    CASS = 9806, ///< \ref Cassini "Cassini-Soldner"
    PST = 9810, ///< Polar Stereographic
    POL = 9818, ///< Polyconic
    LAEA = 1027, ///< Lambert Azimuthal Equal Area
    LEA = 9820, ///< Lambert Equal Area (not implemented)
    ORT = 9840,
  };

  /// Return projection method
  virtual int Id () const = 0;

  /// Convert from XY to geographical coordinates
  virtual mlib::erc XYGeo (double x, double y, double *lat, double *lon) const = 0;

  /// Convert from geographical to XY coordinates
  virtual mlib::erc GeoXY (double *x, double *y, double lat, double lon) const = 0;

  /// Return average scale factor at a given point
  virtual mlib::erc Scale (double x, double y, double *scale);

  /// Return scale factor along the latitude
  virtual double k (double lat, double lon) const;

  /// Return scale factor along meridian
  virtual double h (double lat, double lon) const;

  /// Return convergence value
  virtual double Convergence (double lat, double lon) const;

  const Ellipsoid& ellipsoid () const;
  double unit () const;
  double k0 () const;
  double ref_longitude () const;
  double ref_latitude () const;
  double north_parallel () const;
  double south_parallel () const;
  double false_east () const;
  double false_north () const;
  double skew_azimuth () const;

  static Projection* create (Method proj, const Projection::Params& pp);

protected:
  Projection (const Params& pars);

  double lon_adjust (double lon) const;
  void fwd_adj (double& x, double& y) const;
  void inv_adj (double& x, double& y) const;

  Params par_;

};

/*==================== INLINE FUNCTIONS ===========================*/
/// Set scale factor at origin
inline
Projection::Params& Projection::Params::k0 (double k)
{
  val[Code::SCALE_FACTOR] = k;
  return *this;
}

/// Set reference latitude
inline
Projection::Params& Projection::Params::ref_latitude (double phi)
{
  assert (-M_PI / 2 <= phi && phi <= M_PI / 2);
  val[Code::REF_LATITUDE] = phi;
  return *this;
}

/// Set reference longitude (central meridian)
inline
Projection::Params& Projection::Params::ref_longitude (double lambda)
{
  assert (-M_PI <= lambda && lambda <= M_PI);
  val[Code::REF_LONGITUDE] = lambda;
  return *this;
}

/// Set azimuth of central line
inline
Projection::Params& Projection::Params::skew_azimuth (double alpha)
{
  assert (-M_PI <= alpha && alpha <= M_PI);
  val[Code::SKEW_AZIMUTH] = alpha;
  return *this;
}

/// Set false easting
inline
Projection::Params& Projection::Params::false_east (double x)
{
  val[Code::FALSE_EASTING] = x;
  return *this;
}

/// Set unit conversion factor to meters
inline 
Projection::Params& Projection::Params::unit (double factor)
{
  val[Code::UNIT_FACTOR] = factor;
  return *this;
}

/// Set false northing
inline
Projection::Params& Projection::Params::false_north (double y)
{
  val[Code::FALSE_NORTHING] = y;
  return *this;
}

inline
Projection::Params& Projection::Params::north_latitude (double phin)
{
  assert (-M_PI / 2 <= phin && phin <= M_PI / 2);
  val[Code::NORTH_PARALLEL] = phin;
  return *this;
}

inline
Projection::Params& Projection::Params::south_latitude (double phis)
{
  assert (-M_PI / 2 <= phis && phis <= M_PI / 2);
  val[Code::SOUTH_PARALLEL] = phis;
  return *this;
}

inline
bool Projection::Params::contains (Code p) const
{
  return val.contains (p);
}

//-----------------------------------------------------------------------------

/// Return ellipsoid associated with this projection
inline
const Ellipsoid& Projection::ellipsoid () const
{
  return par_.ellip_;
}

/// Return conversion factor from distance units to meters
inline 
double Projection::unit () const
{
  if (par_.val.contains (Params::Code::UNIT_FACTOR))
    return par_.val.at (Params::Code::UNIT_FACTOR);

  return 1.0;
}

/// Return scale factor at origin
inline 
double Projection::k0() const 
{ 
  if (par_.val.contains (Params::Code::SCALE_FACTOR))
    return par_.val.at (Params::Code::SCALE_FACTOR);

  return 1.; //default value
}

/// Return central meridian
inline 
double Projection::ref_longitude() const 
{ 
  if (par_.val.contains (Params::Code::REF_LONGITUDE))
    return par_.val.at(Params::Code::REF_LONGITUDE);

  return 0.; //default value
}


/// Return reference latitude
inline 
double Projection::ref_latitude() const 
{
  if (par_.val.contains (Params::Code::REF_LATITUDE))
    return par_.val.at(Params::Code::REF_LATITUDE);

  return 0.; //default value
}

/// Return north parallel
inline 
double Projection::north_parallel() const 
{ 
  if (par_.val.contains (Params::Code::NORTH_PARALLEL))
    return par_.val.at (Params::Code::NORTH_PARALLEL);

  return 0.; //default value
}

/// Return south parallel
inline
double Projection::south_parallel() const 
{ 
  if (par_.val.contains (Params::Code::SOUTH_PARALLEL))
    return par_.val.at (Params::Code::SOUTH_PARALLEL);

  return 0.; //default value
};

/// Return X (easting) value at origin
inline 
double Projection::false_east() const 
{ 
  if (par_.val.contains (Params::Code::FALSE_EASTING))
    return par_.val.at (Params::Code::FALSE_EASTING);

  return 0.; //default value
};

/// Return Y (northing) value at origin
inline 
double Projection::false_north() const 
{ 
  if (par_.val.contains (Params::Code::FALSE_NORTHING))
    return par_.val.at (Params::Code::FALSE_NORTHING);

  return 0.; //default value
};

/// Return azimuth of skew
inline
double Projection::skew_azimuth() const 
{ 
  if (par_.val.contains (Params::Code::SKEW_AZIMUTH))
    return par_.val.at (Params::Code::SKEW_AZIMUTH);

  return 0.; //default value
};

/*!
  Default implementation for scale along the meridian returns the same value as
  the scale along the latitude k. This is true only for conformal projections
  and any derived projection that is not conformal would have to re-implement
  this function.
*/
inline
double Projection::h (double lat, double lon) const {return k(lat,lon);};

inline
double Projection::k (double lat, double lon) const {return 1.;};

///  Return an adjusted longitude between -M_PI and M_PI.
inline
double Projection::lon_adjust( double lon ) const
{
  return ((lon>=0)?1:-1)*(fmod(fabs(lon)+M_PI, 2*M_PI)-M_PI);
}

/*!
  Adjust X/Y coordinates for forward transformations

  Coordinates are converted to work units and added false easting and northing
*/
inline
void Projection::fwd_adj (double& x, double& y) const
{
  x = x / unit () + false_east ();
  y = y / unit () + false_north ();
}

/*!
  Adjust X/Y coordinates for inverse transformations

  False easting and northing are subtracted from input coordinates and
  resulting values are converted to meters. 
*/
inline
void Projection::inv_adj (double& x, double& y) const
{
  x = (x - false_east ()) * unit ();
  y = (y - false_north ()) * unit ();
}

/*!
  Default implementation for convergence function returns 0.
  This is true only for Mercator and a few other projections. Derived classes
  should re-implement it.
*/
inline
double Projection::Convergence (double lat, double lon) const
{
  return 0.;
}
