/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <libgeo/ellip.h>
#include <assert.h>

///  \file ELLIP.CPP Ellipsoid class implementation

/*!
  \class Ellipsoid
  Define an ellipsoid using the semi major axis and flattening.
  All other significant values are calculated based on a and f.
  Note that "f" is defined as (a-b)/a
*/


/// Known ellipsoids
static struct ell_data
{
  Ellipsoid::well_known wk;
  const char* name;
  double a;
  double f_1;
} ellip[] = {
  { Ellipsoid::well_known::WGS_84,                "WGS-84",                   6378137.,     298.257223563 },
  { Ellipsoid::well_known::AIRY,                  "Airy",                     6377563.396,  299.3249646 },
  { Ellipsoid::well_known::ATS_77,                "ATS-77",                   6378135.0,    298.257 },
  { Ellipsoid::well_known::AUSTRALIAN,            "Australian National",      6378160.,     298.25 },
  { Ellipsoid::well_known::AIRY_MODIFIED,         "Airy Modified",            6377340.189,  299.3249646 },
  { Ellipsoid::well_known::CLARKE_1858,           "Clarke 1858",              6378293.645,  294.2606764 },
  { Ellipsoid::well_known::BESSEL_1841,           "Bessel 1841",              6377397.155,  299.1528128 },
  { Ellipsoid::well_known::CLARKE_1866,           "Clarke 1866",              6378206.4,    294.9786982 },
  { Ellipsoid::well_known::CLARKE_1880,           "Clarke 1880",              6378249.145,  293.465 },
  { Ellipsoid::well_known::CLARKE_1880_ARC,       "Clarke-1880 (Arc)",        6378249.145,  293.4663077 },
  { Ellipsoid::well_known::CLARKE_1880_IGN,       "Clarke 1880 (IGN)",        6378249.200,  293.466021294, },
  { Ellipsoid::well_known::EVEREST,               "Everest",                  6377276.345,  300.8017 },
  { Ellipsoid::well_known::EVEREST_MALAY_SING,    "Everest (Malay. & Sing)",  6377304.063,  300.8017 },
  { Ellipsoid::well_known::EVEREST_INDIA_1956,    "Everest (India 1956)",     6377301.243,  300.8017 },
  { Ellipsoid::well_known::EVEREST_SABAH_SARWAK,  "Everest (Sabah Sarawak)",  6377298.556,  300.8017 },
  { Ellipsoid::well_known::EVEREST_MALAYSIA_1969, "Everest (Malaysia 1969)",  6377295.664,  300.8017 },
  { Ellipsoid::well_known::EVEREST_PAKISTAN,      "Everest (Pakistan)",       6377309.613,  300.8017 },
  { Ellipsoid::well_known::FISCHER_1960,          "Fischer 1960 (Mercury)",   6378166.,     298.3 },
  { Ellipsoid::well_known::FISCHER_1960_MODIFIED, "Fischer 1960 Modified",    6378155.,     298.3 },
  { Ellipsoid::well_known::FISCHER_1968,          "Fischer 1968",             6378150.,     298.3 },
  { Ellipsoid::well_known::GRS_1967,              "GRS-1967",                 6378160.,     298.247167427 },
  { Ellipsoid::well_known::GRS_1980,              "GRS-1980",                 6378137.,     298.257222101 },
  { Ellipsoid::well_known::HELMERT_1906,          "Helmert 1906",             6378200.,     298.3 },
  { Ellipsoid::well_known::HOUGH,                 "Hough",                    6378270.,     297. },
  { Ellipsoid::well_known::INDONESIAN_1974,       "Indonesian 1974",          6378160.,     298.247 },
  { Ellipsoid::well_known::INTERNATIONAL,         "International",            6378388.,     297. },
  { Ellipsoid::well_known::KRASSOVSKY,            "Krassovsky",               6378245.,     298.3 },
  { Ellipsoid::well_known::SOUTH_AMERICAN_1969,   "SA-1969",                  6378160.,     298.25 },
  { Ellipsoid::well_known::WGS_60,                "WGS-60",                   6378165.,     298.3 },
  { Ellipsoid::well_known::WGS_66,                "WGS-66",                   6378145.,     298.25 },
  { Ellipsoid::well_known::WGS_72,                "WGS-72",                   6378135.,     298.26 },
};

constexpr int num_ellip = sizeof (ellip) / sizeof (ell_data);
static int find_wkidx (double a, double f);


/// Create an ellipsoid object with the given a and f values.
Ellipsoid::Ellipsoid (double a, double f) :
  a_ (a),
  f_ (f)
{
  e2_ = 2. * f_ - f_ * f_;
  e_ = sqrt (e2_);
  int idx;
  if ((idx = find_wkidx (a_, f_)) >= 0)
    name_ = ellip[idx].name;
}

/// Create a well-known ellipsoid
Ellipsoid::Ellipsoid (well_known wk)
{
  assert (wk < num_ellip);
  name_ = known (wk, &a_, &f_);
  e2_ = 2. * f_ - f_ * f_;
  e_ = sqrt (e2_);
}

/// Copy Constructor
Ellipsoid::Ellipsoid (const Ellipsoid& ell) :
  a_ (ell.a_),
  f_ (ell.f_),
  e2_ (ell.e2_),
  e_ (ell.e_),
  name_ (ell.name_)
{
}

/// Assignment operator
Ellipsoid& Ellipsoid::operator= (const Ellipsoid& rhs)
{
  a_ = rhs.a_;
  f_ = rhs.f_;
  e2_ = rhs.e2_;
  e_ = rhs.e_;
  name_ = rhs.name_;
  return *this;
}

///  Return radius of curvature in the prime vertical
///  \f$ rn = \frac a{\sqrt{1-e^2*\sin(lat)^2}} \f$
double Ellipsoid::rn (double lat) const
{
  double slat;
  return a_ / sqrt (1 - e2_ * (slat = sin (lat))*slat);
}

/// Return radius of curvature in the prime meridian
/// \f$ rm = \frac{a*(1-e^2)}{\sqrt{[1-e^2*\sin(lat)^2)]^3}} \f$
double Ellipsoid::rm (double lat) const
{
  double slat, root;
  return a_ * (1 - e2_) / ((root = sqrt (1 - e2_ * (slat = sin (lat))*slat))*root*root);
}

/// Return length from Equator along meridian
double Ellipsoid::lm (double lat) const
{
  if (!c)
  {
    const_cast<Ellipsoid*>(this)->c = std::make_unique<double[] > (4);
    c[0] = ((-0.01953125 * e2_ - 0.046875) * e2_ - 0.250) * e2_ + 1;
    c[1] = ((0.0439453125 * e2_ + 0.093750) * e2_ + 0.375) * e2_;
    c[2] = (0.0439453125 * e2_ + 0.05859375) * e2_ * e2_;
    c[3] =  - 0.01139322916667 * e2_ * e2_ * e2_;
  }
  return a_ * (c[0] * lat - c[1] * sin (2 * lat) + c[2] * sin (4 * lat) + c[3] * sin (6 * lat));
}

/*!
  Return parameters (name, a and f ) of a known ellipsoid

  \param index    index in table of known ellipsoids
  \param pa       pointer to semi-major axis
  \param pf       pointer to inverse of flattening
  \return         ellipsoid name
*/
const char *Ellipsoid::known (well_known wk, double* pa, double* pf)
{
  int i;
  assert (wk < num_ellip);

  //better safe than sorry - ellip table and well_known lists might not be in the same order
  for (i = 0; i < num_ellip; i++)
  {
    if (ellip[i].wk == wk)
      break;
  }
  if (pa)
    *pa = ellip[i].a;
  if (pf)
    *pf = 1. / ellip[i].f_1;
  return ellip[i].name;
}

/*!
  Search the table of well-known ellipsoids for a matching ellipsoid and return its index.
  If not found, the function returns a negative value.
*/
static int find_wkidx (double a, double f)
{
  const int max_index = sizeof (ellip) / sizeof (ell_data) - 1;
  for (int i = 0; i < max_index; i++)
  {
    if (fabs (a - ellip[i].a) < 1e-4 && fabs (1./f - ellip[i].f_1) < 1e-8)
      return i;
  }
  return -1;
}

/*!
  Conversion from geographical coordinates to geocentric.

  Formulas from "GPS Satellite Surveying" by Alfred Leick page 184
*/
void Ellipsoid::geo_ECEF (double lat, double lon, double height, double *x, double *y, double *z) const
{
  double slat = sin (lat);
  double clat = cos (lat);
  double r = rn (lat);
  *x = (r + height)*clat*cos (lon);
  *y = (r + height)*clat*sin (lon);
  *z = ((1 - e2_)*r + height)*slat;
}

/*!
  Conversion from geocentric to geographical coordinates.
  
  Formulas from "GPS Satellite Surveying" by Alfred Leick page 184-185.
  If both x and y are 0 longitude is set to 0.
*/
void Ellipsoid::ECEF_geo (double *lat, double *lon, double *height, double x, double y, double z) const
{
  double rho = hypot (x, y);
  if (rho == 0.)
  {  //deal with North and South Pole cases
    *lat = (z > 0) ? M_PI / 2 : -M_PI / 2;
    *lon = 0.;
    *height = fabs (z) - b ();
  }
  else
  {
    double lat0, lat1;
    double rn;
    lat0 = atan (1 / (1 - e2_) * z / rho);
    lat1 = lat0 + 1;
    while (fabs (lat0 - lat1) > 1e-9) /*approx 2 mm*/
    {
      double slat = sin (lat0);
      rn = a_ / sqrt (1. - e2_ * slat*slat);
      lat1 = lat0;
      lat0 = atan ((z + e2_ * rn*slat) / rho);
    }
    *lon = atan2 (y, x);
    *lat = lat0;
    *height = rho / cos (lat0) - rn;
  }
}

/*!
  Calculate geodetic distance between 2 points using Gauss mid-latitude solution. 
  
  Formulas from "GPS Satellite Surveying" by Alfred Leick page 281-282
*/
double Ellipsoid::gcirc (double lat1, double lon1, double lat2, double lon2, double *az) const
{
  double mid_phi = (lat1 + lat2) / 2.;
  double dphi = lat2 - lat1;
  double dlam = lon2 - lon1;
  double smp = sin (mid_phi);
  double cmp = cos (mid_phi);

  double eta2 = ep2 ()*cmp*cmp;
  double V4 = (1 - eta2)*(1 - eta2);
  double t2 = tan (mid_phi)*tan (mid_phi);

  double C[8] { 
    rm (mid_phi),
    rn (mid_phi), 
    1. / 24., 
    (1 + 9.*eta2*eta2*t2) / (24.*V4),
    (1 - 2 * eta2) / 24.,
    eta2*(1 - t2) / (8 * V4),
    (1 + eta2) / 12.,
    (3 + 8 * eta2) / (24 * V4)
  };

  double ssin = C[1] * dlam*cmp*(1 - C[2] * dlam*dlam*smp*smp + C[3] * dphi*dphi);
  double scos = C[0] * cos (dlam / 2)*dphi*(1 + C[4] * dlam*dlam*cmp*cmp + C[5] * dphi*dphi);
  double dalf = dlam * smp*(1 + C[6] * dlam*dlam*cmp*cmp + C[7] * dphi*dphi);
  if (az)
    *az = atan2 (ssin, scos) - dalf / 2;
  return hypot (ssin, scos);
}

/*!
  Calculate rhumb line distance between 2 points.
*/
double Ellipsoid::rhumb (double lat1, double lon1, double lat2, double lon2, double *az) const
{
  double dist;

  double sl1 = sin (lat1);
  double sl2 = sin (lat2);

  double azm = atan2 (lon2 - lon1,
    log ((1 + sl2) / (1 - sl2)*(1 - sl1) / (1 + sl1) *
      pow ((1 - e_ * sl2) / (1 + e_ * sl2)*(1 + e_ * sl1) / (1 - e_ * sl1), e_))
  );

  if (az)
    *az = azm;

  if (lat1 == lat2)
    dist = fabs (rn (lat1)* cos (lat1) * (lon1 - lon2));
  else
    dist = fabs ((lm (lat2) - lm (lat1)) / cos (azm));
  return dist;
}

