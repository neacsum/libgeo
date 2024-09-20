/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

/// \file UNITS.H Measurement units and conversion functions

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#endif

#define D2R       0.01745329251994          ///< Degrees to radians conversion factor
#define MAS       (M_PI / (180 * 3600000.)) ///< milli-arcsecond
#define A_WGS84   6378137.000               ///< Semimajor axis of WGS84 ellipsoid
#define F_WGS84   0.003352810664747         ///< Flattening of WGS84 ellipsoid
#define F1_WGS84  298.257223563             ///< Inverse of flattening for WGS84 ellipsoid

#define NM2M      1852.                     ///< Nautical mile to meters conversion factor

/*
       Conversion factors from customary units to meter
*/
// International
#define FT2M      0.3048                    ///< International foot to meter conversion factor
#define YD2M      0.9144                    ///< International yard - 1 yd = 3 ft
#define CH2M      20.1168                   ///< International chain - 1 ch = 22 yd
#define LK2M      0.201168                  ///< International link  - 1 lk = 1/100 ch

// US survey
#define FT_US2M   (1200 / 3937.)            ///< US Survey foot to meter conversion factor (1200/3937)
#define CH_US2M   (79200 / 3937.)           ///< US chain  - 1 chUS = 66 ftUS
#define LK_US2M   (792 / 3937.)             ///< US link - 1 lkUS = 1/100 chUS

// Clarke 1865
#define FT_CLA2M   0.3047972654             ///< Foot (Clarke) to meter conversion factor
#define YD_CLA2M   (3 * 0.3047972654)       ///< Yard (Clarke) - 1 ydCla = 3 ftCla
#define CH_CLA2M   (66 * 0.3047972654)      ///< Chain (Clarke) - 1chCla = 22 ydCla
#define LK_CLA2M   0.201166195164           ///< Link (Clarke) - 1 chCla = 100 lkCla

// Benoit 1895
#define FT_BEN2M  (12 / 39.370113)          ///< Foot (Benoit) to meter
#define YD_BEN2M  (36 / 39.370113)          ///< Yard (Benoit) - 1 ydBen = 3 ftBen
#define CH_BEN2M  (792/ 39.370113)          ///< Chain (Benoit) - 1 chBen = 22 ydBen
#define LK_BEN2M  (7.92 / 39.370113)        ///< Link (Benoit) - 1lkBen = 1/100 chBen

// Sears 1922
#define FT_SEA2M  (12 / 39.370147)          ///< Foot (Sears) to meter
#define YD_SEA2M  (36 / 39.370147)          ///< Yard (Sears) - 1 ydSea = 3 ftSea
#define CH_SEA2M  (792 / 39.370147)         ///< Chain (Sears) - 1 chSea = 22 ydSea
#define LK_SEA2M  (7.92 / 39.370147)        ///< Link (Sears) - 1 lkSea = 1/100 chSea


/// Milli-arcseconds literal operator converts a value to radians
///@{
constexpr double operator ""_mas (long double val)
{
  return val * MAS;
}

constexpr double operator ""_mas (unsigned long long val)
{
  return val * MAS;
}
///@}

// ---------- Literal operators for international units ----------------------

/// International foot literal operator converts a value to meters
///@{
constexpr double operator ""_ft (long double val)
{
  return val * FT2M;
}

constexpr double operator ""_ft (unsigned long long val)
{
  return val * FT2M;
}
///@}

/// Yard to meter literal operator converts a value to meters
///@{
constexpr double operator ""_yd (long double val)
{
  return val * YD2M;
}
constexpr double operator ""_yd (unsigned long long val)
{
  return val * YD2M;
}
///@}

/// International chain literal operator converts a value to meters
///@{
constexpr double operator ""_ch (long double val)
{
  return val * CH2M;
}

constexpr double operator ""_ch (unsigned long long val)
{
  return val * CH2M;
}
///@}

/// International link literal operator converts a value to meters
///@{
constexpr double operator ""_lk (long double val)
{
  return val * LK2M;
}

constexpr double operator ""_lk (unsigned long long val)
{
  return val * LK2M;
}
///@}

// ---------- Literal operators for US survey units ----------------------

/// US survey foot literal operator converts a value to meters
///@{
constexpr double operator ""_ftUS (long double val)
{
  return val * FT_US2M;
}

constexpr double operator ""_ftUS (unsigned long long val)
{
  return val * FT_US2M;
}
///@}

/// US chain literal operator converts a value to meters
///@{
constexpr double operator ""_chUS (long double val)
{
  return val * CH_US2M;
}

constexpr double operator ""_chUS (unsigned long long val)
{
  return val * CH_US2M;
}
///@}

/// US survey link literal operator converts a value to meters
///@{
constexpr double operator ""_lkUS (long double val)
{
  return val * LK_US2M;
}

constexpr double operator ""_lkUS (unsigned long long val)
{
  return val * LK_US2M;
}
///@}

// ---------- Literal operators for Clarke's units ----------------------

/// Clarke's foot literal operator converts a value to meters
///@{
constexpr double operator ""_ftCla (long double val)
{
  return val * FT_CLA2M;
}

constexpr double operator ""_ftCla (unsigned long long val)
{
  return val * FT_CLA2M;
}
///@}

/// Clarke's yard to meter literal operator converts a value to meters
///@{
constexpr double operator ""_ydCla (long double val)
{
  return val * YD_CLA2M;
}
constexpr double operator ""_ydCla (unsigned long long val)
{
  return val * YD_CLA2M;
}
///@}

/// Clarke's chain literal operator converts a value to meters
///@{
constexpr double operator ""_chCla (long double val)
{
  return val * CH_CLA2M;
}

constexpr double operator ""_chCla (unsigned long long val)
{
  return val * CH_CLA2M;
}
///@}

/// Clarke's link literal operator converts a value to meters
///@{
constexpr double operator ""_lkCla (long double val)
{
  return val * LK_CLA2M;
}

constexpr double operator ""_lkCla (unsigned long long val)
{
  return val * LK_CLA2M;
}
///@}

// ---------- Literal operators for Benoit units ----------------------

/// Foot (Benoit) literal operator converts a value to meters
///@{
constexpr double operator ""_ftBen (long double val)
{
  return val * FT_BEN2M;
}

constexpr double operator ""_ftBen (unsigned long long val)
{
  return val * FT_BEN2M;
}
///@}

/// Yard (Benoit) to meter literal operator converts a value to meters
///@{
constexpr double operator ""_ydBen (long double val)
{
  return val * YD_BEN2M;
}
constexpr double operator ""_ydBen (unsigned long long val)
{
  return val * YD_BEN2M;
}
///@}

/// Chain (Benoit) literal operator converts a value to meters
///@{
constexpr double operator ""_chBen (long double val)
{
  return val * CH_BEN2M;
}

constexpr double operator ""_chBen (unsigned long long val)
{
  return val * CH_BEN2M;
}
///@}

/// Link (Benoit) literal operator converts a value to meters
///@{
constexpr double operator ""_lkBen (long double val)
{
  return val * LK_BEN2M;
}

constexpr double operator ""_lkBen (unsigned long long val)
{
  return val * LK_BEN2M;
}
///@}

// ---------- Literal operators for Sears units ----------------------

/// Foot (Sears) literal operator converts a value to meters
///@{
constexpr double operator ""_ftSea (long double val)
{
  return val * FT_SEA2M;
}

constexpr double operator ""_ftSea (unsigned long long val)
{
  return val * FT_SEA2M;
}
///@}

/// Yard (Sears) to meter literal operator converts a value to meters
///@{
constexpr double operator ""_ydSea (long double val)
{
  return val * YD_SEA2M;
}
constexpr double operator ""_ydSea (unsigned long long val)
{
  return val * YD_SEA2M;
}
///@}

/// Chain (Sears) literal operator converts a value to meters
///@{
constexpr double operator ""_chSea (long double val)
{
  return val * CH_SEA2M;
}

constexpr double operator ""_chSea (unsigned long long val)
{
  return val * CH_SEA2M;
}
///@}

/// Link (Sears) literal operator converts a value to meters
///@{
constexpr double operator ""_lkSea (long double val)
{
  return val * LK_SEA2M;
}

constexpr double operator ""_lkSea (unsigned long long val)
{
  return val * LK_SEA2M;
}
///@}



/// Nautical miles literal operator  converts a value to meters
///@{
constexpr double operator ""_nmi (long double val)
{
  return val * NM2M;
}

constexpr double operator ""_nmi (unsigned long long val)
{
  return val * NM2M;
}
///@}

// -------------------- Degrees conversion functions --------------------------
/// Convert decimal degrees to radians
constexpr double DEG (double dd)
{
  return dd * D2R;
}

/// Conversion from degrees to radians
constexpr double D2rad (double val)
{
  return val * D2R;
}

/// Degrees literal operator converts a value to radians
///@{
constexpr double operator ""_deg (long double val)
{
  return val * D2R;
}

constexpr double operator ""_deg (unsigned long long val)
{
  return val * D2R;
}
///@}

/// Convert degrees, minutes to radians
constexpr double DM (double dd, double mm)
{
  return (dd + mm / 60.) * D2R;
}

/// Conversion to decimal degrees from DDMM.mmm
constexpr double DM2deg (double value)
{
  int sign = (value >= 0) ? 1 : -1;
  if (value < 0)
    value = -value;
  int deg = (int)(value / 100.);
  value -= (double)deg * 100.;
  return sign * (deg + value / 60.);
}

/// Conversion from degrees, minutes (DDMM.mmm) to radians
constexpr double DM2rad (double val)
{
  return DM2deg (val) * D2R;
}

/// Conversion from decimal degrees to degrees, minutes (DDMM.mmm)
constexpr double deg2DM (double value)
{
  int deg = (int)value;
  return (value - deg) * 60. + deg * 100.;
}

/// Conversion from radians to degrees, minutes (DDMM.mmm)
constexpr double rad2DM (double val)
{
  return deg2DM (val / D2R);
}

/// Degrees-minutes literal operator converts a value to radians
///@{
constexpr double operator ""_dm (long double val)
{
  return DM2rad (val);
}

constexpr double operator ""_dm (unsigned long long val)
{
  return DM2rad ((double)val);
}
///@}


/// Convert degrees, minutes seconds to radians
constexpr double DMS (double dd, double mm, double ss)
{
  return (dd + mm / 60. + ss / 3600.) * D2R;
}

/// Conversion to decimal degrees from DDMMSS.ssss
constexpr double DMS2deg (double value)
{
  int sign = (value >= 0) ? 1 : -1;
  if (value < 0)
    value = -value;
  int deg = (int)(value / 10000.);
  value -= (double)deg * 10000;
  int min = (int)(value / 100.);
  value -= (double)min * 100.;
  return sign * (deg + min / 60. + value / 3600.);
}

/// Conversion from degrees, minutes, seconds (DDMMSS.sss) to radians
constexpr double DMS2rad (double val)
{
  return DMS2deg (val) * D2R;
}

/// Degrees-minutes-seconds literal operator converts a value to radians
///@{
constexpr double operator ""_dms (long double val)
{
  return DMS2rad (val);
}

constexpr double operator ""_dms (unsigned long long val)
{
  return DMS2rad ((double)val);
}
///@}
