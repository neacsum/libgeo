/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#pragma once

/// \file ERRORS.H Error handling and error codes

#include <mlib/errorcode.h>

mlib::errfac& GeoErrors ();

enum GeoError {
  OK,           ///< No error
  PARAM,        ///< Invalid parameters
  RANGE,        ///< Invalid value range
  NONCONV       ///< Non-convergence
};
