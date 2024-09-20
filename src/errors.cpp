/*
    Copyright (c) Mircea Neacsu (2020-2024)
    This is part of libgeo project and licensed under MIT license.
*/
#include <libgeo/errors.h>
/// \file ERRORS.CPP Implementation of error facility class 

/// Error facility 
class Errors : public mlib::errfac
{
public:
  Errors () : errfac ("LIBGEO error") {};
  std::string message (const mlib::erc& e) const;
};


std::string Errors::message (const mlib::erc& e) const
{
  static const char* errtab[] = {
    "No error",
    "Invalid parameter",
    "Invalid value range",
    "Non convergence",
    "Singularity",
  };
  int c = e.code ();
  if (0 <= c && c < _countof (errtab))
    return name () + ": " + errtab[c];
  else
    return errfac::message (e);
}

Errors g;

mlib::errfac& GeoErrors () {
  return g;
};
