#ifndef CRS_SNOMADR_H
#define CRS_SNOMADR_H

#include <streambuf>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

// R defines length() as a macro, which collides with C++ headers.
#ifdef length
#undef length
#endif

#ifdef error
#undef error
#endif

class Routbuf : public std::streambuf {
private:
  int overflow(int c) override {
    if (c != EOF) {
      Rprintf("%.1s", reinterpret_cast<char*>(&c));
    }
    return c;
  }
};

#endif
