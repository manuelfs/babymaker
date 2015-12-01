// MC cross sections

#ifndef H_CROSS_SECTIONS
#define H_CROSS_SECTIONS

#include "TString.h"

namespace xsec{

  float crossSection(const TString &file);
  void signalCrossSection(int glu_mass, double &xsec, double &xsec_unc);
  float fractionNegWeights(const TString &file);
}

#endif
