// Common utilities

#ifndef H_UTILITIES
#define H_UTILITIES

#include <vector>
#include <string>
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <fastjet/PseudoJet.hh>
#include "TString.h"

typedef std::vector<const reco::Candidate*> vCands;

namespace utilities{

  bool greaterPt(const reco::Candidate *a, const reco::Candidate *b);
  bool greaterM(const fastjet::PseudoJet &a, const fastjet::PseudoJet &b);
  std::string execute(const std::string &cmd);
  TString roundNumber(double num, int decimals, double denom=1.);
  float crossSection(const TString &file);

}

#endif
