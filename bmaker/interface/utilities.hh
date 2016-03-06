// Common utilities

#ifndef H_UTILITIES
#define H_UTILITIES

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <fastjet/PseudoJet.hh>
#include "TString.h"

#define ERROR(x) do{throw std::runtime_error(string("Error in file ")+__FILE__+" at line "+to_string(__LINE__)+" (in "+__func__+"): "+x);}while(false)
#define DBG(x) do{std::cerr << "In " << __FILE__ << " at line " << __LINE__ << " (in function " << __func__ << "): " << x << std::endl;}while(false)

typedef std::vector<const reco::Candidate*> vCands;
typedef reco::Candidate::LorentzVector LVector;

namespace utilities{

  enum SysEnum{kSysJER, kSysJECUp, kSysJECDn, kSysLast};

  float dPhi(float phi1, float phi2);
  float dR(float phi1, float phi2, float eta1, float eta2);
  bool greaterPt(const reco::Candidate *a, const reco::Candidate *b);
  bool greaterM(const fastjet::PseudoJet &a, const fastjet::PseudoJet &b);
  float getMT(float pt1, float phi1, float pt2, float phi2);
  float getMT2(float pt1, float phi1, float pt2, float phi2, float met, float met_phi);
  std::string execute(const std::string &cmd);
  TString roundNumber(double num, int decimals, double denom=1.);
  TString addCommas(double num);

}

#endif
