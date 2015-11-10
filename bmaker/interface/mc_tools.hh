// MC_TOOLS: Functions that deal with MC particles

#ifndef H_MC_TOOLS
#define H_MC_TOOLS

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "babymaker/bmaker/interface/utilities.hh"

class mc_tools{

public:

  bool hasDaughter(const reco::GenParticle &mc, size_t id);
  bool isLast(const reco::GenParticle &mc, size_t id);
  bool decaysTo(const reco::GenParticle &mc, size_t id, const reco::GenParticle *&mcDau);
  void printParticle(const reco::GenParticle &mc);
  bool fromWOrWTau(const reco::GenParticle &mc);
  int mom(const reco::GenParticle &mc, const reco::GenParticle *&mcMom);
  bool isFromGSP(const reco::Candidate *c);
  bool isHardProcess(const int status);

  mc_tools();
  ~mc_tools();
};

#endif
