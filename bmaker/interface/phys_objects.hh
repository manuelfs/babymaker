// Common physics objects definitions

#ifndef H_PHYS_OBJECTS
#define H_PHYS_OBJECTS

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "babymaker/bmaker/interface/utilities.hh"

namespace phys_objects{

  ///////////////// CUTS ///////////////////////
  float JetPtCut		= 30.0;
  float JetEtaCut               = 2.4;
  float JetMHTEtaCut            = 5.0;
  float JetHLTPtCut		= 40.0;
  float JetHLTEtaCut            = 3.0;
  float CSVLoose                = 0.605;
  float CSVMedium               = 0.890;
  float CSVTight                = 0.970;

  bool isInJSON(std::string type, int run, int lumiblock);

  bool leptonInJet(const pat::Jet &jet, vCands leptons);
  bool idJet(const pat::Jet &jet);

  bool hasGoodPV(edm::Handle<reco::VertexCollection> vtx);
}

#endif
