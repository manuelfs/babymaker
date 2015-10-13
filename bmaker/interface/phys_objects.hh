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

  bool isInJSON(std::string type, int run, int lumiblock);
  bool hasGoodPV(edm::Handle<reco::VertexCollection> vtx);
  float getMT(float pt1, float phi1, float pt2, float phi2);
}

#endif
