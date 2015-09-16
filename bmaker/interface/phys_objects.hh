// Common physics objects definitions

#ifndef H_PHYS_OBJECTS
#define H_PHYS_OBJECTS

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

namespace phys_objects{
  enum CutLevel{kVeto, kLoose, kMedium, kTight};
  template<class T>
  T chooseVal(CutLevel threshold, T valVeto, T valLoose, T valMedium, T valTight){
   switch(threshold){
   default:
   case kVeto:
     return valVeto;
   case kLoose:
     return valLoose;
   case kMedium:
     return valMedium;
   case kTight:
     return valTight;
   }
   return valVeto;
 }


  double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			const reco::Candidate* ptcl,  
			double r_iso_min, double r_iso_max, double kt_scale,
			bool charged_only);
  bool IdElectron(const pat::Electron &lep, CutLevel threshold, edm::Handle<reco::VertexCollection> vtx, bool do_iso);
}

#endif
