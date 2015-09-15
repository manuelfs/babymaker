// Common physics objects definitions

#ifndef H_PHYS_OBJECTS
#define H_PHYS_OBJECTS

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
		      const reco::Candidate* ptcl,  
		      double r_iso_min, double r_iso_max, double kt_scale,
		      bool charged_only);


#endif
