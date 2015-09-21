// Common utilities

#ifndef H_UTILITIES
#define H_UTILITIES

#include <vector>
#include "DataFormats/PatCandidates/interface/Electron.h"


typedef std::vector<const reco::Candidate*> vCands;

bool greaterPt(const reco::Candidate *a, const reco::Candidate *b);

#endif
