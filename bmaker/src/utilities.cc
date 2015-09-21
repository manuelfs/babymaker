// Common utilities
 
#include <cmath>
#include "babymaker/bmaker/interface/utilities.hh"

bool greaterPt(const reco::Candidate *a, const reco::Candidate *b) {
  return a->pt() > b->pt();
}

