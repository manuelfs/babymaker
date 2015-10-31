#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "babymaker/bmaker/interface/weight_tools.hh"
#include <vector>

weight_tools::weight_tools() {}
weight_tools::~weight_tools() {}

float weight_tools::weight(weight_tools::variationType variation)
{
  if(weights.size()!=0) {
    return weights.at(variation).wgt/weights.at(nominal).wgt;
  }
  else return 1.0;
}

void weight_tools::getWeights(const edm::Event& iEvent)
{
  weights.clear();
  if(!iEvent.isRealData()) { 
    edm::Handle<LHEEventProduct> wLHEEventProduct;
    iEvent.getByLabel("externalLHEProducer", wLHEEventProduct);
    // these event weights are only defined in samples
    // generated from LHE files so an "else" is not desired
    if(wLHEEventProduct.isValid()) {
      weights = wLHEEventProduct->weights();
    }
  }
}
