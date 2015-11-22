#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "babymaker/bmaker/interface/weight_tools.hh"
#include <vector>

weight_tools::weight_tools(std::vector<double> puWeights):
  pileupWeights(puWeights)
{

}
weight_tools::~weight_tools() {}

float weight_tools::pileupWeight(unsigned int ntrupv)
{
  if(pileupWeights.size()!=0 && ntrupv<pileupWeights.size()) return static_cast<float>(pileupWeights.at(ntrupv));
  else return 1.0;
}


float weight_tools::topPtWeight(float toppt1, float toppt2){
  float pt1,pt2;
  if(toppt1>400) pt1=400;
  else pt1=toppt1;
  if(toppt2>400) pt2=400;
  else pt2=toppt2;
  return sqrt(exp(0.156-0.00137*pt1)*exp(0.156-0.00137*pt2));


}

float weight_tools::isrWeight(float isrpt){
  if(isrpt<400) return 0.;
  else if(isrpt<600) return 0.15;
  else return 0.3;


}

float weight_tools::theoryWeight(weight_tools::variationType variation)
{
  if(theoryWeights.size()!=0) {
    return theoryWeights.at(variation).wgt/theoryWeights.at(nominal).wgt;
  }
  else return 1.0;
}

void weight_tools::getTheoryWeights(const edm::Event& iEvent)
{
  theoryWeights.clear();
  if(!iEvent.isRealData()) { 
    edm::Handle<LHEEventProduct> wLHEEventProduct;
    iEvent.getByLabel("externalLHEProducer", wLHEEventProduct);
    // these event weights are only defined in samples
    // generated from LHE files so an "else" is not desired
    if(wLHEEventProduct.isValid()) {
      theoryWeights = wLHEEventProduct->weights();
    }
  }
}
