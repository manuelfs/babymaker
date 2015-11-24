#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "babymaker/bmaker/interface/weight_tools.hh"
#include <vector>

using namespace std;

weight_tools::weight_tools(){
  // Pile up weights assuming sigma = 69 mb
  pileupWeights = vector<double>({
    217.107, 133.498, 123.917, 34.7894, 19.9406, 
    3.8807, 2.22042, 2.76796, 3.42611, 3.16958, 
    2.93076, 2.7112, 2.2135, 1.50491, 0.846353, 
    0.397117, 0.16497, 0.0706261, 0.036932, 0.021373, 
    0.0113103, 0.00498929, 0.00181823, 0.000573713, 0.000173074, 
    5.65199e-05, 2.08184e-05, 8.10237e-06, 3.15537e-06, 1.16908e-06, 
    3.75369e-07, 9.58782e-08, 1.98111e-08, 3.08475e-09, 4.32512e-10, 
    5.26588e-11, 5.69234e-12, 5.60035e-13, 5.37918e-14, 4.09235e-15, 
    0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0
	});
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
