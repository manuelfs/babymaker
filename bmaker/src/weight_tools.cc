#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "babymaker/bmaker/interface/weight_tools.hh"
#include <vector>

using namespace std;

weight_tools::weight_tools(){
  // Pile up weights assuming sigma = 69 mb
  pileupWeights = vector<double>({178.037, 137.543, 90.6595, 23.9371, 11.4287,
	1.57655, 0.963025, 1.97081, 3.05875, 3.09653,
	3.07308, 2.94232, 2.38953, 1.56449, 0.822713,
	0.348372, 0.1212, 0.035636, 0.00948505, 0.00253883,
	0.000780246, 0.000305642, 0.000159835, 0.000105483, 7.88995e-05,
	6.23615e-05, 5.06906e-05, 4.16566e-05, 3.46035e-05, 2.82562e-05,
	2.07778e-05, 1.26429e-05, 6.47437e-06, 2.59927e-06, 9.77605e-07,
	3.32172e-07, 1.04265e-07, 3.09878e-08, 9.38906e-09, 2.68163e-09,
	7.26119e-10, 1.71182e-10, 4.61733e-11, 1.07015e-11, 2.49233e-12,
	5.48282e-13, 2.57415e-13, 7.73396e-14, 6.6361e-14, 0,
	0, 0});
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
