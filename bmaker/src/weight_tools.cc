#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "babymaker/bmaker/interface/weight_tools.hh"
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

weight_tools::weight_tools(){
  theoryWeights.clear();

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


float weight_tools::triggerEfficiency(int &nmus, int &nels, vector<float> &sys_trig){
  sys_trig.resize(2,1.);
  int nleps(nmus+nels);

  if(nleps == 0 || nleps > 2) return 1.;

  float eff_trig(1.);
  if(nleps == 1){
    if(nels==1){
      // eff = 0.941 +0.6-0.7 +-1.0 = 0.941 +0.012 - 0.012
      eff_trig = 0.941;
      sys_trig[0] = (eff_trig+0.012)/eff_trig;
      sys_trig[1] = (eff_trig-0.012)/eff_trig;
    }
    if(nmus==1){
      // eff = 0.951 +0.5-0.5 +-1.0 = 0.951 +0.011 - 0.011
      eff_trig = 0.951;
      sys_trig[0] = (eff_trig+0.011)/eff_trig;
      sys_trig[1] = (eff_trig-0.011)/eff_trig;
    } 
  } // nleps == 1
  if(nleps == 2){
    if(nels==2){
      // eff = 0.983 +0.5-0.7 +-0.5 = 0.983 +0.007 - 0.009
      eff_trig = 0.983;
      sys_trig[0] = (eff_trig+0.007)/eff_trig;
      sys_trig[1] = (eff_trig-0.009)/eff_trig;
    }
    if(nmus==2){
      // eff = 0.991 +0.3-0.5 +-0.5 = 0.991 +0.006 - 0.007
      eff_trig = 0.991;
      sys_trig[0] = (eff_trig+0.006)/eff_trig;
      sys_trig[1] = (eff_trig-0.007)/eff_trig;
    }
    if(nels==1&&nmus==1){
      // eff = 0.989 +0.5-0.9 +-0.5 = 0.989 +0.007 - 0.010
      eff_trig = 0.989;
      sys_trig[0] = (eff_trig+0.007)/eff_trig;
      sys_trig[1] = (eff_trig-0.010)/eff_trig;
    }
  } // nleps == 2
  return eff_trig;
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

float weight_tools::theoryWeight(weight_tools::variationType variation){
  if(theoryWeights.size()!=0) {
    return theoryWeights.at(variation).wgt/theoryWeights.at(nominal).wgt;
  }
  else return 1.0;
}

void weight_tools::getTheoryWeights(edm::Handle<LHEEventProduct> lhe_info){
  theoryWeights.clear();
  theoryWeights = lhe_info->weights();
}

void weight_tools::getPDFWeights(vector<float> &sys_pdf, vector<float> &w_pdf){
  if(theoryWeights.size()!=0) {
    unsigned ind(10), nweights(100); //index of the first pdf weight and number of replicas
    vector<double> pdfwgt = vector<double>(nweights,1.);
    for (unsigned i(0); i<nweights; i++){
      double ipdfw = theoryWeights[i+ind].wgt/theoryWeights[nominal].wgt;
      w_pdf.push_back(ipdfw);
      //if (ipdfw > 1e-3) pdfwgt[i] = ipdfw;
      pdfwgt[i] = ipdfw;
    }
    auto result = minmax_element(pdfwgt.begin(), pdfwgt.end());
    sys_pdf.push_back(*result.second); //max
    sys_pdf.push_back(*result.first);  //min
  } 
}
