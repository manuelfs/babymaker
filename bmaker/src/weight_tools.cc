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

  // the three following sets of weights are for the RPV analysis
  // calculated with
  // pileupCalc.py -i Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt --inputLumiJSON pileup_latest.txt
  // --calcMode true --minBiasXsec 65550 --maxPileupBin 52 --numPileupBins 52 pileup.root
  // and +/- 1 sigma variations of the cross section
  pileupWeightsPlus1Sigma = vector<double>({
    169.879, 122.327, 79.4771, 20.5781, 9.55734,
    1.07177, 0.532029, 0.976134, 2.02173, 2.62484,
    2.85851, 2.82154, 2.44658, 1.81064, 1.12149,
    0.570837, 0.239247, 0.0834243, 0.0252802, 0.00717547,
    0.00213205, 0.000732846, 0.000317964, 0.000182008, 0.000130839,
    0.000107646, 9.48789e-05, 8.61069e-05, 7.97246e-05, 7.30151e-05,
    6.05437e-05, 4.17578e-05, 2.43637e-05, 1.12016e-05, 4.84957e-06,
    1.90652e-06, 6.95957e-07, 2.41781e-07, 8.60729e-08, 2.90321e-08,
    9.33136e-09, 2.62465e-09, 8.4899e-10, 2.37178e-10, 6.69196e-11,
    1.79286e-11, 1.0299e-11, 3.86921e-12, 3.21265e-12});

  pileupWeightsCentral = vector<double>({
    199.584, 132.843, 82.5031, 23.4287, 10.512,
    1.37636, 0.943598, 1.98103, 3.29511, 3.5125,
    3.30143, 2.91134, 2.23547, 1.41984, 0.733959,
    0.306295, 0.10473, 0.0301507, 0.00784197, 0.00205396,
    0.000620814, 0.00024088, 0.000125514, 8.27626e-05, 6.18968e-05,
    4.8922e-05, 3.97662e-05, 3.26791e-05, 2.7146e-05, 2.21666e-05,
    1.63e-05, 9.91824e-06, 5.07907e-06, 2.0391e-06, 7.6692e-07,
    2.60585e-07, 8.1795e-08, 2.43096e-08, 7.36562e-09, 2.10371e-09,
    5.69632e-10, 1.3429e-10, 3.62225e-11, 8.39524e-12, 1.9552e-12,
    4.30121e-13, 2.01939e-13, 6.06721e-14, 5.20595e-14});

  pileupWeightsMinus1Sigma = vector<double>({
    233.452, 144, 87.8541, 26.7559, 11.6959,
    1.94776, 1.9294, 3.84918, 4.98223, 4.39835,
    3.59634, 2.80039, 1.84149, 0.967219, 0.403401,
    0.134488, 0.0370678, 0.00894529, 0.00211198, 0.000561504,
    0.000194316, 9.22467e-05, 5.54701e-05, 3.75286e-05, 2.65067e-05,
    1.90239e-05, 1.38035e-05, 1.00341e-05, 7.32377e-06, 5.22265e-06,
    3.3338e-06, 1.75049e-06, 7.68941e-07, 2.63236e-07, 8.39209e-08,
    2.4027e-08, 6.31717e-09, 1.56329e-09, 3.92062e-10, 9.21377e-11,
    2.04068e-11, 3.91181e-12, 8.52847e-13, 1.58856e-13, 2.96178e-14,
    5.18136e-15, 1.82869e-15, 0, 0});


}
weight_tools::~weight_tools() {}

float weight_tools::pileupWeight(unsigned int ntrupv)
{
  if(pileupWeights.size()!=0 && ntrupv<pileupWeights.size()) return static_cast<float>(pileupWeights.at(ntrupv));
  else return 1.0;
}

float weight_tools::pileupWeightRPV(unsigned int ntrupv, int type)
{
  if(type==-1) {
    if(pileupWeightsMinus1Sigma.size()!=0 && ntrupv<pileupWeightsMinus1Sigma.size()) return static_cast<float>(pileupWeightsMinus1Sigma.at(ntrupv));
  }
  else if(type==1) {
    if(pileupWeightsPlus1Sigma.size()!=0 && ntrupv<pileupWeightsPlus1Sigma.size()) return static_cast<float>(pileupWeightsPlus1Sigma.at(ntrupv));
  }
  else {
    if(pileupWeightsCentral.size()!=0 && ntrupv<pileupWeightsCentral.size()) return static_cast<float>(pileupWeightsCentral.at(ntrupv));
  }
  return 1.0;
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
      pdfwgt[i] = ipdfw;
    }
    auto result = minmax_element(pdfwgt.begin(), pdfwgt.end());
    sys_pdf.push_back(*result.second); //max
    sys_pdf.push_back(*result.first);  //min
  } 
}
