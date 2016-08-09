#include "FWCore/Framework/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "babymaker/bmaker/interface/weight_tools.hh"
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

weight_tools::weight_tools(){
  theoryWeights.clear();

  // obtained with scripts/get_puweights.py
  w_pu_up = vector<double>({1.207e-04, 8.232e-03, 1.234e-02, 1.980e-02, 
    3.640e-02, 2.697e-02, 2.949e-02, 4.402e-02, 1.032e-01, 2.170e-01, 3.557e-01, 
    5.363e-01, 7.440e-01, 9.311e-01, 1.058e+00, 1.040e+00, 9.937e-01, 1.084e+00, 
    1.084e+00, 1.200e+00, 1.098e+00, 1.050e+00, 1.095e+00, 1.171e+00, 1.165e+00, 
    1.303e+00, 1.214e+00, 1.271e+00, 1.128e+00, 1.083e+00, 9.853e-01, 8.872e-01, 
    1.039e+00, 1.103e+00, 1.620e+00, 2.816e+00, 6.427e+00, 9.780e+00, 0.000e+00, 
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00});

  w_pu_nom  = vector<double>({1.954e-04, 9.979e-03, 1.403e-02, 2.365e-02, 
    3.932e-02, 3.098e-02, 3.363e-02, 6.952e-02, 1.646e-01, 3.172e-01, 4.964e-01, 
    7.208e-01, 9.723e-01, 1.159e+00, 1.265e+00, 1.220e+00, 1.142e+00, 1.213e+00, 
    1.181e+00, 1.274e+00, 1.136e+00, 1.057e+00, 1.070e+00, 1.111e+00, 1.070e+00, 
    1.150e+00, 1.018e+00, 1.005e+00, 8.373e-01, 7.528e-01, 6.388e-01, 5.333e-01, 
    5.749e-01, 5.576e-01, 7.453e-01, 1.176e+00, 2.430e+00, 3.343e+00, 0.000e+00, 
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00});

  w_pu_down = vector<double>({3.109e-04, 1.224e-02, 1.587e-02, 2.821e-02, 
    4.299e-02, 3.549e-02, 4.146e-02, 1.192e-01, 2.551e-01, 4.628e-01, 6.873e-01, 
    9.691e-01, 1.250e+00, 1.416e+00, 1.503e+00, 1.420e+00, 1.293e+00, 1.334e+00, 
    1.262e+00, 1.323e+00, 1.145e+00, 1.031e+00, 1.011e+00, 1.012e+00, 9.314e-01, 
    9.452e-01, 7.836e-01, 7.209e-01, 5.582e-01, 4.641e-01, 3.613e-01, 2.744e-01, 
    2.673e-01, 2.334e-01, 2.801e-01, 3.957e-01, 7.322e-01, 9.037e-01, 0.000e+00, 
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 
    0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00});

}
weight_tools::~weight_tools() {}

float weight_tools::pileupWeight(unsigned int ntrupv_mean, int type)
{
  if(type==-1) {
    if(w_pu_down.size()!=0 && ntrupv_mean<w_pu_down.size()) return static_cast<float>(w_pu_down.at(ntrupv_mean));
  }
  else if(type==1) {
    if(w_pu_up.size()!=0 && ntrupv_mean<w_pu_up.size()) return static_cast<float>(w_pu_up.at(ntrupv_mean));
  }
  else {
    if(w_pu_nom.size()!=0 && ntrupv_mean<w_pu_nom.size()) return static_cast<float>(w_pu_nom.at(ntrupv_mean));
  }
  return 1.0;
}


float weight_tools::triggerEfficiency(int &nmus, int &nels, float &met, vector<float> &sys_trig){
  sys_trig.resize(2,1.);
  int nleps(nmus+nels);

  if(nleps == 0 || nleps > 2) return 1.;

  float eff_trig(1.);
  if(nleps == 1){
    if(nels==1){
      if(met>  0&&met<= 50) eff_trig = 0.897;
      if(met> 50&&met<=100) eff_trig = 0.923;
      if(met>100&&met<=150) eff_trig = 0.958;
      if(met>150&&met<=200) eff_trig = 0.988;
      // 1% systematic for met < 200
      sys_trig[0] = eff_trig+0.01;
      sys_trig[1] = min(1., eff_trig-0.01);
      if(met>200) {
        eff_trig = 0.998;
        sys_trig[0] = eff_trig+0.005;
        sys_trig[1] = min(1., eff_trig+0.005);
      }
    }
    if(nmus==1){
      if(met>  0&&met<= 50) eff_trig = 0.938;
      if(met> 50&&met<=100) eff_trig = 0.943;
      if(met>100&&met<=150) eff_trig = 0.965;
      if(met>150&&met<=200) eff_trig = 0.994;
      // 1% systematic for met < 200
      sys_trig[0] = eff_trig+0.01;
      sys_trig[1] = min(1., eff_trig-0.01);
      if(met>200){
        eff_trig = 0.997;
        sys_trig[0] = eff_trig+0.005;
        sys_trig[1] = min(1., eff_trig+0.005);
      }
    }
  } // nleps == 1

  if(nleps == 2){
    eff_trig = 1.;
    sys_trig[0] = 0.995;
    sys_trig[1] = 1.;
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
