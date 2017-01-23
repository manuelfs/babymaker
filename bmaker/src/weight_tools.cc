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
  w_pu_down  = vector<double>({3.510e-01, 1.057e+00, 1.197e+00, 1.063e+00, 1.231e+00, 1.306e+00, 9.260e-01, 
    7.621e-01, 1.089e+00, 1.333e+00, 1.486e+00, 1.532e+00, 1.496e+00, 1.501e+00, 1.498e+00, 1.443e+00, 
    1.366e+00, 1.296e+00, 1.224e+00, 1.163e+00, 1.122e+00, 1.094e+00, 1.065e+00, 1.039e+00, 1.019e+00, 
    1.008e+00, 9.968e-01, 9.833e-01, 9.749e-01, 9.592e-01, 9.145e-01, 8.753e-01, 8.075e-01, 7.335e-01, 
    6.515e-01, 5.645e-01, 4.674e-01, 3.758e-01, 2.889e-01, 2.153e-01, 1.504e-01, 1.005e-01, 6.458e-02, 
    3.955e-02, 2.399e-02, 1.386e-02, 7.572e-03, 4.198e-03, 2.220e-03, 1.189e-03, 6.461e-04, 3.844e-04, 
    2.714e-04, 2.435e-04, 2.917e-04, 4.006e-04, 5.506e-04, 7.688e-04, 1.132e-03, 1.466e-03, 2.504e-03, 
    3.381e-03, 3.440e-03, 4.134e-03, 4.689e-03, 3.913e-03, 3.534e-03, 2.909e-03, 2.829e-03, 2.439e-03, 
    1.932e-03, 1.705e-03, 1.486e-03, 1.236e-03, 8.924e-04});

  w_pu_nom  = vector<double>({3.388e-01, 8.277e-01, 1.139e+00, 9.310e-01, 1.105e+00, 1.187e+00, 8.008e-01, 
    4.921e-01, 7.396e-01, 8.757e-01, 9.640e-01, 1.075e+00, 1.124e+00, 1.176e+00, 1.203e+00, 1.207e+00, 
    1.199e+00, 1.180e+00, 1.141e+00, 1.094e+00, 1.062e+00, 1.053e+00, 1.052e+00, 1.049e+00, 1.049e+00, 
    1.060e+00, 1.072e+00, 1.081e+00, 1.098e+00, 1.111e+00, 1.094e+00, 1.086e+00, 1.042e+00, 9.848e-01, 
    9.114e-01, 8.254e-01, 7.184e-01, 6.120e-01, 5.038e-01, 4.064e-01, 3.105e-01, 2.293e-01, 1.643e-01, 
    1.131e-01, 7.764e-02, 5.108e-02, 3.192e-02, 2.030e-02, 1.229e-02, 7.443e-03, 4.404e-03, 2.610e-03, 
    1.558e-03, 9.718e-04, 7.345e-04, 6.810e-04, 7.393e-04, 9.318e-04, 1.332e-03, 1.731e-03, 3.006e-03, 
    4.147e-03, 4.321e-03, 5.327e-03, 6.203e-03, 5.317e-03, 4.938e-03, 4.182e-03, 4.186e-03, 3.719e-03, 
    3.035e-03, 2.763e-03, 2.486e-03, 2.135e-03, 1.594e-03});

  w_pu_up  = vector<double>({3.301e-01, 6.518e-01, 1.078e+00, 8.173e-01, 1.005e+00, 1.069e+00, 7.303e-01, 
    3.452e-01, 4.988e-01, 6.007e-01, 6.320e-01, 7.344e-01, 8.268e-01, 9.128e-01, 9.610e-01, 9.885e-01, 
    1.024e+00, 1.051e+00, 1.048e+00, 1.025e+00, 1.002e+00, 9.996e-01, 1.015e+00, 1.036e+00, 1.057e+00, 
    1.087e+00, 1.120e+00, 1.153e+00, 1.195e+00, 1.235e+00, 1.246e+00, 1.272e+00, 1.260e+00, 1.233e+00, 
    1.182e+00, 1.111e+00, 1.005e+00, 8.924e-01, 7.702e-01, 6.564e-01, 5.348e-01, 4.251e-01, 3.310e-01, 
    2.497e-01, 1.894e-01, 1.386e-01, 9.691e-02, 6.935e-02, 4.744e-02, 3.256e-02, 2.182e-02, 1.454e-02, 
    9.534e-03, 6.171e-03, 4.309e-03, 3.109e-03, 2.269e-03, 1.888e-03, 1.988e-03, 2.191e-03, 3.561e-03, 
    4.837e-03, 5.075e-03, 6.358e-03, 7.554e-03, 6.620e-03, 6.292e-03, 5.458e-03, 5.599e-03, 5.100e-03, 
    4.272e-03, 3.993e-03, 3.691e-03, 3.258e-03, 2.501e-03});
}
weight_tools::~weight_tools() {}

float weight_tools::pileupWeight(unsigned int ntrupv_mean, int type)
{
  if (ntrupv_mean>=w_pu_down.size()) ntrupv_mean = w_pu_down.size()-1;

  if(type==-1) return static_cast<float>(w_pu_down.at(ntrupv_mean));
  else if(type==1) return static_cast<float>(w_pu_up.at(ntrupv_mean));
  else return static_cast<float>(w_pu_nom.at(ntrupv_mean));
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
      sys_trig[0] = min(1., eff_trig+0.01);
      sys_trig[1] = eff_trig-0.01;
      if(met>200) {
        eff_trig = 0.998;
        sys_trig[0] = min(1., eff_trig+0.005);
        sys_trig[1] = eff_trig-0.005;
      }
    }
    if(nmus==1){
      if(met>  0&&met<= 50) eff_trig = 0.938;
      if(met> 50&&met<=100) eff_trig = 0.943;
      if(met>100&&met<=150) eff_trig = 0.965;
      if(met>150&&met<=200) eff_trig = 0.994;
      // 1% systematic for met < 200
      sys_trig[0] = min(1., eff_trig+0.01);
      sys_trig[1] = eff_trig-0.01;
      if(met>200){
        eff_trig = 0.997;
        sys_trig[0] = min(1., eff_trig+0.005);
        sys_trig[1] = eff_trig-0.005;
      }
    }
  } // nleps == 1

  if(nleps == 2){
    eff_trig = 1.;
    sys_trig[0] = 1.;
    sys_trig[1] = 0.995;
  } // nleps == 2
  return eff_trig;
}

float weight_tools::topPtWeight(float toppt1, float toppt2){
  // 8TeV
  // float pt1,pt2;
  // if(toppt1>400) pt1=400;
  // else pt1=toppt1;
  // if(toppt2>400) pt2=400;
  // else pt2=toppt2;
  //  return sqrt(exp(0.156-0.00137*pt1)*exp(0.156-0.00137*pt2));
  // 13 TeV
  return sqrt(exp(0.0615-0.0005*toppt1)*exp(0.0615-0.0005*toppt2));
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
