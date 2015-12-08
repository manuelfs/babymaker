// WEIGHT_TOOLS: Functions that deal with systematic weights

#ifndef H_WEIGHT_TOOLS
#define H_WEIGHT_TOOLS

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

class weight_tools{

private:
  std::vector<gen::WeightsInfo> theoryWeights;
  std::vector<double> pileupWeights;

public:
  // the enum index corresponds to the index of the variation
  enum variationType {
    nominal=0,
    muFup=1,
    muFdown=2,
    muRup=3,
    muRup_muFup=4,
    muRup_muFdown=5,
    muRdown=6,
    muRdown_muFup=7,
    muRdown_muFdown=8
  };

  float theoryWeight(variationType variation);
  void getTheoryWeights(edm::Handle<LHEEventProduct> lhe_info);
  float pileupWeight(unsigned int ntrupv);
  float triggerEfficiency(int &nmus, int &nels, std::vector<float> &sys_trig);
  float topPtWeight(float top_pt1,float top_pt2);
  float isrWeight(float isrpt);
  void getPDFWeights(std::vector<float> &sys_pdf, std::vector<float> &w_pdf);
  weight_tools();
  ~weight_tools();
};

#endif
