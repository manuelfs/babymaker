// WEIGHT_TOOLS: Functions that deal with systematic weights

#ifndef H_WEIGHT_TOOLS
#define H_WEIGHT_TOOLS

class weight_tools{

private:
  std::vector<gen::WeightsInfo> weights;

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

  float weight(variationType variation);
  void getWeights(const edm::Event& iEvent);    

  weight_tools();
  ~weight_tools();
};

#endif
