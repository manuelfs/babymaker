// jet_met_tools: Functions related to jets, MET, and JECs

#ifndef H_JET_MET_TOOLS
#define H_JET_MET_TOOLS

// System include files
#include <vector>

// FW physics include files
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
// B-tag scale factors                                                                                                                                                                 
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

// ROOT include files
#include "TString.h"
#include "TH3.h"

// User include files
#include "babymaker/bmaker/interface/utilities.hh"

class jet_met_tools{
private:
  float getMCTagEfficiency(int pdgId, float pT, float eta);

public:

  ///////////////// JET CUTS ///////////////////////
  const float JetPtCut	   = 30.0;
  const float JetEtaCut    = 2.4;
  const float JetMHTEtaCut = 5.0;
  const float JetHLTPtCut  = 40.0;
  const float JetHLTEtaCut = 3.0;
  const float CSVLoose     = 0.605;
  const float CSVMedium    = 0.890;
  const float CSVTight     = 0.970;

  const float sizeJet      = 0.4;

  enum CutLevel{kLoose, kTight, kPBNR};

  TString jecName;
  bool doJEC;
  double rhoEvent_;
  edm::Handle<pat::JetCollection> alljets_;
  FactorizedJetCorrectorCalculator *jetCorrector;
  FactorizedJetCorrectorCalculator::VariableValues jetValues;
  std::vector<float> jetTotCorrections, jetL1Corrections;
  std::vector<LVector> corrJet;

  BTagCalibration *calib;
  BTagCalibrationReader *readerBC;
  BTagCalibrationReader *readerUDSG;
  std::string variationTypeBC, variationTypeUDSG, btagEfficiencyFile;
  TH3F *btagEfficiencyParameterization;

  bool leptonInJet(const pat::Jet &jet, vCands leptons);
  bool jetMatched(const pat::Jet &jet, vCands objects);
  bool idJet(const pat::Jet &jet, CutLevel cut);
  bool isLowDphi(vCands jets, float mht_phi, float &dphi1, float &dphi2, float &dphi3, float &dphi4);

  float mismeasurement(const pat::Jet &jet, edm::Handle<edm::View <reco::GenJet> > genjets);
  float trueHT(edm::Handle<edm::View <reco::GenJet> > genjets);

  void getMETRaw(edm::Handle<pat::METCollection> mets, float &metRaw, float &metRawPhi);
  void getMETWithJEC(edm::Handle<pat::METCollection> mets, float &met, float &metPhi);
  void getJetCorrections(edm::Handle<pat::JetCollection> alljets, double rhoEvent);
  
  float jetBTagWeight(const pat::Jet &jet, const LVector &jetp4, bool isBTagged);

  jet_met_tools(TString ijecName, std::string btag_label_BC, std::string btag_label_UDSG, std::string btagEfficiency);
  ~jet_met_tools();
};



#endif
