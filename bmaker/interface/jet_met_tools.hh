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
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
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
  float getMCTagEfficiency(int pdgId, float pT, float eta, bool loose);
  
public:

  ///////////////// JET CUTS ///////////////////////
  const float JetPtCut	    = 30.0;
  const float JetPtCutLoose = 20.0;
  const float JetEtaCut    = 2.4;
  const float JetMHTEtaCut = 5.0;
  const float JetHLTPtCut  = 30.0;
  const float JetHLTEtaCut = 3.0;
  // Set in constructor depending on CMSSW release used
  float CSVLoose  = 999;
  float CSVMedium = 999;
  float CSVTight  = 999;

  const float sizeJet = 0.4;

  enum CutLevel{kLoose, kTight, kPBNR};
  enum btagVariation { kBTagCentral, kBTagUp, kBTagDown, kBTagCentralLoose, kBTagUpLoose, kBTagDownLoose};

  TString jecName;
  bool doSystematics;
  bool isFastSim;
  bool doJEC;
  bool isData;
  double rhoEvent_;
  edm::Handle<pat::JetCollection> alljets_;
  FactorizedJetCorrectorCalculator *jetCorrector;
  FactorizedJetCorrectorCalculator::VariableValues jetValues;
  JetCorrectionUncertainty *jecUncProvider;
  std::vector<float> jetTotCorrections, jetL1Corrections;
  std::vector<LVector> corrJet;
  std::vector<float> jerUnc, jecUnc;
  std::vector<float> genJetPt;

  BTagCalibration *calib;
  BTagCalibration *calibFS;
  std::vector<BTagCalibrationReader*> readersBC;
  std::vector<BTagCalibrationReader*> readersUDSG;
  std::vector<BTagCalibrationReader*> readersBC_fs;
  std::vector<BTagCalibrationReader*> readersUDSG_fs;
  TH3F *btagEfficiencyParameterization;
  TH3F *btagEfficiencyParameterizationLoose;

  bool leptonInJet(const pat::Jet &jet, vCands leptons);
  bool jetMatched(const pat::Jet &jet, vCands objects);
  bool idJet(const pat::Jet &jet, CutLevel cut);
  bool isLowDphi(vCands jets, float mht_phi, float &dphi1, float &dphi2, float &dphi3, float &dphi4);

  float getGenPt(const pat::Jet &jet, edm::Handle<edm::View <reco::GenJet> > genjets);
  float trueHT(edm::Handle<edm::View <reco::GenJet> > genjets);

  void getMETRaw(edm::Handle<pat::METCollection> mets, float &metRaw, float &metRawPhi);
  void getMETWithJEC(edm::Handle<pat::METCollection> mets, float &met, float &metPhi, unsigned isys);

  void getJetCorrections(edm::Handle<edm::View <reco::GenJet> > genjets, edm::Handle<pat::JetCollection> alljets, double rhoEvent);
  void setJetUncertainties(edm::Handle<edm::View <reco::GenJet> > genjets);

  float getJetResolutionSF(float jet_eta);
  float jetBTagWeight(const pat::Jet &jet, const LVector &jetp4, bool isBTagged, 
                      btagVariation readerTypeBC, btagVariation readerTypeUDSG,
                      btagVariation readerTypeBC_fs = kBTagCentral, btagVariation readerTypeUDSG_fs = kBTagCentral);
  void getDeltaRbb(std::vector<float> & deltaRbb, const std::vector<LVector> &jets, const std::vector<float> &jets_csv, const std::vector<bool> &jets_islep);
  static std::vector<size_t> getBRanking(const std::vector<LVector> &momentum, const std::vector<float> &csv,
					 const std::vector<bool> &is_lep);
  static float getDeltaRbb(const std::vector<LVector> &momentum, const std::vector<size_t> &brank, size_t nb);
  static float getDeltaRbbMax(const std::vector<LVector> &momentum, const std::vector<size_t> &brank,
			      size_t nb);
  static float getDeltaRbbMin(const std::vector<LVector> &momentum, const std::vector<size_t> &brank,
			      size_t nb);
  static float getDeltaPhibb(const std::vector<LVector> &momentum, const std::vector<size_t> &brank, size_t nb);
  static float getDeltaPhibbMax(const std::vector<LVector> &momentum, const std::vector<size_t> &brank,
				size_t nb);
  static float getDeltaPhibbMin(const std::vector<LVector> &momentum, const std::vector<size_t> &brank,
				size_t nb);
  static float getMbb(const std::vector<LVector> &momentum, const std::vector<size_t> &brank, size_t nb);
  static float getMbbMax(const std::vector<LVector> &momentum, const std::vector<size_t> &brank,
			 size_t nb);
  static float getMbbMin(const std::vector<LVector> &momentum, const std::vector<size_t> &brank,
			 size_t nb);
  static float getMblepMax2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			    size_t nb, const LVector &lep);
  static float getMblepMin2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			    size_t nb, const LVector &lep);
  static float getMblepMax(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			   size_t nb, const LVector &lep);
  static float getMblepMin(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			   size_t nb, const LVector &lep);
  static float getDeltaRblepMax2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				 size_t nb, const LVector &lep);
  static float getDeltaRblepMin2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				 size_t nb, const LVector &lep);
  static float getDeltaRblepMax(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				size_t nb, const LVector &lep);
  static float getDeltaRblepMin(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				size_t nb, const LVector &lep);
  static float getDeltaPhiblepMax2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				   size_t nb, const LVector &lep);
  static float getDeltaPhiblepMin2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				   size_t nb, const LVector &lep);
  static float getDeltaPhiblepMax(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				  size_t nb, const LVector &lep);
  static float getDeltaPhiblepMin(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				  size_t nb, const LVector &lep);
  static float getMTbmetMax2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			     size_t nb, float met, float met_phi);
  static float getMTbmetMin2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			     size_t nb, float met, float met_phi);
  static float getMTbmetMax(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			    size_t nb, float met, float met_phi);
  static float getMTbmetMin(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
			    size_t nb, float met, float met_phi);
  static float getDeltaPhibmetMax2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				   size_t nb, float met_phi);
  static float getDeltaPhibmetMin2(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				   size_t nb, float met_phi);
  static float getDeltaPhibmetMax(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				  size_t nb, float met_phi);
  static float getDeltaPhibmetMin(const std::vector<LVector> &jets, const std::vector<size_t> &brank,
				  size_t nb, float met_phi);

  void clusterFatJets(int &nfjets, float &mj,
          std::vector<float> &fjets_pt, 
          std::vector<float> &fjets_eta,
          std::vector<float> &fjets_phi, 
          std::vector<float> &fjets_m,
          std::vector<int> &fjets_nconst,
          std::vector<float> &fjets_sumcsv,
          std::vector<float> &fjets_poscsv,
          std::vector<int> &fjets_btags,
          std::vector<int> &jets_fjet_index,
          double radius,
          std::vector<LVector> &jets,
          std::vector<float> &jets_csv);
  double getSysMJ(double radius, std::vector<LVector> &jets);  
  jet_met_tools(TString ijecName, bool doSys, bool isFastSim);
  ~jet_met_tools();
};



#endif
