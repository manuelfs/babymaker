// jet_met_tools: Functions related to jets, MET, and JECs

#ifndef H_JET_MET_TOOLS
#define H_JET_MET_TOOLS

// System include files
#include <vector>

// FW physics include files
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

// ROOT include files
#include "TString.h"

// User include files
#include "babymaker/bmaker/interface/utilities.hh"

class jet_met_tools{
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

  TString jecName;
  bool doJEC;
  double rhoEvent_;
  edm::Handle<pat::JetCollection> alljets_;
  FactorizedJetCorrectorCalculator *jetCorrector;
  FactorizedJetCorrectorCalculator::VariableValues jetValues;
  std::vector<float> jetTotCorrections, jetL1Corrections;
  std::vector<LVector> corrJet;

  bool leptonInJet(const pat::Jet &jet, vCands leptons);
  bool idJet(const pat::Jet &jet, bool doRa2=false);

  void getMETRaw(edm::Handle<pat::METCollection> mets, float &metRaw, float &metRawPhi);
  void getMETWithJEC(edm::Handle<pat::METCollection> mets, float &met, float &metPhi);
  void getJetCorrections(edm::Handle<pat::JetCollection> alljets, double rhoEvent);

  jet_met_tools(TString ijecName);
  ~jet_met_tools();
};



#endif
