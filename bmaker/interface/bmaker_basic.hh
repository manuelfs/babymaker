// Creates baby tree with basic branches

#ifndef H_BMAKER_BASIC
#define H_BMAKER_BASIC

// System include files
#include <memory>
#include <vector>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// FW physics include files
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

// ROOT include files
#include "TTree.h"
#include "TString.h"

// User include files
#include "babymaker/bmaker/interface/baby_basic.hh"
#include "babymaker/bmaker/interface/utilities.hh"

// Class declaration
class bmaker_basic : public edm::EDAnalyzer {
public:
  explicit bmaker_basic(const edm::ParameterSet&);
  ~bmaker_basic();

  TFile *outfile;
  baby_basic baby;
  bool isData;

  // Functions that do the branch writing
  void writeJets(baby_basic &baby, edm::Handle<pat::JetCollection> jets, vCands leptons);
  vCands writeMuons(baby_basic &baby, edm::Handle<pat::MuonCollection> muons, 
		    edm::Handle<pat::PackedCandidateCollection> pfcands, 
		    edm::Handle<reco::VertexCollection> vtx);
  vCands writeElectrons(baby_basic &baby, edm::Handle<pat::ElectronCollection> electrons, 
			edm::Handle<pat::PackedCandidateCollection> pfcands, 
			edm::Handle<reco::VertexCollection> vtx);
  void writeLeptons(baby_basic &baby, vCands leptons); 
  void writeTriggers(baby_basic &baby, const edm::TriggerNames &names, 
                     edm::Handle<edm::TriggerResults> triggerBits,
                     edm::Handle<pat::PackedTriggerPrescales> triggerPrescales);
  void writeHLTObjects(baby_basic &baby, const edm::TriggerNames &names, 
                       edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects);
  void writeFilters(baby_basic &baby, const edm::TriggerNames &fnames,
                    edm::Handle<edm::TriggerResults> filterBits,
                    edm::Handle<reco::VertexCollection> vtx);
  void writeVertices(baby_basic &baby, edm::Handle<reco::VertexCollection> vtx,
		     edm::Handle<std::vector< PileupSummaryInfo > >  pu_info);  

  std::vector<TString> trig_name;

  // Input parameters
  TString outname;
  edm::InputTag met_label;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

#endif
