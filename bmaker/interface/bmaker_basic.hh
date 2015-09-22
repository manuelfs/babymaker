// Creates baby tree with basic branches

#ifndef H_BMAKER_BASIC
#define H_BMAKER_BASIC

// System include files
#include <memory>
#include <vector>
#include <ctime>

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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

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
  time_t startTime;

  // Functions that do the branch writing
  vCands writeJets(edm::Handle<pat::JetCollection> alljets, vCands &leptons);
  void writeFatJets(vCands &jets);
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
		      vCands &jets);
  vCands writeMuons(edm::Handle<pat::MuonCollection> muons, 
		    edm::Handle<pat::PackedCandidateCollection> pfcands, 
		    edm::Handle<reco::VertexCollection> vtx);
  vCands writeElectrons(edm::Handle<pat::ElectronCollection> electrons, 
			edm::Handle<pat::PackedCandidateCollection> pfcands, 
			edm::Handle<reco::VertexCollection> vtx);
  void writeLeptons(vCands &leptons); 
  void writeTriggers(const edm::TriggerNames &names, 
                     edm::Handle<edm::TriggerResults> triggerBits,
                     edm::Handle<pat::PackedTriggerPrescales> triggerPrescales);
  void writeHLTObjects(const edm::TriggerNames &names, 
                       edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects);
  void writeFilters(const edm::TriggerNames &fnames,
                    edm::Handle<edm::TriggerResults> filterBits,
                    edm::Handle<reco::VertexCollection> vtx);
  void writeVertices(edm::Handle<reco::VertexCollection> vtx,
		     edm::Handle<std::vector< PileupSummaryInfo > >  pu_info);  
  void writeGenInfo(edm::Handle<LHEEventProduct> lhe_info);

  std::vector<TString> trig_name;

  // Input parameters
  TString outname;
  edm::InputTag met_label;
  unsigned int nevents_sample;
  unsigned int nevents;
  float xsec;

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
