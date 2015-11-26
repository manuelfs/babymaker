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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
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
#include "babymaker/bmaker/interface/lepton_tools.hh"
#include "babymaker/bmaker/interface/photon_tools.hh"
#include "babymaker/bmaker/interface/event_tools.hh"
#include "babymaker/bmaker/interface/jet_met_tools.hh"
#include "babymaker/bmaker/interface/mc_tools.hh"
#include "babymaker/bmaker/interface/weight_tools.hh"
#include "babymaker/bmaker/interface/utilities.hh"

typedef float& (baby_base::*baby_float)();
typedef std::vector<float>& (baby_base::*baby_vfloat)();
typedef std::vector<bool>& (baby_base::*baby_vbool)();

// Class declaration
class bmaker_basic : public edm::EDAnalyzer {
public:
  explicit bmaker_basic(const edm::ParameterSet&);
  ~bmaker_basic();

  TFile *outfile;
  baby_basic baby;
  bool isData;
  time_t startTime;

  //object classes
  lepton_tools *lepTool;
  jet_met_tools *jetTool;
  photon_tools *photonTool;
  event_tools *eventTool;
  mc_tools *mcTool;
  weight_tools *weightTool;

  // Functions that do the branch writing
  void writeMET(edm::Handle<pat::METCollection> mets, edm::Handle<pat::METCollection> mets_nohf);
  void writeJets(edm::Handle<pat::JetCollection> alljets, 
                 edm::Handle<edm::View <reco::GenJet> > genjets, 
                 vCands &sig_leps, vCands &veto_leps, 
                 vCands &photons, vCands &tks, 
                 std::vector<LVector> &jets, 
                 std::vector<std::vector<LVector> > &sysjets);
  void writeFatJets(std::vector<LVector> &jets);
  vCands writeMuons(edm::Handle<pat::MuonCollection> muons, 
		    edm::Handle<pat::PackedCandidateCollection> pfcands, 
		    edm::Handle<reco::VertexCollection> vtx,
		    vCands &veto_mus, vCands &all_mus, double rhoEventCentral);
  vCands writeElectrons(edm::Handle<pat::ElectronCollection> electrons, 
			edm::Handle<pat::PackedCandidateCollection> pfcands, 
			edm::Handle<reco::VertexCollection> vtx,
			vCands &veto_els, vCands &all_els, double rhoEventCentral);
  void writeDiLep(vCands &sig_mus, vCands &sig_els, vCands &veto_mus, vCands &veto_els);
  void setDiLepMass(vCands leptons, baby_float ll_m, baby_float ll_pt1, baby_float ll_pt2, baby_float ll_pt, 
		    baby_float ll_eta, baby_float ll_phi, baby_vfloat l_pt, baby_vbool l_inz,
		    baby_float ll_w);
  void setElMuMass(vCands leptons1, vCands leptons2, baby_float ll_m, baby_float ll_pt1, baby_float ll_pt2, 
		   baby_float ll_pt, baby_float ll_eta, baby_float ll_phi,
		   baby_float ll_w);
  void writeLeptons(vCands &leptons); 

  vCands writePhotons(edm::Handle<pat::PhotonCollection> allphotons, edm::Handle<std::vector<pat::Electron> > &electrons,
		      edm::Handle<reco::ConversionCollection> &conversions, edm::Handle<reco::BeamSpot> &beamspot, double rho);

  bool writeTriggers(const edm::TriggerNames &names, 
                     edm::Handle<edm::TriggerResults> triggerBits,
                     edm::Handle<pat::PackedTriggerPrescales> triggerPrescales);
  void writeHLTObjects(const edm::TriggerNames &names, 
                       edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects,
		       vCands &all_mus, vCands &all_els, vCands &photons);
  void writeFilters(const edm::TriggerNames &fnames,
                    edm::Handle<edm::TriggerResults> filterBits,
                    edm::Handle<reco::VertexCollection> vtx);
  void writeVertices(edm::Handle<reco::VertexCollection> vtx,
		     edm::Handle<std::vector< PileupSummaryInfo > >  pu_info);  
  void writeGenInfo(edm::Handle<LHEEventProduct> lhe_info);
  void writeMC(edm::Handle<reco::GenParticleCollection> genParticles, vCands &all_mus, vCands &all_els, vCands &photons);

  void reportTime(const edm::Event& iEvent);

  // functions for calculating rebalanced MET
  void rebalancedMET(double& MET, double& METPhi);
  double calculateRescalingFactor(unsigned int jetIdx);
  double calculateRebalancedMET(unsigned int jetIdx, double mu, double& METPhi);

  // for filling additional event weights
  void fillWeights(const vCands &sig_leps);

  std::vector<TString> trig_name;

  // Input parameters
  TString outname;
  std::vector<std::string> inputfiles;
  TString jsonfile;
  TString condor_subtime;
  TString jec_label;
  std::string btag_label_BC;
  std::string btag_label_UDSG;
  edm::InputTag met_label;
  edm::InputTag met_nohf_label;
  edm::InputTag jets_label;
  unsigned int nevents_sample;
  unsigned int nevents;
  bool doMetRebalancing;
  float xsec;

  bool addBTagWeights;
  bool isFastSim;
  bool doSystematics;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  // virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

#endif
