// Creates baby tree with basic branches

#ifndef H_BMAKER_FULL
#define H_BMAKER_FULL

// System include files
#include <memory>
#include <vector>
#include <ctime>
#include "babymaker/bmaker/interface/release.hh"

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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#ifdef POST_7_4
  #include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#endif

// User include files
#include "babymaker/bmaker/interface/baby_full.hh"
#include "babymaker/bmaker/interface/lepton_tools.hh"
#include "babymaker/bmaker/interface/photon_tools.hh"
#include "babymaker/bmaker/interface/event_tools.hh"
#include "babymaker/bmaker/interface/jet_met_tools.hh"
#include "babymaker/bmaker/interface/mc_tools.hh"
#include "babymaker/bmaker/interface/weight_tools.hh"
#include "babymaker/bmaker/interface/utilities.hh"

// ROOT include files
#include "TTree.h"
#include "TString.h"

typedef float& (baby_base::*baby_float)();
typedef std::vector<float>& (baby_base::*baby_vfloat)();
typedef std::vector<bool>& (baby_base::*baby_vbool)();

// Class declaration
class bmaker_full : public edm::EDAnalyzer {
public:
  explicit bmaker_full(const edm::ParameterSet&);
  ~bmaker_full();

  TFile *outfile;
  baby_full baby;
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
  void writeMET(edm::Handle<pat::METCollection> mets, edm::Handle<pat::METCollection> mets_nohf, edm::Handle<pat::METCollection> mets_uncorr, edm::Handle<pat::METCollection> mets_egclean, edm::Handle<pat::METCollection> mets_muclean);
  std::vector<reco::Candidate::LorentzVector> writeJets(edm::Handle<pat::JetCollection> alljets,
							std::vector<unsigned> &all_baby_jets_idx,
                 edm::Handle<edm::View <reco::GenJet> > genjets, 
                 vCands &sig_leps, vCands &veto_leps, 
                 vCands &photons, vCands &tks, 
                 std::vector<std::vector<LVector> > &sysjets,
                 std::vector<double> &jetsMuonEnergyFrac);
  void writeBTagWeights(edm::Handle<pat::JetCollection> alljets,
			std::vector<reco::Candidate::LorentzVector>  &all_baby_jets,
			std::vector<unsigned> &all_baby_jet_idx);
  void writeHiggVars(std::vector<LVector> &baby_jets_p4, std::vector<float> &baby_jets_csv, 
                     std::vector<bool> &baby_jets_h1, std::vector<bool> &baby_jets_h2, 
                     std::vector<bool> &baby_jets_islep, int &baby_nbl, int &baby_nbm, int &baby_nbt,
                     float &baby_hig_am, float &baby_hig_dm, float &baby_hig_drmax, 
                     int &baby_hig_bin, float &baby_mct, bool isSystemtatic = false);
  void writeBBVars(std::vector<reco::Candidate::LorentzVector>  &all_baby_jets, vCands &sig_leps);

  void writeFatJets();
  vCands writeMuons(edm::Handle<pat::MuonCollection> muons, 
		    edm::Handle<pat::PackedCandidateCollection> pfcands, 
		    edm::Handle<reco::VertexCollection> vtx,
		    vCands &veto_mus, vCands &all_mus, double rhoEventCentral);
  vCands writeElectrons(edm::Handle<pat::ElectronCollection> electrons, 
			edm::Handle<pat::ElectronCollection> electrons_pre_gs_fix,
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
  void writeTks(edm::Handle<pat::PackedCandidateCollection> pfcands,edm::Handle<reco::VertexCollection> vtx, double rhoEventCentral);
  vCands writePhotons(edm::Handle<pat::PhotonCollection> allphotons,edm::Handle<pat::PhotonCollection> allphotons_pre_gs_fix, edm::Handle<std::vector<pat::Electron> > &electrons,
		      edm::Handle<reco::ConversionCollection> &conversions, edm::Handle<reco::BeamSpot> &beamspot, double rho);

  bool writeTriggers(const edm::TriggerNames &names, 
                     edm::Handle<edm::TriggerResults> triggerBits,
                     edm::Handle<pat::PackedTriggerPrescales> triggerPrescales);
  void writeHLTObjects(const edm::TriggerNames &names, 
                       edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects,
		       vCands &all_mus, vCands &all_els, vCands &photons);
  void writeFilters(const edm::TriggerNames &fnames,
                    edm::Handle<edm::TriggerResults> filterBits,
                    edm::Handle<reco::VertexCollection> vtx,
                    std::vector<double> jetsMuonEnergyFrac);
  void writeVertices(edm::Handle<reco::VertexCollection> vtx,
		     edm::Handle<std::vector< PileupSummaryInfo > >  pu_info);  
  void writeGenInfo(edm::Handle<LHEEventProduct> lhe_info);
  void writeIFSR(edm::Handle<reco::GenParticleCollection> genParticles, 
                 std::vector<reco::Candidate::LorentzVector> &jets);
  void writeMC(edm::Handle<reco::GenParticleCollection> genParticles, vCands &all_mus, vCands &all_els, vCands &photons);

  void reportTime(const edm::Event& iEvent);

  // functions for calculating rebalanced MET
  void rebalancedMET(double& MET, double& METPhi);
  double calculateRescalingFactor(unsigned int jetIdx);
  double calculateRebalancedMET(unsigned int jetIdx, double mu, double& METPhi);

  // for filling additional event weights
  void writeWeights(const vCands &sig_leps, edm::Handle<GenEventInfoProduct> &gen_event_info, 
		    edm::Handle<LHEEventProduct> &lhe_info);

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
  bool debug;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  int mprod_;
  int mlsp_;

  // Tokens
  edm::EDGetTokenT<edm::TriggerResults> tok_trigResults_hlt_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> tok_patTrig_;
  edm::EDGetTokenT<reco::VertexCollection> tok_primVertex_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > > tok_addPileup_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo > > tok_slimAddPileup_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> tok_packPFCands_;
  edm::EDGetTokenT<double> tok_rhoFastJet_centralNeutral_;
  edm::EDGetTokenT<pat::MuonCollection> tok_muons_;
  edm::EDGetTokenT<pat::ElectronCollection> tok_electrons_;
  edm::EDGetTokenT<pat::ElectronCollection> tok_electrons_before_gsfix_;
  edm::EDGetTokenT<double> tok_rhoFastJet_all_;
  edm::EDGetTokenT<reco::BeamSpot> tok_offBeamSpot_;
  edm::EDGetTokenT<pat::PhotonCollection> tok_photons_;
  edm::EDGetTokenT<pat::PhotonCollection> tok_photons_before_gsfix_;
  edm::EDGetTokenT<std::vector<reco::Conversion> > tok_reducedEgamma_conver_;
  edm::EDGetTokenT<pat::JetCollection> tok_jets_;
  edm::EDGetTokenT<edm::View<reco::GenJet> > tok_genJets_;
  edm::EDGetTokenT<pat::METCollection> tok_met_;
  edm::EDGetTokenT<pat::METCollection> tok_met_noHF_;
  edm::EDGetTokenT<pat::METCollection> tok_met_uncorr_;
  edm::EDGetTokenT<pat::METCollection> tok_met_MuEGClean_;
  edm::EDGetTokenT<pat::METCollection> tok_met_EGClean_;
  edm::EDGetTokenT<bool> tok_HBHENoiseFilter_;
  edm::EDGetTokenT<bool> tok_HBHEIsoNoiseFilter_;
  edm::EDGetTokenT<edm::TriggerResults> tok_trigResults_reco_;
  edm::EDGetTokenT<edm::TriggerResults> tok_trigResults_pat_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> tok_selectedPatTrig_;
  edm::EDGetTokenT<reco::GenParticleCollection> tok_pruneGenPart_;
  edm::EDGetTokenT<LHEEventProduct> tok_extLHEProducer_;
  edm::EDGetTokenT<LHEEventProduct> tok_source_;
  edm::EDGetTokenT<GenEventInfoProduct> tok_generator_;
  edm::EDGetTokenT<bool> tok_badChCandFilter_;
  edm::EDGetTokenT<bool> tok_badPFMuonFilter_;
  #ifdef POST_7_4
    edm::EDGetTokenT<GenLumiInfoHeader> tok_genlumiheader_;
  #endif

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  // virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

#endif
