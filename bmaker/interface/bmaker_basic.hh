// Creates baby tree with basic branches

#ifndef H_BMAKER_BASIC
#define H_BMAKER_BASIC

// System include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// FW physics include files
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

// ROOT include files
#include "TTree.h"

// User include files
#include "babymaker/bmaker/interface/baby_basic.hh"

// Class declaration

class bmaker_basic : public edm::EDAnalyzer {
public:
  explicit bmaker_basic(const edm::ParameterSet&);
  ~bmaker_basic();

  TFile * outfile;
  baby_basic baby;
  void WriteMuons(baby_basic &baby, edm::Handle<pat::MuonCollection> muons, 
		  edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<reco::VertexCollection> vtx);


  static float MinSignalLeptonPt, MinVetoLeptonPt;

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
