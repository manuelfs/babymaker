// Common physics objects definitions

#ifndef H_PHYS_OBJECTS
#define H_PHYS_OBJECTS

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "babymaker/bmaker/interface/utilities.hh"

namespace phys_objects{

  ///////////////// CUTS ///////////////////////
  float JetPtCut		= 30.0;
  float JetEtaCut               = 2.4;
  float CSVLoose                = 0.605;
  float CSVMedium               = 0.890;
  float CSVTight                = 0.970;

  float SignalLeptonPtCut	= 20.0;
  float VetoLeptonPtCut		= 10.0;
  float MuonEtaCut		= 2.4;
  float ElectronEtaCut		= 2.5;
  float MuonMiniIsoCut		= 0.2;
  float ElectronMiniIsoCut	= 0.1;


  enum CutLevel{kVeto, kLoose, kMedium, kTight};
  template<class T>
  T chooseVal(CutLevel threshold, T valVeto, T valLoose, T valMedium, T valTight){
   switch(threshold){
   default:
   case kVeto:
     return valVeto;
   case kLoose:
     return valLoose;
   case kMedium:
     return valMedium;
   case kTight:
     return valTight;
   }
   return valVeto;
 }

  bool isInJSON(std::string type, int run, int lumiblock);

  bool leptonInJet(const pat::Jet &jet, vCands leptons);
  bool idJet(const pat::Jet &jet);

  bool isVetoMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso);
  bool isSignalMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso);
  bool idMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold);
  bool vertexMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0);

  bool isVetoElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso);
  bool isSignalElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso);
  bool idElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold, bool do_iso=false);
  bool vertexElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0);

  double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			const reco::Candidate* ptcl,  
			double r_iso_min, double r_iso_max, double kt_scale,
			bool charged_only);

  bool hasGoodPV(edm::Handle<reco::VertexCollection> vtx);
}

#endif
