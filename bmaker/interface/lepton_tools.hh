// LEPTON_TOOLS: Functions to select analysis leptons

#ifndef H_LEPTON_TOOLS
#define H_LEPTON_TOOLS

#include <utility>

#include "TH2D.h"
#include "TH2F.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "babymaker/bmaker/interface/utilities.hh"

class lepton_tools{

public:
  lepton_tools();
  ~lepton_tools();

  ///////////////// LEPTON CUTS ///////////////////////
  const float SignalLeptonPtCut  = 20.0;
  const float VetoLeptonPtCut    = 10.0;
  const float MuonEtaCut         = 2.4;
  const float ElectronEtaCut     = 2.5;
  const float MuonMiniIsoCut     = 0.2;
  const float ElectronMiniIsoCut = 0.1;

  enum CutLevel{kVeto, kLoose, kMedium, kMediumICHEP, kTight};
  template<class T>
  T chooseVal(CutLevel threshold, T valVeto, T valLoose, T valMedium, T valTight){
    switch(threshold){
    default:
    case kMediumICHEP:
      cms::Exception("BadElectronID")
	<< "CutLevel " << static_cast<int>(threshold) << " is not accepted for electrons." << std::endl;
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

  bool isVetoMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso);
  bool isSignalMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso);
  bool idMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold);
  bool vertexMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0);
  double getEffAreaMuon(double eta);
  double getRelIsolation(const pat::Muon &lep, double rho);

  bool isVetoElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso);
  bool isSignalElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso);
  bool idElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold, bool doIso=false);
  bool vertexElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0);
  double getEffAreaElectron(double eta);
  double getRelIsolation(const pat::Electron &lep, double rho);

  static std::pair<double, double> getScaleFactor(const reco::Candidate &cand);
  static std::pair<double, double> getScaleFactor(const vCands &sig_leps);

  static std::pair<double, double> getScaleFactorFs(const reco::Candidate &cand);
  static std::pair<double, double> getScaleFactorFs(const vCands &sig_leps);

  double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                        const reco::Candidate* ptcl,
                        double r_iso_min, double r_iso_max, double kt_scale,
                        double rho, bool charged_only);

  vCands getIsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi);
  vCands getRA4IsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi,double rhoEventCentral,std::vector<float> &isos,  std::vector<float> &relisos, int primary_pdg);

private:
  static const TH2F sf_full_muon_medium, sf_full_muon_iso, sf_full_electron_medium, sf_full_electron_iso, sf_full_electron_tracking;
  static const TH2D sf_full_muon_tracking, sf_fast_muon_medium, sf_fast_muon_iso, sf_fast_electron_mediumiso;

  static std::pair<double, double> getScaleFactor(const reco::Muon &lep);
  static std::pair<double, double> getScaleFactor(const pat::Electron &lep);
  
  static std::pair<double, double> getScaleFactorFs(const reco::Muon &lep);
  static std::pair<double, double> getScaleFactorFs(const pat::Electron &lep);
};

#endif
