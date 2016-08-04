// LEPTON_TOOLS: Functions to select analysis leptons

#ifndef H_LEPTON_TOOLS
#define H_LEPTON_TOOLS

#include "TH2D.h"
#include "TH3D.h"

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

  static double getScaleFactor(const reco::Candidate &cand);
  static double getScaleFactorUncertainty(const reco::Candidate &cand);

  static double getScaleFactor(const vCands &sig_leps);
  static double getScaleFactorUncertainty(const vCands &sig_leps);

  static double getScaleFactorFs(const reco::Candidate &cand, int npv);
  static double getScaleFactorUncertaintyFs(const reco::Candidate &cand, int npv);

  static double getScaleFactorFs(const vCands &sig_leps, int npv);
  static double getScaleFactorUncertaintyFs(const vCands &sig_leps, int npv);

  double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                        const reco::Candidate* ptcl,
                        double r_iso_min, double r_iso_max, double kt_scale,
                        double rho, bool charged_only);

  vCands getIsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi);
  vCands getRA4IsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi,double rhoEventCentral,std::vector<float> &isos,  std::vector<float> &relisos, int primary_pdg);

private:
  static const TH2D muon_id_sf, muon_iso_sf, electron_id_sf, electron_iso_sf;
  static const TH3D muon_idiso_fs_sf,electron_idiso_fs_sf;
  static double getScaleFactor(const reco::Muon &lep);
  static double getScaleFactorUncertainty(const reco::Muon &lep);
  static double getScaleFactor(const pat::Electron &lep);
  static double getScaleFactorUncertainty(const pat::Electron &lep);
  
  static double getScaleFactorFs(const reco::Muon &lep,  int npv);
  static double getScaleFactorUncertaintyFs(const reco::Muon &lep, int npv);
  static double getScaleFactorFs(const pat::Electron &lep,  int npv);
  static double getScaleFactorUncertaintyFs(const pat::Electron &lep, int npv);
};

#endif
