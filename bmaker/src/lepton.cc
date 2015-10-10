//// LEPTON: Lepton selection and isolation
//// Function names follow the first-lowercase, following words-uppercase. No underscores

// System include files
#include <algorithm>

// FW include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

// User include files
#include "babymaker/bmaker/interface/lepton.hh"

using namespace std;
using namespace utilities;

lepton::lepton():
  SignalLeptonPtCut(20.0),
  VetoLeptonPtCut(10.0),
  MuonEtaCut(2.4),
  ElectronEtaCut(2.5),
  MuonMiniIsoCut(0.2),
  ElectronMiniIsoCut(0.1){
}

lepton::~lepton(){ 
}

//////////////////// Muons
bool lepton::isSignalMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
  // pT, eta cuts
  if(lep.pt() <= SignalLeptonPtCut) return false;
  if(fabs(lep.eta()) > MuonEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idMuon(lep, vtx, kMedium)) return false;
  // Isolation cuts
  if(lep_iso >= 0 && lep_iso > MuonMiniIsoCut) return false;

  return true;
}

bool lepton::isVetoMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
  // pT, eta cuts
  if(lep.pt() <= VetoLeptonPtCut) return false;
  if(fabs(lep.eta()) > MuonEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  //if(!idMuon(lep, vtx, kLoose)) return false;
  if(!idMuon(lep, vtx, kMedium)) return false; // Using RA2/b muons
  // Isolation cuts
  if(lep_iso >= 0 && lep_iso > MuonMiniIsoCut) return false;

  return true;
}

bool lepton::idMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold) {
  double dz(0.), d0(0.);
  if(!vertexMuon(lep, vtx, dz, d0)) return false;

  switch(threshold){
  default:
  case kVeto:
  case kLoose:
    return lep.isLooseMuon();
  case kMedium:
    return lep.isMediumMuon();
  case kTight:
    return lep.isTightMuon(vtx->at(0));
  }
}
  
bool lepton::vertexMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
  dz = 0.; d0 = 0.;
  if(lep.track().isAvailable()){ // If the track is not available we probably don't want the muon
    dz = lep.track()->vz()-vtx->at(0).z();
    d0 = lep.track()->d0()-vtx->at(0).x()*sin(lep.track()->phi())+vtx->at(0).y()*cos(lep.track()->phi());
  } 
  if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

  return true;
}

double lepton::getEffAreaMuon(double eta){
  double abseta = fabs(eta);
  if (abseta < 0.8) return 0.0735;
  else if (abseta < 1.3) return 0.0619;
  else if (abseta < 2.0) return 0.0465;
  else if (abseta < 2.2) return 0.0433;
  else if (abseta < 2.5) return 0.0577;
  else return 0;
}

double lepton::getRelIsolation(const pat::Muon &lep, double rho){
  double ch_iso(lep.pfIsolationR04().sumChargedHadronPt);
  double neu_iso(max(0., lep.pfIsolationR04().sumNeutralHadronEt + lep.pfIsolationR04().sumPhotonEt
	       -rho*getEffAreaMuon(lep.eta())));

  return (ch_iso + neu_iso) / lep.pt();
}


//////////////////// Electrons
bool lepton::isSignalElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
  // pT, eta cuts
  if(lep.pt() <= SignalLeptonPtCut) return false;
  if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idElectron(lep, vtx, kMedium)) return false;
  // Isolation cuts
  if(lep_iso >= 0 && lep_iso > ElectronMiniIsoCut) return false;

  return true;
}

bool lepton::isVetoElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
  // pT, eta cuts
  if(lep.pt() <= VetoLeptonPtCut) return false;
  if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idElectron(lep, vtx, kVeto)) return false;
  // Isolation cuts
  if(lep_iso >= 0 && lep_iso > ElectronMiniIsoCut) return false;

  return true;
}

bool lepton::idElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold, bool do_iso) {

  bool barrel(lep.isEB());
  double deta_cut, dphi_cut, ieta_cut, hovere_cut, d0_cut, dz_cut,
    ooeminusoop_cut, reliso_cut, misshits_cut;
  bool req_conv_veto;

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Working_points_for_Spring15_MC_s
  // Last updated October 8th
  if(barrel){
    ieta_cut        = chooseVal(threshold	,0.0114,  0.0103,  0.0101,  0.0101);
    deta_cut        = chooseVal(threshold	,0.0152,  0.0105,  0.0103,  0.00926);
    dphi_cut        = chooseVal(threshold	,0.216,   0.115,   0.0336,  0.0336);
    hovere_cut      = chooseVal(threshold	,0.181,   0.104,   0.0876,  0.0597);
    reliso_cut      = chooseVal(threshold	,0.126,   0.0893,  0.0766,  0.0354);
    ooeminusoop_cut = chooseVal(threshold	,0.207,   0.102,   0.0174,  0.012);
    d0_cut          = chooseVal(threshold	,0.0564,  0.0261,  0.0118,  0.0111);
    dz_cut          = chooseVal(threshold	,0.472,   0.41,  0.373,   0.0466);
    misshits_cut    = chooseVal(threshold	,2,   2,   2,   2);
    req_conv_veto   = chooseVal(threshold	,true		,  true		,  true		,  true );
  } else {
    ieta_cut        = chooseVal(threshold	,0.0352 , 0.0301 , 0.0283 , 0.0279);
    deta_cut        = chooseVal(threshold	,0.0113 , 0.00814 , 0.00733 , 0.00724);
    dphi_cut        = chooseVal(threshold	,0.237 , 0.182 , 0.114 , 0.0918);
    hovere_cut      = chooseVal(threshold	,0.116 , 0.0897 , 0.0678 , 0.0615);
    reliso_cut      = chooseVal(threshold	,0.144 , 0.121 , 0.0678 , 0.0646);
    ooeminusoop_cut = chooseVal(threshold	,0.174 , 0.126 , 0.0898 , 0.00999);
    d0_cut          = chooseVal(threshold	,0.222 , 0.118 , 0.0739 , 0.0351);
    dz_cut          = chooseVal(threshold	,0.921 , 0.822 , 0.602 , 0.417);
    misshits_cut    = chooseVal(threshold	,3, 1, 1, 1);
    req_conv_veto   = chooseVal(threshold	,true   ,  true   ,  true   ,  true );
  }


  double dz(0.), d0(0.);
  vertexElectron(lep, vtx, dz, d0);

  int mhits(0);
  if(lep.gsfTrack().isAvailable()){
    mhits = lep.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);;
  } 
  const double sigietaieta(lep.full5x5_sigmaIetaIeta());

  return deta_cut > fabs(lep.deltaEtaSuperClusterTrackAtVtx())
    && dphi_cut > fabs(lep.deltaPhiSuperClusterTrackAtVtx())
    && ieta_cut > sigietaieta
    && hovere_cut > lep.hadronicOverEm()
    && d0_cut > fabs(d0)
    && dz_cut > fabs(dz)
    && ooeminusoop_cut > fabs((1.0-lep.eSuperClusterOverP())/lep.ecalEnergy())
    && (!do_iso || reliso_cut > 0) // To be implemented if we want reliso
    && (!req_conv_veto || lep.passConversionVeto())
    && (misshits_cut >= mhits);
}

bool lepton::vertexElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
  dz = 0.; d0 = 0.;
  if(lep.gsfTrack().isAvailable()){ // If the track is not available we probably don't want the electron
    dz = lep.gsfTrack()->vz()-vtx->at(0).z();
    d0 = lep.gsfTrack()->d0()-vtx->at(0).x()*sin(lep.gsfTrack()->phi())+vtx->at(0).y()*cos(lep.gsfTrack()->phi());
  } 
  if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

  return true;
}

double lepton::getEffAreaElectron(double eta){
  double abseta = fabs(eta);
  if (abseta < 1) return 0.1752;
  else if (abseta < 1.479) return 0.1862;
  else if (abseta < 2.0) return 0.1411;
  else if (abseta < 2.2) return 0.1534;
  else if (abseta < 2.3) return 0.1903;
  else if (abseta < 2.4) return 0.2243;
  else if (abseta < 2.5) return 0.2687;
  else return 0;
}

double lepton::getRelIsolation(const pat::Electron &lep, double rho){
  double ch_iso(lep.pfIsolationVariables().sumChargedHadronPt);
  double neu_iso(max(0., lep.pfIsolationVariables().sumNeutralHadronEt + lep.pfIsolationVariables().sumPhotonEt
	       -rho*getEffAreaElectron(lep.eta())));

  return (ch_iso + neu_iso) / lep.pt();
}

double lepton::getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                      const reco::Candidate* ptcl,  
                      double r_iso_min, double r_iso_max, double kt_scale,
                      double rho, bool charged_only) {
  if (ptcl->pt()<5.) return 99999.;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  if(ptcl->isElectron()) {
    if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  } else if(ptcl->isMuon()) {
    deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
  } else {
    //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
  }
  double iso_nh(0.), iso_ch(0.), iso_ph(0.), iso_pu(0.);
  double ptThresh(0.5), r_mini(kt_scale/ptcl->pt());
  if(ptcl->isElectron()) ptThresh = 0;
  double r_iso(max(r_iso_min, min(r_iso_max, r_mini)));

  for (const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;
    double dr = deltaR(pfc, *ptcl);
    if (dr > r_iso) continue;
    if (pfc.charge()==0){ //neutrals
      if (pfc.pt()>ptThresh) {
        if (abs(pfc.pdgId())==22) { //photons
          if(dr < deadcone_ph) continue;
          iso_ph += pfc.pt();
        } else if (abs(pfc.pdgId())==130) { //neutral hadrons
          if(dr < deadcone_nh) continue;
          iso_nh += pfc.pt();
        }
      }
    } else if (pfc.fromPV()>1){ //charged from PV
      if (abs(pfc.pdgId())==211) {
        if(dr < deadcone_ch) continue;
        iso_ch += pfc.pt();
      }
    } else {
      if (pfc.pt()>ptThresh){ //charged from PU
        if(dr < deadcone_pu) continue;
        iso_pu += pfc.pt();
      }
    }
  } // Loop over pf cands

  double effarea = ptcl->isElectron() ? getEffAreaElectron(ptcl->eta()) : getEffAreaMuon(ptcl->eta());
  double pu_corr = rho*effarea*pow(r_mini,2)/(0.3*0.3);
  double iso(0.);
  if (charged_only) iso = iso_ch;
  else iso = iso_ch + max(0.,iso_ph + iso_nh - pu_corr);

  return iso/ptcl->pt();
}
