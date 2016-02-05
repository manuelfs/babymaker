//// LEPTON_TOOLS: Lepton selection and isolation
//// Function names follow the first-lowercase, following words-uppercase. No underscores

// System include files
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <stdexcept>
#include <string>

//ROOT include files
#include "TFile.h"

// FW include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

// User include files
#include "babymaker/bmaker/interface/lepton_tools.hh"

using namespace std;
using namespace utilities;

//////////////////// Scale Factor loading
const TH2D lepton_tools::muon_id_sf = *static_cast<TH2D*>(TFile((string(getenv("CMSSW_BASE"))+"/src/babymaker/bmaker/data/lepton_sf/muon_medium_id.root").c_str(),"read").Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_tag_IsoMu20_pass"));
const TH2D lepton_tools::muon_iso_sf = *static_cast<TH2D*>(TFile((string(getenv("CMSSW_BASE"))+"/src/babymaker/bmaker/data/lepton_sf/muon_mini_iso_0p2.root").c_str(),"read").Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_tag_combRelIsoPF04dBeta_bin0_&_tag_pt_bin0_&_PF_pass_&_tag_IsoMu20_pass"));
const TH2D lepton_tools::electron_id_sf = *static_cast<TH2D*>(TFile((string(getenv("CMSSW_BASE"))+"/src/babymaker/bmaker/data/lepton_sf/electron.root").c_str(),"read").Get("CutBasedMedium"));
const TH2D lepton_tools::electron_iso_sf = *static_cast<TH2D*>(TFile((string(getenv("CMSSW_BASE"))+"/src/babymaker/bmaker/data/lepton_sf/electron.root").c_str(),"read").Get("MiniIso0p1_vs_AbsEta"));
const TH3D lepton_tools::muon_idiso_fs_sf = *static_cast<TH3D*>(TFile((string(getenv("CMSSW_BASE"))+"/src/babymaker/bmaker/data/lepton_sf/sf_mu_mediumID_mini02.root").c_str(),"read").Get("histo3D"));
const TH3D lepton_tools::electron_idiso_fs_sf = *static_cast<TH3D*>(TFile((string(getenv("CMSSW_BASE"))+"/src/babymaker/bmaker/data/lepton_sf/sf_el_mediumCB_mini01.root").c_str(),"read").Get("histo3D"));

//////////////////// Muons
bool lepton_tools::isSignalMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= SignalLeptonPtCut) return false;
  if(fabs(lep.eta()) > MuonEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idMuon(lep, vtx, kMedium)) return false;
  // Isolation cuts
  if(lepIso >= 0 && lepIso > MuonMiniIsoCut) return false;

  return true;
}

bool lepton_tools::isVetoMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= VetoLeptonPtCut) return false;
  if(fabs(lep.eta()) > MuonEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  //if(!idMuon(lep, vtx, kLoose)) return false;
  if(!idMuon(lep, vtx, kMedium)) return false; // Using RA2/b muons
  // Isolation cuts
  if(lepIso >= 0 && lepIso > MuonMiniIsoCut) return false;

  return true;
}

bool lepton_tools::idMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold) {
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

bool lepton_tools::vertexMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
  dz = 0.; d0 = 0.;
  if(lep.track().isAvailable()){ // If the track is not available we probably don't want the muon
    dz = lep.track()->vz()-vtx->at(0).z();
    d0 = lep.track()->d0()-vtx->at(0).x()*sin(lep.track()->phi())+vtx->at(0).y()*cos(lep.track()->phi());
  }
  if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

  return true;
}

double lepton_tools::getEffAreaMuon(double eta){
  double abseta = fabs(eta);
  if (abseta < 0.8) return 0.0735;
  else if (abseta < 1.3) return 0.0619;
  else if (abseta < 2.0) return 0.0465;
  else if (abseta < 2.2) return 0.0433;
  else if (abseta < 2.5) return 0.0577;
  else return 0;
}

double lepton_tools::getRelIsolation(const pat::Muon &lep, double rho){
  double ch_iso(lep.pfIsolationR04().sumChargedHadronPt);
  double neu_iso(max(0., lep.pfIsolationR04().sumNeutralHadronEt + lep.pfIsolationR04().sumPhotonEt
                     -rho*getEffAreaMuon(lep.eta())));

  return (ch_iso + neu_iso) / lep.pt();
}

//////////////////// Electrons
bool lepton_tools::isSignalElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= SignalLeptonPtCut) return false;
  if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idElectron(lep, vtx, kMedium)) return false;
  // Isolation cuts
  if(lepIso >= 0 && lepIso > ElectronMiniIsoCut) return false;

  return true;
}

bool lepton_tools::isVetoElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= VetoLeptonPtCut) return false;
  if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idElectron(lep, vtx, kVeto)) return false;
  // Isolation cuts
  if(lepIso >= 0 && lepIso > ElectronMiniIsoCut) return false;

  return true;
}

bool lepton_tools::idElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold, bool doIso) {

  bool barrel(lep.isEB());
  double deta_cut, dphi_cut, ieta_cut, hovere_cut, d0_cut, dz_cut,
    ooeminusoop_cut, reliso_cut, misshits_cut;
  bool req_conv_veto;

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Working_points_for_Spring15_MC_s
  // Last updated October 8th
  if(barrel){
    ieta_cut        = chooseVal(threshold       ,0.0114,  0.0103,  0.0101,  0.0101);
    deta_cut        = chooseVal(threshold       ,0.0152,  0.0105,  0.0103,  0.00926);
    dphi_cut        = chooseVal(threshold       ,0.216,   0.115,   0.0336,  0.0336);
    hovere_cut      = chooseVal(threshold       ,0.181,   0.104,   0.0876,  0.0597);
    reliso_cut      = chooseVal(threshold       ,0.126,   0.0893,  0.0766,  0.0354);
    ooeminusoop_cut = chooseVal(threshold       ,0.207,   0.102,   0.0174,  0.012);
    d0_cut          = chooseVal(threshold       ,0.0564,  0.0261,  0.0118,  0.0111);
    dz_cut          = chooseVal(threshold       ,0.472,   0.41,  0.373,   0.0466);
    misshits_cut    = chooseVal(threshold       ,2,   2,   2,   2);
    req_conv_veto   = chooseVal(threshold       ,true           ,  true         ,  true         ,  true );
  } else {
    ieta_cut        = chooseVal(threshold       ,0.0352 , 0.0301 , 0.0283 , 0.0279);
    deta_cut        = chooseVal(threshold       ,0.0113 , 0.00814 , 0.00733 , 0.00724);
    dphi_cut        = chooseVal(threshold       ,0.237 , 0.182 , 0.114 , 0.0918);
    hovere_cut      = chooseVal(threshold       ,0.116 , 0.0897 , 0.0678 , 0.0615);
    reliso_cut      = chooseVal(threshold       ,0.144 , 0.121 , 0.0678 , 0.0646);
    ooeminusoop_cut = chooseVal(threshold       ,0.174 , 0.126 , 0.0898 , 0.00999);
    d0_cut          = chooseVal(threshold       ,0.222 , 0.118 , 0.0739 , 0.0351);
    dz_cut          = chooseVal(threshold       ,0.921 , 0.822 , 0.602 , 0.417);
    misshits_cut    = chooseVal(threshold       ,3, 1, 1, 1);
    req_conv_veto   = chooseVal(threshold       ,true   ,  true   ,  true   ,  true );
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
    && (!doIso || reliso_cut > 0) // To be implemented if we want reliso
    && (!req_conv_veto || lep.passConversionVeto())
    && (misshits_cut >= mhits);
}

bool lepton_tools::vertexElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
  dz = 0.; d0 = 0.;
  if(lep.gsfTrack().isAvailable()){ // If the track is not available we probably don't want the electron
    dz = lep.gsfTrack()->vz()-vtx->at(0).z();
    d0 = lep.gsfTrack()->d0()-vtx->at(0).x()*sin(lep.gsfTrack()->phi())+vtx->at(0).y()*cos(lep.gsfTrack()->phi());
  }
  if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

  return true;
}

double lepton_tools::getEffAreaElectron(double eta){
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

double lepton_tools::getRelIsolation(const pat::Electron &lep, double rho){
  double ch_iso(lep.pfIsolationVariables().sumChargedHadronPt);
  double neu_iso(max(0., lep.pfIsolationVariables().sumNeutralHadronEt + lep.pfIsolationVariables().sumPhotonEt
                     -rho*getEffAreaElectron(lep.eta())));

  return (ch_iso + neu_iso) / lep.pt();
}

double lepton_tools::getScaleFactor(const reco::Candidate &cand){
  if(const pat::Electron * ele_ptr = dynamic_cast<const pat::Electron*>(&cand)){
    return getScaleFactor(*ele_ptr);
  }else if(const reco::Muon * mu_ptr = dynamic_cast<const reco::Muon*>(&cand)){
    return getScaleFactor(*mu_ptr);
  }else{
    throw runtime_error(string("Cannot get scale factor for type ")+typeid(cand).name());
  }
  return 1.;
}

double lepton_tools::getScaleFactorUncertainty(const reco::Candidate &cand){
  if(const pat::Electron * ele_ptr = dynamic_cast<const pat::Electron*>(&cand)){
    return getScaleFactorUncertainty(*ele_ptr);
  }else if(const reco::Muon * mu_ptr = dynamic_cast<const reco::Muon*>(&cand)){
    return getScaleFactorUncertainty(*mu_ptr);
  }else{
    throw runtime_error(string("Cannot get scale factor uncertainty for type ")+typeid(cand).name());
  }
  return 0.;
}

double lepton_tools::getScaleFactor(const vCands &sig_leps){
  double scale_factor = 1.;
  for(const auto &lep: sig_leps){
    if(lep == nullptr) throw runtime_error("sig_leps contains a nullptr in lepton_tools::getScaleFactor");
    scale_factor *= getScaleFactor(*lep);
  }
  return scale_factor;
}

double lepton_tools::getScaleFactorUncertainty(const vCands &sig_leps){
  //Crashes if scale factor == 0...
  double uncertainty = 0., totsf(1.);
  for(const auto &lep: sig_leps){
    if(lep == nullptr) throw runtime_error("sig_leps contains a nullptr in lepton_tools::getScaleFactorUncertainty");
    double scale_factor = getScaleFactor(*lep);
    if(scale_factor <= 0) throw runtime_error("One lepton FSSF is 0 in lepton_tools::getScaleFactorUncertainty");
    uncertainty = hypot(uncertainty, getScaleFactorUncertainty(*lep)/scale_factor);
    totsf *= scale_factor;
  }
  return uncertainty*totsf;
}

double lepton_tools::getScaleFactorFs(const reco::Candidate &cand, int npv){
  if(const pat::Electron * ele_ptr = dynamic_cast<const pat::Electron*>(&cand)){
    return getScaleFactorFs(*ele_ptr,npv);
  }else if(const reco::Muon * mu_ptr = dynamic_cast<const reco::Muon*>(&cand)){
    return getScaleFactorFs(*mu_ptr,npv);
  }else{
    throw runtime_error(string("Cannot get scale factor for type ")+typeid(cand).name());
  }
  return 1.;
}

double lepton_tools::getScaleFactorUncertaintyFs(const reco::Candidate &cand, int npv){
  if(const pat::Electron * ele_ptr = dynamic_cast<const pat::Electron*>(&cand)){
    return getScaleFactorUncertaintyFs(*ele_ptr,npv);
  }else if(const reco::Muon * mu_ptr = dynamic_cast<const reco::Muon*>(&cand)){
    return getScaleFactorUncertaintyFs(*mu_ptr,npv);
  }else{
    throw runtime_error(string("Cannot get scale factor uncertainty for type ")+typeid(cand).name());
  }
  return 0.;
}

double lepton_tools::getScaleFactorFs(const vCands &sig_leps, int npv){
  double scale_factor = 1.;
  for(const auto &lep: sig_leps){
    if(lep == nullptr) throw runtime_error("sig_leps contains a nullptr in lepton_tools::getScaleFactor");
    scale_factor *= getScaleFactorFs(*lep,npv);
  }
  return scale_factor;
}

double lepton_tools::getScaleFactorUncertaintyFs(const vCands &sig_leps, int npv){
  //Crashes if scale factor == 0...
  double uncertainty = 0., totsf(1.);
  for(const auto &lep: sig_leps){
    if(lep == nullptr) throw runtime_error("sig_leps contains a nullptr in lepton_tools::getScaleFactorUncertainty");
    double scale_factor = getScaleFactorFs(*lep, npv);
    if(scale_factor <= 0) throw runtime_error("One lepton FSSF is 0 in lepton_tools::getScaleFactorUncertainty");
    uncertainty = hypot(uncertainty, getScaleFactorUncertaintyFs(*lep,npv)/scale_factor);
    totsf *= scale_factor;
  }
  return uncertainty*totsf;
}


double lepton_tools::getScaleFactor(const reco::Muon &lep){
  auto id_bin = muon_id_sf.FindFixBin(lep.pt(), fabs(lep.eta()));
  auto iso_bin = muon_iso_sf.FindFixBin(lep.pt(), fabs(lep.eta()));
  auto id_overflow = muon_id_sf.IsBinOverflow(id_bin);
  auto iso_overflow = muon_iso_sf.IsBinOverflow(id_bin);
  auto id_val = id_overflow ? 1. : muon_id_sf.GetBinContent(id_bin);
  auto iso_val = iso_overflow ? 1. : muon_iso_sf.GetBinContent(iso_bin);
  return id_val * iso_val;
}

double lepton_tools::getScaleFactorUncertainty(const reco::Muon &lep){
  auto id_bin = muon_id_sf.FindFixBin(lep.pt(), fabs(lep.eta()));
  auto iso_bin = muon_iso_sf.FindFixBin(lep.pt(), fabs(lep.eta()));
  auto id_overflow = muon_id_sf.IsBinOverflow(id_bin);
  auto iso_overflow = muon_iso_sf.IsBinOverflow(id_bin);
  auto id_val = id_overflow ? 1. : muon_id_sf.GetBinContent(id_bin);
  auto iso_val = iso_overflow ? 1. : muon_iso_sf.GetBinContent(iso_bin);
  auto id_err = id_overflow ? 0. : muon_id_sf.GetBinError(id_bin);
  auto iso_err = iso_overflow ? 0. : muon_iso_sf.GetBinError(iso_bin);
  auto full_val = id_val*iso_val;
  // Adding relative uncertainty of ISO, ID, and two 1% systematics
  return full_val * hypot(hypot(hypot(iso_err/iso_val, id_err/id_val),0.01/full_val),0.01/full_val);
}

double lepton_tools::getScaleFactor(const pat::Electron &lep){
  auto id_bin = electron_id_sf.FindFixBin(lep.superCluster()->energy()*sin(lep.superClusterPosition().theta()), fabs(lep.superCluster()->eta()));
  auto iso_bin = electron_iso_sf.FindFixBin(lep.superCluster()->energy()*sin(lep.superClusterPosition().theta()), fabs(lep.superCluster()->eta()));
  auto id_val = electron_id_sf.GetBinContent(id_bin);
  auto iso_val = electron_iso_sf.GetBinContent(iso_bin);
  return id_val * iso_val;
}

double lepton_tools::getScaleFactorUncertainty(const pat::Electron &lep){
  auto id_bin = electron_id_sf.FindFixBin(lep.pt(), fabs(lep.superCluster()->eta()));
  auto iso_bin = electron_iso_sf.FindFixBin(lep.pt(), fabs(lep.superCluster()->eta()));
  auto id_val = electron_id_sf.GetBinContent(id_bin);
  auto iso_val = electron_iso_sf.GetBinContent(iso_bin);
  auto id_err = electron_id_sf.GetBinError(id_bin);
  auto iso_err = electron_iso_sf.GetBinError(iso_bin);
  // Adding relative uncertainty of ISO, ID
  auto full_val = id_val*iso_val;
  return full_val * hypot(iso_err/iso_val, id_err/id_val);
}

double lepton_tools::getScaleFactorFs(const reco::Muon &lep, int npv){
  auto bin = muon_idiso_fs_sf.FindFixBin(lep.pt(), fabs(lep.eta()), npv);
  auto overflow = muon_idiso_fs_sf.IsBinOverflow(bin);
  auto val = overflow ? 1. : muon_idiso_fs_sf.GetBinContent(bin); 
  return val;
}

double lepton_tools::getScaleFactorUncertaintyFs(const reco::Muon &lep, int npv){
  auto bin = muon_idiso_fs_sf.FindFixBin(lep.pt(), fabs(lep.eta()), npv);
  auto overflow = muon_idiso_fs_sf.IsBinOverflow(bin);
  auto val = overflow ? 1. : muon_idiso_fs_sf.GetBinContent(bin); 

  // Systematics : https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSFMC#Recommendations
  float syst=0;
  if(lep.pt()>10 && lep.pt()<=20)	syst = val*0.03;
  else if(lep.pt()>20 && lep.pt()<=30)	syst = val*0.01;
  else if(lep.pt()>30)			syst = val*0.01; 
  return syst;
}

double lepton_tools::getScaleFactorFs(const pat::Electron &lep, int npv){
  auto bin = electron_idiso_fs_sf.FindFixBin(lep.pt(), fabs(lep.eta()), npv);
  auto overflow = electron_idiso_fs_sf.IsBinOverflow(bin);
  auto val = overflow ? 1. : electron_idiso_fs_sf.GetBinContent(bin);
  return val;
}

double lepton_tools::getScaleFactorUncertaintyFs(const pat::Electron &lep, int npv){
  auto bin = electron_idiso_fs_sf.FindFixBin(lep.pt(), fabs(lep.eta()), npv);
  auto overflow = electron_idiso_fs_sf.IsBinOverflow(bin);
  auto val = overflow ? 1. : electron_idiso_fs_sf.GetBinContent(bin);

  // Systematics : https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSFMC#Recommendations
  float syst=0;
  if(lep.pt()>10 && lep.pt()<=20)	syst=val*0.10;
  else if(lep.pt()>20 && lep.pt()<=30)	syst=val*0.08;
  else if(lep.pt()>30)			syst=val*0.05; 
  return syst;
}

double lepton_tools::getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                                    const reco::Candidate* ptcl,
                                    double r_iso_min, double r_iso_max, double kt_scale,
                                    double rho, bool charged_only) {
  if (ptcl->pt()<1.) return 99999.;

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
    if (&pfc == ptcl) continue;
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
  double pu_corr = rho*effarea*pow(r_iso/0.3, 2);
  double iso(0.);
  if (charged_only) iso = iso_ch;
  else iso = iso_ch + max(0.,iso_ph + iso_nh - pu_corr);

  return iso/ptcl->pt();
}

vCands lepton_tools::getIsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi){

  vCands tks;
  //common parameters
  float mt_max = 100.;
  float eta_max = 2.5;
  float dz_max = 0.1;
  float cone_size = 0.3;

  for (size_t i(0); i < pfcands->size(); i++) {
    const pat::PackedCandidate &tk = (*pfcands)[i];
    unsigned int id = abs(tk.pdgId());

    //id-specific parameters
    float pt_min = id==211 ? 10. : 5.;
    float iso_max = id==211 ? 0.1 : 0.2;

    // track selection
    if (tk.charge()==0) continue;
    if (id!=11 && id!=13 && id!=211) continue;
    if (tk.pt() < pt_min) continue;
    if (fabs(tk.eta()) > eta_max) continue;
    if (fabs(tk.dz()) > dz_max) continue;
    if (mt_max>0.01 && getMT(met, met_phi,  tk.pt(), tk.phi())>mt_max) continue;

    // calculate track isolation
    double iso = 0.;
    for (size_t j(0); j < pfcands->size(); j++) {
      if (i==j) continue;
      const pat::PackedCandidate &pfc = (*pfcands)[j];

      if (pfc.charge()==0) continue;
      if (deltaR(tk,pfc) > cone_size) continue;
      if (fabs(pfc.dz()) > dz_max) continue;
      iso += pfc.pt();
    }
    if (iso/tk.pt()>iso_max) continue;

    tks.push_back(dynamic_cast<const reco::Candidate *>(&tk));
  }

  return tks;
}

vCands lepton_tools::getRA4IsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi, double rhoEventCentral, vector<float> &isos, int primary_pdg){

  vCands tks;
  //common parameters
  // float eta_max = 2.5;
  //float dz_max = 0.1;
  // float cone_size = 0.3;
  float eta_max = 5;
  float dz_max = 5;

  for (size_t i(0); i < pfcands->size(); i++) {
    const pat::PackedCandidate &tk = (*pfcands)[i];
    unsigned int id = abs(tk.pdgId());

    //id-specific parameters
    float pt_min = id==211 ? 5. : 5.;
    // float iso_max = id==211 ? 0.1 : 0.2;

    // track selection
    if (tk.charge()==0 || tk.charge()==(-primary_pdg)) continue;
    if (id!=11 && id!=13 && id!=211) continue;
    if (tk.pt() < pt_min) continue;
    if (fabs(tk.eta()) > eta_max) continue;
    if (fabs(tk.dz()) > dz_max) continue;

    //if (mt_max>0.01 && getMT(met, met_phi,  tk.pt(), tk.phi())>mt_max) continue;

    // calculate track isolation
    double iso = 0.;
    if(id!=211)
      iso = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&tk), 0.05, 0.2, 10., rhoEventCentral, false);
    else
      iso = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&tk), 0.05, 0.2, 10., rhoEventCentral, true);

    if(iso>1.0) continue;
    isos.push_back(iso);
    tks.push_back(dynamic_cast<const reco::Candidate *>(&tk));
  }

  return tks;
}

lepton_tools::lepton_tools(){
}

lepton_tools::~lepton_tools(){
}
