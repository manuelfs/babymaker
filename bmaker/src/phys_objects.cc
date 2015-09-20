//// PHYS_OBJECTS: Common physics objects definitions
//// Function names follow the first-lowercase, following words-uppercase. No underscores

// System include files
#include <algorithm>

// FW include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

// User include files
#include "babymaker/bmaker/interface/phys_objects.hh"

using namespace std;
using namespace phys_objects;


namespace phys_objects{

  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// CUTS /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  // float JetPtCut		= 30.0;
  // float JetEtaCut               = 2.4;

  float SignalLeptonPtCut	= 20.0;
  float VetoLeptonPtCut		= 10.0;
  float MuonEtaCut		= 2.4;
  float ElectronEtaCut		= 2.5;
  float MuonMiniIsoCut		= 0.2;
  float ElectronMiniIsoCut	= 0.1;


  ////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// LEPTONS ///////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  //////////////////// Muons
  bool isSignalMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
    // pT, eta cuts
    if(lep.pt() <= VetoLeptonPtCut) return false;
    if(fabs(lep.eta()) > MuonEtaCut) return false;
    // ID cuts (includes dz/dz0 cuts)
    if(!idMuon(lep, vtx, kMedium)) return false;
    // Isolation cuts
    if(lep_iso >= 0 && lep_iso > MuonMiniIsoCut) return false;

    return true;
  }

  bool isVetoMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
    // pT, eta cuts
    if(lep.pt() <= VetoLeptonPtCut) return false;
    if(fabs(lep.eta()) > MuonEtaCut) return false;
    // ID cuts (includes dz/dz0 cuts)
    if(!idMuon(lep, vtx, kLoose)) return false;
    // Isolation cuts
    if(lep_iso >= 0 && lep_iso > MuonMiniIsoCut) return false;

    return true;
  }

  bool idMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold) {
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
    
  bool vertexMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
    dz = 0.; d0 = 0.;
    if(lep.track().isAvailable()){ // If the track is not available we probably don't want the muon
      dz = lep.track()->vz()-vtx->at(0).z();
      d0 = lep.track()->d0()-vtx->at(0).x()*sin(lep.track()->phi())+vtx->at(0).y()*cos(lep.track()->phi());
    } 
    if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

    return true;
  }

  
  //////////////////// Electrons
  bool isSignalElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
    // pT, eta cuts
    if(lep.pt() <= VetoLeptonPtCut) return false;
    if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
    // ID cuts (includes dz/dz0 cuts)
    if(!idElectron(lep, vtx, kMedium)) return false;
    // Isolation cuts
    if(lep_iso >= 0 && lep_iso > ElectronMiniIsoCut) return false;

    return true;
  }

  bool isVetoElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lep_iso){
    // pT, eta cuts
    if(lep.pt() <= VetoLeptonPtCut) return false;
    if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
    // ID cuts (includes dz/dz0 cuts)
    if(!idElectron(lep, vtx, kVeto)) return false;
    // Isolation cuts
    if(lep_iso >= 0 && lep_iso > ElectronMiniIsoCut) return false;

    return true;
  }

  bool idElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold, bool do_iso) {

    bool barrel(lep.isEB());
    double deta_cut, dphi_cut, ieta_cut, hovere_cut, d0_cut, dz_cut,
      ooeminusoop_cut, reliso_cut, misshits_cut;
    bool req_conv_veto;

    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Working_points_for_Spring15_MC_s
    // Values from September 16th
    if(barrel){
      ieta_cut        = chooseVal(threshold	, 0.012		,  0.0105	,  0.0101	, 0.0101 );
      deta_cut        = chooseVal(threshold	, 0.0126	,  0.00976	,  0.0094	, 0.0095 );
      dphi_cut        = chooseVal(threshold	, 0.107		,  0.0929	,  0.0296	, 0.0291 );
      hovere_cut      = chooseVal(threshold	, 0.186		,  0.0765	,  0.0372	, 0.0372 );
      reliso_cut      = chooseVal(threshold	, 0.161		, 0.118		,  0.0987	, 0.0468 );
      ooeminusoop_cut = chooseVal(threshold	, 0.239		, 0.184		,  0.118	, 0.0174 );
      d0_cut          = chooseVal(threshold	, 0.0621	,  0.0227	,  0.0151	, 0.0144 );
      dz_cut          = chooseVal(threshold	, 0.613		, 0.379		,  0.238	, 0.323  );
      misshits_cut    = chooseVal(threshold	, 2		,  2		,  2		,  2	 );
      req_conv_veto   = chooseVal(threshold	, true		,  true		,  true		,  true );
    } else {
      ieta_cut        = chooseVal(threshold	, 0.0339	,  0.0318	,  0.0287	, 0.0287 );
      deta_cut        = chooseVal(threshold	, 0.0109	,  0.00952	,  0.00773	, 0.00762);
      dphi_cut        = chooseVal(threshold	, 0.219		,  0.181	,  0.148	, 0.0439 );
      hovere_cut      = chooseVal(threshold	, 0.0962	,  0.0824	,  0.0546	, 0.0544 );
      reliso_cut      = chooseVal(threshold	, 0.193		,  0.118	,  0.0902	, 0.0759 );
      ooeminusoop_cut = chooseVal(threshold	, 0.141		,  0.125	,  0.104	, 0.01   );
      d0_cut          = chooseVal(threshold	, 0.279		,  0.242	,  0.0535	, 0.0377 );
      dz_cut          = chooseVal(threshold	, 0.947		,  0.921	,  0.572	, 0.571	 );
      misshits_cut    = chooseVal(threshold	, 3		,  1		,  1		, 1	 );
      req_conv_veto   = chooseVal(threshold	, true		,  true		,  true		,  true	 );
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

  bool vertexElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
    dz = 0.; d0 = 0.;
    if(lep.gsfTrack().isAvailable()){ // If the track is not available we probably don't want the electron
      dz = lep.gsfTrack()->vz()-vtx->at(0).z();
      d0 = lep.gsfTrack()->d0()-vtx->at(0).x()*sin(lep.gsfTrack()->phi())+vtx->at(0).y()*cos(lep.gsfTrack()->phi());
    } 
    if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

    return true;
  }

  double getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			const reco::Candidate* ptcl,  
			double r_iso_min, double r_iso_max, double kt_scale,
			bool charged_only) {

    if (ptcl->pt()<5.) return 99999.;

    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    double iso_nh(0.); double iso_ch(0.); 
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5), r_mini(kt_scale/ptcl->pt());
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = max(r_iso_min, min(r_iso_max, r_mini));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      double dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
	if (pfc.pt()>ptThresh) {
	  /////////// PHOTONS ////////////
	  if (abs(pfc.pdgId())==22) {
	    if(dr < deadcone_ph) continue;
	    iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
	  } else if (abs(pfc.pdgId())==130) {
	    if(dr < deadcone_nh) continue;
	    iso_nh += pfc.pt();
	  }
	}
	//////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
	if (abs(pfc.pdgId())==211) {
	  if(dr < deadcone_ch) continue;
	  iso_ch += pfc.pt();
	}
	//////////////////  CHARGED from PU  /////////////////////////
      } else {
	if (pfc.pt()>ptThresh){
	  if(dr < deadcone_pu) continue;
	  iso_pu += pfc.pt();
	}
      }
    }
    double iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;
  }


}
