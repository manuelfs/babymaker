//// BMAKER_FULL: Creates baby tree with full branches
//// Function names follow the first-lowercase, following words-uppercase. No underscores

// System include files
#include <cmath>
#include <memory>
#include <iostream>
#include <stdlib.h>
#include <iomanip>  // put_time

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// FW physics include files
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#ifdef POST_7_4
  #include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#endif
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"

// ROOT include files
#include "TFile.h"
#include "TROOT.h"

// User include files
#include "babymaker/bmaker/interface/bmaker_full.hh"
#include "babymaker/bmaker/interface/baby_full.hh"
#include "babymaker/bmaker/interface/cross_sections.hh"

using namespace std;
using namespace utilities;

///////////////////////// analyze: METHOD CALLED EACH EVENT ///////////////////////////
void bmaker_full::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  nevents++;
  isData = iEvent.isRealData();
  baby.Clear();

  ////////////////////// Event info /////////////////////
  baby.run() = iEvent.id().run();
  baby.event() = iEvent.id().event();
  baby.lumiblock() = iEvent.luminosityBlock();
  // We are applying the golden JSON with lumisToProcess in bmaker_full_cfg.py
  if(isData){
    if (debug) cout<<"INFO: Checking JSON..."<<endl;
    baby.nonblind() = eventTool->isInJSON("nonblind", baby.run(), baby.lumiblock());
    baby.json2p6() = eventTool->isInJSON("json2p6", baby.run(), baby.lumiblock());
  } else {
    baby.nonblind() = true;
    baby.json2p6() = true;
  }

  ////////////////////// Trigger /////////////////////
  if (debug) cout<<"INFO: Processing trigger info..."<<endl;
  bool triggerFired;
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(tok_trigResults_hlt_,triggerBits);  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(tok_patTrig_,triggerPrescales);  

  if(triggerBits.isValid() && triggerPrescales.isValid()){
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    // Not saving data events that don't have triggers we care about
    triggerFired = writeTriggers(names, triggerBits, triggerPrescales);
  }
  else
    triggerFired = true;
  if(!triggerFired && isData) {
    reportTime(iEvent);
    return;
  }

  ////////////////////// Primary vertices /////////////////////
  if (debug) cout<<"INFO: Writing vertices..."<<endl;
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByToken(tok_primVertex_, vtx);
  edm::Handle<std::vector< PileupSummaryInfo > >  pu_info;
  if(!isData) {
    iEvent.getByToken(tok_addPileup_, pu_info);
    if(!pu_info.isValid()) iEvent.getByToken(tok_slimAddPileup_, pu_info);
  }
  writeVertices(vtx, pu_info);

  ////////////////////// Leptons /////////////////////
  if (debug) cout<<"INFO: Writing leptons..."<<endl;
  // pfcands, to be used in calculating isolation
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(tok_packPFCands_, pfcands);

  // Event energy density, to be used in calculating isolation and JECs
  edm::Handle<double> rhoEventCentral_h;
  iEvent.getByToken(tok_rhoFastJet_centralNeutral_, rhoEventCentral_h);
  double rhoEventCentral = static_cast<double>(*rhoEventCentral_h);

  vCands sig_leps, veto_leps, sig_mus, veto_mus, sig_els, veto_els;
  vCands all_mus, all_els;
  edm::Handle<pat::MuonCollection> allmuons;
  iEvent.getByToken(tok_muons_, allmuons);
  sig_mus = writeMuons(allmuons, pfcands, vtx, veto_mus, all_mus, rhoEventCentral);
  edm::Handle<pat::ElectronCollection> allelectrons;
  iEvent.getByToken(tok_electrons_, allelectrons);
  sig_els = writeElectrons(allelectrons, pfcands, vtx, veto_els, all_els, rhoEventCentral);
  
  writeDiLep(sig_mus, sig_els, veto_mus, veto_els);

  // Putting muons and electrons together
  sig_leps = sig_mus;
  sig_leps.insert(sig_leps.end(), sig_els.begin(), sig_els.end());
  writeLeptons(sig_leps);
  // if (baby.nleps()<1) return;
  veto_leps = veto_mus;
  veto_leps.insert(veto_leps.end(), veto_els.begin(), veto_els.end());

  ///////////////////////////////// Photons ////////////////////////////////
  if (debug) cout<<"INFO: Writing photons..."<<endl;
  vCands photons;
  edm::Handle<double> rhoEvent_h;
  iEvent.getByToken(tok_rhoFastJet_all_, rhoEvent_h);
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByToken(tok_offBeamSpot_, beamspot);
  edm::Handle<pat::PhotonCollection> allphotons;
  iEvent.getByToken(tok_photons_, allphotons);
  edm::Handle<vector<reco::Conversion> > conversions;
  iEvent.getByToken(tok_reducedEgamma_conver_, conversions);
  photons = writePhotons(allphotons, allelectrons, conversions, beamspot, *rhoEvent_h);

  //////////////////////////// MET/JETs with JECs ///////////////////////////
  if (debug) cout<<"INFO: Applying JECs..."<<endl;
  edm::Handle<pat::JetCollection> alljets;
  iEvent.getByToken(tok_jets_, alljets);
  edm::Handle<edm::View<reco::GenJet> > genjets;
  iEvent.getByToken(tok_genJets_, genjets) ;
  jetTool->getJetCorrections(genjets, alljets, *rhoEvent_h);

  /// MET
  if (debug) cout<<"INFO: Writing MET..."<<endl;
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(tok_met_, mets);
  edm::Handle<pat::METCollection> mets_nohf;
  iEvent.getByToken(tok_met_noHF_, mets_nohf);
  writeMET(mets, mets_nohf);

  /// isolated tracks need to be after MET calculation and before jets cleaning
  if (debug) cout<<"INFO: Calculating track veto..."<<endl;
  vCands tks,ra4tks;
  if (eventTool->hasGoodPV(vtx)){
    // RA2/b track veto
    tks = lepTool->getIsoTracks(pfcands, baby.met(), baby.met_phi());
    baby.ntks() = tks.size();

    // RA4 track veto
    writeTks(pfcands,vtx,rhoEventCentral);
  }

  /// Jets
  if (debug) cout<<"INFO: Writing jets..."<<endl;
  vector<LVector> jets, clean_jets;
  vector<double> jetsMuonEnergyFrac;
  vector<vector<LVector> > sys_jets;
  if (!doSystematics) {
    clean_jets = writeJets(alljets, genjets, sig_leps, veto_leps, photons, tks, jets, sys_jets, jetsMuonEnergyFrac);
  } else {
    clean_jets = writeJets(alljets, genjets, sig_leps, veto_leps, photons, tks, jets, sys_jets, jetsMuonEnergyFrac);
    for (unsigned isys(0); isys<kSysLast; isys++){
      baby.sys_mj().push_back(jetTool->getSysMJ(1.2, sys_jets[isys]));
      baby.sys_mj08().push_back(jetTool->getSysMJ(0.8, sys_jets[isys]));
      baby.sys_mj13().push_back(jetTool->getSysMJ(1.3, sys_jets[isys]));
      baby.sys_mj14().push_back(jetTool->getSysMJ(1.4, sys_jets[isys]));
      baby.sys_mj15().push_back(jetTool->getSysMJ(1.5, sys_jets[isys]));
      baby.sys_mj16().push_back(jetTool->getSysMJ(1.6, sys_jets[isys]));
    }
  }
  writeFatJets(jets);

  ////////////////////// mT, dphi /////////////////////
  // It requires previous storing of MET
  if (debug) cout<<"INFO: Calculating mT..."<<endl;
  if(sig_leps.size()>0){
    float wx = baby.met()*cos(baby.met_phi()) + sig_leps[0]->px();
    float wy = baby.met()*sin(baby.met_phi()) + sig_leps[0]->py();
    float wphi = atan2(wy, wx);
    baby.dphi_wlep() = deltaPhi(wphi, sig_leps[0]->phi());

    baby.mt() = getMT(baby.met(), baby.met_phi(),  sig_leps[0]->pt(), sig_leps[0]->phi());
    baby.mt_nohf() = getMT(baby.met_nohf(), baby.met_nohf_phi(),  sig_leps[0]->pt(), sig_leps[0]->phi());
  }   
  // calculating met with systematics and corresponding mT
  if (doSystematics) {
    baby.sys_met().resize(kSysLast,-999.);
    baby.sys_mt().resize(kSysLast,-999.);
    for (unsigned isys(0); isys<kSysLast; isys++) {
      float sys_met_phi(0.);
      jetTool->getMETWithJEC(mets, baby.sys_met()[isys], sys_met_phi, isys);
      if (sig_leps.size()>0) 
        baby.sys_mt()[isys] = getMT(baby.sys_met()[isys], sys_met_phi, sig_leps[0]->pt(), sig_leps[0]->phi());                              
    }
  }

  ///////////////////// Filters ///////////////////////
  // the HBHE noise filter needs to be recomputed in early 2015 data
  if (debug) cout<<"INFO: Writing filters..."<<endl; 
/* now these are taken from miniAOD 
  edm::Handle<bool> filter_hbhe;
  if(iEvent.getByToken(tok_HBHENoiseFilter_,filter_hbhe)) { 
    if(filter_hbhe.isValid()) 
      baby.pass_hbhe() = *filter_hbhe;
    else
      baby.pass_hbhe() = true;

    iEvent.getByToken(tok_HBHEIsoNoiseFilter_,filter_hbhe);
    if(filter_hbhe.isValid()) 
      baby.pass_hbheiso() = *filter_hbhe;
    else
      baby.pass_hbheiso() = true;
  }
*/
  edm::Handle<edm::TriggerResults> filterBits;
  if(isData) iEvent.getByToken(tok_trigResults_reco_,filterBits);  
  else iEvent.getByToken(tok_trigResults_pat_,filterBits);  
  // re-recoed data will have the process label "PAT" rather than "RECO";
  // if the attempt to find data with "RECO" process fails, try "PAT"
  if(!filterBits.isValid() && isData) 
    iEvent.getByToken(tok_trigResults_pat_,filterBits);  
  if(filterBits.isValid()){
    const edm::TriggerNames &fnames = iEvent.triggerNames(*filterBits);
    //this method uses baby.pass_jets(), so call only after writeJets()!!
    writeFilters(fnames, filterBits, vtx, jetsMuonEnergyFrac);
  }

  //////////////// HLT objects //////////////////
  if (debug) cout<<"INFO: Writing HLT objects..."<<endl;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(tok_selectedPatTrig_,triggerObjects);  
  // Requires having called writeMuons and writeElectrons for truth-matching
  if(triggerBits.isValid() && triggerPrescales.isValid()){
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    writeHLTObjects(names, triggerObjects, all_mus, all_els, photons);
  }

  ////////////////// MC particles and Truth-matching //////////////////
  if (!isData) {
    if (debug) cout<<"INFO: Writing MC particles..."<<endl;
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(tok_pruneGenPart_, genParticles);
    writeMC(genParticles, all_mus, all_els, photons);
    nisrMatch(genParticles, clean_jets);
  }

  ////////////////// resolution-corrected MET /////////////////////////
  if (debug) cout<<"INFO: Calculating MET rebalancing..."<<endl;
  baby.jetmismeas() = false;
  if(doMetRebalancing && sig_leps.size()==1) {
    double temp_met = baby.met();
    double temp_met_phi = baby.met_phi();
    rebalancedMET(temp_met, temp_met_phi);
    baby.met_rebal() = temp_met;
    baby.mt_rebal() = getMT(temp_met, temp_met_phi, sig_leps[0]->pt(), sig_leps[0]->phi());
    if(baby.met_rebal()/baby.met()<0.2 && baby.mt_rebal()<150) baby.jetmismeas()=true;
  } else {
    // use default values for events that do not have exactly one lepton
    baby.mt_rebal() = baby.mt();
    baby.met_rebal() = baby.met();
  }

  ///////////////////// MC hard scatter info ///////////////////////
  if (debug) cout<<"INFO: Retrieving hard scatter info..."<<endl;
  edm::Handle<LHEEventProduct> lhe_info;
  baby.stitch() = true;
  if (outname.Contains("SMS-") && outname.Contains("PUSpring16Fast")) {
    baby.mgluino() = mprod_;
    baby.mlsp() = mlsp_;
  } else if (!isData) {
    iEvent.getByToken(tok_extLHEProducer_, lhe_info);
    if(!lhe_info.isValid()) iEvent.getByToken(tok_source_, lhe_info);
    if(lhe_info.isValid()) writeGenInfo(lhe_info);
    if((outname.Contains("TTJets") && (outname.Contains("Lept") || outname.Contains("TTJets_Tune")) && baby.ht_isr_me()>600)
       || (outname.Contains("DYJetsToLL_M-50_TuneCUETP8M")  && baby.ht_isr_me()>100)
       || (outname.Contains("WJetsToLNu_TuneCUETP8M1")  && baby.ht_isr_me()>100))
      baby.stitch() = false;
    if(outname.Contains("TTJets_Tune") && baby.ntruleps()!=0) baby.stitch()=false; 
  } // if it is not data
         
  //////////////// Weights and systematics //////////////////
  // w_btag calculated in writeJets
  // w_toppt and sys_isr calculated in writeMC
  edm::Handle<GenEventInfoProduct> gen_event_info;
  iEvent.getByToken(tok_generator_, gen_event_info);
  writeWeights(sig_leps, gen_event_info, lhe_info);

  ////////////////// Filling the tree //////////////////
  baby.Fill();

  reportTime(iEvent);

}


/*
  ______                      _                     _ _   _             
  | ___ \                    | |                   (_) | (_)            
  | |_/ /_ __ __ _ _ __   ___| |__   __      ___ __ _| |_ _ _ __   __ _ 
  | ___ \ '__/ _` | '_ \ / __| '_ \  \ \ /\ / / '__| | __| | '_ \ / _` |
  | |_/ / | | (_| | | | | (__| | | |  \ V  V /| |  | | |_| | | | | (_| |
  \____/|_|  \__,_|_| |_|\___|_| |_|   \_/\_/ |_|  |_|\__|_|_| |_|\__, |
  __/ |
  |___/ 
*/

// Requires having called jetTool->getJetCorrections(alljets, rhoEvent_) beforehand
void bmaker_full::writeMET(edm::Handle<pat::METCollection> mets, edm::Handle<pat::METCollection> mets_nohf){
  jetTool->getMETWithJEC(mets, baby.met(), baby.met_phi(), kSysLast);
  jetTool->getMETRaw(mets, baby.met_raw(), baby.met_raw_phi());
  baby.met_mini() = mets->at(0).pt();
  baby.met_mini_phi() = mets->at(0).phi();
  baby.met_calo() = mets->at(0).caloMETPt();
  baby.met_calo_phi() = mets->at(0).caloMETPhi();
  if(mets_nohf.isValid()){
    baby.met_nohf() = mets_nohf->at(0).pt();
    baby.met_nohf_phi() = mets_nohf->at(0).phi();
  }
  if(!isData){
    baby.met_tru() = mets->at(0).genMET()->pt();
    baby.met_tru_phi() = mets->at(0).genMET()->phi();
  }
  if(isnan(baby.met())) {
    cout<<"MET is NaN. Perhaps you are running on old miniAOD with a new release. Setting MET to zero."<<endl;
    baby.met() = 0;
  }
}

// Requires having called jetTool->getJetCorrections(alljets, rhoEvent_) beforehand
vector<LVector> bmaker_full::writeJets(edm::Handle<pat::JetCollection> alljets, 
                            edm::Handle<edm::View <reco::GenJet> > genjets,
                            vCands &sig_leps, vCands &veto_leps, vCands &photons, vCands &tks,
                            vector<LVector> &jets, 
                            vector<vector<LVector> > &sys_jets,
                            vector<double> &jetsMuonEnergyFrac){
  vector<int> hi_csv(5,-1); // Indices of the 5 jets with highest CSV
  vector<LVector> clean_jets;
  vCands jets_ra2;  vCands jets20_cands; vector<LVector> jets20;
  LVector jetsys_p4, jetsys_nob_p4;
  baby.njets() = 0; baby.nbl() = 0; baby.nbm() = 0;  baby.nbt() = 0;
  baby.njets20() = 0; baby.nbl20() = 0; baby.nbm20() = 0;  baby.nbt20() = 0;  
  baby.ht() = 0.; baby.ht_hlt() = 0.;
  baby.njets_ra2() = 0; baby.njets_clean() = 0; baby.nbm_ra2() = 0; baby.ht_ra2() = 0.; baby.ht_clean() = 0.; 
  baby.pass_jets() = true; baby.pass_jets_nohf() = true; baby.pass_jets_tight() = true; 
  baby.pass_jets_ra2() = true; baby.pass_jets_tight_ra2() = true; 
  baby.sys_bctag().resize(2, 1.); baby.sys_udsgtag().resize(2, 1.);
  baby.sys_bctag_loose().resize(2, 1.); baby.sys_udsgtag_loose().resize(2, 1.);
  baby.w_btag() = baby.w_btag_loose() = 1.;
  baby.hig_dphi() = 9999.;
  if (isFastSim){ baby.sys_fs_bctag().resize(2, 1.); baby.sys_fs_udsgtag().resize(2, 1.);}
  if (doSystematics) {
    baby.sys_njets().resize(kSysLast, 0); baby.sys_nbm().resize(kSysLast, 0); 
    baby.sys_pass().resize(kSysLast, true); baby.sys_ht().resize(kSysLast, 0.); 
    sys_jets.resize(kSysLast, vector<LVector>());
  }
  float mht_px(0.), mht_py(0.), mht_clean_px(0.), mht_clean_py(0.);
  for (size_t ijet(0); ijet < alljets->size(); ijet++) {
    const pat::Jet &jet = (*alljets)[ijet];
    if(fabs(jet.eta()) > 5) continue;

    LVector jetp4(jetTool->corrJet[ijet]);
    float csv(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    bool isLep = jetTool->leptonInJet(jet, sig_leps);
    bool looseID = jetTool->idJet(jet, jetTool->kLoose);
    bool tightID = jetTool->idJet(jet, jetTool->kTight);
    bool goodPtEtaLoose = jetp4.pt() > jetTool->JetPtCutLoose && fabs(jet.eta()) <= jetTool->JetEtaCut;
    bool goodPtEta = jetp4.pt() > jetTool->JetPtCut && fabs(jet.eta()) <= jetTool->JetEtaCut;
    
    //   RA4 Jet Quality filters
    //--------------------------------
    if(jetp4.pt() > jetTool->JetPtCut && !isLep) {
      if(goodPtEta && !looseID) baby.pass_jets_nohf() = false;
      if(!looseID) baby.pass_jets() = false;
      if(!tightID) baby.pass_jets_tight() = false;
    }

    //    RA4 Jet Variables
    //----------------------------
    if((looseID && goodPtEtaLoose) || isLep) {
      baby.jets_pt().push_back(jetp4.pt());
      baby.jets_eta().push_back(jet.eta());
      baby.jets_phi().push_back(jet.phi());
      baby.jets_m().push_back(jetp4.mass());
      baby.jets_islep().push_back(isLep);
      jetsMuonEnergyFrac.push_back(jet.muonEnergyFraction());
      if(!isData && jetTool->genJetPt[ijet]>0.) baby.jets_pt_res().push_back(jetp4.pt()/jetTool->genJetPt[ijet]);
      else baby.jets_pt_res().push_back(-99999.);
      baby.jets_hflavor().push_back(jet.hadronFlavour());
      baby.jets_csv().push_back(csv);
      if (goodPtEta && !isLep) clean_jets.push_back(jetp4);
      if (goodPtEta || isLep) jets.push_back(jetp4);
      jets20_cands.push_back(dynamic_cast<const reco::Candidate *>(&jet));
      jets20.push_back(jetp4);
      baby.jets_h1().push_back(false);
      baby.jets_h2().push_back(false);
      float dphi = reco::deltaPhi(jet.phi(), baby.met_phi());
      if(dphi < baby.hig_dphi()) baby.hig_dphi() = dphi;

      // Finding the N jets with highest CSV values
      for(size_t ind(0); ind<hi_csv.size(); ind++){
        int icsv = hi_csv[ind];
        if(icsv==-1 || csv > baby.jets_csv()[icsv]){
          for(size_t ind2(hi_csv.size()-1); ind2>=ind+1; ind2--) hi_csv[ind2] = hi_csv[ind2-1];
          hi_csv[ind] = baby.jets_csv().size()-1;
          break;
        }
      } // Loop over highest CSV jets


      if(!isLep){
        if(addBTagWeights && goodPtEta) {
          bool btag(csv > jetTool->CSVMedium);
          bool btagLoose(csv > jetTool->CSVLoose);
          //central weight for fastsim taken into account together with the fullsim inside jetTool->jetBTagWeight()
          baby.w_btag() *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM,
						  "central", "central");
          // now, vary only the full sim scale factor, regardless of whether we run on Fast or Full sim
          // this is necessary since uncertainties for FastSim and FullSim are uncorrelated 
          baby.sys_bctag()[0]   *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM, "up",      "central");
          baby.sys_bctag()[1]   *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM, "down",    "central");
          baby.sys_udsgtag()[0] *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM, "central", "up");
          baby.sys_udsgtag()[1] *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM, "central", "down");
          baby.w_btag_loose()         *= jetTool->jetBTagWeight(jet, jetp4, btagLoose, BTagEntry::OP_LOOSE, "central", "central");
          baby.sys_bctag_loose()[0]   *= jetTool->jetBTagWeight(jet, jetp4, btagLoose, BTagEntry::OP_LOOSE, "up",      "central");
          baby.sys_bctag_loose()[1]   *= jetTool->jetBTagWeight(jet, jetp4, btagLoose, BTagEntry::OP_LOOSE, "down",    "central");
          baby.sys_udsgtag_loose()[0] *= jetTool->jetBTagWeight(jet, jetp4, btagLoose, BTagEntry::OP_LOOSE, "central", "up");
          baby.sys_udsgtag_loose()[1] *= jetTool->jetBTagWeight(jet, jetp4, btagLoose, BTagEntry::OP_LOOSE, "central", "down");
          if (isFastSim) { 
            // now we vary only the FastSim SF
            baby.sys_fs_bctag()[0]   *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM,
							       "central", "central", "up",      "central");
            baby.sys_fs_bctag()[1]   *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM,
							       "central", "central", "down",    "central");
            baby.sys_fs_udsgtag()[0] *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM,
							       "central", "central", "central", "up");
            baby.sys_fs_udsgtag()[1] *= jetTool->jetBTagWeight(jet, jetp4, btag, BTagEntry::OP_MEDIUM,
							       "central", "central", "central", "down");
	  }
	}
	if (goodPtEta){
	  jetsys_p4 += jet.p4();
	  baby.njets()++;
	  baby.ht() += jetp4.pt();
	  if(csv > jetTool->CSVLoose)  baby.nbl()++;
	  if(csv > jetTool->CSVMedium) baby.nbm()++;
	  else jetsys_nob_p4 += jet.p4();
	  if(csv > jetTool->CSVTight)  baby.nbt()++;
	} 
	baby.njets20()++;
	if(csv > jetTool->CSVLoose)  baby.nbl20()++;
	if(csv > jetTool->CSVMedium) baby.nbm20()++;
	if(csv > jetTool->CSVTight)  baby.nbt20()++;
	
      }
    }

    //    HLT HT definition
    //----------------------------
    if(jetp4.pt() > jetTool->JetHLTPtCut && fabs(jet.eta()) <= jetTool->JetHLTEtaCut) baby.ht_hlt() += jetp4.pt();

    //    RA2 Jet Variables
    //----------------------------
    bool isLep_ra2 = jetTool->jetMatched(jet, veto_leps); // Uses RA2/b's loose matching, dpt/pt < 100%, dR < 0.4
    bool isPhoton = jetTool->jetMatched(jet, photons);    // Uses RA2/b's loose matching, dpt/pt < 100%, dR < 0.4
    bool isIsoTrack = jetTool->jetMatched(jet, tks);      // Uses RA2/b's loose matching, dpt/pt < 100%, dR < 0.4
    bool applyId_ra2 = !isLep_ra2 && !isPhoton && !isIsoTrack; // Only check ID if jet not matched
    bool goodJet_ra2 = (looseID || !applyId_ra2);
    bool tightJet_ra2 = (tightID || !applyId_ra2);

    if(jetp4.pt() > jetTool->JetPtCut) { // Check jet ID on 30 GeV jets
      if(!goodJet_ra2) baby.pass_jets_ra2() = false;
      if(!tightJet_ra2) baby.pass_jets_tight_ra2() = false;
    }

    if(goodPtEta && goodJet_ra2) {
      baby.njets_ra2()++;
      baby.ht_ra2() += jetp4.pt();
      jets_ra2.push_back(dynamic_cast<const reco::Candidate *>(&jet)); 
      if(csv > jetTool->CSVMedium) baby.nbm_ra2()++;
      if(!isLep_ra2 && !isPhoton) {
        baby.njets_clean()++;
        baby.ht_clean() += jetp4.pt();
      }
    }
    if(goodJet_ra2 && jetp4.pt() > jetTool->JetPtCut && fabs(jet.eta()) <= jetTool->JetMHTEtaCut){
      mht_px -= jet.px();
      mht_py -= jet.py();
      if(!isLep_ra2 && !isPhoton) {
        mht_clean_px -= jet.px();
        mht_clean_py -= jet.py();
      }
    }

    //   Systematic variations
    //----------------------------
    if (doSystematics){
      for (unsigned isys(0); isys<kSysLast; isys++){
        jetp4 = jetTool->corrJet[ijet];
        if (isys == kSysJER) jetp4 *= (1 + jetTool->jerUnc[ijet]);
        else if (isys == kSysJECUp) jetp4 *= (1 + jetTool->jecUnc[ijet]);
        else if (isys == kSysJECDn) jetp4 *= (1 - jetTool->jecUnc[ijet]);
        //for now store pass_jets in the baby.sys_pass
        //variable will be updated with rest of filters in writeFilters()
        if(jetp4.pt() > jetTool->JetPtCut && !isLep && !looseID) baby.sys_pass()[isys] = false;

        if(looseID && jetp4.pt() > jetTool->JetPtCut && fabs(jet.eta()) <= jetTool->JetEtaCut) {
          sys_jets[isys].push_back(jetp4);
          if(!isLep){
            baby.sys_njets()[isys]++;
            baby.sys_ht()[isys] += jetp4.pt();
            if(csv > jetTool->CSVMedium) baby.sys_nbm()[isys]++;
          }
        }
      } // loop over systematics
    } // do systematics
  } // Loop over jets  

  if(!isData) baby.ht_tru() = jetTool->trueHT(genjets);

  baby.mht() = hypot(mht_px, mht_py);
  baby.mht_phi() = atan2(mht_py, mht_px);
  baby.mht_clean() = hypot(mht_clean_px, mht_clean_py);
  baby.mht_clean_phi() = atan2(mht_clean_py, mht_clean_px);
  baby.low_dphi() = jetTool->isLowDphi(jets_ra2, baby.mht_phi(), baby.dphi1(), baby.dphi2(), baby.dphi3(), baby.dphi4());
  // saving the dphi values for 20 GeV jets instead of for RA2b
  baby.low_dphi20() = jetTool->isLowDphi(jets20_cands, baby.met_phi(), baby.dphi1(), baby.dphi2(), baby.dphi3(), baby.dphi4());

  for (unsigned i(0); i<2; i++){
    baby.sys_bctag()[i] /= baby.w_btag();
    baby.sys_udsgtag()[i] /= baby.w_btag();
    baby.sys_bctag_loose()[i] /= baby.w_btag_loose();
    baby.sys_udsgtag_loose()[i] /= baby.w_btag_loose();
    if (isFastSim) { 
      baby.sys_fs_bctag()[i] /= baby.w_btag();
      baby.sys_fs_udsgtag()[i] /= baby.w_btag();
    }
  }
  // ISR system for Z->ll and tt->llbb configurations
  baby.jetsys_pt()  = jetsys_p4.pt();
  baby.jetsys_eta() = jetsys_p4.eta();
  baby.jetsys_phi() = jetsys_p4.phi();
  baby.jetsys_m()   = jetsys_p4.mass();
  baby.jetsys_nob_pt()  = jetsys_nob_p4.pt();
  baby.jetsys_nob_eta() = jetsys_nob_p4.eta();
  baby.jetsys_nob_phi() = jetsys_nob_p4.phi();
  baby.jetsys_nob_m()   = jetsys_nob_p4.mass();

  // write deltaR between csvm jets
  jetTool->getDeltaRbb(baby.dr_bb(), jets, baby.jets_csv(), baby.jets_islep());

  vector<size_t> branks = jet_met_tools::getBRanking(jets, baby.jets_csv(), baby.jets_islep());
  baby.dr_bb_2() = jet_met_tools::getDeltaRbb(jets, branks, baby.nbm());
  baby.dr_bb_max() = jet_met_tools::getDeltaRbbMax(jets, branks, baby.nbm());
  baby.dr_bb_min() = jet_met_tools::getDeltaRbbMin(jets, branks, baby.nbm());
  baby.dphi_bb_2() = jet_met_tools::getDeltaPhibb(jets, branks, baby.nbm());
  baby.dphi_bb_max() = jet_met_tools::getDeltaPhibbMax(jets, branks, baby.nbm());
  baby.dphi_bb_min() = jet_met_tools::getDeltaPhibbMin(jets, branks, baby.nbm());
  baby.m_bb_2() = jet_met_tools::getMbb(jets, branks, baby.nbm());
  baby.m_bb_max() = jet_met_tools::getMbbMax(jets, branks, baby.nbm());
  baby.m_bb_min() = jet_met_tools::getMbbMin(jets, branks, baby.nbm());

  if(sig_leps.size() > 0){
    baby.dr_bl_min() = jet_met_tools::getDeltaRblepMin(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dr_bl_max() = jet_met_tools::getDeltaRblepMax(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dr_bl_min2() = jet_met_tools::getDeltaRblepMin2(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dr_bl_max2() = jet_met_tools::getDeltaRblepMax2(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dphi_bl_min() = jet_met_tools::getDeltaPhiblepMin(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dphi_bl_max() = jet_met_tools::getDeltaPhiblepMax(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dphi_bl_min2() = jet_met_tools::getDeltaPhiblepMin2(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dphi_bl_max2() = jet_met_tools::getDeltaPhiblepMax2(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.m_bl_min() = jet_met_tools::getMblepMin(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.m_bl_max() = jet_met_tools::getMblepMax(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.m_bl_min2() = jet_met_tools::getMblepMin2(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.m_bl_max2() = jet_met_tools::getMblepMax2(jets, branks, baby.nbm(), sig_leps.at(0)->p4());
    baby.dphi_lmet() = fabs(reco::deltaPhi(sig_leps.at(0)->phi(), baby.met_phi()));
  }else{
    baby.dr_bl_min() = -1.;
    baby.dr_bl_max() = -1.;
    baby.dr_bl_min2() = -1.;
    baby.dr_bl_max2() = -1.;
    baby.dphi_bl_min() = -1.;
    baby.dphi_bl_max() = -1.;
    baby.dphi_bl_min2() = -1.;
    baby.dphi_bl_max2() = -1.;
    baby.m_bl_min() = -1.;
    baby.m_bl_max() = -1.;
    baby.m_bl_min2() = -1.;
    baby.m_bl_max2() = -1.;
    baby.dphi_lmet() = -1.;
  }

  baby.dphi_bmet_min() = jet_met_tools::getDeltaPhibmetMin(jets, branks, baby.nbm(), baby.met_phi());
  baby.dphi_bmet_max() = jet_met_tools::getDeltaPhibmetMax(jets, branks, baby.nbm(), baby.met_phi());
  baby.dphi_bmet_min2() = jet_met_tools::getDeltaPhibmetMin2(jets, branks, baby.nbm(), baby.met_phi());
  baby.dphi_bmet_max2() = jet_met_tools::getDeltaPhibmetMax2(jets, branks, baby.nbm(), baby.met_phi());
  baby.mt_bmet_min() = jet_met_tools::getMTbmetMin(jets, branks, baby.nbm(), baby.met(), baby.met_phi());
  baby.mt_bmet_max() = jet_met_tools::getMTbmetMax(jets, branks, baby.nbm(), baby.met(), baby.met_phi());
  baby.mt_bmet_min2() = jet_met_tools::getMTbmetMin2(jets, branks, baby.nbm(), baby.met(), baby.met_phi());
  baby.mt_bmet_max2() = jet_met_tools::getMTbmetMax2(jets, branks, baby.nbm(), baby.met(), baby.met_phi());


  if(baby.nbm() >= 2 && sig_leps.size()>0){
    const auto &jet1 = jets.at(branks.at(0));
    const auto &jet2 = jets.at(branks.at(1));
    double px = baby.met()*cos(baby.met_phi())+sig_leps.at(0)->px();
    double py = baby.met()*sin(baby.met_phi())+sig_leps.at(0)->py();
    baby.mt2() = getMT2(jet1.mass(), jet1.pt(), jet1.phi(),
                        jet2.mass(), jet2.pt(), jet2.phi(),
                        hypot(px, py), atan2(py, px));
    baby.mt2_0mass() = getMT2(0., jet1.pt(), jet1.phi(),
			      0., jet2.pt(), jet2.phi(),
			      hypot(px, py), atan2(py, px));
  }else{
    baby.mt2() = -1.;
  }

  //// Variables for the Higgsino analysis
  if(hi_csv[3]>=0){
    // hig_p4 has the p4 of the jet if row==col, and if not the sum of the p4 for the row-th and col-th jets 
    vector<vector<LVector> > hig_p4;
    for(size_t row=0; row<hi_csv.size(); row++){
      hig_p4.push_back(vector<LVector>());
      for(size_t col=0; col<=row; col++){
        LVector jetp4(jets20[hi_csv[row]]);
        hig_p4.back().push_back(jetp4);
        if(row!=col) hig_p4.back().back() += hig_p4[col][col];
      } // Loop over columns in hig_p4
    } // Loop over rows in hig_p4
    // Loop over all possible Higgs combination in the first nCSVs jets of highest CSV
    size_t nCSVs = hi_csv.size();
    nCSVs = 4;
    vector<int> hig_ind(4,0);
    float minDm(9999.);
    for(size_t ind0=0; ind0<nCSVs; ind0++){
      for(size_t ind1=0; ind1<ind0; ind1++){
        for(size_t ind2=ind0+1; ind2<nCSVs; ind2++){
          if(ind2==ind1) continue;
          for(size_t ind3=0; ind3<ind2; ind3++){
            if(ind3==ind0 || ind3==ind1) continue;
            float thisDm = fabs(hig_p4[ind0][ind1].mass() - hig_p4[ind2][ind3].mass());
            if(thisDm < minDm) {
              hig_ind[0] = ind0; hig_ind[1] = ind1; hig_ind[2] = ind2; hig_ind[3] = ind3;
              minDm = thisDm;
            }
          } // ind3
        } // ind2
      } // ind1
    } // ind0
    baby.jets_h1()[hi_csv[hig_ind[0]]] = true; baby.jets_h1()[hi_csv[hig_ind[1]]] = true; 
    baby.jets_h2()[hi_csv[hig_ind[2]]] = true; baby.jets_h2()[hi_csv[hig_ind[3]]] = true; 

    LVector hig1 = hig_p4[hig_ind[0]][hig_ind[1]], hig2 = hig_p4[hig_ind[2]][hig_ind[3]];
    baby.hig_dm()   = fabs(hig1.mass() - hig2.mass());
    baby.hig_am()   = (hig1.mass() + hig2.mass())/2.;
    baby.hig_drmax() = max(deltaR(jets20[hi_csv[hig_ind[0]]], jets20[hi_csv[hig_ind[1]]]),
                           deltaR(jets20[hi_csv[hig_ind[2]]], jets20[hi_csv[hig_ind[3]]]));

    baby.hig1_pt()  = hig1.pt();
    baby.hig1_eta() = hig1.eta();
    baby.hig1_phi() = hig1.phi();
    baby.hig1_m()   = hig1.mass();
    baby.hig2_pt()  = hig2.pt();
    baby.hig2_eta() = hig2.eta();
    baby.hig2_phi() = hig2.phi();
    baby.hig2_m()   = hig2.mass();

    // Setting up the ABCD bin: 
    // 2 -> SIG, 1 -> SB, 0 -> in between, not used
    if(baby.hig_dm()<40 && baby.hig_am()>100 && baby.hig_am()<140) baby.hig_bin() = 2;
    else if(baby.hig_dm()>40 || baby.hig_am()<100 || baby.hig_am()>140) baby.hig_bin() = 1;
    else baby.hig_bin() = 0;
    // 20 -> 2b, 30 -> 3b, 40 -> 4b
    if(baby.nbt20()>=2) {
      baby.hig_bin() += 20;
      if(baby.nbm20()>=3) {
        baby.hig_bin() += 10;
        if(baby.nbl20()>=4) baby.hig_bin() += 10;
      }
    }

    if(hi_csv[4]>=0) nCSVs = 5; // If there are 5 jets, find all combinations of those 5
    minDm = 9999.;
    for(size_t ind0=0; ind0<nCSVs; ind0++){
      for(size_t ind1=0; ind1<ind0; ind1++){
        for(size_t ind2=ind0+1; ind2<nCSVs; ind2++){
          if(ind2==ind1) continue;
          for(size_t ind3=0; ind3<ind2; ind3++){
            if(ind3==ind0 || ind3==ind1) continue;
            float thisDm = fabs(hig_p4[ind0][ind1].mass() - hig_p4[ind2][ind3].mass());
            //cout<<ind0<<"-"<<ind1<<", "<<ind2<<"-"<<ind3<<", thisDM "<<thisDm<<", minDm "<<minDm<<endl;
            if(thisDm < minDm) {
              hig_ind[0] = ind0; hig_ind[1] = ind1; hig_ind[2] = ind2; hig_ind[3] = ind3;
              minDm = thisDm;
            }
          } // ind3
        } // ind2
      } // ind1
    } // ind0
    hig1 = hig_p4[hig_ind[0]][hig_ind[1]]; hig2 = hig_p4[hig_ind[2]][hig_ind[3]];
    baby.hig5_dm()   = fabs(hig1.mass() - hig2.mass());
    baby.hig5_am()   = (hig1.mass() + hig2.mass())/2.;
    baby.hig5_drmax() = max(deltaR(jets20[hi_csv[hig_ind[0]]], jets20[hi_csv[hig_ind[1]]]),
                            deltaR(jets20[hi_csv[hig_ind[2]]], jets20[hi_csv[hig_ind[3]]]));


  } // if njets >= 4
  return clean_jets;
} // writeJets

void bmaker_full::writeFatJets(vector<LVector> &jets){
  jetTool->clusterFatJets(baby.nfjets(), baby.mj(),
                          baby.fjets_pt(), baby.fjets_eta(),
                          baby.fjets_phi(), baby.fjets_m(),
                          baby.fjets_nconst(),
                          baby.fjets_sumcsv(), baby.fjets_poscsv(),
                          baby.fjets_btags(), baby.jets_fjet_index(),
                          1.2, jets, baby.jets_csv());

  jetTool->clusterFatJets(baby.nfjets08(), baby.mj08(),
                          baby.fjets08_pt(), baby.fjets08_eta(),
                          baby.fjets08_phi(), baby.fjets08_m(),
                          baby.fjets08_nconst(),
                          baby.fjets08_sumcsv(), baby.fjets08_poscsv(),
                          baby.fjets08_btags(), baby.jets_fjet08_index(),
                          0.8, jets, baby.jets_csv());

  jetTool->clusterFatJets(baby.nfjets13(), baby.mj13(),
                          baby.fjets13_pt(), baby.fjets13_eta(),
                          baby.fjets13_phi(), baby.fjets13_m(),
                          baby.fjets13_nconst(),
                          baby.fjets13_sumcsv(), baby.fjets13_poscsv(),
                          baby.fjets13_btags(), baby.jets_fjet13_index(),
                          1.3, jets, baby.jets_csv());

  jetTool->clusterFatJets(baby.nfjets14(), baby.mj14(),
                          baby.fjets14_pt(), baby.fjets14_eta(),
                          baby.fjets14_phi(), baby.fjets14_m(),
                          baby.fjets14_nconst(),
                          baby.fjets14_sumcsv(), baby.fjets14_poscsv(),
                          baby.fjets14_btags(), baby.jets_fjet14_index(),
                          1.4, jets, baby.jets_csv());

  jetTool->clusterFatJets(baby.nfjets15(), baby.mj15(),
                          baby.fjets15_pt(), baby.fjets15_eta(),
                          baby.fjets15_phi(), baby.fjets15_m(),
                          baby.fjets15_nconst(),
                          baby.fjets15_sumcsv(), baby.fjets15_poscsv(),
                          baby.fjets15_btags(), baby.jets_fjet15_index(),
                          1.5, jets, baby.jets_csv());

  jetTool->clusterFatJets(baby.nfjets16(), baby.mj16(),
                          baby.fjets16_pt(), baby.fjets16_eta(),
                          baby.fjets16_phi(), baby.fjets16_m(),
                          baby.fjets16_nconst(),
                          baby.fjets16_sumcsv(), baby.fjets16_poscsv(),
                          baby.fjets16_btags(), baby.jets_fjet16_index(),
                          1.6, jets, baby.jets_csv());

}

vCands bmaker_full::writeMuons(edm::Handle<pat::MuonCollection> muons, 
                               edm::Handle<pat::PackedCandidateCollection> pfcands, 
                               edm::Handle<reco::VertexCollection> vtx,
                               vCands &veto_mus, vCands &all_mus, double rhoEventCentral){
  vCands sig_mus; 
  veto_mus.clear(); all_mus.clear();
  baby.nmus() = 0; baby.nvmus() = 0;
  for (size_t ilep(0); ilep < muons->size(); ilep++) {
    const pat::Muon &lep = (*muons)[ilep];    
    if(!lepTool->isVetoMuon(lep, vtx, -99.)) continue; // Storing leptons that pass all veto cuts except for iso

    double lep_iso(lepTool->getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., rhoEventCentral, false));
    double lep_reliso(lepTool->getRelIsolation(lep, rhoEventCentral));
    double dz(0.), d0(0.);
    lepTool->vertexMuon(lep, vtx, dz, d0); // Calculating dz and d0

    baby.mus_pt().push_back(lep.pt());
    baby.mus_eta().push_back(lep.eta());
    baby.mus_phi().push_back(lep.phi());
    baby.mus_dz().push_back(dz);
    baby.mus_d0().push_back(d0);
    baby.mus_charge().push_back(lep.charge());
    baby.mus_sigid().push_back(lepTool->idMuon(lep, vtx, lepTool->kMedium));
    baby.mus_tight().push_back(lepTool->idMuon(lep, vtx, lepTool->kTight));
    baby.mus_miniso().push_back(lep_iso);
    baby.mus_reliso().push_back(lep_reliso);
    baby.mus_tm().push_back(false);        // Filled in writeMC
    baby.mus_inz().push_back(false);       // Filled in writeDiLep
    baby.mus_inzv().push_back(false);      // Filled in writeDiLep
    baby.mus_vvvl().push_back(false);      // Filled in writeHLTObjects
    baby.mus_isomu18().push_back(false);   // Filled in writeHLTObjects
    baby.mus_mu50().push_back(false);      // Filled in writeHLTObjects
    baby.mus_mu8().push_back(false);       // Filled in writeHLTObjects
    if (lep.track().isNonnull()){
      baby.mus_trk_quality().push_back(lep.innerTrack()->quality(reco::TrackBase::highPurity));
      baby.mus_pterr().push_back(lep.innerTrack()->ptError());
      baby.mus_trk_nholes_in().push_back(lep.innerTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      baby.mus_trk_nholes_out().push_back(lep.innerTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS));
      baby.mus_trk_algo().push_back(lep.innerTrack()->algo());
    } else {
      baby.mus_trk_quality().push_back(-99999);
      baby.mus_pterr().push_back(-99999);
      baby.mus_trk_nholes_in().push_back(-99999);
      baby.mus_trk_nholes_out().push_back(-99999);
      baby.mus_trk_algo().push_back(-99999);
    }
    baby.mus_em_e().push_back(lep.calEnergy().em);
    baby.mus_had_e().push_back(lep.calEnergy().had);
    all_mus.push_back(dynamic_cast<const reco::Candidate *>(&lep)); // For truth-matching in writeMC

    if(lepTool->isVetoMuon(lep, vtx, lep_iso)) {
      baby.nvmus()++;
      veto_mus.push_back(dynamic_cast<const reco::Candidate *>(&lep));
    }
    if(lepTool->isSignalMuon(lep, vtx, lep_iso)) {
      baby.nmus()++;
      sig_mus.push_back(dynamic_cast<const reco::Candidate *>(&lep));
      baby.mus_sig().push_back(true); 
    } else baby.mus_sig().push_back(false); 
  } // Loop over muons
  
  return sig_mus;
}


vCands bmaker_full::writeElectrons(edm::Handle<pat::ElectronCollection> electrons, 
                                   edm::Handle<pat::PackedCandidateCollection> pfcands, 
                                   edm::Handle<reco::VertexCollection> vtx,
                                   vCands &veto_els, vCands &all_els, double rhoEventCentral){
  vCands sig_els; 
  veto_els.clear(); all_els.clear();
  baby.nels() = 0; baby.nvels() = 0;
  for (size_t ilep(0); ilep < electrons->size(); ilep++) {
    const pat::Electron &lep = (*electrons)[ilep];    
    if(!lepTool->isVetoElectron(lep, vtx, -99.)) continue; // Storing leptons that pass all veto cuts except for iso

    double lep_iso(lepTool->getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., rhoEventCentral, false));
    double lep_reliso(lepTool->getRelIsolation(lep, rhoEventCentral));
    double dz(0.), d0(0.);
    lepTool->vertexElectron(lep, vtx, dz, d0); // Calculating dz and d0

    baby.els_pt().push_back(lep.pt());
    baby.els_sceta().push_back(lep.superCluster()->position().eta());
    baby.els_eta().push_back(lep.eta());
    baby.els_phi().push_back(lep.phi());
    baby.els_dz().push_back(dz);
    baby.els_d0().push_back(d0);
    baby.els_ip3d().push_back(lep.ip3d());
    baby.els_charge().push_back(lep.charge());
    baby.els_sigid().push_back(lepTool->idElectron(lep, vtx, lepTool->kMedium));
    baby.els_ispf().push_back(lep.numberOfSourceCandidatePtrs()==2 && abs(lep.sourceCandidatePtr(1)->pdgId())==11);
    baby.els_tight().push_back(lepTool->idElectron(lep, vtx, lepTool->kTight));
    baby.els_miniso().push_back(lep_iso);
    baby.els_reliso().push_back(lep_reliso);
    baby.els_tm().push_back(false);       // Filled in writeMC
    baby.els_inz().push_back(false);      // Filled in writeDiLep
    baby.els_inzv().push_back(false);     // Filled in writeDiLep
    baby.els_vvvl().push_back(false);     // Filled in writeHLTObjects
    baby.els_ele23().push_back(false);    // Filled in writeHLTObjects
    baby.els_ele105().push_back(false);   // Filled in writeHLTObjects
    baby.els_ele8().push_back(false);     // Filled in writeHLTObjects
    baby.els_hovere().push_back(lep.hadronicOverEm());
    baby.els_eoverp().push_back(lep.eSuperClusterOverP());
    baby.els_em_e().push_back(lep.ecalEnergy());
    baby.els_deta_sctrk().push_back(lep.deltaEtaSuperClusterTrackAtVtx());
    baby.els_dphi_sctrk().push_back(lep.deltaPhiSuperClusterTrackAtVtx());
    if (lep.gsfTrack().isNonnull()){
      baby.els_trk_pt().push_back(lep.gsfTrack()->pt());
      baby.els_trk_pterr().push_back(lep.gsfTrack()->ptError());
      baby.els_trk_nholes().push_back(lep.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    } else {
      baby.els_trk_pt().push_back(-99999);
      baby.els_trk_pterr().push_back(-99999);
      baby.els_trk_nholes().push_back(-99999);
    }
    all_els.push_back(dynamic_cast<const reco::Candidate *>(&lep)); // For truth-matching in writeMC

    if(lepTool->isVetoElectron(lep, vtx, lep_iso)){
      baby.nvels()++;
      veto_els.push_back(dynamic_cast<const reco::Candidate *>(&lep));
    }
    if(lepTool->isSignalElectron(lep, vtx, lep_iso)) {
      baby.nels()++;
      sig_els.push_back(dynamic_cast<const reco::Candidate *>(&lep));
      baby.els_sig().push_back(true); 
    } else baby.els_sig().push_back(false); 
  } // Loop over electrons

  return sig_els;
}

void bmaker_full::writeLeptons(vCands &leptons){ 
  baby.nleps() = baby.nmus() + baby.nels();
  baby.nvleps() = baby.nvmus() + baby.nvels();
  sort(leptons.begin(), leptons.end(), greaterPt);
  for(size_t ind(0); ind < leptons.size(); ind++) {
    baby.leps_pt().push_back(leptons[ind]->pt());
    baby.leps_eta().push_back(leptons[ind]->eta());
    baby.leps_phi().push_back(leptons[ind]->phi());
    baby.leps_id().push_back(leptons[ind]->pdgId());
  } // Loop over leptons
}

void bmaker_full::writeTks(edm::Handle<pat::PackedCandidateCollection> pfcands,edm::Handle<reco::VertexCollection> vtx, double rhoEventCentral ){
  if(baby.leps_id().size()>0){    
    vector<float> isos;
    vector<float> relisos;
    vCands ra4tks;
    ra4tks = lepTool->getRA4IsoTracks(pfcands, baby.met(), baby.met_phi(),rhoEventCentral,isos,relisos,baby.leps_id().at(0));
     
    int nveto=0;

    for(unsigned i=0;i<ra4tks.size();i++){
      baby.tks_pt().push_back(ra4tks.at(i)->pt());
      baby.tks_eta().push_back(ra4tks.at(i)->eta());
      baby.tks_phi().push_back(ra4tks.at(i)->phi());
      baby.tks_d0().push_back( sqrt(pow(ra4tks.at(i)->vx()-vtx->at(0).x(),2) + pow(vtx->at(0).y()-ra4tks.at(i)->vy(),2)));
      baby.tks_dz().push_back(ra4tks.at(i)->vz()-vtx->at(0).z());
      baby.tks_pdg().push_back(ra4tks.at(i)->pdgId());
      baby.tks_miniso().push_back(isos.at(i));
      baby.tks_reliso().push_back(relisos.at(i));
      baby.tks_tm().push_back(false); //filled in writeMC
                                                        
      baby.tks_mt2().push_back(getMT2(baby.leps_pt().at(0),baby.leps_phi().at(0),baby.tks_pt().back(),baby.tks_phi().back(),baby.met(),baby.met_phi()));
      baby.tks_mt().push_back(getMT(baby.tks_pt().back(),baby.tks_phi().back(),baby.met(),baby.met_phi()));
       
      if(fabs(baby.tks_pdg().back())==211  && baby.tks_pt().back()>15. && baby.tks_miniso().back()<0.1 && baby.tks_mt2().back()<60 && baby.tks_dz().back()<0.07 && baby.tks_d0().back()<0.05 ) nveto++;
      else if (fabs(baby.tks_pdg().back())==13 && baby.tks_pt().back()>10. && baby.tks_miniso().back()<0.2 && baby.tks_mt2().back()<80 && baby.tks_dz().back()<0.07 && baby.tks_d0().back()<0.05) nveto++;
      else if (fabs(baby.tks_pdg().back())==11 && baby.tks_pt().back()>10. && baby.tks_miniso().back()<0.2 && baby.tks_mt2().back()<80 && baby.tks_dz().back()<0.07 && baby.tks_d0().back()<0.05) nveto++;

    }
 
    baby.nveto()=nveto;

  }  
} // if goodPV

 
void bmaker_full::writeDiLep(vCands &sig_mus, vCands &sig_els, vCands &veto_mus, vCands &veto_els){
  setDiLepMass(sig_mus,  &baby_base::mumu_m,  &baby_base::mumu_pt1,  &baby_base::mumu_pt2,  &baby_base::mumu_pt,
               &baby_base::mumu_eta,  &baby_base::mumu_phi, &baby_base::mus_pt, &baby_base::mus_inz,
               &baby_base::mumu_w);
  setDiLepMass(veto_mus, &baby_base::mumuv_m, &baby_base::mumuv_pt1, &baby_base::mumuv_pt2, &baby_base::mumuv_pt,
               &baby_base::mumuv_eta,  &baby_base::mumuv_phi, &baby_base::mus_pt, &baby_base::mus_inzv,
               &baby_base::mumuv_w);
  setDiLepMass(sig_els,  &baby_base::elel_m,  &baby_base::elel_pt1,  &baby_base::elel_pt2,  &baby_base::elel_pt,
               &baby_base::elel_eta,  &baby_base::elel_phi, &baby_base::els_pt, &baby_base::els_inz,
               &baby_base::elel_w);
  setDiLepMass(veto_els, &baby_base::elelv_m, &baby_base::elelv_pt1, &baby_base::elelv_pt2, &baby_base::elelv_pt,
               &baby_base::elelv_eta,  &baby_base::elelv_phi, &baby_base::els_pt, &baby_base::els_inzv,
               &baby_base::elelv_w);
  setElMuMass(sig_els, sig_mus, &baby_base::elmu_m, &baby_base::elmu_pt1, &baby_base::elmu_pt2, &baby_base::elmu_pt,
              &baby_base::elmu_eta,  &baby_base::elmu_phi,
              &baby_base::elmu_w);
  // setElMuMass(veto_els, veto_mus, &baby_base::elmuv_m, &baby_base::elmuv_pt1, &baby_base::elmuv_pt2, &baby_base::elmuv_pt,
  //          &baby_base::elmuv_eta,  &baby_base::elmuv_phi);
}

void bmaker_full::setDiLepMass(vCands leptons, baby_float ll_m, baby_float ll_pt1, baby_float ll_pt2, 
                               baby_float ll_pt, baby_float ll_eta, baby_float ll_phi, baby_vfloat l_pt, baby_vbool l_inz,
                               baby_float ll_w){
  for(size_t lep1(0); lep1 < leptons.size(); lep1++){
    for(size_t lep2(lep1+1); lep2 < leptons.size(); lep2++){
      if(leptons[lep1]->charge()*leptons[lep2]->charge()<0){
        LVector z_p4(leptons[lep1]->p4()); 
        z_p4 += leptons[lep2]->p4();
        (baby.*ll_m)()   = z_p4.mass();
        (baby.*ll_pt)()  = z_p4.pt();
        (baby.*ll_eta)() = z_p4.eta();
        (baby.*ll_phi)() = z_p4.phi();
        float pt1(leptons[lep1]->pt()), pt2(leptons[lep2]->pt());
        (baby.*ll_pt1)() = max(pt1, pt2); 
        (baby.*ll_pt2)() = min(pt1, pt2);
        for(size_t ilep(0); ilep < (baby.*l_pt)().size(); ilep++){
          if(fabs(pt1 - (baby.*l_pt)()[ilep]) < 1e-7) (baby.*l_inz)()[ilep] = true;
          if(fabs(pt2 - (baby.*l_pt)()[ilep]) < 1e-7) (baby.*l_inz)()[ilep] = true;
        }
        (baby.*ll_w)() = lepton_tools::getScaleFactor({leptons[lep1], leptons[lep2]}).first;
        return; // We only set it with the first good ll combination
      }
    } // Loop over lep2
  } // Loop over lep1
}

void bmaker_full::setElMuMass(vCands leptons1, vCands leptons2, baby_float ll_m, baby_float ll_pt1, baby_float ll_pt2, 
                              baby_float ll_pt, baby_float ll_eta, baby_float ll_phi,
                              baby_float ll_w){
  for(size_t lep1(0); lep1 < leptons1.size(); lep1++){
    for(size_t lep2(0); lep2 < leptons2.size(); lep2++){
      if(leptons1[lep1]->charge()*leptons2[lep2]->charge()<0){
        LVector z_p4(leptons1[lep1]->p4()); 
        z_p4 += leptons2[lep2]->p4();
        (baby.*ll_m)()   = z_p4.mass();
        (baby.*ll_pt)()  = z_p4.pt();
        (baby.*ll_eta)() = z_p4.eta();
        (baby.*ll_phi)() = z_p4.phi();
        float pt1(leptons1[lep1]->pt()), pt2(leptons2[lep2]->pt());
        (baby.*ll_pt1)() = pt1; 
        (baby.*ll_pt2)() = pt2;
        (baby.*ll_w)() = lepton_tools::getScaleFactor({leptons1[lep1], leptons2[lep2]}).first;
        return; // We only set it with the first good ll combination
      }
    } // Loop over lep2
  } // Loop over lep1
}


vCands bmaker_full::writePhotons(edm::Handle<pat::PhotonCollection> allphotons, 
                                 edm::Handle<std::vector<pat::Electron> > &electrons,
                                 edm::Handle<reco::ConversionCollection> &conversions,
                                 edm::Handle<reco::BeamSpot> &beamspot, double rho){
  vCands photons;
  baby.nph() = 0;

  for (size_t ind(0); ind < allphotons->size(); ind++) {
    const pat::Photon &photon = (*allphotons)[ind];    

    if(photon.pt() < 50) continue;
    if(!photonTool->idPhoton(photon, electrons, conversions, beamspot, rho)) continue;

    if(photon.pt() > photonTool->PhotonPtCut) baby.nph()++;
    baby.ph_pt().push_back(photon.pt());
    baby.ph_eta().push_back(photon.eta());
    baby.ph_phi().push_back(photon.phi());
    baby.ph_tm().push_back(false);          // Filled in writeMC
    baby.ph_ph90().push_back(false);    // Filled in writeHLTObjects    

    photons.push_back(dynamic_cast<const reco::Candidate *>(&photon)); 
  } // Loop over allphotons

  return photons;
} // writePhotons

bool bmaker_full::writeTriggers(const edm::TriggerNames &names, 
                                edm::Handle<edm::TriggerResults> triggerBits, 
                                edm::Handle<pat::PackedTriggerPrescales> triggerPrescales){
  bool keep(false);
  baby.trig().resize(trig_name.size(), false);
  baby.trig_prescale().resize(trig_name.size(), -1.);
  for (size_t itrig(0); itrig < triggerBits->size(); itrig++) {
    for(size_t itn(0); itn < trig_name.size(); itn++){
      if(names.triggerName(itrig).find(trig_name[itn])!=string::npos){
        baby.trig()[itn] = triggerBits->accept(itrig);
        baby.trig_prescale()[itn] = triggerPrescales->getPrescaleForIndex(itrig);
        if(baby.trig()[itn]) keep = true;
      }
    }
  }
  // trig_ra4 = (Lep15_HT350 ||  Lep15_HT400 || MET100 || METNoMu100)
  vector<size_t> trig_or = {3,4,7,8,13,33};
  baby.trig_ra4() = false;
  for(size_t itrig=0; itrig<trig_or.size(); itrig++)
    baby.trig_ra4() = (baby.trig_ra4() || baby.trig()[trig_or[itrig]]);
  
  return keep;
}

void bmaker_full::writeHLTObjects(const edm::TriggerNames &names, 
                                  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects, 
                                  vCands &all_mus, vCands &all_els, vCands &photons){
  baby.nmus_vvvl() = 0; baby.nmus_isomu18() = 0; 
  baby.nels_vvvl() = 0; baby.nels_ele23() = 0; 
  const float relptThreshold(1), drThreshold(0.3);      
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackPathNames(names);
    TString name(obj.collection());
    float objpt(obj.pt());
    if(name=="hltPFMETProducer::HLT") baby.onmet() = objpt;
    if(name=="hltPFHT::HLT") baby.onht() = objpt; // There's 2 of these, and we want the last one

    if(name=="hltL3MuonCandidates::HLT"){
      bool vvvl(obj.hasFilterLabel("hltL3MuVVVLIsoFIlter"));
      bool isomu18(obj.hasFilterLabel("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09") );
      bool mu50(obj.hasFilterLabel("hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered50Q"));
      bool mu8(obj.hasFilterLabel("hltDoubleMu8Mass8L3Filtered"));
      if(vvvl) {
        baby.nmus_vvvl()++;
        baby.mus_vvvl_pt().push_back(objpt);
        baby.mus_vvvl_eta().push_back(obj.eta());
        baby.mus_vvvl_phi().push_back(obj.phi());
      }
      if(isomu18) baby.nmus_isomu18()++;
      if(vvvl && baby.onmu_vvvl()<objpt) baby.onmu_vvvl() = objpt;
      if(isomu18 && baby.onmu_isomu18()<objpt) baby.onmu_isomu18() = objpt;
      if(mu50 && baby.onmu_mu50()<objpt) baby.onmu_mu50() = objpt;
      if(mu8 && baby.onmu_mu8()<objpt) baby.onmu_mu8() = objpt;
      double mindr(999.);
      int minind(-1);
      if(vvvl || isomu18 || mu50 || mu8){
        for(size_t ind(0); ind < all_mus.size(); ind++) {
          double dr(deltaR(obj, *(all_mus[ind])));
          //double drelpt(fabs((all_mus[ind]->pt() - objpt)/objpt));
          //if(dr > drThreshold || drelpt > relptThreshold) continue;
          if(dr > drThreshold) continue;
          if(dr < mindr){
            mindr = dr;
            minind = ind;
          }
        } // Loop over reco muons
        if(minind>=0){
          baby.mus_vvvl()[minind] = vvvl;
          baby.mus_isomu18()[minind] = isomu18;
          baby.mus_mu50()[minind] = mu50;
          baby.mus_mu8()[minind] = mu8;
        }
      } // At least one match
    }
    if(name=="hltEgammaCandidates::HLT" || name=="hltDoubleEle8HLTPixelMatchElectronProducer::HLT"){
      bool vvvl(obj.hasFilterLabel("hltEle15VVVLGsfTrackIsoFilter"));
      bool ele23(obj.hasFilterLabel("hltEle23WPLooseGsfTrackIsoFilter") );
      bool ele105(obj.hasFilterLabel("hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter"));
      bool ele8(obj.hasFilterLabel("hltDoubleEle8Mass8Filter"));
      if(vvvl) {
        baby.nels_vvvl()++;
        baby.els_vvvl_pt().push_back(objpt);
        baby.els_vvvl_eta().push_back(obj.eta());
        baby.els_vvvl_phi().push_back(obj.phi());
      }
      if(ele23) baby.nels_ele23()++;
      if(vvvl && baby.onel_vvvl()<objpt) baby.onel_vvvl() = objpt;
      if(ele23 && baby.onel_ele23()<objpt) baby.onel_ele23() = objpt;
      if(ele105 && baby.onel_ele105()<objpt) baby.onel_ele105() = objpt;
      if(ele8 && baby.onel_ele8()<objpt) baby.onel_ele8() = objpt;
      if(vvvl || ele23 || ele105 || ele8){
        double mindr(999.);
        int minind(-1);
        for(size_t ind(0); ind < all_els.size(); ind++) {
          double dr(deltaR(obj, *(all_els[ind])));
          double drelpt(fabs((all_els[ind]->pt() - objpt)/objpt));
          if(dr > drThreshold || drelpt > relptThreshold) continue;
          if(dr < mindr){
            mindr = dr;
            minind = ind;
          }
        } // Loop over reco elecrons
        if(minind>=0){
          baby.els_vvvl()[minind] = vvvl;
          baby.els_ele23()[minind] = ele23;
          baby.els_ele105()[minind] = ele105;
          baby.els_ele8()[minind] = ele8;
        }
      } // At least one electron match
      bool ph90(obj.hasFilterLabel("hltEG90L1SingleEG40HEFilter"));
      if(ph90 && baby.onph_ph90()<objpt) baby.onph_ph90() = objpt;
      if(ph90){
        double mindr(999.);
        int minind(-1);
        for(size_t ind(0); ind < photons.size(); ind++) {
          double dr(deltaR(obj, *(photons[ind])));
          double drelpt(fabs((photons[ind]->pt() - objpt)/objpt));
          if(dr > drThreshold || drelpt > relptThreshold) continue;
          if(dr < mindr){
            mindr = dr;
            minind = ind;
          }
        } // Loop over reco photons
        if(minind>=0){
          baby.ph_ph90()[minind] = ph90;
        }
      } // At least one photon match
    }
  } // Loop over trigger objects
}

// From https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
void bmaker_full::writeFilters(const edm::TriggerNames &fnames,
                               edm::Handle<edm::TriggerResults> filterBits,
                               edm::Handle<reco::VertexCollection> vtx,
                               vector<double> jetsMuonEnergyFrac){
  baby.pass_goodv() = true; baby.pass_cschalo() = true; baby.pass_eebadsc() = true;
  baby.pass_ecaldeadcell() = true; baby.pass_hbhe() = true; baby.pass_hbheiso() = true;
  baby.pass_ra2_badmu() = true;
  for (size_t i(0); i < filterBits->size(); ++i) {
    string name = fnames.triggerName(i);
    bool pass = static_cast<bool>(filterBits->accept(i));
    if (name=="Flag_goodVertices") baby.pass_goodv() = pass;
    //else if (name=="Flag_CSCTightHaloFilter") baby.pass_cschalo() = pass; // Requires reading it from txt file
    else if (name=="Flag_CSCTightHalo2015Filter") baby.pass_cschalo() = pass;
    else if (name=="Flag_eeBadScFilter") baby.pass_eebadsc() = pass;
    else if (name=="Flag_EcalDeadCellTriggerPrimitiveFilter") baby.pass_ecaldeadcell() = pass;
    else if (name=="Flag_HBHENoiseFilter") baby.pass_hbhe() = pass; 
    else if (name=="Flag_HBHENoiseIsoFilter") baby.pass_hbheiso() = pass; 
  }

  //baby.pass_goodv() &= eventTool->hasGoodPV(vtx); // We needed to re-run it for Run2015B
  //baby.pass_cschalo() = eventTool->passBeamHalo(baby.run(), baby.event()); // now taken from miniAOD

  baby.pass() = baby.pass_goodv() && baby.pass_eebadsc() && baby.pass_cschalo() && baby.pass_hbhe() && baby.pass_hbheiso() && baby.pass_ecaldeadcell() 
    && baby.pass_jets();
  baby.pass_ra2() = baby.pass_goodv() && baby.pass_eebadsc() && baby.pass_cschalo() && baby.pass_hbhe() && baby.pass_hbheiso() && baby.pass_ecaldeadcell()  
    && baby.pass_jets_ra2();
  baby.pass_nohf() = baby.pass_goodv() && baby.pass_eebadsc() && baby.pass_cschalo() && baby.pass_hbhe() && baby.pass_hbheiso() && baby.pass_ecaldeadcell() 
    && baby.pass_jets_nohf();

  for (size_t ijet(0); ijet < baby.jets_pt().size(); ijet++){
    if (abs(baby.jets_eta()[ijet])>2.4 || baby.jets_pt()[ijet]<200.) continue;
    if (baby.jets_islep()[ijet]) continue;
    if (jetsMuonEnergyFrac[ijet]<=0.5) continue;
    if (abs(dPhi(baby.jets_phi()[ijet],baby.met_phi()))<(3.14159-0.4)) continue;
    baby.pass_ra2_badmu() = false;
    break;
  }

  if (doSystematics){
    for (unsigned isys(0); isys<kSysLast; isys++){
      // sys_pass_jets already stored in the value of this variable in the baby
      baby.sys_pass()[isys] = baby.sys_pass()[isys] && baby.pass_goodv() && baby.pass_eebadsc() && baby.pass_cschalo() && baby.pass_hbhe() && baby.pass_hbheiso() && baby.pass_ecaldeadcell();
    }
  }
}

void bmaker_full::writeVertices(edm::Handle<reco::VertexCollection> vtx,
                                edm::Handle<std::vector< PileupSummaryInfo > >  pu_info){  
  baby.npv() = vtx->size();
  if(pu_info.isValid()){
    for(size_t bc(0); bc<pu_info->size(); ++bc){
      if(pu_info->at(bc).getBunchCrossing()==0){
        baby.ntrupv() = pu_info->at(bc).getPU_NumInteractions();
        baby.ntrupv_mean() = pu_info->at(bc).getTrueNumInteractions();
        break;
      }
    } // Loop over true PVs
  } // if pu_info is valid
}

void bmaker_full::writeGenInfo(edm::Handle<LHEEventProduct> lhe_info){
  baby.nisr_me()=0; baby.ht_isr_me()=0.; 
  for ( unsigned int icount = 0 ; icount < (unsigned int)lhe_info->hepeup().NUP; icount++ ) {
    unsigned int pdgid = abs(lhe_info->hepeup().IDUP[icount]);
    int status = lhe_info->hepeup().ISTUP[icount];
    int mom1id = abs(lhe_info->hepeup().IDUP[lhe_info->hepeup().MOTHUP[icount].first-1]);
    int mom2id = abs(lhe_info->hepeup().IDUP[lhe_info->hepeup().MOTHUP[icount].second-1]);
    float px = (lhe_info->hepeup().PUP[icount])[0];
    float py = (lhe_info->hepeup().PUP[icount])[1];
    float pt = sqrt(px*px+py*py);

    if(status==1 && (pdgid<6 || pdgid==21) && mom1id!=6 && mom2id!=6 && mom1id!=24 && mom2id!=24 
       && mom1id!=23 && mom2id!=23 && mom1id!=25 && mom2id!=25) {
      baby.nisr_me()++;
      baby.ht_isr_me() += pt;
    }

  } // Loop over generator particles
  
  if(outname.Contains("SMS")){ //Get mgluino and mlsp
    
    typedef std::vector<std::string>::const_iterator comments_const_iterator;
    
    comments_const_iterator c_begin = lhe_info->comments_begin();
    comments_const_iterator c_end = lhe_info->comments_end();
    
    TString model_params;
    for(comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
      size_t found = (*cit).find("model");
      if(found != std::string::npos)   {
        //    std::cout <<"BABYMAKER: "<< *cit <<"end"<< std::endl;  
        model_params = *cit;
      }
    }
    mcTool->getMassPoints(model_params,baby.mgluino(),baby.mlsp());
  }
} // writeGenInfo

void bmaker_full::nisrMatch(edm::Handle<reco::GenParticleCollection> genParticles, 
                            vector<reco::Candidate::LorentzVector> &clean_jets){
  bool verbose = false;
  int nisr(0);
  for (size_t ijet(0); ijet<clean_jets.size(); ijet++){

    bool matched=false;
    for (size_t imc(0); imc < genParticles->size(); imc++) {
      if (matched) break;
      const reco::GenParticle &mc = (*genParticles)[imc];
      if (mc.status()!=23 || abs(mc.pdgId())>5) continue;
      int momid = abs(mc.mother()->pdgId());
      if(!(momid==6 || momid==23 || momid==24 || momid==25 || momid>1e6)) continue; 
      //check against daughter in case of hard initial splitting
      for (size_t idau(0); idau < mc.numberOfDaughters(); idau++) {
        float dR = deltaR(clean_jets[ijet], mc.daughter(idau)->p4());
        if(dR<0.3){
          if (verbose) cout<<"Jet: ("<<clean_jets[ijet].pt()<<", "<<clean_jets[ijet].eta()<<", "<<clean_jets[ijet].phi()
            <<"), MC: ("<<mc.daughter(idau)->pt()<<", "<<mc.daughter(idau)->eta()<<", "<<mc.daughter(idau)->phi()<<"), ID "<<mc.daughter(idau)->pdgId()<<". dR "<<dR <<endl;
            matched = true;
            break;
        }
      }
    } // Loop over MC particles
    if(!matched) {
      nisr++;
    }
  } // Loop over jets

  // int nisr=0;
  // set<int> unmatched;
  // for (size_t ijet(0); ijet<baby.jets_pt().size(); ijet++){
  //   bool goodjet = baby.jets_pt()[ijet] > jetTool->JetPtCut && fabs(baby.jets_eta()[ijet]) <= jetTool->JetEtaCut && !baby.jets_islep()[ijet];
  //   if(!goodjet) continue;
  //   bool matched=false;
  //   for (size_t imc(0); imc<baby.mc_pt().size(); imc++){
  //     if(baby.mc_status()[imc]!=23 || abs(baby.mc_id()[imc])>5) continue;
  //     // In our ntuples where all taus come from W
  //     if(!(abs(baby.mc_mom()[imc])==6 || abs(baby.mc_mom()[imc])==23 || 
  //         abs(baby.mc_mom()[imc])==24 || abs(baby.mc_mom()[imc])==15 || abs(baby.mc_mom()[imc])>1e6)) continue; 
  //     float dR = deltaR(baby.jets_eta()[ijet], baby.jets_phi()[ijet], baby.mc_eta()[imc], baby.mc_phi()[imc]);
  //     if(dR<0.4){
  //       if (verbose) cout<<"Jet: ("<<baby.jets_pt()[ijet]<<", "
  //         <<baby.jets_eta()[ijet]<<", "<<baby.jets_phi()[ijet]
  //         <<"), MC: ("<<baby.mc_pt()[imc]<<", "<<baby.mc_eta()[imc]<<", "
  //         <<baby.mc_phi()[imc]<<"), ID "<<baby.mc_id()[imc]<<". dR "<<dR <<endl;
  //         matched = true;
  //         break;
  //     }
  //   } // Loop over MC particles
  //   if(!matched) {
  //     nisr++;
  //     unmatched.insert(ijet);
  //   }
  // } // Loop over jets
  // if (verbose && nisr>0){
  //   for (size_t ijet(0); ijet<baby.jets_pt().size(); ijet++){
  //     if (unmatched.find(ijet)==unmatched.end()) continue;
  //     cout<<"Jet: ("<<baby.jets_pt()[ijet]<<", "<<baby.jets_eta()[ijet]
  //                      <<", "<<baby.jets_phi()[ijet]<<" not matched"<<endl;
  //     for (size_t imc(0); imc < genParticles->size(); imc++) {
  //       const reco::GenParticle &mc = (*genParticles)[imc];
  //       float dR = deltaR(baby.jets_eta()[ijet], baby.jets_phi()[ijet], mc.eta(), mc.phi());
  //       if (dR<0.6){
  //         cout<<"MC particle with dR = "<<dR<<endl;
  //         mcTool->printParticle(mc);
  //       }
  //     }
  //   }
    // for (size_t imc(0); imc<baby.mc_pt().size(); imc++){
    //   if(baby.mc_status()[imc]!=23 || abs(baby.mc_id()[imc])>5) continue;
    //   if(!(abs(baby.mc_mom()[imc])==6 || abs(baby.mc_mom()[imc])==23 || 
    //       abs(baby.mc_mom()[imc])==24 || abs(baby.mc_mom()[imc])==15 || abs(baby.mc_mom()[imc])>1e6)) continue; 
    //   cout<<" MC: ("<<baby.mc_pt()[imc]<<", "<<baby.mc_eta()[imc]<<", "
    //       <<baby.mc_phi()[imc]<<"), ID "<<baby.mc_id()[imc]<<endl;
    // }
  //   cout<<" ======== New event: njets "<<baby.njets()<<", nisr "<<nisr<<endl<<endl;
  // }
  baby.nisr() = nisr;
}

void bmaker_full::writeMC(edm::Handle<reco::GenParticleCollection> genParticles, 
                          vCands &all_mus, vCands &all_els, vCands &photons){
  LVector isr_p4;
  float metw_tru_x(0.), metw_tru_y(0.);
  float lep_tru_pt(0.), lep_tru_phi(0.);
  baby.ntruleps()=0; baby.ntrumus()=0; baby.ntruels()=0; baby.ntrutaush()=0; baby.ntrutausl()=0;
  baby.nleps_tm()=0;
  baby.fromGS()=false;
  baby.m_tt()=0;
  vector<float> top_pt;
  int topIndex=-1;
  int antitopIndex=-1;
  const size_t bsmid(1000000);
  for (size_t imc(0); imc < genParticles->size(); imc++) {
    const reco::GenParticle &mc = (*genParticles)[imc];
    size_t id = abs(mc.pdgId());
    bool isLast = mcTool->isLast(mc, id);
    bool isFirst = true;
    const reco::GenParticle *mcMom = static_cast<const reco::GenParticle *>(mc.mother());
    if(mcMom){
      if(mcMom->pdgId() == mc.pdgId()) isFirst = false;
    }
    // mcTool->printParticle(mc); // Prints various properties of the MC particle

    const reco::GenParticle *mom = nullptr;
    size_t momid = abs(mcTool->mom(mc, mom));
    bool isTop(id==6);
    bool isNewPhysics(id>=bsmid);
    bool isGluino(id==1000021);
    bool isZ(id==23);
    bool isW(id==24);
    bool isH(id==25);
    bool bTopOrBSM(id==5 && (momid==6 || momid>=bsmid));
    bool nuFromZ((id==12 || id==14 || id==16) && momid==23);
    bool eFromTopZ(id==11 && (momid==24 || momid==23));
    bool muFromTopZ(id==13 && (momid==24 || momid==23));
    bool tauFromTopZ(id==15 && (momid==24 || momid==23));
    bool fromTau(mcTool->fromTau(mc));
    bool fromWOrWTau(mcTool->fromWOrWTau(mc));
    bool chgPionFromTau(id==211 && momid==15 && fromWOrWTau);
    bool eFromTopZtau(id==11 && (momid==24 || momid==23 || fromTau ));
    bool muFromTopZtau(id==13 && (momid==24 || momid==23 || fromTau ));
    bool fromHiggs(momid==25);
    bool fromTop(momid==6);
    bool mg_me = mc.status()==22 || mc.status()==23;

    //////// Saving interesting true particles
    if((isLast && (isTop || isNewPhysics || isZ || isH))
       || (isFirst && (fromHiggs || isW || bTopOrBSM || eFromTopZ || muFromTopZ 
                       || tauFromTopZ || nuFromZ || fromWOrWTau || fromTop))
       || mg_me){
      baby.mc_status().push_back(mc.status());
      baby.mc_id().push_back(mc.pdgId());
      baby.mc_pt().push_back(mc.pt());
      baby.mc_eta().push_back(mc.eta());
      baby.mc_phi().push_back(mc.phi());
      baby.mc_mass().push_back(mc.mass());
      baby.mc_mom().push_back(mcTool->mom(mc,mom));
    }

    if(isLast){
      if(isTop) mc.pdgId()>0 ? topIndex=imc : antitopIndex=imc;

      //////// Finding p4 of ME ISR system
      if((isTop && outname.Contains("TTJets"))
         || (isGluino && (outname.Contains("SMS") || outname.Contains("RPV")))
         || (isZ && outname.Contains("DY"))) isr_p4 -= mc.p4();

      if(isTop && (outname.Contains("TTJets") || outname.Contains("TT_"))){
        top_pt.push_back(mc.pt());
      }

      //////// Counting true leptons
      if(muFromTopZ) baby.ntrumus()++;
      if(eFromTopZ)  baby.ntruels()++;
      if(tauFromTopZ){
        const reco::GenParticle *tauDaughter(0);
        if(mcTool->decaysTo(mc, 11, tauDaughter) || mcTool->decaysTo(mc, 13, tauDaughter)){
          baby.ntrutausl()++;
        } else baby.ntrutaush()++;
      }
    }

    //Finding lost leptons
    const float relptThreshold(0.3), drThreshold(0.1); 
    float relptThres(2.), drThres(0.15); 
    //if(mcTool->fromWOrWTau(mc) && (laste||lastmu||lasttau|| (mcTool->fromTau(mc) && abs(mc.pdgId())==211))){
    if(eFromTopZtau|| muFromTopZtau || chgPionFromTau){
      double mindr(999.);
      int minind(-1);
      
      for(size_t ind(0); ind < baby.leps_pt().size(); ind++) {
        double dr = dR(baby.leps_phi()[ind],mc.phi(),baby.leps_eta()[ind],mc.eta());
        //      cout<<"mc lep eta phi = "<<mc.eta()<<" "<<mc.phi()<<", reco lep eta phi "<<baby.leps_eta()[ind]<<" "<<baby.leps_phi()[ind]<<", dr = "<<dr<<endl;
        double drelpt(fabs((baby.leps_pt()[ind] - mc.pt())/mc.pt()));
        if(dr > drThres || drelpt > relptThres) continue;
        if(dr < mindr){
          mindr = dr;
          minind = ind;
        }
      }
      
    
      if(minind<0){ //Lepton is lost
        //Try to match to veto track collection
        double min_dr(999.);
        int minindex(-1);
        for(size_t ind(0); ind < baby.tks_pt().size(); ind++) {
          double dr = dR(baby.tks_phi()[ind],mc.phi(),baby.tks_eta()[ind],mc.eta());
          double drelpt(fabs((baby.tks_pt()[ind] - mc.pt())/mc.pt()));
          if(dr > drThreshold || drelpt > relptThres) continue; //use loose pt comparison in case lost lepton is very poorly measured
          if(dr < min_dr){
            min_dr = dr;
            minindex = ind;
          }
        } // Loop over tks
        if(minindex >= 0) {
          baby.tks_tm()[minindex] = true;
        }


      }
    }
    //////// Finding truth-matched leptons
    //const float relptThreshold(0.3), drThreshold(0.1);  needed earlier    
    if(id==11 && fromWOrWTau){
      double mindr(999.);
      int minind(-1);
      for(size_t ind(0); ind < all_els.size(); ind++) {
        double dr(deltaR(mc, *(all_els[ind])));
        double drelpt(fabs((all_els[ind]->pt() - mc.pt())/mc.pt()));
        if(dr > drThreshold || drelpt > relptThreshold) continue;
        if(dr < mindr){
          mindr = dr;
          minind = ind;
        }
      } // Loop over all_els
      if(minind >= 0) {
        baby.els_tm()[minind] = true;
        if(baby.els_sig()[minind]) baby.nleps_tm()++;
      }
      if(lep_tru_pt < mc.pt()){
        lep_tru_pt = mc.pt();
        lep_tru_phi = mc.phi();
      } // Lepton pt to find mt_tru
    } // If it is an electron
    if(id==13 && fromWOrWTau){
      double mindr(999.);
      int minind(-1);
      for(size_t ind(0); ind < all_mus.size(); ind++) {
        double dr(deltaR(mc, *(all_mus[ind])));
        double drelpt(fabs((all_mus[ind]->pt() - mc.pt())/mc.pt()));
        if(dr > drThreshold || drelpt > relptThreshold) continue;
        if(dr < mindr){
          mindr = dr;
          minind = ind;
        }
      } // Loop over all_mus
      if(minind >= 0) {
        baby.mus_tm()[minind] = true;
        if(baby.mus_sig()[minind]) baby.nleps_tm()++;
      }
      if(lep_tru_pt < mc.pt()){
        lep_tru_pt = mc.pt();
        lep_tru_phi = mc.phi();
      } // Lepton pt to find mt_tru
    } // If it is a muon
    if(id==22){
      double mindr(999.);
      int minind(-1);
      for(size_t ind(0); ind < photons.size(); ind++) {
        double dr(deltaR(mc, *(photons[ind])));
        double drelpt(fabs((photons[ind]->pt() - mc.pt())/mc.pt()));
        if(dr > drThreshold || drelpt > relptThreshold) continue;
        if(dr < mindr){
          mindr = dr;
          minind = ind;
        }
      } // Loop over photons
      if(minind >= 0) baby.ph_tm()[minind] = true;
    } // If it is a photon

    //////// Finding true MET
    if((id==12 || id==14 || id==16 || id==18 || id==1000012 || id==1000014 || id==1000016
        || id==1000022 || id==1000023 || id==1000025 || id==1000035 || id==1000039) &&
       id != momid){ // neutrinos decay to themselves
      if(fromWOrWTau) {
        metw_tru_x += mc.px();
        metw_tru_y += mc.py();
      }
    } // If undetected neutral particle

    // don't need to check for gluon splitting if flag is already set
    if(!baby.fromGS()) baby.fromGS()|=mcTool->isFromGSP(dynamic_cast<const reco::Candidate*>(&mc));
  } // Loop over genParticles
    // calculate invariant mass of ttbar pair
  if(topIndex>=0 && antitopIndex>=0) {
    reco::Candidate::LorentzVector topP4 = genParticles->at(topIndex).p4();
    reco::Candidate::LorentzVector antitopP4 = genParticles->at(antitopIndex).p4();
    reco::Candidate::LorentzVector ttbarP4 = topP4+antitopP4;
    baby.m_tt()=ttbarP4.mass();
  }

  baby.ntruleps() = baby.ntrumus()+baby.ntruels()+baby.ntrutaush()+baby.ntrutausl();
  baby.isr_tru_pt() = isr_p4.pt();
  baby.isr_tru_eta() = isr_p4.eta();
  baby.isr_tru_phi() = isr_p4.phi();

  vector<float> isr_sys;
  if(outname.Contains("SMS") || outname.Contains("RPV")){
    isr_sys.push_back(1. + weightTool->isrWeight(baby.isr_tru_pt()));
    isr_sys.push_back(1. - weightTool->isrWeight(baby.isr_tru_pt()));
  }
  else{ isr_sys.push_back(1.); isr_sys.push_back(1.);}
  baby.sys_isr()=isr_sys;

  if((outname.Contains("TTJets") || outname.Contains("TT_")) && top_pt.size() == 2) baby.w_toppt() = weightTool->topPtWeight(top_pt.at(0),top_pt.at(1));
  else baby.w_toppt() = 1.;

  baby.met_tru_nuw() = hypot(metw_tru_x, metw_tru_y);
  baby.met_tru_nuw_phi() = atan2(metw_tru_y, metw_tru_x);

  baby.mt_tru()     = getMT(baby.met_tru(),     baby.met_tru_phi(),     lep_tru_pt, lep_tru_phi);
  baby.mt_tru_nuw() = getMT(baby.met_tru_nuw(), baby.met_tru_nuw_phi(), lep_tru_pt, lep_tru_phi);
} // writeMC

// Finds the jet that minimizes the MET when a variation is performed
void bmaker_full::rebalancedMET( double& minMET, double& minMETPhi)
{
  for(unsigned int iJet=0; iJet<baby.jets_pt().size(); iJet++) {
    // calculate best rescaling factor for this jet
    double rescalingFactor=calculateRescalingFactor(iJet);
    double newMETPhi=0;
    double newMET=calculateRebalancedMET(iJet, rescalingFactor, newMETPhi);
    if(newMET<minMET) {
      minMET=newMET;
      minMETPhi=newMETPhi;
    }
  }
}

// calculate a rebalancing of the jet momentum that minimizes MET
double bmaker_full::calculateRescalingFactor(unsigned int jetIdx)
{
  
  // don't allow jet pt to be scaled by more than this factor
  const double scaleCutoff=1;
  
  TVector3 jet, metVector;
  jet.SetPtEtaPhi(baby.jets_pt().at(jetIdx), baby.jets_eta().at(jetIdx), baby.jets_phi().at(jetIdx));
  metVector.SetPtEtaPhi(baby.met(), 0, baby.met_phi());
  
  double denominator = -jet.Px()*jet.Px()-jet.Py()*jet.Py();
  double numerator = jet.Px()*metVector.Px()+jet.Py()*metVector.Py();
  
  double rescalingFactor=1e6;
  if(denominator!=0) rescalingFactor = numerator/denominator;
  if(fabs(rescalingFactor)>scaleCutoff) rescalingFactor=scaleCutoff*rescalingFactor/fabs(rescalingFactor);
  // the resolution tail is on the _low_ side, not the high side
  // so we always need to subtract pT
  if(rescalingFactor>0) rescalingFactor=0;
  
  return rescalingFactor;
}

double bmaker_full::calculateRebalancedMET(unsigned int jetIdx, double mu, double& METPhi)
{
  TVector3 jet, metVector;
  jet.SetPtEtaPhi(baby.jets_pt().at(jetIdx), baby.jets_eta().at(jetIdx), baby.jets_phi().at(jetIdx));
  metVector.SetPtEtaPhi(baby.met(), 0, baby.met_phi());
 
  double sumPx = metVector.Px()+mu*jet.Px();
  double sumPy = metVector.Py()+mu*jet.Py();

  METPhi=atan(sumPy/sumPx);

  return sqrt(sumPx*sumPx+sumPy*sumPy);
}

void bmaker_full::writeWeights(const vCands &sig_leps, edm::Handle<GenEventInfoProduct> &gen_event_info, 
                               edm::Handle<LHEEventProduct> &lhe_info){
  if (debug) cout<<"INFO: Filling weights..."<<endl;

  // Initializing weights
  if(isData) {
    baby.eff_trig() = baby.w_btag() = baby.w_btag_loose() = baby.w_pu() = baby.w_pu() = baby.w_lep() = baby.w_fs_lep() = baby.w_toppt() = 1.;
    baby.eff_jetid() = baby.w_lumi() = baby.weight() = baby.weight_rpv() = 1.;
    return;
  }

  // Luminosity weight
  const float luminosity = 1000., fneg(xsec::fractionNegWeights(outname));
  baby.w_lumi() = xsec*luminosity / (static_cast<double>(nevents_sample)*(1-2*fneg));
  if (gen_event_info->weight() < 0) baby.w_lumi() *= -1.;
  
  // Pile-up weight
  baby.w_pu() = weightTool->pileupWeight(baby.ntrupv_mean(),0);
  baby.sys_pu().resize(2, 1.);
  baby.sys_pu()[0] = weightTool->pileupWeight(baby.ntrupv_mean(), 1)/baby.w_pu();
  baby.sys_pu()[1] = weightTool->pileupWeight(baby.ntrupv_mean(), -1)/baby.w_pu();

  // Lepton SFs
  pair<double, double> sf_err = lepTool->getScaleFactor(sig_leps);
  double sf  = sf_err.first;
  double unc = sf_err.second/sf_err.first;
  baby.w_lep() = sf;
  baby.sys_lep().resize(2, 1.);
  baby.sys_lep().at(0) = 1.+unc;
  baby.sys_lep().at(1) = 1.-unc; 

  // Lepton SFs in FastSim
  baby.sys_fs_lep().resize(2, 1.);
  if(isFastSim){ 
    pair<double, double> sf_err_fs = lepTool->getScaleFactorFs(sig_leps);
    double sf_fs  = sf_err_fs.first;
    double unc_fs = sf_err_fs.second/sf_err_fs.first;
    baby.w_fs_lep() = sf_fs;
    baby.sys_fs_lep()[0] = 1.+unc_fs;
    baby.sys_fs_lep()[1] = 1.-unc_fs; 
  } else baby.w_fs_lep() = 1.;

  // VVVL trigger efficiency
  baby.eff_trig() = weightTool->triggerEfficiency(baby.nmus(), baby.nels(), baby.sys_trig());
  
  // In FastSim the JetID is broken, so we just apply 0.99 +- 0.01
  if(isFastSim) baby.eff_jetid() = 0.99;
  else baby.eff_jetid() = 1.;
  ////////////  Total weight  ////////////
  // w_btag calculated in writeJets
  // w_toppt and sys_isr calculated in writeMC
  baby.weight() = baby.w_lumi() * baby.w_lep() * baby.w_fs_lep() * baby.w_toppt() * baby.w_btag() 
    * baby.eff_jetid();
  baby.weight_rpv() = baby.w_lumi() * baby.w_lep() * baby.w_fs_lep() * baby.w_toppt() * baby.w_btag()
    * baby.w_pu() * baby.eff_jetid();

  /////// Systematics that do not change central value /////////
  if(lhe_info.isValid()) weightTool->getTheoryWeights(lhe_info);
  // Renormalization and Factorization scales
  baby.sys_mur().push_back(weightTool->theoryWeight(weight_tools::muRup));
  baby.sys_mur().push_back(weightTool->theoryWeight(weight_tools::muRdown));
  baby.sys_muf().push_back(weightTool->theoryWeight(weight_tools::muFup));
  baby.sys_muf().push_back(weightTool->theoryWeight(weight_tools::muFdown));
  baby.sys_murf().push_back(weightTool->theoryWeight(weight_tools::muRup_muFup));
  baby.sys_murf().push_back(weightTool->theoryWeight(weight_tools::muRdown_muFdown));
  // PDF variations
  weightTool->getPDFWeights(baby.sys_pdf(), baby.w_pdf());

}

/*
  _____                 _                   _                 
  /  __ \               | |                 | |                
  | /  \/ ___  _ __  ___| |_ _ __ _   _  ___| |_ ___  _ __ ___ 
  | |    / _ \| '_ \/ __| __| '__| | | |/ __| __/ _ \| '__/ __|
  | \__/\ (_) | | | \__ \ |_| |  | |_| | (__| || (_) | |  \__ \
  \____/\___/|_| |_|___/\__|_|   \__,_|\___|\__\___/|_|  |___/
*/
//Constructor (this line searchable)
bmaker_full::bmaker_full(const edm::ParameterSet& iConfig):
  outname(TString(iConfig.getParameter<string>("outputFile"))),
  inputfiles(iConfig.getParameter<vector<string> >("inputFiles")),
  jsonfile(iConfig.getParameter<string>("json")),
  condor_subtime(iConfig.getParameter<string>("condor_subtime")),
  jec_label(iConfig.getParameter<string>("jec")),
  met_label(iConfig.getParameter<edm::InputTag>("met")),
  met_nohf_label(iConfig.getParameter<edm::InputTag>("met_nohf")),
  jets_label(iConfig.getParameter<edm::InputTag>("jets")),
  nevents_sample(iConfig.getParameter<unsigned int>("nEventsSample")),
  nevents(0),
  doMetRebalancing(iConfig.getParameter<bool>("doMetRebalancing")),
  addBTagWeights(iConfig.getParameter<bool>("addBTagWeights")),
  isFastSim(iConfig.getParameter<bool>("isFastSim")),
  doSystematics(iConfig.getParameter<bool>("doSystematics")),
  debug(iConfig.getParameter<bool>("debugMode")),

  // Initialize tokens
  tok_trigResults_hlt_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))),
  tok_patTrig_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),
  tok_primVertex_(consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"))),
  tok_addPileup_(consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("addPileupInfo"))),
  tok_slimAddPileup_(consumes<std::vector< PileupSummaryInfo > >(edm::InputTag("slimmedAddPileupInfo"))),
  tok_packPFCands_(consumes<pat::PackedCandidateCollection>(edm::InputTag("packedPFCandidates"))),
  tok_rhoFastJet_centralNeutral_(consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"))),
  tok_muons_(consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"))),
  tok_electrons_(consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"))),
  tok_rhoFastJet_all_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
  tok_offBeamSpot_(consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"))),
  tok_photons_(consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"))),
  tok_reducedEgamma_conver_(consumes<vector<reco::Conversion> >(edm::InputTag("reducedEgamma","reducedConversions"))),
  tok_jets_(consumes<pat::JetCollection>(jets_label)),
  tok_genJets_(consumes<edm::View<reco::GenJet> >(edm::InputTag("slimmedGenJets"))),
  tok_met_(consumes<pat::METCollection>(met_label)),
  tok_met_noHF_(consumes<pat::METCollection>(met_nohf_label)),
  tok_HBHENoiseFilter_(consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"))),
  tok_HBHEIsoNoiseFilter_(consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"))),
  tok_trigResults_reco_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","RECO"))),
  tok_trigResults_pat_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","PAT"))),
  tok_selectedPatTrig_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("selectedPatTrigger"))),
  tok_pruneGenPart_(consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"))),
  tok_extLHEProducer_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"))),
  tok_source_(consumes<LHEEventProduct>(edm::InputTag("source"))),
  tok_generator_(consumes<GenEventInfoProduct>(edm::InputTag("generator")))
  #ifdef POST_7_4
  ,
    tok_genlumiheader_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator")))
  #endif
{
  time(&startTime);


  lepTool    = new lepton_tools();
  jetTool    = new jet_met_tools(jec_label, doSystematics, isFastSim);
  photonTool = new photon_tools();
  mcTool     = new mc_tools();
  weightTool = new weight_tools();
  eventTool  = new event_tools(outname);
  

  outfile = new TFile(outname, "recreate");
  outfile->cd();
  baby.tree_.SetDirectory(outfile);

  xsec = xsec::crossSection(outname);
  // if(xsec<=0) {
  //   cout<<"BABYMAKER: Cross section not found, aborting"<<endl<<endl;
  //   exit(1);
  // }

  trig_name = vector<TString>();
  if(outname.Contains("Run201")){ // Would like to define isData, but need iEvent?
    if(outname.Contains("Run2016")){
      trig_name.push_back("HLT_PFHT300_PFMET100_v");                            // 0 
      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v");                // 1 
      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");                        // 2
      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT400_v");                        // 3
      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_v");                        // 4 
      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v");               // 5 
      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");                       // 6
      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT400_v");                       // 7
      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_v");                       // 8 
      trig_name.push_back("HLT_DoubleMu8_Mass8_PFHT300_v");                     // 9
      trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v");   // 10

      trig_name.push_back("HLT_PFHT475_v");                                     // 11
      trig_name.push_back("HLT_PFHT800_v");                                     // 12
      trig_name.push_back("HLT_PFMET100_PFMHT100_IDTight_v");                   // 13
      trig_name.push_back("HLT_PFMET110_PFMHT110_IDTight_v");                   // 14
      trig_name.push_back("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v");           // 15
      trig_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");		// 16
      trig_name.push_back("HLT_Mu45_eta2p1_v");                                 // 17
      trig_name.push_back("HLT_IsoMu18_v");                                     // 18
      trig_name.push_back("HLT_IsoMu24_v");					// 19
      trig_name.push_back("HLT_IsoMu27_v");                                     // 20

      trig_name.push_back("HLT_Mu50_v");                                        // 21
      trig_name.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");                    // 22
      trig_name.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v");                    // 23
      trig_name.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");                   // 24
      trig_name.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");       // 25
      trig_name.push_back("HLT_Photon175_v");					// 26
      trig_name.push_back("HLT_Photon90_CaloIdL_PFHT500_v");                    // 27
      trig_name.push_back("HLT_PFMET90_PFMHT90_IDTight_v");			// 28
      trig_name.push_back("HLT_Ele23_WPLoose_Gsf_v");			        // 29
      trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_v");                   // 30

      trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");           // 31
      trig_name.push_back("HLT_IsoMu22_v");					// 32
      trig_name.push_back("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v");           // 33
      trig_name.push_back("HLT_Mu50_IsoVVVL_PFHT400_v");                        // 34
      trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_v");           // 35
      trig_name.push_back("HLT_Ele50_IsoVVVL_PFHT400_v");                       // 36
      trig_name.push_back("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_v");          // 37
    } else {
      trig_name.push_back("HLT_PFHT350_PFMET100_");                               // 0 
      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v");                  // 1 
      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");                          // 2
      trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_v");              // 3
      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_v");                          // 4 
      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v");                 // 5 
      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");                         // 6
      trig_name.push_back("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_v");             // 7
      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_v");                         // 8 
      trig_name.push_back("HLT_DoubleMu8_Mass8_PFHT300_v");                       // 9
      trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v");     // 10
      trig_name.push_back("HLT_PFHT475_v");                                       // 11
      trig_name.push_back("HLT_PFHT800_v");                                       // 12
      trig_name.push_back("HLT_PFMET120_JetIdCleaned_Mu5_v");                     // 13
      trig_name.push_back("HLT_PFMET170_JetIdCleaned_v");                         // 14
      trig_name.push_back("HLT_DoubleIsoMu17_eta2p1_v");                          // 15
      trig_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");               // 16
      trig_name.push_back("HLT_IsoMu20_v");                                       // 17
      trig_name.push_back("HLT_IsoMu18_v");                                       // 18
      trig_name.push_back("HLT_IsoMu24_eta2p1_v");                                // 19
      trig_name.push_back("HLT_IsoMu27_v");                                       // 20
      trig_name.push_back("HLT_Mu50_v");                                          // 21
      trig_name.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");                      // 22
      trig_name.push_back("HLT_Ele23_WPLoose_Gsf_v");                             // 23
      trig_name.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");                     // 24
      trig_name.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v");                // 25
      trig_name.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v");             // 26
      trig_name.push_back("HLT_Photon90_CaloIdL_PFHT500_v");                      // 27
      trig_name.push_back("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v");  // 28
    }
  } else {
    trig_name.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_v");                 // 0 
    trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT400_PFMET70_v");                  // 1 
    trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");                          // 2
    trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV07_PFHT400_v");                // 3
    trig_name.push_back("HLT_Mu15_PFHT300_v");                                  // 4 
    trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT400_PFMET70_v");                 // 5 
    trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");                         // 6
    trig_name.push_back("HLT_Ele15_IsoVVVL_BTagtop8CSV07_PFHT400_v");           // 7
    trig_name.push_back("HLT_Ele15_PFHT300_v");                                 // 8 
    trig_name.push_back("HLT_DoubleMu8_Mass8_PFHT300_v");                       // 9
    trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v");     // 10
    trig_name.push_back("HLT_PFHT350_v");                                       // 11
    trig_name.push_back("HLT_PFHT900_v");                                       // 12
    trig_name.push_back("HLT_PFMET120_NoiseCleaned_Mu5_v");                     // 13
    trig_name.push_back("HLT_PFMET170_NoiseCleaned_v");                         // 14
    trig_name.push_back("HLT_DoubleIsoMu17_eta2p1_v");                          // 15
    trig_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");               // 16
    trig_name.push_back("HLT_IsoMu20_v");                                       // 17
    trig_name.push_back("HLT_IsoMu17_eta2p1_v");                                // 18
    trig_name.push_back("HLT_IsoMu24_eta2p1_v");                                // 19
    trig_name.push_back("HLT_IsoMu27_v");                                       // 20
    trig_name.push_back("HLT_Mu50_v");                                          // 21
    trig_name.push_back("HLT_Ele27_eta2p1_WP75_Gsf_v");                         // 22
    trig_name.push_back("HLT_Ele22_eta2p1_WP75_Gsf_v");                         // 23
    trig_name.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");                     // 24
    trig_name.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v");                // 25
    trig_name.push_back("HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v");                // 26
    trig_name.push_back("HLT_Photon90_CaloIdL_PFHT500_v");                      // 27
    trig_name.push_back("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v");  // 28
  }

}


bmaker_full::~bmaker_full(){
  outfile->cd();
  baby.tree_.SetDirectory(outfile);
  baby.Write();

  string commit_s = execute("git rev-parse HEAD");
  while(!commit_s.empty() && commit_s.at(commit_s.length()-1) == '\n') commit_s.erase(commit_s.length()-1);
  TString commit = commit_s;
  TString type = baby.Type();
  TString root_version = gROOT->GetVersion();
  TString root_tutorial_dir = gROOT->GetTutorialsDir();
  TString user(getenv("ORIGIN_USER"));
  if (user=="") user = getenv("USER");
  TString cmssw(getenv("CMSSW_BASE"));
  time_t curTime;
  time(&curTime);
  char time_c[100];
  struct tm * timeinfo = localtime(&curTime);
  strftime(time_c,100,"%Y-%m-%d %H:%M:%S",timeinfo);
  TString date(time_c);
  int seconds(floor(difftime(curTime,startTime)+0.5));

  vector<TString> sys_names;
  sys_names.resize(kSysLast,"");
  sys_names[kSysJER] = "jer";
  sys_names[kSysJECUp] = "jec_up";
  sys_names[kSysJECDn] = "jec_dn";

  TTree treeglobal("treeglobal", "treeglobal");
  treeglobal.Branch("nev_sample", &nevents_sample);
  treeglobal.Branch("nev_file", &nevents);
  treeglobal.Branch("runtime_seconds", &seconds);
  treeglobal.Branch("git_commit", &commit);
  // treeglobal.Branch("model", &model);
  treeglobal.Branch("baby_type", &type);
  treeglobal.Branch("root_version", &root_version);
  treeglobal.Branch("root_tutorial_dir", &root_tutorial_dir);
  treeglobal.Branch("trig_names", &trig_name);
  treeglobal.Branch("sys_names", &sys_names);
  treeglobal.Branch("xsec", &xsec);
  treeglobal.Branch("user", &user);
  treeglobal.Branch("cmssw", &cmssw);
  treeglobal.Branch("jec", &jec_label);
  treeglobal.Branch("json", &jsonfile);
  treeglobal.Branch("date", &date);
  treeglobal.Branch("inputfiles", &inputfiles);
  treeglobal.Branch("condor_subtime", &condor_subtime);
  treeglobal.Fill();
  treeglobal.SetDirectory(outfile);
  treeglobal.Write();
  
  outfile->Close();

  int minutes((seconds/60)%60), hours(seconds/3600);
  TString runtime("");
  if(hours<10) runtime += "0";
  runtime += hours; runtime += ":";
  if(minutes<10) runtime += "0";
  runtime += minutes; runtime += ":";
  if((seconds%60)<10) runtime += "0";
  runtime += seconds%60; 
  float hertz(nevents); hertz /= seconds;
  cout<<endl<<"BABYMAKER: Written "<<nevents<<" events in "<<outname<<". It took "<<seconds<<" seconds to run ("<<runtime<<"), "
      <<roundNumber(hertz,1)<<" Hz, "<<roundNumber(1000,2,hertz)<<" ms per event"<<endl<<endl;
  cout<<"BABYMAKER: *********** List of input files ***********"<<endl;
  for(size_t ifile(0); ifile < inputfiles.size(); ifile++)
    cout<<"BABYMAKER: "<<inputfiles[ifile].c_str()<<endl;
  cout<<endl;

  delete outfile;

  delete lepTool;
  delete photonTool;
  delete jetTool;
  delete mcTool;
  delete weightTool;
}

void bmaker_full::reportTime(const edm::Event& iEvent){
  // Time reporting
  if(nevents==1) {
    time_t curTime;
    time(&curTime);
    cout<<endl<<"BABYMAKER: Took "<<roundNumber(difftime(curTime,startTime),1)<<" seconds for set up and run first event"
        <<endl<<endl;
    time(&startTime);
  }
  if(debug || (nevents<100&&nevents%10==0) || (nevents<1000&&nevents%100==0) 
     || (nevents<10000&&nevents%1000==0) || nevents%10000==0) {
    time_t curTime;
    time(&curTime);
    float seconds(difftime(curTime,startTime));
    cout<<"BABYMAKER: Run "<<iEvent.id().run()<<", Event "<< setw(8)<<iEvent.id().event()
        <<", LumiSection "<< setw(5)<< iEvent.luminosityBlock()
        <<". Ran "<<setw(7)<<nevents<<" events in "<<setw(7)<<seconds<<" seconds -> "
        <<setw(5)<<roundNumber(nevents-1,1,seconds)<<" Hz, "
        <<setw(5)<<roundNumber(seconds*1000,2,nevents-1)<<" ms per event"<<endl;
  }
}

// ------------ method called once each job just before starting event loop  ------------
void bmaker_full::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void bmaker_full::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  bmaker_full::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  bmaker_full::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------

void bmaker_full::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iEvent){
#ifdef POST_7_4
  if (outname.Contains("PUSpring16Fast") && outname.Contains("SMS-")){
    edm::Handle<GenLumiInfoHeader> gen_header;
    iLumi.getByToken(tok_genlumiheader_, gen_header);  
    string model = gen_header->configDescription();
    mcTool->getMassPoints(model, mprod_, mlsp_);
  }
#endif
}


// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  bmaker_full::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bmaker_full::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bmaker_full);
