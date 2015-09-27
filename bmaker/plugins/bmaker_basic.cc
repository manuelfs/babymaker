//// BMAKER_BASIC: Creates baby tree with basic branches
//// Function names follow the first-lowercase, following words-uppercase. No underscores


// System include files
#include <cmath>
#include <memory>
#include <iostream>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// FW physics include files
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>

// ROOT include files
#include "TFile.h"
#include "TROOT.h"

// User include files
#include "babymaker/bmaker/interface/bmaker_basic.hh"
#include "babymaker/bmaker/interface/baby_basic.hh"
#include "babymaker/bmaker/interface/phys_objects.hh"

using namespace std;
using namespace phys_objects;
using namespace utilities;

///////////////////////// analyze: METHOD CALLED EACH EVENT ///////////////////////////
void bmaker_basic::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  isData = iEvent.isRealData();
  nevents++;

  ////////////////////// Event info /////////////////////
  baby.run() = iEvent.id().run();
  baby.event() = iEvent.id().event();
  baby.lumiblock() = iEvent.luminosityBlock();
  if(isData){
    if(!isInJSON("dcs", baby.run(), baby.lumiblock())) return;
    baby.json() = isInJSON("golden", baby.run(), baby.lumiblock());
  } else baby.json() = true;

  const float luminosity = 1000.;

  
  baby.weight() = 1.;
  if(!isData) {
    edm::Handle<GenEventInfoProduct> gen_event_info;
    iEvent.getByLabel("generator", gen_event_info);
    if (gen_event_info->weight() < 0) baby.weight() *= -1.;
    baby.weight() *= xsec*luminosity / static_cast<double>(nevents_sample);

  }

  //////////////////////////// MET ///////////////////////////
  edm::Handle<pat::METCollection> mets;
  iEvent.getByLabel(met_label, mets);
  baby.met() = mets->at(0).pt();
  baby.met_phi() = mets->at(0).phi();

  ////////////////////// Primary vertices /////////////////////
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtx);
  edm::Handle<std::vector< PileupSummaryInfo > >  pu_info;
  if(!isData) iEvent.getByLabel("addPileupInfo", pu_info);
  writeVertices(vtx, pu_info);

  ////////////////////// Trigger /////////////////////
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"),triggerBits);  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByLabel("patTrigger",triggerPrescales);  
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  writeTriggers(names, triggerBits, triggerPrescales);

  //////////////// HLT objects //////////////////
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByLabel("selectedPatTrigger",triggerObjects);  
  writeHLTObjects(names, triggerObjects);

  //////////////////// pfcands  //////////////////////
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByLabel("packedPFCandidates", pfcands);

  ////////////////////// Leptons /////////////////////
  vCands leptons, muons, electrons;
  edm::Handle<pat::MuonCollection> allmuons;
  iEvent.getByLabel("slimmedMuons", allmuons);
  muons = writeMuons(allmuons, pfcands, vtx);
  edm::Handle<pat::ElectronCollection> allelectrons;
  iEvent.getByLabel("slimmedElectrons", allelectrons);
  electrons = writeElectrons(allelectrons, pfcands, vtx);

  // Putting muons and electrons together
  leptons = muons;
  leptons.insert(leptons.end(), electrons.begin(), electrons.end());
  writeLeptons(leptons);

  ////////////////////// mT, dphi /////////////////////
  if(leptons.size()>0){
    float wx = baby.met()*cos(baby.met_phi()) + leptons[0]->px();
    float wy = baby.met()*sin(baby.met_phi()) + leptons[0]->py();
    float wphi = atan2(wy, wx);

    baby.dphi_wlep() = deltaPhi(wphi, leptons[0]->phi());
    baby.mt() = sqrt(2.*baby.met()*leptons[0]->pt()*(1.-cos(leptons[0]->phi()-baby.met_phi())));
  }   
    

  ////////////////////// Jets /////////////////////
  vCands jets;
  edm::Handle<pat::JetCollection> alljets;
  iEvent.getByLabel(jets_label, alljets);
  jets = writeJets(alljets, leptons);
  writeFatJets(jets);

  ///////////////////// Filters ///////////////////////
  edm::Handle<edm::TriggerResults> filterBits;
  std::string processLabel = isData ? "RECO":"PAT"; // prompt reco runs in the "RECO" process
  iEvent.getByLabel(edm::InputTag("TriggerResults","",processLabel),filterBits);  
  // re-recoed data will have the process label "PAT" rather than "RECO";
  // if the attempt to find data with "RECO" process fails, try "PAT"
  if(!filterBits.isValid() && isData) 
      iEvent.getByLabel(edm::InputTag("TriggerResults", "", "PAT"),filterBits);  
  const edm::TriggerNames &fnames = iEvent.triggerNames(*filterBits);
  // the HBHE noise filter needs to be recomputed in early 2015 data
  edm::Handle<bool> filter_hbhe;
  if(isData && iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult",filter_hbhe)) { 
    if(*filter_hbhe) baby.pass_hbhe() = true;
    else baby.pass_hbhe() = false;
  }
  writeFilters(fnames, filterBits, vtx);

  ///////////////////// Truth Info ///////////////////////
  if (!isData) {
    edm::Handle<LHEEventProduct> lhe_info;
    iEvent.getByLabel( "externalLHEProducer", lhe_info);
    writeGenInfo(lhe_info);
  }

  ////////////////// Filling the tree //////////////////
  baby.Fill();

  // Time reporting
  if(nevents==1) {
    time_t curTime;
    time(&curTime);
    cout<<endl<<"Took "<<roundNumber(difftime(curTime,startTime),1)<<" seconds for set up and run first event"<<endl<<endl;
    time(&startTime);
  }
  if(nevents<10 || (nevents<100&&nevents%10==0) || (nevents<1000&&nevents%100==0) 
     || (nevents<10000&&nevents%1000==0) || nevents%10000==0) {
    time_t curTime;
    time(&curTime);
    float seconds(difftime(curTime,startTime));
    cout<<"Run "<<iEvent.id().run()<<", Events "<< iEvent.id().event()<<", LumiSection "<< iEvent.luminosityBlock()
	<<". Ran "<<nevents<<" events in "<<seconds<<" seconds -> "<<roundNumber(nevents-1,1,seconds)<<" Hz and "
	<<roundNumber(seconds*1000,2,nevents-1)<<" ms"<<endl;
  }
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
vCands bmaker_basic::writeJets(edm::Handle<pat::JetCollection> alljets, vCands &leptons){
  vCands jets;
  baby.njets() = 0; baby.nbl() = 0; baby.nbm() = 0;  baby.nbt() = 0;  
  baby.ht() = 0.;
  baby.pass_jets() = true; 
  bool goodPtEta, isLep, goodID, goodJet;
  for (size_t ijet(0); ijet < alljets->size(); ijet++) {
    const pat::Jet &jet = (*alljets)[ijet];

    // Saving good jets and jets corresponding to signal leptons
    isLep = leptonInJet(jet, leptons);
    goodID = idJet(jet);
    goodPtEta = jet.pt() > JetPtCut && fabs(jet.eta()) <= JetEtaCut;
    goodJet = (!isLep) && goodID && goodPtEta;
    if(goodPtEta && !isLep && !goodID) baby.pass_jets_nohf() = false;
    if(jet.pt() > JetPtCut && !isLep && !goodID) baby.pass_jets() = false;
    if(!isLep && !goodJet) continue;

    // Filling branches
    baby.jets_pt().push_back(jet.pt());
    baby.jets_eta().push_back(jet.eta());
    baby.jets_phi().push_back(jet.phi());
    baby.jets_m().push_back(jet.mass());
    baby.jets_islep().push_back(isLep);
    
    float csv(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    baby.jets_csv().push_back(csv);

    jets.push_back(dynamic_cast<const reco::Candidate *>(&jet));

    if(goodJet){
      baby.njets()++;
      baby.ht() += jet.pt();
      if(csv>CSVLoose)  baby.nbl()++;
      if(csv>CSVMedium) baby.nbm()++;
      if(csv>CSVTight)  baby.nbt()++;
    }
  } // Loop over jets  

  return jets; // Returning jets that pass acceptance and ID, regardless of whether they're leptons
} 

void bmaker_basic::writeFatJets(vCands &jets){
  clusterFatJets(baby.nfjets(), baby.mj(),
		 baby.fjets_pt(), baby.fjets_eta(),
		 baby.fjets_phi(), baby.fjets_m(),
		 baby.fjets_nconst(),
		 baby.fjets_sumcsv(), baby.fjets_poscsv(),
		 baby.fjets_btags(), baby.jets_fjet_index(),
		 1.2, jets);
  clusterFatJets(baby.nfjets08(), baby.mj08(),
                 baby.fjets08_pt(), baby.fjets08_eta(),
                 baby.fjets08_phi(), baby.fjets08_m(),
                 baby.fjets08_nconst(),
                 baby.fjets08_sumcsv(), baby.fjets08_poscsv(),
                 baby.fjets08_btags(), baby.jets_fjet08_index(),
                 0.8, jets);

}
void bmaker_basic::clusterFatJets(int &nfjets, float &mj,
				  vector<float> &fjets_pt, 
				  vector<float> &fjets_eta,
				  vector<float> &fjets_phi, 
				  vector<float> &fjets_m,
				  vector<int> &fjets_nconst,
				  vector<float> &fjets_sumcsv,
				  vector<float> &fjets_poscsv,
				  vector<int> &fjets_btags,
				  vector<int> &jets_fjet_index,
				  double radius,
				  vCands &jets){

  vector<fastjet::PseudoJet> sjets(0);
  vector<float> csvs(0);
  for(size_t ijet(0); ijet < jets.size(); ijet++){
    const fastjet::PseudoJet this_pj(jets[ijet]->px(), jets[ijet]->py(), jets[ijet]->pz(), jets[ijet]->energy());
    sjets.push_back(this_pj);
    csvs.push_back(baby.jets_csv()[ijet]); // This require to have the same jets in baby and jets
  } // Loop over skinny jets
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, radius);
  fastjet::ClusterSequence cs(sjets, jet_def);
  vector<fastjet::PseudoJet> fjets = cs.inclusive_jets();
  sort(fjets.begin(), fjets.end(), greaterM);
  nfjets = 0;
  mj = 0.;
  fjets_pt.resize(fjets.size());
  fjets_eta.resize(fjets.size());
  fjets_phi.resize(fjets.size());
  fjets_m.resize(fjets.size());
  fjets_nconst.resize(fjets.size());
  fjets_sumcsv.resize(fjets.size());
  fjets_poscsv.resize(fjets.size());
  fjets_btags.resize(fjets.size());
  jets_fjet_index.resize(jets.size());

  for(size_t ipj = 0; ipj < fjets.size(); ++ipj){
    const fastjet::PseudoJet &pj = fjets.at(ipj);
    fjets_pt.at(ipj) = pj.pt();
    fjets_eta.at(ipj) = pj.eta();
    fjets_phi.at(ipj) = pj.phi_std();
    fjets_m.at(ipj) = pj.m();
    const vector<fastjet::PseudoJet> &cjets = pj.constituents();
    fjets_nconst.at(ipj) = cjets.size();
    mj += pj.m();
    ++nfjets;
    fjets_btags.at(ipj) = 0;
    fjets_sumcsv.at(ipj) = 0.;
    fjets_poscsv.at(ipj) = 0.;
    for(size_t ijet = 0; ijet < jets.size(); ++ijet){
      for(size_t cjet = 0; cjet < cjets.size(); ++ cjet){
        if((cjets.at(cjet) - sjets.at(ijet)).pt() < 0.0001){
          jets_fjet_index.at(ijet) = ipj;
          fjets_sumcsv.at(ipj) += csvs.at(ijet);
          if(csvs.at(ijet) > 0.){
            fjets_poscsv.at(ipj) += csvs.at(ijet);
          }
          if(csvs.at(ijet) > CSVMedium){
            ++(fjets_btags.at(ipj));
          }
        }
      } // Loop over fat jet constituents
    } // Loop over skinny jets
  } // Loop over fat jets

}


vCands bmaker_basic::writeMuons(edm::Handle<pat::MuonCollection> muons, 
				edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<reco::VertexCollection> vtx){
  vCands signalmuons; 
  baby.nmus() = 0; baby.nvmus() = 0;
  for (size_t ilep(0); ilep < muons->size(); ilep++) {
    const pat::Muon &lep = (*muons)[ilep];    
    if(!isVetoMuon(lep, vtx, -99.)) continue; // Storing leptons that pass all veto cuts except for iso

    double lep_iso(getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., false));
    double dz(0.), d0(0.);
    vertexMuon(lep, vtx, dz, d0); // Calculating dz and d0

    baby.mus_pt().push_back(lep.pt());
    baby.mus_eta().push_back(lep.eta());
    baby.mus_phi().push_back(lep.phi());
    baby.mus_dz().push_back(dz);
    baby.mus_d0().push_back(d0);
    baby.mus_charge().push_back(lep.charge());
    baby.mus_medium().push_back(idMuon(lep, vtx, kMedium));
    baby.mus_tight().push_back(idMuon(lep, vtx, kTight));
    baby.mus_miniso().push_back(lep_iso);

    if(isVetoMuon(lep, vtx, lep_iso))   baby.nvmus()++;
    if(isSignalMuon(lep, vtx, lep_iso)) {
      baby.nmus()++;
      signalmuons.push_back(dynamic_cast<const reco::Candidate *>(&lep));
    }
  } // Loop over muons
  
  return signalmuons;
}


vCands bmaker_basic::writeElectrons(edm::Handle<pat::ElectronCollection> electrons, 
				    edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<reco::VertexCollection> vtx){
  vCands signalelectrons; 
  baby.nels() = 0; baby.nvels() = 0;
  for (size_t ilep(0); ilep < electrons->size(); ilep++) {
    const pat::Electron &lep = (*electrons)[ilep];    
    if(!isVetoElectron(lep, vtx, -99.)) continue; // Storing leptons that pass all veto cuts except for iso

    double lep_iso(getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., false));
    double dz(0.), d0(0.);
    vertexElectron(lep, vtx, dz, d0); // Calculating dz and d0

    baby.els_pt().push_back(lep.pt());
    baby.els_sceta().push_back(lep.superCluster()->position().eta());
    baby.els_eta().push_back(lep.eta());
    baby.els_phi().push_back(lep.phi());
    baby.els_dz().push_back(dz);
    baby.els_d0().push_back(d0);
    baby.els_charge().push_back(lep.charge());
    baby.els_medium().push_back(idElectron(lep, vtx, kMedium));
    baby.els_ispf().push_back(lep.numberOfSourceCandidatePtrs()==2 && abs(lep.sourceCandidatePtr(1)->pdgId())==11);
    baby.els_tight().push_back(idElectron(lep, vtx, kTight));
    baby.els_miniso().push_back(lep_iso);

    if(isVetoElectron(lep, vtx, lep_iso))   baby.nvels()++;
    if(isSignalElectron(lep, vtx, lep_iso)) {
      baby.nels()++;
      signalelectrons.push_back(dynamic_cast<const reco::Candidate *>(&lep));
    }
  } // Loop over electrons

  return signalelectrons;
}

void bmaker_basic::writeLeptons(vCands &leptons){ 
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


void bmaker_basic::writeTriggers(const edm::TriggerNames &names, 
                                 edm::Handle<edm::TriggerResults> triggerBits, 
                                 edm::Handle<pat::PackedTriggerPrescales> triggerPrescales){
  baby.trig().resize(trig_name.size(), false);
  baby.trig_prescale().resize(trig_name.size(), -1.);
  for (size_t itrig(0); itrig < triggerBits->size(); itrig++) {
    for(size_t itn(0); itn < trig_name.size(); itn++){
      if(names.triggerName(itrig).find(trig_name[itn])!=string::npos){
        baby.trig()[itn] = triggerBits->accept(itrig);
        baby.trig_prescale()[itn] = triggerPrescales->getPrescaleForIndex(itrig);
      }
    }
  }
}

void bmaker_basic::writeHLTObjects(const edm::TriggerNames &names, 
                                   edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects){
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackPathNames(names);
    TString name(obj.collection());
    float objpt(obj.pt());
    if(name=="hltPFMETProducer::HLT") baby.onmet() = objpt;
    else if(name=="hltPFHT::HLT") baby.onht() = objpt; // There's 2 of these, and we want the last one
    else if(name=="hltL3MuonCandidates::HLT" && baby.onmaxmu()<objpt) baby.onmaxmu() = objpt;
    else if(name=="hltEgammaCandidates::HLT" && baby.onmaxel()<objpt) baby.onmaxel() = objpt;
  }
}

void bmaker_basic::writeFilters(const edm::TriggerNames &fnames,
                                  edm::Handle<edm::TriggerResults> filterBits,
                                  edm::Handle<reco::VertexCollection> vtx){
  for (size_t i(0); i < filterBits->size(); ++i) {
    string name = fnames.triggerName(i);
    bool pass = static_cast<bool>(filterBits->accept(i));
    if (name=="Flag_goodVertices") baby.pass_goodv() = pass;
    else if (name=="Flag_CSCTightHaloFilter") baby.pass_cschalo() = pass;
    else if (name=="Flag_eeBadScFilter") baby.pass_eebadsc() = pass;
    //else if (name=="Flag_HBHENoiseFilter") baby.pass_hbhe() = pass;
  }

  baby.pass_goodv() &= hasGoodPV(vtx);

  baby.pass() = baby.pass_hbhe() && baby.pass_goodv() && baby.pass_cschalo() && baby.pass_eebadsc() && baby.pass_jets();
  baby.pass_nohf() = baby.pass_hbhe() && baby.pass_goodv() && baby.pass_cschalo() && baby.pass_eebadsc() && baby.pass_jets_nohf();
}

void bmaker_basic::writeVertices(edm::Handle<reco::VertexCollection> vtx,
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

void bmaker_basic::writeGenInfo(edm::Handle<LHEEventProduct> lhe_info){
  baby.nisr_me()=0; baby.ht_isr_me()=0.; baby.ntruleps()=0;
  for ( unsigned int icount = 0 ; icount < (unsigned int)lhe_info->hepeup().NUP; icount++ ) {
    unsigned int pdgid = abs(lhe_info->hepeup().IDUP[icount]);
    int status = lhe_info->hepeup().ISTUP[icount];
    int mom1id = abs(lhe_info->hepeup().IDUP[lhe_info->hepeup().MOTHUP[icount].first-1]);
    int mom2id = abs(lhe_info->hepeup().IDUP[lhe_info->hepeup().MOTHUP[icount].second-1]);
    float px = (lhe_info->hepeup().PUP[icount])[0];
    float py = (lhe_info->hepeup().PUP[icount])[1];
    float pt = sqrt(px*px+py*py);

    if(status==1 && (pdgid<6 || pdgid==21) && mom1id!=6 && mom1id!=24 && mom2id!=6 && mom2id!=24) {
       baby.nisr_me()++;
       baby.ht_isr_me() += pt;
    }

    if (status==1 && (pdgid==11 || pdgid==13 || pdgid==15)) baby.ntruleps()++;
  } // Loop over generator particles
}

/*
 _____                 _                   _                 
/  __ \               | |                 | |                
| /  \/ ___  _ __  ___| |_ _ __ _   _  ___| |_ ___  _ __ ___ 
| |    / _ \| '_ \/ __| __| '__| | | |/ __| __/ _ \| '__/ __|
| \__/\ (_) | | | \__ \ |_| |  | |_| | (__| || (_) | |  \__ \
 \____/\___/|_| |_|___/\__|_|   \__,_|\___|\__\___/|_|  |___/
*/

bmaker_basic::bmaker_basic(const edm::ParameterSet& iConfig):
  outname(TString(iConfig.getParameter<string>("outputFile"))),
  met_label(iConfig.getParameter<edm::InputTag>("met")),
  jets_label(iConfig.getParameter<edm::InputTag>("jets")),
  nevents_sample(iConfig.getParameter<unsigned int>("nEventsSample")),
  nevents(0){

  time(&startTime);

  outfile = new TFile(outname, "recreate");
  outfile->cd();
  baby.tree_.SetDirectory(outfile);

  xsec = crossSection(outname);

  trig_name = vector<TString>();
  trig_name.push_back("HLT_PFHT350_PFMET100_JetIdCleaned_v");			// 0 
  trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v");			// 1 
  trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");				// 2
  trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_v");		// 3
  trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_v");				// 4 
  trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v");			// 5 
  trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");				// 6
  trig_name.push_back("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_v");		// 7
  trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_v");				// 8 
  trig_name.push_back("HLT_DoubleMu8_Mass8_PFHT300_v");				// 9
  trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v");	// 10
  trig_name.push_back("HLT_PFHT475_v");						// 11
  trig_name.push_back("HLT_PFHT800_v");						// 12
  trig_name.push_back("HLT_PFMET120_NoiseCleaned_Mu5_v");			// 13
  trig_name.push_back("HLT_PFMET170_NoiseCleaned_v");				// 14
  trig_name.push_back("HLT_DoubleIsoMu17_eta2p1_v");				// 15
  trig_name.push_back("HLT_Mu17_TrkIsoVVL_v");					// 16
  trig_name.push_back("HLT_IsoMu17_eta2p1_v");					// 17
  trig_name.push_back("HLT_IsoMu20_v");						// 18
  trig_name.push_back("HLT_IsoMu24_eta2p1_v");					// 19
  trig_name.push_back("HLT_IsoMu27_v");						// 20
  trig_name.push_back("HLT_Mu50_v");						// 21
  trig_name.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");			// 22
  trig_name.push_back("HLT_Ele32_eta2p1_WPLoose_Gsf_v");			// 23
  trig_name.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");			// 24
  trig_name.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v");		// 25
  trig_name.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v");		// 26

}


bmaker_basic::~bmaker_basic(){
  outfile->cd();
  baby.tree_.SetDirectory(outfile);
  baby.Write();

  string commit_s = execute("git rev-parse HEAD");
  while(!commit_s.empty() && commit_s.at(commit_s.length()-1) == '\n') commit_s.erase(commit_s.length()-1);
  TString commit = commit_s;
  TString type = baby.Type();
  TString root_version = gROOT->GetVersion();
  TString root_tutorial_dir = gROOT->GetTutorialsDir();

  time_t curTime;
  time(&curTime);
  int seconds(floor(difftime(curTime,startTime)+0.5));
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
  treeglobal.Branch("xsec", &xsec);
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
  cout<<endl<<"Written "<<nevents<<" events in "<<outname<<". It took "<<seconds<<" seconds to run ("<<runtime<<"), "
      <<roundNumber(hertz,1)<<" Hz and "<<roundNumber(1000,2,hertz)<<" ms."<<endl<<endl;
  delete outfile;
}


// ------------ method called once each job just before starting event loop  ------------
void bmaker_basic::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void bmaker_basic::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void 
bmaker_basic::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
bmaker_basic::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
bmaker_basic::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
bmaker_basic::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bmaker_basic::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bmaker_basic);
