//// BMAKER_BASIC: Creates baby tree with basic branches
//// Function names follow the first-lowercase, following words-uppercase. No underscores


// System include files
#include <memory>
#include <iostream>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// FW physics include files

// ROOT include files
#include "TFile.h"

// User include files
#include "babymaker/bmaker/interface/bmaker_basic.hh"
#include "babymaker/bmaker/interface/baby_basic.hh"
#include "babymaker/bmaker/interface/phys_objects.hh"

using namespace std;
using namespace phys_objects;

///////////////////////// analyze: METHOD CALLED EACH EVENT ///////////////////////////
void bmaker_basic::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  ////////////////////// Event info /////////////////////
  baby.run() = iEvent.id().run();
  baby.event() = iEvent.id().event();
  baby.lumiblock() = iEvent.luminosityBlock();

  ////////////////////// Primary vertices /////////////////////
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtx);

  ////////////////////// Trigger /////////////////////
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"),triggerBits);  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByLabel("patTrigger",triggerPrescales);  
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  writeTriggers(baby, names, triggerBits, triggerPrescales);

  //////////////// HLT objects //////////////////
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByLabel("selectedPatTrigger",triggerObjects);  
  writeHLTObjects(baby, names, triggerObjects);

  ///////////////////// Filters ///////////////////////
  edm::Handle<edm::TriggerResults> filterBits;
  std::string processLabel = iEvent.isRealData() ? "RECO":"PAT"; // prompt reco runs in the "RECO" process
  iEvent.getByLabel(edm::InputTag("TriggerResults","",processLabel),filterBits);  
  // re-recoed data will have the process label "PAT" rather than "RECO";
  // if the attempt to find data with "RECO" process fails, try "PAT"
  if(!filterBits.isValid() && iEvent.isRealData()) 
      iEvent.getByLabel(edm::InputTag("TriggerResults", "", "PAT"),filterBits);  
  const edm::TriggerNames &fnames = iEvent.triggerNames(*filterBits);
  writeFilters(baby, fnames, filterBits, vtx);
  // the HBHE noise filter needs to be recomputed in early 2015 data
  edm::Handle<bool> filter_hbhe;
  if(iEvent.isRealData() && iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult",filter_hbhe)) { 
    if(*filter_hbhe) baby.pass_hbhe() = true;
    else baby.pass_hbhe() = false;
  }

  //////////////////// pfcands  //////////////////////
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByLabel("packedPFCandidates", pfcands);

  ////////////////////// Leptons /////////////////////
  vCands leptons, muons, electrons;
  edm::Handle<pat::MuonCollection> allmuons;
  iEvent.getByLabel("slimmedMuons", allmuons);
  muons = writeMuons(baby, allmuons, pfcands, vtx);
  edm::Handle<pat::ElectronCollection> allelectrons;
  iEvent.getByLabel("slimmedElectrons", allelectrons);
  electrons = writeElectrons(baby, allelectrons, pfcands, vtx);

  // Putting muons and electrons together
  baby.nleps() = baby.nmus() + baby.nels();
  baby.nvleps() = baby.nvmus() + baby.nvels();
  leptons = muons;
  leptons.insert(leptons.end(), electrons.begin(), electrons.end());

  ////////////////////// Jets /////////////////////
  edm::Handle<pat::JetCollection> alljets;
  iEvent.getByLabel("patJetsReapplyJEC", alljets);
  writeJets(baby, alljets, leptons);

  ////////////////// Filling the tree //////////////////
  baby.Fill();
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
void bmaker_basic::writeJets(baby_basic &baby, edm::Handle<pat::JetCollection> jets, vCands leptons){
  baby.njets() = 0; baby.nbl() = 0; baby.nbm() = 0;  baby.nbt() = 0;  
  baby.ht() = 0.;
  for (unsigned int ijet(0); ijet < jets->size(); ijet++) {
    const pat::Jet &jet = (*jets)[ijet];
    if(!isGoodJet(jet, JetPtCut, JetEtaCut, leptons)) continue;

    baby.jets_pt().push_back(jet.pt());
    baby.jets_eta().push_back(jet.eta());
    baby.jets_phi().push_back(jet.phi());
    
    float csv(jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    baby.jets_csv().push_back(csv);
    if(csv>CSVLoose)  baby.nbl()++;
    if(csv>CSVMedium) baby.nbm()++;
    if(csv>CSVTight)  baby.nbt()++;

    baby.njets()++;
    baby.ht() += jet.pt();
  } // Loop over jets  
} 

vCands bmaker_basic::writeMuons(baby_basic &baby, edm::Handle<pat::MuonCollection> muons, 
				edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<reco::VertexCollection> vtx){
  vCands signalmuons; 
  baby.nmus() = 0; baby.nvmus() = 0;
  for (unsigned int ilep(0); ilep < muons->size(); ilep++) {
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


vCands bmaker_basic::writeElectrons(baby_basic &baby, edm::Handle<pat::ElectronCollection> electrons, 
				    edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<reco::VertexCollection> vtx){
  vCands signalelectrons; 
  baby.nels() = 0; baby.nvels() = 0;
  for (unsigned int ilep(0); ilep < electrons->size(); ilep++) {
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

void bmaker_basic::writeTriggers(baby_basic &baby, const edm::TriggerNames &names, 
                                 edm::Handle<edm::TriggerResults> triggerBits, 
                                 edm::Handle<pat::PackedTriggerPrescales> triggerPrescales){
  baby.trig().resize(trig_name.size(), false);
  baby.trig_prescale().resize(trig_name.size(), -1.);
  for (unsigned int itrig(0); itrig < triggerBits->size(); itrig++) {
    for(unsigned itn(0); itn < trig_name.size(); itn++){
      if(names.triggerName(itrig).find(trig_name[itn])!=string::npos){
        baby.trig()[itn] = triggerBits->accept(itrig);
        baby.trig_prescale()[itn] = triggerPrescales->getPrescaleForIndex(itrig);
      }
    }
  }
}

void bmaker_basic::writeHLTObjects(baby_basic &baby, const edm::TriggerNames &names, 
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

void bmaker_basic::writeFilters(baby_basic &baby, const edm::TriggerNames &fnames,
                                  edm::Handle<edm::TriggerResults> filterBits,
                                  edm::Handle<reco::VertexCollection> vtx){
  for (unsigned int i(0); i < filterBits->size(); ++i) {
    string name = fnames.triggerName(i);
    bool pass = static_cast<bool>(filterBits->accept(i));
    if (name=="Flag_goodVertices") baby.pass_goodv() = pass;
    else if (name=="Flag_CSCTightHaloFilter") baby.pass_cschalo() = pass;
    else if (name=="Flag_eeBadScFilter") baby.pass_eebadsc() = pass;
    else if (name=="Flag_HBHENoiseFilter") baby.pass_hbhe() = pass;
  }

  baby.pass_goodv() &= hasGoodPV(vtx);
  baby.pass_jets() = true; //FIXME

  baby.pass() = baby.pass_hbhe() && baby.pass_goodv() && baby.pass_cschalo() && baby.pass_eebadsc() && baby.pass_jets();
}

/*
 _____                 _                   _                 
/  __ \               | |                 | |                
| /  \/ ___  _ __  ___| |_ _ __ _   _  ___| |_ ___  _ __ ___ 
| |    / _ \| '_ \/ __| __| '__| | | |/ __| __/ _ \| '__/ __|
| \__/\ (_) | | | \__ \ |_| |  | |_| | (__| || (_) | |  \__ \
 \____/\___/|_| |_|___/\__|_|   \__,_|\___|\__\___/|_|  |___/
*/

bmaker_basic::bmaker_basic(const edm::ParameterSet& iConfig){
  outname = TString(iConfig.getParameter<string>("outputFile"));
  outfile = new TFile(outname, "recreate");
  outfile->cd();
  baby.tree_.SetDirectory(outfile);

  trig_name = vector<TString>();
  trig_name.push_back("HLT_PFHT350_PFMET100_NoiseCleaned_v");     // 0
  trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET70_v");      // 1
  trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");        // 2
  trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_v");    // 3
  trig_name.push_back("HLT_Mu15_PFHT300_v");          // 4
  trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET70_v");     // 5
  trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");       // 6
  trig_name.push_back("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_v");   // 7
  trig_name.push_back("HLT_Ele15_PFHT300_v");         // 8
  trig_name.push_back("HLT_DoubleMu8_Mass8_PFHT300_v");       // 9
  trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v"); // 10
  trig_name.push_back("HLT_PFHT475_v");           // 11
  trig_name.push_back("HLT_PFHT800_v");           // 12
  trig_name.push_back("HLT_PFMET120_NoiseCleaned_Mu5_v");     // 13
  trig_name.push_back("HLT_PFMET170_NoiseCleaned_v");       // 14
  trig_name.push_back("HLT_DoubleIsoMu17_eta2p1_v");        // 15
  trig_name.push_back("HLT_Mu17_TrkIsoVVL_v");          // 16
  trig_name.push_back("HLT_IsoMu17_eta2p1_v");          // 17
  trig_name.push_back("HLT_IsoMu20_v");           // 18
  trig_name.push_back("HLT_IsoMu24_eta2p1_v");          // 19
  trig_name.push_back("HLT_IsoMu27_v");           // 20
  trig_name.push_back("HLT_Mu50_v");            // 21
  trig_name.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");      // 22
  trig_name.push_back("HLT_Ele32_eta2p1_WPLoose_Gsf_v");      // 23
  trig_name.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");     // 24
  trig_name.push_back("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v");   // 25
  trig_name.push_back("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v");   // 26

}


bmaker_basic::~bmaker_basic(){
  outfile->cd();
  baby.tree_.SetDirectory(outfile);
  baby.Write();

  TTree treeglobal("treeglobal", "treeglobal");
  // treeglobal.Branch("nev_file", &num_entries);
  // treeglobal.Branch("nev_sample", &num_total_entries);
  // treeglobal.Branch("commit", &commit);
  // treeglobal.Branch("model", &model);
  // treeglobal.Branch("type", &type);
  // treeglobal.Branch("root_version", &root_version);
  // treeglobal.Branch("root_tutorial_dir", &root_tutorial_dir);
  treeglobal.Branch("trig_name", &trig_name);
  treeglobal.Fill();
  treeglobal.SetDirectory(outfile);
  treeglobal.Write();
  
  outfile->Close();
  cout<<endl<<"Written baby in "<<outname<<endl<<endl;
  delete outfile;
}


// ------------ method called once each job just before starting event loop  ------------
void bmaker_basic::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bmaker_basic::endJob() 
{
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
