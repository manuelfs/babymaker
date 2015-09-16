// Creates baby tree with basic branches


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
#include "DataFormats/PatCandidates/interface/Jet.h"

// ROOT include files
#include "TFile.h"

// User include files
#include "babymaker/bmaker/interface/bmaker_basic.hh"
#include "babymaker/bmaker/interface/baby_basic.hh"
#include "babymaker/bmaker/interface/phys_objects.hh"

using namespace std;
using namespace phys_objects;

float bmaker_basic::MinSignalLeptonPt = 20.0;
float bmaker_basic::MinVetoLeptonPt = 10.0;
float bmaker_basic::MuMiniIsoCut = 0.2;
float bmaker_basic::ElMiniIsoCut = 0.1;

///////////////////////// analyze: METHOD CALLED EACH EVENT ///////////////////////////
void bmaker_basic::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  ////////////////////// Event info /////////////////////
  baby.run() = iEvent.id().run();
  baby.event() = iEvent.id().event();
  baby.lumiblock() = iEvent.luminosityBlock();

  ////////////////////// Primary vertices /////////////////////
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtx);

  //////////////////// pfcands  //////////////////////
  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByLabel("packedPFCandidates", pfcands);

  ////////////////////// Leptons /////////////////////
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByLabel("slimmedMuons", muons);
  WriteMuons(baby, muons, pfcands, vtx);
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByLabel("slimmedElectrons", electrons);
  WriteElectrons(baby, electrons, pfcands, vtx);

  baby.nleps() = baby.nmus() + baby.nels();
  baby.nvleps() = baby.nvmus() + baby.nvels();

  ////////////////////// Jets /////////////////////
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel("slimmedJets", jets);
  baby.njets() = jets->size();

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
void bmaker_basic::WriteMuons(baby_basic &baby, edm::Handle<pat::MuonCollection> muons, 
			      edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<reco::VertexCollection> vtx){
  baby.nmus() = 0; baby.nvmus() = 0;
  for (unsigned int ilep(0); ilep < muons->size(); ilep++) {
    const pat::Muon &lep = (*muons)[ilep];    
    double dz(0.), d0(0.);
    if(lep.track().isAvailable()){
      dz = lep.track()->vz()-vtx->at(0).z();
      d0 = lep.track()->d0()-vtx->at(0).x()*sin(lep.track()->phi())+vtx->at(0).y()*cos(lep.track()->phi());
    } //else cout<<ilep<<": (pt,eta,phi) = ("<<lep.pt()<<", "<<lep.eta()<<", "<<lep.phi()<<"). Is loose "<<lep.isLooseMuon()
    //       <<", is tight "<<lep.isTightMuon(vtx->at(0))<<endl;
    if(!lep.isLooseMuon() || lep.pt() <= MinVetoLeptonPt || fabs(lep.eta()) > 2.4 || 
       fabs(dz) > 0.5 || fabs(d0) > 0.2) continue;

    baby.mus_pt().push_back(lep.pt());
    baby.mus_eta().push_back(lep.eta());
    baby.mus_phi().push_back(lep.phi());
    baby.mus_dz().push_back(dz);
    baby.mus_d0().push_back(d0);
    baby.mus_charge().push_back(lep.charge());
    baby.mus_medium().push_back(lep.isMediumMuon());
    baby.mus_tight().push_back(lep.isTightMuon(vtx->at(0)));
    baby.mus_miniso().push_back(getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., false));

    if(baby.mus_miniso().back() < MuMiniIsoCut){
      baby.nvmus()++;
      if(baby.mus_medium().back() && lep.pt() > MinSignalLeptonPt) baby.nmus()++;
    }
  } // Loop over muons

}


void bmaker_basic::WriteElectrons(baby_basic &baby, edm::Handle<pat::ElectronCollection> electrons, 
				  edm::Handle<pat::PackedCandidateCollection> pfcands, edm::Handle<reco::VertexCollection> vtx){
  baby.nels() = 0; baby.nvels() = 0;
  for (unsigned int ilep(0); ilep < electrons->size(); ilep++) {
    const pat::Electron &lep = (*electrons)[ilep];    
    double dz(0.), d0(0.);
    if(lep.gsfTrack().isAvailable()){
      dz = lep.gsfTrack()->vz()-vtx->at(0).z();
      d0 = lep.gsfTrack()->d0()-vtx->at(0).x()*sin(lep.gsfTrack()->phi())+vtx->at(0).y()*cos(lep.gsfTrack()->phi());
    } 
    if(!IdElectron(lep, kVeto, vtx, false) || lep.pt() <= MinVetoLeptonPt || fabs(lep.eta()) > 2.5) continue;

    baby.els_pt().push_back(lep.pt());
    baby.els_eta().push_back(lep.eta());
    baby.els_phi().push_back(lep.phi());
    baby.els_dz().push_back(dz);
    baby.els_d0().push_back(d0);
    baby.els_charge().push_back(lep.charge());
    baby.els_medium().push_back(IdElectron(lep, kMedium, vtx, false));
    baby.els_tight().push_back(IdElectron(lep, kTight, vtx, false));
    baby.els_miniso().push_back(getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&lep), 0.05, 0.2, 10., false));

    if(baby.els_miniso().back() < ElMiniIsoCut){
      baby.nvels()++;
      if(baby.els_medium().back() && lep.pt() > MinSignalLeptonPt) baby.nels()++;
    }
  } // Loop over electrons

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
  outfile = new TFile(TString(iConfig.getParameter<std::string>("outputFile")), "recreate");
  outfile->cd();
  baby.tree_.SetDirectory(outfile);
}


bmaker_basic::~bmaker_basic(){
  baby.Write();
  outfile->Close();

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
