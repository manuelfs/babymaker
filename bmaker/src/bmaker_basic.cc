// Original Author:  Manuel Franco Sevilla
//         Created:  Mon, 14 Sep 2015 08:56:08 GMT


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

using namespace std;

// ------------ method called for each event  ------------
void bmaker_basic::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //////////////// Jets //////////////////
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel("slimmedJets", jets);
  baby.njets() = jets->size();
  cout<<"Found "<<baby.njets()<<" jets"<<endl;

  // Filling the tree
  baby.Fill();
}


// Constructors and destructor
bmaker_basic::bmaker_basic(const edm::ParameterSet& iConfig){
  outfile = new TFile("baby.root", "recreate");
  outfile->cd();
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
