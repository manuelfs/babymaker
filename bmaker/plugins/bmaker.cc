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

// User include files
#include "babymaker/bmaker/interface/bmaker.hh"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


// Constructors and destructor
bmaker::bmaker(const edm::ParameterSet& iConfig){
  // edm::Service<TFileService> fs;      
  // tree_=fs->make<TTree>("tree_global","tree_global tree");
  outfile = new TFile("baby.root", "recreate");
  outfile->cd();
  tree_ = new TTree("tree_global","tree_global");
  tree_->Branch("nevent", &nevent, "nevent/I");
  nevent = 0;
}


bmaker::~bmaker(){
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  tree_->Write();
  outfile->Close();

  delete tree_;
  delete outfile;
}


// ------------ method called for each event  ------------
void bmaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   nevent++;
   std::cout<<nevent<<": Hola. Header file!"<<std::endl;
   tree_->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void bmaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void 
bmaker::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
bmaker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
bmaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
bmaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
bmaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bmaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(bmaker);
