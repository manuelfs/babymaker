// Original Author:  Manuel Franco Sevilla
//         Created:  Mon, 14 Sep 2015 08:56:08 GMT


#ifndef H_BMAKER
#define H_BMAKER

// System include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ROOT include files
#include "TTree.h"

// User include files
#include "babymaker/bmaker/interface/baby_basic.hh"

// Class declaration

class bmaker : public edm::EDAnalyzer {
public:
  explicit bmaker(const edm::ParameterSet&);
  ~bmaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  TFile * outfile;
  baby_basic baby;
  long nevent;


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
};

#endif
