//// PHYS_OBJECTS: Common physics objects definitions
//// Function names follow the first-lowercase, following words-uppercase. No underscores

// System include files
#include <algorithm>

// FW include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

// User include files
#include "babymaker/bmaker/interface/phys_objects.hh"
#include "babymaker/bmaker/interface/in_json.hh"

using namespace std;

namespace phys_objects{

  const std::vector<std::vector<int> > VRunLumi2015golden(MakeVRunLumi("2015golden"));
  const std::vector<std::vector<int> > VRunLumi2015nohfgolden(MakeVRunLumi("2015nohfgolden"));
  const std::vector<std::vector<int> > VRunLumi2015dcs(MakeVRunLumi("2015dcs"));

  bool isInJSON(string type, int run, int lumiblock){
    if(type=="golden") return inJSON(VRunLumi2015golden, run, lumiblock);
    if(type=="nohf_golden") return inJSON(VRunLumi2015nohfgolden, run, lumiblock);
    if(type=="dcs") return inJSON(VRunLumi2015dcs, run, lumiblock);

    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// JETS /////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  bool leptonInJet(const pat::Jet &jet, vCands leptons){
    for(unsigned ilep(0); ilep < leptons.size(); ilep++){
      int indpf(-1);
      unsigned npflep(leptons[ilep]->numberOfSourceCandidatePtrs());
      if(leptons[ilep]->isMuon() && npflep==1) indpf = 0;
      if(leptons[ilep]->isElectron() && npflep==2) indpf = 1; // Electrons have a missing reference at 0
      if(indpf>=0){ // The lepton is PF -> looping over PF cands in jet
	for (unsigned ijet(0); ijet < jet.numberOfSourceCandidatePtrs(); ijet++) 
	  if(jet.sourceCandidatePtr(ijet) == leptons[ilep]->sourceCandidatePtr(indpf))
	    return true;
      } else { // The lepton is not PF, matching with deltaR
	if(deltaR(jet, *leptons[ilep]) < 0.4) return true;
      }
    } // Loop over leptons

    return false;
  }

  bool idJet(const pat::Jet &jet){
    //LooseID from https://twiki.cern.ch/twiki/bin/view/CMS/JetID
    double eta = jet.eta();
    double NHF = jet.neutralHadronEnergyFraction();
    double NEMF = jet.neutralEmEnergyFraction();
    double CHF = jet.chargedHadronEnergyFraction();
    double CEMF = jet.chargedEmEnergyFraction();
    double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    double NumNeutralParticles =jet.neutralMultiplicity();
    double CHM = jet.chargedMultiplicity(); 
    //double MUF = jet.muonEnergyFraction(); //Only used in TightID
    
    bool eta_leq_3 = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(eta)>2.4);
    bool eta_g_3 = NEMF<0.90 && NumNeutralParticles>10;

    return  (eta_leq_3 && fabs(eta)<=3.) || (eta_g_3 && fabs(eta)>3.);
  }

  bool hasGoodPV(edm::Handle<reco::VertexCollection> vtx){
    bool one_good_pv(false);
    for(unsigned ipv(0); ipv < vtx->size(); ipv++){
      const double pv_rho(sqrt(pow(vtx->at(ipv).x(),2) + pow(vtx->at(ipv).y(),2)));
      if(vtx->at(ipv).ndof()>4 && fabs(vtx->at(ipv).z())<24. && pv_rho<2.0 && vtx->at(ipv).isFake()==0){
        one_good_pv = true;
        break;
      }
    } // Loop over vertices
    return one_good_pv;
  }
}
