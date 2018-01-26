// photon_tools: Functions related to RA2/b photons
// Copied from RA2/b's TreeMaker
// https://github.com/TreeMaker/TreeMaker/blob/74848c2a9fd16acaaf8d711842c6e2dc1d5bd56b/Utils/src/PhotonIDisoProducer.cc

// System include files
#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>

// FW physics include files
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// User include files
#include "babymaker/bmaker/interface/release.hh"
#include "babymaker/bmaker/interface/photon_tools.hh"


using namespace std;

bool photon_tools::idPhoton(const pat::Photon &photon, edm::Handle<std::vector<pat::Electron> > &electrons,
			    edm::Handle<reco::ConversionCollection> &conversions,
			    edm::Handle<reco::BeamSpot> &beamspot, double rho){
  double eta(fabs(photon.eta()));
  bool barrel(eta<1.4442), endcap(eta>1.566 && eta<2.5);

  // Barrel cuts
  double hadTowOverEmCut(0.028), sigmaIetaIetaCut(0.0107);
  double chIsoCut(2.67), nuIsoCut(7.23 + exp(0.0028*(photon.pt()+0.5408))), gamIsoCut(2.11 + 0.0014*(photon.pt()));
  if(endcap) { // Endcap cuts
    hadTowOverEmCut  = 0.093;
    sigmaIetaIetaCut = 0.0272;
    chIsoCut         = 1.79;
    nuIsoCut         = 8.89 + 0.01725*(photon.pt());
    gamIsoCut        = 3.09 + 0.0091*(photon.pt());
  }

  ////// Check acceptance
  if(!barrel && !endcap) return false;
  ////// Check ID
  if(photon.hadTowOverEm() >= hadTowOverEmCut || photon.sigmaIetaIeta() >= sigmaIetaIetaCut) return false;
  if(electronMatch(photon.superCluster(), electrons, conversions, beamspot->position())) return false;
  ////// Check iso
  double chIso  = rhoCorrectedIso(pfCh , photon.chargedHadronIso(), photon.eta(), rho); 
  double nuIso  = rhoCorrectedIso(pfNu , photon.neutralHadronIso(), photon.eta(), rho); 
  double gamIso = rhoCorrectedIso(pfGam, photon.photonIso()       , photon.eta(), rho); 
  if(chIso >= chIsoCut || nuIso >= nuIsoCut || gamIso >= gamIsoCut) return false;

  return true;
}

// From https://github.com/RazorCMS/SUSYBSMAnalysis-RazorTuplizer/blob/6072ffb43bbeb3f6b34cf8a96426c7f104c5b902/plugins/RazorAux.cc#L127
//check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
//and not matching any conversion in the collection passing the quality cuts
bool photon_tools::electronMatch(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<pat::Electron> > &electrons,
				 const edm::Handle<reco::ConversionCollection> &conversions, const math::XYZPoint &beamspot,
				 double lxyMin, double probMin, unsigned int nHitsBeforeVtxMax) {

  if (sc.isNull()) return false;
  for (std::vector<pat::Electron>::const_iterator it = electrons->begin(); it!=electrons->end(); ++it) {
    //match electron to supercluster
    if (it->superCluster()!=sc) continue;
    //check expected inner hits
#ifdef CMSSW_9_4
    if (it->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
#else 
    if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
#endif
    //check if electron is matching to a conversion
    if (ConversionTools::hasMatchedConversion(*it,conversions,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
    return true;
  }
  return false;
}

double photon_tools::rhoCorrectedIso(isoType isoType_, double isoVal, double eta, double rho){
  int iEta = -1;
  for(unsigned int i = 0; i < etaHigh.size(); i++)   
    if(fabs(eta) < etaHigh[i] && fabs(eta) > etaLow[i]){
      iEta = i;
      break;
    }

  if(iEta < 0){
    std::cout << "AHHHHH: couldn't match eta region... eta = " << fabs(eta) << std::endl;
    return 99999. ;
  } else return max(isoVal - effA[iEta][isoType_]*rho , 0.);
};

void photon_tools::addEffA(double etaLow_,double etaHigh_,double effA_pfCh_,double effA_pfNu_,double effA_pfGa_){
  etaHigh.push_back( etaHigh_ );
  etaLow.push_back( etaLow_ );
  vector<double> temp;
  temp.push_back(effA_pfCh_);
  temp.push_back(effA_pfNu_);
  temp.push_back(effA_pfGa_);
  effA.push_back( temp );

};

photon_tools::photon_tools(){
  addEffA( 0.0, 1.0, 0.0234, 0.0053, 0.078 );
  addEffA( 1.0, 1.479, 0.0189, 0.0130, 0.0629 );
  addEffA( 1.479, 2.0, 0.0171, 0.0057, 0.0264 );
  addEffA( 2.0, 2.2, 0.0129, 0.0070, 0.0462 );
  addEffA( 2.2, 2.3, 0.0110, 0.0152, 0.0740 );
  addEffA( 2.3, 2.4, 0.0074, 0.0232, 0.0924 );
  addEffA( 2.4, 99., 0.0035, 0.1709, 0.1484 );
}

photon_tools::~photon_tools(){
}

