// jet_met_tools: Functions related to jets, MET, and JECs

// System include files
#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

// User include files
#include "babymaker/bmaker/interface/jet_met_tools.hh"
#include "babymaker/bmaker/interface/release.hh"

using namespace std;
using namespace utilities;


bool jet_met_tools::leptonInJet(const pat::Jet &jet, vCands leptons){
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

bool jet_met_tools::idJet(const pat::Jet &jet){
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


void jet_met_tools::getJetCorrections(edm::Handle<pat::JetCollection> alljets, double rhoEvent){
  jetTotCorrections.resize(alljets->size(), 1.);
  jetL1Corrections.resize(alljets->size(), 1.);
  corrJet.clear();
  if(!doJEC) {
    for (size_t ijet(0); ijet < alljets->size(); ijet++) {
      const pat::Jet &jet = (*alljets)[ijet];
      corrJet.push_back(jet.p4());
    }
    return;
  }
  
  rhoEvent_ = rhoEvent;
  alljets_ = alljets;
  for (size_t ijet(0); ijet < alljets->size(); ijet++) {
    const pat::Jet &jet = (*alljets)[ijet];
    float rawFactor(jet.jecFactor("Uncorrected"));
    jetValues.setJetPt(jet.pt()*rawFactor);
    jetValues.setJetEta(jet.eta());
    jetValues.setJetA(jet.jetArea());
    jetValues.setRho(rhoEvent);
    vector<float> corr_vals = jetCorrector->getSubCorrections(jetValues);
    jetTotCorrections[ijet] = corr_vals.at(corr_vals.size()-1);      // All corrections
    jetL1Corrections[ijet] = corr_vals.at(0);                        // L1 PU correction (offset)
    corrJet.push_back(jet.p4()*rawFactor*jetTotCorrections[ijet]);   // LorentzVecor with all corrections * raw factor
  } // Loop over alljets
    
}

void jet_met_tools::getMETRaw(edm::Handle<pat::METCollection> mets, float &metRaw, float &metRawPhi){
#ifdef PRE_7_4_12
  metRaw = mets->at(0).uncorrectedPt();
  metRawPhi = mets->at(0).uncorrectedPhi();
#else
  metRaw = mets->at(0).uncorPt();
  metRawPhi = mets->at(0).uncorPhi();
#endif

}

void jet_met_tools::getMETWithJEC(edm::Handle<pat::METCollection> mets, float &met, float &metPhi){
  if(!doJEC) {
    met = mets->at(0).pt();
    metPhi = mets->at(0).phi();
    return;
  }

  float metRaw, metRawPhi;
  getMETRaw(mets, metRaw, metRawPhi);
  float metx(metRaw*cos(metRawPhi)), mety(metRaw*sin(metRawPhi));

  // Code to skip muons and EM from
  // https://github.com/cms-sw/cmssw/blob/8d582ad580a446865fc58675d16f6cdf2dae3605/JetMETCorrections/Type1MET/interface/PFJetMETcorrInputProducerT.h#L173-L184
  StringCutObjectSelector<reco::Candidate> skipMuonSelection("isGlobalMuon | isStandAloneMuon",true);
  for (size_t ijet(0); ijet < alljets_->size(); ijet++) {
    const pat::Jet &jet = (*alljets_)[ijet];

    double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
    if(emEnergyFraction > 0.90 || fabs(jet.eta()) > 9.9) continue;
    
    reco::Candidate::LorentzVector rawJetP4 = jet.p4()*jet.jecFactor("Uncorrected");
    float totCorr(jetTotCorrections[ijet]), l1Corr(jetL1Corrections[ijet]);

    const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
    for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
	  cand != cands.end(); ++cand ) {
      const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
      const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
      if ( mu != 0 && skipMuonSelection(*mu) ) {
	reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
	rawJetP4 -= muonP4;
	jetValues.setJetPt(rawJetP4.pt());
	jetValues.setJetEta(rawJetP4.eta());
	jetValues.setJetA(jet.jetArea());
	jetValues.setRho(rhoEvent_);
	vector<float> corr_vals = jetCorrector->getSubCorrections(jetValues);
	totCorr = corr_vals.at(corr_vals.size()-1); // All corrections
	l1Corr = corr_vals.at(0);      // L1 PU correction (offset)
      }
    }

    if((rawJetP4.pt()*totCorr) <= 15.) continue;
    metx -= rawJetP4.px()*(totCorr - l1Corr);
    mety -= rawJetP4.py()*(totCorr - l1Corr);
  } // Loop over alljets_

  met = hypot(metx,mety);
  metPhi = atan2(mety,metx);
}

jet_met_tools::jet_met_tools(TString ijecName):
  jecName(ijecName){

  doJEC = !jecName.Contains("miniAOD");
  if(doJEC) {
    vector<JetCorrectorParameters> jecFiles;
    string basename(getenv("CMSSW_BASE")); basename += "/src/babymaker/txt/jec/"; basename += jecName.Data();
    cout<<endl<<"BABYMAKER: jet_met_tools: Applying JECs on-the-fly with files "<<basename.c_str()<<"*.txt"<<endl<<endl;
    jecFiles.push_back(JetCorrectorParameters(basename+"_L1FastJet_AK4PFchs.txt"));
    jecFiles.push_back(JetCorrectorParameters(basename+"_L2Relative_AK4PFchs.txt"));
    jecFiles.push_back(JetCorrectorParameters(basename+"_L3Absolute_AK4PFchs.txt"));
    if(jecName.Contains("DATA")) jecFiles.push_back(JetCorrectorParameters(basename+"_L2L3Residual_AK4PFchs.txt"));
    jetCorrector = new FactorizedJetCorrectorCalculator(jecFiles);
  }
  }

jet_met_tools::~jet_met_tools(){
  delete jetCorrector;
}

