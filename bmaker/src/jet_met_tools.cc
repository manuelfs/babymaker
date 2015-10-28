// jet_met_tools: Functions related to jets, MET, and JECs

// System include files
#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <DataFormats/Math/interface/deltaR.h>

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
      if(deltaR(jet, *leptons[ilep]) < sizeJet) return true;
    }
  } // Loop over leptons

  return false;
}

// Loose jet matching from RA2/b just for cleaning
bool jet_met_tools::jetMatched(const pat::Jet &jet, vCands objects){
  int nph(0); // Just cleaning the first photon with 100 GeV
  for(unsigned ind(0); ind < objects.size(); ind++){
    double dr(deltaR(jet, *(objects[ind])));
    double drelpt(fabs((objects[ind]->pt() - jet.pt())/objects[ind]->pt()));
    if(objects[ind]->pdgId()==22) {
      if(objects[ind]->pt()>100){
	nph++;
	if(nph>1) continue;
      } else continue;
    } // If it is a photon
    if(drelpt < 1. && dr < sizeJet) return true;
  } // Loop over objects
  return false;
}

float jet_met_tools::getGenPt(const pat::Jet &jet, edm::Handle<edm::View <reco::GenJet> > genjets){
  for (size_t ijet(0); ijet < genjets->size(); ijet++) {
    const reco::GenJet &genjet = (*genjets)[ijet];
    double dr(deltaR(jet, genjet));
    if(dr < 0.2) return genjet.pt();
  }
  return -99999.;    
}

bool jet_met_tools::isLowDphi(vCands jets, float mht_phi, float &dphi1, float &dphi2, float &dphi3, float &dphi4){
  dphi1 = 10.; dphi2 = 10.; dphi3 = 10.; dphi4 = 10.; 
  float *dphi[] = {&dphi1, &dphi2, &dphi3, &dphi4}; 
  for(unsigned ind(0); ind < jets.size() && ind < 4; ind++)
    *(dphi[ind]) = abs(reco::deltaPhi(jets[ind]->phi(), mht_phi));
  return (dphi1<0.5 || dphi2<0.5 || dphi3<0.3 || dphi4<0.3);
}

float jet_met_tools::trueHT(edm::Handle<edm::View <reco::GenJet> > genjets){
  float ht(0.);
  for (size_t ijet(0); ijet < genjets->size(); ijet++) {
    const reco::GenJet &jet = (*genjets)[ijet];
    if(jet.pt() > JetHLTPtCut && fabs(jet.eta()) <= JetHLTEtaCut) ht += jet.pt();
  }
  return ht;
}

bool jet_met_tools::idJet(const pat::Jet &jet, CutLevel cut){
  // From https://twiki.cern.ch/twiki/bin/view/CMS/JetID
  double eta = jet.eta();
  double NHF = jet.neutralHadronEnergyFraction();
  double NEMF = jet.neutralEmEnergyFraction();
  double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
  double NumNeutralParticles =jet.neutralMultiplicity();

  double CHF = jet.chargedHadronEnergyFraction();
  double CHM = jet.chargedMultiplicity(); 
  double CEMF = jet.chargedEmEnergyFraction();
  
  if(cut == kPBNR){  // RA2/b's PBNR and old Jet ID
    bool eta_l_2p4 =  NumConst>=2 && NHF<0.9 && NEMF<0.95 && CHM>0 && CHF>0 && CEMF<0.99;
    bool eta_geq_2p4 =  NHF<0.9 && NEMF<0.95 && NumConst>=2;
    return (eta_l_2p4 && fabs(eta)<2.4) || (fabs(eta)>=2.4 && eta_geq_2p4);   
  }

  double NHFCut, NEMFCut, NumConstCut, CHFCut, CHMCut, CEMFCut, NEMFCut_HF, NumNeuCut;
  switch(cut){
  case kTight:
    NHFCut	= 0.9;
    NEMFCut	= 0.9;
    NumConstCut = 1;

    CHFCut	= 0;
    CHMCut	= 0;
    CEMFCut	= 0.99;

    NEMFCut_HF	= 0.9;
    NumNeuCut	= 10;
    break;
  case kLoose:
  default:
    NHFCut	= 0.99;
    NEMFCut	= 0.99;
    NumConstCut = 1;

    CHFCut	= 0;
    CHMCut	= 0;
    CEMFCut	= 0.99;

    NEMFCut_HF	= 0.9;
    NumNeuCut	= 10;
    break;
  }
    
  bool eta_leq_3 = (NHF<NHFCut && NEMF<NEMFCut && NumConst>NumConstCut) && 
    ((fabs(eta)<=2.4 && CHF>CHFCut && CHM>CHMCut && CEMF<CEMFCut) || fabs(eta)>2.4);
  bool eta_g_3 = NEMF<NEMFCut_HF && NumNeutralParticles>NumNeuCut;

  return  (eta_leq_3 && fabs(eta)<=3.) || (eta_g_3 && fabs(eta)>3.);  // Official recommendation
}


void jet_met_tools::getJetCorrections(edm::Handle<edm::View <reco::GenJet> > genjets, edm::Handle<pat::JetCollection> alljets, double rhoEvent){
  jetTotCorrections.resize(alljets->size(), 1.);
  jetL1Corrections.resize(alljets->size(), 1.);
  corrJet.clear();
  genJetPt.clear();
  if(!doJEC) {
    for (size_t ijet(0); ijet < alljets->size(); ijet++) {
      const pat::Jet &jet = (*alljets)[ijet];
      genJetPt.push_back(getGenPt(jet, genjets));
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
    
    //smear jets
    bool doSmearJets = true;
    genJetPt.push_back(getGenPt(jet, genjets));
    if (doSmearJets){
      if (genJetPt[ijet]>0.) {
        float corr_pt = jet.p4().pt()*rawFactor*jetTotCorrections[ijet];
        float smeared_pt = genJetPt[ijet] + getJetResolutionSF(jet.eta())*(corr_pt - genJetPt[ijet]);
        if (smeared_pt < 0.) smeared_pt = 0.;
        jetTotCorrections[ijet] *= smeared_pt/corr_pt;
      }
    }  
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


// from 8TeV dijet measurement with an extra 50% 
// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#2012
float jet_met_tools::getJetResolutionSF(float jet_eta){
  double abseta = fabs(jet_eta);
  if (abseta < 0.5) return (1. + 1.5*0.079);
  else if (abseta < 1.1) return (1. + 1.5*0.099);
  else if (abseta < 1.7) return (1. + 1.5*0.121);
  else if (abseta < 2.3) return (1. + 1.5*0.208);
  else if (abseta < 2.8) return (1. + 1.5*0.254);
  else if (abseta < 3.2) return (1. + 1.5*0.395);
  else if (abseta < 5.0) return (1. + 1.5*0.056);
  else return 1.;
}

jet_met_tools::~jet_met_tools(){
  if(doJEC) delete jetCorrector;
}

