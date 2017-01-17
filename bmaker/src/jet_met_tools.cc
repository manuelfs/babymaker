// jet_met_tools: Functions related to jets, MET, and JECs

// System include files
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <DataFormats/Math/interface/deltaR.h>
#include <DataFormats/Math/interface/deltaPhi.h>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>

// User include files
#include "babymaker/bmaker/interface/jet_met_tools.hh"
#include "babymaker/bmaker/interface/release.hh"

// ROOT include files
#include "TFile.h"

using namespace std;
using namespace utilities;

const vector<BTagEntry::OperatingPoint> jet_met_tools::op_pts_{BTagEntry::OP_MEDIUM, BTagEntry::OP_LOOSE, BTagEntry::OP_TIGHT};
const vector<BTagEntry::JetFlavor> jet_met_tools::flavors_{BTagEntry::FLAV_B, BTagEntry::FLAV_C, BTagEntry::FLAV_UDSG};

namespace{
  template<typename T, typename... Args>
    std::unique_ptr<T> MakeUnique(Args&&... args){
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }
}

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
  if(!genjets.isValid()) return -99999.;
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
    NHFCut      = 0.9;
    NEMFCut     = 0.9;
    NumConstCut = 1;

    CHFCut      = 0;
    CHMCut      = 0;
    CEMFCut     = 0.99;

    NEMFCut_HF  = 0.9;
    NumNeuCut   = 10;
    break;
  case kLoose:
  default:
    NHFCut      = 0.99;
    NEMFCut     = 0.99;
    NumConstCut = 1;

    CHFCut      = 0;
    CHMCut      = 0;
    CEMFCut     = 0.99;

    NEMFCut_HF  = 0.9;
    NumNeuCut   = 10;
    break;
  }
    
  bool eta_leq_3 = (NHF<NHFCut && NEMF<NEMFCut && NumConst>NumConstCut) && 
    ((fabs(eta)<=2.4 && CHF>CHFCut && CHM>CHMCut && CEMF<CEMFCut) || fabs(eta)>2.4);
  bool eta_g_3 = NEMF<NEMFCut_HF && NumNeutralParticles>NumNeuCut;

  return  (eta_leq_3 && fabs(eta)<=3.) || (eta_g_3 && fabs(eta)>3.);  // Official recommendation
}


void jet_met_tools::getJetCorrections(edm::Handle<edm::View <reco::GenJet> > genjets, edm::Handle<pat::JetCollection> alljets, double rhoEvent){
  rhoEvent_ = rhoEvent;
  alljets_ = alljets;
  jetTotCorrections.resize(alljets->size(), 1.);
  jetL1Corrections.resize(alljets->size(), 1.);
  jecUnc.resize(alljets->size(), 0.);
  jerUnc.resize(alljets->size(), 0.);
  corrJet.clear();
  genJetPt.clear();
  if(!doJEC) {
    for (size_t ijet(0); ijet < alljets->size(); ijet++) {
      const pat::Jet &jet = (*alljets)[ijet];
      corrJet.push_back(jet.p4());
      genJetPt.push_back(getGenPt(jet, genjets));
    }
    if (doSystematics) setJetUncertainties(genjets);
    return;
  }
  
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
    
    //the genJets should be obtained even if we are not running systematics, since we save the pT resolution in the babies
    genJetPt.push_back(getGenPt(jet, genjets));
    if (doSystematics){
      //set JECs uncertainty values
      jecUncProvider->setJetEta(jet.eta());
      jecUncProvider->setJetPt(corrJet[ijet].pt());
      jecUnc[ijet] = jecUncProvider->getUncertainty(true);
      if (isData) jecUnc[ijet] = sqrt(pow(jecUnc[ijet],2) + pow((corr_vals.at(corr_vals.size()-1)/corr_vals.at(corr_vals.size()-2)-1.),2));
      //set JER uncertainty values

      float smearedJetPt(0.);
      if (genJetPt[ijet]>0) smearedJetPt = genJetPt[ijet] + getJetResolutionSF(jet.eta())*(corrJet[ijet].pt() - genJetPt[ijet]);
      jerUnc[ijet] = smearedJetPt/corrJet[ijet].pt() - 1.;
      if (smearedJetPt < 0.01) jerUnc[ijet] = 0.; // so data will not have resolution unc. 
    }
  } // Loop over alljets
    
}

void jet_met_tools::setJetUncertainties(edm::Handle<edm::View <reco::GenJet> > genjets){
  // if we are doing systematics but not applying JECs, call this function to load uncertainties independently of JECs
  for (size_t ijet(0); ijet < alljets_->size(); ijet++) {
    const pat::Jet &jet = (*alljets_)[ijet];
    //set JECs uncertainty values
    jecUncProvider->setJetEta(jet.eta());
    jecUncProvider->setJetPt(corrJet[ijet].pt());
    jecUnc[ijet] = jecUncProvider->getUncertainty(true);
    if (isData) {
      float rawFactor(jet.jecFactor("Uncorrected"));
      jetValues.setJetPt(jet.pt()*rawFactor);
      jetValues.setJetEta(jet.eta());
      jetValues.setJetA(jet.jetArea());
      jetValues.setRho(rhoEvent_);
      vector<float> corr_vals = jetCorrector->getSubCorrections(jetValues);
      jecUnc[ijet] = sqrt(pow(jecUnc[ijet],2) + pow((corr_vals.at(corr_vals.size()-1)/corr_vals.at(corr_vals.size()-2)-1.),2));
    }
    //set JER uncertainty values
    float smearedJetPt(0.);
    if (genJetPt[ijet]>0) smearedJetPt = genJetPt[ijet] + getJetResolutionSF(jet.eta())*(corrJet[ijet].pt() - genJetPt[ijet]);
    jerUnc[ijet] = smearedJetPt/corrJet[ijet].pt() - 1.;
    if (smearedJetPt < 0.01) jerUnc[ijet] = 0.; // so data will not have resolution unc. 
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

void jet_met_tools::getMETWithJEC(edm::Handle<pat::METCollection> mets, float &met, float &metPhi, unsigned isys){
  if(!doJEC && (!doSystematics || isys==kSysLast)) {
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
    
    if (isys == kSysJER) totCorr *= 1+jerUnc[ijet];
    else if (isys == kSysJECUp) totCorr *= 1+jecUnc[ijet];
    else if (isys == kSysJECDn) totCorr *= 1-jecUnc[ijet];

    if((rawJetP4.pt()*totCorr) <= 15.) continue;
    metx -= rawJetP4.px()*(totCorr - l1Corr);
    mety -= rawJetP4.py()*(totCorr - l1Corr);
  } // Loop over alljets_

  met = hypot(metx,mety);
  metPhi = atan2(mety,metx);

}

// the jetp4 and isBTaggged are technically redundant but avoid recalculating information
float jet_met_tools::jetBTagWeight(const pat::Jet &jet, const LVector &jetp4, bool isBTagged,
				   BTagEntry::OperatingPoint op,
				   const string &bc_full_syst, const string &udsg_full_syst) const{
  return jetBTagWeight(jet, jetp4, isBTagged, op, bc_full_syst, udsg_full_syst, "central", "central");				     
}

// the jetp4 and isBTaggged are technically redundant but avoid recalculating information
float jet_met_tools::jetBTagWeight(const pat::Jet &jet, const LVector &jetp4, bool isBTagged,
				   BTagEntry::OperatingPoint op,
				   const string &bc_full_syst, const string &udsg_full_syst,
				   const string &bc_fast_syst, const string &udsg_fast_syst) const{
  int hadronFlavour = abs(jet.hadronFlavour());
  BTagEntry::JetFlavor flav;
  string full_syst, fast_syst;
  switch(hadronFlavour){
  case 5: flav = BTagEntry::FLAV_B; break;
  case 4: flav = BTagEntry::FLAV_C; break;
  default: flav = BTagEntry::FLAV_UDSG; break;
  }

  double eff = getMCTagEfficiency(hadronFlavour, jetp4.pt(), jetp4.eta(), op);

  switch(flav){
  case BTagEntry::FLAV_B:
  case BTagEntry::FLAV_C:
    full_syst = bc_full_syst;
    fast_syst = bc_fast_syst;
    break;
  case BTagEntry::FLAV_UDSG:
    full_syst = udsg_full_syst;
    fast_syst = udsg_fast_syst;
    break;
  default:
    ERROR("Did not recognize BTagEntry::JetFlavor " << static_cast<int>(flav));
  }

  double sf = readers_full_.at(op)->eval_auto_bounds(full_syst, flav, jetp4.eta(), jetp4.pt());
  double sf_fs = isFastSim ? readers_fast_.at(op)->eval_auto_bounds(fast_syst, flav, jetp4.eta(), jetp4.pt()) : 1.;

  // procedure from https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
  double jet_scalefactor;
  if(isFastSim){
    double eff_fs = eff/sf_fs;
    jet_scalefactor = isBTagged ? sf*sf_fs : (1-sf*sf_fs*eff_fs)/(1-eff_fs);
  }else{
    jet_scalefactor = isBTagged ? sf : (1-sf*eff)/(1-eff);
  }
  return jet_scalefactor;
}

float jet_met_tools::getMCTagEfficiency(int pdgId, float pT, float eta, BTagEntry::OperatingPoint op) const{
  size_t rdr_idx = distance(op_pts_.cbegin(), find(op_pts_.cbegin(), op_pts_.cend(), op));
  pdgId = abs(pdgId);
  if(pdgId != 4 && pdgId != 5){
    // in the ghost clustering scheme to determine flavor, there are only b, c and other (id=0) flavors
    pdgId = 0;
  }
  int bin = btag_efficiencies_.at(rdr_idx)->FindFixBin(fabs(eta), pT, pdgId);
  float eff = btag_efficiencies_.at(rdr_idx)->GetBinContent(bin);
  return eff;
}

// get deltaR between two b-tagged CSVM jets
void jet_met_tools::fillDeltaRbb(vector<float> & deltaRbb, vector<float> & bb_pt,vector<float> & bb_m,vector<int> & bb_jet_idx1, vector<int> & bb_jet_idx2,vector<int> & bb_gs_idx,vector<int> & bb_gs_flavor, const vector<LVector> &jets, const vector<float> &jets_csv, const vector<bool> &jets_islep, const vector<float> &jets_pt, const vector<size_t> &branks, int & highcsv_index)
{
  for(size_t ijet(0); ijet<jets.size(); ijet++) {
    for(size_t jjet=ijet+1; jjet<jets.size(); jjet++) {
      if(jets_pt[ijet]<JetPtCut || jets_pt[jjet]<JetPtCut) continue;
      if(jets_csv[ijet]<CSVMedium || jets_csv[jjet]<CSVMedium) continue;
      if(jets_islep[ijet] || jets_islep[jjet]) continue;
      deltaRbb.push_back(deltaR(jets.at(ijet), jets.at(jjet)));
      bb_pt.push_back(sumPt(jets.at(ijet), jets.at(jjet)));
      bb_m.push_back(sumMass(jets.at(ijet), jets.at(jjet)));
      bb_jet_idx1.push_back(ijet);
      bb_jet_idx2.push_back(jjet);
      
      //Save index of pair with highest CSV
      if((ijet==branks.at(0) && jjet==branks.at(1)) || (ijet==branks.at(1) && jjet==branks.at(0))) highcsv_index = deltaRbb.size()-1;
      
      bb_gs_idx.push_back(-1); //Filled in writeMC
      bb_gs_flavor.push_back(0); //Filled in writeMC
    }
  }

}

// ranks jets by CSV, breaking ties by pt and mass, with jets less than 30 GeV or matched to leptons last
vector<size_t> jet_met_tools::getBRanking(const vector<LVector> &momentum, const vector<float> &csv,
                                          const vector<bool> &is_lep){
  typedef tuple<bool, float, float, float> Bness;
  vector<pair<Bness, long> > jets;
  for(size_t i = 0; i < momentum.size(); ++i){
    jets.push_back(pair<Bness, size_t>(Bness(!(is_lep.at(i)||momentum.at(i).pt()<30.0), csv.at(i), momentum.at(i).pt(), momentum.at(i).M()), -i));
  }
  sort(jets.rbegin(), jets.rend());
  vector<size_t> indices(jets.size());
  for(size_t i = 0; i < jets.size(); ++i){
    indices.at(i) = -jets.at(i).second;
  }
  return indices;
}

float jet_met_tools::getDeltaRbbHighCSV(const vector<LVector> &momentum, const vector<size_t> &brank, size_t nb){
  return (brank.size() > 1 && nb > 1) ? deltaR(momentum.at(brank.at(0)), momentum.at(brank.at(1))) : -1.;
}

float jet_met_tools::getDeltaRbbMax(const vector<LVector> &momentum, const vector<size_t> &brank,
                                    size_t nb){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    for(size_t j = i + 1; j < nb && j < brank.size(); ++j){
      float x = deltaR(momentum.at(brank.at(i)), momentum.at(brank.at(j)));
      if(x > max) max = x;
    }
  }
  return max;
}

float jet_met_tools::getDeltaRbbMin(const vector<LVector> &momentum, const vector<size_t> &brank,
                                    size_t nb){
  float min = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    for(size_t j = i + 1; j < nb && j < brank.size(); ++j){
      float x = deltaR(momentum.at(brank.at(i)), momentum.at(brank.at(j)));
      if(x < min || min < 0.) min = x;
    }
  }
  return min;
}

float jet_met_tools::getDeltaPhibb(const vector<LVector> &momentum, const vector<size_t> &brank, size_t nb){
  return (brank.size() > 1 && nb > 1) ? fabs(reco::deltaPhi(momentum.at(brank.at(0)).phi(), momentum.at(brank.at(1)).phi())) : -1.;
}

float jet_met_tools::getDeltaPhibbMax(const vector<LVector> &momentum, const vector<size_t> &brank,
                                      size_t nb){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    for(size_t j = i + 1; j < nb && j < brank.size(); ++j){
      float x = fabs(reco::deltaPhi(momentum.at(brank.at(i)).phi(), momentum.at(brank.at(j)).phi()));
      if(x > max) max = x;
    }
  }
  return max;
}

float jet_met_tools::getDeltaPhibbMin(const vector<LVector> &momentum, const vector<size_t> &brank,
                                      size_t nb){
  float min = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    for(size_t j = i + 1; j < nb && j < brank.size(); ++j){
      float x = fabs(reco::deltaPhi(momentum.at(brank.at(i)).phi(), momentum.at(brank.at(j)).phi()));
      if(x < min || min < 0.) min = x;
    }
  }
  return min;
}

float jet_met_tools::getMbb(const vector<LVector> &momentum, const vector<size_t> &brank, size_t nb){
  return (brank.size() > 1 && nb > 1) ? sumMass(momentum.at(brank.at(0)), momentum.at(brank.at(1))) : -1.;
}

float jet_met_tools::getMbbMax(const vector<LVector> &momentum, const vector<size_t> &brank,
                               size_t nb){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    for(size_t j = i + 1; j < nb && j < brank.size(); ++j){
      float x = sumMass(momentum.at(brank.at(i)), momentum.at(brank.at(j)));
      if(x > max) max = x;
    }
  }
  return max;
}

float jet_met_tools::getMbbMin(const vector<LVector> &momentum, const vector<size_t> &brank,
                               size_t nb){
  float min = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    for(size_t j = i + 1; j < nb && j < brank.size(); ++j){
      float x = sumMass(momentum.at(brank.at(i)), momentum.at(brank.at(j)));
      if(x < min || min < 0.) min = x;
    }
  }
  return min;
}

float jet_met_tools::getMblepMax2(const vector<LVector> &jets, const vector<size_t> &brank,
                                  size_t nb, const LVector &lep){
  switch(nb){
  case 0: return -1.;
  case 1: return sumMass(jets.at(brank.at(0)), lep);
  default:
    float a = sumMass(jets.at(brank.at(0)), lep);
    float b = sumMass(jets.at(brank.at(1)), lep);
    return max(a,b);
  }
}

float jet_met_tools::getMblepMin2(const vector<LVector> &jets, const vector<size_t> &brank,
                                  size_t nb, const LVector &lep){
  switch(nb){
  case 0: return -1.;
  case 1: return sumMass(jets.at(brank.at(0)), lep);
  default:
    float a = sumMass(jets.at(brank.at(0)), lep);
    float b = sumMass(jets.at(brank.at(1)), lep);
    return min(a,b);
  }
}

float jet_met_tools::getMblepMax(const vector<LVector> &jets, const vector<size_t> &brank,
                                 size_t nb, const LVector &lep){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = sumMass(jets.at(brank.at(i)), lep);
    if(m>max) max = m;
  }
  return max;
}

float jet_met_tools::getMblepMin(const vector<LVector> &jets, const vector<size_t> &brank,
                                 size_t nb, const LVector &lep){
  float min = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = sumMass(jets.at(brank.at(i)), lep);
    if(m<min || min<0.) min = m;
  }
  return min;
}

float jet_met_tools::getDeltaRblepMax2(const vector<LVector> &jets, const vector<size_t> &brank,
                                       size_t nb, const LVector &lep){
  switch(nb){
  case 0: return -1.;
  case 1: return deltaR(jets.at(brank.at(0)), lep);
  default:
    float a = deltaR(jets.at(brank.at(0)), lep);
    float b = deltaR(jets.at(brank.at(1)), lep);
    return max(a,b);
  }
}

float jet_met_tools::getDeltaRblepMin2(const vector<LVector> &jets, const vector<size_t> &brank,
                                       size_t nb, const LVector &lep){
  switch(nb){
  case 0: return -1.;
  case 1: return deltaR(jets.at(brank.at(0)), lep);
  default:
    float a = deltaR(jets.at(brank.at(0)), lep);
    float b = deltaR(jets.at(brank.at(1)), lep);
    return min(a,b);
  }
}

float jet_met_tools::getDeltaRblepMax(const vector<LVector> &jets, const vector<size_t> &brank,
                                      size_t nb, const LVector &lep){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = deltaR(jets.at(brank.at(i)), lep);
    if(m>max) max = m;
  }
  return max;
}

float jet_met_tools::getDeltaRblepMin(const vector<LVector> &jets, const vector<size_t> &brank,
                                      size_t nb, const LVector &lep){
  float min = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = deltaR(jets.at(brank.at(i)), lep);
    if(m<min || min < 0.) min = m;
  }
  return min;
}

float jet_met_tools::getDeltaPhiblepMax2(const vector<LVector> &jets, const vector<size_t> &brank,
                                         size_t nb, const LVector &lep){
  switch(nb){
  case 0: return -1.;
  case 1: return fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), lep.phi()));
  default:
    float a = fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), lep.phi()));
    float b = fabs(reco::deltaPhi(jets.at(brank.at(1)).phi(), lep.phi()));
    return max(a,b);
  }
}

float jet_met_tools::getDeltaPhiblepMin2(const vector<LVector> &jets, const vector<size_t> &brank,
                                         size_t nb, const LVector &lep){
  switch(nb){
  case 0: return -1.;
  case 1: return fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), lep.phi()));
  default:
    float a = fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), lep.phi()));
    float b = fabs(reco::deltaPhi(jets.at(brank.at(1)).phi(), lep.phi()));
    return min(a,b);
  }
}

float jet_met_tools::getDeltaPhiblepMax(const vector<LVector> &jets, const vector<size_t> &brank,
                                        size_t nb, const LVector &lep){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = fabs(reco::deltaPhi(jets.at(brank.at(i)).phi(), lep.phi()));
    if(m>max) max = m;
  }
  return max;
}

float jet_met_tools::getDeltaPhiblepMin(const vector<LVector> &jets, const vector<size_t> &brank,
                                        size_t nb, const LVector &lep){
  float min = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = fabs(reco::deltaPhi(jets.at(brank.at(i)).phi(), lep.phi()));
    if(m<min || min < 0.) min = m;
  }
  return min;
}

float jet_met_tools::getMTbmetMax2(const vector<LVector> &jets, const vector<size_t> &brank,
                                   size_t nb, float met, float met_phi){
  if(nb == 0){
    return -1.;
  }else if(nb == 1){
    const auto &j = jets.at(brank.at(0));
    return getMT(j.M(), j.pt(), j.phi(), 0., met, met_phi);
  }else{
    const auto &ja = jets.at(brank.at(0));
    const auto &jb = jets.at(brank.at(1));
    float a = getMT(ja.M(), ja.pt(), ja.phi(), 0., met, met_phi);
    float b = getMT(jb.M(), jb.pt(), jb.phi(), 0., met, met_phi);
    return max(a,b);
  }
}

float jet_met_tools::getMTbmetMin2(const vector<LVector> &jets, const vector<size_t> &brank,
                                   size_t nb, float met, float met_phi){
  if(nb == 0){
    return -1.;
  }else if(nb == 1){
    const auto &j = jets.at(brank.at(0));
    return getMT(j.M(), j.pt(), j.phi(), 0., met, met_phi);
  }else{
    const auto &ja = jets.at(brank.at(0));
    const auto &jb = jets.at(brank.at(1));
    float a = getMT(ja.M(), ja.pt(), ja.phi(), 0., met, met_phi);
    float b = getMT(jb.M(), jb.pt(), jb.phi(), 0., met, met_phi);
    return min(a,b);
  }
}

float jet_met_tools::getMTbmetMax(const vector<LVector> &jets, const vector<size_t> &brank,
                                  size_t nb, float met, float met_phi){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    const auto &jet = jets.at(brank.at(i));
    float m = getMT(jet.M(), jet.pt(), jet.phi(), 0., met, met_phi);
    if(m>max) m = max;
  }
  return max;
}

float jet_met_tools::getMTbmetMin(const vector<LVector> &jets, const vector<size_t> &brank,
                                  size_t nb, float met, float met_phi){
  float min = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    const auto &jet = jets.at(brank.at(i));
    float m = getMT(jet.M(), jet.pt(), jet.phi(), 0., met, met_phi);
    if(m<min || min<0.) m = min;
  }
  return min;
}

float jet_met_tools::getDeltaPhibmetMax2(const vector<LVector> &jets, const vector<size_t> &brank,
                                         size_t nb, float met_phi){
  switch(nb){
  case 0: return -1.;
  case 1: return fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), met_phi));
  default:
    float a = fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), met_phi));
    float b = fabs(reco::deltaPhi(jets.at(brank.at(1)).phi(), met_phi));
    return max(a,b);
  }
}

float jet_met_tools::getDeltaPhibmetMin2(const vector<LVector> &jets, const vector<size_t> &brank,
                                         size_t nb, float met_phi){
  switch(nb){
  case 0: return -1.;
  case 1: return fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), met_phi));
  default:
    float a = fabs(reco::deltaPhi(jets.at(brank.at(0)).phi(), met_phi));
    float b = fabs(reco::deltaPhi(jets.at(brank.at(1)).phi(), met_phi));
    return min(a,b);
  }
}

float jet_met_tools::getDeltaPhibmetMax(const vector<LVector> &jets, const vector<size_t> &brank,
                                        size_t nb, float met_phi){
  float max = -1.;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = fabs(reco::deltaPhi(jets.at(brank.at(i)).phi(), met_phi));
    if(m>max) max = m;
  }
  return max;
}

float jet_met_tools::getDeltaPhibmetMin(const vector<LVector> &jets, const vector<size_t> &brank,
                                        size_t nb, float met_phi){
  float min = -1;
  for(size_t i = 0; i < nb && i < brank.size(); ++i){
    float m = fabs(reco::deltaPhi(jets.at(brank.at(i)).phi(), met_phi));
    if(m<min || min < 0.) min = m;
  }
  return min;
}

void jet_met_tools::clusterFatJets(int &nfjets, float &mj,
                                   std::vector<float> &fjets_pt, 
                                   std::vector<float> &fjets_eta,
                                   std::vector<float> &fjets_phi, 
                                   std::vector<float> &fjets_m,
                                   std::vector<int> &fjets_nconst,
                                   std::vector<int> &jets_fjet_index,
                                   baby_full &baby,
                                   double radius,
                                   double min_jets_pt,
                                   bool cluster_leps){

  vector<fastjet::PseudoJet> sjets(0);
  for(size_t ijet(0); ijet < baby.jets_pt().size(); ijet++){
    if (!cluster_leps && baby.jets_islep()[ijet]) continue;
    if (baby.jets_pt()[ijet]>min_jets_pt || (cluster_leps && baby.jets_islep()[ijet])) {
      TLorentzVector vjet; vjet.SetPtEtaPhiM(baby.jets_pt()[ijet], baby.jets_eta()[ijet], baby.jets_phi()[ijet], baby.jets_m()[ijet]);
      const fastjet::PseudoJet this_pj(vjet.Px(), vjet.Py(), vjet.Pz(), vjet.E());
      sjets.push_back(this_pj);
    }
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
  jets_fjet_index.resize(baby.jets_csv().size(), -1);

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
    for(size_t ijet = 0; ijet < baby.jets_pt().size(); ++ijet){
      for(size_t cjet = 0; cjet < cjets.size(); ++ cjet){
        if(fabs(cjets.at(cjet).pt() - baby.jets_pt()[ijet]) < 0.0001){
          jets_fjet_index.at(ijet) = ipj;
        }
      } // Loop over fat jet constituents
    } // Loop over skinny jets
  } // Loop over fat jets

}

double jet_met_tools::getSysMJ(double radius, vector<LVector> &jets, vector<bool> &jets_islep, double min_jets_pt, bool cluster_leps){

  double mj(0.);
  vector<fastjet::PseudoJet> sjets(0);
  for(size_t ijet(0); ijet < jets.size(); ijet++){
    if (!cluster_leps && jets_islep[ijet]) continue;
    if (jets[ijet].pt()>min_jets_pt || (cluster_leps && jets_islep[ijet])) {
      const fastjet::PseudoJet this_pj(jets[ijet].px(), jets[ijet].py(), jets[ijet].pz(), jets[ijet].energy());
      sjets.push_back(this_pj);
    }
  } // Loop over skinny jets
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, radius);
  fastjet::ClusterSequence cs(sjets, jet_def);
  vector<fastjet::PseudoJet> fjets = cs.inclusive_jets();
  sort(fjets.begin(), fjets.end(), greaterM);
  for(size_t ipj = 0; ipj < fjets.size(); ++ipj) mj += fjets.at(ipj).m();

  return mj;
}
jet_met_tools::jet_met_tools(TString ijecName, bool doSys, bool fastSim):
  jecName(ijecName),
  doSystematics(doSys),
  isFastSim(fastSim),
  calib_full_(),
  calib_fast_(),
  readers_full_(),
  readers_fast_(),
  btag_efficiencies_(op_pts_.size()){
  if (jecName.Contains("DATA")) isData = true;
  else if (jecName.Contains("MC")) isData = false;
  else cout<<endl<<"BABYMAKER: jet_met_tools: The jecLabel string must contain either 'DATA' or 'MC'. "
                 <<"Currently running with "<<jecName<<endl<<endl;

  doJEC = jecName.Contains("onthefly"); 
  jecName.ReplaceAll("onthefly_","");
  string basename(getenv("CMSSW_BASE"));
  basename += "/src/babymaker/data/jec/";
  basename += jecName.Data();

  if(doJEC || doSystematics) {
    vector<JetCorrectorParameters> jecFiles;
    if (doJEC) cout<<endl<<"BABYMAKER: jet_met_tools: Applying JECs on-the-fly with files "
      <<basename.c_str()<<"*.txt"<<endl<<endl;
    //these files need to be loaded even if doing only systematics
    jecFiles.push_back(JetCorrectorParameters(basename+"_L1FastJet_AK4PFchs.txt"));
    jecFiles.push_back(JetCorrectorParameters(basename+"_L2Relative_AK4PFchs.txt"));
    jecFiles.push_back(JetCorrectorParameters(basename+"_L3Absolute_AK4PFchs.txt"));
    if(jecName.Contains("DATA")) 
      jecFiles.push_back(JetCorrectorParameters(basename+"_L2L3Residual_AK4PFchs.txt"));
    jetCorrector.reset(new FactorizedJetCorrectorCalculator(jecFiles));
  }
  if (doSystematics){
    cout<<endl<<"BABYMAKER: jet_met_tools: Using JEC uncertainties from file "
              <<basename.c_str()<<"_Uncertainty_AK4PFchs.txt"<<endl<<endl;
    jecUncProvider.reset(new JetCorrectionUncertainty(basename+"_Uncertainty_AK4PFchs.txt"));
  }

  // set btag working points
  TString cmssw(getenv("CMSSW_VERSION"));
  if(cmssw.Contains("CMSSW_7_4")){
    CSVLoose  = 0.605;
    CSVMedium = 0.890;
    CSVTight  = 0.970;
  }
  else if(cmssw.Contains("CMSSW_8_0")){
    CSVLoose  = 0.460;
    CSVMedium = 0.800;
    CSVTight  = 0.935;
  }

  // only add b-tagging weights if requested
  string scaleFactorFile(getenv("CMSSW_BASE"));
  scaleFactorFile+="/src/babymaker/bmaker/data/CSVv2.csv";
  calib_full_.reset(new BTagCalibration("csvv2", scaleFactorFile));
  for(const auto &op: op_pts_){
    readers_full_[op] = MakeUnique<BTagCalibrationReader>(op, "central", vector<string>{"up", "down"});
    readers_full_.at(op)->load(*calib_full_, BTagEntry::FLAV_UDSG, "incl");
    readers_full_.at(op)->load(*calib_full_, BTagEntry::FLAV_C, "comb");
    readers_full_.at(op)->load(*calib_full_, BTagEntry::FLAV_B, "comb");
  }

  if (isFastSim){
    string scaleFactorFileFastSim(getenv("CMSSW_BASE"));
    scaleFactorFileFastSim+="/src/babymaker/bmaker/data/CSV_13TEV_Combined_14_7_2016.csv";
    calib_fast_.reset(new BTagCalibration("csvv2", scaleFactorFileFastSim));
    for(const auto &op: op_pts_){
      readers_fast_[op] = MakeUnique<BTagCalibrationReader>(op, "central", vector<string>{"up", "down"});
      readers_fast_.at(op)->load(*calib_fast_, BTagEntry::FLAV_UDSG, "fastsim");
      readers_fast_.at(op)->load(*calib_fast_, BTagEntry::FLAV_C, "fastsim");
      readers_fast_.at(op)->load(*calib_fast_, BTagEntry::FLAV_B, "fastsim");
    }
  }

  string filename(getenv("CMSSW_BASE"));
  filename+="/src/babymaker/bmaker/data/btagEfficiency.root";
  TFile *efficiencyFile = TFile::Open(filename.c_str());
  if(efficiencyFile->IsOpen()) {
    for(size_t i = 0; i < op_pts_.size(); ++i){
      string hist_name;
      switch(op_pts_.at(i)){
      case BTagEntry::OP_LOOSE: hist_name = "btagEfficiency_loose"; break;
      case BTagEntry::OP_MEDIUM: hist_name = "btagEfficiency_medium"; break;
      case BTagEntry::OP_TIGHT: hist_name = "btagEfficiency_tight"; break;
      case BTagEntry::OP_RESHAPING: hist_name = "btagEfficiency_reshaping"; break;
      default: hist_name = "btagEfficiency"; break;
      }
      btag_efficiencies_.at(i) = static_cast<const TH3D*>(efficiencyFile->Get(hist_name.c_str()));
      if(!btag_efficiencies_.at(i)){
	ERROR("Could not find efficiency parametrization " << hist_name);
      }
    }
  }else{
    ERROR("Could not find efficiency file " << filename);
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
