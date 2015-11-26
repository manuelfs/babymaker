// jet_met_tools: Functions related to jets, MET, and JECs

// System include files
#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <DataFormats/Math/interface/deltaR.h>
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
				   btagVariation readerTypeBC, btagVariation readerTypeUDSG,
				   btagVariation readerTypeBC_fs, btagVariation readerTypeUDSG_fs) {
  double jet_scalefactor = 1.0;
  int hadronFlavour = abs(jet.hadronFlavour());
  // only apply weights if readers are initialized
  if(readersBC.size()>0 && readersUDSG.size()>0) {
    double jetpttemp = jetp4.pt();
    // maximum pt in the parameterizations is 670 GeV
    if(jetpttemp>670) jetpttemp=669.99;
    try {
      double eff = getMCTagEfficiency(abs(hadronFlavour), jetpttemp, fabs(jetp4.eta()));
      double sf(1.), sf_fs(1.); 
      // procedure from https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
      switch ( abs(hadronFlavour) ) {
      case 5:
        sf = readersBC.at(readerTypeBC)->eval(BTagEntry::FLAV_B, jetp4.eta(), jetpttemp);
        if (isFastSim) sf_fs = readersBC_fs.at(readerTypeBC_fs)->eval(BTagEntry::FLAV_B, jetp4.eta(), jetpttemp);
        break;
      case 4:
        sf = readersBC.at(readerTypeBC)->eval(BTagEntry::FLAV_C, jetp4.eta(), jetpttemp);
        if (isFastSim) sf_fs = readersBC_fs.at(readerTypeBC_fs)->eval(BTagEntry::FLAV_C, jetp4.eta(), jetpttemp);
        break;
      default:
        sf = readersUDSG.at(readerTypeUDSG)->eval(BTagEntry::FLAV_UDSG, jetp4.eta(), jetpttemp);
        if (isFastSim) sf_fs = readersUDSG_fs.at(readerTypeUDSG_fs)->eval(BTagEntry::FLAV_UDSG, jetp4.eta(), jetpttemp);
        break;
      }
      if (isFastSim) {
        double eff_fs = eff/sf_fs;
        jet_scalefactor = isBTagged ? sf*sf_fs : (1-sf*sf_fs*eff_fs)/(1-eff_fs);
      } else {
        jet_scalefactor = isBTagged ? sf : (1-sf*eff)/(1-eff);
      }
    }
    catch (std::exception &e) {
      std::cout << "Caught exception: " << e.what() << " for a jet with (pt,eta,hadron flavor)=(" 
                << jetp4.pt() << "," << jetp4.eta() << "," << hadronFlavour << ")" << std::endl;
    }
  }
  return jet_scalefactor;
}

float jet_met_tools::getMCTagEfficiency(int pdgId, float pT, float eta)
{
  // for testing purposes, just use a plausible efficiency
  // need to add code for getting efficiencies
  if(abs(pdgId)==4 || abs(pdgId)==5) {
    int bin = btagEfficiencyParameterization->FindBin(eta, pT, pdgId);
    return btagEfficiencyParameterization->GetBinContent(bin);
  }
  else {
    // in the ghost clustering scheme to determine flavor, there are only b, c and other (id=0) flavors
    int bin = btagEfficiencyParameterization->FindBin(eta, pT, 0);
    return btagEfficiencyParameterization->GetBinContent(bin);
  }
}

void jet_met_tools::clusterFatJets(int &nfjets, float &mj,
                                  vector<float> &fjets_pt, 
                                  vector<float> &fjets_eta,
                                  vector<float> &fjets_phi, 
                                  vector<float> &fjets_m,
                                  vector<int> &fjets_nconst,
                                  vector<float> &fjets_sumcsv,
                                  vector<float> &fjets_poscsv,
                                  vector<int> &fjets_btags,
                                  vector<int> &jets_fjet_index,
                                  double radius,
                                  vector<LVector> &jets,
                                  vector<float> &jets_csv){

  vector<fastjet::PseudoJet> sjets(0);
  vector<float> csvs(0);
  for(size_t ijet(0); ijet < jets.size(); ijet++){
    const fastjet::PseudoJet this_pj(jets[ijet].px(), jets[ijet].py(), jets[ijet].pz(), jets[ijet].energy());
    sjets.push_back(this_pj);
    csvs.push_back(jets_csv[ijet]); // This require to have the same jets in baby and jets
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
  fjets_sumcsv.resize(fjets.size());
  fjets_poscsv.resize(fjets.size());
  fjets_btags.resize(fjets.size());
  jets_fjet_index.resize(jets.size());

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
    fjets_btags.at(ipj) = 0;
    fjets_sumcsv.at(ipj) = 0.;
    fjets_poscsv.at(ipj) = 0.;
    for(size_t ijet = 0; ijet < jets.size(); ++ijet){
      for(size_t cjet = 0; cjet < cjets.size(); ++ cjet){
        if((cjets.at(cjet) - sjets.at(ijet)).pt() < 0.0001){
          jets_fjet_index.at(ijet) = ipj;
          fjets_sumcsv.at(ipj) += csvs.at(ijet);
          if(csvs.at(ijet) > 0.){
            fjets_poscsv.at(ipj) += csvs.at(ijet);
          }
          if(csvs.at(ijet) > CSVMedium){
            ++(fjets_btags.at(ipj));
          }
        }
      } // Loop over fat jet constituents
    } // Loop over skinny jets
  } // Loop over fat jets

}

double jet_met_tools::getSysMJ(double radius, vector<LVector> &jets){

  double mj(0.);
  vector<fastjet::PseudoJet> sjets(0);
  for(size_t ijet(0); ijet < jets.size(); ijet++){
    const fastjet::PseudoJet this_pj(jets[ijet].px(), jets[ijet].py(), jets[ijet].pz(), jets[ijet].energy());
    sjets.push_back(this_pj);
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
  calib(0),
  calibFS(0){

  if (jecName.Contains("DATA")) isData = true;
  else if (jecName.Contains("MC")) isData = false;
  else cout<<endl<<"BABYMAKER: jet_met_tools: The jecLabel string must contain either 'DATA' or 'MC'. Currently running with "<<jecName<<endl<<endl;

  doJEC = !jecName.Contains("miniAOD");
  jecName.ReplaceAll("miniAOD_","");
  string basename(getenv("CMSSW_BASE"));
  basename += "/src/babymaker/data/jec/";
  basename += jecName.Data();

  if(doJEC || doSystematics) {
    vector<JetCorrectorParameters> jecFiles;
    if (doJEC) cout<<endl<<"BABYMAKER: jet_met_tools: Applying JECs on-the-fly with files "<<basename.c_str()<<"*.txt"<<endl<<endl;
    else       cout<<endl<<"BABYMAKER: jet_met_tools: Reading JECs needed to determine "
                         <<"residuals uncertainty with files "<<basename.c_str()<<"*.txt"<<endl<<endl;
    jecFiles.push_back(JetCorrectorParameters(basename+"_L1FastJet_AK4PFchs.txt"));
    jecFiles.push_back(JetCorrectorParameters(basename+"_L2Relative_AK4PFchs.txt"));
    jecFiles.push_back(JetCorrectorParameters(basename+"_L3Absolute_AK4PFchs.txt"));
    if(jecName.Contains("DATA")) jecFiles.push_back(JetCorrectorParameters(basename+"_L2L3Residual_AK4PFchs.txt"));
    jetCorrector = new FactorizedJetCorrectorCalculator(jecFiles);
  }
  if (doSystematics){
    cout<<endl<<"BABYMAKER: jet_met_tools: Using JEC uncertainties from file "<<basename.c_str()<<"_Uncertainty_AK4PFchs.txt"<<endl<<endl;
    jecUncProvider = new JetCorrectionUncertainty(basename+"_Uncertainty_AK4PFchs.txt");
  }

  // only add b-tagging weights if requested

  std::string scaleFactorFile(getenv("CMSSW_BASE"));
  scaleFactorFile+="/src/babymaker/bmaker/data/CSVv2.csv";
  calib   = new BTagCalibration("csvv1", scaleFactorFile);
  std::vector<std::string> variationTypes = {"central", "up", "down"};
  for(auto itype : variationTypes) {
    // BC full sim
    readersBC.push_back(new BTagCalibrationReader(calib, // calibration instance 
					    BTagEntry::OP_MEDIUM, // operating point
					    "mujets", // measurement type ("comb" or "mujets")
					    itype)); // systematics type ("central", "up", or "down")
    // UDSG full sim
    readersUDSG.push_back(new BTagCalibrationReader(calib, // calibration instance 
					      BTagEntry::OP_MEDIUM, // operating point
					      "comb", // measurement type ("comb" or "mujets")
					      itype)); // systematics type ("central", "up", or "down")
  }
  if (isFastSim){
    std::string scaleFactorFileFastSim(getenv("CMSSW_BASE"));
    scaleFactorFileFastSim+="/src/babymaker/bmaker/data/CSV_13TEV_Combined_20_11_2015.csv";
    calibFS = new BTagCalibration("csvv1", scaleFactorFileFastSim);
    for(auto itype : variationTypes) {
      // BC fastsim
      readersBC_fs.push_back(new BTagCalibrationReader(calibFS, // calibration instance 
  					    BTagEntry::OP_MEDIUM, // operating point
  					    "fastsim", // measurement type ("comb" or "mujets")
  					    itype)); // systematics type ("central", "up", or "down")
      // UDSG fastsim
      readersUDSG_fs.push_back(new BTagCalibrationReader(calibFS, // calibration instance 
  					      BTagEntry::OP_MEDIUM, // operating point
  					      "fastsim", // measurement type ("comb" or "mujets")
  					      itype)); // systematics type ("central", "up", or "down")
    }
  }

  std::string filename(getenv("CMSSW_BASE"));
  filename+="/src/babymaker/bmaker/data/btagEfficiency.root";
  TFile *efficiencyFile = TFile::Open(filename.c_str());
  if(efficiencyFile->IsOpen()) {
    btagEfficiencyParameterization = static_cast<TH3F*>(efficiencyFile->Get("btagEfficiency"));
  }
  else {
    throw cms::Exception("FileNotFound") 
      << "Could not find efficiency file " << filename << "." << std::endl;
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
  if(calib !=0 ) delete calib;
  if(calibFS !=0 ) delete calibFS;
  for(auto ireader : readersBC) delete ireader;
  for(auto ireader : readersUDSG) delete ireader;
}

