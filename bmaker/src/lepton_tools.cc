//// LEPTON_TOOLS: Lepton selection and isolation
//// Function names follow the first-lowercase, following words-uppercase. No underscores

// System include files
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <numeric>

//ROOT include files
#include "TFile.h"
#include "TGraphAsymmErrors.h"

// FW include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

// User include files
#include "babymaker/bmaker/interface/lepton_tools.hh"

using namespace std;
using namespace utilities;

namespace{
  template<typename T>
    T GetSF(const string &file_name, const string &item_name){
    string path = string(getenv("CMSSW_BASE"))+"/src/babymaker/bmaker/data/lepton_sf/"+file_name;
    TFile f(path.c_str(), "read");
    if(!f.IsOpen()) ERROR("Could not open " << file_name);
    T* item = static_cast<T*>(f.Get(item_name.c_str()));
    if(!item) ERROR("Could not find " << item_name << " in " << file_name);
    return *item;
  }

  template<typename T>
    pair<double, double> GetSF(const T &h, double x, double y, bool ignore_error = false){
    pair<double, double> sf;
    auto bin = h.FindFixBin(x, y);
    if((h.IsBinOverflow(bin) || h.IsBinUnderflow(bin))
       && h.GetBinContent(bin) == 0. && h.GetBinError(bin) == 0.){
      auto bin_x = h.GetXaxis()->FindFixBin(x);
      auto bin_y = h.GetYaxis()->FindFixBin(y);
      if(bin_x <= 0) bin_x = 1;
      if(bin_x > h.GetNbinsX()) bin_x = h.GetNbinsX();
      if(bin_y <= 0) bin_y = 1;
      if(bin_y > h.GetNbinsY()) bin_y = h.GetNbinsY();
      sf = {h.GetBinContent(bin_x, bin_y), h.GetBinError(bin_x, bin_y)};
    }else{
      sf = {h.GetBinContent(bin), h.GetBinError(bin)};
    }
    if(ignore_error) sf.second = 0.;
    return sf;
  }

  pair<double, double> MergeSF(pair<double, double> a,
                               pair<double, double> b){
    double sf = a.first * b.first;
    double err = hypot(a.first*b.second, b.first*a.second);
    return {sf, err};
  }

  TH2D GraphToHist(const TGraphAsymmErrors &g){
    struct Point{
      double xl, xh, y, e;
      Point(double xl_in, double xh_in, double y_in, double e_in):
        xl(xl_in),
        xh(xh_in),
        y(y_in),
        e(e_in){
      }
      bool operator<(const Point &p) const{
        return make_tuple(xl, xh, fabs(log(fabs(y))), fabs(e))
          <make_tuple(p.xl, p.xh, fabs(log(fabs(p.y))), fabs(p.e));
      }
    };
    vector<Point> bins;
    Double_t *x = g.GetX();
    Double_t *xl = g.GetEXlow();
    Double_t *xh = g.GetEXhigh();
    Double_t *y = g.GetY();
    Double_t *yl = g.GetEYlow();
    Double_t *yh = g.GetEYhigh();
    for(int i = 0; i < g.GetN(); ++i){
      bins.emplace_back(x[i]-fabs(xl[i]), x[i]+fabs(xh[i]),
                        y[i], max(fabs(yl[i]), fabs(yh[i])));
    }
    bool problems = true;
    while(problems){
      stable_sort(bins.begin(), bins.end());
      problems = false;
      for(auto low = bins.begin(); !problems && low != bins.end(); ++low){
        auto high = low;
        ++high;
        if(high == bins.end()) break;
        double new_y = sqrt(low->y * high->y);
        double top = max(low->y+low->e, high->y+high->e);
        double bot = min(low->y-low->e, high->y-high->e);
        double new_e = max(top-new_y, new_y-bot);
        if(low->xh < high->xl){
          //Gap
          bins.insert(high, Point(low->xh, high->xl, new_y, new_e));
        }else if(low->xh > high->xl){
          //Overlap
          problems = true;
          if(low->xh < high->xh){
            //Plain overlap
            Point new_low(low->xl, high->xl, low->y, low->e);
            Point new_mid(high->xl, low->xh, new_y, new_e);
            Point new_high(low->xh, high->xh, high->y, high->e);
            *low = new_low;
            *high = new_high;
            bins.insert(high, new_mid);
          }else if(low->xh == high->xh){
            //Subset -> 2 bins
            Point new_low(low->xl, high->xl, low->y, low->e);
            Point new_high(high->xl, high->xh, new_y, new_e);
            *low = new_low;
            *high = new_high;
          }else{
            //Subset -> 3 bins
            Point new_low(low->xl, high->xl, low->y, low->e);
            Point new_mid(high->xl, high->xh, new_y, new_e);
            Point new_high(high->xh, low->xh, low->y, low->e);
            *low = new_low;
            *high = new_high;
            bins.insert(high, new_mid);
          }
        }
      }
    }
    vector<double> bin_edges(bins.size()+1);
    for(size_t i = 0; i < bins.size(); ++i){
      bin_edges.at(i) = bins.at(i).xl;
    }
    bin_edges.back() = bins.back().xh;
    TH2D h(g.GetName(), (string(g.GetTitle())+";"+g.GetXaxis()->GetTitle()+";"+g.GetYaxis()->GetTitle()).c_str(),
           1, 0., 1.e4, bin_edges.size()-1, &bin_edges.at(0));
    for(int ix = 0; ix <= 2; ++ix){
      h.SetBinContent(ix, 0, 1.);
      h.SetBinError(ix, 0, 1.);
      h.SetBinContent(ix, h.GetNbinsY()+1, 1.);
      h.SetBinError(ix, h.GetNbinsY()+1, 1.);
      for(int iy = 1; iy <= h.GetNbinsY(); ++iy){
        h.SetBinContent(ix, iy, bins.at(iy-1).y);
        h.SetBinError(ix ,iy, bins.at(iy-1).e);
      }
    }
    return h;
  }
}

//////////////////// Scale Factor loading
const TH2F lepton_tools::sf_full_muon_medium = GetSF<TH2F>("sf_full_muon_medium.root",
                                                           "pt_abseta_PLOT_pair_probeMultiplicity_bin0");
const TH2F lepton_tools::sf_full_muon_iso = GetSF<TH2F>("sf_full_muon_iso.root",
                                                        "pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass");
const TH2D lepton_tools::sf_full_muon_tracking = GraphToHist(GetSF<TGraphAsymmErrors>("sf_full_muon_tracking.root",
                                                                                      "ratio_eta"));
const TH2F lepton_tools::sf_full_electron_medium = GetSF<TH2F>("sf_full_electron_idiso.root",
                                                               "GsfElectronToMedium");
const TH2F lepton_tools::sf_full_electron_iso = GetSF<TH2F>("sf_full_electron_idiso.root",
                                                            "MVAVLooseElectronToMini");
const TH2F lepton_tools::sf_full_electron_tracking = GetSF<TH2F>("sf_full_electron_tracking.root",
                                                                 "EGamma_SF2D");
const TH2D lepton_tools::sf_fast_muon_medium = GetSF<TH2D>("sf_fast_muon_medium.root",
                                                           "histo2D");
const TH2D lepton_tools::sf_fast_muon_iso = GetSF<TH2D>("sf_fast_muon_iso.root",
                                                        "histo2D");
const TH2D lepton_tools::sf_fast_electron_mediumiso = GetSF<TH2D>("sf_fast_electron_mediumiso.root",
                                                                  "histo2D");

//////////////////// Muons
bool lepton_tools::isSignalMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= SignalLeptonPtCut) return false;
  if(fabs(lep.eta()) > MuonEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idMuon(lep, vtx, kMedium)) return false;
  // Isolation cuts
  if(lepIso >= 0 && lepIso > MuonMiniIsoCut) return false;

  return true;
}

bool lepton_tools::isVetoMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= VetoLeptonPtCut) return false;
  if(fabs(lep.eta()) > MuonEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idMuon(lep, vtx, kMedium)) return false; // Using RA2/b muons
  // Isolation cuts
  if(lepIso >= 0 && lepIso > MuonMiniIsoCut) return false;

  return true;
}

bool lepton_tools::idMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold) {
  double dz(0.), d0(0.);
  if(!vertexMuon(lep, vtx, dz, d0)) return false;

  bool good_global;
  switch(threshold){
  default:
  case kVeto:
  case kLoose:
    return lep.isLooseMuon();
  case kMedium:
    return lep.isMediumMuon();
  case kMediumICHEP:
    //From https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Short_Term_Medium_Muon_Definitio
    good_global = lep.isGlobalMuon()
      && lep.globalTrack()->normalizedChi2() < 3.
      && lep.combinedQuality().chi2LocalPosition < 12.
      && lep.combinedQuality().trkKink < 20.;
    return  muon::isLooseMuon(lep)
      && lep.innerTrack()->validFraction() > 0.49
      && muon::segmentCompatibility(lep) > (good_global ? 0.303 : 0.451);
  case kTight:
    return lep.isTightMuon(vtx->at(0));
  }
}

bool lepton_tools::vertexMuon(const pat::Muon &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
  dz = 0.; d0 = 0.;
  if(lep.track().isAvailable()){ // If the track is not available we probably don't want the muon
    dz = lep.track()->vz()-vtx->at(0).z();
    d0 = lep.track()->d0()-vtx->at(0).x()*sin(lep.track()->phi())+vtx->at(0).y()*cos(lep.track()->phi());
  }
  if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

  return true;
}

double lepton_tools::getEffAreaMuon(double eta){
  double abseta = fabs(eta);
  if (abseta < 0.8) return 0.0735;
  else if (abseta < 1.3) return 0.0619;
  else if (abseta < 2.0) return 0.0465;
  else if (abseta < 2.2) return 0.0433;
  else if (abseta < 2.5) return 0.0577;
  else return 0;
}

double lepton_tools::getRelIsolation(const pat::Muon &lep, double rho){
  double ch_iso(lep.pfIsolationR04().sumChargedHadronPt);
  double neu_iso(max(0., lep.pfIsolationR04().sumNeutralHadronEt + lep.pfIsolationR04().sumPhotonEt
                     -rho*getEffAreaMuon(lep.eta())));

  return (ch_iso + neu_iso) / lep.pt();
}

//////////////////// Electrons
bool lepton_tools::isSignalElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= SignalLeptonPtCut) return false;
  if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idElectron(lep, vtx, kMedium)) return false;
  // Isolation cuts
  if(lepIso >= 0 && lepIso > ElectronMiniIsoCut) return false;

  return true;
}

bool lepton_tools::isVetoElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double lepIso){
  // pT, eta cuts
  if(lep.pt() <= VetoLeptonPtCut) return false;
  if(fabs(lep.superCluster()->position().eta()) > ElectronEtaCut) return false;
  // ID cuts (includes dz/dz0 cuts)
  if(!idElectron(lep, vtx, kVeto)) return false;
  // Isolation cuts
  if(lepIso >= 0 && lepIso > ElectronMiniIsoCut) return false;

  return true;
}

bool lepton_tools::idElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, CutLevel threshold, bool doIso) {

  bool barrel(lep.isEB());
  double deta_cut, dphi_cut, ieta_cut, hovere_cut, d0_cut, dz_cut,
    ooeminusoop_cut, reliso_cut, misshits_cut;
  bool req_conv_veto;
  
  // if(barrel){
  //   ieta_cut        = chooseVal(threshold       ,0.0115,  0.011,   0.00998, 0.00998);
  //   deta_cut        = chooseVal(threshold       ,0.00749, 0.00477, 0.00311, 0.00308);
  //   dphi_cut        = chooseVal(threshold       ,0.228,   0.222,   0.103,   0.0816);
  //   hovere_cut      = chooseVal(threshold       ,0.356,   0.298,   0.253,   0.0414);
  //   reliso_cut      = chooseVal(threshold       ,0.175,   0.0994,  0.0695,  0.0588);
  //   ooeminusoop_cut = chooseVal(threshold       ,0.299,   0.241,   0.134,   0.0129);
  //   d0_cut          = chooseVal(threshold       ,0.05,    0.05,    0.05,    0.05);
  //   dz_cut          = chooseVal(threshold       ,0.10 ,   0.10,    0.10,    0.10);
  //   misshits_cut    = chooseVal(threshold       ,2,   1,   1,   1);
  //   req_conv_veto   = chooseVal(threshold       ,true           ,  true         ,  true         ,  true );
  // } else {
  //   ieta_cut        = chooseVal(threshold       ,0.037,   0.0314,  0.0298,  0.0292);
  //   deta_cut        = chooseVal(threshold       ,0.00895, 0.00868, 0.00609, 0.00605);
  //   dphi_cut        = chooseVal(threshold       ,0.213,   0.213,   0.045,   0.0394);
  //   hovere_cut      = chooseVal(threshold       ,0.211,   0.101,   0.0878,  0.0641);
  //   reliso_cut      = chooseVal(threshold       ,0.159,   0.107,   0.0821,  0.0571);
  //   ooeminusoop_cut = chooseVal(threshold       ,0.15,    0.14,    0.13,    0.0129);
  //   d0_cut          = chooseVal(threshold       ,0.10 ,   0.10,    0.10,    0.10);
  //   dz_cut          = chooseVal(threshold       ,0.20 ,   0.20,    0.20,    0.20);
  //   misshits_cut    = chooseVal(threshold       ,3, 1, 1, 1);
  //   req_conv_veto   = chooseVal(threshold       ,true   ,  true   ,  true   ,  true );
  // }

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Working_points_for_Spring15_MC_s
  // Last updated October 8th
  if(barrel){
    ieta_cut        = chooseVal(threshold       ,0.0114,  0.0103,  0.0101,  0.0101);
    deta_cut        = chooseVal(threshold       ,0.0152,  0.0105,  0.0103,  0.00926);
    dphi_cut        = chooseVal(threshold       ,0.216,   0.115,   0.0336,  0.0336);
    hovere_cut      = chooseVal(threshold       ,0.181,   0.104,   0.0876,  0.0597);
    reliso_cut      = chooseVal(threshold       ,0.126,   0.0893,  0.0766,  0.0354);
    ooeminusoop_cut = chooseVal(threshold       ,0.207,   0.102,   0.0174,  0.012);
    d0_cut          = chooseVal(threshold       ,0.0564,  0.0261,  0.0118,  0.0111);
    dz_cut          = chooseVal(threshold       ,0.472,   0.41,  0.373,   0.0466);
    misshits_cut    = chooseVal(threshold       ,2,   2,   2,   2);
    req_conv_veto   = chooseVal(threshold       ,true           ,  true         ,  true         ,  true );
  } else {
    ieta_cut        = chooseVal(threshold       ,0.0352 , 0.0301 , 0.0283 , 0.0279);
    deta_cut        = chooseVal(threshold       ,0.0113 , 0.00814 , 0.00733 , 0.00724);
    dphi_cut        = chooseVal(threshold       ,0.237 , 0.182 , 0.114 , 0.0918);
    hovere_cut      = chooseVal(threshold       ,0.116 , 0.0897 , 0.0678 , 0.0615);
    reliso_cut      = chooseVal(threshold       ,0.144 , 0.121 , 0.0678 , 0.0646);
    ooeminusoop_cut = chooseVal(threshold       ,0.174 , 0.126 , 0.0898 , 0.00999);
    d0_cut          = chooseVal(threshold       ,0.222 , 0.118 , 0.0739 , 0.0351);
    dz_cut          = chooseVal(threshold       ,0.921 , 0.822 , 0.602 , 0.417);
    misshits_cut    = chooseVal(threshold       ,3, 1, 1, 1);
    req_conv_veto   = chooseVal(threshold       ,true   ,  true   ,  true   ,  true );
  }

  double dz(0.), d0(0.);
  vertexElectron(lep, vtx, dz, d0);

  int mhits(0);
  if(lep.gsfTrack().isAvailable()){
    mhits = lep.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);;
  }
  const double sigietaieta(lep.full5x5_sigmaIetaIeta());

  return deta_cut > fabs(lep.deltaEtaSuperClusterTrackAtVtx())
    && dphi_cut > fabs(lep.deltaPhiSuperClusterTrackAtVtx())
    && ieta_cut > sigietaieta
    && hovere_cut > lep.hadronicOverEm()
    && d0_cut > fabs(d0)
    && dz_cut > fabs(dz)
    && ooeminusoop_cut > fabs((1.0-lep.eSuperClusterOverP())/lep.ecalEnergy())
    && (!doIso || reliso_cut > 0) // To be implemented if we want reliso
    && (!req_conv_veto || lep.passConversionVeto())
    && (misshits_cut >= mhits);
}

bool lepton_tools::vertexElectron(const pat::Electron &lep, edm::Handle<reco::VertexCollection> vtx, double &dz, double &d0){
  dz = 0.; d0 = 0.;
  if(lep.gsfTrack().isAvailable()){ // If the track is not available we probably don't want the electron
    dz = lep.gsfTrack()->vz()-vtx->at(0).z();
    d0 = lep.gsfTrack()->d0()-vtx->at(0).x()*sin(lep.gsfTrack()->phi())+vtx->at(0).y()*cos(lep.gsfTrack()->phi());
  }
  if(fabs(dz) > 0.5 || fabs(d0) > 0.2) return false;

  return true;
}

double lepton_tools::getEffAreaElectron(double eta){
  double abseta = fabs(eta);
  if (abseta < 1) return 0.1752;
  else if (abseta < 1.479) return 0.1862;
  else if (abseta < 2.0) return 0.1411;
  else if (abseta < 2.2) return 0.1534;
  else if (abseta < 2.3) return 0.1903;
  else if (abseta < 2.4) return 0.2243;
  else if (abseta < 2.5) return 0.2687;
  else return 0;
}

double lepton_tools::getRelIsolation(const pat::Electron &lep, double rho){
  double ch_iso(lep.pfIsolationVariables().sumChargedHadronPt);
  double neu_iso(max(0., lep.pfIsolationVariables().sumNeutralHadronEt + lep.pfIsolationVariables().sumPhotonEt
                     -rho*getEffAreaElectron(lep.eta())));

  return (ch_iso + neu_iso) / lep.pt();
}

pair<double, double> lepton_tools::getScaleFactor(const reco::Candidate &cand){
  if(const pat::Electron * ele_ptr = dynamic_cast<const pat::Electron*>(&cand)){
    return getScaleFactor(*ele_ptr);
  }else if(const reco::Muon * mu_ptr = dynamic_cast<const reco::Muon*>(&cand)){
    return getScaleFactor(*mu_ptr);
  }else{
    ERROR("Cannot get scale factor for type " << typeid(cand).name());
  }
  return {1., 0.};
}

pair<double, double> lepton_tools::getScaleFactor(const vCands &sig_leps){
  return accumulate(sig_leps.cbegin(), sig_leps.cend(), make_pair(1., 0.),
                    [](pair<double, double> sf, const reco::Candidate* cand){
                      if(cand == nullptr) ERROR("Dereferencing nullptr");
                      return MergeSF(sf, getScaleFactor(*cand));
                    });
}

pair<double, double> lepton_tools::getScaleFactorFs(const reco::Candidate &cand){
  if(const pat::Electron * ele_ptr = dynamic_cast<const pat::Electron*>(&cand)){
    return getScaleFactorFs(*ele_ptr);
  }else if(const reco::Muon * mu_ptr = dynamic_cast<const reco::Muon*>(&cand)){
    return getScaleFactorFs(*mu_ptr);
  }else{
    ERROR("Cannot get scale factor for type " << typeid(cand).name());
  }
  return {1., 0.};
}

pair<double, double> lepton_tools::getScaleFactorFs(const vCands &sig_leps){
  return accumulate(sig_leps.cbegin(), sig_leps.cend(), make_pair(1., 0.),
                    [](pair<double, double> sf, const reco::Candidate* cand){
                      if(cand == nullptr) ERROR("Dereferencing nullptr");
                      return MergeSF(sf, getScaleFactorFs(*cand));
                    });
}

pair<double, double> lepton_tools::getScaleFactor(const reco::Muon &lep){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Data_leading_order_FullSim_MC_co
  //ID, iso, tracking SFs applied
  //No stat error, 3% systematic from ID, iso
  double pt = lep.pt();
  double eta = lep.eta();
  double abseta = fabs(eta);
  vector<pair<double, double> > sfs{
    GetSF(sf_full_muon_medium, pt, abseta, false),
      make_pair(1., 0.03),//Systematic uncertainty
      GetSF(sf_full_muon_iso, pt, abseta, false),
      make_pair(1., 0.03),//Systematic uncertainty
      GetSF(sf_full_muon_tracking, pt, eta)//Asymmetric in eta
      };
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

pair<double, double> lepton_tools::getScaleFactor(const pat::Electron &lep){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#Data_leading_order_FullSim_M_AN1
  //ID, iso, tracking SFs applied
  //ID iso systematics built-in
  //Tracking SFs from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
  //3% tracking systematic below 20 GeV
  double pt = lep.superCluster()->energy()*sin(lep.superClusterPosition().theta());
  double eta = lep.superCluster()->eta();
  double abseta = fabs(eta);
  vector<pair<double, double> > sfs{
    GetSF(sf_full_electron_medium, pt, abseta),
      GetSF(sf_full_electron_iso, pt, abseta),
      GetSF(sf_full_electron_tracking, eta, pt),//Axes swapped, asymmetric in eta
      make_pair(1., pt<20. ? 0.03 : 0.)//Systematic uncertainty
      };
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

pair<double, double> lepton_tools::getScaleFactorFs(const reco::Muon &lep){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_compari
  //ID, iso SFs applied
  //No stat error, 2% systematic from ID, iso
  double pt = lep.pt();
  double abseta = fabs(lep.eta());
  vector<pair<double, double> > sfs{
    GetSF(sf_fast_muon_medium, pt, abseta, false),
      make_pair(1., 0.02),
      GetSF(sf_fast_muon_iso, pt, abseta, false),
      make_pair(1., 0.02),
      };
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

pair<double, double> lepton_tools::getScaleFactorFs(const pat::Electron &lep){
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#FullSim_FastSim_TTBar_MC_com_AN1
  //ID, iso SFs applied
  //No stat error, 2% systematic from ID, iso
  double pt = lep.superCluster()->energy()*sin(lep.superClusterPosition().theta());
  double abseta = fabs(lep.superCluster()->eta());
  vector<pair<double, double> > sfs{
    GetSF(sf_fast_electron_mediumiso, pt, abseta, false),
      make_pair(1., 0.02),//Systematic uncertainty
      make_pair(1., 0.02)//Systematic uncertainty
      };
  return accumulate(sfs.cbegin(), sfs.cend(), make_pair(1., 0.), MergeSF);
}

double lepton_tools::getPFIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
                                    const reco::Candidate* ptcl,
                                    double r_iso_min, double r_iso_max, double kt_scale,
                                    double rho, bool charged_only) {
  if (ptcl->pt()<1.) return 99999.;

  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  if(ptcl->isElectron()) {
    if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  } else if(ptcl->isMuon()) {
    deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
  } else {
    //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
  }
  double iso_nh(0.), iso_ch(0.), iso_ph(0.), iso_pu(0.);
  double ptThresh(0.5), r_mini(kt_scale/ptcl->pt());
  if(ptcl->isElectron()) ptThresh = 0;
  double r_iso(max(r_iso_min, min(r_iso_max, r_mini)));

  for (const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;
    double dr = deltaR(pfc, *ptcl);
    if (dr > r_iso) continue;
    if (&pfc == ptcl) continue;
    if (pfc.charge()==0){ //neutrals
      if (pfc.pt()>ptThresh) {
        if (abs(pfc.pdgId())==22) { //photons
          if(dr < deadcone_ph) continue;
          iso_ph += pfc.pt();
        } else if (abs(pfc.pdgId())==130) { //neutral hadrons
          if(dr < deadcone_nh) continue;
          iso_nh += pfc.pt();
        }
      }
    } else if (pfc.fromPV()>1){ //charged from PV
      if (abs(pfc.pdgId())==211) {
        if(dr < deadcone_ch) continue;
        iso_ch += pfc.pt();
      }
    } else {
      if (pfc.pt()>ptThresh){ //charged from PU
        if(dr < deadcone_pu) continue;
        iso_pu += pfc.pt();
      }
    }
  } // Loop over pf cands

  double effarea = ptcl->isElectron() ? getEffAreaElectron(ptcl->eta()) : getEffAreaMuon(ptcl->eta());
  double pu_corr = rho*effarea*pow(r_iso/0.3, 2);
  double iso(0.);
  if (charged_only) iso = iso_ch;
  else iso = iso_ch + max(0.,iso_ph + iso_nh - pu_corr);

  return iso/ptcl->pt();
}

vCands lepton_tools::getIsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi){

  vCands tks;
  //common parameters
  float mt_max = 100.;
  float eta_max = 2.5;
  float dz_max = 0.1;
  float cone_size = 0.3;

  for (size_t i(0); i < pfcands->size(); i++) {
    const pat::PackedCandidate &tk = (*pfcands)[i];
    unsigned int id = abs(tk.pdgId());

    //id-specific parameters
    float pt_min = id==211 ? 10. : 5.;
    float iso_max = id==211 ? 0.1 : 0.2;

    // track selection
    if (tk.charge()==0) continue;
    if (id!=11 && id!=13 && id!=211) continue;
    if (tk.pt() < pt_min) continue;
    if (fabs(tk.eta()) > eta_max) continue;
    if (fabs(tk.dz()) > dz_max) continue;
    if (mt_max>0.01 && getMT(met, met_phi,  tk.pt(), tk.phi())>mt_max) continue;

    // calculate track isolation
    double iso = 0.;
    for (size_t j(0); j < pfcands->size(); j++) {
      if (i==j) continue;
      const pat::PackedCandidate &pfc = (*pfcands)[j];

      if (pfc.charge()==0) continue;
      if (deltaR(tk,pfc) > cone_size) continue;
      if (fabs(pfc.dz()) > dz_max) continue;
      iso += pfc.pt();
    }
    if (iso/tk.pt()>iso_max) continue;

    tks.push_back(dynamic_cast<const reco::Candidate *>(&tk));
  }

  return tks;
}

vCands lepton_tools::getRA4IsoTracks(edm::Handle<pat::PackedCandidateCollection> pfcands, double met, double met_phi, double rhoEventCentral, vector<float> &isos, vector<float> &relisos, int primary_pdg){

  vCands tks;
  //Very loose cuts to have flexibility, not the veto definition
  float dz_max = 0.5;
  float pt_min = 5;
  float iso_max = 1.;

  for (size_t i(0); i < pfcands->size(); i++) {
    const pat::PackedCandidate &tk = (*pfcands)[i];
    unsigned int id = abs(tk.pdgId());

    // track selection
    if (tk.charge()==0) continue;
    if (id!=11 && id!=13 && id!=211) continue;
    if (tk.pt() < pt_min) continue;
    if (fabs(tk.dz()) > dz_max) continue;

    //Only save opposite sign tracks
    if(id == 211){
      if(tk.pdgId() * primary_pdg < 0 ) continue;
    }
    else{
      if(tk.pdgId() * primary_pdg > 0 ) continue;
    }
     

    // calculate track isolation
    double iso = 0.;
    double reliso=0.;
    if(id!=211){ //Use charged+neutral for e's and mu's
      iso = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&tk), 0.05, 0.2, 10., rhoEventCentral, false);
      reliso = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&tk), 0.3, 0.3, 10., rhoEventCentral, false);
    }
    else{ //Use charged only for hadronic tracks
      iso = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&tk), 0.05, 0.2, 10., rhoEventCentral, true);
      reliso =  getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&tk), 0.3, 0.3, 10., rhoEventCentral, true);
    }
    
    if(iso>iso_max && reliso>iso_max) continue; //Save tracks that pass a very loose isolation cut, for either isolation
    isos.push_back(iso);
    relisos.push_back(reliso);
    tks.push_back(dynamic_cast<const reco::Candidate *>(&tk));
  }

  return tks;
}

// Bad Global Muon filter
bool lepton_tools::outInOnly(const reco::Muon &mu) const {
  const reco::Track &tk = *mu.innerTrack();
  return tk.algoMask().count() == 1 && tk.isAlgoInMask(reco::Track::muonSeededStepOutIn);
}
bool lepton_tools::preselection(const reco::Muon &mu, bool selectClones) const { 
  return (!selectClones || outInOnly(mu));
}
bool lepton_tools::tighterId(const reco::Muon &mu) const { 
  return muon::isMediumMuon(mu) && mu.numberOfMatchedStations() >= 2; 
}
bool lepton_tools::tightGlobal(const reco::Muon &mu) const {
  return mu.isGlobalMuon() && (mu.globalTrack()->hitPattern().muonStationsWithValidHits() >= 3 && mu.globalTrack()->normalizedChi2() <= 20);
}
bool lepton_tools::safeId(const reco::Muon &mu) const { 
  if (mu.muonBestTrack()->ptError() > 0.2 * mu.muonBestTrack()->pt()) { return false; }
  return mu.numberOfMatchedStations() >= 1 || tightGlobal(mu);
}
bool lepton_tools::partnerId(const reco::Muon &mu) const {
  return mu.pt() >= 10 && mu.numberOfMatchedStations() >= 1;
}

set<unsigned> lepton_tools::badGlobalMuonSelector(edm::Handle<reco::VertexCollection> vtx, 
                                         edm::Handle<pat::MuonCollection> muptr, bool selectClones) {
    using namespace edm;
    float ptCut_ = 20.;
    bool verbose_ = false;
    assert(vtx->size() >= 1);
    const auto &PV = vtx->front().position();

    set<unsigned> badmus; 
    std::vector<int> goodMuon;
    const pat::MuonCollection &muons = *muptr;
    for (auto & mu : muons) {
        if (!mu.isPFMuon() || mu.innerTrack().isNull()) {
            goodMuon.push_back(-1); // bad but we don't care
            continue;
        } 
        if (preselection(mu, selectClones)) {
            float dxypv = std::abs(mu.innerTrack()->dxy(PV));
            float dzpv  = std::abs(mu.innerTrack()->dz(PV));
            if (tighterId(mu)) {
                bool ipLoose = ((dxypv < 0.5 && dzpv < 2.0) || mu.innerTrack()->hitPattern().pixelLayersWithMeasurement() >= 2);
                goodMuon.push_back(ipLoose || (!selectClones && tightGlobal(mu)));
            } else if (safeId(mu)) {
                bool ipTight = (dxypv < 0.2 && dzpv < 0.5);
                goodMuon.push_back(ipTight);
           } else {
                goodMuon.push_back(0);
            }
        } else {
            goodMuon.push_back(3); // maybe good, maybe bad, but we don't care
        }
    }

    for (unsigned int i = 0, n = muons.size(); i < n; ++i) {
        if (muons[i].pt() < ptCut_ || goodMuon[i] != 0) continue;
        if (verbose_) printf("potentially bad muon %d of pt %.1f eta %+.3f phi %+.3f\n", int(i+1), muons[i].pt(), muons[i].eta(), muons[i].phi());
        bool bad = true;
        if (selectClones) {
            bad = false; // unless proven otherwise
            unsigned int n1 = muons[i].numberOfMatches(reco::Muon::SegmentArbitration);
            for (unsigned int j = 0; j < n; ++j) {
                if (j == i || goodMuon[j] <= 0 || !partnerId(muons[j])) continue;
                unsigned int n2 = muons[j].numberOfMatches(reco::Muon::SegmentArbitration);
                if (deltaR2(muons[i],muons[j]) < 0.16 || (n1 > 0 && n2 > 0 && muon::sharedSegments(muons[i],muons[j]) >= 0.5*std::min(n1,n2))) {
                    if (verbose_) printf("     tagged as clone of muon %d of pt %.1f eta %+.3f phi %+.3f\n", int(j+1), muons[j].pt(), muons[j].eta(), muons[j].phi());
                    bad = true;
                    break;
                } 
            }
        }
        if (bad) badmus.insert(i);
    }

    return badmus;
}


lepton_tools::lepton_tools(){
}

lepton_tools::~lepton_tools(){
}
