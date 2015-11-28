// baby_base: base class to handle reduce tree ntuples
//File generated with generate_baby.exe

#include "baby_base.hh"

#include <stdexcept>
#include <string>
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
using namespace std;

bool baby_base::VectorLoader::loaded_ = false;

baby_base::VectorLoader baby_base::vl_ = baby_base::VectorLoader();

baby_base::VectorLoader::VectorLoader(){
  if(!loaded_){
    gROOT->ProcessLine("#include <vector>");
    loaded_ = true;
  }
}

const double baby_base::bad_val_ = -999.;

baby_base::baby_base():
  chain_("junk", "junk"),
  tree_("tree", "tree"),
  entry_(0),
  read_only_(false),
  event_(static_cast<Long64_t >(bad_val_)),
  b_event_(tree_.Branch("event", &event_)),
  c_event_(false),
  fromGS_(static_cast<bool >(bad_val_)),
  b_fromGS_(tree_.Branch("fromGS", &fromGS_)),
  c_fromGS_(false),
  jetmismeas_(static_cast<bool >(bad_val_)),
  b_jetmismeas_(tree_.Branch("jetmismeas", &jetmismeas_)),
  c_jetmismeas_(false),
  low_dphi_(static_cast<bool >(bad_val_)),
  b_low_dphi_(tree_.Branch("low_dphi", &low_dphi_)),
  c_low_dphi_(false),
  nonblind_(static_cast<bool >(bad_val_)),
  b_nonblind_(tree_.Branch("nonblind", &nonblind_)),
  c_nonblind_(false),
  pass_(static_cast<bool >(bad_val_)),
  b_pass_(tree_.Branch("pass", &pass_)),
  c_pass_(false),
  pass_cschalo_(static_cast<bool >(bad_val_)),
  b_pass_cschalo_(tree_.Branch("pass_cschalo", &pass_cschalo_)),
  c_pass_cschalo_(false),
  pass_eebadsc_(static_cast<bool >(bad_val_)),
  b_pass_eebadsc_(tree_.Branch("pass_eebadsc", &pass_eebadsc_)),
  c_pass_eebadsc_(false),
  pass_goodv_(static_cast<bool >(bad_val_)),
  b_pass_goodv_(tree_.Branch("pass_goodv", &pass_goodv_)),
  c_pass_goodv_(false),
  pass_hbhe_(static_cast<bool >(bad_val_)),
  b_pass_hbhe_(tree_.Branch("pass_hbhe", &pass_hbhe_)),
  c_pass_hbhe_(false),
  pass_hbheiso_(static_cast<bool >(bad_val_)),
  b_pass_hbheiso_(tree_.Branch("pass_hbheiso", &pass_hbheiso_)),
  c_pass_hbheiso_(false),
  pass_jets_(static_cast<bool >(bad_val_)),
  b_pass_jets_(tree_.Branch("pass_jets", &pass_jets_)),
  c_pass_jets_(false),
  pass_jets_nohf_(static_cast<bool >(bad_val_)),
  b_pass_jets_nohf_(tree_.Branch("pass_jets_nohf", &pass_jets_nohf_)),
  c_pass_jets_nohf_(false),
  pass_jets_ra2_(static_cast<bool >(bad_val_)),
  b_pass_jets_ra2_(tree_.Branch("pass_jets_ra2", &pass_jets_ra2_)),
  c_pass_jets_ra2_(false),
  pass_jets_tight_(static_cast<bool >(bad_val_)),
  b_pass_jets_tight_(tree_.Branch("pass_jets_tight", &pass_jets_tight_)),
  c_pass_jets_tight_(false),
  pass_jets_tight_ra2_(static_cast<bool >(bad_val_)),
  b_pass_jets_tight_ra2_(tree_.Branch("pass_jets_tight_ra2", &pass_jets_tight_ra2_)),
  c_pass_jets_tight_ra2_(false),
  pass_nohf_(static_cast<bool >(bad_val_)),
  b_pass_nohf_(tree_.Branch("pass_nohf", &pass_nohf_)),
  c_pass_nohf_(false),
  pass_ra2_(static_cast<bool >(bad_val_)),
  b_pass_ra2_(tree_.Branch("pass_ra2", &pass_ra2_)),
  c_pass_ra2_(false),
  stitch_(static_cast<bool >(bad_val_)),
  b_stitch_(tree_.Branch("stitch", &stitch_)),
  c_stitch_(false),
  dphi1_(static_cast<float >(bad_val_)),
  b_dphi1_(tree_.Branch("dphi1", &dphi1_)),
  c_dphi1_(false),
  dphi2_(static_cast<float >(bad_val_)),
  b_dphi2_(tree_.Branch("dphi2", &dphi2_)),
  c_dphi2_(false),
  dphi3_(static_cast<float >(bad_val_)),
  b_dphi3_(tree_.Branch("dphi3", &dphi3_)),
  c_dphi3_(false),
  dphi4_(static_cast<float >(bad_val_)),
  b_dphi4_(tree_.Branch("dphi4", &dphi4_)),
  c_dphi4_(false),
  dphi_wlep_(static_cast<float >(bad_val_)),
  b_dphi_wlep_(tree_.Branch("dphi_wlep", &dphi_wlep_)),
  c_dphi_wlep_(false),
  eff_trig_(static_cast<float >(bad_val_)),
  b_eff_trig_(tree_.Branch("eff_trig", &eff_trig_)),
  c_eff_trig_(false),
  elel_eta_(static_cast<float >(bad_val_)),
  b_elel_eta_(tree_.Branch("elel_eta", &elel_eta_)),
  c_elel_eta_(false),
  elel_m_(static_cast<float >(bad_val_)),
  b_elel_m_(tree_.Branch("elel_m", &elel_m_)),
  c_elel_m_(false),
  elel_phi_(static_cast<float >(bad_val_)),
  b_elel_phi_(tree_.Branch("elel_phi", &elel_phi_)),
  c_elel_phi_(false),
  elel_pt_(static_cast<float >(bad_val_)),
  b_elel_pt_(tree_.Branch("elel_pt", &elel_pt_)),
  c_elel_pt_(false),
  elel_pt1_(static_cast<float >(bad_val_)),
  b_elel_pt1_(tree_.Branch("elel_pt1", &elel_pt1_)),
  c_elel_pt1_(false),
  elel_pt2_(static_cast<float >(bad_val_)),
  b_elel_pt2_(tree_.Branch("elel_pt2", &elel_pt2_)),
  c_elel_pt2_(false),
  elel_w_(static_cast<float >(bad_val_)),
  b_elel_w_(tree_.Branch("elel_w", &elel_w_)),
  c_elel_w_(false),
  elelv_eta_(static_cast<float >(bad_val_)),
  b_elelv_eta_(tree_.Branch("elelv_eta", &elelv_eta_)),
  c_elelv_eta_(false),
  elelv_m_(static_cast<float >(bad_val_)),
  b_elelv_m_(tree_.Branch("elelv_m", &elelv_m_)),
  c_elelv_m_(false),
  elelv_phi_(static_cast<float >(bad_val_)),
  b_elelv_phi_(tree_.Branch("elelv_phi", &elelv_phi_)),
  c_elelv_phi_(false),
  elelv_pt_(static_cast<float >(bad_val_)),
  b_elelv_pt_(tree_.Branch("elelv_pt", &elelv_pt_)),
  c_elelv_pt_(false),
  elelv_pt1_(static_cast<float >(bad_val_)),
  b_elelv_pt1_(tree_.Branch("elelv_pt1", &elelv_pt1_)),
  c_elelv_pt1_(false),
  elelv_pt2_(static_cast<float >(bad_val_)),
  b_elelv_pt2_(tree_.Branch("elelv_pt2", &elelv_pt2_)),
  c_elelv_pt2_(false),
  elelv_w_(static_cast<float >(bad_val_)),
  b_elelv_w_(tree_.Branch("elelv_w", &elelv_w_)),
  c_elelv_w_(false),
  elmu_eta_(static_cast<float >(bad_val_)),
  b_elmu_eta_(tree_.Branch("elmu_eta", &elmu_eta_)),
  c_elmu_eta_(false),
  elmu_m_(static_cast<float >(bad_val_)),
  b_elmu_m_(tree_.Branch("elmu_m", &elmu_m_)),
  c_elmu_m_(false),
  elmu_phi_(static_cast<float >(bad_val_)),
  b_elmu_phi_(tree_.Branch("elmu_phi", &elmu_phi_)),
  c_elmu_phi_(false),
  elmu_pt_(static_cast<float >(bad_val_)),
  b_elmu_pt_(tree_.Branch("elmu_pt", &elmu_pt_)),
  c_elmu_pt_(false),
  elmu_pt1_(static_cast<float >(bad_val_)),
  b_elmu_pt1_(tree_.Branch("elmu_pt1", &elmu_pt1_)),
  c_elmu_pt1_(false),
  elmu_pt2_(static_cast<float >(bad_val_)),
  b_elmu_pt2_(tree_.Branch("elmu_pt2", &elmu_pt2_)),
  c_elmu_pt2_(false),
  elmu_w_(static_cast<float >(bad_val_)),
  b_elmu_w_(tree_.Branch("elmu_w", &elmu_w_)),
  c_elmu_w_(false),
  ht_(static_cast<float >(bad_val_)),
  b_ht_(tree_.Branch("ht", &ht_)),
  c_ht_(false),
  ht_clean_(static_cast<float >(bad_val_)),
  b_ht_clean_(tree_.Branch("ht_clean", &ht_clean_)),
  c_ht_clean_(false),
  ht_hlt_(static_cast<float >(bad_val_)),
  b_ht_hlt_(tree_.Branch("ht_hlt", &ht_hlt_)),
  c_ht_hlt_(false),
  ht_isr_me_(static_cast<float >(bad_val_)),
  b_ht_isr_me_(tree_.Branch("ht_isr_me", &ht_isr_me_)),
  c_ht_isr_me_(false),
  ht_ra2_(static_cast<float >(bad_val_)),
  b_ht_ra2_(tree_.Branch("ht_ra2", &ht_ra2_)),
  c_ht_ra2_(false),
  ht_tru_(static_cast<float >(bad_val_)),
  b_ht_tru_(tree_.Branch("ht_tru", &ht_tru_)),
  c_ht_tru_(false),
  isr_tru_eta_(static_cast<float >(bad_val_)),
  b_isr_tru_eta_(tree_.Branch("isr_tru_eta", &isr_tru_eta_)),
  c_isr_tru_eta_(false),
  isr_tru_phi_(static_cast<float >(bad_val_)),
  b_isr_tru_phi_(tree_.Branch("isr_tru_phi", &isr_tru_phi_)),
  c_isr_tru_phi_(false),
  isr_tru_pt_(static_cast<float >(bad_val_)),
  b_isr_tru_pt_(tree_.Branch("isr_tru_pt", &isr_tru_pt_)),
  c_isr_tru_pt_(false),
  jetsys_eta_(static_cast<float >(bad_val_)),
  b_jetsys_eta_(tree_.Branch("jetsys_eta", &jetsys_eta_)),
  c_jetsys_eta_(false),
  jetsys_nob_eta_(static_cast<float >(bad_val_)),
  b_jetsys_nob_eta_(tree_.Branch("jetsys_nob_eta", &jetsys_nob_eta_)),
  c_jetsys_nob_eta_(false),
  jetsys_nob_phi_(static_cast<float >(bad_val_)),
  b_jetsys_nob_phi_(tree_.Branch("jetsys_nob_phi", &jetsys_nob_phi_)),
  c_jetsys_nob_phi_(false),
  jetsys_nob_pt_(static_cast<float >(bad_val_)),
  b_jetsys_nob_pt_(tree_.Branch("jetsys_nob_pt", &jetsys_nob_pt_)),
  c_jetsys_nob_pt_(false),
  jetsys_phi_(static_cast<float >(bad_val_)),
  b_jetsys_phi_(tree_.Branch("jetsys_phi", &jetsys_phi_)),
  c_jetsys_phi_(false),
  jetsys_pt_(static_cast<float >(bad_val_)),
  b_jetsys_pt_(tree_.Branch("jetsys_pt", &jetsys_pt_)),
  c_jetsys_pt_(false),
  met_(static_cast<float >(bad_val_)),
  b_met_(tree_.Branch("met", &met_)),
  c_met_(false),
  met_calo_(static_cast<float >(bad_val_)),
  b_met_calo_(tree_.Branch("met_calo", &met_calo_)),
  c_met_calo_(false),
  met_calo_phi_(static_cast<float >(bad_val_)),
  b_met_calo_phi_(tree_.Branch("met_calo_phi", &met_calo_phi_)),
  c_met_calo_phi_(false),
  met_mini_(static_cast<float >(bad_val_)),
  b_met_mini_(tree_.Branch("met_mini", &met_mini_)),
  c_met_mini_(false),
  met_mini_phi_(static_cast<float >(bad_val_)),
  b_met_mini_phi_(tree_.Branch("met_mini_phi", &met_mini_phi_)),
  c_met_mini_phi_(false),
  met_nohf_(static_cast<float >(bad_val_)),
  b_met_nohf_(tree_.Branch("met_nohf", &met_nohf_)),
  c_met_nohf_(false),
  met_nohf_phi_(static_cast<float >(bad_val_)),
  b_met_nohf_phi_(tree_.Branch("met_nohf_phi", &met_nohf_phi_)),
  c_met_nohf_phi_(false),
  met_phi_(static_cast<float >(bad_val_)),
  b_met_phi_(tree_.Branch("met_phi", &met_phi_)),
  c_met_phi_(false),
  met_raw_(static_cast<float >(bad_val_)),
  b_met_raw_(tree_.Branch("met_raw", &met_raw_)),
  c_met_raw_(false),
  met_raw_phi_(static_cast<float >(bad_val_)),
  b_met_raw_phi_(tree_.Branch("met_raw_phi", &met_raw_phi_)),
  c_met_raw_phi_(false),
  met_rebal_(static_cast<float >(bad_val_)),
  b_met_rebal_(tree_.Branch("met_rebal", &met_rebal_)),
  c_met_rebal_(false),
  met_tru_(static_cast<float >(bad_val_)),
  b_met_tru_(tree_.Branch("met_tru", &met_tru_)),
  c_met_tru_(false),
  met_tru_nuw_(static_cast<float >(bad_val_)),
  b_met_tru_nuw_(tree_.Branch("met_tru_nuw", &met_tru_nuw_)),
  c_met_tru_nuw_(false),
  met_tru_nuw_phi_(static_cast<float >(bad_val_)),
  b_met_tru_nuw_phi_(tree_.Branch("met_tru_nuw_phi", &met_tru_nuw_phi_)),
  c_met_tru_nuw_phi_(false),
  met_tru_phi_(static_cast<float >(bad_val_)),
  b_met_tru_phi_(tree_.Branch("met_tru_phi", &met_tru_phi_)),
  c_met_tru_phi_(false),
  mht_(static_cast<float >(bad_val_)),
  b_mht_(tree_.Branch("mht", &mht_)),
  c_mht_(false),
  mht_clean_(static_cast<float >(bad_val_)),
  b_mht_clean_(tree_.Branch("mht_clean", &mht_clean_)),
  c_mht_clean_(false),
  mht_clean_phi_(static_cast<float >(bad_val_)),
  b_mht_clean_phi_(tree_.Branch("mht_clean_phi", &mht_clean_phi_)),
  c_mht_clean_phi_(false),
  mht_phi_(static_cast<float >(bad_val_)),
  b_mht_phi_(tree_.Branch("mht_phi", &mht_phi_)),
  c_mht_phi_(false),
  mj_(static_cast<float >(bad_val_)),
  b_mj_(tree_.Branch("mj", &mj_)),
  c_mj_(false),
  mj08_(static_cast<float >(bad_val_)),
  b_mj08_(tree_.Branch("mj08", &mj08_)),
  c_mj08_(false),
  mt_(static_cast<float >(bad_val_)),
  b_mt_(tree_.Branch("mt", &mt_)),
  c_mt_(false),
  mt_nohf_(static_cast<float >(bad_val_)),
  b_mt_nohf_(tree_.Branch("mt_nohf", &mt_nohf_)),
  c_mt_nohf_(false),
  mt_rebal_(static_cast<float >(bad_val_)),
  b_mt_rebal_(tree_.Branch("mt_rebal", &mt_rebal_)),
  c_mt_rebal_(false),
  mt_tru_(static_cast<float >(bad_val_)),
  b_mt_tru_(tree_.Branch("mt_tru", &mt_tru_)),
  c_mt_tru_(false),
  mt_tru_nuw_(static_cast<float >(bad_val_)),
  b_mt_tru_nuw_(tree_.Branch("mt_tru_nuw", &mt_tru_nuw_)),
  c_mt_tru_nuw_(false),
  mumu_eta_(static_cast<float >(bad_val_)),
  b_mumu_eta_(tree_.Branch("mumu_eta", &mumu_eta_)),
  c_mumu_eta_(false),
  mumu_m_(static_cast<float >(bad_val_)),
  b_mumu_m_(tree_.Branch("mumu_m", &mumu_m_)),
  c_mumu_m_(false),
  mumu_phi_(static_cast<float >(bad_val_)),
  b_mumu_phi_(tree_.Branch("mumu_phi", &mumu_phi_)),
  c_mumu_phi_(false),
  mumu_pt_(static_cast<float >(bad_val_)),
  b_mumu_pt_(tree_.Branch("mumu_pt", &mumu_pt_)),
  c_mumu_pt_(false),
  mumu_pt1_(static_cast<float >(bad_val_)),
  b_mumu_pt1_(tree_.Branch("mumu_pt1", &mumu_pt1_)),
  c_mumu_pt1_(false),
  mumu_pt2_(static_cast<float >(bad_val_)),
  b_mumu_pt2_(tree_.Branch("mumu_pt2", &mumu_pt2_)),
  c_mumu_pt2_(false),
  mumu_w_(static_cast<float >(bad_val_)),
  b_mumu_w_(tree_.Branch("mumu_w", &mumu_w_)),
  c_mumu_w_(false),
  mumuv_eta_(static_cast<float >(bad_val_)),
  b_mumuv_eta_(tree_.Branch("mumuv_eta", &mumuv_eta_)),
  c_mumuv_eta_(false),
  mumuv_m_(static_cast<float >(bad_val_)),
  b_mumuv_m_(tree_.Branch("mumuv_m", &mumuv_m_)),
  c_mumuv_m_(false),
  mumuv_phi_(static_cast<float >(bad_val_)),
  b_mumuv_phi_(tree_.Branch("mumuv_phi", &mumuv_phi_)),
  c_mumuv_phi_(false),
  mumuv_pt_(static_cast<float >(bad_val_)),
  b_mumuv_pt_(tree_.Branch("mumuv_pt", &mumuv_pt_)),
  c_mumuv_pt_(false),
  mumuv_pt1_(static_cast<float >(bad_val_)),
  b_mumuv_pt1_(tree_.Branch("mumuv_pt1", &mumuv_pt1_)),
  c_mumuv_pt1_(false),
  mumuv_pt2_(static_cast<float >(bad_val_)),
  b_mumuv_pt2_(tree_.Branch("mumuv_pt2", &mumuv_pt2_)),
  c_mumuv_pt2_(false),
  mumuv_w_(static_cast<float >(bad_val_)),
  b_mumuv_w_(tree_.Branch("mumuv_w", &mumuv_w_)),
  c_mumuv_w_(false),
  ntrupv_mean_(static_cast<float >(bad_val_)),
  b_ntrupv_mean_(tree_.Branch("ntrupv_mean", &ntrupv_mean_)),
  c_ntrupv_mean_(false),
  onel_ele105_(static_cast<float >(bad_val_)),
  b_onel_ele105_(tree_.Branch("onel_ele105", &onel_ele105_)),
  c_onel_ele105_(false),
  onel_ele23_(static_cast<float >(bad_val_)),
  b_onel_ele23_(tree_.Branch("onel_ele23", &onel_ele23_)),
  c_onel_ele23_(false),
  onel_ele8_(static_cast<float >(bad_val_)),
  b_onel_ele8_(tree_.Branch("onel_ele8", &onel_ele8_)),
  c_onel_ele8_(false),
  onel_vvvl_(static_cast<float >(bad_val_)),
  b_onel_vvvl_(tree_.Branch("onel_vvvl", &onel_vvvl_)),
  c_onel_vvvl_(false),
  onht_(static_cast<float >(bad_val_)),
  b_onht_(tree_.Branch("onht", &onht_)),
  c_onht_(false),
  onmet_(static_cast<float >(bad_val_)),
  b_onmet_(tree_.Branch("onmet", &onmet_)),
  c_onmet_(false),
  onmu_isomu18_(static_cast<float >(bad_val_)),
  b_onmu_isomu18_(tree_.Branch("onmu_isomu18", &onmu_isomu18_)),
  c_onmu_isomu18_(false),
  onmu_mu50_(static_cast<float >(bad_val_)),
  b_onmu_mu50_(tree_.Branch("onmu_mu50", &onmu_mu50_)),
  c_onmu_mu50_(false),
  onmu_mu8_(static_cast<float >(bad_val_)),
  b_onmu_mu8_(tree_.Branch("onmu_mu8", &onmu_mu8_)),
  c_onmu_mu8_(false),
  onmu_vvvl_(static_cast<float >(bad_val_)),
  b_onmu_vvvl_(tree_.Branch("onmu_vvvl", &onmu_vvvl_)),
  c_onmu_vvvl_(false),
  onph_ph90_(static_cast<float >(bad_val_)),
  b_onph_ph90_(tree_.Branch("onph_ph90", &onph_ph90_)),
  c_onph_ph90_(false),
  w_btag_(static_cast<float >(bad_val_)),
  b_w_btag_(tree_.Branch("w_btag", &w_btag_)),
  c_w_btag_(false),
  w_fs_lep_(static_cast<float >(bad_val_)),
  b_w_fs_lep_(tree_.Branch("w_fs_lep", &w_fs_lep_)),
  c_w_fs_lep_(false),
  w_lep_(static_cast<float >(bad_val_)),
  b_w_lep_(tree_.Branch("w_lep", &w_lep_)),
  c_w_lep_(false),
  w_lumi_(static_cast<float >(bad_val_)),
  b_w_lumi_(tree_.Branch("w_lumi", &w_lumi_)),
  c_w_lumi_(false),
  w_pu_(static_cast<float >(bad_val_)),
  b_w_pu_(tree_.Branch("w_pu", &w_pu_)),
  c_w_pu_(false),
  w_toppt_(static_cast<float >(bad_val_)),
  b_w_toppt_(tree_.Branch("w_toppt", &w_toppt_)),
  c_w_toppt_(false),
  weight_(static_cast<float >(bad_val_)),
  b_weight_(tree_.Branch("weight", &weight_)),
  c_weight_(false),
  lumiblock_(static_cast<int >(bad_val_)),
  b_lumiblock_(tree_.Branch("lumiblock", &lumiblock_)),
  c_lumiblock_(false),
  mgluino_(static_cast<int >(bad_val_)),
  b_mgluino_(tree_.Branch("mgluino", &mgluino_)),
  c_mgluino_(false),
  mlsp_(static_cast<int >(bad_val_)),
  b_mlsp_(tree_.Branch("mlsp", &mlsp_)),
  c_mlsp_(false),
  nbl_(static_cast<int >(bad_val_)),
  b_nbl_(tree_.Branch("nbl", &nbl_)),
  c_nbl_(false),
  nbm_(static_cast<int >(bad_val_)),
  b_nbm_(tree_.Branch("nbm", &nbm_)),
  c_nbm_(false),
  nbm_ra2_(static_cast<int >(bad_val_)),
  b_nbm_ra2_(tree_.Branch("nbm_ra2", &nbm_ra2_)),
  c_nbm_ra2_(false),
  nbt_(static_cast<int >(bad_val_)),
  b_nbt_(tree_.Branch("nbt", &nbt_)),
  c_nbt_(false),
  nels_(static_cast<int >(bad_val_)),
  b_nels_(tree_.Branch("nels", &nels_)),
  c_nels_(false),
  nels_ele23_(static_cast<int >(bad_val_)),
  b_nels_ele23_(tree_.Branch("nels_ele23", &nels_ele23_)),
  c_nels_ele23_(false),
  nels_vvvl_(static_cast<int >(bad_val_)),
  b_nels_vvvl_(tree_.Branch("nels_vvvl", &nels_vvvl_)),
  c_nels_vvvl_(false),
  nfjets_(static_cast<int >(bad_val_)),
  b_nfjets_(tree_.Branch("nfjets", &nfjets_)),
  c_nfjets_(false),
  nfjets08_(static_cast<int >(bad_val_)),
  b_nfjets08_(tree_.Branch("nfjets08", &nfjets08_)),
  c_nfjets08_(false),
  nisr_me_(static_cast<int >(bad_val_)),
  b_nisr_me_(tree_.Branch("nisr_me", &nisr_me_)),
  c_nisr_me_(false),
  njets_(static_cast<int >(bad_val_)),
  b_njets_(tree_.Branch("njets", &njets_)),
  c_njets_(false),
  njets_clean_(static_cast<int >(bad_val_)),
  b_njets_clean_(tree_.Branch("njets_clean", &njets_clean_)),
  c_njets_clean_(false),
  njets_ra2_(static_cast<int >(bad_val_)),
  b_njets_ra2_(tree_.Branch("njets_ra2", &njets_ra2_)),
  c_njets_ra2_(false),
  nleps_(static_cast<int >(bad_val_)),
  b_nleps_(tree_.Branch("nleps", &nleps_)),
  c_nleps_(false),
  nleps_tm_(static_cast<int >(bad_val_)),
  b_nleps_tm_(tree_.Branch("nleps_tm", &nleps_tm_)),
  c_nleps_tm_(false),
  nmus_(static_cast<int >(bad_val_)),
  b_nmus_(tree_.Branch("nmus", &nmus_)),
  c_nmus_(false),
  nmus_isomu18_(static_cast<int >(bad_val_)),
  b_nmus_isomu18_(tree_.Branch("nmus_isomu18", &nmus_isomu18_)),
  c_nmus_isomu18_(false),
  nmus_vvvl_(static_cast<int >(bad_val_)),
  b_nmus_vvvl_(tree_.Branch("nmus_vvvl", &nmus_vvvl_)),
  c_nmus_vvvl_(false),
  nph_(static_cast<int >(bad_val_)),
  b_nph_(tree_.Branch("nph", &nph_)),
  c_nph_(false),
  npv_(static_cast<int >(bad_val_)),
  b_npv_(tree_.Branch("npv", &npv_)),
  c_npv_(false),
  ntks_(static_cast<int >(bad_val_)),
  b_ntks_(tree_.Branch("ntks", &ntks_)),
  c_ntks_(false),
  ntruels_(static_cast<int >(bad_val_)),
  b_ntruels_(tree_.Branch("ntruels", &ntruels_)),
  c_ntruels_(false),
  ntruleps_(static_cast<int >(bad_val_)),
  b_ntruleps_(tree_.Branch("ntruleps", &ntruleps_)),
  c_ntruleps_(false),
  ntrumus_(static_cast<int >(bad_val_)),
  b_ntrumus_(tree_.Branch("ntrumus", &ntrumus_)),
  c_ntrumus_(false),
  ntrupv_(static_cast<int >(bad_val_)),
  b_ntrupv_(tree_.Branch("ntrupv", &ntrupv_)),
  c_ntrupv_(false),
  ntrutaush_(static_cast<int >(bad_val_)),
  b_ntrutaush_(tree_.Branch("ntrutaush", &ntrutaush_)),
  c_ntrutaush_(false),
  ntrutausl_(static_cast<int >(bad_val_)),
  b_ntrutausl_(tree_.Branch("ntrutausl", &ntrutausl_)),
  c_ntrutausl_(false),
  nvels_(static_cast<int >(bad_val_)),
  b_nvels_(tree_.Branch("nvels", &nvels_)),
  c_nvels_(false),
  nveto_(static_cast<int >(bad_val_)),
  b_nveto_(tree_.Branch("nveto", &nveto_)),
  c_nveto_(false),
  nvleps_(static_cast<int >(bad_val_)),
  b_nvleps_(tree_.Branch("nvleps", &nvleps_)),
  c_nvleps_(false),
  nvmus_(static_cast<int >(bad_val_)),
  b_nvmus_(tree_.Branch("nvmus", &nvmus_)),
  c_nvmus_(false),
  run_(static_cast<int >(bad_val_)),
  b_run_(tree_.Branch("run", &run_)),
  c_run_(false),
  els_ele105_(0),
  p_els_ele105_(&els_ele105_),
  b_els_ele105_(tree_.Branch("els_ele105", &p_els_ele105_)),
  c_els_ele105_(false),
  els_ele23_(0),
  p_els_ele23_(&els_ele23_),
  b_els_ele23_(tree_.Branch("els_ele23", &p_els_ele23_)),
  c_els_ele23_(false),
  els_ele8_(0),
  p_els_ele8_(&els_ele8_),
  b_els_ele8_(tree_.Branch("els_ele8", &p_els_ele8_)),
  c_els_ele8_(false),
  els_inz_(0),
  p_els_inz_(&els_inz_),
  b_els_inz_(tree_.Branch("els_inz", &p_els_inz_)),
  c_els_inz_(false),
  els_inzv_(0),
  p_els_inzv_(&els_inzv_),
  b_els_inzv_(tree_.Branch("els_inzv", &p_els_inzv_)),
  c_els_inzv_(false),
  els_ispf_(0),
  p_els_ispf_(&els_ispf_),
  b_els_ispf_(tree_.Branch("els_ispf", &p_els_ispf_)),
  c_els_ispf_(false),
  els_sig_(0),
  p_els_sig_(&els_sig_),
  b_els_sig_(tree_.Branch("els_sig", &p_els_sig_)),
  c_els_sig_(false),
  els_sigid_(0),
  p_els_sigid_(&els_sigid_),
  b_els_sigid_(tree_.Branch("els_sigid", &p_els_sigid_)),
  c_els_sigid_(false),
  els_tight_(0),
  p_els_tight_(&els_tight_),
  b_els_tight_(tree_.Branch("els_tight", &p_els_tight_)),
  c_els_tight_(false),
  els_tm_(0),
  p_els_tm_(&els_tm_),
  b_els_tm_(tree_.Branch("els_tm", &p_els_tm_)),
  c_els_tm_(false),
  els_vvvl_(0),
  p_els_vvvl_(&els_vvvl_),
  b_els_vvvl_(tree_.Branch("els_vvvl", &p_els_vvvl_)),
  c_els_vvvl_(false),
  jets_islep_(0),
  p_jets_islep_(&jets_islep_),
  b_jets_islep_(tree_.Branch("jets_islep", &p_jets_islep_)),
  c_jets_islep_(false),
  mus_inz_(0),
  p_mus_inz_(&mus_inz_),
  b_mus_inz_(tree_.Branch("mus_inz", &p_mus_inz_)),
  c_mus_inz_(false),
  mus_inzv_(0),
  p_mus_inzv_(&mus_inzv_),
  b_mus_inzv_(tree_.Branch("mus_inzv", &p_mus_inzv_)),
  c_mus_inzv_(false),
  mus_isomu18_(0),
  p_mus_isomu18_(&mus_isomu18_),
  b_mus_isomu18_(tree_.Branch("mus_isomu18", &p_mus_isomu18_)),
  c_mus_isomu18_(false),
  mus_mu50_(0),
  p_mus_mu50_(&mus_mu50_),
  b_mus_mu50_(tree_.Branch("mus_mu50", &p_mus_mu50_)),
  c_mus_mu50_(false),
  mus_mu8_(0),
  p_mus_mu8_(&mus_mu8_),
  b_mus_mu8_(tree_.Branch("mus_mu8", &p_mus_mu8_)),
  c_mus_mu8_(false),
  mus_sig_(0),
  p_mus_sig_(&mus_sig_),
  b_mus_sig_(tree_.Branch("mus_sig", &p_mus_sig_)),
  c_mus_sig_(false),
  mus_sigid_(0),
  p_mus_sigid_(&mus_sigid_),
  b_mus_sigid_(tree_.Branch("mus_sigid", &p_mus_sigid_)),
  c_mus_sigid_(false),
  mus_tight_(0),
  p_mus_tight_(&mus_tight_),
  b_mus_tight_(tree_.Branch("mus_tight", &p_mus_tight_)),
  c_mus_tight_(false),
  mus_tm_(0),
  p_mus_tm_(&mus_tm_),
  b_mus_tm_(tree_.Branch("mus_tm", &p_mus_tm_)),
  c_mus_tm_(false),
  mus_vvvl_(0),
  p_mus_vvvl_(&mus_vvvl_),
  b_mus_vvvl_(tree_.Branch("mus_vvvl", &p_mus_vvvl_)),
  c_mus_vvvl_(false),
  ph_ph90_(0),
  p_ph_ph90_(&ph_ph90_),
  b_ph_ph90_(tree_.Branch("ph_ph90", &p_ph_ph90_)),
  c_ph_ph90_(false),
  ph_tm_(0),
  p_ph_tm_(&ph_tm_),
  b_ph_tm_(tree_.Branch("ph_tm", &p_ph_tm_)),
  c_ph_tm_(false),
  sys_pass_(0),
  p_sys_pass_(&sys_pass_),
  b_sys_pass_(tree_.Branch("sys_pass", &p_sys_pass_)),
  c_sys_pass_(false),
  trig_(0),
  p_trig_(&trig_),
  b_trig_(tree_.Branch("trig", &p_trig_)),
  c_trig_(false),
  els_d0_(0),
  p_els_d0_(&els_d0_),
  b_els_d0_(tree_.Branch("els_d0", &p_els_d0_)),
  c_els_d0_(false),
  els_dz_(0),
  p_els_dz_(&els_dz_),
  b_els_dz_(tree_.Branch("els_dz", &p_els_dz_)),
  c_els_dz_(false),
  els_eta_(0),
  p_els_eta_(&els_eta_),
  b_els_eta_(tree_.Branch("els_eta", &p_els_eta_)),
  c_els_eta_(false),
  els_miniso_(0),
  p_els_miniso_(&els_miniso_),
  b_els_miniso_(tree_.Branch("els_miniso", &p_els_miniso_)),
  c_els_miniso_(false),
  els_phi_(0),
  p_els_phi_(&els_phi_),
  b_els_phi_(tree_.Branch("els_phi", &p_els_phi_)),
  c_els_phi_(false),
  els_pt_(0),
  p_els_pt_(&els_pt_),
  b_els_pt_(tree_.Branch("els_pt", &p_els_pt_)),
  c_els_pt_(false),
  els_reliso_(0),
  p_els_reliso_(&els_reliso_),
  b_els_reliso_(tree_.Branch("els_reliso", &p_els_reliso_)),
  c_els_reliso_(false),
  els_sceta_(0),
  p_els_sceta_(&els_sceta_),
  b_els_sceta_(tree_.Branch("els_sceta", &p_els_sceta_)),
  c_els_sceta_(false),
  els_vvvl_eta_(0),
  p_els_vvvl_eta_(&els_vvvl_eta_),
  b_els_vvvl_eta_(tree_.Branch("els_vvvl_eta", &p_els_vvvl_eta_)),
  c_els_vvvl_eta_(false),
  els_vvvl_phi_(0),
  p_els_vvvl_phi_(&els_vvvl_phi_),
  b_els_vvvl_phi_(tree_.Branch("els_vvvl_phi", &p_els_vvvl_phi_)),
  c_els_vvvl_phi_(false),
  els_vvvl_pt_(0),
  p_els_vvvl_pt_(&els_vvvl_pt_),
  b_els_vvvl_pt_(tree_.Branch("els_vvvl_pt", &p_els_vvvl_pt_)),
  c_els_vvvl_pt_(false),
  fjets08_eta_(0),
  p_fjets08_eta_(&fjets08_eta_),
  b_fjets08_eta_(tree_.Branch("fjets08_eta", &p_fjets08_eta_)),
  c_fjets08_eta_(false),
  fjets08_m_(0),
  p_fjets08_m_(&fjets08_m_),
  b_fjets08_m_(tree_.Branch("fjets08_m", &p_fjets08_m_)),
  c_fjets08_m_(false),
  fjets08_phi_(0),
  p_fjets08_phi_(&fjets08_phi_),
  b_fjets08_phi_(tree_.Branch("fjets08_phi", &p_fjets08_phi_)),
  c_fjets08_phi_(false),
  fjets08_poscsv_(0),
  p_fjets08_poscsv_(&fjets08_poscsv_),
  b_fjets08_poscsv_(tree_.Branch("fjets08_poscsv", &p_fjets08_poscsv_)),
  c_fjets08_poscsv_(false),
  fjets08_pt_(0),
  p_fjets08_pt_(&fjets08_pt_),
  b_fjets08_pt_(tree_.Branch("fjets08_pt", &p_fjets08_pt_)),
  c_fjets08_pt_(false),
  fjets08_sumcsv_(0),
  p_fjets08_sumcsv_(&fjets08_sumcsv_),
  b_fjets08_sumcsv_(tree_.Branch("fjets08_sumcsv", &p_fjets08_sumcsv_)),
  c_fjets08_sumcsv_(false),
  fjets_eta_(0),
  p_fjets_eta_(&fjets_eta_),
  b_fjets_eta_(tree_.Branch("fjets_eta", &p_fjets_eta_)),
  c_fjets_eta_(false),
  fjets_m_(0),
  p_fjets_m_(&fjets_m_),
  b_fjets_m_(tree_.Branch("fjets_m", &p_fjets_m_)),
  c_fjets_m_(false),
  fjets_phi_(0),
  p_fjets_phi_(&fjets_phi_),
  b_fjets_phi_(tree_.Branch("fjets_phi", &p_fjets_phi_)),
  c_fjets_phi_(false),
  fjets_poscsv_(0),
  p_fjets_poscsv_(&fjets_poscsv_),
  b_fjets_poscsv_(tree_.Branch("fjets_poscsv", &p_fjets_poscsv_)),
  c_fjets_poscsv_(false),
  fjets_pt_(0),
  p_fjets_pt_(&fjets_pt_),
  b_fjets_pt_(tree_.Branch("fjets_pt", &p_fjets_pt_)),
  c_fjets_pt_(false),
  fjets_sumcsv_(0),
  p_fjets_sumcsv_(&fjets_sumcsv_),
  b_fjets_sumcsv_(tree_.Branch("fjets_sumcsv", &p_fjets_sumcsv_)),
  c_fjets_sumcsv_(false),
  jets_csv_(0),
  p_jets_csv_(&jets_csv_),
  b_jets_csv_(tree_.Branch("jets_csv", &p_jets_csv_)),
  c_jets_csv_(false),
  jets_eta_(0),
  p_jets_eta_(&jets_eta_),
  b_jets_eta_(tree_.Branch("jets_eta", &p_jets_eta_)),
  c_jets_eta_(false),
  jets_m_(0),
  p_jets_m_(&jets_m_),
  b_jets_m_(tree_.Branch("jets_m", &p_jets_m_)),
  c_jets_m_(false),
  jets_phi_(0),
  p_jets_phi_(&jets_phi_),
  b_jets_phi_(tree_.Branch("jets_phi", &p_jets_phi_)),
  c_jets_phi_(false),
  jets_pt_(0),
  p_jets_pt_(&jets_pt_),
  b_jets_pt_(tree_.Branch("jets_pt", &p_jets_pt_)),
  c_jets_pt_(false),
  jets_pt_res_(0),
  p_jets_pt_res_(&jets_pt_res_),
  b_jets_pt_res_(tree_.Branch("jets_pt_res", &p_jets_pt_res_)),
  c_jets_pt_res_(false),
  leps_eta_(0),
  p_leps_eta_(&leps_eta_),
  b_leps_eta_(tree_.Branch("leps_eta", &p_leps_eta_)),
  c_leps_eta_(false),
  leps_id_(0),
  p_leps_id_(&leps_id_),
  b_leps_id_(tree_.Branch("leps_id", &p_leps_id_)),
  c_leps_id_(false),
  leps_phi_(0),
  p_leps_phi_(&leps_phi_),
  b_leps_phi_(tree_.Branch("leps_phi", &p_leps_phi_)),
  c_leps_phi_(false),
  leps_pt_(0),
  p_leps_pt_(&leps_pt_),
  b_leps_pt_(tree_.Branch("leps_pt", &p_leps_pt_)),
  c_leps_pt_(false),
  mc_eta_(0),
  p_mc_eta_(&mc_eta_),
  b_mc_eta_(tree_.Branch("mc_eta", &p_mc_eta_)),
  c_mc_eta_(false),
  mc_mass_(0),
  p_mc_mass_(&mc_mass_),
  b_mc_mass_(tree_.Branch("mc_mass", &p_mc_mass_)),
  c_mc_mass_(false),
  mc_phi_(0),
  p_mc_phi_(&mc_phi_),
  b_mc_phi_(tree_.Branch("mc_phi", &p_mc_phi_)),
  c_mc_phi_(false),
  mc_pt_(0),
  p_mc_pt_(&mc_pt_),
  b_mc_pt_(tree_.Branch("mc_pt", &p_mc_pt_)),
  c_mc_pt_(false),
  mus_d0_(0),
  p_mus_d0_(&mus_d0_),
  b_mus_d0_(tree_.Branch("mus_d0", &p_mus_d0_)),
  c_mus_d0_(false),
  mus_dz_(0),
  p_mus_dz_(&mus_dz_),
  b_mus_dz_(tree_.Branch("mus_dz", &p_mus_dz_)),
  c_mus_dz_(false),
  mus_eta_(0),
  p_mus_eta_(&mus_eta_),
  b_mus_eta_(tree_.Branch("mus_eta", &p_mus_eta_)),
  c_mus_eta_(false),
  mus_miniso_(0),
  p_mus_miniso_(&mus_miniso_),
  b_mus_miniso_(tree_.Branch("mus_miniso", &p_mus_miniso_)),
  c_mus_miniso_(false),
  mus_phi_(0),
  p_mus_phi_(&mus_phi_),
  b_mus_phi_(tree_.Branch("mus_phi", &p_mus_phi_)),
  c_mus_phi_(false),
  mus_pt_(0),
  p_mus_pt_(&mus_pt_),
  b_mus_pt_(tree_.Branch("mus_pt", &p_mus_pt_)),
  c_mus_pt_(false),
  mus_reliso_(0),
  p_mus_reliso_(&mus_reliso_),
  b_mus_reliso_(tree_.Branch("mus_reliso", &p_mus_reliso_)),
  c_mus_reliso_(false),
  mus_vvvl_eta_(0),
  p_mus_vvvl_eta_(&mus_vvvl_eta_),
  b_mus_vvvl_eta_(tree_.Branch("mus_vvvl_eta", &p_mus_vvvl_eta_)),
  c_mus_vvvl_eta_(false),
  mus_vvvl_phi_(0),
  p_mus_vvvl_phi_(&mus_vvvl_phi_),
  b_mus_vvvl_phi_(tree_.Branch("mus_vvvl_phi", &p_mus_vvvl_phi_)),
  c_mus_vvvl_phi_(false),
  mus_vvvl_pt_(0),
  p_mus_vvvl_pt_(&mus_vvvl_pt_),
  b_mus_vvvl_pt_(tree_.Branch("mus_vvvl_pt", &p_mus_vvvl_pt_)),
  c_mus_vvvl_pt_(false),
  ph_eta_(0),
  p_ph_eta_(&ph_eta_),
  b_ph_eta_(tree_.Branch("ph_eta", &p_ph_eta_)),
  c_ph_eta_(false),
  ph_phi_(0),
  p_ph_phi_(&ph_phi_),
  b_ph_phi_(tree_.Branch("ph_phi", &p_ph_phi_)),
  c_ph_phi_(false),
  ph_pt_(0),
  p_ph_pt_(&ph_pt_),
  b_ph_pt_(tree_.Branch("ph_pt", &p_ph_pt_)),
  c_ph_pt_(false),
  sys_bctag_(0),
  p_sys_bctag_(&sys_bctag_),
  b_sys_bctag_(tree_.Branch("sys_bctag", &p_sys_bctag_)),
  c_sys_bctag_(false),
  sys_fs_bctag_(0),
  p_sys_fs_bctag_(&sys_fs_bctag_),
  b_sys_fs_bctag_(tree_.Branch("sys_fs_bctag", &p_sys_fs_bctag_)),
  c_sys_fs_bctag_(false),
  sys_fs_lep_(0),
  p_sys_fs_lep_(&sys_fs_lep_),
  b_sys_fs_lep_(tree_.Branch("sys_fs_lep", &p_sys_fs_lep_)),
  c_sys_fs_lep_(false),
  sys_fs_udsgtag_(0),
  p_sys_fs_udsgtag_(&sys_fs_udsgtag_),
  b_sys_fs_udsgtag_(tree_.Branch("sys_fs_udsgtag", &p_sys_fs_udsgtag_)),
  c_sys_fs_udsgtag_(false),
  sys_ht_(0),
  p_sys_ht_(&sys_ht_),
  b_sys_ht_(tree_.Branch("sys_ht", &p_sys_ht_)),
  c_sys_ht_(false),
  sys_isr_(0),
  p_sys_isr_(&sys_isr_),
  b_sys_isr_(tree_.Branch("sys_isr", &p_sys_isr_)),
  c_sys_isr_(false),
  sys_lep_(0),
  p_sys_lep_(&sys_lep_),
  b_sys_lep_(tree_.Branch("sys_lep", &p_sys_lep_)),
  c_sys_lep_(false),
  sys_met_(0),
  p_sys_met_(&sys_met_),
  b_sys_met_(tree_.Branch("sys_met", &p_sys_met_)),
  c_sys_met_(false),
  sys_mj_(0),
  p_sys_mj_(&sys_mj_),
  b_sys_mj_(tree_.Branch("sys_mj", &p_sys_mj_)),
  c_sys_mj_(false),
  sys_mt_(0),
  p_sys_mt_(&sys_mt_),
  b_sys_mt_(tree_.Branch("sys_mt", &p_sys_mt_)),
  c_sys_mt_(false),
  sys_muf_(0),
  p_sys_muf_(&sys_muf_),
  b_sys_muf_(tree_.Branch("sys_muf", &p_sys_muf_)),
  c_sys_muf_(false),
  sys_mur_(0),
  p_sys_mur_(&sys_mur_),
  b_sys_mur_(tree_.Branch("sys_mur", &p_sys_mur_)),
  c_sys_mur_(false),
  sys_murf_(0),
  p_sys_murf_(&sys_murf_),
  b_sys_murf_(tree_.Branch("sys_murf", &p_sys_murf_)),
  c_sys_murf_(false),
  sys_pdf_(0),
  p_sys_pdf_(&sys_pdf_),
  b_sys_pdf_(tree_.Branch("sys_pdf", &p_sys_pdf_)),
  c_sys_pdf_(false),
  sys_trig_(0),
  p_sys_trig_(&sys_trig_),
  b_sys_trig_(tree_.Branch("sys_trig", &p_sys_trig_)),
  c_sys_trig_(false),
  sys_udsgtag_(0),
  p_sys_udsgtag_(&sys_udsgtag_),
  b_sys_udsgtag_(tree_.Branch("sys_udsgtag", &p_sys_udsgtag_)),
  c_sys_udsgtag_(false),
  tks_eta_(0),
  p_tks_eta_(&tks_eta_),
  b_tks_eta_(tree_.Branch("tks_eta", &p_tks_eta_)),
  c_tks_eta_(false),
  tks_miniso_(0),
  p_tks_miniso_(&tks_miniso_),
  b_tks_miniso_(tree_.Branch("tks_miniso", &p_tks_miniso_)),
  c_tks_miniso_(false),
  tks_mt_(0),
  p_tks_mt_(&tks_mt_),
  b_tks_mt_(tree_.Branch("tks_mt", &p_tks_mt_)),
  c_tks_mt_(false),
  tks_mt2_(0),
  p_tks_mt2_(&tks_mt2_),
  b_tks_mt2_(tree_.Branch("tks_mt2", &p_tks_mt2_)),
  c_tks_mt2_(false),
  tks_phi_(0),
  p_tks_phi_(&tks_phi_),
  b_tks_phi_(tree_.Branch("tks_phi", &p_tks_phi_)),
  c_tks_phi_(false),
  tks_pt_(0),
  p_tks_pt_(&tks_pt_),
  b_tks_pt_(tree_.Branch("tks_pt", &p_tks_pt_)),
  c_tks_pt_(false),
  trig_prescale_(0),
  p_trig_prescale_(&trig_prescale_),
  b_trig_prescale_(tree_.Branch("trig_prescale", &p_trig_prescale_)),
  c_trig_prescale_(false),
  w_pdf_(0),
  p_w_pdf_(&w_pdf_),
  b_w_pdf_(tree_.Branch("w_pdf", &p_w_pdf_)),
  c_w_pdf_(false),
  els_charge_(0),
  p_els_charge_(&els_charge_),
  b_els_charge_(tree_.Branch("els_charge", &p_els_charge_)),
  c_els_charge_(false),
  fjets08_btags_(0),
  p_fjets08_btags_(&fjets08_btags_),
  b_fjets08_btags_(tree_.Branch("fjets08_btags", &p_fjets08_btags_)),
  c_fjets08_btags_(false),
  fjets08_nconst_(0),
  p_fjets08_nconst_(&fjets08_nconst_),
  b_fjets08_nconst_(tree_.Branch("fjets08_nconst", &p_fjets08_nconst_)),
  c_fjets08_nconst_(false),
  fjets_btags_(0),
  p_fjets_btags_(&fjets_btags_),
  b_fjets_btags_(tree_.Branch("fjets_btags", &p_fjets_btags_)),
  c_fjets_btags_(false),
  fjets_nconst_(0),
  p_fjets_nconst_(&fjets_nconst_),
  b_fjets_nconst_(tree_.Branch("fjets_nconst", &p_fjets_nconst_)),
  c_fjets_nconst_(false),
  jets_fjet08_index_(0),
  p_jets_fjet08_index_(&jets_fjet08_index_),
  b_jets_fjet08_index_(tree_.Branch("jets_fjet08_index", &p_jets_fjet08_index_)),
  c_jets_fjet08_index_(false),
  jets_fjet_index_(0),
  p_jets_fjet_index_(&jets_fjet_index_),
  b_jets_fjet_index_(tree_.Branch("jets_fjet_index", &p_jets_fjet_index_)),
  c_jets_fjet_index_(false),
  jets_hflavor_(0),
  p_jets_hflavor_(&jets_hflavor_),
  b_jets_hflavor_(tree_.Branch("jets_hflavor", &p_jets_hflavor_)),
  c_jets_hflavor_(false),
  mc_id_(0),
  p_mc_id_(&mc_id_),
  b_mc_id_(tree_.Branch("mc_id", &p_mc_id_)),
  c_mc_id_(false),
  mc_mom_(0),
  p_mc_mom_(&mc_mom_),
  b_mc_mom_(tree_.Branch("mc_mom", &p_mc_mom_)),
  c_mc_mom_(false),
  mus_charge_(0),
  p_mus_charge_(&mus_charge_),
  b_mus_charge_(tree_.Branch("mus_charge", &p_mus_charge_)),
  c_mus_charge_(false),
  sys_nbm_(0),
  p_sys_nbm_(&sys_nbm_),
  b_sys_nbm_(tree_.Branch("sys_nbm", &p_sys_nbm_)),
  c_sys_nbm_(false),
  sys_njets_(0),
  p_sys_njets_(&sys_njets_),
  b_sys_njets_(tree_.Branch("sys_njets", &p_sys_njets_)),
  c_sys_njets_(false),
  tks_pdg_(0),
  b_tks_pdg_(tree_.Branch("tks_pdg", &tks_pdg_)),
  c_tks_pdg_(false){
}

baby_base::baby_base(const string &filename):
  chain_("tree","tree"),
  tree_("junk","junk"),
  entry_(0),
  read_only_(true),
  event_(static_cast<Long64_t >(bad_val_)),
  b_event_(NULL),
  c_event_(false),
  fromGS_(static_cast<bool >(bad_val_)),
  b_fromGS_(NULL),
  c_fromGS_(false),
  jetmismeas_(static_cast<bool >(bad_val_)),
  b_jetmismeas_(NULL),
  c_jetmismeas_(false),
  low_dphi_(static_cast<bool >(bad_val_)),
  b_low_dphi_(NULL),
  c_low_dphi_(false),
  nonblind_(static_cast<bool >(bad_val_)),
  b_nonblind_(NULL),
  c_nonblind_(false),
  pass_(static_cast<bool >(bad_val_)),
  b_pass_(NULL),
  c_pass_(false),
  pass_cschalo_(static_cast<bool >(bad_val_)),
  b_pass_cschalo_(NULL),
  c_pass_cschalo_(false),
  pass_eebadsc_(static_cast<bool >(bad_val_)),
  b_pass_eebadsc_(NULL),
  c_pass_eebadsc_(false),
  pass_goodv_(static_cast<bool >(bad_val_)),
  b_pass_goodv_(NULL),
  c_pass_goodv_(false),
  pass_hbhe_(static_cast<bool >(bad_val_)),
  b_pass_hbhe_(NULL),
  c_pass_hbhe_(false),
  pass_hbheiso_(static_cast<bool >(bad_val_)),
  b_pass_hbheiso_(NULL),
  c_pass_hbheiso_(false),
  pass_jets_(static_cast<bool >(bad_val_)),
  b_pass_jets_(NULL),
  c_pass_jets_(false),
  pass_jets_nohf_(static_cast<bool >(bad_val_)),
  b_pass_jets_nohf_(NULL),
  c_pass_jets_nohf_(false),
  pass_jets_ra2_(static_cast<bool >(bad_val_)),
  b_pass_jets_ra2_(NULL),
  c_pass_jets_ra2_(false),
  pass_jets_tight_(static_cast<bool >(bad_val_)),
  b_pass_jets_tight_(NULL),
  c_pass_jets_tight_(false),
  pass_jets_tight_ra2_(static_cast<bool >(bad_val_)),
  b_pass_jets_tight_ra2_(NULL),
  c_pass_jets_tight_ra2_(false),
  pass_nohf_(static_cast<bool >(bad_val_)),
  b_pass_nohf_(NULL),
  c_pass_nohf_(false),
  pass_ra2_(static_cast<bool >(bad_val_)),
  b_pass_ra2_(NULL),
  c_pass_ra2_(false),
  stitch_(static_cast<bool >(bad_val_)),
  b_stitch_(NULL),
  c_stitch_(false),
  dphi1_(static_cast<float >(bad_val_)),
  b_dphi1_(NULL),
  c_dphi1_(false),
  dphi2_(static_cast<float >(bad_val_)),
  b_dphi2_(NULL),
  c_dphi2_(false),
  dphi3_(static_cast<float >(bad_val_)),
  b_dphi3_(NULL),
  c_dphi3_(false),
  dphi4_(static_cast<float >(bad_val_)),
  b_dphi4_(NULL),
  c_dphi4_(false),
  dphi_wlep_(static_cast<float >(bad_val_)),
  b_dphi_wlep_(NULL),
  c_dphi_wlep_(false),
  eff_trig_(static_cast<float >(bad_val_)),
  b_eff_trig_(NULL),
  c_eff_trig_(false),
  elel_eta_(static_cast<float >(bad_val_)),
  b_elel_eta_(NULL),
  c_elel_eta_(false),
  elel_m_(static_cast<float >(bad_val_)),
  b_elel_m_(NULL),
  c_elel_m_(false),
  elel_phi_(static_cast<float >(bad_val_)),
  b_elel_phi_(NULL),
  c_elel_phi_(false),
  elel_pt_(static_cast<float >(bad_val_)),
  b_elel_pt_(NULL),
  c_elel_pt_(false),
  elel_pt1_(static_cast<float >(bad_val_)),
  b_elel_pt1_(NULL),
  c_elel_pt1_(false),
  elel_pt2_(static_cast<float >(bad_val_)),
  b_elel_pt2_(NULL),
  c_elel_pt2_(false),
  elel_w_(static_cast<float >(bad_val_)),
  b_elel_w_(NULL),
  c_elel_w_(false),
  elelv_eta_(static_cast<float >(bad_val_)),
  b_elelv_eta_(NULL),
  c_elelv_eta_(false),
  elelv_m_(static_cast<float >(bad_val_)),
  b_elelv_m_(NULL),
  c_elelv_m_(false),
  elelv_phi_(static_cast<float >(bad_val_)),
  b_elelv_phi_(NULL),
  c_elelv_phi_(false),
  elelv_pt_(static_cast<float >(bad_val_)),
  b_elelv_pt_(NULL),
  c_elelv_pt_(false),
  elelv_pt1_(static_cast<float >(bad_val_)),
  b_elelv_pt1_(NULL),
  c_elelv_pt1_(false),
  elelv_pt2_(static_cast<float >(bad_val_)),
  b_elelv_pt2_(NULL),
  c_elelv_pt2_(false),
  elelv_w_(static_cast<float >(bad_val_)),
  b_elelv_w_(NULL),
  c_elelv_w_(false),
  elmu_eta_(static_cast<float >(bad_val_)),
  b_elmu_eta_(NULL),
  c_elmu_eta_(false),
  elmu_m_(static_cast<float >(bad_val_)),
  b_elmu_m_(NULL),
  c_elmu_m_(false),
  elmu_phi_(static_cast<float >(bad_val_)),
  b_elmu_phi_(NULL),
  c_elmu_phi_(false),
  elmu_pt_(static_cast<float >(bad_val_)),
  b_elmu_pt_(NULL),
  c_elmu_pt_(false),
  elmu_pt1_(static_cast<float >(bad_val_)),
  b_elmu_pt1_(NULL),
  c_elmu_pt1_(false),
  elmu_pt2_(static_cast<float >(bad_val_)),
  b_elmu_pt2_(NULL),
  c_elmu_pt2_(false),
  elmu_w_(static_cast<float >(bad_val_)),
  b_elmu_w_(NULL),
  c_elmu_w_(false),
  ht_(static_cast<float >(bad_val_)),
  b_ht_(NULL),
  c_ht_(false),
  ht_clean_(static_cast<float >(bad_val_)),
  b_ht_clean_(NULL),
  c_ht_clean_(false),
  ht_hlt_(static_cast<float >(bad_val_)),
  b_ht_hlt_(NULL),
  c_ht_hlt_(false),
  ht_isr_me_(static_cast<float >(bad_val_)),
  b_ht_isr_me_(NULL),
  c_ht_isr_me_(false),
  ht_ra2_(static_cast<float >(bad_val_)),
  b_ht_ra2_(NULL),
  c_ht_ra2_(false),
  ht_tru_(static_cast<float >(bad_val_)),
  b_ht_tru_(NULL),
  c_ht_tru_(false),
  isr_tru_eta_(static_cast<float >(bad_val_)),
  b_isr_tru_eta_(NULL),
  c_isr_tru_eta_(false),
  isr_tru_phi_(static_cast<float >(bad_val_)),
  b_isr_tru_phi_(NULL),
  c_isr_tru_phi_(false),
  isr_tru_pt_(static_cast<float >(bad_val_)),
  b_isr_tru_pt_(NULL),
  c_isr_tru_pt_(false),
  jetsys_eta_(static_cast<float >(bad_val_)),
  b_jetsys_eta_(NULL),
  c_jetsys_eta_(false),
  jetsys_nob_eta_(static_cast<float >(bad_val_)),
  b_jetsys_nob_eta_(NULL),
  c_jetsys_nob_eta_(false),
  jetsys_nob_phi_(static_cast<float >(bad_val_)),
  b_jetsys_nob_phi_(NULL),
  c_jetsys_nob_phi_(false),
  jetsys_nob_pt_(static_cast<float >(bad_val_)),
  b_jetsys_nob_pt_(NULL),
  c_jetsys_nob_pt_(false),
  jetsys_phi_(static_cast<float >(bad_val_)),
  b_jetsys_phi_(NULL),
  c_jetsys_phi_(false),
  jetsys_pt_(static_cast<float >(bad_val_)),
  b_jetsys_pt_(NULL),
  c_jetsys_pt_(false),
  met_(static_cast<float >(bad_val_)),
  b_met_(NULL),
  c_met_(false),
  met_calo_(static_cast<float >(bad_val_)),
  b_met_calo_(NULL),
  c_met_calo_(false),
  met_calo_phi_(static_cast<float >(bad_val_)),
  b_met_calo_phi_(NULL),
  c_met_calo_phi_(false),
  met_mini_(static_cast<float >(bad_val_)),
  b_met_mini_(NULL),
  c_met_mini_(false),
  met_mini_phi_(static_cast<float >(bad_val_)),
  b_met_mini_phi_(NULL),
  c_met_mini_phi_(false),
  met_nohf_(static_cast<float >(bad_val_)),
  b_met_nohf_(NULL),
  c_met_nohf_(false),
  met_nohf_phi_(static_cast<float >(bad_val_)),
  b_met_nohf_phi_(NULL),
  c_met_nohf_phi_(false),
  met_phi_(static_cast<float >(bad_val_)),
  b_met_phi_(NULL),
  c_met_phi_(false),
  met_raw_(static_cast<float >(bad_val_)),
  b_met_raw_(NULL),
  c_met_raw_(false),
  met_raw_phi_(static_cast<float >(bad_val_)),
  b_met_raw_phi_(NULL),
  c_met_raw_phi_(false),
  met_rebal_(static_cast<float >(bad_val_)),
  b_met_rebal_(NULL),
  c_met_rebal_(false),
  met_tru_(static_cast<float >(bad_val_)),
  b_met_tru_(NULL),
  c_met_tru_(false),
  met_tru_nuw_(static_cast<float >(bad_val_)),
  b_met_tru_nuw_(NULL),
  c_met_tru_nuw_(false),
  met_tru_nuw_phi_(static_cast<float >(bad_val_)),
  b_met_tru_nuw_phi_(NULL),
  c_met_tru_nuw_phi_(false),
  met_tru_phi_(static_cast<float >(bad_val_)),
  b_met_tru_phi_(NULL),
  c_met_tru_phi_(false),
  mht_(static_cast<float >(bad_val_)),
  b_mht_(NULL),
  c_mht_(false),
  mht_clean_(static_cast<float >(bad_val_)),
  b_mht_clean_(NULL),
  c_mht_clean_(false),
  mht_clean_phi_(static_cast<float >(bad_val_)),
  b_mht_clean_phi_(NULL),
  c_mht_clean_phi_(false),
  mht_phi_(static_cast<float >(bad_val_)),
  b_mht_phi_(NULL),
  c_mht_phi_(false),
  mj_(static_cast<float >(bad_val_)),
  b_mj_(NULL),
  c_mj_(false),
  mj08_(static_cast<float >(bad_val_)),
  b_mj08_(NULL),
  c_mj08_(false),
  mt_(static_cast<float >(bad_val_)),
  b_mt_(NULL),
  c_mt_(false),
  mt_nohf_(static_cast<float >(bad_val_)),
  b_mt_nohf_(NULL),
  c_mt_nohf_(false),
  mt_rebal_(static_cast<float >(bad_val_)),
  b_mt_rebal_(NULL),
  c_mt_rebal_(false),
  mt_tru_(static_cast<float >(bad_val_)),
  b_mt_tru_(NULL),
  c_mt_tru_(false),
  mt_tru_nuw_(static_cast<float >(bad_val_)),
  b_mt_tru_nuw_(NULL),
  c_mt_tru_nuw_(false),
  mumu_eta_(static_cast<float >(bad_val_)),
  b_mumu_eta_(NULL),
  c_mumu_eta_(false),
  mumu_m_(static_cast<float >(bad_val_)),
  b_mumu_m_(NULL),
  c_mumu_m_(false),
  mumu_phi_(static_cast<float >(bad_val_)),
  b_mumu_phi_(NULL),
  c_mumu_phi_(false),
  mumu_pt_(static_cast<float >(bad_val_)),
  b_mumu_pt_(NULL),
  c_mumu_pt_(false),
  mumu_pt1_(static_cast<float >(bad_val_)),
  b_mumu_pt1_(NULL),
  c_mumu_pt1_(false),
  mumu_pt2_(static_cast<float >(bad_val_)),
  b_mumu_pt2_(NULL),
  c_mumu_pt2_(false),
  mumu_w_(static_cast<float >(bad_val_)),
  b_mumu_w_(NULL),
  c_mumu_w_(false),
  mumuv_eta_(static_cast<float >(bad_val_)),
  b_mumuv_eta_(NULL),
  c_mumuv_eta_(false),
  mumuv_m_(static_cast<float >(bad_val_)),
  b_mumuv_m_(NULL),
  c_mumuv_m_(false),
  mumuv_phi_(static_cast<float >(bad_val_)),
  b_mumuv_phi_(NULL),
  c_mumuv_phi_(false),
  mumuv_pt_(static_cast<float >(bad_val_)),
  b_mumuv_pt_(NULL),
  c_mumuv_pt_(false),
  mumuv_pt1_(static_cast<float >(bad_val_)),
  b_mumuv_pt1_(NULL),
  c_mumuv_pt1_(false),
  mumuv_pt2_(static_cast<float >(bad_val_)),
  b_mumuv_pt2_(NULL),
  c_mumuv_pt2_(false),
  mumuv_w_(static_cast<float >(bad_val_)),
  b_mumuv_w_(NULL),
  c_mumuv_w_(false),
  ntrupv_mean_(static_cast<float >(bad_val_)),
  b_ntrupv_mean_(NULL),
  c_ntrupv_mean_(false),
  onel_ele105_(static_cast<float >(bad_val_)),
  b_onel_ele105_(NULL),
  c_onel_ele105_(false),
  onel_ele23_(static_cast<float >(bad_val_)),
  b_onel_ele23_(NULL),
  c_onel_ele23_(false),
  onel_ele8_(static_cast<float >(bad_val_)),
  b_onel_ele8_(NULL),
  c_onel_ele8_(false),
  onel_vvvl_(static_cast<float >(bad_val_)),
  b_onel_vvvl_(NULL),
  c_onel_vvvl_(false),
  onht_(static_cast<float >(bad_val_)),
  b_onht_(NULL),
  c_onht_(false),
  onmet_(static_cast<float >(bad_val_)),
  b_onmet_(NULL),
  c_onmet_(false),
  onmu_isomu18_(static_cast<float >(bad_val_)),
  b_onmu_isomu18_(NULL),
  c_onmu_isomu18_(false),
  onmu_mu50_(static_cast<float >(bad_val_)),
  b_onmu_mu50_(NULL),
  c_onmu_mu50_(false),
  onmu_mu8_(static_cast<float >(bad_val_)),
  b_onmu_mu8_(NULL),
  c_onmu_mu8_(false),
  onmu_vvvl_(static_cast<float >(bad_val_)),
  b_onmu_vvvl_(NULL),
  c_onmu_vvvl_(false),
  onph_ph90_(static_cast<float >(bad_val_)),
  b_onph_ph90_(NULL),
  c_onph_ph90_(false),
  w_btag_(static_cast<float >(bad_val_)),
  b_w_btag_(NULL),
  c_w_btag_(false),
  w_fs_lep_(static_cast<float >(bad_val_)),
  b_w_fs_lep_(NULL),
  c_w_fs_lep_(false),
  w_lep_(static_cast<float >(bad_val_)),
  b_w_lep_(NULL),
  c_w_lep_(false),
  w_lumi_(static_cast<float >(bad_val_)),
  b_w_lumi_(NULL),
  c_w_lumi_(false),
  w_pu_(static_cast<float >(bad_val_)),
  b_w_pu_(NULL),
  c_w_pu_(false),
  w_toppt_(static_cast<float >(bad_val_)),
  b_w_toppt_(NULL),
  c_w_toppt_(false),
  weight_(static_cast<float >(bad_val_)),
  b_weight_(NULL),
  c_weight_(false),
  lumiblock_(static_cast<int >(bad_val_)),
  b_lumiblock_(NULL),
  c_lumiblock_(false),
  mgluino_(static_cast<int >(bad_val_)),
  b_mgluino_(NULL),
  c_mgluino_(false),
  mlsp_(static_cast<int >(bad_val_)),
  b_mlsp_(NULL),
  c_mlsp_(false),
  nbl_(static_cast<int >(bad_val_)),
  b_nbl_(NULL),
  c_nbl_(false),
  nbm_(static_cast<int >(bad_val_)),
  b_nbm_(NULL),
  c_nbm_(false),
  nbm_ra2_(static_cast<int >(bad_val_)),
  b_nbm_ra2_(NULL),
  c_nbm_ra2_(false),
  nbt_(static_cast<int >(bad_val_)),
  b_nbt_(NULL),
  c_nbt_(false),
  nels_(static_cast<int >(bad_val_)),
  b_nels_(NULL),
  c_nels_(false),
  nels_ele23_(static_cast<int >(bad_val_)),
  b_nels_ele23_(NULL),
  c_nels_ele23_(false),
  nels_vvvl_(static_cast<int >(bad_val_)),
  b_nels_vvvl_(NULL),
  c_nels_vvvl_(false),
  nfjets_(static_cast<int >(bad_val_)),
  b_nfjets_(NULL),
  c_nfjets_(false),
  nfjets08_(static_cast<int >(bad_val_)),
  b_nfjets08_(NULL),
  c_nfjets08_(false),
  nisr_me_(static_cast<int >(bad_val_)),
  b_nisr_me_(NULL),
  c_nisr_me_(false),
  njets_(static_cast<int >(bad_val_)),
  b_njets_(NULL),
  c_njets_(false),
  njets_clean_(static_cast<int >(bad_val_)),
  b_njets_clean_(NULL),
  c_njets_clean_(false),
  njets_ra2_(static_cast<int >(bad_val_)),
  b_njets_ra2_(NULL),
  c_njets_ra2_(false),
  nleps_(static_cast<int >(bad_val_)),
  b_nleps_(NULL),
  c_nleps_(false),
  nleps_tm_(static_cast<int >(bad_val_)),
  b_nleps_tm_(NULL),
  c_nleps_tm_(false),
  nmus_(static_cast<int >(bad_val_)),
  b_nmus_(NULL),
  c_nmus_(false),
  nmus_isomu18_(static_cast<int >(bad_val_)),
  b_nmus_isomu18_(NULL),
  c_nmus_isomu18_(false),
  nmus_vvvl_(static_cast<int >(bad_val_)),
  b_nmus_vvvl_(NULL),
  c_nmus_vvvl_(false),
  nph_(static_cast<int >(bad_val_)),
  b_nph_(NULL),
  c_nph_(false),
  npv_(static_cast<int >(bad_val_)),
  b_npv_(NULL),
  c_npv_(false),
  ntks_(static_cast<int >(bad_val_)),
  b_ntks_(NULL),
  c_ntks_(false),
  ntruels_(static_cast<int >(bad_val_)),
  b_ntruels_(NULL),
  c_ntruels_(false),
  ntruleps_(static_cast<int >(bad_val_)),
  b_ntruleps_(NULL),
  c_ntruleps_(false),
  ntrumus_(static_cast<int >(bad_val_)),
  b_ntrumus_(NULL),
  c_ntrumus_(false),
  ntrupv_(static_cast<int >(bad_val_)),
  b_ntrupv_(NULL),
  c_ntrupv_(false),
  ntrutaush_(static_cast<int >(bad_val_)),
  b_ntrutaush_(NULL),
  c_ntrutaush_(false),
  ntrutausl_(static_cast<int >(bad_val_)),
  b_ntrutausl_(NULL),
  c_ntrutausl_(false),
  nvels_(static_cast<int >(bad_val_)),
  b_nvels_(NULL),
  c_nvels_(false),
  nveto_(static_cast<int >(bad_val_)),
  b_nveto_(NULL),
  c_nveto_(false),
  nvleps_(static_cast<int >(bad_val_)),
  b_nvleps_(NULL),
  c_nvleps_(false),
  nvmus_(static_cast<int >(bad_val_)),
  b_nvmus_(NULL),
  c_nvmus_(false),
  run_(static_cast<int >(bad_val_)),
  b_run_(NULL),
  c_run_(false),
  els_ele105_(0),
  p_els_ele105_(&els_ele105_),
  b_els_ele105_(NULL),
  c_els_ele105_(false),
  els_ele23_(0),
  p_els_ele23_(&els_ele23_),
  b_els_ele23_(NULL),
  c_els_ele23_(false),
  els_ele8_(0),
  p_els_ele8_(&els_ele8_),
  b_els_ele8_(NULL),
  c_els_ele8_(false),
  els_inz_(0),
  p_els_inz_(&els_inz_),
  b_els_inz_(NULL),
  c_els_inz_(false),
  els_inzv_(0),
  p_els_inzv_(&els_inzv_),
  b_els_inzv_(NULL),
  c_els_inzv_(false),
  els_ispf_(0),
  p_els_ispf_(&els_ispf_),
  b_els_ispf_(NULL),
  c_els_ispf_(false),
  els_sig_(0),
  p_els_sig_(&els_sig_),
  b_els_sig_(NULL),
  c_els_sig_(false),
  els_sigid_(0),
  p_els_sigid_(&els_sigid_),
  b_els_sigid_(NULL),
  c_els_sigid_(false),
  els_tight_(0),
  p_els_tight_(&els_tight_),
  b_els_tight_(NULL),
  c_els_tight_(false),
  els_tm_(0),
  p_els_tm_(&els_tm_),
  b_els_tm_(NULL),
  c_els_tm_(false),
  els_vvvl_(0),
  p_els_vvvl_(&els_vvvl_),
  b_els_vvvl_(NULL),
  c_els_vvvl_(false),
  jets_islep_(0),
  p_jets_islep_(&jets_islep_),
  b_jets_islep_(NULL),
  c_jets_islep_(false),
  mus_inz_(0),
  p_mus_inz_(&mus_inz_),
  b_mus_inz_(NULL),
  c_mus_inz_(false),
  mus_inzv_(0),
  p_mus_inzv_(&mus_inzv_),
  b_mus_inzv_(NULL),
  c_mus_inzv_(false),
  mus_isomu18_(0),
  p_mus_isomu18_(&mus_isomu18_),
  b_mus_isomu18_(NULL),
  c_mus_isomu18_(false),
  mus_mu50_(0),
  p_mus_mu50_(&mus_mu50_),
  b_mus_mu50_(NULL),
  c_mus_mu50_(false),
  mus_mu8_(0),
  p_mus_mu8_(&mus_mu8_),
  b_mus_mu8_(NULL),
  c_mus_mu8_(false),
  mus_sig_(0),
  p_mus_sig_(&mus_sig_),
  b_mus_sig_(NULL),
  c_mus_sig_(false),
  mus_sigid_(0),
  p_mus_sigid_(&mus_sigid_),
  b_mus_sigid_(NULL),
  c_mus_sigid_(false),
  mus_tight_(0),
  p_mus_tight_(&mus_tight_),
  b_mus_tight_(NULL),
  c_mus_tight_(false),
  mus_tm_(0),
  p_mus_tm_(&mus_tm_),
  b_mus_tm_(NULL),
  c_mus_tm_(false),
  mus_vvvl_(0),
  p_mus_vvvl_(&mus_vvvl_),
  b_mus_vvvl_(NULL),
  c_mus_vvvl_(false),
  ph_ph90_(0),
  p_ph_ph90_(&ph_ph90_),
  b_ph_ph90_(NULL),
  c_ph_ph90_(false),
  ph_tm_(0),
  p_ph_tm_(&ph_tm_),
  b_ph_tm_(NULL),
  c_ph_tm_(false),
  sys_pass_(0),
  p_sys_pass_(&sys_pass_),
  b_sys_pass_(NULL),
  c_sys_pass_(false),
  trig_(0),
  p_trig_(&trig_),
  b_trig_(NULL),
  c_trig_(false),
  els_d0_(0),
  p_els_d0_(&els_d0_),
  b_els_d0_(NULL),
  c_els_d0_(false),
  els_dz_(0),
  p_els_dz_(&els_dz_),
  b_els_dz_(NULL),
  c_els_dz_(false),
  els_eta_(0),
  p_els_eta_(&els_eta_),
  b_els_eta_(NULL),
  c_els_eta_(false),
  els_miniso_(0),
  p_els_miniso_(&els_miniso_),
  b_els_miniso_(NULL),
  c_els_miniso_(false),
  els_phi_(0),
  p_els_phi_(&els_phi_),
  b_els_phi_(NULL),
  c_els_phi_(false),
  els_pt_(0),
  p_els_pt_(&els_pt_),
  b_els_pt_(NULL),
  c_els_pt_(false),
  els_reliso_(0),
  p_els_reliso_(&els_reliso_),
  b_els_reliso_(NULL),
  c_els_reliso_(false),
  els_sceta_(0),
  p_els_sceta_(&els_sceta_),
  b_els_sceta_(NULL),
  c_els_sceta_(false),
  els_vvvl_eta_(0),
  p_els_vvvl_eta_(&els_vvvl_eta_),
  b_els_vvvl_eta_(NULL),
  c_els_vvvl_eta_(false),
  els_vvvl_phi_(0),
  p_els_vvvl_phi_(&els_vvvl_phi_),
  b_els_vvvl_phi_(NULL),
  c_els_vvvl_phi_(false),
  els_vvvl_pt_(0),
  p_els_vvvl_pt_(&els_vvvl_pt_),
  b_els_vvvl_pt_(NULL),
  c_els_vvvl_pt_(false),
  fjets08_eta_(0),
  p_fjets08_eta_(&fjets08_eta_),
  b_fjets08_eta_(NULL),
  c_fjets08_eta_(false),
  fjets08_m_(0),
  p_fjets08_m_(&fjets08_m_),
  b_fjets08_m_(NULL),
  c_fjets08_m_(false),
  fjets08_phi_(0),
  p_fjets08_phi_(&fjets08_phi_),
  b_fjets08_phi_(NULL),
  c_fjets08_phi_(false),
  fjets08_poscsv_(0),
  p_fjets08_poscsv_(&fjets08_poscsv_),
  b_fjets08_poscsv_(NULL),
  c_fjets08_poscsv_(false),
  fjets08_pt_(0),
  p_fjets08_pt_(&fjets08_pt_),
  b_fjets08_pt_(NULL),
  c_fjets08_pt_(false),
  fjets08_sumcsv_(0),
  p_fjets08_sumcsv_(&fjets08_sumcsv_),
  b_fjets08_sumcsv_(NULL),
  c_fjets08_sumcsv_(false),
  fjets_eta_(0),
  p_fjets_eta_(&fjets_eta_),
  b_fjets_eta_(NULL),
  c_fjets_eta_(false),
  fjets_m_(0),
  p_fjets_m_(&fjets_m_),
  b_fjets_m_(NULL),
  c_fjets_m_(false),
  fjets_phi_(0),
  p_fjets_phi_(&fjets_phi_),
  b_fjets_phi_(NULL),
  c_fjets_phi_(false),
  fjets_poscsv_(0),
  p_fjets_poscsv_(&fjets_poscsv_),
  b_fjets_poscsv_(NULL),
  c_fjets_poscsv_(false),
  fjets_pt_(0),
  p_fjets_pt_(&fjets_pt_),
  b_fjets_pt_(NULL),
  c_fjets_pt_(false),
  fjets_sumcsv_(0),
  p_fjets_sumcsv_(&fjets_sumcsv_),
  b_fjets_sumcsv_(NULL),
  c_fjets_sumcsv_(false),
  jets_csv_(0),
  p_jets_csv_(&jets_csv_),
  b_jets_csv_(NULL),
  c_jets_csv_(false),
  jets_eta_(0),
  p_jets_eta_(&jets_eta_),
  b_jets_eta_(NULL),
  c_jets_eta_(false),
  jets_m_(0),
  p_jets_m_(&jets_m_),
  b_jets_m_(NULL),
  c_jets_m_(false),
  jets_phi_(0),
  p_jets_phi_(&jets_phi_),
  b_jets_phi_(NULL),
  c_jets_phi_(false),
  jets_pt_(0),
  p_jets_pt_(&jets_pt_),
  b_jets_pt_(NULL),
  c_jets_pt_(false),
  jets_pt_res_(0),
  p_jets_pt_res_(&jets_pt_res_),
  b_jets_pt_res_(NULL),
  c_jets_pt_res_(false),
  leps_eta_(0),
  p_leps_eta_(&leps_eta_),
  b_leps_eta_(NULL),
  c_leps_eta_(false),
  leps_id_(0),
  p_leps_id_(&leps_id_),
  b_leps_id_(NULL),
  c_leps_id_(false),
  leps_phi_(0),
  p_leps_phi_(&leps_phi_),
  b_leps_phi_(NULL),
  c_leps_phi_(false),
  leps_pt_(0),
  p_leps_pt_(&leps_pt_),
  b_leps_pt_(NULL),
  c_leps_pt_(false),
  mc_eta_(0),
  p_mc_eta_(&mc_eta_),
  b_mc_eta_(NULL),
  c_mc_eta_(false),
  mc_mass_(0),
  p_mc_mass_(&mc_mass_),
  b_mc_mass_(NULL),
  c_mc_mass_(false),
  mc_phi_(0),
  p_mc_phi_(&mc_phi_),
  b_mc_phi_(NULL),
  c_mc_phi_(false),
  mc_pt_(0),
  p_mc_pt_(&mc_pt_),
  b_mc_pt_(NULL),
  c_mc_pt_(false),
  mus_d0_(0),
  p_mus_d0_(&mus_d0_),
  b_mus_d0_(NULL),
  c_mus_d0_(false),
  mus_dz_(0),
  p_mus_dz_(&mus_dz_),
  b_mus_dz_(NULL),
  c_mus_dz_(false),
  mus_eta_(0),
  p_mus_eta_(&mus_eta_),
  b_mus_eta_(NULL),
  c_mus_eta_(false),
  mus_miniso_(0),
  p_mus_miniso_(&mus_miniso_),
  b_mus_miniso_(NULL),
  c_mus_miniso_(false),
  mus_phi_(0),
  p_mus_phi_(&mus_phi_),
  b_mus_phi_(NULL),
  c_mus_phi_(false),
  mus_pt_(0),
  p_mus_pt_(&mus_pt_),
  b_mus_pt_(NULL),
  c_mus_pt_(false),
  mus_reliso_(0),
  p_mus_reliso_(&mus_reliso_),
  b_mus_reliso_(NULL),
  c_mus_reliso_(false),
  mus_vvvl_eta_(0),
  p_mus_vvvl_eta_(&mus_vvvl_eta_),
  b_mus_vvvl_eta_(NULL),
  c_mus_vvvl_eta_(false),
  mus_vvvl_phi_(0),
  p_mus_vvvl_phi_(&mus_vvvl_phi_),
  b_mus_vvvl_phi_(NULL),
  c_mus_vvvl_phi_(false),
  mus_vvvl_pt_(0),
  p_mus_vvvl_pt_(&mus_vvvl_pt_),
  b_mus_vvvl_pt_(NULL),
  c_mus_vvvl_pt_(false),
  ph_eta_(0),
  p_ph_eta_(&ph_eta_),
  b_ph_eta_(NULL),
  c_ph_eta_(false),
  ph_phi_(0),
  p_ph_phi_(&ph_phi_),
  b_ph_phi_(NULL),
  c_ph_phi_(false),
  ph_pt_(0),
  p_ph_pt_(&ph_pt_),
  b_ph_pt_(NULL),
  c_ph_pt_(false),
  sys_bctag_(0),
  p_sys_bctag_(&sys_bctag_),
  b_sys_bctag_(NULL),
  c_sys_bctag_(false),
  sys_fs_bctag_(0),
  p_sys_fs_bctag_(&sys_fs_bctag_),
  b_sys_fs_bctag_(NULL),
  c_sys_fs_bctag_(false),
  sys_fs_lep_(0),
  p_sys_fs_lep_(&sys_fs_lep_),
  b_sys_fs_lep_(NULL),
  c_sys_fs_lep_(false),
  sys_fs_udsgtag_(0),
  p_sys_fs_udsgtag_(&sys_fs_udsgtag_),
  b_sys_fs_udsgtag_(NULL),
  c_sys_fs_udsgtag_(false),
  sys_ht_(0),
  p_sys_ht_(&sys_ht_),
  b_sys_ht_(NULL),
  c_sys_ht_(false),
  sys_isr_(0),
  p_sys_isr_(&sys_isr_),
  b_sys_isr_(NULL),
  c_sys_isr_(false),
  sys_lep_(0),
  p_sys_lep_(&sys_lep_),
  b_sys_lep_(NULL),
  c_sys_lep_(false),
  sys_met_(0),
  p_sys_met_(&sys_met_),
  b_sys_met_(NULL),
  c_sys_met_(false),
  sys_mj_(0),
  p_sys_mj_(&sys_mj_),
  b_sys_mj_(NULL),
  c_sys_mj_(false),
  sys_mt_(0),
  p_sys_mt_(&sys_mt_),
  b_sys_mt_(NULL),
  c_sys_mt_(false),
  sys_muf_(0),
  p_sys_muf_(&sys_muf_),
  b_sys_muf_(NULL),
  c_sys_muf_(false),
  sys_mur_(0),
  p_sys_mur_(&sys_mur_),
  b_sys_mur_(NULL),
  c_sys_mur_(false),
  sys_murf_(0),
  p_sys_murf_(&sys_murf_),
  b_sys_murf_(NULL),
  c_sys_murf_(false),
  sys_pdf_(0),
  p_sys_pdf_(&sys_pdf_),
  b_sys_pdf_(NULL),
  c_sys_pdf_(false),
  sys_trig_(0),
  p_sys_trig_(&sys_trig_),
  b_sys_trig_(NULL),
  c_sys_trig_(false),
  sys_udsgtag_(0),
  p_sys_udsgtag_(&sys_udsgtag_),
  b_sys_udsgtag_(NULL),
  c_sys_udsgtag_(false),
  tks_eta_(0),
  p_tks_eta_(&tks_eta_),
  b_tks_eta_(NULL),
  c_tks_eta_(false),
  tks_miniso_(0),
  p_tks_miniso_(&tks_miniso_),
  b_tks_miniso_(NULL),
  c_tks_miniso_(false),
  tks_mt_(0),
  p_tks_mt_(&tks_mt_),
  b_tks_mt_(NULL),
  c_tks_mt_(false),
  tks_mt2_(0),
  p_tks_mt2_(&tks_mt2_),
  b_tks_mt2_(NULL),
  c_tks_mt2_(false),
  tks_phi_(0),
  p_tks_phi_(&tks_phi_),
  b_tks_phi_(NULL),
  c_tks_phi_(false),
  tks_pt_(0),
  p_tks_pt_(&tks_pt_),
  b_tks_pt_(NULL),
  c_tks_pt_(false),
  trig_prescale_(0),
  p_trig_prescale_(&trig_prescale_),
  b_trig_prescale_(NULL),
  c_trig_prescale_(false),
  w_pdf_(0),
  p_w_pdf_(&w_pdf_),
  b_w_pdf_(NULL),
  c_w_pdf_(false),
  els_charge_(0),
  p_els_charge_(&els_charge_),
  b_els_charge_(NULL),
  c_els_charge_(false),
  fjets08_btags_(0),
  p_fjets08_btags_(&fjets08_btags_),
  b_fjets08_btags_(NULL),
  c_fjets08_btags_(false),
  fjets08_nconst_(0),
  p_fjets08_nconst_(&fjets08_nconst_),
  b_fjets08_nconst_(NULL),
  c_fjets08_nconst_(false),
  fjets_btags_(0),
  p_fjets_btags_(&fjets_btags_),
  b_fjets_btags_(NULL),
  c_fjets_btags_(false),
  fjets_nconst_(0),
  p_fjets_nconst_(&fjets_nconst_),
  b_fjets_nconst_(NULL),
  c_fjets_nconst_(false),
  jets_fjet08_index_(0),
  p_jets_fjet08_index_(&jets_fjet08_index_),
  b_jets_fjet08_index_(NULL),
  c_jets_fjet08_index_(false),
  jets_fjet_index_(0),
  p_jets_fjet_index_(&jets_fjet_index_),
  b_jets_fjet_index_(NULL),
  c_jets_fjet_index_(false),
  jets_hflavor_(0),
  p_jets_hflavor_(&jets_hflavor_),
  b_jets_hflavor_(NULL),
  c_jets_hflavor_(false),
  mc_id_(0),
  p_mc_id_(&mc_id_),
  b_mc_id_(NULL),
  c_mc_id_(false),
  mc_mom_(0),
  p_mc_mom_(&mc_mom_),
  b_mc_mom_(NULL),
  c_mc_mom_(false),
  mus_charge_(0),
  p_mus_charge_(&mus_charge_),
  b_mus_charge_(NULL),
  c_mus_charge_(false),
  sys_nbm_(0),
  p_sys_nbm_(&sys_nbm_),
  b_sys_nbm_(NULL),
  c_sys_nbm_(false),
  sys_njets_(0),
  p_sys_njets_(&sys_njets_),
  b_sys_njets_(NULL),
  c_sys_njets_(false),
  tks_pdg_(0),
  p_tks_pdg_(&tks_pdg_),
  b_tks_pdg_(NULL),
  c_tks_pdg_(false){
  chain_.Add(filename.c_str());
  chain_.SetBranchAddress("event", &event_, &b_event_);
  chain_.SetBranchAddress("fromGS", &fromGS_, &b_fromGS_);
  chain_.SetBranchAddress("jetmismeas", &jetmismeas_, &b_jetmismeas_);
  chain_.SetBranchAddress("low_dphi", &low_dphi_, &b_low_dphi_);
  chain_.SetBranchAddress("nonblind", &nonblind_, &b_nonblind_);
  chain_.SetBranchAddress("pass", &pass_, &b_pass_);
  chain_.SetBranchAddress("pass_cschalo", &pass_cschalo_, &b_pass_cschalo_);
  chain_.SetBranchAddress("pass_eebadsc", &pass_eebadsc_, &b_pass_eebadsc_);
  chain_.SetBranchAddress("pass_goodv", &pass_goodv_, &b_pass_goodv_);
  chain_.SetBranchAddress("pass_hbhe", &pass_hbhe_, &b_pass_hbhe_);
  chain_.SetBranchAddress("pass_hbheiso", &pass_hbheiso_, &b_pass_hbheiso_);
  chain_.SetBranchAddress("pass_jets", &pass_jets_, &b_pass_jets_);
  chain_.SetBranchAddress("pass_jets_nohf", &pass_jets_nohf_, &b_pass_jets_nohf_);
  chain_.SetBranchAddress("pass_jets_ra2", &pass_jets_ra2_, &b_pass_jets_ra2_);
  chain_.SetBranchAddress("pass_jets_tight", &pass_jets_tight_, &b_pass_jets_tight_);
  chain_.SetBranchAddress("pass_jets_tight_ra2", &pass_jets_tight_ra2_, &b_pass_jets_tight_ra2_);
  chain_.SetBranchAddress("pass_nohf", &pass_nohf_, &b_pass_nohf_);
  chain_.SetBranchAddress("pass_ra2", &pass_ra2_, &b_pass_ra2_);
  chain_.SetBranchAddress("stitch", &stitch_, &b_stitch_);
  chain_.SetBranchAddress("dphi1", &dphi1_, &b_dphi1_);
  chain_.SetBranchAddress("dphi2", &dphi2_, &b_dphi2_);
  chain_.SetBranchAddress("dphi3", &dphi3_, &b_dphi3_);
  chain_.SetBranchAddress("dphi4", &dphi4_, &b_dphi4_);
  chain_.SetBranchAddress("dphi_wlep", &dphi_wlep_, &b_dphi_wlep_);
  chain_.SetBranchAddress("eff_trig", &eff_trig_, &b_eff_trig_);
  chain_.SetBranchAddress("elel_eta", &elel_eta_, &b_elel_eta_);
  chain_.SetBranchAddress("elel_m", &elel_m_, &b_elel_m_);
  chain_.SetBranchAddress("elel_phi", &elel_phi_, &b_elel_phi_);
  chain_.SetBranchAddress("elel_pt", &elel_pt_, &b_elel_pt_);
  chain_.SetBranchAddress("elel_pt1", &elel_pt1_, &b_elel_pt1_);
  chain_.SetBranchAddress("elel_pt2", &elel_pt2_, &b_elel_pt2_);
  chain_.SetBranchAddress("elel_w", &elel_w_, &b_elel_w_);
  chain_.SetBranchAddress("elelv_eta", &elelv_eta_, &b_elelv_eta_);
  chain_.SetBranchAddress("elelv_m", &elelv_m_, &b_elelv_m_);
  chain_.SetBranchAddress("elelv_phi", &elelv_phi_, &b_elelv_phi_);
  chain_.SetBranchAddress("elelv_pt", &elelv_pt_, &b_elelv_pt_);
  chain_.SetBranchAddress("elelv_pt1", &elelv_pt1_, &b_elelv_pt1_);
  chain_.SetBranchAddress("elelv_pt2", &elelv_pt2_, &b_elelv_pt2_);
  chain_.SetBranchAddress("elelv_w", &elelv_w_, &b_elelv_w_);
  chain_.SetBranchAddress("elmu_eta", &elmu_eta_, &b_elmu_eta_);
  chain_.SetBranchAddress("elmu_m", &elmu_m_, &b_elmu_m_);
  chain_.SetBranchAddress("elmu_phi", &elmu_phi_, &b_elmu_phi_);
  chain_.SetBranchAddress("elmu_pt", &elmu_pt_, &b_elmu_pt_);
  chain_.SetBranchAddress("elmu_pt1", &elmu_pt1_, &b_elmu_pt1_);
  chain_.SetBranchAddress("elmu_pt2", &elmu_pt2_, &b_elmu_pt2_);
  chain_.SetBranchAddress("elmu_w", &elmu_w_, &b_elmu_w_);
  chain_.SetBranchAddress("ht", &ht_, &b_ht_);
  chain_.SetBranchAddress("ht_clean", &ht_clean_, &b_ht_clean_);
  chain_.SetBranchAddress("ht_hlt", &ht_hlt_, &b_ht_hlt_);
  chain_.SetBranchAddress("ht_isr_me", &ht_isr_me_, &b_ht_isr_me_);
  chain_.SetBranchAddress("ht_ra2", &ht_ra2_, &b_ht_ra2_);
  chain_.SetBranchAddress("ht_tru", &ht_tru_, &b_ht_tru_);
  chain_.SetBranchAddress("isr_tru_eta", &isr_tru_eta_, &b_isr_tru_eta_);
  chain_.SetBranchAddress("isr_tru_phi", &isr_tru_phi_, &b_isr_tru_phi_);
  chain_.SetBranchAddress("isr_tru_pt", &isr_tru_pt_, &b_isr_tru_pt_);
  chain_.SetBranchAddress("jetsys_eta", &jetsys_eta_, &b_jetsys_eta_);
  chain_.SetBranchAddress("jetsys_nob_eta", &jetsys_nob_eta_, &b_jetsys_nob_eta_);
  chain_.SetBranchAddress("jetsys_nob_phi", &jetsys_nob_phi_, &b_jetsys_nob_phi_);
  chain_.SetBranchAddress("jetsys_nob_pt", &jetsys_nob_pt_, &b_jetsys_nob_pt_);
  chain_.SetBranchAddress("jetsys_phi", &jetsys_phi_, &b_jetsys_phi_);
  chain_.SetBranchAddress("jetsys_pt", &jetsys_pt_, &b_jetsys_pt_);
  chain_.SetBranchAddress("met", &met_, &b_met_);
  chain_.SetBranchAddress("met_calo", &met_calo_, &b_met_calo_);
  chain_.SetBranchAddress("met_calo_phi", &met_calo_phi_, &b_met_calo_phi_);
  chain_.SetBranchAddress("met_mini", &met_mini_, &b_met_mini_);
  chain_.SetBranchAddress("met_mini_phi", &met_mini_phi_, &b_met_mini_phi_);
  chain_.SetBranchAddress("met_nohf", &met_nohf_, &b_met_nohf_);
  chain_.SetBranchAddress("met_nohf_phi", &met_nohf_phi_, &b_met_nohf_phi_);
  chain_.SetBranchAddress("met_phi", &met_phi_, &b_met_phi_);
  chain_.SetBranchAddress("met_raw", &met_raw_, &b_met_raw_);
  chain_.SetBranchAddress("met_raw_phi", &met_raw_phi_, &b_met_raw_phi_);
  chain_.SetBranchAddress("met_rebal", &met_rebal_, &b_met_rebal_);
  chain_.SetBranchAddress("met_tru", &met_tru_, &b_met_tru_);
  chain_.SetBranchAddress("met_tru_nuw", &met_tru_nuw_, &b_met_tru_nuw_);
  chain_.SetBranchAddress("met_tru_nuw_phi", &met_tru_nuw_phi_, &b_met_tru_nuw_phi_);
  chain_.SetBranchAddress("met_tru_phi", &met_tru_phi_, &b_met_tru_phi_);
  chain_.SetBranchAddress("mht", &mht_, &b_mht_);
  chain_.SetBranchAddress("mht_clean", &mht_clean_, &b_mht_clean_);
  chain_.SetBranchAddress("mht_clean_phi", &mht_clean_phi_, &b_mht_clean_phi_);
  chain_.SetBranchAddress("mht_phi", &mht_phi_, &b_mht_phi_);
  chain_.SetBranchAddress("mj", &mj_, &b_mj_);
  chain_.SetBranchAddress("mj08", &mj08_, &b_mj08_);
  chain_.SetBranchAddress("mt", &mt_, &b_mt_);
  chain_.SetBranchAddress("mt_nohf", &mt_nohf_, &b_mt_nohf_);
  chain_.SetBranchAddress("mt_rebal", &mt_rebal_, &b_mt_rebal_);
  chain_.SetBranchAddress("mt_tru", &mt_tru_, &b_mt_tru_);
  chain_.SetBranchAddress("mt_tru_nuw", &mt_tru_nuw_, &b_mt_tru_nuw_);
  chain_.SetBranchAddress("mumu_eta", &mumu_eta_, &b_mumu_eta_);
  chain_.SetBranchAddress("mumu_m", &mumu_m_, &b_mumu_m_);
  chain_.SetBranchAddress("mumu_phi", &mumu_phi_, &b_mumu_phi_);
  chain_.SetBranchAddress("mumu_pt", &mumu_pt_, &b_mumu_pt_);
  chain_.SetBranchAddress("mumu_pt1", &mumu_pt1_, &b_mumu_pt1_);
  chain_.SetBranchAddress("mumu_pt2", &mumu_pt2_, &b_mumu_pt2_);
  chain_.SetBranchAddress("mumu_w", &mumu_w_, &b_mumu_w_);
  chain_.SetBranchAddress("mumuv_eta", &mumuv_eta_, &b_mumuv_eta_);
  chain_.SetBranchAddress("mumuv_m", &mumuv_m_, &b_mumuv_m_);
  chain_.SetBranchAddress("mumuv_phi", &mumuv_phi_, &b_mumuv_phi_);
  chain_.SetBranchAddress("mumuv_pt", &mumuv_pt_, &b_mumuv_pt_);
  chain_.SetBranchAddress("mumuv_pt1", &mumuv_pt1_, &b_mumuv_pt1_);
  chain_.SetBranchAddress("mumuv_pt2", &mumuv_pt2_, &b_mumuv_pt2_);
  chain_.SetBranchAddress("mumuv_w", &mumuv_w_, &b_mumuv_w_);
  chain_.SetBranchAddress("ntrupv_mean", &ntrupv_mean_, &b_ntrupv_mean_);
  chain_.SetBranchAddress("onel_ele105", &onel_ele105_, &b_onel_ele105_);
  chain_.SetBranchAddress("onel_ele23", &onel_ele23_, &b_onel_ele23_);
  chain_.SetBranchAddress("onel_ele8", &onel_ele8_, &b_onel_ele8_);
  chain_.SetBranchAddress("onel_vvvl", &onel_vvvl_, &b_onel_vvvl_);
  chain_.SetBranchAddress("onht", &onht_, &b_onht_);
  chain_.SetBranchAddress("onmet", &onmet_, &b_onmet_);
  chain_.SetBranchAddress("onmu_isomu18", &onmu_isomu18_, &b_onmu_isomu18_);
  chain_.SetBranchAddress("onmu_mu50", &onmu_mu50_, &b_onmu_mu50_);
  chain_.SetBranchAddress("onmu_mu8", &onmu_mu8_, &b_onmu_mu8_);
  chain_.SetBranchAddress("onmu_vvvl", &onmu_vvvl_, &b_onmu_vvvl_);
  chain_.SetBranchAddress("onph_ph90", &onph_ph90_, &b_onph_ph90_);
  chain_.SetBranchAddress("w_btag", &w_btag_, &b_w_btag_);
  chain_.SetBranchAddress("w_fs_lep", &w_fs_lep_, &b_w_fs_lep_);
  chain_.SetBranchAddress("w_lep", &w_lep_, &b_w_lep_);
  chain_.SetBranchAddress("w_lumi", &w_lumi_, &b_w_lumi_);
  chain_.SetBranchAddress("w_pu", &w_pu_, &b_w_pu_);
  chain_.SetBranchAddress("w_toppt", &w_toppt_, &b_w_toppt_);
  chain_.SetBranchAddress("weight", &weight_, &b_weight_);
  chain_.SetBranchAddress("lumiblock", &lumiblock_, &b_lumiblock_);
  chain_.SetBranchAddress("mgluino", &mgluino_, &b_mgluino_);
  chain_.SetBranchAddress("mlsp", &mlsp_, &b_mlsp_);
  chain_.SetBranchAddress("nbl", &nbl_, &b_nbl_);
  chain_.SetBranchAddress("nbm", &nbm_, &b_nbm_);
  chain_.SetBranchAddress("nbm_ra2", &nbm_ra2_, &b_nbm_ra2_);
  chain_.SetBranchAddress("nbt", &nbt_, &b_nbt_);
  chain_.SetBranchAddress("nels", &nels_, &b_nels_);
  chain_.SetBranchAddress("nels_ele23", &nels_ele23_, &b_nels_ele23_);
  chain_.SetBranchAddress("nels_vvvl", &nels_vvvl_, &b_nels_vvvl_);
  chain_.SetBranchAddress("nfjets", &nfjets_, &b_nfjets_);
  chain_.SetBranchAddress("nfjets08", &nfjets08_, &b_nfjets08_);
  chain_.SetBranchAddress("nisr_me", &nisr_me_, &b_nisr_me_);
  chain_.SetBranchAddress("njets", &njets_, &b_njets_);
  chain_.SetBranchAddress("njets_clean", &njets_clean_, &b_njets_clean_);
  chain_.SetBranchAddress("njets_ra2", &njets_ra2_, &b_njets_ra2_);
  chain_.SetBranchAddress("nleps", &nleps_, &b_nleps_);
  chain_.SetBranchAddress("nleps_tm", &nleps_tm_, &b_nleps_tm_);
  chain_.SetBranchAddress("nmus", &nmus_, &b_nmus_);
  chain_.SetBranchAddress("nmus_isomu18", &nmus_isomu18_, &b_nmus_isomu18_);
  chain_.SetBranchAddress("nmus_vvvl", &nmus_vvvl_, &b_nmus_vvvl_);
  chain_.SetBranchAddress("nph", &nph_, &b_nph_);
  chain_.SetBranchAddress("npv", &npv_, &b_npv_);
  chain_.SetBranchAddress("ntks", &ntks_, &b_ntks_);
  chain_.SetBranchAddress("ntruels", &ntruels_, &b_ntruels_);
  chain_.SetBranchAddress("ntruleps", &ntruleps_, &b_ntruleps_);
  chain_.SetBranchAddress("ntrumus", &ntrumus_, &b_ntrumus_);
  chain_.SetBranchAddress("ntrupv", &ntrupv_, &b_ntrupv_);
  chain_.SetBranchAddress("ntrutaush", &ntrutaush_, &b_ntrutaush_);
  chain_.SetBranchAddress("ntrutausl", &ntrutausl_, &b_ntrutausl_);
  chain_.SetBranchAddress("nvels", &nvels_, &b_nvels_);
  chain_.SetBranchAddress("nveto", &nveto_, &b_nveto_);
  chain_.SetBranchAddress("nvleps", &nvleps_, &b_nvleps_);
  chain_.SetBranchAddress("nvmus", &nvmus_, &b_nvmus_);
  chain_.SetBranchAddress("run", &run_, &b_run_);
  chain_.SetBranchAddress("els_ele105", &p_els_ele105_, &b_els_ele105_);
  chain_.SetBranchAddress("els_ele23", &p_els_ele23_, &b_els_ele23_);
  chain_.SetBranchAddress("els_ele8", &p_els_ele8_, &b_els_ele8_);
  chain_.SetBranchAddress("els_inz", &p_els_inz_, &b_els_inz_);
  chain_.SetBranchAddress("els_inzv", &p_els_inzv_, &b_els_inzv_);
  chain_.SetBranchAddress("els_ispf", &p_els_ispf_, &b_els_ispf_);
  chain_.SetBranchAddress("els_sig", &p_els_sig_, &b_els_sig_);
  chain_.SetBranchAddress("els_sigid", &p_els_sigid_, &b_els_sigid_);
  chain_.SetBranchAddress("els_tight", &p_els_tight_, &b_els_tight_);
  chain_.SetBranchAddress("els_tm", &p_els_tm_, &b_els_tm_);
  chain_.SetBranchAddress("els_vvvl", &p_els_vvvl_, &b_els_vvvl_);
  chain_.SetBranchAddress("jets_islep", &p_jets_islep_, &b_jets_islep_);
  chain_.SetBranchAddress("mus_inz", &p_mus_inz_, &b_mus_inz_);
  chain_.SetBranchAddress("mus_inzv", &p_mus_inzv_, &b_mus_inzv_);
  chain_.SetBranchAddress("mus_isomu18", &p_mus_isomu18_, &b_mus_isomu18_);
  chain_.SetBranchAddress("mus_mu50", &p_mus_mu50_, &b_mus_mu50_);
  chain_.SetBranchAddress("mus_mu8", &p_mus_mu8_, &b_mus_mu8_);
  chain_.SetBranchAddress("mus_sig", &p_mus_sig_, &b_mus_sig_);
  chain_.SetBranchAddress("mus_sigid", &p_mus_sigid_, &b_mus_sigid_);
  chain_.SetBranchAddress("mus_tight", &p_mus_tight_, &b_mus_tight_);
  chain_.SetBranchAddress("mus_tm", &p_mus_tm_, &b_mus_tm_);
  chain_.SetBranchAddress("mus_vvvl", &p_mus_vvvl_, &b_mus_vvvl_);
  chain_.SetBranchAddress("ph_ph90", &p_ph_ph90_, &b_ph_ph90_);
  chain_.SetBranchAddress("ph_tm", &p_ph_tm_, &b_ph_tm_);
  chain_.SetBranchAddress("sys_pass", &p_sys_pass_, &b_sys_pass_);
  chain_.SetBranchAddress("trig", &p_trig_, &b_trig_);
  chain_.SetBranchAddress("els_d0", &p_els_d0_, &b_els_d0_);
  chain_.SetBranchAddress("els_dz", &p_els_dz_, &b_els_dz_);
  chain_.SetBranchAddress("els_eta", &p_els_eta_, &b_els_eta_);
  chain_.SetBranchAddress("els_miniso", &p_els_miniso_, &b_els_miniso_);
  chain_.SetBranchAddress("els_phi", &p_els_phi_, &b_els_phi_);
  chain_.SetBranchAddress("els_pt", &p_els_pt_, &b_els_pt_);
  chain_.SetBranchAddress("els_reliso", &p_els_reliso_, &b_els_reliso_);
  chain_.SetBranchAddress("els_sceta", &p_els_sceta_, &b_els_sceta_);
  chain_.SetBranchAddress("els_vvvl_eta", &p_els_vvvl_eta_, &b_els_vvvl_eta_);
  chain_.SetBranchAddress("els_vvvl_phi", &p_els_vvvl_phi_, &b_els_vvvl_phi_);
  chain_.SetBranchAddress("els_vvvl_pt", &p_els_vvvl_pt_, &b_els_vvvl_pt_);
  chain_.SetBranchAddress("fjets08_eta", &p_fjets08_eta_, &b_fjets08_eta_);
  chain_.SetBranchAddress("fjets08_m", &p_fjets08_m_, &b_fjets08_m_);
  chain_.SetBranchAddress("fjets08_phi", &p_fjets08_phi_, &b_fjets08_phi_);
  chain_.SetBranchAddress("fjets08_poscsv", &p_fjets08_poscsv_, &b_fjets08_poscsv_);
  chain_.SetBranchAddress("fjets08_pt", &p_fjets08_pt_, &b_fjets08_pt_);
  chain_.SetBranchAddress("fjets08_sumcsv", &p_fjets08_sumcsv_, &b_fjets08_sumcsv_);
  chain_.SetBranchAddress("fjets_eta", &p_fjets_eta_, &b_fjets_eta_);
  chain_.SetBranchAddress("fjets_m", &p_fjets_m_, &b_fjets_m_);
  chain_.SetBranchAddress("fjets_phi", &p_fjets_phi_, &b_fjets_phi_);
  chain_.SetBranchAddress("fjets_poscsv", &p_fjets_poscsv_, &b_fjets_poscsv_);
  chain_.SetBranchAddress("fjets_pt", &p_fjets_pt_, &b_fjets_pt_);
  chain_.SetBranchAddress("fjets_sumcsv", &p_fjets_sumcsv_, &b_fjets_sumcsv_);
  chain_.SetBranchAddress("jets_csv", &p_jets_csv_, &b_jets_csv_);
  chain_.SetBranchAddress("jets_eta", &p_jets_eta_, &b_jets_eta_);
  chain_.SetBranchAddress("jets_m", &p_jets_m_, &b_jets_m_);
  chain_.SetBranchAddress("jets_phi", &p_jets_phi_, &b_jets_phi_);
  chain_.SetBranchAddress("jets_pt", &p_jets_pt_, &b_jets_pt_);
  chain_.SetBranchAddress("jets_pt_res", &p_jets_pt_res_, &b_jets_pt_res_);
  chain_.SetBranchAddress("leps_eta", &p_leps_eta_, &b_leps_eta_);
  chain_.SetBranchAddress("leps_id", &p_leps_id_, &b_leps_id_);
  chain_.SetBranchAddress("leps_phi", &p_leps_phi_, &b_leps_phi_);
  chain_.SetBranchAddress("leps_pt", &p_leps_pt_, &b_leps_pt_);
  chain_.SetBranchAddress("mc_eta", &p_mc_eta_, &b_mc_eta_);
  chain_.SetBranchAddress("mc_mass", &p_mc_mass_, &b_mc_mass_);
  chain_.SetBranchAddress("mc_phi", &p_mc_phi_, &b_mc_phi_);
  chain_.SetBranchAddress("mc_pt", &p_mc_pt_, &b_mc_pt_);
  chain_.SetBranchAddress("mus_d0", &p_mus_d0_, &b_mus_d0_);
  chain_.SetBranchAddress("mus_dz", &p_mus_dz_, &b_mus_dz_);
  chain_.SetBranchAddress("mus_eta", &p_mus_eta_, &b_mus_eta_);
  chain_.SetBranchAddress("mus_miniso", &p_mus_miniso_, &b_mus_miniso_);
  chain_.SetBranchAddress("mus_phi", &p_mus_phi_, &b_mus_phi_);
  chain_.SetBranchAddress("mus_pt", &p_mus_pt_, &b_mus_pt_);
  chain_.SetBranchAddress("mus_reliso", &p_mus_reliso_, &b_mus_reliso_);
  chain_.SetBranchAddress("mus_vvvl_eta", &p_mus_vvvl_eta_, &b_mus_vvvl_eta_);
  chain_.SetBranchAddress("mus_vvvl_phi", &p_mus_vvvl_phi_, &b_mus_vvvl_phi_);
  chain_.SetBranchAddress("mus_vvvl_pt", &p_mus_vvvl_pt_, &b_mus_vvvl_pt_);
  chain_.SetBranchAddress("ph_eta", &p_ph_eta_, &b_ph_eta_);
  chain_.SetBranchAddress("ph_phi", &p_ph_phi_, &b_ph_phi_);
  chain_.SetBranchAddress("ph_pt", &p_ph_pt_, &b_ph_pt_);
  chain_.SetBranchAddress("sys_bctag", &p_sys_bctag_, &b_sys_bctag_);
  chain_.SetBranchAddress("sys_fs_bctag", &p_sys_fs_bctag_, &b_sys_fs_bctag_);
  chain_.SetBranchAddress("sys_fs_lep", &p_sys_fs_lep_, &b_sys_fs_lep_);
  chain_.SetBranchAddress("sys_fs_udsgtag", &p_sys_fs_udsgtag_, &b_sys_fs_udsgtag_);
  chain_.SetBranchAddress("sys_ht", &p_sys_ht_, &b_sys_ht_);
  chain_.SetBranchAddress("sys_isr", &p_sys_isr_, &b_sys_isr_);
  chain_.SetBranchAddress("sys_lep", &p_sys_lep_, &b_sys_lep_);
  chain_.SetBranchAddress("sys_met", &p_sys_met_, &b_sys_met_);
  chain_.SetBranchAddress("sys_mj", &p_sys_mj_, &b_sys_mj_);
  chain_.SetBranchAddress("sys_mt", &p_sys_mt_, &b_sys_mt_);
  chain_.SetBranchAddress("sys_muf", &p_sys_muf_, &b_sys_muf_);
  chain_.SetBranchAddress("sys_mur", &p_sys_mur_, &b_sys_mur_);
  chain_.SetBranchAddress("sys_murf", &p_sys_murf_, &b_sys_murf_);
  chain_.SetBranchAddress("sys_pdf", &p_sys_pdf_, &b_sys_pdf_);
  chain_.SetBranchAddress("sys_trig", &p_sys_trig_, &b_sys_trig_);
  chain_.SetBranchAddress("sys_udsgtag", &p_sys_udsgtag_, &b_sys_udsgtag_);
  chain_.SetBranchAddress("tks_eta", &p_tks_eta_, &b_tks_eta_);
  chain_.SetBranchAddress("tks_miniso", &p_tks_miniso_, &b_tks_miniso_);
  chain_.SetBranchAddress("tks_mt", &p_tks_mt_, &b_tks_mt_);
  chain_.SetBranchAddress("tks_mt2", &p_tks_mt2_, &b_tks_mt2_);
  chain_.SetBranchAddress("tks_phi", &p_tks_phi_, &b_tks_phi_);
  chain_.SetBranchAddress("tks_pt", &p_tks_pt_, &b_tks_pt_);
  chain_.SetBranchAddress("trig_prescale", &p_trig_prescale_, &b_trig_prescale_);
  chain_.SetBranchAddress("w_pdf", &p_w_pdf_, &b_w_pdf_);
  chain_.SetBranchAddress("els_charge", &p_els_charge_, &b_els_charge_);
  chain_.SetBranchAddress("fjets08_btags", &p_fjets08_btags_, &b_fjets08_btags_);
  chain_.SetBranchAddress("fjets08_nconst", &p_fjets08_nconst_, &b_fjets08_nconst_);
  chain_.SetBranchAddress("fjets_btags", &p_fjets_btags_, &b_fjets_btags_);
  chain_.SetBranchAddress("fjets_nconst", &p_fjets_nconst_, &b_fjets_nconst_);
  chain_.SetBranchAddress("jets_fjet08_index", &p_jets_fjet08_index_, &b_jets_fjet08_index_);
  chain_.SetBranchAddress("jets_fjet_index", &p_jets_fjet_index_, &b_jets_fjet_index_);
  chain_.SetBranchAddress("jets_hflavor", &p_jets_hflavor_, &b_jets_hflavor_);
  chain_.SetBranchAddress("mc_id", &p_mc_id_, &b_mc_id_);
  chain_.SetBranchAddress("mc_mom", &p_mc_mom_, &b_mc_mom_);
  chain_.SetBranchAddress("mus_charge", &p_mus_charge_, &b_mus_charge_);
  chain_.SetBranchAddress("sys_nbm", &p_sys_nbm_, &b_sys_nbm_);
  chain_.SetBranchAddress("sys_njets", &p_sys_njets_, &b_sys_njets_);
  chain_.SetBranchAddress("tks_pdg", &p_tks_pdg_, &b_tks_pdg_);
}

void baby_base::Fill(){
  if(read_only_){
    throw std::logic_error("Trying to write to read-only tree");
  }else{
    tree_.Fill();
  }

}

void baby_base::Clear(){
  if(read_only_){
    throw std::logic_error("Trying to write to read-only tree");
  }else{
  //Resetting variables
  event_ = static_cast<Long64_t >(bad_val_);
  fromGS_ = static_cast<bool >(bad_val_);
  jetmismeas_ = static_cast<bool >(bad_val_);
  low_dphi_ = static_cast<bool >(bad_val_);
  nonblind_ = static_cast<bool >(bad_val_);
  pass_ = static_cast<bool >(bad_val_);
  pass_cschalo_ = static_cast<bool >(bad_val_);
  pass_eebadsc_ = static_cast<bool >(bad_val_);
  pass_goodv_ = static_cast<bool >(bad_val_);
  pass_hbhe_ = static_cast<bool >(bad_val_);
  pass_hbheiso_ = static_cast<bool >(bad_val_);
  pass_jets_ = static_cast<bool >(bad_val_);
  pass_jets_nohf_ = static_cast<bool >(bad_val_);
  pass_jets_ra2_ = static_cast<bool >(bad_val_);
  pass_jets_tight_ = static_cast<bool >(bad_val_);
  pass_jets_tight_ra2_ = static_cast<bool >(bad_val_);
  pass_nohf_ = static_cast<bool >(bad_val_);
  pass_ra2_ = static_cast<bool >(bad_val_);
  stitch_ = static_cast<bool >(bad_val_);
  dphi1_ = static_cast<float >(bad_val_);
  dphi2_ = static_cast<float >(bad_val_);
  dphi3_ = static_cast<float >(bad_val_);
  dphi4_ = static_cast<float >(bad_val_);
  dphi_wlep_ = static_cast<float >(bad_val_);
  eff_trig_ = static_cast<float >(bad_val_);
  elel_eta_ = static_cast<float >(bad_val_);
  elel_m_ = static_cast<float >(bad_val_);
  elel_phi_ = static_cast<float >(bad_val_);
  elel_pt_ = static_cast<float >(bad_val_);
  elel_pt1_ = static_cast<float >(bad_val_);
  elel_pt2_ = static_cast<float >(bad_val_);
  elel_w_ = static_cast<float >(bad_val_);
  elelv_eta_ = static_cast<float >(bad_val_);
  elelv_m_ = static_cast<float >(bad_val_);
  elelv_phi_ = static_cast<float >(bad_val_);
  elelv_pt_ = static_cast<float >(bad_val_);
  elelv_pt1_ = static_cast<float >(bad_val_);
  elelv_pt2_ = static_cast<float >(bad_val_);
  elelv_w_ = static_cast<float >(bad_val_);
  elmu_eta_ = static_cast<float >(bad_val_);
  elmu_m_ = static_cast<float >(bad_val_);
  elmu_phi_ = static_cast<float >(bad_val_);
  elmu_pt_ = static_cast<float >(bad_val_);
  elmu_pt1_ = static_cast<float >(bad_val_);
  elmu_pt2_ = static_cast<float >(bad_val_);
  elmu_w_ = static_cast<float >(bad_val_);
  ht_ = static_cast<float >(bad_val_);
  ht_clean_ = static_cast<float >(bad_val_);
  ht_hlt_ = static_cast<float >(bad_val_);
  ht_isr_me_ = static_cast<float >(bad_val_);
  ht_ra2_ = static_cast<float >(bad_val_);
  ht_tru_ = static_cast<float >(bad_val_);
  isr_tru_eta_ = static_cast<float >(bad_val_);
  isr_tru_phi_ = static_cast<float >(bad_val_);
  isr_tru_pt_ = static_cast<float >(bad_val_);
  jetsys_eta_ = static_cast<float >(bad_val_);
  jetsys_nob_eta_ = static_cast<float >(bad_val_);
  jetsys_nob_phi_ = static_cast<float >(bad_val_);
  jetsys_nob_pt_ = static_cast<float >(bad_val_);
  jetsys_phi_ = static_cast<float >(bad_val_);
  jetsys_pt_ = static_cast<float >(bad_val_);
  met_ = static_cast<float >(bad_val_);
  met_calo_ = static_cast<float >(bad_val_);
  met_calo_phi_ = static_cast<float >(bad_val_);
  met_mini_ = static_cast<float >(bad_val_);
  met_mini_phi_ = static_cast<float >(bad_val_);
  met_nohf_ = static_cast<float >(bad_val_);
  met_nohf_phi_ = static_cast<float >(bad_val_);
  met_phi_ = static_cast<float >(bad_val_);
  met_raw_ = static_cast<float >(bad_val_);
  met_raw_phi_ = static_cast<float >(bad_val_);
  met_rebal_ = static_cast<float >(bad_val_);
  met_tru_ = static_cast<float >(bad_val_);
  met_tru_nuw_ = static_cast<float >(bad_val_);
  met_tru_nuw_phi_ = static_cast<float >(bad_val_);
  met_tru_phi_ = static_cast<float >(bad_val_);
  mht_ = static_cast<float >(bad_val_);
  mht_clean_ = static_cast<float >(bad_val_);
  mht_clean_phi_ = static_cast<float >(bad_val_);
  mht_phi_ = static_cast<float >(bad_val_);
  mj_ = static_cast<float >(bad_val_);
  mj08_ = static_cast<float >(bad_val_);
  mt_ = static_cast<float >(bad_val_);
  mt_nohf_ = static_cast<float >(bad_val_);
  mt_rebal_ = static_cast<float >(bad_val_);
  mt_tru_ = static_cast<float >(bad_val_);
  mt_tru_nuw_ = static_cast<float >(bad_val_);
  mumu_eta_ = static_cast<float >(bad_val_);
  mumu_m_ = static_cast<float >(bad_val_);
  mumu_phi_ = static_cast<float >(bad_val_);
  mumu_pt_ = static_cast<float >(bad_val_);
  mumu_pt1_ = static_cast<float >(bad_val_);
  mumu_pt2_ = static_cast<float >(bad_val_);
  mumu_w_ = static_cast<float >(bad_val_);
  mumuv_eta_ = static_cast<float >(bad_val_);
  mumuv_m_ = static_cast<float >(bad_val_);
  mumuv_phi_ = static_cast<float >(bad_val_);
  mumuv_pt_ = static_cast<float >(bad_val_);
  mumuv_pt1_ = static_cast<float >(bad_val_);
  mumuv_pt2_ = static_cast<float >(bad_val_);
  mumuv_w_ = static_cast<float >(bad_val_);
  ntrupv_mean_ = static_cast<float >(bad_val_);
  onel_ele105_ = static_cast<float >(bad_val_);
  onel_ele23_ = static_cast<float >(bad_val_);
  onel_ele8_ = static_cast<float >(bad_val_);
  onel_vvvl_ = static_cast<float >(bad_val_);
  onht_ = static_cast<float >(bad_val_);
  onmet_ = static_cast<float >(bad_val_);
  onmu_isomu18_ = static_cast<float >(bad_val_);
  onmu_mu50_ = static_cast<float >(bad_val_);
  onmu_mu8_ = static_cast<float >(bad_val_);
  onmu_vvvl_ = static_cast<float >(bad_val_);
  onph_ph90_ = static_cast<float >(bad_val_);
  w_btag_ = static_cast<float >(bad_val_);
  w_fs_lep_ = static_cast<float >(bad_val_);
  w_lep_ = static_cast<float >(bad_val_);
  w_lumi_ = static_cast<float >(bad_val_);
  w_pu_ = static_cast<float >(bad_val_);
  w_toppt_ = static_cast<float >(bad_val_);
  weight_ = static_cast<float >(bad_val_);
  lumiblock_ = static_cast<int >(bad_val_);
  mgluino_ = static_cast<int >(bad_val_);
  mlsp_ = static_cast<int >(bad_val_);
  nbl_ = static_cast<int >(bad_val_);
  nbm_ = static_cast<int >(bad_val_);
  nbm_ra2_ = static_cast<int >(bad_val_);
  nbt_ = static_cast<int >(bad_val_);
  nels_ = static_cast<int >(bad_val_);
  nels_ele23_ = static_cast<int >(bad_val_);
  nels_vvvl_ = static_cast<int >(bad_val_);
  nfjets_ = static_cast<int >(bad_val_);
  nfjets08_ = static_cast<int >(bad_val_);
  nisr_me_ = static_cast<int >(bad_val_);
  njets_ = static_cast<int >(bad_val_);
  njets_clean_ = static_cast<int >(bad_val_);
  njets_ra2_ = static_cast<int >(bad_val_);
  nleps_ = static_cast<int >(bad_val_);
  nleps_tm_ = static_cast<int >(bad_val_);
  nmus_ = static_cast<int >(bad_val_);
  nmus_isomu18_ = static_cast<int >(bad_val_);
  nmus_vvvl_ = static_cast<int >(bad_val_);
  nph_ = static_cast<int >(bad_val_);
  npv_ = static_cast<int >(bad_val_);
  ntks_ = static_cast<int >(bad_val_);
  ntruels_ = static_cast<int >(bad_val_);
  ntruleps_ = static_cast<int >(bad_val_);
  ntrumus_ = static_cast<int >(bad_val_);
  ntrupv_ = static_cast<int >(bad_val_);
  ntrutaush_ = static_cast<int >(bad_val_);
  ntrutausl_ = static_cast<int >(bad_val_);
  nvels_ = static_cast<int >(bad_val_);
  nveto_ = static_cast<int >(bad_val_);
  nvleps_ = static_cast<int >(bad_val_);
  nvmus_ = static_cast<int >(bad_val_);
  run_ = static_cast<int >(bad_val_);
  els_ele105_.clear();
  els_ele23_.clear();
  els_ele8_.clear();
  els_inz_.clear();
  els_inzv_.clear();
  els_ispf_.clear();
  els_sig_.clear();
  els_sigid_.clear();
  els_tight_.clear();
  els_tm_.clear();
  els_vvvl_.clear();
  jets_islep_.clear();
  mus_inz_.clear();
  mus_inzv_.clear();
  mus_isomu18_.clear();
  mus_mu50_.clear();
  mus_mu8_.clear();
  mus_sig_.clear();
  mus_sigid_.clear();
  mus_tight_.clear();
  mus_tm_.clear();
  mus_vvvl_.clear();
  ph_ph90_.clear();
  ph_tm_.clear();
  sys_pass_.clear();
  trig_.clear();
  els_d0_.clear();
  els_dz_.clear();
  els_eta_.clear();
  els_miniso_.clear();
  els_phi_.clear();
  els_pt_.clear();
  els_reliso_.clear();
  els_sceta_.clear();
  els_vvvl_eta_.clear();
  els_vvvl_phi_.clear();
  els_vvvl_pt_.clear();
  fjets08_eta_.clear();
  fjets08_m_.clear();
  fjets08_phi_.clear();
  fjets08_poscsv_.clear();
  fjets08_pt_.clear();
  fjets08_sumcsv_.clear();
  fjets_eta_.clear();
  fjets_m_.clear();
  fjets_phi_.clear();
  fjets_poscsv_.clear();
  fjets_pt_.clear();
  fjets_sumcsv_.clear();
  jets_csv_.clear();
  jets_eta_.clear();
  jets_m_.clear();
  jets_phi_.clear();
  jets_pt_.clear();
  jets_pt_res_.clear();
  leps_eta_.clear();
  leps_id_.clear();
  leps_phi_.clear();
  leps_pt_.clear();
  mc_eta_.clear();
  mc_mass_.clear();
  mc_phi_.clear();
  mc_pt_.clear();
  mus_d0_.clear();
  mus_dz_.clear();
  mus_eta_.clear();
  mus_miniso_.clear();
  mus_phi_.clear();
  mus_pt_.clear();
  mus_reliso_.clear();
  mus_vvvl_eta_.clear();
  mus_vvvl_phi_.clear();
  mus_vvvl_pt_.clear();
  ph_eta_.clear();
  ph_phi_.clear();
  ph_pt_.clear();
  sys_bctag_.clear();
  sys_fs_bctag_.clear();
  sys_fs_lep_.clear();
  sys_fs_udsgtag_.clear();
  sys_ht_.clear();
  sys_isr_.clear();
  sys_lep_.clear();
  sys_met_.clear();
  sys_mj_.clear();
  sys_mt_.clear();
  sys_muf_.clear();
  sys_mur_.clear();
  sys_murf_.clear();
  sys_pdf_.clear();
  sys_trig_.clear();
  sys_udsgtag_.clear();
  tks_eta_.clear();
  tks_miniso_.clear();
  tks_mt_.clear();
  tks_mt2_.clear();
  tks_phi_.clear();
  tks_pt_.clear();
  trig_prescale_.clear();
  w_pdf_.clear();
  els_charge_.clear();
  fjets08_btags_.clear();
  fjets08_nconst_.clear();
  fjets_btags_.clear();
  fjets_nconst_.clear();
  jets_fjet08_index_.clear();
  jets_fjet_index_.clear();
  jets_hflavor_.clear();
  mc_id_.clear();
  mc_mom_.clear();
  mus_charge_.clear();
  sys_nbm_.clear();
  sys_njets_.clear();
  tks_pdg_.clear();
  }

}

void baby_base::Write(){
  if(read_only_){
    throw std::logic_error("Trying to write to read-only tree.");
  }else{
    tree_.Write();
  }
}

string baby_base::Type() const{
  return "";
}

baby_base::~baby_base(){
}

int baby_base::Add(const std::string &filename){
  if(!read_only_){
    throw std::logic_error("Trying to add files to tree opened for writing.");
  }
  return chain_.Add(filename.c_str());
}

bool baby_base::PassString(TString cut){
  bool result = (cut=="dummy");
 return result;
}

long baby_base::GetEntries() const{
  if(read_only_){
    return chain_.GetEntries();
  }else{
    return tree_.GetEntries();
  }
}

void baby_base::GetEntry(const long entry){
  if(!read_only_){
    throw std::logic_error("Trying to read from write-only tree.");
  }

  c_event_ = false;
  c_fromGS_ = false;
  c_jetmismeas_ = false;
  c_low_dphi_ = false;
  c_nonblind_ = false;
  c_pass_ = false;
  c_pass_cschalo_ = false;
  c_pass_eebadsc_ = false;
  c_pass_goodv_ = false;
  c_pass_hbhe_ = false;
  c_pass_hbheiso_ = false;
  c_pass_jets_ = false;
  c_pass_jets_nohf_ = false;
  c_pass_jets_ra2_ = false;
  c_pass_jets_tight_ = false;
  c_pass_jets_tight_ra2_ = false;
  c_pass_nohf_ = false;
  c_pass_ra2_ = false;
  c_stitch_ = false;
  c_dphi1_ = false;
  c_dphi2_ = false;
  c_dphi3_ = false;
  c_dphi4_ = false;
  c_dphi_wlep_ = false;
  c_eff_trig_ = false;
  c_elel_eta_ = false;
  c_elel_m_ = false;
  c_elel_phi_ = false;
  c_elel_pt_ = false;
  c_elel_pt1_ = false;
  c_elel_pt2_ = false;
  c_elel_w_ = false;
  c_elelv_eta_ = false;
  c_elelv_m_ = false;
  c_elelv_phi_ = false;
  c_elelv_pt_ = false;
  c_elelv_pt1_ = false;
  c_elelv_pt2_ = false;
  c_elelv_w_ = false;
  c_elmu_eta_ = false;
  c_elmu_m_ = false;
  c_elmu_phi_ = false;
  c_elmu_pt_ = false;
  c_elmu_pt1_ = false;
  c_elmu_pt2_ = false;
  c_elmu_w_ = false;
  c_ht_ = false;
  c_ht_clean_ = false;
  c_ht_hlt_ = false;
  c_ht_isr_me_ = false;
  c_ht_ra2_ = false;
  c_ht_tru_ = false;
  c_isr_tru_eta_ = false;
  c_isr_tru_phi_ = false;
  c_isr_tru_pt_ = false;
  c_jetsys_eta_ = false;
  c_jetsys_nob_eta_ = false;
  c_jetsys_nob_phi_ = false;
  c_jetsys_nob_pt_ = false;
  c_jetsys_phi_ = false;
  c_jetsys_pt_ = false;
  c_met_ = false;
  c_met_calo_ = false;
  c_met_calo_phi_ = false;
  c_met_mini_ = false;
  c_met_mini_phi_ = false;
  c_met_nohf_ = false;
  c_met_nohf_phi_ = false;
  c_met_phi_ = false;
  c_met_raw_ = false;
  c_met_raw_phi_ = false;
  c_met_rebal_ = false;
  c_met_tru_ = false;
  c_met_tru_nuw_ = false;
  c_met_tru_nuw_phi_ = false;
  c_met_tru_phi_ = false;
  c_mht_ = false;
  c_mht_clean_ = false;
  c_mht_clean_phi_ = false;
  c_mht_phi_ = false;
  c_mj_ = false;
  c_mj08_ = false;
  c_mt_ = false;
  c_mt_nohf_ = false;
  c_mt_rebal_ = false;
  c_mt_tru_ = false;
  c_mt_tru_nuw_ = false;
  c_mumu_eta_ = false;
  c_mumu_m_ = false;
  c_mumu_phi_ = false;
  c_mumu_pt_ = false;
  c_mumu_pt1_ = false;
  c_mumu_pt2_ = false;
  c_mumu_w_ = false;
  c_mumuv_eta_ = false;
  c_mumuv_m_ = false;
  c_mumuv_phi_ = false;
  c_mumuv_pt_ = false;
  c_mumuv_pt1_ = false;
  c_mumuv_pt2_ = false;
  c_mumuv_w_ = false;
  c_ntrupv_mean_ = false;
  c_onel_ele105_ = false;
  c_onel_ele23_ = false;
  c_onel_ele8_ = false;
  c_onel_vvvl_ = false;
  c_onht_ = false;
  c_onmet_ = false;
  c_onmu_isomu18_ = false;
  c_onmu_mu50_ = false;
  c_onmu_mu8_ = false;
  c_onmu_vvvl_ = false;
  c_onph_ph90_ = false;
  c_w_btag_ = false;
  c_w_fs_lep_ = false;
  c_w_lep_ = false;
  c_w_lumi_ = false;
  c_w_pu_ = false;
  c_w_toppt_ = false;
  c_weight_ = false;
  c_lumiblock_ = false;
  c_mgluino_ = false;
  c_mlsp_ = false;
  c_nbl_ = false;
  c_nbm_ = false;
  c_nbm_ra2_ = false;
  c_nbt_ = false;
  c_nels_ = false;
  c_nels_ele23_ = false;
  c_nels_vvvl_ = false;
  c_nfjets_ = false;
  c_nfjets08_ = false;
  c_nisr_me_ = false;
  c_njets_ = false;
  c_njets_clean_ = false;
  c_njets_ra2_ = false;
  c_nleps_ = false;
  c_nleps_tm_ = false;
  c_nmus_ = false;
  c_nmus_isomu18_ = false;
  c_nmus_vvvl_ = false;
  c_nph_ = false;
  c_npv_ = false;
  c_ntks_ = false;
  c_ntruels_ = false;
  c_ntruleps_ = false;
  c_ntrumus_ = false;
  c_ntrupv_ = false;
  c_ntrutaush_ = false;
  c_ntrutausl_ = false;
  c_nvels_ = false;
  c_nveto_ = false;
  c_nvleps_ = false;
  c_nvmus_ = false;
  c_run_ = false;
  c_els_ele105_ = false;
  c_els_ele23_ = false;
  c_els_ele8_ = false;
  c_els_inz_ = false;
  c_els_inzv_ = false;
  c_els_ispf_ = false;
  c_els_sig_ = false;
  c_els_sigid_ = false;
  c_els_tight_ = false;
  c_els_tm_ = false;
  c_els_vvvl_ = false;
  c_jets_islep_ = false;
  c_mus_inz_ = false;
  c_mus_inzv_ = false;
  c_mus_isomu18_ = false;
  c_mus_mu50_ = false;
  c_mus_mu8_ = false;
  c_mus_sig_ = false;
  c_mus_sigid_ = false;
  c_mus_tight_ = false;
  c_mus_tm_ = false;
  c_mus_vvvl_ = false;
  c_ph_ph90_ = false;
  c_ph_tm_ = false;
  c_sys_pass_ = false;
  c_trig_ = false;
  c_els_d0_ = false;
  c_els_dz_ = false;
  c_els_eta_ = false;
  c_els_miniso_ = false;
  c_els_phi_ = false;
  c_els_pt_ = false;
  c_els_reliso_ = false;
  c_els_sceta_ = false;
  c_els_vvvl_eta_ = false;
  c_els_vvvl_phi_ = false;
  c_els_vvvl_pt_ = false;
  c_fjets08_eta_ = false;
  c_fjets08_m_ = false;
  c_fjets08_phi_ = false;
  c_fjets08_poscsv_ = false;
  c_fjets08_pt_ = false;
  c_fjets08_sumcsv_ = false;
  c_fjets_eta_ = false;
  c_fjets_m_ = false;
  c_fjets_phi_ = false;
  c_fjets_poscsv_ = false;
  c_fjets_pt_ = false;
  c_fjets_sumcsv_ = false;
  c_jets_csv_ = false;
  c_jets_eta_ = false;
  c_jets_m_ = false;
  c_jets_phi_ = false;
  c_jets_pt_ = false;
  c_jets_pt_res_ = false;
  c_leps_eta_ = false;
  c_leps_id_ = false;
  c_leps_phi_ = false;
  c_leps_pt_ = false;
  c_mc_eta_ = false;
  c_mc_mass_ = false;
  c_mc_phi_ = false;
  c_mc_pt_ = false;
  c_mus_d0_ = false;
  c_mus_dz_ = false;
  c_mus_eta_ = false;
  c_mus_miniso_ = false;
  c_mus_phi_ = false;
  c_mus_pt_ = false;
  c_mus_reliso_ = false;
  c_mus_vvvl_eta_ = false;
  c_mus_vvvl_phi_ = false;
  c_mus_vvvl_pt_ = false;
  c_ph_eta_ = false;
  c_ph_phi_ = false;
  c_ph_pt_ = false;
  c_sys_bctag_ = false;
  c_sys_fs_bctag_ = false;
  c_sys_fs_lep_ = false;
  c_sys_fs_udsgtag_ = false;
  c_sys_ht_ = false;
  c_sys_isr_ = false;
  c_sys_lep_ = false;
  c_sys_met_ = false;
  c_sys_mj_ = false;
  c_sys_mt_ = false;
  c_sys_muf_ = false;
  c_sys_mur_ = false;
  c_sys_murf_ = false;
  c_sys_pdf_ = false;
  c_sys_trig_ = false;
  c_sys_udsgtag_ = false;
  c_tks_eta_ = false;
  c_tks_miniso_ = false;
  c_tks_mt_ = false;
  c_tks_mt2_ = false;
  c_tks_phi_ = false;
  c_tks_pt_ = false;
  c_trig_prescale_ = false;
  c_w_pdf_ = false;
  c_els_charge_ = false;
  c_fjets08_btags_ = false;
  c_fjets08_nconst_ = false;
  c_fjets_btags_ = false;
  c_fjets_nconst_ = false;
  c_jets_fjet08_index_ = false;
  c_jets_fjet_index_ = false;
  c_jets_hflavor_ = false;
  c_mc_id_ = false;
  c_mc_mom_ = false;
  c_mus_charge_ = false;
  c_sys_nbm_ = false;
  c_sys_njets_ = false;
  c_tks_pdg_ = false;
  entry_ = chain_.LoadTree(entry);
}

Long64_t  const & baby_base::event() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_event_ && b_event_){
    b_event_->GetEntry(entry_);
    c_event_ = true;
  }
  return event_;
}

bool  const & baby_base::fromGS() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fromGS_ && b_fromGS_){
    b_fromGS_->GetEntry(entry_);
    c_fromGS_ = true;
  }
  return fromGS_;
}

bool  const & baby_base::jetmismeas() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jetmismeas_ && b_jetmismeas_){
    b_jetmismeas_->GetEntry(entry_);
    c_jetmismeas_ = true;
  }
  return jetmismeas_;
}

bool  const & baby_base::low_dphi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_low_dphi_ && b_low_dphi_){
    b_low_dphi_->GetEntry(entry_);
    c_low_dphi_ = true;
  }
  return low_dphi_;
}

bool  const & baby_base::nonblind() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nonblind_ && b_nonblind_){
    b_nonblind_->GetEntry(entry_);
    c_nonblind_ = true;
  }
  return nonblind_;
}

bool  const & baby_base::pass() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_ && b_pass_){
    b_pass_->GetEntry(entry_);
    c_pass_ = true;
  }
  return pass_;
}

bool  const & baby_base::pass_cschalo() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_cschalo_ && b_pass_cschalo_){
    b_pass_cschalo_->GetEntry(entry_);
    c_pass_cschalo_ = true;
  }
  return pass_cschalo_;
}

bool  const & baby_base::pass_eebadsc() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_eebadsc_ && b_pass_eebadsc_){
    b_pass_eebadsc_->GetEntry(entry_);
    c_pass_eebadsc_ = true;
  }
  return pass_eebadsc_;
}

bool  const & baby_base::pass_goodv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_goodv_ && b_pass_goodv_){
    b_pass_goodv_->GetEntry(entry_);
    c_pass_goodv_ = true;
  }
  return pass_goodv_;
}

bool  const & baby_base::pass_hbhe() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_hbhe_ && b_pass_hbhe_){
    b_pass_hbhe_->GetEntry(entry_);
    c_pass_hbhe_ = true;
  }
  return pass_hbhe_;
}

bool  const & baby_base::pass_hbheiso() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_hbheiso_ && b_pass_hbheiso_){
    b_pass_hbheiso_->GetEntry(entry_);
    c_pass_hbheiso_ = true;
  }
  return pass_hbheiso_;
}

bool  const & baby_base::pass_jets() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_jets_ && b_pass_jets_){
    b_pass_jets_->GetEntry(entry_);
    c_pass_jets_ = true;
  }
  return pass_jets_;
}

bool  const & baby_base::pass_jets_nohf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_jets_nohf_ && b_pass_jets_nohf_){
    b_pass_jets_nohf_->GetEntry(entry_);
    c_pass_jets_nohf_ = true;
  }
  return pass_jets_nohf_;
}

bool  const & baby_base::pass_jets_ra2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_jets_ra2_ && b_pass_jets_ra2_){
    b_pass_jets_ra2_->GetEntry(entry_);
    c_pass_jets_ra2_ = true;
  }
  return pass_jets_ra2_;
}

bool  const & baby_base::pass_jets_tight() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_jets_tight_ && b_pass_jets_tight_){
    b_pass_jets_tight_->GetEntry(entry_);
    c_pass_jets_tight_ = true;
  }
  return pass_jets_tight_;
}

bool  const & baby_base::pass_jets_tight_ra2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_jets_tight_ra2_ && b_pass_jets_tight_ra2_){
    b_pass_jets_tight_ra2_->GetEntry(entry_);
    c_pass_jets_tight_ra2_ = true;
  }
  return pass_jets_tight_ra2_;
}

bool  const & baby_base::pass_nohf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_nohf_ && b_pass_nohf_){
    b_pass_nohf_->GetEntry(entry_);
    c_pass_nohf_ = true;
  }
  return pass_nohf_;
}

bool  const & baby_base::pass_ra2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_pass_ra2_ && b_pass_ra2_){
    b_pass_ra2_->GetEntry(entry_);
    c_pass_ra2_ = true;
  }
  return pass_ra2_;
}

bool  const & baby_base::stitch() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_stitch_ && b_stitch_){
    b_stitch_->GetEntry(entry_);
    c_stitch_ = true;
  }
  return stitch_;
}

float  const & baby_base::dphi1() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_dphi1_ && b_dphi1_){
    b_dphi1_->GetEntry(entry_);
    c_dphi1_ = true;
  }
  return dphi1_;
}

float  const & baby_base::dphi2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_dphi2_ && b_dphi2_){
    b_dphi2_->GetEntry(entry_);
    c_dphi2_ = true;
  }
  return dphi2_;
}

float  const & baby_base::dphi3() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_dphi3_ && b_dphi3_){
    b_dphi3_->GetEntry(entry_);
    c_dphi3_ = true;
  }
  return dphi3_;
}

float  const & baby_base::dphi4() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_dphi4_ && b_dphi4_){
    b_dphi4_->GetEntry(entry_);
    c_dphi4_ = true;
  }
  return dphi4_;
}

float  const & baby_base::dphi_wlep() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_dphi_wlep_ && b_dphi_wlep_){
    b_dphi_wlep_->GetEntry(entry_);
    c_dphi_wlep_ = true;
  }
  return dphi_wlep_;
}

float  const & baby_base::eff_trig() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_eff_trig_ && b_eff_trig_){
    b_eff_trig_->GetEntry(entry_);
    c_eff_trig_ = true;
  }
  return eff_trig_;
}

float  const & baby_base::elel_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elel_eta_ && b_elel_eta_){
    b_elel_eta_->GetEntry(entry_);
    c_elel_eta_ = true;
  }
  return elel_eta_;
}

float  const & baby_base::elel_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elel_m_ && b_elel_m_){
    b_elel_m_->GetEntry(entry_);
    c_elel_m_ = true;
  }
  return elel_m_;
}

float  const & baby_base::elel_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elel_phi_ && b_elel_phi_){
    b_elel_phi_->GetEntry(entry_);
    c_elel_phi_ = true;
  }
  return elel_phi_;
}

float  const & baby_base::elel_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elel_pt_ && b_elel_pt_){
    b_elel_pt_->GetEntry(entry_);
    c_elel_pt_ = true;
  }
  return elel_pt_;
}

float  const & baby_base::elel_pt1() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elel_pt1_ && b_elel_pt1_){
    b_elel_pt1_->GetEntry(entry_);
    c_elel_pt1_ = true;
  }
  return elel_pt1_;
}

float  const & baby_base::elel_pt2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elel_pt2_ && b_elel_pt2_){
    b_elel_pt2_->GetEntry(entry_);
    c_elel_pt2_ = true;
  }
  return elel_pt2_;
}

float  const & baby_base::elel_w() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elel_w_ && b_elel_w_){
    b_elel_w_->GetEntry(entry_);
    c_elel_w_ = true;
  }
  return elel_w_;
}

float  const & baby_base::elelv_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elelv_eta_ && b_elelv_eta_){
    b_elelv_eta_->GetEntry(entry_);
    c_elelv_eta_ = true;
  }
  return elelv_eta_;
}

float  const & baby_base::elelv_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elelv_m_ && b_elelv_m_){
    b_elelv_m_->GetEntry(entry_);
    c_elelv_m_ = true;
  }
  return elelv_m_;
}

float  const & baby_base::elelv_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elelv_phi_ && b_elelv_phi_){
    b_elelv_phi_->GetEntry(entry_);
    c_elelv_phi_ = true;
  }
  return elelv_phi_;
}

float  const & baby_base::elelv_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elelv_pt_ && b_elelv_pt_){
    b_elelv_pt_->GetEntry(entry_);
    c_elelv_pt_ = true;
  }
  return elelv_pt_;
}

float  const & baby_base::elelv_pt1() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elelv_pt1_ && b_elelv_pt1_){
    b_elelv_pt1_->GetEntry(entry_);
    c_elelv_pt1_ = true;
  }
  return elelv_pt1_;
}

float  const & baby_base::elelv_pt2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elelv_pt2_ && b_elelv_pt2_){
    b_elelv_pt2_->GetEntry(entry_);
    c_elelv_pt2_ = true;
  }
  return elelv_pt2_;
}

float  const & baby_base::elelv_w() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elelv_w_ && b_elelv_w_){
    b_elelv_w_->GetEntry(entry_);
    c_elelv_w_ = true;
  }
  return elelv_w_;
}

float  const & baby_base::elmu_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elmu_eta_ && b_elmu_eta_){
    b_elmu_eta_->GetEntry(entry_);
    c_elmu_eta_ = true;
  }
  return elmu_eta_;
}

float  const & baby_base::elmu_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elmu_m_ && b_elmu_m_){
    b_elmu_m_->GetEntry(entry_);
    c_elmu_m_ = true;
  }
  return elmu_m_;
}

float  const & baby_base::elmu_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elmu_phi_ && b_elmu_phi_){
    b_elmu_phi_->GetEntry(entry_);
    c_elmu_phi_ = true;
  }
  return elmu_phi_;
}

float  const & baby_base::elmu_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elmu_pt_ && b_elmu_pt_){
    b_elmu_pt_->GetEntry(entry_);
    c_elmu_pt_ = true;
  }
  return elmu_pt_;
}

float  const & baby_base::elmu_pt1() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elmu_pt1_ && b_elmu_pt1_){
    b_elmu_pt1_->GetEntry(entry_);
    c_elmu_pt1_ = true;
  }
  return elmu_pt1_;
}

float  const & baby_base::elmu_pt2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elmu_pt2_ && b_elmu_pt2_){
    b_elmu_pt2_->GetEntry(entry_);
    c_elmu_pt2_ = true;
  }
  return elmu_pt2_;
}

float  const & baby_base::elmu_w() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_elmu_w_ && b_elmu_w_){
    b_elmu_w_->GetEntry(entry_);
    c_elmu_w_ = true;
  }
  return elmu_w_;
}

float  const & baby_base::ht() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ht_ && b_ht_){
    b_ht_->GetEntry(entry_);
    c_ht_ = true;
  }
  return ht_;
}

float  const & baby_base::ht_clean() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ht_clean_ && b_ht_clean_){
    b_ht_clean_->GetEntry(entry_);
    c_ht_clean_ = true;
  }
  return ht_clean_;
}

float  const & baby_base::ht_hlt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ht_hlt_ && b_ht_hlt_){
    b_ht_hlt_->GetEntry(entry_);
    c_ht_hlt_ = true;
  }
  return ht_hlt_;
}

float  const & baby_base::ht_isr_me() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ht_isr_me_ && b_ht_isr_me_){
    b_ht_isr_me_->GetEntry(entry_);
    c_ht_isr_me_ = true;
  }
  return ht_isr_me_;
}

float  const & baby_base::ht_ra2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ht_ra2_ && b_ht_ra2_){
    b_ht_ra2_->GetEntry(entry_);
    c_ht_ra2_ = true;
  }
  return ht_ra2_;
}

float  const & baby_base::ht_tru() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ht_tru_ && b_ht_tru_){
    b_ht_tru_->GetEntry(entry_);
    c_ht_tru_ = true;
  }
  return ht_tru_;
}

float  const & baby_base::isr_tru_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_isr_tru_eta_ && b_isr_tru_eta_){
    b_isr_tru_eta_->GetEntry(entry_);
    c_isr_tru_eta_ = true;
  }
  return isr_tru_eta_;
}

float  const & baby_base::isr_tru_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_isr_tru_phi_ && b_isr_tru_phi_){
    b_isr_tru_phi_->GetEntry(entry_);
    c_isr_tru_phi_ = true;
  }
  return isr_tru_phi_;
}

float  const & baby_base::isr_tru_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_isr_tru_pt_ && b_isr_tru_pt_){
    b_isr_tru_pt_->GetEntry(entry_);
    c_isr_tru_pt_ = true;
  }
  return isr_tru_pt_;
}

float  const & baby_base::jetsys_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jetsys_eta_ && b_jetsys_eta_){
    b_jetsys_eta_->GetEntry(entry_);
    c_jetsys_eta_ = true;
  }
  return jetsys_eta_;
}

float  const & baby_base::jetsys_nob_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jetsys_nob_eta_ && b_jetsys_nob_eta_){
    b_jetsys_nob_eta_->GetEntry(entry_);
    c_jetsys_nob_eta_ = true;
  }
  return jetsys_nob_eta_;
}

float  const & baby_base::jetsys_nob_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jetsys_nob_phi_ && b_jetsys_nob_phi_){
    b_jetsys_nob_phi_->GetEntry(entry_);
    c_jetsys_nob_phi_ = true;
  }
  return jetsys_nob_phi_;
}

float  const & baby_base::jetsys_nob_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jetsys_nob_pt_ && b_jetsys_nob_pt_){
    b_jetsys_nob_pt_->GetEntry(entry_);
    c_jetsys_nob_pt_ = true;
  }
  return jetsys_nob_pt_;
}

float  const & baby_base::jetsys_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jetsys_phi_ && b_jetsys_phi_){
    b_jetsys_phi_->GetEntry(entry_);
    c_jetsys_phi_ = true;
  }
  return jetsys_phi_;
}

float  const & baby_base::jetsys_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jetsys_pt_ && b_jetsys_pt_){
    b_jetsys_pt_->GetEntry(entry_);
    c_jetsys_pt_ = true;
  }
  return jetsys_pt_;
}

float  const & baby_base::met() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_ && b_met_){
    b_met_->GetEntry(entry_);
    c_met_ = true;
  }
  return met_;
}

float  const & baby_base::met_calo() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_calo_ && b_met_calo_){
    b_met_calo_->GetEntry(entry_);
    c_met_calo_ = true;
  }
  return met_calo_;
}

float  const & baby_base::met_calo_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_calo_phi_ && b_met_calo_phi_){
    b_met_calo_phi_->GetEntry(entry_);
    c_met_calo_phi_ = true;
  }
  return met_calo_phi_;
}

float  const & baby_base::met_mini() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_mini_ && b_met_mini_){
    b_met_mini_->GetEntry(entry_);
    c_met_mini_ = true;
  }
  return met_mini_;
}

float  const & baby_base::met_mini_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_mini_phi_ && b_met_mini_phi_){
    b_met_mini_phi_->GetEntry(entry_);
    c_met_mini_phi_ = true;
  }
  return met_mini_phi_;
}

float  const & baby_base::met_nohf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_nohf_ && b_met_nohf_){
    b_met_nohf_->GetEntry(entry_);
    c_met_nohf_ = true;
  }
  return met_nohf_;
}

float  const & baby_base::met_nohf_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_nohf_phi_ && b_met_nohf_phi_){
    b_met_nohf_phi_->GetEntry(entry_);
    c_met_nohf_phi_ = true;
  }
  return met_nohf_phi_;
}

float  const & baby_base::met_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_phi_ && b_met_phi_){
    b_met_phi_->GetEntry(entry_);
    c_met_phi_ = true;
  }
  return met_phi_;
}

float  const & baby_base::met_raw() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_raw_ && b_met_raw_){
    b_met_raw_->GetEntry(entry_);
    c_met_raw_ = true;
  }
  return met_raw_;
}

float  const & baby_base::met_raw_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_raw_phi_ && b_met_raw_phi_){
    b_met_raw_phi_->GetEntry(entry_);
    c_met_raw_phi_ = true;
  }
  return met_raw_phi_;
}

float  const & baby_base::met_rebal() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_rebal_ && b_met_rebal_){
    b_met_rebal_->GetEntry(entry_);
    c_met_rebal_ = true;
  }
  return met_rebal_;
}

float  const & baby_base::met_tru() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_tru_ && b_met_tru_){
    b_met_tru_->GetEntry(entry_);
    c_met_tru_ = true;
  }
  return met_tru_;
}

float  const & baby_base::met_tru_nuw() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_tru_nuw_ && b_met_tru_nuw_){
    b_met_tru_nuw_->GetEntry(entry_);
    c_met_tru_nuw_ = true;
  }
  return met_tru_nuw_;
}

float  const & baby_base::met_tru_nuw_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_tru_nuw_phi_ && b_met_tru_nuw_phi_){
    b_met_tru_nuw_phi_->GetEntry(entry_);
    c_met_tru_nuw_phi_ = true;
  }
  return met_tru_nuw_phi_;
}

float  const & baby_base::met_tru_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_met_tru_phi_ && b_met_tru_phi_){
    b_met_tru_phi_->GetEntry(entry_);
    c_met_tru_phi_ = true;
  }
  return met_tru_phi_;
}

float  const & baby_base::mht() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mht_ && b_mht_){
    b_mht_->GetEntry(entry_);
    c_mht_ = true;
  }
  return mht_;
}

float  const & baby_base::mht_clean() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mht_clean_ && b_mht_clean_){
    b_mht_clean_->GetEntry(entry_);
    c_mht_clean_ = true;
  }
  return mht_clean_;
}

float  const & baby_base::mht_clean_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mht_clean_phi_ && b_mht_clean_phi_){
    b_mht_clean_phi_->GetEntry(entry_);
    c_mht_clean_phi_ = true;
  }
  return mht_clean_phi_;
}

float  const & baby_base::mht_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mht_phi_ && b_mht_phi_){
    b_mht_phi_->GetEntry(entry_);
    c_mht_phi_ = true;
  }
  return mht_phi_;
}

float  const & baby_base::mj() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mj_ && b_mj_){
    b_mj_->GetEntry(entry_);
    c_mj_ = true;
  }
  return mj_;
}

float  const & baby_base::mj08() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mj08_ && b_mj08_){
    b_mj08_->GetEntry(entry_);
    c_mj08_ = true;
  }
  return mj08_;
}

float  const & baby_base::mt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mt_ && b_mt_){
    b_mt_->GetEntry(entry_);
    c_mt_ = true;
  }
  return mt_;
}

float  const & baby_base::mt_nohf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mt_nohf_ && b_mt_nohf_){
    b_mt_nohf_->GetEntry(entry_);
    c_mt_nohf_ = true;
  }
  return mt_nohf_;
}

float  const & baby_base::mt_rebal() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mt_rebal_ && b_mt_rebal_){
    b_mt_rebal_->GetEntry(entry_);
    c_mt_rebal_ = true;
  }
  return mt_rebal_;
}

float  const & baby_base::mt_tru() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mt_tru_ && b_mt_tru_){
    b_mt_tru_->GetEntry(entry_);
    c_mt_tru_ = true;
  }
  return mt_tru_;
}

float  const & baby_base::mt_tru_nuw() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mt_tru_nuw_ && b_mt_tru_nuw_){
    b_mt_tru_nuw_->GetEntry(entry_);
    c_mt_tru_nuw_ = true;
  }
  return mt_tru_nuw_;
}

float  const & baby_base::mumu_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumu_eta_ && b_mumu_eta_){
    b_mumu_eta_->GetEntry(entry_);
    c_mumu_eta_ = true;
  }
  return mumu_eta_;
}

float  const & baby_base::mumu_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumu_m_ && b_mumu_m_){
    b_mumu_m_->GetEntry(entry_);
    c_mumu_m_ = true;
  }
  return mumu_m_;
}

float  const & baby_base::mumu_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumu_phi_ && b_mumu_phi_){
    b_mumu_phi_->GetEntry(entry_);
    c_mumu_phi_ = true;
  }
  return mumu_phi_;
}

float  const & baby_base::mumu_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumu_pt_ && b_mumu_pt_){
    b_mumu_pt_->GetEntry(entry_);
    c_mumu_pt_ = true;
  }
  return mumu_pt_;
}

float  const & baby_base::mumu_pt1() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumu_pt1_ && b_mumu_pt1_){
    b_mumu_pt1_->GetEntry(entry_);
    c_mumu_pt1_ = true;
  }
  return mumu_pt1_;
}

float  const & baby_base::mumu_pt2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumu_pt2_ && b_mumu_pt2_){
    b_mumu_pt2_->GetEntry(entry_);
    c_mumu_pt2_ = true;
  }
  return mumu_pt2_;
}

float  const & baby_base::mumu_w() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumu_w_ && b_mumu_w_){
    b_mumu_w_->GetEntry(entry_);
    c_mumu_w_ = true;
  }
  return mumu_w_;
}

float  const & baby_base::mumuv_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumuv_eta_ && b_mumuv_eta_){
    b_mumuv_eta_->GetEntry(entry_);
    c_mumuv_eta_ = true;
  }
  return mumuv_eta_;
}

float  const & baby_base::mumuv_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumuv_m_ && b_mumuv_m_){
    b_mumuv_m_->GetEntry(entry_);
    c_mumuv_m_ = true;
  }
  return mumuv_m_;
}

float  const & baby_base::mumuv_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumuv_phi_ && b_mumuv_phi_){
    b_mumuv_phi_->GetEntry(entry_);
    c_mumuv_phi_ = true;
  }
  return mumuv_phi_;
}

float  const & baby_base::mumuv_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumuv_pt_ && b_mumuv_pt_){
    b_mumuv_pt_->GetEntry(entry_);
    c_mumuv_pt_ = true;
  }
  return mumuv_pt_;
}

float  const & baby_base::mumuv_pt1() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumuv_pt1_ && b_mumuv_pt1_){
    b_mumuv_pt1_->GetEntry(entry_);
    c_mumuv_pt1_ = true;
  }
  return mumuv_pt1_;
}

float  const & baby_base::mumuv_pt2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumuv_pt2_ && b_mumuv_pt2_){
    b_mumuv_pt2_->GetEntry(entry_);
    c_mumuv_pt2_ = true;
  }
  return mumuv_pt2_;
}

float  const & baby_base::mumuv_w() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mumuv_w_ && b_mumuv_w_){
    b_mumuv_w_->GetEntry(entry_);
    c_mumuv_w_ = true;
  }
  return mumuv_w_;
}

float  const & baby_base::ntrupv_mean() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntrupv_mean_ && b_ntrupv_mean_){
    b_ntrupv_mean_->GetEntry(entry_);
    c_ntrupv_mean_ = true;
  }
  return ntrupv_mean_;
}

float  const & baby_base::onel_ele105() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onel_ele105_ && b_onel_ele105_){
    b_onel_ele105_->GetEntry(entry_);
    c_onel_ele105_ = true;
  }
  return onel_ele105_;
}

float  const & baby_base::onel_ele23() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onel_ele23_ && b_onel_ele23_){
    b_onel_ele23_->GetEntry(entry_);
    c_onel_ele23_ = true;
  }
  return onel_ele23_;
}

float  const & baby_base::onel_ele8() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onel_ele8_ && b_onel_ele8_){
    b_onel_ele8_->GetEntry(entry_);
    c_onel_ele8_ = true;
  }
  return onel_ele8_;
}

float  const & baby_base::onel_vvvl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onel_vvvl_ && b_onel_vvvl_){
    b_onel_vvvl_->GetEntry(entry_);
    c_onel_vvvl_ = true;
  }
  return onel_vvvl_;
}

float  const & baby_base::onht() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onht_ && b_onht_){
    b_onht_->GetEntry(entry_);
    c_onht_ = true;
  }
  return onht_;
}

float  const & baby_base::onmet() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onmet_ && b_onmet_){
    b_onmet_->GetEntry(entry_);
    c_onmet_ = true;
  }
  return onmet_;
}

float  const & baby_base::onmu_isomu18() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onmu_isomu18_ && b_onmu_isomu18_){
    b_onmu_isomu18_->GetEntry(entry_);
    c_onmu_isomu18_ = true;
  }
  return onmu_isomu18_;
}

float  const & baby_base::onmu_mu50() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onmu_mu50_ && b_onmu_mu50_){
    b_onmu_mu50_->GetEntry(entry_);
    c_onmu_mu50_ = true;
  }
  return onmu_mu50_;
}

float  const & baby_base::onmu_mu8() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onmu_mu8_ && b_onmu_mu8_){
    b_onmu_mu8_->GetEntry(entry_);
    c_onmu_mu8_ = true;
  }
  return onmu_mu8_;
}

float  const & baby_base::onmu_vvvl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onmu_vvvl_ && b_onmu_vvvl_){
    b_onmu_vvvl_->GetEntry(entry_);
    c_onmu_vvvl_ = true;
  }
  return onmu_vvvl_;
}

float  const & baby_base::onph_ph90() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_onph_ph90_ && b_onph_ph90_){
    b_onph_ph90_->GetEntry(entry_);
    c_onph_ph90_ = true;
  }
  return onph_ph90_;
}

float  const & baby_base::w_btag() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_w_btag_ && b_w_btag_){
    b_w_btag_->GetEntry(entry_);
    c_w_btag_ = true;
  }
  return w_btag_;
}

float  const & baby_base::w_fs_lep() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_w_fs_lep_ && b_w_fs_lep_){
    b_w_fs_lep_->GetEntry(entry_);
    c_w_fs_lep_ = true;
  }
  return w_fs_lep_;
}

float  const & baby_base::w_lep() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_w_lep_ && b_w_lep_){
    b_w_lep_->GetEntry(entry_);
    c_w_lep_ = true;
  }
  return w_lep_;
}

float  const & baby_base::w_lumi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_w_lumi_ && b_w_lumi_){
    b_w_lumi_->GetEntry(entry_);
    c_w_lumi_ = true;
  }
  return w_lumi_;
}

float  const & baby_base::w_pu() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_w_pu_ && b_w_pu_){
    b_w_pu_->GetEntry(entry_);
    c_w_pu_ = true;
  }
  return w_pu_;
}

float  const & baby_base::w_toppt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_w_toppt_ && b_w_toppt_){
    b_w_toppt_->GetEntry(entry_);
    c_w_toppt_ = true;
  }
  return w_toppt_;
}

float  const & baby_base::weight() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_weight_ && b_weight_){
    b_weight_->GetEntry(entry_);
    c_weight_ = true;
  }
  return weight_;
}

int  const & baby_base::lumiblock() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_lumiblock_ && b_lumiblock_){
    b_lumiblock_->GetEntry(entry_);
    c_lumiblock_ = true;
  }
  return lumiblock_;
}

int  const & baby_base::mgluino() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mgluino_ && b_mgluino_){
    b_mgluino_->GetEntry(entry_);
    c_mgluino_ = true;
  }
  return mgluino_;
}

int  const & baby_base::mlsp() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mlsp_ && b_mlsp_){
    b_mlsp_->GetEntry(entry_);
    c_mlsp_ = true;
  }
  return mlsp_;
}

int  const & baby_base::nbl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nbl_ && b_nbl_){
    b_nbl_->GetEntry(entry_);
    c_nbl_ = true;
  }
  return nbl_;
}

int  const & baby_base::nbm() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nbm_ && b_nbm_){
    b_nbm_->GetEntry(entry_);
    c_nbm_ = true;
  }
  return nbm_;
}

int  const & baby_base::nbm_ra2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nbm_ra2_ && b_nbm_ra2_){
    b_nbm_ra2_->GetEntry(entry_);
    c_nbm_ra2_ = true;
  }
  return nbm_ra2_;
}

int  const & baby_base::nbt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nbt_ && b_nbt_){
    b_nbt_->GetEntry(entry_);
    c_nbt_ = true;
  }
  return nbt_;
}

int  const & baby_base::nels() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nels_ && b_nels_){
    b_nels_->GetEntry(entry_);
    c_nels_ = true;
  }
  return nels_;
}

int  const & baby_base::nels_ele23() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nels_ele23_ && b_nels_ele23_){
    b_nels_ele23_->GetEntry(entry_);
    c_nels_ele23_ = true;
  }
  return nels_ele23_;
}

int  const & baby_base::nels_vvvl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nels_vvvl_ && b_nels_vvvl_){
    b_nels_vvvl_->GetEntry(entry_);
    c_nels_vvvl_ = true;
  }
  return nels_vvvl_;
}

int  const & baby_base::nfjets() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nfjets_ && b_nfjets_){
    b_nfjets_->GetEntry(entry_);
    c_nfjets_ = true;
  }
  return nfjets_;
}

int  const & baby_base::nfjets08() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nfjets08_ && b_nfjets08_){
    b_nfjets08_->GetEntry(entry_);
    c_nfjets08_ = true;
  }
  return nfjets08_;
}

int  const & baby_base::nisr_me() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nisr_me_ && b_nisr_me_){
    b_nisr_me_->GetEntry(entry_);
    c_nisr_me_ = true;
  }
  return nisr_me_;
}

int  const & baby_base::njets() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_njets_ && b_njets_){
    b_njets_->GetEntry(entry_);
    c_njets_ = true;
  }
  return njets_;
}

int  const & baby_base::njets_clean() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_njets_clean_ && b_njets_clean_){
    b_njets_clean_->GetEntry(entry_);
    c_njets_clean_ = true;
  }
  return njets_clean_;
}

int  const & baby_base::njets_ra2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_njets_ra2_ && b_njets_ra2_){
    b_njets_ra2_->GetEntry(entry_);
    c_njets_ra2_ = true;
  }
  return njets_ra2_;
}

int  const & baby_base::nleps() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nleps_ && b_nleps_){
    b_nleps_->GetEntry(entry_);
    c_nleps_ = true;
  }
  return nleps_;
}

int  const & baby_base::nleps_tm() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nleps_tm_ && b_nleps_tm_){
    b_nleps_tm_->GetEntry(entry_);
    c_nleps_tm_ = true;
  }
  return nleps_tm_;
}

int  const & baby_base::nmus() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nmus_ && b_nmus_){
    b_nmus_->GetEntry(entry_);
    c_nmus_ = true;
  }
  return nmus_;
}

int  const & baby_base::nmus_isomu18() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nmus_isomu18_ && b_nmus_isomu18_){
    b_nmus_isomu18_->GetEntry(entry_);
    c_nmus_isomu18_ = true;
  }
  return nmus_isomu18_;
}

int  const & baby_base::nmus_vvvl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nmus_vvvl_ && b_nmus_vvvl_){
    b_nmus_vvvl_->GetEntry(entry_);
    c_nmus_vvvl_ = true;
  }
  return nmus_vvvl_;
}

int  const & baby_base::nph() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nph_ && b_nph_){
    b_nph_->GetEntry(entry_);
    c_nph_ = true;
  }
  return nph_;
}

int  const & baby_base::npv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_npv_ && b_npv_){
    b_npv_->GetEntry(entry_);
    c_npv_ = true;
  }
  return npv_;
}

int  const & baby_base::ntks() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntks_ && b_ntks_){
    b_ntks_->GetEntry(entry_);
    c_ntks_ = true;
  }
  return ntks_;
}

int  const & baby_base::ntruels() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntruels_ && b_ntruels_){
    b_ntruels_->GetEntry(entry_);
    c_ntruels_ = true;
  }
  return ntruels_;
}

int  const & baby_base::ntruleps() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntruleps_ && b_ntruleps_){
    b_ntruleps_->GetEntry(entry_);
    c_ntruleps_ = true;
  }
  return ntruleps_;
}

int  const & baby_base::ntrumus() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntrumus_ && b_ntrumus_){
    b_ntrumus_->GetEntry(entry_);
    c_ntrumus_ = true;
  }
  return ntrumus_;
}

int  const & baby_base::ntrupv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntrupv_ && b_ntrupv_){
    b_ntrupv_->GetEntry(entry_);
    c_ntrupv_ = true;
  }
  return ntrupv_;
}

int  const & baby_base::ntrutaush() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntrutaush_ && b_ntrutaush_){
    b_ntrutaush_->GetEntry(entry_);
    c_ntrutaush_ = true;
  }
  return ntrutaush_;
}

int  const & baby_base::ntrutausl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ntrutausl_ && b_ntrutausl_){
    b_ntrutausl_->GetEntry(entry_);
    c_ntrutausl_ = true;
  }
  return ntrutausl_;
}

int  const & baby_base::nvels() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nvels_ && b_nvels_){
    b_nvels_->GetEntry(entry_);
    c_nvels_ = true;
  }
  return nvels_;
}

int  const & baby_base::nveto() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nveto_ && b_nveto_){
    b_nveto_->GetEntry(entry_);
    c_nveto_ = true;
  }
  return nveto_;
}

int  const & baby_base::nvleps() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nvleps_ && b_nvleps_){
    b_nvleps_->GetEntry(entry_);
    c_nvleps_ = true;
  }
  return nvleps_;
}

int  const & baby_base::nvmus() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nvmus_ && b_nvmus_){
    b_nvmus_->GetEntry(entry_);
    c_nvmus_ = true;
  }
  return nvmus_;
}

int  const & baby_base::run() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_run_ && b_run_){
    b_run_->GetEntry(entry_);
    c_run_ = true;
  }
  return run_;
}

std::vector<bool>  const & baby_base::els_ele105() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_ele105_ && b_els_ele105_){
    b_els_ele105_->GetEntry(entry_);
    c_els_ele105_ = true;
  }
  return els_ele105_;
}

std::vector<bool>  const & baby_base::els_ele23() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_ele23_ && b_els_ele23_){
    b_els_ele23_->GetEntry(entry_);
    c_els_ele23_ = true;
  }
  return els_ele23_;
}

std::vector<bool>  const & baby_base::els_ele8() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_ele8_ && b_els_ele8_){
    b_els_ele8_->GetEntry(entry_);
    c_els_ele8_ = true;
  }
  return els_ele8_;
}

std::vector<bool>  const & baby_base::els_inz() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_inz_ && b_els_inz_){
    b_els_inz_->GetEntry(entry_);
    c_els_inz_ = true;
  }
  return els_inz_;
}

std::vector<bool>  const & baby_base::els_inzv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_inzv_ && b_els_inzv_){
    b_els_inzv_->GetEntry(entry_);
    c_els_inzv_ = true;
  }
  return els_inzv_;
}

std::vector<bool>  const & baby_base::els_ispf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_ispf_ && b_els_ispf_){
    b_els_ispf_->GetEntry(entry_);
    c_els_ispf_ = true;
  }
  return els_ispf_;
}

std::vector<bool>  const & baby_base::els_sig() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_sig_ && b_els_sig_){
    b_els_sig_->GetEntry(entry_);
    c_els_sig_ = true;
  }
  return els_sig_;
}

std::vector<bool>  const & baby_base::els_sigid() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_sigid_ && b_els_sigid_){
    b_els_sigid_->GetEntry(entry_);
    c_els_sigid_ = true;
  }
  return els_sigid_;
}

std::vector<bool>  const & baby_base::els_tight() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_tight_ && b_els_tight_){
    b_els_tight_->GetEntry(entry_);
    c_els_tight_ = true;
  }
  return els_tight_;
}

std::vector<bool>  const & baby_base::els_tm() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_tm_ && b_els_tm_){
    b_els_tm_->GetEntry(entry_);
    c_els_tm_ = true;
  }
  return els_tm_;
}

std::vector<bool>  const & baby_base::els_vvvl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_vvvl_ && b_els_vvvl_){
    b_els_vvvl_->GetEntry(entry_);
    c_els_vvvl_ = true;
  }
  return els_vvvl_;
}

std::vector<bool>  const & baby_base::jets_islep() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_islep_ && b_jets_islep_){
    b_jets_islep_->GetEntry(entry_);
    c_jets_islep_ = true;
  }
  return jets_islep_;
}

std::vector<bool>  const & baby_base::mus_inz() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_inz_ && b_mus_inz_){
    b_mus_inz_->GetEntry(entry_);
    c_mus_inz_ = true;
  }
  return mus_inz_;
}

std::vector<bool>  const & baby_base::mus_inzv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_inzv_ && b_mus_inzv_){
    b_mus_inzv_->GetEntry(entry_);
    c_mus_inzv_ = true;
  }
  return mus_inzv_;
}

std::vector<bool>  const & baby_base::mus_isomu18() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_isomu18_ && b_mus_isomu18_){
    b_mus_isomu18_->GetEntry(entry_);
    c_mus_isomu18_ = true;
  }
  return mus_isomu18_;
}

std::vector<bool>  const & baby_base::mus_mu50() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_mu50_ && b_mus_mu50_){
    b_mus_mu50_->GetEntry(entry_);
    c_mus_mu50_ = true;
  }
  return mus_mu50_;
}

std::vector<bool>  const & baby_base::mus_mu8() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_mu8_ && b_mus_mu8_){
    b_mus_mu8_->GetEntry(entry_);
    c_mus_mu8_ = true;
  }
  return mus_mu8_;
}

std::vector<bool>  const & baby_base::mus_sig() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_sig_ && b_mus_sig_){
    b_mus_sig_->GetEntry(entry_);
    c_mus_sig_ = true;
  }
  return mus_sig_;
}

std::vector<bool>  const & baby_base::mus_sigid() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_sigid_ && b_mus_sigid_){
    b_mus_sigid_->GetEntry(entry_);
    c_mus_sigid_ = true;
  }
  return mus_sigid_;
}

std::vector<bool>  const & baby_base::mus_tight() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_tight_ && b_mus_tight_){
    b_mus_tight_->GetEntry(entry_);
    c_mus_tight_ = true;
  }
  return mus_tight_;
}

std::vector<bool>  const & baby_base::mus_tm() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_tm_ && b_mus_tm_){
    b_mus_tm_->GetEntry(entry_);
    c_mus_tm_ = true;
  }
  return mus_tm_;
}

std::vector<bool>  const & baby_base::mus_vvvl() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_vvvl_ && b_mus_vvvl_){
    b_mus_vvvl_->GetEntry(entry_);
    c_mus_vvvl_ = true;
  }
  return mus_vvvl_;
}

std::vector<bool>  const & baby_base::ph_ph90() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ph_ph90_ && b_ph_ph90_){
    b_ph_ph90_->GetEntry(entry_);
    c_ph_ph90_ = true;
  }
  return ph_ph90_;
}

std::vector<bool>  const & baby_base::ph_tm() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ph_tm_ && b_ph_tm_){
    b_ph_tm_->GetEntry(entry_);
    c_ph_tm_ = true;
  }
  return ph_tm_;
}

std::vector<bool>  const & baby_base::sys_pass() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_pass_ && b_sys_pass_){
    b_sys_pass_->GetEntry(entry_);
    c_sys_pass_ = true;
  }
  return sys_pass_;
}

std::vector<bool>  const & baby_base::trig() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_trig_ && b_trig_){
    b_trig_->GetEntry(entry_);
    c_trig_ = true;
  }
  return trig_;
}

std::vector<float>  const & baby_base::els_d0() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_d0_ && b_els_d0_){
    b_els_d0_->GetEntry(entry_);
    c_els_d0_ = true;
  }
  return els_d0_;
}

std::vector<float>  const & baby_base::els_dz() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_dz_ && b_els_dz_){
    b_els_dz_->GetEntry(entry_);
    c_els_dz_ = true;
  }
  return els_dz_;
}

std::vector<float>  const & baby_base::els_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_eta_ && b_els_eta_){
    b_els_eta_->GetEntry(entry_);
    c_els_eta_ = true;
  }
  return els_eta_;
}

std::vector<float>  const & baby_base::els_miniso() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_miniso_ && b_els_miniso_){
    b_els_miniso_->GetEntry(entry_);
    c_els_miniso_ = true;
  }
  return els_miniso_;
}

std::vector<float>  const & baby_base::els_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_phi_ && b_els_phi_){
    b_els_phi_->GetEntry(entry_);
    c_els_phi_ = true;
  }
  return els_phi_;
}

std::vector<float>  const & baby_base::els_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_pt_ && b_els_pt_){
    b_els_pt_->GetEntry(entry_);
    c_els_pt_ = true;
  }
  return els_pt_;
}

std::vector<float>  const & baby_base::els_reliso() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_reliso_ && b_els_reliso_){
    b_els_reliso_->GetEntry(entry_);
    c_els_reliso_ = true;
  }
  return els_reliso_;
}

std::vector<float>  const & baby_base::els_sceta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_sceta_ && b_els_sceta_){
    b_els_sceta_->GetEntry(entry_);
    c_els_sceta_ = true;
  }
  return els_sceta_;
}

std::vector<float>  const & baby_base::els_vvvl_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_vvvl_eta_ && b_els_vvvl_eta_){
    b_els_vvvl_eta_->GetEntry(entry_);
    c_els_vvvl_eta_ = true;
  }
  return els_vvvl_eta_;
}

std::vector<float>  const & baby_base::els_vvvl_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_vvvl_phi_ && b_els_vvvl_phi_){
    b_els_vvvl_phi_->GetEntry(entry_);
    c_els_vvvl_phi_ = true;
  }
  return els_vvvl_phi_;
}

std::vector<float>  const & baby_base::els_vvvl_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_vvvl_pt_ && b_els_vvvl_pt_){
    b_els_vvvl_pt_->GetEntry(entry_);
    c_els_vvvl_pt_ = true;
  }
  return els_vvvl_pt_;
}

std::vector<float>  const & baby_base::fjets08_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_eta_ && b_fjets08_eta_){
    b_fjets08_eta_->GetEntry(entry_);
    c_fjets08_eta_ = true;
  }
  return fjets08_eta_;
}

std::vector<float>  const & baby_base::fjets08_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_m_ && b_fjets08_m_){
    b_fjets08_m_->GetEntry(entry_);
    c_fjets08_m_ = true;
  }
  return fjets08_m_;
}

std::vector<float>  const & baby_base::fjets08_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_phi_ && b_fjets08_phi_){
    b_fjets08_phi_->GetEntry(entry_);
    c_fjets08_phi_ = true;
  }
  return fjets08_phi_;
}

std::vector<float>  const & baby_base::fjets08_poscsv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_poscsv_ && b_fjets08_poscsv_){
    b_fjets08_poscsv_->GetEntry(entry_);
    c_fjets08_poscsv_ = true;
  }
  return fjets08_poscsv_;
}

std::vector<float>  const & baby_base::fjets08_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_pt_ && b_fjets08_pt_){
    b_fjets08_pt_->GetEntry(entry_);
    c_fjets08_pt_ = true;
  }
  return fjets08_pt_;
}

std::vector<float>  const & baby_base::fjets08_sumcsv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_sumcsv_ && b_fjets08_sumcsv_){
    b_fjets08_sumcsv_->GetEntry(entry_);
    c_fjets08_sumcsv_ = true;
  }
  return fjets08_sumcsv_;
}

std::vector<float>  const & baby_base::fjets_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_eta_ && b_fjets_eta_){
    b_fjets_eta_->GetEntry(entry_);
    c_fjets_eta_ = true;
  }
  return fjets_eta_;
}

std::vector<float>  const & baby_base::fjets_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_m_ && b_fjets_m_){
    b_fjets_m_->GetEntry(entry_);
    c_fjets_m_ = true;
  }
  return fjets_m_;
}

std::vector<float>  const & baby_base::fjets_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_phi_ && b_fjets_phi_){
    b_fjets_phi_->GetEntry(entry_);
    c_fjets_phi_ = true;
  }
  return fjets_phi_;
}

std::vector<float>  const & baby_base::fjets_poscsv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_poscsv_ && b_fjets_poscsv_){
    b_fjets_poscsv_->GetEntry(entry_);
    c_fjets_poscsv_ = true;
  }
  return fjets_poscsv_;
}

std::vector<float>  const & baby_base::fjets_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_pt_ && b_fjets_pt_){
    b_fjets_pt_->GetEntry(entry_);
    c_fjets_pt_ = true;
  }
  return fjets_pt_;
}

std::vector<float>  const & baby_base::fjets_sumcsv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_sumcsv_ && b_fjets_sumcsv_){
    b_fjets_sumcsv_->GetEntry(entry_);
    c_fjets_sumcsv_ = true;
  }
  return fjets_sumcsv_;
}

std::vector<float>  const & baby_base::jets_csv() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_csv_ && b_jets_csv_){
    b_jets_csv_->GetEntry(entry_);
    c_jets_csv_ = true;
  }
  return jets_csv_;
}

std::vector<float>  const & baby_base::jets_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_eta_ && b_jets_eta_){
    b_jets_eta_->GetEntry(entry_);
    c_jets_eta_ = true;
  }
  return jets_eta_;
}

std::vector<float>  const & baby_base::jets_m() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_m_ && b_jets_m_){
    b_jets_m_->GetEntry(entry_);
    c_jets_m_ = true;
  }
  return jets_m_;
}

std::vector<float>  const & baby_base::jets_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_phi_ && b_jets_phi_){
    b_jets_phi_->GetEntry(entry_);
    c_jets_phi_ = true;
  }
  return jets_phi_;
}

std::vector<float>  const & baby_base::jets_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_pt_ && b_jets_pt_){
    b_jets_pt_->GetEntry(entry_);
    c_jets_pt_ = true;
  }
  return jets_pt_;
}

std::vector<float>  const & baby_base::jets_pt_res() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_pt_res_ && b_jets_pt_res_){
    b_jets_pt_res_->GetEntry(entry_);
    c_jets_pt_res_ = true;
  }
  return jets_pt_res_;
}

std::vector<float>  const & baby_base::leps_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_leps_eta_ && b_leps_eta_){
    b_leps_eta_->GetEntry(entry_);
    c_leps_eta_ = true;
  }
  return leps_eta_;
}

std::vector<float>  const & baby_base::leps_id() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_leps_id_ && b_leps_id_){
    b_leps_id_->GetEntry(entry_);
    c_leps_id_ = true;
  }
  return leps_id_;
}

std::vector<float>  const & baby_base::leps_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_leps_phi_ && b_leps_phi_){
    b_leps_phi_->GetEntry(entry_);
    c_leps_phi_ = true;
  }
  return leps_phi_;
}

std::vector<float>  const & baby_base::leps_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_leps_pt_ && b_leps_pt_){
    b_leps_pt_->GetEntry(entry_);
    c_leps_pt_ = true;
  }
  return leps_pt_;
}

std::vector<float>  const & baby_base::mc_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mc_eta_ && b_mc_eta_){
    b_mc_eta_->GetEntry(entry_);
    c_mc_eta_ = true;
  }
  return mc_eta_;
}

std::vector<float>  const & baby_base::mc_mass() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mc_mass_ && b_mc_mass_){
    b_mc_mass_->GetEntry(entry_);
    c_mc_mass_ = true;
  }
  return mc_mass_;
}

std::vector<float>  const & baby_base::mc_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mc_phi_ && b_mc_phi_){
    b_mc_phi_->GetEntry(entry_);
    c_mc_phi_ = true;
  }
  return mc_phi_;
}

std::vector<float>  const & baby_base::mc_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mc_pt_ && b_mc_pt_){
    b_mc_pt_->GetEntry(entry_);
    c_mc_pt_ = true;
  }
  return mc_pt_;
}

std::vector<float>  const & baby_base::mus_d0() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_d0_ && b_mus_d0_){
    b_mus_d0_->GetEntry(entry_);
    c_mus_d0_ = true;
  }
  return mus_d0_;
}

std::vector<float>  const & baby_base::mus_dz() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_dz_ && b_mus_dz_){
    b_mus_dz_->GetEntry(entry_);
    c_mus_dz_ = true;
  }
  return mus_dz_;
}

std::vector<float>  const & baby_base::mus_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_eta_ && b_mus_eta_){
    b_mus_eta_->GetEntry(entry_);
    c_mus_eta_ = true;
  }
  return mus_eta_;
}

std::vector<float>  const & baby_base::mus_miniso() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_miniso_ && b_mus_miniso_){
    b_mus_miniso_->GetEntry(entry_);
    c_mus_miniso_ = true;
  }
  return mus_miniso_;
}

std::vector<float>  const & baby_base::mus_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_phi_ && b_mus_phi_){
    b_mus_phi_->GetEntry(entry_);
    c_mus_phi_ = true;
  }
  return mus_phi_;
}

std::vector<float>  const & baby_base::mus_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_pt_ && b_mus_pt_){
    b_mus_pt_->GetEntry(entry_);
    c_mus_pt_ = true;
  }
  return mus_pt_;
}

std::vector<float>  const & baby_base::mus_reliso() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_reliso_ && b_mus_reliso_){
    b_mus_reliso_->GetEntry(entry_);
    c_mus_reliso_ = true;
  }
  return mus_reliso_;
}

std::vector<float>  const & baby_base::mus_vvvl_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_vvvl_eta_ && b_mus_vvvl_eta_){
    b_mus_vvvl_eta_->GetEntry(entry_);
    c_mus_vvvl_eta_ = true;
  }
  return mus_vvvl_eta_;
}

std::vector<float>  const & baby_base::mus_vvvl_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_vvvl_phi_ && b_mus_vvvl_phi_){
    b_mus_vvvl_phi_->GetEntry(entry_);
    c_mus_vvvl_phi_ = true;
  }
  return mus_vvvl_phi_;
}

std::vector<float>  const & baby_base::mus_vvvl_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_vvvl_pt_ && b_mus_vvvl_pt_){
    b_mus_vvvl_pt_->GetEntry(entry_);
    c_mus_vvvl_pt_ = true;
  }
  return mus_vvvl_pt_;
}

std::vector<float>  const & baby_base::ph_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ph_eta_ && b_ph_eta_){
    b_ph_eta_->GetEntry(entry_);
    c_ph_eta_ = true;
  }
  return ph_eta_;
}

std::vector<float>  const & baby_base::ph_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ph_phi_ && b_ph_phi_){
    b_ph_phi_->GetEntry(entry_);
    c_ph_phi_ = true;
  }
  return ph_phi_;
}

std::vector<float>  const & baby_base::ph_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_ph_pt_ && b_ph_pt_){
    b_ph_pt_->GetEntry(entry_);
    c_ph_pt_ = true;
  }
  return ph_pt_;
}

std::vector<float>  const & baby_base::sys_bctag() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_bctag_ && b_sys_bctag_){
    b_sys_bctag_->GetEntry(entry_);
    c_sys_bctag_ = true;
  }
  return sys_bctag_;
}

std::vector<float>  const & baby_base::sys_fs_bctag() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_fs_bctag_ && b_sys_fs_bctag_){
    b_sys_fs_bctag_->GetEntry(entry_);
    c_sys_fs_bctag_ = true;
  }
  return sys_fs_bctag_;
}

std::vector<float>  const & baby_base::sys_fs_lep() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_fs_lep_ && b_sys_fs_lep_){
    b_sys_fs_lep_->GetEntry(entry_);
    c_sys_fs_lep_ = true;
  }
  return sys_fs_lep_;
}

std::vector<float>  const & baby_base::sys_fs_udsgtag() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_fs_udsgtag_ && b_sys_fs_udsgtag_){
    b_sys_fs_udsgtag_->GetEntry(entry_);
    c_sys_fs_udsgtag_ = true;
  }
  return sys_fs_udsgtag_;
}

std::vector<float>  const & baby_base::sys_ht() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_ht_ && b_sys_ht_){
    b_sys_ht_->GetEntry(entry_);
    c_sys_ht_ = true;
  }
  return sys_ht_;
}

std::vector<float>  const & baby_base::sys_isr() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_isr_ && b_sys_isr_){
    b_sys_isr_->GetEntry(entry_);
    c_sys_isr_ = true;
  }
  return sys_isr_;
}

std::vector<float>  const & baby_base::sys_lep() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_lep_ && b_sys_lep_){
    b_sys_lep_->GetEntry(entry_);
    c_sys_lep_ = true;
  }
  return sys_lep_;
}

std::vector<float>  const & baby_base::sys_met() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_met_ && b_sys_met_){
    b_sys_met_->GetEntry(entry_);
    c_sys_met_ = true;
  }
  return sys_met_;
}

std::vector<float>  const & baby_base::sys_mj() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_mj_ && b_sys_mj_){
    b_sys_mj_->GetEntry(entry_);
    c_sys_mj_ = true;
  }
  return sys_mj_;
}

std::vector<float>  const & baby_base::sys_mt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_mt_ && b_sys_mt_){
    b_sys_mt_->GetEntry(entry_);
    c_sys_mt_ = true;
  }
  return sys_mt_;
}

std::vector<float>  const & baby_base::sys_muf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_muf_ && b_sys_muf_){
    b_sys_muf_->GetEntry(entry_);
    c_sys_muf_ = true;
  }
  return sys_muf_;
}

std::vector<float>  const & baby_base::sys_mur() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_mur_ && b_sys_mur_){
    b_sys_mur_->GetEntry(entry_);
    c_sys_mur_ = true;
  }
  return sys_mur_;
}

std::vector<float>  const & baby_base::sys_murf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_murf_ && b_sys_murf_){
    b_sys_murf_->GetEntry(entry_);
    c_sys_murf_ = true;
  }
  return sys_murf_;
}

std::vector<float>  const & baby_base::sys_pdf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_pdf_ && b_sys_pdf_){
    b_sys_pdf_->GetEntry(entry_);
    c_sys_pdf_ = true;
  }
  return sys_pdf_;
}

std::vector<float>  const & baby_base::sys_trig() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_trig_ && b_sys_trig_){
    b_sys_trig_->GetEntry(entry_);
    c_sys_trig_ = true;
  }
  return sys_trig_;
}

std::vector<float>  const & baby_base::sys_udsgtag() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_udsgtag_ && b_sys_udsgtag_){
    b_sys_udsgtag_->GetEntry(entry_);
    c_sys_udsgtag_ = true;
  }
  return sys_udsgtag_;
}

std::vector<float>  const & baby_base::tks_eta() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_tks_eta_ && b_tks_eta_){
    b_tks_eta_->GetEntry(entry_);
    c_tks_eta_ = true;
  }
  return tks_eta_;
}

std::vector<float>  const & baby_base::tks_miniso() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_tks_miniso_ && b_tks_miniso_){
    b_tks_miniso_->GetEntry(entry_);
    c_tks_miniso_ = true;
  }
  return tks_miniso_;
}

std::vector<float>  const & baby_base::tks_mt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_tks_mt_ && b_tks_mt_){
    b_tks_mt_->GetEntry(entry_);
    c_tks_mt_ = true;
  }
  return tks_mt_;
}

std::vector<float>  const & baby_base::tks_mt2() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_tks_mt2_ && b_tks_mt2_){
    b_tks_mt2_->GetEntry(entry_);
    c_tks_mt2_ = true;
  }
  return tks_mt2_;
}

std::vector<float>  const & baby_base::tks_phi() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_tks_phi_ && b_tks_phi_){
    b_tks_phi_->GetEntry(entry_);
    c_tks_phi_ = true;
  }
  return tks_phi_;
}

std::vector<float>  const & baby_base::tks_pt() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_tks_pt_ && b_tks_pt_){
    b_tks_pt_->GetEntry(entry_);
    c_tks_pt_ = true;
  }
  return tks_pt_;
}

std::vector<float>  const & baby_base::trig_prescale() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_trig_prescale_ && b_trig_prescale_){
    b_trig_prescale_->GetEntry(entry_);
    c_trig_prescale_ = true;
  }
  return trig_prescale_;
}

std::vector<float>  const & baby_base::w_pdf() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_w_pdf_ && b_w_pdf_){
    b_w_pdf_->GetEntry(entry_);
    c_w_pdf_ = true;
  }
  return w_pdf_;
}

std::vector<int>  const & baby_base::els_charge() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_els_charge_ && b_els_charge_){
    b_els_charge_->GetEntry(entry_);
    c_els_charge_ = true;
  }
  return els_charge_;
}

std::vector<int>  const & baby_base::fjets08_btags() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_btags_ && b_fjets08_btags_){
    b_fjets08_btags_->GetEntry(entry_);
    c_fjets08_btags_ = true;
  }
  return fjets08_btags_;
}

std::vector<int>  const & baby_base::fjets08_nconst() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets08_nconst_ && b_fjets08_nconst_){
    b_fjets08_nconst_->GetEntry(entry_);
    c_fjets08_nconst_ = true;
  }
  return fjets08_nconst_;
}

std::vector<int>  const & baby_base::fjets_btags() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_btags_ && b_fjets_btags_){
    b_fjets_btags_->GetEntry(entry_);
    c_fjets_btags_ = true;
  }
  return fjets_btags_;
}

std::vector<int>  const & baby_base::fjets_nconst() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_fjets_nconst_ && b_fjets_nconst_){
    b_fjets_nconst_->GetEntry(entry_);
    c_fjets_nconst_ = true;
  }
  return fjets_nconst_;
}

std::vector<int>  const & baby_base::jets_fjet08_index() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_fjet08_index_ && b_jets_fjet08_index_){
    b_jets_fjet08_index_->GetEntry(entry_);
    c_jets_fjet08_index_ = true;
  }
  return jets_fjet08_index_;
}

std::vector<int>  const & baby_base::jets_fjet_index() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_fjet_index_ && b_jets_fjet_index_){
    b_jets_fjet_index_->GetEntry(entry_);
    c_jets_fjet_index_ = true;
  }
  return jets_fjet_index_;
}

std::vector<int>  const & baby_base::jets_hflavor() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_jets_hflavor_ && b_jets_hflavor_){
    b_jets_hflavor_->GetEntry(entry_);
    c_jets_hflavor_ = true;
  }
  return jets_hflavor_;
}

std::vector<int>  const & baby_base::mc_id() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mc_id_ && b_mc_id_){
    b_mc_id_->GetEntry(entry_);
    c_mc_id_ = true;
  }
  return mc_id_;
}

std::vector<int>  const & baby_base::mc_mom() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mc_mom_ && b_mc_mom_){
    b_mc_mom_->GetEntry(entry_);
    c_mc_mom_ = true;
  }
  return mc_mom_;
}

std::vector<int>  const & baby_base::mus_charge() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_mus_charge_ && b_mus_charge_){
    b_mus_charge_->GetEntry(entry_);
    c_mus_charge_ = true;
  }
  return mus_charge_;
}

std::vector<int>  const & baby_base::sys_nbm() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_nbm_ && b_sys_nbm_){
    b_sys_nbm_->GetEntry(entry_);
    c_sys_nbm_ = true;
  }
  return sys_nbm_;
}

std::vector<int>  const & baby_base::sys_njets() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_sys_njets_ && b_sys_njets_){
    b_sys_njets_->GetEntry(entry_);
    c_sys_njets_ = true;
  }
  return sys_njets_;
}

std::vector<int>  const & baby_base::tks_pdg() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_tks_pdg_ && b_tks_pdg_){
    b_tks_pdg_->GetEntry(entry_);
    c_tks_pdg_ = true;
  }
  return tks_pdg_;
}

Long64_t  & baby_base::event(){
  if(read_only_ && !c_event_ && b_event_){
    b_event_->GetEntry(entry_);
    c_event_ = true;
  }
  return event_;
}

bool  & baby_base::fromGS(){
  if(read_only_ && !c_fromGS_ && b_fromGS_){
    b_fromGS_->GetEntry(entry_);
    c_fromGS_ = true;
  }
  return fromGS_;
}

bool  & baby_base::jetmismeas(){
  if(read_only_ && !c_jetmismeas_ && b_jetmismeas_){
    b_jetmismeas_->GetEntry(entry_);
    c_jetmismeas_ = true;
  }
  return jetmismeas_;
}

bool  & baby_base::low_dphi(){
  if(read_only_ && !c_low_dphi_ && b_low_dphi_){
    b_low_dphi_->GetEntry(entry_);
    c_low_dphi_ = true;
  }
  return low_dphi_;
}

bool  & baby_base::nonblind(){
  if(read_only_ && !c_nonblind_ && b_nonblind_){
    b_nonblind_->GetEntry(entry_);
    c_nonblind_ = true;
  }
  return nonblind_;
}

bool  & baby_base::pass(){
  if(read_only_ && !c_pass_ && b_pass_){
    b_pass_->GetEntry(entry_);
    c_pass_ = true;
  }
  return pass_;
}

bool  & baby_base::pass_cschalo(){
  if(read_only_ && !c_pass_cschalo_ && b_pass_cschalo_){
    b_pass_cschalo_->GetEntry(entry_);
    c_pass_cschalo_ = true;
  }
  return pass_cschalo_;
}

bool  & baby_base::pass_eebadsc(){
  if(read_only_ && !c_pass_eebadsc_ && b_pass_eebadsc_){
    b_pass_eebadsc_->GetEntry(entry_);
    c_pass_eebadsc_ = true;
  }
  return pass_eebadsc_;
}

bool  & baby_base::pass_goodv(){
  if(read_only_ && !c_pass_goodv_ && b_pass_goodv_){
    b_pass_goodv_->GetEntry(entry_);
    c_pass_goodv_ = true;
  }
  return pass_goodv_;
}

bool  & baby_base::pass_hbhe(){
  if(read_only_ && !c_pass_hbhe_ && b_pass_hbhe_){
    b_pass_hbhe_->GetEntry(entry_);
    c_pass_hbhe_ = true;
  }
  return pass_hbhe_;
}

bool  & baby_base::pass_hbheiso(){
  if(read_only_ && !c_pass_hbheiso_ && b_pass_hbheiso_){
    b_pass_hbheiso_->GetEntry(entry_);
    c_pass_hbheiso_ = true;
  }
  return pass_hbheiso_;
}

bool  & baby_base::pass_jets(){
  if(read_only_ && !c_pass_jets_ && b_pass_jets_){
    b_pass_jets_->GetEntry(entry_);
    c_pass_jets_ = true;
  }
  return pass_jets_;
}

bool  & baby_base::pass_jets_nohf(){
  if(read_only_ && !c_pass_jets_nohf_ && b_pass_jets_nohf_){
    b_pass_jets_nohf_->GetEntry(entry_);
    c_pass_jets_nohf_ = true;
  }
  return pass_jets_nohf_;
}

bool  & baby_base::pass_jets_ra2(){
  if(read_only_ && !c_pass_jets_ra2_ && b_pass_jets_ra2_){
    b_pass_jets_ra2_->GetEntry(entry_);
    c_pass_jets_ra2_ = true;
  }
  return pass_jets_ra2_;
}

bool  & baby_base::pass_jets_tight(){
  if(read_only_ && !c_pass_jets_tight_ && b_pass_jets_tight_){
    b_pass_jets_tight_->GetEntry(entry_);
    c_pass_jets_tight_ = true;
  }
  return pass_jets_tight_;
}

bool  & baby_base::pass_jets_tight_ra2(){
  if(read_only_ && !c_pass_jets_tight_ra2_ && b_pass_jets_tight_ra2_){
    b_pass_jets_tight_ra2_->GetEntry(entry_);
    c_pass_jets_tight_ra2_ = true;
  }
  return pass_jets_tight_ra2_;
}

bool  & baby_base::pass_nohf(){
  if(read_only_ && !c_pass_nohf_ && b_pass_nohf_){
    b_pass_nohf_->GetEntry(entry_);
    c_pass_nohf_ = true;
  }
  return pass_nohf_;
}

bool  & baby_base::pass_ra2(){
  if(read_only_ && !c_pass_ra2_ && b_pass_ra2_){
    b_pass_ra2_->GetEntry(entry_);
    c_pass_ra2_ = true;
  }
  return pass_ra2_;
}

bool  & baby_base::stitch(){
  if(read_only_ && !c_stitch_ && b_stitch_){
    b_stitch_->GetEntry(entry_);
    c_stitch_ = true;
  }
  return stitch_;
}

float  & baby_base::dphi1(){
  if(read_only_ && !c_dphi1_ && b_dphi1_){
    b_dphi1_->GetEntry(entry_);
    c_dphi1_ = true;
  }
  return dphi1_;
}

float  & baby_base::dphi2(){
  if(read_only_ && !c_dphi2_ && b_dphi2_){
    b_dphi2_->GetEntry(entry_);
    c_dphi2_ = true;
  }
  return dphi2_;
}

float  & baby_base::dphi3(){
  if(read_only_ && !c_dphi3_ && b_dphi3_){
    b_dphi3_->GetEntry(entry_);
    c_dphi3_ = true;
  }
  return dphi3_;
}

float  & baby_base::dphi4(){
  if(read_only_ && !c_dphi4_ && b_dphi4_){
    b_dphi4_->GetEntry(entry_);
    c_dphi4_ = true;
  }
  return dphi4_;
}

float  & baby_base::dphi_wlep(){
  if(read_only_ && !c_dphi_wlep_ && b_dphi_wlep_){
    b_dphi_wlep_->GetEntry(entry_);
    c_dphi_wlep_ = true;
  }
  return dphi_wlep_;
}

float  & baby_base::eff_trig(){
  if(read_only_ && !c_eff_trig_ && b_eff_trig_){
    b_eff_trig_->GetEntry(entry_);
    c_eff_trig_ = true;
  }
  return eff_trig_;
}

float  & baby_base::elel_eta(){
  if(read_only_ && !c_elel_eta_ && b_elel_eta_){
    b_elel_eta_->GetEntry(entry_);
    c_elel_eta_ = true;
  }
  return elel_eta_;
}

float  & baby_base::elel_m(){
  if(read_only_ && !c_elel_m_ && b_elel_m_){
    b_elel_m_->GetEntry(entry_);
    c_elel_m_ = true;
  }
  return elel_m_;
}

float  & baby_base::elel_phi(){
  if(read_only_ && !c_elel_phi_ && b_elel_phi_){
    b_elel_phi_->GetEntry(entry_);
    c_elel_phi_ = true;
  }
  return elel_phi_;
}

float  & baby_base::elel_pt(){
  if(read_only_ && !c_elel_pt_ && b_elel_pt_){
    b_elel_pt_->GetEntry(entry_);
    c_elel_pt_ = true;
  }
  return elel_pt_;
}

float  & baby_base::elel_pt1(){
  if(read_only_ && !c_elel_pt1_ && b_elel_pt1_){
    b_elel_pt1_->GetEntry(entry_);
    c_elel_pt1_ = true;
  }
  return elel_pt1_;
}

float  & baby_base::elel_pt2(){
  if(read_only_ && !c_elel_pt2_ && b_elel_pt2_){
    b_elel_pt2_->GetEntry(entry_);
    c_elel_pt2_ = true;
  }
  return elel_pt2_;
}

float  & baby_base::elel_w(){
  if(read_only_ && !c_elel_w_ && b_elel_w_){
    b_elel_w_->GetEntry(entry_);
    c_elel_w_ = true;
  }
  return elel_w_;
}

float  & baby_base::elelv_eta(){
  if(read_only_ && !c_elelv_eta_ && b_elelv_eta_){
    b_elelv_eta_->GetEntry(entry_);
    c_elelv_eta_ = true;
  }
  return elelv_eta_;
}

float  & baby_base::elelv_m(){
  if(read_only_ && !c_elelv_m_ && b_elelv_m_){
    b_elelv_m_->GetEntry(entry_);
    c_elelv_m_ = true;
  }
  return elelv_m_;
}

float  & baby_base::elelv_phi(){
  if(read_only_ && !c_elelv_phi_ && b_elelv_phi_){
    b_elelv_phi_->GetEntry(entry_);
    c_elelv_phi_ = true;
  }
  return elelv_phi_;
}

float  & baby_base::elelv_pt(){
  if(read_only_ && !c_elelv_pt_ && b_elelv_pt_){
    b_elelv_pt_->GetEntry(entry_);
    c_elelv_pt_ = true;
  }
  return elelv_pt_;
}

float  & baby_base::elelv_pt1(){
  if(read_only_ && !c_elelv_pt1_ && b_elelv_pt1_){
    b_elelv_pt1_->GetEntry(entry_);
    c_elelv_pt1_ = true;
  }
  return elelv_pt1_;
}

float  & baby_base::elelv_pt2(){
  if(read_only_ && !c_elelv_pt2_ && b_elelv_pt2_){
    b_elelv_pt2_->GetEntry(entry_);
    c_elelv_pt2_ = true;
  }
  return elelv_pt2_;
}

float  & baby_base::elelv_w(){
  if(read_only_ && !c_elelv_w_ && b_elelv_w_){
    b_elelv_w_->GetEntry(entry_);
    c_elelv_w_ = true;
  }
  return elelv_w_;
}

float  & baby_base::elmu_eta(){
  if(read_only_ && !c_elmu_eta_ && b_elmu_eta_){
    b_elmu_eta_->GetEntry(entry_);
    c_elmu_eta_ = true;
  }
  return elmu_eta_;
}

float  & baby_base::elmu_m(){
  if(read_only_ && !c_elmu_m_ && b_elmu_m_){
    b_elmu_m_->GetEntry(entry_);
    c_elmu_m_ = true;
  }
  return elmu_m_;
}

float  & baby_base::elmu_phi(){
  if(read_only_ && !c_elmu_phi_ && b_elmu_phi_){
    b_elmu_phi_->GetEntry(entry_);
    c_elmu_phi_ = true;
  }
  return elmu_phi_;
}

float  & baby_base::elmu_pt(){
  if(read_only_ && !c_elmu_pt_ && b_elmu_pt_){
    b_elmu_pt_->GetEntry(entry_);
    c_elmu_pt_ = true;
  }
  return elmu_pt_;
}

float  & baby_base::elmu_pt1(){
  if(read_only_ && !c_elmu_pt1_ && b_elmu_pt1_){
    b_elmu_pt1_->GetEntry(entry_);
    c_elmu_pt1_ = true;
  }
  return elmu_pt1_;
}

float  & baby_base::elmu_pt2(){
  if(read_only_ && !c_elmu_pt2_ && b_elmu_pt2_){
    b_elmu_pt2_->GetEntry(entry_);
    c_elmu_pt2_ = true;
  }
  return elmu_pt2_;
}

float  & baby_base::elmu_w(){
  if(read_only_ && !c_elmu_w_ && b_elmu_w_){
    b_elmu_w_->GetEntry(entry_);
    c_elmu_w_ = true;
  }
  return elmu_w_;
}

float  & baby_base::ht(){
  if(read_only_ && !c_ht_ && b_ht_){
    b_ht_->GetEntry(entry_);
    c_ht_ = true;
  }
  return ht_;
}

float  & baby_base::ht_clean(){
  if(read_only_ && !c_ht_clean_ && b_ht_clean_){
    b_ht_clean_->GetEntry(entry_);
    c_ht_clean_ = true;
  }
  return ht_clean_;
}

float  & baby_base::ht_hlt(){
  if(read_only_ && !c_ht_hlt_ && b_ht_hlt_){
    b_ht_hlt_->GetEntry(entry_);
    c_ht_hlt_ = true;
  }
  return ht_hlt_;
}

float  & baby_base::ht_isr_me(){
  if(read_only_ && !c_ht_isr_me_ && b_ht_isr_me_){
    b_ht_isr_me_->GetEntry(entry_);
    c_ht_isr_me_ = true;
  }
  return ht_isr_me_;
}

float  & baby_base::ht_ra2(){
  if(read_only_ && !c_ht_ra2_ && b_ht_ra2_){
    b_ht_ra2_->GetEntry(entry_);
    c_ht_ra2_ = true;
  }
  return ht_ra2_;
}

float  & baby_base::ht_tru(){
  if(read_only_ && !c_ht_tru_ && b_ht_tru_){
    b_ht_tru_->GetEntry(entry_);
    c_ht_tru_ = true;
  }
  return ht_tru_;
}

float  & baby_base::isr_tru_eta(){
  if(read_only_ && !c_isr_tru_eta_ && b_isr_tru_eta_){
    b_isr_tru_eta_->GetEntry(entry_);
    c_isr_tru_eta_ = true;
  }
  return isr_tru_eta_;
}

float  & baby_base::isr_tru_phi(){
  if(read_only_ && !c_isr_tru_phi_ && b_isr_tru_phi_){
    b_isr_tru_phi_->GetEntry(entry_);
    c_isr_tru_phi_ = true;
  }
  return isr_tru_phi_;
}

float  & baby_base::isr_tru_pt(){
  if(read_only_ && !c_isr_tru_pt_ && b_isr_tru_pt_){
    b_isr_tru_pt_->GetEntry(entry_);
    c_isr_tru_pt_ = true;
  }
  return isr_tru_pt_;
}

float  & baby_base::jetsys_eta(){
  if(read_only_ && !c_jetsys_eta_ && b_jetsys_eta_){
    b_jetsys_eta_->GetEntry(entry_);
    c_jetsys_eta_ = true;
  }
  return jetsys_eta_;
}

float  & baby_base::jetsys_nob_eta(){
  if(read_only_ && !c_jetsys_nob_eta_ && b_jetsys_nob_eta_){
    b_jetsys_nob_eta_->GetEntry(entry_);
    c_jetsys_nob_eta_ = true;
  }
  return jetsys_nob_eta_;
}

float  & baby_base::jetsys_nob_phi(){
  if(read_only_ && !c_jetsys_nob_phi_ && b_jetsys_nob_phi_){
    b_jetsys_nob_phi_->GetEntry(entry_);
    c_jetsys_nob_phi_ = true;
  }
  return jetsys_nob_phi_;
}

float  & baby_base::jetsys_nob_pt(){
  if(read_only_ && !c_jetsys_nob_pt_ && b_jetsys_nob_pt_){
    b_jetsys_nob_pt_->GetEntry(entry_);
    c_jetsys_nob_pt_ = true;
  }
  return jetsys_nob_pt_;
}

float  & baby_base::jetsys_phi(){
  if(read_only_ && !c_jetsys_phi_ && b_jetsys_phi_){
    b_jetsys_phi_->GetEntry(entry_);
    c_jetsys_phi_ = true;
  }
  return jetsys_phi_;
}

float  & baby_base::jetsys_pt(){
  if(read_only_ && !c_jetsys_pt_ && b_jetsys_pt_){
    b_jetsys_pt_->GetEntry(entry_);
    c_jetsys_pt_ = true;
  }
  return jetsys_pt_;
}

float  & baby_base::met(){
  if(read_only_ && !c_met_ && b_met_){
    b_met_->GetEntry(entry_);
    c_met_ = true;
  }
  return met_;
}

float  & baby_base::met_calo(){
  if(read_only_ && !c_met_calo_ && b_met_calo_){
    b_met_calo_->GetEntry(entry_);
    c_met_calo_ = true;
  }
  return met_calo_;
}

float  & baby_base::met_calo_phi(){
  if(read_only_ && !c_met_calo_phi_ && b_met_calo_phi_){
    b_met_calo_phi_->GetEntry(entry_);
    c_met_calo_phi_ = true;
  }
  return met_calo_phi_;
}

float  & baby_base::met_mini(){
  if(read_only_ && !c_met_mini_ && b_met_mini_){
    b_met_mini_->GetEntry(entry_);
    c_met_mini_ = true;
  }
  return met_mini_;
}

float  & baby_base::met_mini_phi(){
  if(read_only_ && !c_met_mini_phi_ && b_met_mini_phi_){
    b_met_mini_phi_->GetEntry(entry_);
    c_met_mini_phi_ = true;
  }
  return met_mini_phi_;
}

float  & baby_base::met_nohf(){
  if(read_only_ && !c_met_nohf_ && b_met_nohf_){
    b_met_nohf_->GetEntry(entry_);
    c_met_nohf_ = true;
  }
  return met_nohf_;
}

float  & baby_base::met_nohf_phi(){
  if(read_only_ && !c_met_nohf_phi_ && b_met_nohf_phi_){
    b_met_nohf_phi_->GetEntry(entry_);
    c_met_nohf_phi_ = true;
  }
  return met_nohf_phi_;
}

float  & baby_base::met_phi(){
  if(read_only_ && !c_met_phi_ && b_met_phi_){
    b_met_phi_->GetEntry(entry_);
    c_met_phi_ = true;
  }
  return met_phi_;
}

float  & baby_base::met_raw(){
  if(read_only_ && !c_met_raw_ && b_met_raw_){
    b_met_raw_->GetEntry(entry_);
    c_met_raw_ = true;
  }
  return met_raw_;
}

float  & baby_base::met_raw_phi(){
  if(read_only_ && !c_met_raw_phi_ && b_met_raw_phi_){
    b_met_raw_phi_->GetEntry(entry_);
    c_met_raw_phi_ = true;
  }
  return met_raw_phi_;
}

float  & baby_base::met_rebal(){
  if(read_only_ && !c_met_rebal_ && b_met_rebal_){
    b_met_rebal_->GetEntry(entry_);
    c_met_rebal_ = true;
  }
  return met_rebal_;
}

float  & baby_base::met_tru(){
  if(read_only_ && !c_met_tru_ && b_met_tru_){
    b_met_tru_->GetEntry(entry_);
    c_met_tru_ = true;
  }
  return met_tru_;
}

float  & baby_base::met_tru_nuw(){
  if(read_only_ && !c_met_tru_nuw_ && b_met_tru_nuw_){
    b_met_tru_nuw_->GetEntry(entry_);
    c_met_tru_nuw_ = true;
  }
  return met_tru_nuw_;
}

float  & baby_base::met_tru_nuw_phi(){
  if(read_only_ && !c_met_tru_nuw_phi_ && b_met_tru_nuw_phi_){
    b_met_tru_nuw_phi_->GetEntry(entry_);
    c_met_tru_nuw_phi_ = true;
  }
  return met_tru_nuw_phi_;
}

float  & baby_base::met_tru_phi(){
  if(read_only_ && !c_met_tru_phi_ && b_met_tru_phi_){
    b_met_tru_phi_->GetEntry(entry_);
    c_met_tru_phi_ = true;
  }
  return met_tru_phi_;
}

float  & baby_base::mht(){
  if(read_only_ && !c_mht_ && b_mht_){
    b_mht_->GetEntry(entry_);
    c_mht_ = true;
  }
  return mht_;
}

float  & baby_base::mht_clean(){
  if(read_only_ && !c_mht_clean_ && b_mht_clean_){
    b_mht_clean_->GetEntry(entry_);
    c_mht_clean_ = true;
  }
  return mht_clean_;
}

float  & baby_base::mht_clean_phi(){
  if(read_only_ && !c_mht_clean_phi_ && b_mht_clean_phi_){
    b_mht_clean_phi_->GetEntry(entry_);
    c_mht_clean_phi_ = true;
  }
  return mht_clean_phi_;
}

float  & baby_base::mht_phi(){
  if(read_only_ && !c_mht_phi_ && b_mht_phi_){
    b_mht_phi_->GetEntry(entry_);
    c_mht_phi_ = true;
  }
  return mht_phi_;
}

float  & baby_base::mj(){
  if(read_only_ && !c_mj_ && b_mj_){
    b_mj_->GetEntry(entry_);
    c_mj_ = true;
  }
  return mj_;
}

float  & baby_base::mj08(){
  if(read_only_ && !c_mj08_ && b_mj08_){
    b_mj08_->GetEntry(entry_);
    c_mj08_ = true;
  }
  return mj08_;
}

float  & baby_base::mt(){
  if(read_only_ && !c_mt_ && b_mt_){
    b_mt_->GetEntry(entry_);
    c_mt_ = true;
  }
  return mt_;
}

float  & baby_base::mt_nohf(){
  if(read_only_ && !c_mt_nohf_ && b_mt_nohf_){
    b_mt_nohf_->GetEntry(entry_);
    c_mt_nohf_ = true;
  }
  return mt_nohf_;
}

float  & baby_base::mt_rebal(){
  if(read_only_ && !c_mt_rebal_ && b_mt_rebal_){
    b_mt_rebal_->GetEntry(entry_);
    c_mt_rebal_ = true;
  }
  return mt_rebal_;
}

float  & baby_base::mt_tru(){
  if(read_only_ && !c_mt_tru_ && b_mt_tru_){
    b_mt_tru_->GetEntry(entry_);
    c_mt_tru_ = true;
  }
  return mt_tru_;
}

float  & baby_base::mt_tru_nuw(){
  if(read_only_ && !c_mt_tru_nuw_ && b_mt_tru_nuw_){
    b_mt_tru_nuw_->GetEntry(entry_);
    c_mt_tru_nuw_ = true;
  }
  return mt_tru_nuw_;
}

float  & baby_base::mumu_eta(){
  if(read_only_ && !c_mumu_eta_ && b_mumu_eta_){
    b_mumu_eta_->GetEntry(entry_);
    c_mumu_eta_ = true;
  }
  return mumu_eta_;
}

float  & baby_base::mumu_m(){
  if(read_only_ && !c_mumu_m_ && b_mumu_m_){
    b_mumu_m_->GetEntry(entry_);
    c_mumu_m_ = true;
  }
  return mumu_m_;
}

float  & baby_base::mumu_phi(){
  if(read_only_ && !c_mumu_phi_ && b_mumu_phi_){
    b_mumu_phi_->GetEntry(entry_);
    c_mumu_phi_ = true;
  }
  return mumu_phi_;
}

float  & baby_base::mumu_pt(){
  if(read_only_ && !c_mumu_pt_ && b_mumu_pt_){
    b_mumu_pt_->GetEntry(entry_);
    c_mumu_pt_ = true;
  }
  return mumu_pt_;
}

float  & baby_base::mumu_pt1(){
  if(read_only_ && !c_mumu_pt1_ && b_mumu_pt1_){
    b_mumu_pt1_->GetEntry(entry_);
    c_mumu_pt1_ = true;
  }
  return mumu_pt1_;
}

float  & baby_base::mumu_pt2(){
  if(read_only_ && !c_mumu_pt2_ && b_mumu_pt2_){
    b_mumu_pt2_->GetEntry(entry_);
    c_mumu_pt2_ = true;
  }
  return mumu_pt2_;
}

float  & baby_base::mumu_w(){
  if(read_only_ && !c_mumu_w_ && b_mumu_w_){
    b_mumu_w_->GetEntry(entry_);
    c_mumu_w_ = true;
  }
  return mumu_w_;
}

float  & baby_base::mumuv_eta(){
  if(read_only_ && !c_mumuv_eta_ && b_mumuv_eta_){
    b_mumuv_eta_->GetEntry(entry_);
    c_mumuv_eta_ = true;
  }
  return mumuv_eta_;
}

float  & baby_base::mumuv_m(){
  if(read_only_ && !c_mumuv_m_ && b_mumuv_m_){
    b_mumuv_m_->GetEntry(entry_);
    c_mumuv_m_ = true;
  }
  return mumuv_m_;
}

float  & baby_base::mumuv_phi(){
  if(read_only_ && !c_mumuv_phi_ && b_mumuv_phi_){
    b_mumuv_phi_->GetEntry(entry_);
    c_mumuv_phi_ = true;
  }
  return mumuv_phi_;
}

float  & baby_base::mumuv_pt(){
  if(read_only_ && !c_mumuv_pt_ && b_mumuv_pt_){
    b_mumuv_pt_->GetEntry(entry_);
    c_mumuv_pt_ = true;
  }
  return mumuv_pt_;
}

float  & baby_base::mumuv_pt1(){
  if(read_only_ && !c_mumuv_pt1_ && b_mumuv_pt1_){
    b_mumuv_pt1_->GetEntry(entry_);
    c_mumuv_pt1_ = true;
  }
  return mumuv_pt1_;
}

float  & baby_base::mumuv_pt2(){
  if(read_only_ && !c_mumuv_pt2_ && b_mumuv_pt2_){
    b_mumuv_pt2_->GetEntry(entry_);
    c_mumuv_pt2_ = true;
  }
  return mumuv_pt2_;
}

float  & baby_base::mumuv_w(){
  if(read_only_ && !c_mumuv_w_ && b_mumuv_w_){
    b_mumuv_w_->GetEntry(entry_);
    c_mumuv_w_ = true;
  }
  return mumuv_w_;
}

float  & baby_base::ntrupv_mean(){
  if(read_only_ && !c_ntrupv_mean_ && b_ntrupv_mean_){
    b_ntrupv_mean_->GetEntry(entry_);
    c_ntrupv_mean_ = true;
  }
  return ntrupv_mean_;
}

float  & baby_base::onel_ele105(){
  if(read_only_ && !c_onel_ele105_ && b_onel_ele105_){
    b_onel_ele105_->GetEntry(entry_);
    c_onel_ele105_ = true;
  }
  return onel_ele105_;
}

float  & baby_base::onel_ele23(){
  if(read_only_ && !c_onel_ele23_ && b_onel_ele23_){
    b_onel_ele23_->GetEntry(entry_);
    c_onel_ele23_ = true;
  }
  return onel_ele23_;
}

float  & baby_base::onel_ele8(){
  if(read_only_ && !c_onel_ele8_ && b_onel_ele8_){
    b_onel_ele8_->GetEntry(entry_);
    c_onel_ele8_ = true;
  }
  return onel_ele8_;
}

float  & baby_base::onel_vvvl(){
  if(read_only_ && !c_onel_vvvl_ && b_onel_vvvl_){
    b_onel_vvvl_->GetEntry(entry_);
    c_onel_vvvl_ = true;
  }
  return onel_vvvl_;
}

float  & baby_base::onht(){
  if(read_only_ && !c_onht_ && b_onht_){
    b_onht_->GetEntry(entry_);
    c_onht_ = true;
  }
  return onht_;
}

float  & baby_base::onmet(){
  if(read_only_ && !c_onmet_ && b_onmet_){
    b_onmet_->GetEntry(entry_);
    c_onmet_ = true;
  }
  return onmet_;
}

float  & baby_base::onmu_isomu18(){
  if(read_only_ && !c_onmu_isomu18_ && b_onmu_isomu18_){
    b_onmu_isomu18_->GetEntry(entry_);
    c_onmu_isomu18_ = true;
  }
  return onmu_isomu18_;
}

float  & baby_base::onmu_mu50(){
  if(read_only_ && !c_onmu_mu50_ && b_onmu_mu50_){
    b_onmu_mu50_->GetEntry(entry_);
    c_onmu_mu50_ = true;
  }
  return onmu_mu50_;
}

float  & baby_base::onmu_mu8(){
  if(read_only_ && !c_onmu_mu8_ && b_onmu_mu8_){
    b_onmu_mu8_->GetEntry(entry_);
    c_onmu_mu8_ = true;
  }
  return onmu_mu8_;
}

float  & baby_base::onmu_vvvl(){
  if(read_only_ && !c_onmu_vvvl_ && b_onmu_vvvl_){
    b_onmu_vvvl_->GetEntry(entry_);
    c_onmu_vvvl_ = true;
  }
  return onmu_vvvl_;
}

float  & baby_base::onph_ph90(){
  if(read_only_ && !c_onph_ph90_ && b_onph_ph90_){
    b_onph_ph90_->GetEntry(entry_);
    c_onph_ph90_ = true;
  }
  return onph_ph90_;
}

float  & baby_base::w_btag(){
  if(read_only_ && !c_w_btag_ && b_w_btag_){
    b_w_btag_->GetEntry(entry_);
    c_w_btag_ = true;
  }
  return w_btag_;
}

float  & baby_base::w_fs_lep(){
  if(read_only_ && !c_w_fs_lep_ && b_w_fs_lep_){
    b_w_fs_lep_->GetEntry(entry_);
    c_w_fs_lep_ = true;
  }
  return w_fs_lep_;
}

float  & baby_base::w_lep(){
  if(read_only_ && !c_w_lep_ && b_w_lep_){
    b_w_lep_->GetEntry(entry_);
    c_w_lep_ = true;
  }
  return w_lep_;
}

float  & baby_base::w_lumi(){
  if(read_only_ && !c_w_lumi_ && b_w_lumi_){
    b_w_lumi_->GetEntry(entry_);
    c_w_lumi_ = true;
  }
  return w_lumi_;
}

float  & baby_base::w_pu(){
  if(read_only_ && !c_w_pu_ && b_w_pu_){
    b_w_pu_->GetEntry(entry_);
    c_w_pu_ = true;
  }
  return w_pu_;
}

float  & baby_base::w_toppt(){
  if(read_only_ && !c_w_toppt_ && b_w_toppt_){
    b_w_toppt_->GetEntry(entry_);
    c_w_toppt_ = true;
  }
  return w_toppt_;
}

float  & baby_base::weight(){
  if(read_only_ && !c_weight_ && b_weight_){
    b_weight_->GetEntry(entry_);
    c_weight_ = true;
  }
  return weight_;
}

int  & baby_base::lumiblock(){
  if(read_only_ && !c_lumiblock_ && b_lumiblock_){
    b_lumiblock_->GetEntry(entry_);
    c_lumiblock_ = true;
  }
  return lumiblock_;
}

int  & baby_base::mgluino(){
  if(read_only_ && !c_mgluino_ && b_mgluino_){
    b_mgluino_->GetEntry(entry_);
    c_mgluino_ = true;
  }
  return mgluino_;
}

int  & baby_base::mlsp(){
  if(read_only_ && !c_mlsp_ && b_mlsp_){
    b_mlsp_->GetEntry(entry_);
    c_mlsp_ = true;
  }
  return mlsp_;
}

int  & baby_base::nbl(){
  if(read_only_ && !c_nbl_ && b_nbl_){
    b_nbl_->GetEntry(entry_);
    c_nbl_ = true;
  }
  return nbl_;
}

int  & baby_base::nbm(){
  if(read_only_ && !c_nbm_ && b_nbm_){
    b_nbm_->GetEntry(entry_);
    c_nbm_ = true;
  }
  return nbm_;
}

int  & baby_base::nbm_ra2(){
  if(read_only_ && !c_nbm_ra2_ && b_nbm_ra2_){
    b_nbm_ra2_->GetEntry(entry_);
    c_nbm_ra2_ = true;
  }
  return nbm_ra2_;
}

int  & baby_base::nbt(){
  if(read_only_ && !c_nbt_ && b_nbt_){
    b_nbt_->GetEntry(entry_);
    c_nbt_ = true;
  }
  return nbt_;
}

int  & baby_base::nels(){
  if(read_only_ && !c_nels_ && b_nels_){
    b_nels_->GetEntry(entry_);
    c_nels_ = true;
  }
  return nels_;
}

int  & baby_base::nels_ele23(){
  if(read_only_ && !c_nels_ele23_ && b_nels_ele23_){
    b_nels_ele23_->GetEntry(entry_);
    c_nels_ele23_ = true;
  }
  return nels_ele23_;
}

int  & baby_base::nels_vvvl(){
  if(read_only_ && !c_nels_vvvl_ && b_nels_vvvl_){
    b_nels_vvvl_->GetEntry(entry_);
    c_nels_vvvl_ = true;
  }
  return nels_vvvl_;
}

int  & baby_base::nfjets(){
  if(read_only_ && !c_nfjets_ && b_nfjets_){
    b_nfjets_->GetEntry(entry_);
    c_nfjets_ = true;
  }
  return nfjets_;
}

int  & baby_base::nfjets08(){
  if(read_only_ && !c_nfjets08_ && b_nfjets08_){
    b_nfjets08_->GetEntry(entry_);
    c_nfjets08_ = true;
  }
  return nfjets08_;
}

int  & baby_base::nisr_me(){
  if(read_only_ && !c_nisr_me_ && b_nisr_me_){
    b_nisr_me_->GetEntry(entry_);
    c_nisr_me_ = true;
  }
  return nisr_me_;
}

int  & baby_base::njets(){
  if(read_only_ && !c_njets_ && b_njets_){
    b_njets_->GetEntry(entry_);
    c_njets_ = true;
  }
  return njets_;
}

int  & baby_base::njets_clean(){
  if(read_only_ && !c_njets_clean_ && b_njets_clean_){
    b_njets_clean_->GetEntry(entry_);
    c_njets_clean_ = true;
  }
  return njets_clean_;
}

int  & baby_base::njets_ra2(){
  if(read_only_ && !c_njets_ra2_ && b_njets_ra2_){
    b_njets_ra2_->GetEntry(entry_);
    c_njets_ra2_ = true;
  }
  return njets_ra2_;
}

int  & baby_base::nleps(){
  if(read_only_ && !c_nleps_ && b_nleps_){
    b_nleps_->GetEntry(entry_);
    c_nleps_ = true;
  }
  return nleps_;
}

int  & baby_base::nleps_tm(){
  if(read_only_ && !c_nleps_tm_ && b_nleps_tm_){
    b_nleps_tm_->GetEntry(entry_);
    c_nleps_tm_ = true;
  }
  return nleps_tm_;
}

int  & baby_base::nmus(){
  if(read_only_ && !c_nmus_ && b_nmus_){
    b_nmus_->GetEntry(entry_);
    c_nmus_ = true;
  }
  return nmus_;
}

int  & baby_base::nmus_isomu18(){
  if(read_only_ && !c_nmus_isomu18_ && b_nmus_isomu18_){
    b_nmus_isomu18_->GetEntry(entry_);
    c_nmus_isomu18_ = true;
  }
  return nmus_isomu18_;
}

int  & baby_base::nmus_vvvl(){
  if(read_only_ && !c_nmus_vvvl_ && b_nmus_vvvl_){
    b_nmus_vvvl_->GetEntry(entry_);
    c_nmus_vvvl_ = true;
  }
  return nmus_vvvl_;
}

int  & baby_base::nph(){
  if(read_only_ && !c_nph_ && b_nph_){
    b_nph_->GetEntry(entry_);
    c_nph_ = true;
  }
  return nph_;
}

int  & baby_base::npv(){
  if(read_only_ && !c_npv_ && b_npv_){
    b_npv_->GetEntry(entry_);
    c_npv_ = true;
  }
  return npv_;
}

int  & baby_base::ntks(){
  if(read_only_ && !c_ntks_ && b_ntks_){
    b_ntks_->GetEntry(entry_);
    c_ntks_ = true;
  }
  return ntks_;
}

int  & baby_base::ntruels(){
  if(read_only_ && !c_ntruels_ && b_ntruels_){
    b_ntruels_->GetEntry(entry_);
    c_ntruels_ = true;
  }
  return ntruels_;
}

int  & baby_base::ntruleps(){
  if(read_only_ && !c_ntruleps_ && b_ntruleps_){
    b_ntruleps_->GetEntry(entry_);
    c_ntruleps_ = true;
  }
  return ntruleps_;
}

int  & baby_base::ntrumus(){
  if(read_only_ && !c_ntrumus_ && b_ntrumus_){
    b_ntrumus_->GetEntry(entry_);
    c_ntrumus_ = true;
  }
  return ntrumus_;
}

int  & baby_base::ntrupv(){
  if(read_only_ && !c_ntrupv_ && b_ntrupv_){
    b_ntrupv_->GetEntry(entry_);
    c_ntrupv_ = true;
  }
  return ntrupv_;
}

int  & baby_base::ntrutaush(){
  if(read_only_ && !c_ntrutaush_ && b_ntrutaush_){
    b_ntrutaush_->GetEntry(entry_);
    c_ntrutaush_ = true;
  }
  return ntrutaush_;
}

int  & baby_base::ntrutausl(){
  if(read_only_ && !c_ntrutausl_ && b_ntrutausl_){
    b_ntrutausl_->GetEntry(entry_);
    c_ntrutausl_ = true;
  }
  return ntrutausl_;
}

int  & baby_base::nvels(){
  if(read_only_ && !c_nvels_ && b_nvels_){
    b_nvels_->GetEntry(entry_);
    c_nvels_ = true;
  }
  return nvels_;
}

int  & baby_base::nveto(){
  if(read_only_ && !c_nveto_ && b_nveto_){
    b_nveto_->GetEntry(entry_);
    c_nveto_ = true;
  }
  return nveto_;
}

int  & baby_base::nvleps(){
  if(read_only_ && !c_nvleps_ && b_nvleps_){
    b_nvleps_->GetEntry(entry_);
    c_nvleps_ = true;
  }
  return nvleps_;
}

int  & baby_base::nvmus(){
  if(read_only_ && !c_nvmus_ && b_nvmus_){
    b_nvmus_->GetEntry(entry_);
    c_nvmus_ = true;
  }
  return nvmus_;
}

int  & baby_base::run(){
  if(read_only_ && !c_run_ && b_run_){
    b_run_->GetEntry(entry_);
    c_run_ = true;
  }
  return run_;
}

std::vector<bool>  & baby_base::els_ele105(){
  if(read_only_ && !c_els_ele105_ && b_els_ele105_){
    b_els_ele105_->GetEntry(entry_);
    c_els_ele105_ = true;
  }
  return els_ele105_;
}

std::vector<bool>  & baby_base::els_ele23(){
  if(read_only_ && !c_els_ele23_ && b_els_ele23_){
    b_els_ele23_->GetEntry(entry_);
    c_els_ele23_ = true;
  }
  return els_ele23_;
}

std::vector<bool>  & baby_base::els_ele8(){
  if(read_only_ && !c_els_ele8_ && b_els_ele8_){
    b_els_ele8_->GetEntry(entry_);
    c_els_ele8_ = true;
  }
  return els_ele8_;
}

std::vector<bool>  & baby_base::els_inz(){
  if(read_only_ && !c_els_inz_ && b_els_inz_){
    b_els_inz_->GetEntry(entry_);
    c_els_inz_ = true;
  }
  return els_inz_;
}

std::vector<bool>  & baby_base::els_inzv(){
  if(read_only_ && !c_els_inzv_ && b_els_inzv_){
    b_els_inzv_->GetEntry(entry_);
    c_els_inzv_ = true;
  }
  return els_inzv_;
}

std::vector<bool>  & baby_base::els_ispf(){
  if(read_only_ && !c_els_ispf_ && b_els_ispf_){
    b_els_ispf_->GetEntry(entry_);
    c_els_ispf_ = true;
  }
  return els_ispf_;
}

std::vector<bool>  & baby_base::els_sig(){
  if(read_only_ && !c_els_sig_ && b_els_sig_){
    b_els_sig_->GetEntry(entry_);
    c_els_sig_ = true;
  }
  return els_sig_;
}

std::vector<bool>  & baby_base::els_sigid(){
  if(read_only_ && !c_els_sigid_ && b_els_sigid_){
    b_els_sigid_->GetEntry(entry_);
    c_els_sigid_ = true;
  }
  return els_sigid_;
}

std::vector<bool>  & baby_base::els_tight(){
  if(read_only_ && !c_els_tight_ && b_els_tight_){
    b_els_tight_->GetEntry(entry_);
    c_els_tight_ = true;
  }
  return els_tight_;
}

std::vector<bool>  & baby_base::els_tm(){
  if(read_only_ && !c_els_tm_ && b_els_tm_){
    b_els_tm_->GetEntry(entry_);
    c_els_tm_ = true;
  }
  return els_tm_;
}

std::vector<bool>  & baby_base::els_vvvl(){
  if(read_only_ && !c_els_vvvl_ && b_els_vvvl_){
    b_els_vvvl_->GetEntry(entry_);
    c_els_vvvl_ = true;
  }
  return els_vvvl_;
}

std::vector<bool>  & baby_base::jets_islep(){
  if(read_only_ && !c_jets_islep_ && b_jets_islep_){
    b_jets_islep_->GetEntry(entry_);
    c_jets_islep_ = true;
  }
  return jets_islep_;
}

std::vector<bool>  & baby_base::mus_inz(){
  if(read_only_ && !c_mus_inz_ && b_mus_inz_){
    b_mus_inz_->GetEntry(entry_);
    c_mus_inz_ = true;
  }
  return mus_inz_;
}

std::vector<bool>  & baby_base::mus_inzv(){
  if(read_only_ && !c_mus_inzv_ && b_mus_inzv_){
    b_mus_inzv_->GetEntry(entry_);
    c_mus_inzv_ = true;
  }
  return mus_inzv_;
}

std::vector<bool>  & baby_base::mus_isomu18(){
  if(read_only_ && !c_mus_isomu18_ && b_mus_isomu18_){
    b_mus_isomu18_->GetEntry(entry_);
    c_mus_isomu18_ = true;
  }
  return mus_isomu18_;
}

std::vector<bool>  & baby_base::mus_mu50(){
  if(read_only_ && !c_mus_mu50_ && b_mus_mu50_){
    b_mus_mu50_->GetEntry(entry_);
    c_mus_mu50_ = true;
  }
  return mus_mu50_;
}

std::vector<bool>  & baby_base::mus_mu8(){
  if(read_only_ && !c_mus_mu8_ && b_mus_mu8_){
    b_mus_mu8_->GetEntry(entry_);
    c_mus_mu8_ = true;
  }
  return mus_mu8_;
}

std::vector<bool>  & baby_base::mus_sig(){
  if(read_only_ && !c_mus_sig_ && b_mus_sig_){
    b_mus_sig_->GetEntry(entry_);
    c_mus_sig_ = true;
  }
  return mus_sig_;
}

std::vector<bool>  & baby_base::mus_sigid(){
  if(read_only_ && !c_mus_sigid_ && b_mus_sigid_){
    b_mus_sigid_->GetEntry(entry_);
    c_mus_sigid_ = true;
  }
  return mus_sigid_;
}

std::vector<bool>  & baby_base::mus_tight(){
  if(read_only_ && !c_mus_tight_ && b_mus_tight_){
    b_mus_tight_->GetEntry(entry_);
    c_mus_tight_ = true;
  }
  return mus_tight_;
}

std::vector<bool>  & baby_base::mus_tm(){
  if(read_only_ && !c_mus_tm_ && b_mus_tm_){
    b_mus_tm_->GetEntry(entry_);
    c_mus_tm_ = true;
  }
  return mus_tm_;
}

std::vector<bool>  & baby_base::mus_vvvl(){
  if(read_only_ && !c_mus_vvvl_ && b_mus_vvvl_){
    b_mus_vvvl_->GetEntry(entry_);
    c_mus_vvvl_ = true;
  }
  return mus_vvvl_;
}

std::vector<bool>  & baby_base::ph_ph90(){
  if(read_only_ && !c_ph_ph90_ && b_ph_ph90_){
    b_ph_ph90_->GetEntry(entry_);
    c_ph_ph90_ = true;
  }
  return ph_ph90_;
}

std::vector<bool>  & baby_base::ph_tm(){
  if(read_only_ && !c_ph_tm_ && b_ph_tm_){
    b_ph_tm_->GetEntry(entry_);
    c_ph_tm_ = true;
  }
  return ph_tm_;
}

std::vector<bool>  & baby_base::sys_pass(){
  if(read_only_ && !c_sys_pass_ && b_sys_pass_){
    b_sys_pass_->GetEntry(entry_);
    c_sys_pass_ = true;
  }
  return sys_pass_;
}

std::vector<bool>  & baby_base::trig(){
  if(read_only_ && !c_trig_ && b_trig_){
    b_trig_->GetEntry(entry_);
    c_trig_ = true;
  }
  return trig_;
}

std::vector<float>  & baby_base::els_d0(){
  if(read_only_ && !c_els_d0_ && b_els_d0_){
    b_els_d0_->GetEntry(entry_);
    c_els_d0_ = true;
  }
  return els_d0_;
}

std::vector<float>  & baby_base::els_dz(){
  if(read_only_ && !c_els_dz_ && b_els_dz_){
    b_els_dz_->GetEntry(entry_);
    c_els_dz_ = true;
  }
  return els_dz_;
}

std::vector<float>  & baby_base::els_eta(){
  if(read_only_ && !c_els_eta_ && b_els_eta_){
    b_els_eta_->GetEntry(entry_);
    c_els_eta_ = true;
  }
  return els_eta_;
}

std::vector<float>  & baby_base::els_miniso(){
  if(read_only_ && !c_els_miniso_ && b_els_miniso_){
    b_els_miniso_->GetEntry(entry_);
    c_els_miniso_ = true;
  }
  return els_miniso_;
}

std::vector<float>  & baby_base::els_phi(){
  if(read_only_ && !c_els_phi_ && b_els_phi_){
    b_els_phi_->GetEntry(entry_);
    c_els_phi_ = true;
  }
  return els_phi_;
}

std::vector<float>  & baby_base::els_pt(){
  if(read_only_ && !c_els_pt_ && b_els_pt_){
    b_els_pt_->GetEntry(entry_);
    c_els_pt_ = true;
  }
  return els_pt_;
}

std::vector<float>  & baby_base::els_reliso(){
  if(read_only_ && !c_els_reliso_ && b_els_reliso_){
    b_els_reliso_->GetEntry(entry_);
    c_els_reliso_ = true;
  }
  return els_reliso_;
}

std::vector<float>  & baby_base::els_sceta(){
  if(read_only_ && !c_els_sceta_ && b_els_sceta_){
    b_els_sceta_->GetEntry(entry_);
    c_els_sceta_ = true;
  }
  return els_sceta_;
}

std::vector<float>  & baby_base::els_vvvl_eta(){
  if(read_only_ && !c_els_vvvl_eta_ && b_els_vvvl_eta_){
    b_els_vvvl_eta_->GetEntry(entry_);
    c_els_vvvl_eta_ = true;
  }
  return els_vvvl_eta_;
}

std::vector<float>  & baby_base::els_vvvl_phi(){
  if(read_only_ && !c_els_vvvl_phi_ && b_els_vvvl_phi_){
    b_els_vvvl_phi_->GetEntry(entry_);
    c_els_vvvl_phi_ = true;
  }
  return els_vvvl_phi_;
}

std::vector<float>  & baby_base::els_vvvl_pt(){
  if(read_only_ && !c_els_vvvl_pt_ && b_els_vvvl_pt_){
    b_els_vvvl_pt_->GetEntry(entry_);
    c_els_vvvl_pt_ = true;
  }
  return els_vvvl_pt_;
}

std::vector<float>  & baby_base::fjets08_eta(){
  if(read_only_ && !c_fjets08_eta_ && b_fjets08_eta_){
    b_fjets08_eta_->GetEntry(entry_);
    c_fjets08_eta_ = true;
  }
  return fjets08_eta_;
}

std::vector<float>  & baby_base::fjets08_m(){
  if(read_only_ && !c_fjets08_m_ && b_fjets08_m_){
    b_fjets08_m_->GetEntry(entry_);
    c_fjets08_m_ = true;
  }
  return fjets08_m_;
}

std::vector<float>  & baby_base::fjets08_phi(){
  if(read_only_ && !c_fjets08_phi_ && b_fjets08_phi_){
    b_fjets08_phi_->GetEntry(entry_);
    c_fjets08_phi_ = true;
  }
  return fjets08_phi_;
}

std::vector<float>  & baby_base::fjets08_poscsv(){
  if(read_only_ && !c_fjets08_poscsv_ && b_fjets08_poscsv_){
    b_fjets08_poscsv_->GetEntry(entry_);
    c_fjets08_poscsv_ = true;
  }
  return fjets08_poscsv_;
}

std::vector<float>  & baby_base::fjets08_pt(){
  if(read_only_ && !c_fjets08_pt_ && b_fjets08_pt_){
    b_fjets08_pt_->GetEntry(entry_);
    c_fjets08_pt_ = true;
  }
  return fjets08_pt_;
}

std::vector<float>  & baby_base::fjets08_sumcsv(){
  if(read_only_ && !c_fjets08_sumcsv_ && b_fjets08_sumcsv_){
    b_fjets08_sumcsv_->GetEntry(entry_);
    c_fjets08_sumcsv_ = true;
  }
  return fjets08_sumcsv_;
}

std::vector<float>  & baby_base::fjets_eta(){
  if(read_only_ && !c_fjets_eta_ && b_fjets_eta_){
    b_fjets_eta_->GetEntry(entry_);
    c_fjets_eta_ = true;
  }
  return fjets_eta_;
}

std::vector<float>  & baby_base::fjets_m(){
  if(read_only_ && !c_fjets_m_ && b_fjets_m_){
    b_fjets_m_->GetEntry(entry_);
    c_fjets_m_ = true;
  }
  return fjets_m_;
}

std::vector<float>  & baby_base::fjets_phi(){
  if(read_only_ && !c_fjets_phi_ && b_fjets_phi_){
    b_fjets_phi_->GetEntry(entry_);
    c_fjets_phi_ = true;
  }
  return fjets_phi_;
}

std::vector<float>  & baby_base::fjets_poscsv(){
  if(read_only_ && !c_fjets_poscsv_ && b_fjets_poscsv_){
    b_fjets_poscsv_->GetEntry(entry_);
    c_fjets_poscsv_ = true;
  }
  return fjets_poscsv_;
}

std::vector<float>  & baby_base::fjets_pt(){
  if(read_only_ && !c_fjets_pt_ && b_fjets_pt_){
    b_fjets_pt_->GetEntry(entry_);
    c_fjets_pt_ = true;
  }
  return fjets_pt_;
}

std::vector<float>  & baby_base::fjets_sumcsv(){
  if(read_only_ && !c_fjets_sumcsv_ && b_fjets_sumcsv_){
    b_fjets_sumcsv_->GetEntry(entry_);
    c_fjets_sumcsv_ = true;
  }
  return fjets_sumcsv_;
}

std::vector<float>  & baby_base::jets_csv(){
  if(read_only_ && !c_jets_csv_ && b_jets_csv_){
    b_jets_csv_->GetEntry(entry_);
    c_jets_csv_ = true;
  }
  return jets_csv_;
}

std::vector<float>  & baby_base::jets_eta(){
  if(read_only_ && !c_jets_eta_ && b_jets_eta_){
    b_jets_eta_->GetEntry(entry_);
    c_jets_eta_ = true;
  }
  return jets_eta_;
}

std::vector<float>  & baby_base::jets_m(){
  if(read_only_ && !c_jets_m_ && b_jets_m_){
    b_jets_m_->GetEntry(entry_);
    c_jets_m_ = true;
  }
  return jets_m_;
}

std::vector<float>  & baby_base::jets_phi(){
  if(read_only_ && !c_jets_phi_ && b_jets_phi_){
    b_jets_phi_->GetEntry(entry_);
    c_jets_phi_ = true;
  }
  return jets_phi_;
}

std::vector<float>  & baby_base::jets_pt(){
  if(read_only_ && !c_jets_pt_ && b_jets_pt_){
    b_jets_pt_->GetEntry(entry_);
    c_jets_pt_ = true;
  }
  return jets_pt_;
}

std::vector<float>  & baby_base::jets_pt_res(){
  if(read_only_ && !c_jets_pt_res_ && b_jets_pt_res_){
    b_jets_pt_res_->GetEntry(entry_);
    c_jets_pt_res_ = true;
  }
  return jets_pt_res_;
}

std::vector<float>  & baby_base::leps_eta(){
  if(read_only_ && !c_leps_eta_ && b_leps_eta_){
    b_leps_eta_->GetEntry(entry_);
    c_leps_eta_ = true;
  }
  return leps_eta_;
}

std::vector<float>  & baby_base::leps_id(){
  if(read_only_ && !c_leps_id_ && b_leps_id_){
    b_leps_id_->GetEntry(entry_);
    c_leps_id_ = true;
  }
  return leps_id_;
}

std::vector<float>  & baby_base::leps_phi(){
  if(read_only_ && !c_leps_phi_ && b_leps_phi_){
    b_leps_phi_->GetEntry(entry_);
    c_leps_phi_ = true;
  }
  return leps_phi_;
}

std::vector<float>  & baby_base::leps_pt(){
  if(read_only_ && !c_leps_pt_ && b_leps_pt_){
    b_leps_pt_->GetEntry(entry_);
    c_leps_pt_ = true;
  }
  return leps_pt_;
}

std::vector<float>  & baby_base::mc_eta(){
  if(read_only_ && !c_mc_eta_ && b_mc_eta_){
    b_mc_eta_->GetEntry(entry_);
    c_mc_eta_ = true;
  }
  return mc_eta_;
}

std::vector<float>  & baby_base::mc_mass(){
  if(read_only_ && !c_mc_mass_ && b_mc_mass_){
    b_mc_mass_->GetEntry(entry_);
    c_mc_mass_ = true;
  }
  return mc_mass_;
}

std::vector<float>  & baby_base::mc_phi(){
  if(read_only_ && !c_mc_phi_ && b_mc_phi_){
    b_mc_phi_->GetEntry(entry_);
    c_mc_phi_ = true;
  }
  return mc_phi_;
}

std::vector<float>  & baby_base::mc_pt(){
  if(read_only_ && !c_mc_pt_ && b_mc_pt_){
    b_mc_pt_->GetEntry(entry_);
    c_mc_pt_ = true;
  }
  return mc_pt_;
}

std::vector<float>  & baby_base::mus_d0(){
  if(read_only_ && !c_mus_d0_ && b_mus_d0_){
    b_mus_d0_->GetEntry(entry_);
    c_mus_d0_ = true;
  }
  return mus_d0_;
}

std::vector<float>  & baby_base::mus_dz(){
  if(read_only_ && !c_mus_dz_ && b_mus_dz_){
    b_mus_dz_->GetEntry(entry_);
    c_mus_dz_ = true;
  }
  return mus_dz_;
}

std::vector<float>  & baby_base::mus_eta(){
  if(read_only_ && !c_mus_eta_ && b_mus_eta_){
    b_mus_eta_->GetEntry(entry_);
    c_mus_eta_ = true;
  }
  return mus_eta_;
}

std::vector<float>  & baby_base::mus_miniso(){
  if(read_only_ && !c_mus_miniso_ && b_mus_miniso_){
    b_mus_miniso_->GetEntry(entry_);
    c_mus_miniso_ = true;
  }
  return mus_miniso_;
}

std::vector<float>  & baby_base::mus_phi(){
  if(read_only_ && !c_mus_phi_ && b_mus_phi_){
    b_mus_phi_->GetEntry(entry_);
    c_mus_phi_ = true;
  }
  return mus_phi_;
}

std::vector<float>  & baby_base::mus_pt(){
  if(read_only_ && !c_mus_pt_ && b_mus_pt_){
    b_mus_pt_->GetEntry(entry_);
    c_mus_pt_ = true;
  }
  return mus_pt_;
}

std::vector<float>  & baby_base::mus_reliso(){
  if(read_only_ && !c_mus_reliso_ && b_mus_reliso_){
    b_mus_reliso_->GetEntry(entry_);
    c_mus_reliso_ = true;
  }
  return mus_reliso_;
}

std::vector<float>  & baby_base::mus_vvvl_eta(){
  if(read_only_ && !c_mus_vvvl_eta_ && b_mus_vvvl_eta_){
    b_mus_vvvl_eta_->GetEntry(entry_);
    c_mus_vvvl_eta_ = true;
  }
  return mus_vvvl_eta_;
}

std::vector<float>  & baby_base::mus_vvvl_phi(){
  if(read_only_ && !c_mus_vvvl_phi_ && b_mus_vvvl_phi_){
    b_mus_vvvl_phi_->GetEntry(entry_);
    c_mus_vvvl_phi_ = true;
  }
  return mus_vvvl_phi_;
}

std::vector<float>  & baby_base::mus_vvvl_pt(){
  if(read_only_ && !c_mus_vvvl_pt_ && b_mus_vvvl_pt_){
    b_mus_vvvl_pt_->GetEntry(entry_);
    c_mus_vvvl_pt_ = true;
  }
  return mus_vvvl_pt_;
}

std::vector<float>  & baby_base::ph_eta(){
  if(read_only_ && !c_ph_eta_ && b_ph_eta_){
    b_ph_eta_->GetEntry(entry_);
    c_ph_eta_ = true;
  }
  return ph_eta_;
}

std::vector<float>  & baby_base::ph_phi(){
  if(read_only_ && !c_ph_phi_ && b_ph_phi_){
    b_ph_phi_->GetEntry(entry_);
    c_ph_phi_ = true;
  }
  return ph_phi_;
}

std::vector<float>  & baby_base::ph_pt(){
  if(read_only_ && !c_ph_pt_ && b_ph_pt_){
    b_ph_pt_->GetEntry(entry_);
    c_ph_pt_ = true;
  }
  return ph_pt_;
}

std::vector<float>  & baby_base::sys_bctag(){
  if(read_only_ && !c_sys_bctag_ && b_sys_bctag_){
    b_sys_bctag_->GetEntry(entry_);
    c_sys_bctag_ = true;
  }
  return sys_bctag_;
}

std::vector<float>  & baby_base::sys_fs_bctag(){
  if(read_only_ && !c_sys_fs_bctag_ && b_sys_fs_bctag_){
    b_sys_fs_bctag_->GetEntry(entry_);
    c_sys_fs_bctag_ = true;
  }
  return sys_fs_bctag_;
}

std::vector<float>  & baby_base::sys_fs_lep(){
  if(read_only_ && !c_sys_fs_lep_ && b_sys_fs_lep_){
    b_sys_fs_lep_->GetEntry(entry_);
    c_sys_fs_lep_ = true;
  }
  return sys_fs_lep_;
}

std::vector<float>  & baby_base::sys_fs_udsgtag(){
  if(read_only_ && !c_sys_fs_udsgtag_ && b_sys_fs_udsgtag_){
    b_sys_fs_udsgtag_->GetEntry(entry_);
    c_sys_fs_udsgtag_ = true;
  }
  return sys_fs_udsgtag_;
}

std::vector<float>  & baby_base::sys_ht(){
  if(read_only_ && !c_sys_ht_ && b_sys_ht_){
    b_sys_ht_->GetEntry(entry_);
    c_sys_ht_ = true;
  }
  return sys_ht_;
}

std::vector<float>  & baby_base::sys_isr(){
  if(read_only_ && !c_sys_isr_ && b_sys_isr_){
    b_sys_isr_->GetEntry(entry_);
    c_sys_isr_ = true;
  }
  return sys_isr_;
}

std::vector<float>  & baby_base::sys_lep(){
  if(read_only_ && !c_sys_lep_ && b_sys_lep_){
    b_sys_lep_->GetEntry(entry_);
    c_sys_lep_ = true;
  }
  return sys_lep_;
}

std::vector<float>  & baby_base::sys_met(){
  if(read_only_ && !c_sys_met_ && b_sys_met_){
    b_sys_met_->GetEntry(entry_);
    c_sys_met_ = true;
  }
  return sys_met_;
}

std::vector<float>  & baby_base::sys_mj(){
  if(read_only_ && !c_sys_mj_ && b_sys_mj_){
    b_sys_mj_->GetEntry(entry_);
    c_sys_mj_ = true;
  }
  return sys_mj_;
}

std::vector<float>  & baby_base::sys_mt(){
  if(read_only_ && !c_sys_mt_ && b_sys_mt_){
    b_sys_mt_->GetEntry(entry_);
    c_sys_mt_ = true;
  }
  return sys_mt_;
}

std::vector<float>  & baby_base::sys_muf(){
  if(read_only_ && !c_sys_muf_ && b_sys_muf_){
    b_sys_muf_->GetEntry(entry_);
    c_sys_muf_ = true;
  }
  return sys_muf_;
}

std::vector<float>  & baby_base::sys_mur(){
  if(read_only_ && !c_sys_mur_ && b_sys_mur_){
    b_sys_mur_->GetEntry(entry_);
    c_sys_mur_ = true;
  }
  return sys_mur_;
}

std::vector<float>  & baby_base::sys_murf(){
  if(read_only_ && !c_sys_murf_ && b_sys_murf_){
    b_sys_murf_->GetEntry(entry_);
    c_sys_murf_ = true;
  }
  return sys_murf_;
}

std::vector<float>  & baby_base::sys_pdf(){
  if(read_only_ && !c_sys_pdf_ && b_sys_pdf_){
    b_sys_pdf_->GetEntry(entry_);
    c_sys_pdf_ = true;
  }
  return sys_pdf_;
}

std::vector<float>  & baby_base::sys_trig(){
  if(read_only_ && !c_sys_trig_ && b_sys_trig_){
    b_sys_trig_->GetEntry(entry_);
    c_sys_trig_ = true;
  }
  return sys_trig_;
}

std::vector<float>  & baby_base::sys_udsgtag(){
  if(read_only_ && !c_sys_udsgtag_ && b_sys_udsgtag_){
    b_sys_udsgtag_->GetEntry(entry_);
    c_sys_udsgtag_ = true;
  }
  return sys_udsgtag_;
}

std::vector<float>  & baby_base::tks_eta(){
  if(read_only_ && !c_tks_eta_ && b_tks_eta_){
    b_tks_eta_->GetEntry(entry_);
    c_tks_eta_ = true;
  }
  return tks_eta_;
}

std::vector<float>  & baby_base::tks_miniso(){
  if(read_only_ && !c_tks_miniso_ && b_tks_miniso_){
    b_tks_miniso_->GetEntry(entry_);
    c_tks_miniso_ = true;
  }
  return tks_miniso_;
}

std::vector<float>  & baby_base::tks_mt(){
  if(read_only_ && !c_tks_mt_ && b_tks_mt_){
    b_tks_mt_->GetEntry(entry_);
    c_tks_mt_ = true;
  }
  return tks_mt_;
}

std::vector<float>  & baby_base::tks_mt2(){
  if(read_only_ && !c_tks_mt2_ && b_tks_mt2_){
    b_tks_mt2_->GetEntry(entry_);
    c_tks_mt2_ = true;
  }
  return tks_mt2_;
}

std::vector<float>  & baby_base::tks_phi(){
  if(read_only_ && !c_tks_phi_ && b_tks_phi_){
    b_tks_phi_->GetEntry(entry_);
    c_tks_phi_ = true;
  }
  return tks_phi_;
}

std::vector<float>  & baby_base::tks_pt(){
  if(read_only_ && !c_tks_pt_ && b_tks_pt_){
    b_tks_pt_->GetEntry(entry_);
    c_tks_pt_ = true;
  }
  return tks_pt_;
}

std::vector<float>  & baby_base::trig_prescale(){
  if(read_only_ && !c_trig_prescale_ && b_trig_prescale_){
    b_trig_prescale_->GetEntry(entry_);
    c_trig_prescale_ = true;
  }
  return trig_prescale_;
}

std::vector<float>  & baby_base::w_pdf(){
  if(read_only_ && !c_w_pdf_ && b_w_pdf_){
    b_w_pdf_->GetEntry(entry_);
    c_w_pdf_ = true;
  }
  return w_pdf_;
}

std::vector<int>  & baby_base::els_charge(){
  if(read_only_ && !c_els_charge_ && b_els_charge_){
    b_els_charge_->GetEntry(entry_);
    c_els_charge_ = true;
  }
  return els_charge_;
}

std::vector<int>  & baby_base::fjets08_btags(){
  if(read_only_ && !c_fjets08_btags_ && b_fjets08_btags_){
    b_fjets08_btags_->GetEntry(entry_);
    c_fjets08_btags_ = true;
  }
  return fjets08_btags_;
}

std::vector<int>  & baby_base::fjets08_nconst(){
  if(read_only_ && !c_fjets08_nconst_ && b_fjets08_nconst_){
    b_fjets08_nconst_->GetEntry(entry_);
    c_fjets08_nconst_ = true;
  }
  return fjets08_nconst_;
}

std::vector<int>  & baby_base::fjets_btags(){
  if(read_only_ && !c_fjets_btags_ && b_fjets_btags_){
    b_fjets_btags_->GetEntry(entry_);
    c_fjets_btags_ = true;
  }
  return fjets_btags_;
}

std::vector<int>  & baby_base::fjets_nconst(){
  if(read_only_ && !c_fjets_nconst_ && b_fjets_nconst_){
    b_fjets_nconst_->GetEntry(entry_);
    c_fjets_nconst_ = true;
  }
  return fjets_nconst_;
}

std::vector<int>  & baby_base::jets_fjet08_index(){
  if(read_only_ && !c_jets_fjet08_index_ && b_jets_fjet08_index_){
    b_jets_fjet08_index_->GetEntry(entry_);
    c_jets_fjet08_index_ = true;
  }
  return jets_fjet08_index_;
}

std::vector<int>  & baby_base::jets_fjet_index(){
  if(read_only_ && !c_jets_fjet_index_ && b_jets_fjet_index_){
    b_jets_fjet_index_->GetEntry(entry_);
    c_jets_fjet_index_ = true;
  }
  return jets_fjet_index_;
}

std::vector<int>  & baby_base::jets_hflavor(){
  if(read_only_ && !c_jets_hflavor_ && b_jets_hflavor_){
    b_jets_hflavor_->GetEntry(entry_);
    c_jets_hflavor_ = true;
  }
  return jets_hflavor_;
}

std::vector<int>  & baby_base::mc_id(){
  if(read_only_ && !c_mc_id_ && b_mc_id_){
    b_mc_id_->GetEntry(entry_);
    c_mc_id_ = true;
  }
  return mc_id_;
}

std::vector<int>  & baby_base::mc_mom(){
  if(read_only_ && !c_mc_mom_ && b_mc_mom_){
    b_mc_mom_->GetEntry(entry_);
    c_mc_mom_ = true;
  }
  return mc_mom_;
}

std::vector<int>  & baby_base::mus_charge(){
  if(read_only_ && !c_mus_charge_ && b_mus_charge_){
    b_mus_charge_->GetEntry(entry_);
    c_mus_charge_ = true;
  }
  return mus_charge_;
}

std::vector<int>  & baby_base::sys_nbm(){
  if(read_only_ && !c_sys_nbm_ && b_sys_nbm_){
    b_sys_nbm_->GetEntry(entry_);
    c_sys_nbm_ = true;
  }
  return sys_nbm_;
}

std::vector<int>  & baby_base::sys_njets(){
  if(read_only_ && !c_sys_njets_ && b_sys_njets_){
    b_sys_njets_->GetEntry(entry_);
    c_sys_njets_ = true;
  }
  return sys_njets_;
}

std::vector<int>  & baby_base::tks_pdg(){
  if(read_only_ && !c_tks_pdg_ && b_tks_pdg_){
    b_tks_pdg_->GetEntry(entry_);
    c_tks_pdg_ = true;
  }
  return tks_pdg_;
}

#include "baby_basic.hh"
baby_base* NewTree(const std::type_info &type){

  if(type == typeid(baby_base)) return new baby_base;
  else if(type == typeid(baby_basic)) return static_cast<baby_base*>(new baby_basic);
  else return new baby_base;
}

