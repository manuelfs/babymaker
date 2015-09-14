// baby: base class to handle reduced tree ntuples
// File generated with generate_baby.exe

#ifndef H_BABY
#define H_BABY

#include <vector>
#include <string>

#include "TTree.h"
#include "TChain.h"
#include "TTreeFormula.h"

class baby{
public:
  baby(); // Constructor to create tree
  baby(const std::string &filename); // Constructor to read tree

  int Add(const std::string &filename);
  long GetEntries() const;
  virtual void GetEntry(const long entry);
  bool PassString(TString cut);

  virtual void Fill();
  void Write();

  virtual std::string Type() const;

  static const double bad_val_;

  virtual ~baby();


  __attribute__((noreturn)) virtual bool  const & hfjet() const;
  __attribute__((noreturn)) virtual bool  & hfjet();
  __attribute__((noreturn)) virtual bool  const & json_golden() const;
  __attribute__((noreturn)) virtual bool  & json_golden();
  __attribute__((noreturn)) virtual bool  const & pass() const;
  __attribute__((noreturn)) virtual bool  & pass();
  __attribute__((noreturn)) virtual bool  const & pass_cschalo() const;
  __attribute__((noreturn)) virtual bool  & pass_cschalo();
  __attribute__((noreturn)) virtual bool  const & pass_eebadsc() const;
  __attribute__((noreturn)) virtual bool  & pass_eebadsc();
  __attribute__((noreturn)) virtual bool  const & pass_goodv() const;
  __attribute__((noreturn)) virtual bool  & pass_goodv();
  __attribute__((noreturn)) virtual bool  const & pass_hbhe() const;
  __attribute__((noreturn)) virtual bool  & pass_hbhe();
  __attribute__((noreturn)) virtual bool  const & pass_jets() const;
  __attribute__((noreturn)) virtual bool  & pass_jets();
  __attribute__((noreturn)) virtual float  const & dphi_wlep() const;
  __attribute__((noreturn)) virtual float  & dphi_wlep();
  __attribute__((noreturn)) virtual float  const & dphi_wlep_reliso() const;
  __attribute__((noreturn)) virtual float  & dphi_wlep_reliso();
  __attribute__((noreturn)) virtual float  const & elel_m() const;
  __attribute__((noreturn)) virtual float  & elel_m();
  __attribute__((noreturn)) virtual float  const & elel_pt1() const;
  __attribute__((noreturn)) virtual float  & elel_pt1();
  __attribute__((noreturn)) virtual float  const & elel_pt2() const;
  __attribute__((noreturn)) virtual float  & elel_pt2();
  __attribute__((noreturn)) virtual float  const & elel_zpt() const;
  __attribute__((noreturn)) virtual float  & elel_zpt();
  __attribute__((noreturn)) virtual float  const & elelv_m() const;
  __attribute__((noreturn)) virtual float  & elelv_m();
  __attribute__((noreturn)) virtual float  const & elelv_pt1() const;
  __attribute__((noreturn)) virtual float  & elelv_pt1();
  __attribute__((noreturn)) virtual float  const & elelv_pt2() const;
  __attribute__((noreturn)) virtual float  & elelv_pt2();
  __attribute__((noreturn)) virtual float  const & elelv_zpt() const;
  __attribute__((noreturn)) virtual float  & elelv_zpt();
  __attribute__((noreturn)) virtual float  const & ht() const;
  __attribute__((noreturn)) virtual float  & ht();
  __attribute__((noreturn)) virtual float  const & ht40() const;
  __attribute__((noreturn)) virtual float  & ht40();
  __attribute__((noreturn)) virtual float  const & ht_hf() const;
  __attribute__((noreturn)) virtual float  & ht_hf();
  __attribute__((noreturn)) virtual float  const & ht_hlt() const;
  __attribute__((noreturn)) virtual float  & ht_hlt();
  __attribute__((noreturn)) virtual float  const & ht_nohf() const;
  __attribute__((noreturn)) virtual float  & ht_nohf();
  __attribute__((noreturn)) virtual float  const & ht_ra2b() const;
  __attribute__((noreturn)) virtual float  & ht_ra2b();
  __attribute__((noreturn)) virtual float  const & ht_reliso() const;
  __attribute__((noreturn)) virtual float  & ht_reliso();
  __attribute__((noreturn)) virtual float  const & lep_eta() const;
  __attribute__((noreturn)) virtual float  & lep_eta();
  __attribute__((noreturn)) virtual float  const & lep_eta_reliso() const;
  __attribute__((noreturn)) virtual float  & lep_eta_reliso();
  __attribute__((noreturn)) virtual float  const & lep_phi() const;
  __attribute__((noreturn)) virtual float  & lep_phi();
  __attribute__((noreturn)) virtual float  const & lep_phi_reliso() const;
  __attribute__((noreturn)) virtual float  & lep_phi_reliso();
  __attribute__((noreturn)) virtual float  const & lep_pt() const;
  __attribute__((noreturn)) virtual float  & lep_pt();
  __attribute__((noreturn)) virtual float  const & lep_pt_reliso() const;
  __attribute__((noreturn)) virtual float  & lep_pt_reliso();
  __attribute__((noreturn)) virtual float  const & met() const;
  __attribute__((noreturn)) virtual float  & met();
  __attribute__((noreturn)) virtual float  const & met_hf() const;
  __attribute__((noreturn)) virtual float  & met_hf();
  __attribute__((noreturn)) virtual float  const & met_hf_phi() const;
  __attribute__((noreturn)) virtual float  & met_hf_phi();
  __attribute__((noreturn)) virtual float  const & met_mini() const;
  __attribute__((noreturn)) virtual float  & met_mini();
  __attribute__((noreturn)) virtual float  const & met_mini_phi() const;
  __attribute__((noreturn)) virtual float  & met_mini_phi();
  __attribute__((noreturn)) virtual float  const & met_nohf() const;
  __attribute__((noreturn)) virtual float  & met_nohf();
  __attribute__((noreturn)) virtual float  const & met_nohf_phi() const;
  __attribute__((noreturn)) virtual float  & met_nohf_phi();
  __attribute__((noreturn)) virtual float  const & met_nohf_sumEt() const;
  __attribute__((noreturn)) virtual float  & met_nohf_sumEt();
  __attribute__((noreturn)) virtual float  const & met_phi() const;
  __attribute__((noreturn)) virtual float  & met_phi();
  __attribute__((noreturn)) virtual float  const & mht() const;
  __attribute__((noreturn)) virtual float  & mht();
  __attribute__((noreturn)) virtual float  const & mht_ra2b() const;
  __attribute__((noreturn)) virtual float  & mht_ra2b();
  __attribute__((noreturn)) virtual float  const & mindphin_metjet() const;
  __attribute__((noreturn)) virtual float  & mindphin_metjet();
  __attribute__((noreturn)) virtual float  const & mj() const;
  __attribute__((noreturn)) virtual float  & mj();
  __attribute__((noreturn)) virtual float  const & mj08() const;
  __attribute__((noreturn)) virtual float  & mj08();
  __attribute__((noreturn)) virtual float  const & mt() const;
  __attribute__((noreturn)) virtual float  & mt();
  __attribute__((noreturn)) virtual float  const & mt_reliso() const;
  __attribute__((noreturn)) virtual float  & mt_reliso();
  __attribute__((noreturn)) virtual float  const & mumu_m() const;
  __attribute__((noreturn)) virtual float  & mumu_m();
  __attribute__((noreturn)) virtual float  const & mumu_pt1() const;
  __attribute__((noreturn)) virtual float  & mumu_pt1();
  __attribute__((noreturn)) virtual float  const & mumu_pt2() const;
  __attribute__((noreturn)) virtual float  & mumu_pt2();
  __attribute__((noreturn)) virtual float  const & mumu_zpt() const;
  __attribute__((noreturn)) virtual float  & mumu_zpt();
  __attribute__((noreturn)) virtual float  const & mumuv_m() const;
  __attribute__((noreturn)) virtual float  & mumuv_m();
  __attribute__((noreturn)) virtual float  const & mumuv_pt1() const;
  __attribute__((noreturn)) virtual float  & mumuv_pt1();
  __attribute__((noreturn)) virtual float  const & mumuv_pt2() const;
  __attribute__((noreturn)) virtual float  & mumuv_pt2();
  __attribute__((noreturn)) virtual float  const & mumuv_zpt() const;
  __attribute__((noreturn)) virtual float  & mumuv_zpt();
  __attribute__((noreturn)) virtual float  const & ntrupv_mean() const;
  __attribute__((noreturn)) virtual float  & ntrupv_mean();
  __attribute__((noreturn)) virtual float  const & onht() const;
  __attribute__((noreturn)) virtual float  & onht();
  __attribute__((noreturn)) virtual float  const & onmaxel() const;
  __attribute__((noreturn)) virtual float  & onmaxel();
  __attribute__((noreturn)) virtual float  const & onmaxmu() const;
  __attribute__((noreturn)) virtual float  & onmaxmu();
  __attribute__((noreturn)) virtual float  const & onmet() const;
  __attribute__((noreturn)) virtual float  & onmet();
  __attribute__((noreturn)) virtual float  const & st() const;
  __attribute__((noreturn)) virtual float  & st();
  __attribute__((noreturn)) virtual float  const & st_reliso() const;
  __attribute__((noreturn)) virtual float  & st_reliso();
  __attribute__((noreturn)) virtual float  const & trutop1_phi() const;
  __attribute__((noreturn)) virtual float  & trutop1_phi();
  __attribute__((noreturn)) virtual float  const & trutop1_pt() const;
  __attribute__((noreturn)) virtual float  & trutop1_pt();
  __attribute__((noreturn)) virtual float  const & trutop2_phi() const;
  __attribute__((noreturn)) virtual float  & trutop2_phi();
  __attribute__((noreturn)) virtual float  const & trutop2_pt() const;
  __attribute__((noreturn)) virtual float  & trutop2_pt();
  __attribute__((noreturn)) virtual float  const & weight() const;
  __attribute__((noreturn)) virtual float  & weight();
  __attribute__((noreturn)) virtual float  const & wpu() const;
  __attribute__((noreturn)) virtual float  & wpu();
  __attribute__((noreturn)) virtual int  const & event() const;
  __attribute__((noreturn)) virtual int  & event();
  __attribute__((noreturn)) virtual int  const & lep_charge() const;
  __attribute__((noreturn)) virtual int  & lep_charge();
  __attribute__((noreturn)) virtual int  const & lep_charge_reliso() const;
  __attribute__((noreturn)) virtual int  & lep_charge_reliso();
  __attribute__((noreturn)) virtual int  const & lumiblock() const;
  __attribute__((noreturn)) virtual int  & lumiblock();
  __attribute__((noreturn)) virtual int  const & nbl() const;
  __attribute__((noreturn)) virtual int  & nbl();
  __attribute__((noreturn)) virtual int  const & nbl40() const;
  __attribute__((noreturn)) virtual int  & nbl40();
  __attribute__((noreturn)) virtual int  const & nbm() const;
  __attribute__((noreturn)) virtual int  & nbm();
  __attribute__((noreturn)) virtual int  const & nbm40() const;
  __attribute__((noreturn)) virtual int  & nbm40();
  __attribute__((noreturn)) virtual int  const & nbm_ra2b() const;
  __attribute__((noreturn)) virtual int  & nbm_ra2b();
  __attribute__((noreturn)) virtual int  const & nbt() const;
  __attribute__((noreturn)) virtual int  & nbt();
  __attribute__((noreturn)) virtual int  const & nbt40() const;
  __attribute__((noreturn)) virtual int  & nbt40();
  __attribute__((noreturn)) virtual int  const & nels() const;
  __attribute__((noreturn)) virtual int  & nels();
  __attribute__((noreturn)) virtual int  const & nels_reliso() const;
  __attribute__((noreturn)) virtual int  & nels_reliso();
  __attribute__((noreturn)) virtual int  const & nevents() const;
  __attribute__((noreturn)) virtual int  & nevents();
  __attribute__((noreturn)) virtual int  const & nfjets() const;
  __attribute__((noreturn)) virtual int  & nfjets();
  __attribute__((noreturn)) virtual int  const & nfjets08() const;
  __attribute__((noreturn)) virtual int  & nfjets08();
  __attribute__((noreturn)) virtual int  const & njets() const;
  __attribute__((noreturn)) virtual int  & njets();
  __attribute__((noreturn)) virtual int  const & njets40() const;
  __attribute__((noreturn)) virtual int  & njets40();
  __attribute__((noreturn)) virtual int  const & njets_hf() const;
  __attribute__((noreturn)) virtual int  & njets_hf();
  __attribute__((noreturn)) virtual int  const & njets_nohf() const;
  __attribute__((noreturn)) virtual int  & njets_nohf();
  __attribute__((noreturn)) virtual int  const & njets_ra2b() const;
  __attribute__((noreturn)) virtual int  & njets_ra2b();
  __attribute__((noreturn)) virtual int  const & nleps() const;
  __attribute__((noreturn)) virtual int  & nleps();
  __attribute__((noreturn)) virtual int  const & nleps_reliso() const;
  __attribute__((noreturn)) virtual int  & nleps_reliso();
  __attribute__((noreturn)) virtual int  const & nmus() const;
  __attribute__((noreturn)) virtual int  & nmus();
  __attribute__((noreturn)) virtual int  const & nmus_reliso() const;
  __attribute__((noreturn)) virtual int  & nmus_reliso();
  __attribute__((noreturn)) virtual int  const & npv() const;
  __attribute__((noreturn)) virtual int  & npv();
  __attribute__((noreturn)) virtual int  const & ntruels() const;
  __attribute__((noreturn)) virtual int  & ntruels();
  __attribute__((noreturn)) virtual int  const & ntruleps() const;
  __attribute__((noreturn)) virtual int  & ntruleps();
  __attribute__((noreturn)) virtual int  const & ntrumeisr() const;
  __attribute__((noreturn)) virtual int  & ntrumeisr();
  __attribute__((noreturn)) virtual int  const & ntrumus() const;
  __attribute__((noreturn)) virtual int  & ntrumus();
  __attribute__((noreturn)) virtual int  const & ntrupv() const;
  __attribute__((noreturn)) virtual int  & ntrupv();
  __attribute__((noreturn)) virtual int  const & ntrutaush() const;
  __attribute__((noreturn)) virtual int  & ntrutaush();
  __attribute__((noreturn)) virtual int  const & ntrutausl() const;
  __attribute__((noreturn)) virtual int  & ntrutausl();
  __attribute__((noreturn)) virtual int  const & nvels() const;
  __attribute__((noreturn)) virtual int  & nvels();
  __attribute__((noreturn)) virtual int  const & nvels_reliso() const;
  __attribute__((noreturn)) virtual int  & nvels_reliso();
  __attribute__((noreturn)) virtual int  const & nvleps() const;
  __attribute__((noreturn)) virtual int  & nvleps();
  __attribute__((noreturn)) virtual int  const & nvmus() const;
  __attribute__((noreturn)) virtual int  & nvmus();
  __attribute__((noreturn)) virtual int  const & nvmus_reliso() const;
  __attribute__((noreturn)) virtual int  & nvmus_reliso();
  __attribute__((noreturn)) virtual int  const & run() const;
  __attribute__((noreturn)) virtual int  & run();
  __attribute__((noreturn)) virtual std::vector<bool>  const & els_ispf() const;
  __attribute__((noreturn)) virtual std::vector<bool>  & els_ispf();
  __attribute__((noreturn)) virtual std::vector<bool>  const & els_sigid() const;
  __attribute__((noreturn)) virtual std::vector<bool>  & els_sigid();
  __attribute__((noreturn)) virtual std::vector<bool>  const & els_tru_tm() const;
  __attribute__((noreturn)) virtual std::vector<bool>  & els_tru_tm();
  __attribute__((noreturn)) virtual std::vector<bool>  const & jets_islep() const;
  __attribute__((noreturn)) virtual std::vector<bool>  & jets_islep();
  __attribute__((noreturn)) virtual std::vector<bool>  const & mus_sigid() const;
  __attribute__((noreturn)) virtual std::vector<bool>  & mus_sigid();
  __attribute__((noreturn)) virtual std::vector<bool>  const & mus_tru_tm() const;
  __attribute__((noreturn)) virtual std::vector<bool>  & mus_tru_tm();
  __attribute__((noreturn)) virtual std::vector<bool>  const & trig() const;
  __attribute__((noreturn)) virtual std::vector<bool>  & trig();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_d0() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_d0();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_dz() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_dz();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_eta() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_eta();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_miniso() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_miniso();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_mt() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_mt();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_phi() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_phi();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_pt() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_pt();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_reliso() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_reliso();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_sceta() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_sceta();
  __attribute__((noreturn)) virtual std::vector<float>  const & els_tru_dr() const;
  __attribute__((noreturn)) virtual std::vector<float>  & els_tru_dr();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets08_eta() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets08_eta();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets08_m() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets08_m();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets08_phi() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets08_phi();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets08_poscsv() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets08_poscsv();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets08_pt() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets08_pt();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets08_sumcsv() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets08_sumcsv();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets_eta() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets_eta();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets_m() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets_m();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets_phi() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets_phi();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets_poscsv() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets_poscsv();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets_pt() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets_pt();
  __attribute__((noreturn)) virtual std::vector<float>  const & fjets_sumcsv() const;
  __attribute__((noreturn)) virtual std::vector<float>  & fjets_sumcsv();
  __attribute__((noreturn)) virtual std::vector<float>  const & jets_csv() const;
  __attribute__((noreturn)) virtual std::vector<float>  & jets_csv();
  __attribute__((noreturn)) virtual std::vector<float>  const & jets_eta() const;
  __attribute__((noreturn)) virtual std::vector<float>  & jets_eta();
  __attribute__((noreturn)) virtual std::vector<float>  const & jets_id() const;
  __attribute__((noreturn)) virtual std::vector<float>  & jets_id();
  __attribute__((noreturn)) virtual std::vector<float>  const & jets_m() const;
  __attribute__((noreturn)) virtual std::vector<float>  & jets_m();
  __attribute__((noreturn)) virtual std::vector<float>  const & jets_phi() const;
  __attribute__((noreturn)) virtual std::vector<float>  & jets_phi();
  __attribute__((noreturn)) virtual std::vector<float>  const & jets_pt() const;
  __attribute__((noreturn)) virtual std::vector<float>  & jets_pt();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_d0() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_d0();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_dz() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_dz();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_eta() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_eta();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_miniso() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_miniso();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_mt() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_mt();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_phi() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_phi();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_pt() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_pt();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_reliso() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_reliso();
  __attribute__((noreturn)) virtual std::vector<float>  const & mus_tru_dr() const;
  __attribute__((noreturn)) virtual std::vector<float>  & mus_tru_dr();
  __attribute__((noreturn)) virtual std::vector<float>  const & trig_prescale() const;
  __attribute__((noreturn)) virtual std::vector<float>  & trig_prescale();
  __attribute__((noreturn)) virtual std::vector<int>  const & els_charge() const;
  __attribute__((noreturn)) virtual std::vector<int>  & els_charge();
  __attribute__((noreturn)) virtual std::vector<int>  const & els_tru_id() const;
  __attribute__((noreturn)) virtual std::vector<int>  & els_tru_id();
  __attribute__((noreturn)) virtual std::vector<int>  const & els_tru_momid() const;
  __attribute__((noreturn)) virtual std::vector<int>  & els_tru_momid();
  __attribute__((noreturn)) virtual std::vector<int>  const & fjets08_btags() const;
  __attribute__((noreturn)) virtual std::vector<int>  & fjets08_btags();
  __attribute__((noreturn)) virtual std::vector<int>  const & fjets08_nconst() const;
  __attribute__((noreturn)) virtual std::vector<int>  & fjets08_nconst();
  __attribute__((noreturn)) virtual std::vector<int>  const & fjets_btags() const;
  __attribute__((noreturn)) virtual std::vector<int>  & fjets_btags();
  __attribute__((noreturn)) virtual std::vector<int>  const & fjets_nconst() const;
  __attribute__((noreturn)) virtual std::vector<int>  & fjets_nconst();
  __attribute__((noreturn)) virtual std::vector<int>  const & jets_fjet08_index() const;
  __attribute__((noreturn)) virtual std::vector<int>  & jets_fjet08_index();
  __attribute__((noreturn)) virtual std::vector<int>  const & jets_fjet_index() const;
  __attribute__((noreturn)) virtual std::vector<int>  & jets_fjet_index();
  __attribute__((noreturn)) virtual std::vector<int>  const & mus_charge() const;
  __attribute__((noreturn)) virtual std::vector<int>  & mus_charge();
  __attribute__((noreturn)) virtual std::vector<int>  const & mus_tru_id() const;
  __attribute__((noreturn)) virtual std::vector<int>  & mus_tru_id();
  __attribute__((noreturn)) virtual std::vector<int>  const & mus_tru_momid() const;
  __attribute__((noreturn)) virtual std::vector<int>  & mus_tru_momid();

protected:
  TChain chain_;
  TTree tree_;
  long entry_;
  const bool read_only_;

private:
  class VectorLoader{
  public:
    VectorLoader();
  private:
    static bool loaded_;
  };

  static VectorLoader vl_;
};

baby* NewTree(const std::type_info &type);

#include "baby_basic.hpp"
#include "baby_quick.hpp"

#endif
