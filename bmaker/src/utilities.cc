// Common utilities
 
#include <cmath>
#include <iostream>
#include "babymaker/bmaker/interface/utilities.hh"

using namespace std;

namespace utilities{

  bool greaterPt(const reco::Candidate *a, const reco::Candidate *b) {
    return a->pt() > b->pt();
  }

  bool greaterM(const fastjet::PseudoJet &a, const fastjet::PseudoJet &b){
    return a.m() > b.m();
  }

  string execute(const string &cmd){
    FILE *pipe = popen(cmd.c_str(), "r");
    if(!pipe) throw runtime_error("Could not open pipe.");
    const size_t buffer_size = 128;
    char buffer[buffer_size];
    string result = "";
    while(!feof(pipe)){
      if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
    }

    pclose(pipe);
    return result;
  }

  TString roundNumber(double num, int decimals, double denom){
    if(denom==0) return " - ";
    double neg = 1; if(num*denom<0) neg = -1;
    num /= neg*denom; num += 0.5*pow(10.,-decimals);
    long num_int = static_cast<long>(num);
    long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
    TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
    TString result="";
    if(neg<0) result+="-";
    result+= num_int;
    if(decimals>0) {
      result+="."; result+=s_dec;
    }

    TString afterdot = result;
    afterdot.Remove(0,afterdot.First(".")+1);
    for(int i=0; i<decimals-afterdot.Length(); i++)
      result += "0";
    return result;
  }

  float crossSection(const TString &file){
    float xsec(0.);

    if(file.Contains("Run2015"))   xsec = 1.;


    // From https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVgluglu
    if(file.Contains("T1tttt") && file.Contains("825_"))   xsec = 1.2167;
    if(file.Contains("T1tttt") && file.Contains("1025_"))  xsec = 0.272778;
    if(file.Contains("T1tttt") && file.Contains("1150_"))  xsec = 0.117687;
    if(file.Contains("T1tttt") && file.Contains("1200_"))  xsec = 0.0856418;
    if(file.Contains("T1tttt") && file.Contains("1500_"))  xsec = 0.0141903;

    if(file.Contains("T2tt") && file.Contains("650_"))  xsec = 0.107045;
    if(file.Contains("T2tt") && file.Contains("850_"))  xsec = 0.0189612;

    if(file.Contains("SMS-T2tt_2J_mStop-425_mLSP-325"))  xsec = 1.31169;
    if(file.Contains("SMS-T2tt_2J_mStop-500_mLSP-325"))  xsec = 0.51848;
    if(file.Contains("SMS-T1bbbb_2J_mGl-1500_mLSP-100"))  xsec = 0.0141903;
    if(file.Contains("SMS-T1bbbb_2J_mGl-1000_mLSP-900"))  xsec = 0.325388;
    if(file.Contains("SMS-T1qqqq_2J_mGl-1400_mLSP-100"))  xsec = 0.0252977;
    if(file.Contains("SMS-T1qqqq_2J_mGl-1000_mLSP-800"))  xsec = 0.325388;
    if(file.Contains("SMS-T2bb_2J_mStop-600_mLSP-580"))  xsec = 0.174599;
    if(file.Contains("SMS-T2bb_2J_mStop-900_mLSP-100"))  xsec = 0.0128895;

    //  Cross-section taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    // Alternative option: https://twiki.cern.ch/twiki/bin/view/Sandbox/FullNNLOcrossSections#Top_cross_section_for_13_TeV
    if(file.Contains("TTJets_Tune") || file.Contains("TT_"))  xsec = 815.96;
    if(file.Contains("TTJets_HT")){//LO cross sections with k-factor of 1.625 already applied
      if(file.Contains("2500toInf")) xsec = 0.0023234211;
      if(file.Contains("1200to2500")) xsec = 0.194972521;
      if(file.Contains("800to1200")) xsec = 1.07722318;
      if(file.Contains("600to800")) xsec = 2.61537118;

    }

    if(file.Contains("TTJets_DiLept")) xsec = 85.66; // (3*0.108)^2*815.96
    if(file.Contains("TTJets_SingleLept")) xsec = 178.7; //(1- ((1-3*0.108)^2+(3*0.108)^2))*815.96*0.5 per half

                

    
    // From https://cms-pdmv.cern.ch/mcm
    // k-factors from https://mangano.web.cern.ch/mangano/public/MECCA/samples_50ns_miniaod.txt
    // k-factors are ratio of https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    // NLO/NNLO cross-sections to that of an inclusive sample in mcm at lower order (LO/NLO)

    if(file.Contains("WJetsToLNu") && file.Contains("amcatnloFXFX")) xsec=61526.7; //NNLO from Lesya's summary table 

    //cross-section per slice changed due to change in genHT definition
    if(file.Contains("WJetsToLNu_HT-100To200"))  xsec = 1347.*1.21; //updated based on MCM
    if(file.Contains("WJetsToLNu_HT-200To400"))  xsec = 360.*1.21;
    if(file.Contains("WJetsToLNu_HT-400To600"))  xsec = 48.98*1.21;
    if(file.Contains("WJetsToLNu_HT-600ToInf"))  xsec = 18.77*1.21;
    if(file.Contains("WToENu"))   xsec = 16000.0;
    if(file.Contains("WToMuNu"))  xsec = 16100.0;

    if(file.Contains("QCD_HT-100To250_13TeV-madgraph"))  xsec = 28730000.;
    if(file.Contains("QCD_HT_250To500_13TeV-madgraph"))  xsec = 670500.0;
    if(file.Contains("QCD_HT-500To1000_13TeV-madgraph")) xsec = 26740.0;
    if(file.Contains("QCD_HT_1000ToInf_13TeV-madgraph")) xsec = 769.7;

    if(file.Contains("QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"))   xsec = 1735000;
    if(file.Contains("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"))   xsec = 366800;
    if(file.Contains("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"))   xsec = 29370;
    if(file.Contains("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"))  xsec = 6524;
    if(file.Contains("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")) xsec = 1064;
    if(file.Contains("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")) xsec = 121.5;
    if(file.Contains("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"))  xsec = 25.42;

    if(file.Contains("QCD_Pt-1800_")) xsec = 0.1091;

    // only thing that changed for Spring 15 is "QCD_Pt-" --> "QCD_Pt_"
    if (file.Contains("QCD_Pt")){ 
      if(file.Contains("5to10"))      xsec = 80710000000;
      if(file.Contains("10to15"))     xsec = 7528000000;
      if(file.Contains("15to30"))     xsec = 2237000000;
      if(file.Contains("30to50"))     xsec = 161500000;
      if(file.Contains("50to80"))     xsec = 22110000;
      if(file.Contains("80to120"))   xsec = 2762530;// xsec = 3000114.3;
      if(file.Contains("120to170"))   xsec = 471100;//493200;
      if(file.Contains("170to300"))   xsec = 117276;//120300;
      if(file.Contains("300to470"))   xsec = 7823;//7475;
      if(file.Contains("470to600"))   xsec = 648.2;//587.1;
      if(file.Contains("600to800"))   xsec = 186.9;//167;
      if(file.Contains("800to1000"))  xsec = 32.293;//28.25;
      if(file.Contains("1000to1400")) xsec = 9.4183;//8.195;
      if(file.Contains("1400to1800")) xsec = 0.84265;// 0.7346;
      if(file.Contains("1800to2400")) xsec = 0.114943;//0.102;
      if(file.Contains("2400to3200")) xsec = 0.00682981;//0.00644;
      if(file.Contains("3200toInf"))       xsec = 0.000165445;// 0.000163;
    }

    // Cross sections from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
    // multiplied by BF(W->mu,e,tau) = 0.324
    if (file.Contains("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8"))     xsec = 3.34;
    if (file.Contains("ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8")) xsec = 26.23;
    if (file.Contains("ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8"))     xsec = 44.07;
    if (file.Contains("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8"))     xsec = 35.8;
    if (file.Contains("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8"))     xsec = 35.8; 

    if(file.Contains("DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV")) xsec = 18610*1.23;
    if(file.Contains("DYJetsToLL_M-50_TuneCUETP8M1_13TeV"))     xsec = 4895*1.23;

    if(file.Contains("DYJetsToLL_M-50_HT-100to200"))    xsec = 139.4*1.23;
    if(file.Contains("DYJetsToLL_M-50_HT-200to400"))    xsec = 42.75*1.23;
    if(file.Contains("DYJetsToLL_M-50_HT-400to600"))    xsec = 5.497*1.23;
    if(file.Contains("DYJetsToLL_M-50_HT-600toInf"))    xsec = 2.21*1.23;

    if(file.Contains("ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola"))  xsec =372.6*1.27;
    if(file.Contains("ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola"))  xsec =100.8*1.27;
    if(file.Contains("ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola"))  xsec =11.99*1.27;
    if(file.Contains("ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola"))  xsec =4.113*1.27;

    if(file.Contains("TTZJets_Tune4C_13TeV-madgraph-tauola"))    xsec = 0.7598;
    if(file.Contains("TTWJets_Tune4C_13TeV-madgraph-tauola"))    xsec = 0.5662;
    // Calculated at 13 TeV in
    // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV
    // Higgs branching ratios from
    // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
    if(file.Contains("ZH_HToBB_ZToLL_M-125_13TeV_powheg-herwigpp"))    xsec = 0.569*0.033658*0.8696;
    if(file.Contains("ZH_HToBB_ZToNuNu_M-125_13TeV_powheg-herwigpp"))    xsec = 0.569*0.2*0.8696;
    if(file.Contains("WH_HToBB_WToLNu_M-125_13TeV_powheg-herwigpp"))    xsec = 0.569*0.1086*1.380;

    if(xsec<=0) std::cout<<"BABYMAKER: Cross section not found for "<<file<<std::endl;

    return xsec;
  }

}
