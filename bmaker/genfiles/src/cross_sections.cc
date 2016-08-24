// Lists the MC cross sections

#include <iostream>
#include "cross_sections.hh"

using namespace std;

namespace xsec{

  float crossSection(const TString &file){
    float xsec(0.);

    if(file.Contains("Run2015"))   xsec = 1.;


    // From https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVgluglu
    if(file.Contains("T1tttt")) xsec = 0.0856418; //just to avoid crash on scans
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

    // cross sections from mcm
    if(file.Contains("TTG")) xsec = 3.697;                
    if(file.Contains("TTTT_TuneCUETP8M1_13TeV-amcatnlo")) xsec = 0.009103;
    if(file.Contains("WJetsToQQ_HT-600ToInf")) xsec = 95.14;
    if(file.Contains("ZJetsToQQ_HT600toInf")) xsec = 5.67;
    
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

    if(file.Contains("QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")) xsec = 27540000;
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
    if(file.Contains("TTZToQQ"))                xsec = 0.5297;
    if(file.Contains("TTZToLLNuNu_M-10"))       xsec = 0.2529;
    if(file.Contains("TTWJetsToQQ"))            xsec = 0.4062;
    if(file.Contains("TTWJetsToLNu"))           xsec = 0.2043;
   
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Diboson
    if(file.Contains("WWTo2L2Nu_13TeV-powheg"))   xsec = 12.178; //NNLO
    if(file.Contains("WWToLNuQQ_13TeV-powheg"))   xsec = 49.997; //NNLO
    if(file.Contains("ttHJetTobb_M125_13TeV_amcatnloFXFX"))   xsec = 0.2934;

    if(file.Contains("WZTo2L2Q_13TeV_amcatnloFXFX"))   xsec = 5.595;
    if(file.Contains("WZTo3LNu_TuneCUETP8M1_13TeV-powheg-"))   xsec = 4.42965;
    if(file.Contains("VVTo2L2Nu_13TeV_amcatnloFXFX"))   xsec = 11.95;

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

  float fractionNegWeights(const TString &file){
    float fneg(0.);

    // ttH, ttZ, ttW, ttGamma
    if(file.Contains("ttHJetTobb_M125_13TeV_amcatnloFXFX"))     fneg = 0.3515;
    if(file.Contains("TTZToQQ"))                                fneg = 0.2657;
    if(file.Contains("TTZToLLNuNu_M-10"))                       fneg = 0.2676;
    if(file.Contains("TTWJetsToQQ"))                            fneg = 0.2412;
    if(file.Contains("TTWJetsToLNu"))                           fneg = 0.2433;
    if(file.Contains("TTG"))                                    fneg = 0.342; // from MCM

    if(file.Contains("TTTT_TuneCUETP8M1_13TeV-amcatnlo"))       fneg = 0.41; // from MCM
    if(file.Contains("VVTo2L2Nu_13TeV_amcatnloFXFX"))       fneg = 0.20; // from MCM

    // Single top
    if (file.Contains("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8")) fneg = 0.1884;

    return fneg;
  }

  void higgsinoCrossSection(int hig_mass, float &xsec, float &xsec_unc) {
        if(hig_mass == 127) { xsec = .5824*.5824*7602.2/1000; xsec_unc = 0.0393921; return;}
        else if(hig_mass == 150) { xsec = .5824*.5824*3832.31/1000; xsec_unc = 0.0413612; return;}
        else if(hig_mass == 175) { xsec = .5824*.5824*2267.94/1000; xsec_unc = 0.044299;  return;}
        else if(hig_mass == 200) { xsec = .5824*.5824*1335.62/1000; xsec_unc = 0.0474362; return;}
        else if(hig_mass == 225) { xsec = .5824*.5824*860.597/1000; xsec_unc = 0.0504217; return;}
        else if(hig_mass == 250) { xsec = .5824*.5824*577.314/1000; xsec_unc = 0.0532731; return;}
        else if(hig_mass == 275) { xsec = .5824*.5824*400.107/1000; xsec_unc = 0.0560232; return;}
        else if(hig_mass == 300) { xsec = .5824*.5824*284.855/1000; xsec_unc = 0.0586867; return;}
        else if(hig_mass == 325) { xsec = .5824*.5824* 207.36/1000; xsec_unc = 0.0613554; return;}
        else if(hig_mass == 350) { xsec = .5824*.5824*153.841/1000; xsec_unc = 0.0640598; return;}
        else if(hig_mass == 375) { xsec = .5824*.5824*116.006/1000; xsec_unc = 0.066892;  return;}
        else if(hig_mass == 400) { xsec = .5824*.5824*88.7325/1000; xsec_unc = 0.0697517; return;}
        else if(hig_mass == 425) { xsec = .5824*.5824*68.6963/1000; xsec_unc = 0.0723531; return;}
        else if(hig_mass == 450) { xsec = .5824*.5824*53.7702/1000; xsec_unc = 0.0748325; return;}
        else if(hig_mass == 475) { xsec = .5824*.5824*42.4699/1000; xsec_unc = 0.0775146; return;}
        else if(hig_mass == 500) { xsec = .5824*.5824*33.8387/1000; xsec_unc = 0.0802572; return;}
        else if(hig_mass == 525) { xsec = .5824*.5824*27.1867/1000; xsec_unc = 0.0825803; return;}
        else if(hig_mass == 550) { xsec = .5824*.5824*21.9868/1000; xsec_unc = 0.0849278; return;}
        else if(hig_mass == 575) { xsec = .5824*.5824*17.9062/1000; xsec_unc = 0.087561;  return;}
        else if(hig_mass == 600) { xsec = .5824*.5824*14.6677/1000; xsec_unc = 0.0900693; return;}
        else if(hig_mass == 625) { xsec = .5824*.5824* 12.062/1000; xsec_unc = 0.091959;  return;}
        else if(hig_mass == 650) { xsec = .5824*.5824*9.96406/1000; xsec_unc = 0.094065;  return;}
        else if(hig_mass == 675) { xsec = .5824*.5824*8.28246/1000; xsec_unc = 0.0957436; return;}
        else if(hig_mass == 700) { xsec = .5824*.5824*6.89981/1000; xsec_unc = 0.0982894; return;}
        else if(hig_mass == 725) { xsec = .5824*.5824*5.78355/1000; xsec_unc = 0.0999915; return;}
        else if(hig_mass == 750) { xsec = .5824*.5824* 4.8731/1000; xsec_unc = 0.101211;  return;}
        else if(hig_mass == 775) { xsec = .5824*.5824*4.09781/1000; xsec_unc = 0.104646;  return;}
        else if(hig_mass == 800) { xsec = .5824*.5824*3.46143/1000; xsec_unc = 0.107618;  return;}
        else if(hig_mass == 825) { xsec = .5824*.5824* 2.9337/1000; xsec_unc = 0.108353;  return;}
        else if(hig_mass == 850) { xsec = .5824*.5824* 2.4923/1000; xsec_unc = 0.110016;  return;}
        else if(hig_mass == 875) { xsec = .5824*.5824*2.13679/1000; xsec_unc = 0.112636;  return;}
        else if(hig_mass == 900) { xsec = .5824*.5824*1.80616/1000; xsec_unc = 0.1134;    return;}
        else if(hig_mass == 925) { xsec = .5824*.5824*1.55453/1000; xsec_unc = 0.116949;  return;}
        else if(hig_mass == 950) { xsec = .5824*.5824*1.32692/1000; xsec_unc = 0.117027;  return;}
        else if(hig_mass == 975) { xsec = .5824*.5824*1.12975/1000; xsec_unc = 0.121244;  return;}
        else if(hig_mass ==1000) { xsec = .5824*.5824*0.968853/1000; xsec_unc = 0.126209; return;}
        else{ xsec = 0; xsec_unc = 0;}
  }

  void signalCrossSection(int glu_mass, float &xsec, float &xsec_unc){
    if (glu_mass == 595) { xsec = 0.; xsec_unc = 0.; return; } // we shouldn't have these points
    else if (glu_mass == 600) { xsec = 9.20353; xsec_unc = 0.137185; return; }
    else if (glu_mass == 605) { xsec = 8.74315; xsec_unc = 0.137502; return; }
    else if (glu_mass == 610) { xsec = 8.30988; xsec_unc = 0.136818; return; }
    else if (glu_mass == 615) { xsec = 7.9012; xsec_unc = 0.137122; return; }
    else if (glu_mass == 620) { xsec = 7.51811; xsec_unc = 0.136123; return; }
    else if (glu_mass == 625) { xsec = 7.15194; xsec_unc = 0.136874; return; }
    else if (glu_mass == 630) { xsec = 6.80558; xsec_unc = 0.136655; return; }
    else if (glu_mass == 635) { xsec = 6.47541; xsec_unc = 0.136459; return; }
    else if (glu_mass == 640) { xsec = 6.17196; xsec_unc = 0.136449; return; }
    else if (glu_mass == 645) { xsec = 5.87366; xsec_unc = 0.136392; return; }
    else if (glu_mass == 650) { xsec = 5.60048; xsec_unc = 0.136262; return; }
    else if (glu_mass == 655) { xsec = 5.33799; xsec_unc = 0.137031; return; }
    else if (glu_mass == 660) { xsec = 5.09822; xsec_unc = 0.137329; return; }
    else if (glu_mass == 665) { xsec = 4.86409; xsec_unc = 0.137702; return; }
    else if (glu_mass == 670) { xsec = 4.64349; xsec_unc = 0.138022; return; }
    else if (glu_mass == 675) { xsec = 4.43132; xsec_unc = 0.138749; return; }
    else if (glu_mass == 680) { xsec = 4.23046; xsec_unc = 0.139166; return; }
    else if (glu_mass == 685) { xsec = 4.03841; xsec_unc = 0.139934; return; }
    else if (glu_mass == 690) { xsec = 3.85666; xsec_unc = 0.139917; return; }
    else if (glu_mass == 695) { xsec = 3.68567; xsec_unc = 0.140759; return; }
    else if (glu_mass == 700) { xsec = 3.5251; xsec_unc = 0.141034; return; }
    else if (glu_mass == 705) { xsec = 3.3737; xsec_unc = 0.141609; return; }
    else if (glu_mass == 710) { xsec = 3.22336; xsec_unc = 0.141972; return; }
    else if (glu_mass == 715) { xsec = 3.0811; xsec_unc = 0.142311; return; }
    else if (glu_mass == 720) { xsec = 2.9509; xsec_unc = 0.142518; return; }
    else if (glu_mass == 725) { xsec = 2.81957; xsec_unc = 0.143333; return; }
    else if (glu_mass == 730) { xsec = 2.7; xsec_unc = 0.143772; return; }
    else if (glu_mass == 735) { xsec = 2.57737; xsec_unc = 0.144452; return; }
    else if (glu_mass == 740) { xsec = 2.47729; xsec_unc = 0.144485; return; }
    else if (glu_mass == 745) { xsec = 2.3661; xsec_unc = 0.145381; return; }
    else if (glu_mass == 750) { xsec = 2.26585; xsec_unc = 0.145653; return; }
    else if (glu_mass == 755) { xsec = 2.17436; xsec_unc = 0.145861; return; }
    else if (glu_mass == 760) { xsec = 2.08446; xsec_unc = 0.146279; return; }
    else if (glu_mass == 765) { xsec = 1.99341; xsec_unc = 0.147278; return; }
    else if (glu_mass == 770) { xsec = 1.91352; xsec_unc = 0.147424; return; }
    else if (glu_mass == 775) { xsec = 1.83188; xsec_unc = 0.147835; return; }
    else if (glu_mass == 780) { xsec = 1.76145; xsec_unc = 0.148078; return; }
    else if (glu_mass == 785) { xsec = 1.68078; xsec_unc = 0.148956; return; }
    else if (glu_mass == 790) { xsec = 1.62071; xsec_unc = 0.149017; return; }
    else if (glu_mass == 795) { xsec = 1.54896; xsec_unc = 0.149976; return; }
    else if (glu_mass == 800) { xsec = 1.4891; xsec_unc = 0.150167; return; }
    else if (glu_mass == 805) { xsec = 1.42888; xsec_unc = 0.150599; return; }
    else if (glu_mass == 810) { xsec = 1.36759; xsec_unc = 0.151122; return; }
    else if (glu_mass == 815) { xsec = 1.31749; xsec_unc = 0.151184; return; }
    else if (glu_mass == 820) { xsec = 1.26659; xsec_unc = 0.151928; return; }
    else if (glu_mass == 825) { xsec = 1.2167; xsec_unc = 0.152141; return; }
    else if (glu_mass == 830) { xsec = 1.16617; xsec_unc = 0.152437; return; }
    else if (glu_mass == 835) { xsec = 1.12555; xsec_unc = 0.153009; return; }
    else if (glu_mass == 840) { xsec = 1.07523; xsec_unc = 0.15367; return; }
    else if (glu_mass == 845) { xsec = 1.03426; xsec_unc = 0.154018; return; }
    else if (glu_mass == 850) { xsec = 0.996137; xsec_unc = 0.154252; return; }
    else if (glu_mass == 855) { xsec = 0.957975; xsec_unc = 0.154597; return; }
    else if (glu_mass == 860) { xsec = 0.921447; xsec_unc = 0.155362; return; }
    else if (glu_mass == 865) { xsec = 0.885917; xsec_unc = 0.155643; return; }
    else if (glu_mass == 870) { xsec = 0.852433; xsec_unc = 0.156368; return; }
    else if (glu_mass == 875) { xsec = 0.820259; xsec_unc = 0.156742; return; }
    else if (glu_mass == 880) { xsec = 0.788789; xsec_unc = 0.156746; return; }
    else if (glu_mass == 885) { xsec = 0.759346; xsec_unc = 0.157507; return; }
    else if (glu_mass == 890) { xsec = 0.731213; xsec_unc = 0.157879; return; }
    else if (glu_mass == 895) { xsec = 0.703532; xsec_unc = 0.158276; return; }
    else if (glu_mass == 900) { xsec = 0.677478; xsec_unc = 0.158762; return; }
    else if (glu_mass == 905) { xsec = 0.652317; xsec_unc = 0.15914; return; }
    else if (glu_mass == 910) { xsec = 0.627695; xsec_unc = 0.159569; return; }
    else if (glu_mass == 915) { xsec = 0.605596; xsec_unc = 0.159838; return; }
    else if (glu_mass == 920) { xsec = 0.58302; xsec_unc = 0.16029; return; }
    else if (glu_mass == 925) { xsec = 0.561889; xsec_unc = 0.160626; return; }
    else if (glu_mass == 930) { xsec = 0.540533; xsec_unc = 0.161499; return; }
    else if (glu_mass == 935) { xsec = 0.521159; xsec_unc = 0.161607; return; }
    else if (glu_mass == 940) { xsec = 0.501865; xsec_unc = 0.16245; return; }
    else if (glu_mass == 945) { xsec = 0.483546; xsec_unc = 0.162492; return; }
    else if (glu_mass == 950) { xsec = 0.466352; xsec_unc = 0.163378; return; }
    else if (glu_mass == 955) { xsec = 0.45012; xsec_unc = 0.163303; return; }
    else if (glu_mass == 960) { xsec = 0.433842; xsec_unc = 0.164161; return; }
    else if (glu_mass == 965) { xsec = 0.418744; xsec_unc = 0.164473; return; }
    else if (glu_mass == 970) { xsec = 0.403514; xsec_unc = 0.164538; return; }
    else if (glu_mass == 975) { xsec = 0.389266; xsec_unc = 0.165308; return; }
    else if (glu_mass == 980) { xsec = 0.375053; xsec_unc = 0.165398; return; }
    else if (glu_mass == 985) { xsec = 0.36182; xsec_unc = 0.16619; return; }
    else if (glu_mass == 990) { xsec = 0.349764; xsec_unc = 0.166462; return; }
    else if (glu_mass == 995) { xsec = 0.337454; xsec_unc = 0.166888; return; }
    else if (glu_mass == 1000) { xsec = 0.325388; xsec_unc = 0.16758; return; }
    else if (glu_mass == 1005) { xsec = 0.314329; xsec_unc = 0.167865; return; }
    else if (glu_mass == 1010) { xsec = 0.30314; xsec_unc = 0.168766; return; }
    else if (glu_mass == 1015) { xsec = 0.292987; xsec_unc = 0.168793; return; }
    else if (glu_mass == 1020) { xsec = 0.282927; xsec_unc = 0.169098; return; }
    else if (glu_mass == 1025) { xsec = 0.272778; xsec_unc = 0.169917; return; }
    else if (glu_mass == 1030) { xsec = 0.263724; xsec_unc = 0.170244; return; }
    else if (glu_mass == 1035) { xsec = 0.254721; xsec_unc = 0.170758; return; }
    else if (glu_mass == 1040) { xsec = 0.245426; xsec_unc = 0.171325; return; }
    else if (glu_mass == 1045) { xsec = 0.237403; xsec_unc = 0.171542; return; }
    else if (glu_mass == 1050) { xsec = 0.229367; xsec_unc = 0.171975; return; }
    else if (glu_mass == 1055) { xsec = 0.221273; xsec_unc = 0.172482; return; }
    else if (glu_mass == 1060) { xsec = 0.214167; xsec_unc = 0.173167; return; }
    else if (glu_mass == 1065) { xsec = 0.207025; xsec_unc = 0.173211; return; }
    else if (glu_mass == 1070) { xsec = 0.199967; xsec_unc = 0.173603; return; }
    else if (glu_mass == 1075) { xsec = 0.193881; xsec_unc = 0.174329; return; }
    else if (glu_mass == 1080) { xsec = 0.186836; xsec_unc = 0.174816; return; }
    else if (glu_mass == 1085) { xsec = 0.180783; xsec_unc = 0.175245; return; }
    else if (glu_mass == 1090) { xsec = 0.174652; xsec_unc = 0.175336; return; }
    else if (glu_mass == 1095) { xsec = 0.168526; xsec_unc = 0.176231; return; }
    else if (glu_mass == 1100) { xsec = 0.163491; xsec_unc = 0.176402; return; }
    else if (glu_mass == 1105) { xsec = 0.158451; xsec_unc = 0.176564; return; }
    else if (glu_mass == 1110) { xsec = 0.153298; xsec_unc = 0.177266; return; }
    else if (glu_mass == 1115) { xsec = 0.148246; xsec_unc = 0.177755; return; }
    else if (glu_mass == 1120) { xsec = 0.143169; xsec_unc = 0.17813; return; }
    else if (glu_mass == 1125) { xsec = 0.139009; xsec_unc = 0.178569; return; }
    else if (glu_mass == 1130) { xsec = 0.133972; xsec_unc = 0.179205; return; }
    else if (glu_mass == 1135) { xsec = 0.129938; xsec_unc = 0.17938; return; }
    else if (glu_mass == 1140) { xsec = 0.125799; xsec_unc = 0.179658; return; }
    else if (glu_mass == 1145) { xsec = 0.121755; xsec_unc = 0.180222; return; }
    else if (glu_mass == 1150) { xsec = 0.117687; xsec_unc = 0.180655; return; }
    else if (glu_mass == 1155) { xsec = 0.11358; xsec_unc = 0.181327; return; }
    else if (glu_mass == 1160) { xsec = 0.110557; xsec_unc = 0.181465; return; }
    else if (glu_mass == 1165) { xsec = 0.107532; xsec_unc = 0.181655; return; }
    else if (glu_mass == 1170) { xsec = 0.10339; xsec_unc = 0.182421; return; }
    else if (glu_mass == 1175) { xsec = 0.10036; xsec_unc = 0.182686; return; }
    else if (glu_mass == 1180) { xsec = 0.0971485; xsec_unc = 0.183142; return; }
    else if (glu_mass == 1185) { xsec = 0.0942072; xsec_unc = 0.183623; return; }
    else if (glu_mass == 1190) { xsec = 0.0912756; xsec_unc = 0.183957; return; }
    else if (glu_mass == 1195) { xsec = 0.0883712; xsec_unc = 0.184467; return; }
    else if (glu_mass == 1200) { xsec = 0.0856418; xsec_unc = 0.184814; return; }
    else if (glu_mass == 1205) { xsec = 0.0830236; xsec_unc = 0.185276; return; }
    else if (glu_mass == 1210) { xsec = 0.0804313; xsec_unc = 0.185714; return; }
    else if (glu_mass == 1215) { xsec = 0.0779039; xsec_unc = 0.186096; return; }
    else if (glu_mass == 1220) { xsec = 0.0755801; xsec_unc = 0.186429; return; }
    else if (glu_mass == 1225) { xsec = 0.0732255; xsec_unc = 0.187227; return; }
    else if (glu_mass == 1230) { xsec = 0.0709683; xsec_unc = 0.187266; return; }
    else if (glu_mass == 1235) { xsec = 0.0688462; xsec_unc = 0.187544; return; }
    else if (glu_mass == 1240) { xsec = 0.0666928; xsec_unc = 0.188404; return; }
    else if (glu_mass == 1245) { xsec = 0.0646423; xsec_unc = 0.188414; return; }
    else if (glu_mass == 1250) { xsec = 0.0627027; xsec_unc = 0.189328; return; }
    else if (glu_mass == 1255) { xsec = 0.0607803; xsec_unc = 0.189693; return; }
    else if (glu_mass == 1260) { xsec = 0.0589319; xsec_unc = 0.189695; return; }
    else if (glu_mass == 1265) { xsec = 0.0571859; xsec_unc = 0.190561; return; }
    else if (glu_mass == 1270) { xsec = 0.0554225; xsec_unc = 0.191806; return; }
    else if (glu_mass == 1275) { xsec = 0.0536906; xsec_unc = 0.192452; return; }
    else if (glu_mass == 1280) { xsec = 0.052051; xsec_unc = 0.192396; return; }
    else if (glu_mass == 1285) { xsec = 0.0504982; xsec_unc = 0.193577; return; }
    else if (glu_mass == 1290) { xsec = 0.0489404; xsec_unc = 0.194903; return; }
    else if (glu_mass == 1295) { xsec = 0.047474; xsec_unc = 0.195871; return; }
    else if (glu_mass == 1300) { xsec = 0.0460525; xsec_unc = 0.1964; return; }
    else if (glu_mass == 1305) { xsec = 0.0447038; xsec_unc = 0.197627; return; }
    else if (glu_mass == 1310) { xsec = 0.0433373; xsec_unc = 0.198601; return; }
    else if (glu_mass == 1315) { xsec = 0.0420362; xsec_unc = 0.198634; return; }
    else if (glu_mass == 1320) { xsec = 0.0407723; xsec_unc = 0.199586; return; }
    else if (glu_mass == 1325) { xsec = 0.0395728; xsec_unc = 0.19951; return; }
    else if (glu_mass == 1330) { xsec = 0.0383587; xsec_unc = 0.19993; return; }
    else if (glu_mass == 1335) { xsec = 0.0372043; xsec_unc = 0.201012; return; }
    else if (glu_mass == 1340) { xsec = 0.0361694; xsec_unc = 0.202191; return; }
    else if (glu_mass == 1345) { xsec = 0.0350586; xsec_unc = 0.201714; return; }
    else if (glu_mass == 1350) { xsec = 0.0340187; xsec_unc = 0.203088; return; }
    else if (glu_mass == 1355) { xsec = 0.0330251; xsec_unc = 0.202807; return; }
    else if (glu_mass == 1360) { xsec = 0.0320787; xsec_unc = 0.203682; return; }
    else if (glu_mass == 1365) { xsec = 0.0311325; xsec_unc = 0.205466; return; }
    else if (glu_mass == 1370) { xsec = 0.0302294; xsec_unc = 0.204724; return; }
    else if (glu_mass == 1375) { xsec = 0.0292919; xsec_unc = 0.206217; return; }
    else if (glu_mass == 1380) { xsec = 0.0284627; xsec_unc = 0.207773; return; }
    else if (glu_mass == 1385) { xsec = 0.0276679; xsec_unc = 0.206729; return; }
    else if (glu_mass == 1390) { xsec = 0.0268339; xsec_unc = 0.208251; return; }
    else if (glu_mass == 1395) { xsec = 0.0260313; xsec_unc = 0.207488; return; }
    else if (glu_mass == 1400) { xsec = 0.0252977; xsec_unc = 0.209163; return; }
    else if (glu_mass == 1405) { xsec = 0.0245679; xsec_unc = 0.210704; return; }
    else if (glu_mass == 1410) { xsec = 0.0238741; xsec_unc = 0.209586; return; }
    else if (glu_mass == 1415) { xsec = 0.0231433; xsec_unc = 0.211204; return; }
    else if (glu_mass == 1420) { xsec = 0.0225194; xsec_unc = 0.212481; return; }
    else if (glu_mass == 1425) { xsec = 0.0218959; xsec_unc = 0.214183; return; }
    else if (glu_mass == 1430) { xsec = 0.0211928; xsec_unc = 0.21365; return; }
    else if (glu_mass == 1435) { xsec = 0.0206244; xsec_unc = 0.217574; return; }
    else if (glu_mass == 1440) { xsec = 0.0200458; xsec_unc = 0.216629; return; }
    else if (glu_mass == 1445) { xsec = 0.0194648; xsec_unc = 0.215531; return; }
    else if (glu_mass == 1450) { xsec = 0.0188887; xsec_unc = 0.219548; return; }
    else if (glu_mass == 1455) { xsec = 0.018364; xsec_unc = 0.221266; return; }
    else if (glu_mass == 1460) { xsec = 0.0178858; xsec_unc = 0.220054; return; }
    else if (glu_mass == 1465) { xsec = 0.0173622; xsec_unc = 0.221916; return; }
    else if (glu_mass == 1470) { xsec = 0.0168403; xsec_unc = 0.223972; return; }
    else if (glu_mass == 1475) { xsec = 0.0163556; xsec_unc = 0.222173; return; }
    else if (glu_mass == 1480) { xsec = 0.0159386; xsec_unc = 0.223581; return; }
    else if (glu_mass == 1485) { xsec = 0.0154568; xsec_unc = 0.222281; return; }
    else if (glu_mass == 1490) { xsec = 0.0150345; xsec_unc = 0.224111; return; }
    else if (glu_mass == 1495) { xsec = 0.0146102; xsec_unc = 0.225293; return; }
    else if (glu_mass == 1500) { xsec = 0.0141903; xsec_unc = 0.227296; return; }
    else if (glu_mass == 1505) { xsec = 0.01377; xsec_unc = 0.229402; return; }
    else if (glu_mass == 1510) { xsec = 0.0133923; xsec_unc = 0.226528; return; }
    else if (glu_mass == 1515) { xsec = 0.0130286; xsec_unc = 0.232697; return; }
    else if (glu_mass == 1520) { xsec = 0.012649; xsec_unc = 0.230194; return; }
    else if (glu_mass == 1525) { xsec = 0.0123374; xsec_unc = 0.231801; return; }
    else if (glu_mass == 1530) { xsec = 0.0119628; xsec_unc = 0.229449; return; }
    else if (glu_mass == 1535) { xsec = 0.0116378; xsec_unc = 0.231293; return; }
    else if (glu_mass == 1540) { xsec = 0.0113183; xsec_unc = 0.233535; return; }
    else if (glu_mass == 1545) { xsec = 0.0110039; xsec_unc = 0.23456; return; }
    else if (glu_mass == 1550) { xsec = 0.0107027; xsec_unc = 0.234971; return; }
    else if (glu_mass == 1555) { xsec = 0.0103967; xsec_unc = 0.23505; return; }
    else if (glu_mass == 1560) { xsec = 0.0101149; xsec_unc = 0.236723; return; }
    else if (glu_mass == 1565) { xsec = 0.00984079; xsec_unc = 0.237486; return; }
    else if (glu_mass == 1570) { xsec = 0.00956216; xsec_unc = 0.238011; return; }
    else if (glu_mass == 1575) { xsec = 0.00930893; xsec_unc = 0.238712; return; }
    else if (glu_mass == 1580) { xsec = 0.00905112; xsec_unc = 0.239145; return; }
    else if (glu_mass == 1585) { xsec = 0.00880102; xsec_unc = 0.24088; return; }
    else if (glu_mass == 1590) { xsec = 0.00856388; xsec_unc = 0.241033; return; }
    else if (glu_mass == 1595) { xsec = 0.00832287; xsec_unc = 0.242052; return; }
    else if (glu_mass == 1600) { xsec = 0.00810078; xsec_unc = 0.242679; return; }
    else if (glu_mass == 1605) { xsec = 0.0078785; xsec_unc = 0.243322; return; }
    else if (glu_mass == 1610) { xsec = 0.00767087; xsec_unc = 0.244839; return; }
    else if (glu_mass == 1615) { xsec = 0.00745579; xsec_unc = 0.245137; return; }
    else if (glu_mass == 1620) { xsec = 0.00725443; xsec_unc = 0.24569; return; }
    else if (glu_mass == 1625) { xsec = 0.00705942; xsec_unc = 0.246853; return; }
    else if (glu_mass == 1630) { xsec = 0.00687457; xsec_unc = 0.24804; return; }
    else if (glu_mass == 1635) { xsec = 0.00668418; xsec_unc = 0.248672; return; }
    else if (glu_mass == 1640) { xsec = 0.00651001; xsec_unc = 0.249776; return; }
    else if (glu_mass == 1645) { xsec = 0.00633268; xsec_unc = 0.250679; return; }
    else if (glu_mass == 1650) { xsec = 0.00616072; xsec_unc = 0.25138; return; }
    else if (glu_mass == 1655) { xsec = 0.00599673; xsec_unc = 0.252591; return; }
    else if (glu_mass == 1660) { xsec = 0.00583243; xsec_unc = 0.253829; return; }
    else if (glu_mass == 1665) { xsec = 0.00567868; xsec_unc = 0.255006; return; }
    else if (glu_mass == 1670) { xsec = 0.00553066; xsec_unc = 0.255203; return; }
    else if (glu_mass == 1675) { xsec = 0.00538094; xsec_unc = 0.255439; return; }
    else if (glu_mass == 1680) { xsec = 0.00523764; xsec_unc = 0.256602; return; }
    else if (glu_mass == 1685) { xsec = 0.00509647; xsec_unc = 0.258745; return; }
    else if (glu_mass == 1690) { xsec = 0.0049577; xsec_unc = 0.258847; return; }
    else if (glu_mass == 1695) { xsec = 0.00483094; xsec_unc = 0.260944; return; }
    else if (glu_mass == 1700) { xsec = 0.00470323; xsec_unc = 0.261021; return; }
    else if (glu_mass == 1705) { xsec = 0.0045807; xsec_unc = 0.262095; return; }
    else if (glu_mass == 1710) { xsec = 0.00445824; xsec_unc = 0.263238; return; }
    else if (glu_mass == 1715) { xsec = 0.0043369; xsec_unc = 0.263092; return; }
    else if (glu_mass == 1720) { xsec = 0.00422488; xsec_unc = 0.264093; return; }
    else if (glu_mass == 1725) { xsec = 0.00411276; xsec_unc = 0.26513; return; }
    else if (glu_mass == 1730) { xsec = 0.00400698; xsec_unc = 0.267386; return; }
    else if (glu_mass == 1735) { xsec = 0.00389655; xsec_unc = 0.267109; return; }
    else if (glu_mass == 1740) { xsec = 0.00379497; xsec_unc = 0.268072; return; }
    else if (glu_mass == 1745) { xsec = 0.00370003; xsec_unc = 0.2704; return; }
    else if (glu_mass == 1750) { xsec = 0.00359842; xsec_unc = 0.271502; return; }
    else if (glu_mass == 1755) { xsec = 0.00350486; xsec_unc = 0.27229; return; }
    else if (glu_mass == 1760) { xsec = 0.00341375; xsec_unc = 0.273209; return; }
    else if (glu_mass == 1765) { xsec = 0.00332255; xsec_unc = 0.27416; return; }
    else if (glu_mass == 1770) { xsec = 0.00323809; xsec_unc = 0.276458; return; }
    else if (glu_mass == 1775) { xsec = 0.00314866; xsec_unc = 0.275834; return; }
    else if (glu_mass == 1780) { xsec = 0.00306841; xsec_unc = 0.276481; return; }
    else if (glu_mass == 1785) { xsec = 0.00298808; xsec_unc = 0.277145; return; }
    else if (glu_mass == 1790) { xsec = 0.00291365; xsec_unc = 0.279548; return; }
    else if (glu_mass == 1795) { xsec = 0.0028312; xsec_unc = 0.280642; return; }
    else if (glu_mass == 1800) { xsec = 0.00276133; xsec_unc = 0.28108; return; }
    else if (glu_mass == 1805) { xsec = 0.00269156; xsec_unc = 0.281566; return; }
    else if (glu_mass == 1810) { xsec = 0.00262156; xsec_unc = 0.282017; return; }
    else if (glu_mass == 1815) { xsec = 0.00254938; xsec_unc = 0.282755; return; }
    else if (glu_mass == 1820) { xsec = 0.00248581; xsec_unc = 0.285102; return; }
    else if (glu_mass == 1825) { xsec = 0.00241549; xsec_unc = 0.285869; return; }
    else if (glu_mass == 1830) { xsec = 0.00235625; xsec_unc = 0.286103; return; }
    else if (glu_mass == 1835) { xsec = 0.00229576; xsec_unc = 0.28596; return; }
    else if (glu_mass == 1840) { xsec = 0.00223603; xsec_unc = 0.286654; return; }
    else if (glu_mass == 1845) { xsec = 0.00218302; xsec_unc = 0.288949; return; }
    else if (glu_mass == 1850) { xsec = 0.00212345; xsec_unc = 0.289167; return; }
    else if (glu_mass == 1855) { xsec = 0.00207; xsec_unc = 0.291835; return; }
    else if (glu_mass == 1860) { xsec = 0.00200972; xsec_unc = 0.291901; return; }
    else if (glu_mass == 1865) { xsec = 0.00196025; xsec_unc = 0.292103; return; }
    else if (glu_mass == 1870) { xsec = 0.00191132; xsec_unc = 0.291893; return; }
    else if (glu_mass == 1875) { xsec = 0.00185789; xsec_unc = 0.294928; return; }
    else if (glu_mass == 1880) { xsec = 0.00181527; xsec_unc = 0.29723; return; }
    else if (glu_mass == 1885) { xsec = 0.00176658; xsec_unc = 0.297236; return; }
    else if (glu_mass == 1890) { xsec = 0.00172274; xsec_unc = 0.299813; return; }
    else if (glu_mass == 1895) { xsec = 0.00167806; xsec_unc = 0.296455; return; }
    else if (glu_mass == 1900) { xsec = 0.00163547; xsec_unc = 0.299045; return; }
    else if (glu_mass == 1905) { xsec = 0.0015925; xsec_unc = 0.302039; return; }
    else if (glu_mass == 1910) { xsec = 0.00155445; xsec_unc = 0.301015; return; }
    else if (glu_mass == 1915) { xsec = 0.00151503; xsec_unc = 0.300356; return; }
    else if (glu_mass == 1920) { xsec = 0.00147199; xsec_unc = 0.303575; return; }
    else if (glu_mass == 1925) { xsec = 0.0014401; xsec_unc = 0.305951; return; }
    else if (glu_mass == 1930) { xsec = 0.0014016; xsec_unc = 0.305171; return; }
    else if (glu_mass == 1935) { xsec = 0.00136297; xsec_unc = 0.304873; return; }
    else if (glu_mass == 1940) { xsec = 0.001331; xsec_unc = 0.307414; return; }
    else if (glu_mass == 1945) { xsec = 0.001299; xsec_unc = 0.310066; return; }
    else if (glu_mass == 1950) { xsec = 0.0012642; xsec_unc = 0.304581; return; }
    else if (glu_mass == 1955) { xsec = 0.00123087; xsec_unc = 0.308644; return; }
    else if (glu_mass == 1960) { xsec = 0.00120048; xsec_unc = 0.309669; return; }
    else if (glu_mass == 1965) { xsec = 0.00117053; xsec_unc = 0.310216; return; }
    else if (glu_mass == 1970) { xsec = 0.00114051; xsec_unc = 0.310814; return; }
    else if (glu_mass == 1975) { xsec = 0.00111722; xsec_unc = 0.315357; return; }
    else if (glu_mass == 1980) { xsec = 0.00108758; xsec_unc = 0.315568; return; }
    else if (glu_mass == 1985) { xsec = 0.00105813; xsec_unc = 0.315103; return; }
    else if (glu_mass == 1990) { xsec = 0.00102936; xsec_unc = 0.314167; return; }
    else if (glu_mass == 1995) { xsec = 0.00100614; xsec_unc = 0.317628; return; }
    else if (glu_mass == 2000) { xsec = 0.000981077; xsec_unc = 0.318422; return; }
    else {xsec = 0.; xsec_unc = 0.;} 
  }

  void stopCrossSection(int stop_mass, float &xsec, float &xsec_unc){
    
    if(stop_mass == 100) { xsec = 1521.11; xsec_unc = 0.154038; return; }
    else if(stop_mass == 105) { xsec = 1233.18; xsec_unc = 0.154059; return; }
    else if(stop_mass == 110) { xsec = 1013.76; xsec_unc = 0.154088; return; }
    else if(stop_mass == 115) { xsec = 832.656; xsec_unc = 0.151503; return; }
    else if(stop_mass == 120) { xsec = 689.799; xsec_unc = 0.15044; return; }
    else if(stop_mass == 125) { xsec = 574.981; xsec_unc = 0.149895; return; }
    else if(stop_mass == 130) { xsec = 481.397; xsec_unc = 0.148906; return; }
    else if(stop_mass == 135) { xsec = 405.159; xsec_unc = 0.148952; return; }
    else if(stop_mass == 140) { xsec = 342.865; xsec_unc = 0.149119; return; }
    else if(stop_mass == 145) { xsec = 291.752; xsec_unc = 0.148022; return; }
    else if(stop_mass == 150) { xsec = 249.409; xsec_unc = 0.147477; return; }
    else if(stop_mass == 155) { xsec = 214.221; xsec_unc = 0.145928; return; }
    else if(stop_mass == 160) { xsec = 184.623; xsec_unc = 0.145821; return; }
    else if(stop_mass == 165) { xsec = 159.614; xsec_unc = 0.147859; return; }
    else if(stop_mass == 170) { xsec = 139.252; xsec_unc = 0.14547; return; }
    else if(stop_mass == 175) { xsec = 121.416; xsec_unc = 0.146341; return; }
    else if(stop_mass == 180) { xsec = 106.194; xsec_unc = 0.142033; return; }
    else if(stop_mass == 185) { xsec = 93.3347; xsec_unc = 0.144893; return; }
    else if(stop_mass == 190) { xsec = 82.2541; xsec_unc = 0.144677; return; }
    else if(stop_mass == 195) { xsec = 72.7397; xsec_unc = 0.144452; return; }
    else if(stop_mass == 200) { xsec = 64.5085; xsec_unc = 0.144098; return; }
    else if(stop_mass == 205) { xsec = 57.2279; xsec_unc = 0.144191; return; }
    else if(stop_mass == 210) { xsec = 50.9226; xsec_unc = 0.142457; return; }
    else if(stop_mass == 215) { xsec = 45.3761; xsec_unc = 0.14344; return; }
    else if(stop_mass == 220) { xsec = 40.5941; xsec_unc = 0.142634; return; }
    else if(stop_mass == 225) { xsec = 36.3818; xsec_unc = 0.142189; return; }
    else if(stop_mass == 230) { xsec = 32.6679; xsec_unc = 0.141592; return; }
    else if(stop_mass == 235) { xsec = 29.3155; xsec_unc = 0.142233; return; }
    else if(stop_mass == 240) { xsec = 26.4761; xsec_unc = 0.141723; return; }
    else if(stop_mass == 245) { xsec = 23.8853; xsec_unc = 0.139482; return; }
    else if(stop_mass == 250) { xsec = 21.5949; xsec_unc = 0.140595; return; }
    else if(stop_mass == 255) { xsec = 19.5614; xsec_unc = 0.138755; return; }
    else if(stop_mass == 260) { xsec = 17.6836; xsec_unc = 0.139505; return; }
    else if(stop_mass == 265) { xsec = 16.112; xsec_unc = 0.139531; return; }
    else if(stop_mass == 270) { xsec = 14.6459; xsec_unc = 0.139278; return; }
    else if(stop_mass == 275) { xsec = 13.3231; xsec_unc = 0.142549; return; }
    else if(stop_mass == 280) { xsec = 12.1575; xsec_unc = 0.141584; return; }
    else if(stop_mass == 285) { xsec = 11.0925; xsec_unc = 0.140904; return; }
    else if(stop_mass == 290) { xsec = 10.1363; xsec_unc = 0.138967; return; }
    else if(stop_mass == 295) { xsec = 9.29002; xsec_unc = 0.139107; return; }
    else if(stop_mass == 300) { xsec = 8.51615; xsec_unc = 0.139223; return; }
    else if(stop_mass == 305) { xsec = 7.81428; xsec_unc = 0.138996; return; }
    else if(stop_mass == 310) { xsec = 7.17876; xsec_unc = 0.139357; return; }
    else if(stop_mass == 315) { xsec = 6.60266; xsec_unc = 0.139256; return; }
    else if(stop_mass == 320) { xsec = 6.08444; xsec_unc = 0.137957; return; }
    else if(stop_mass == 325) { xsec = 5.60471; xsec_unc = 0.138144; return; }
    else if(stop_mass == 330) { xsec = 5.17188; xsec_unc = 0.136954; return; }
    else if(stop_mass == 335) { xsec = 4.77871; xsec_unc = 0.137554; return; }
    else if(stop_mass == 340) { xsec = 4.41629; xsec_unc = 0.137945; return; }
    else if(stop_mass == 345) { xsec = 4.08881; xsec_unc = 0.137075; return; }
    else if(stop_mass == 350) { xsec = 3.78661; xsec_unc = 0.136877; return; }
    else if(stop_mass == 355) { xsec = 3.50911; xsec_unc = 0.138089; return; }
    else if(stop_mass == 360) { xsec = 3.25619; xsec_unc = 0.138002; return; }
    else if(stop_mass == 365) { xsec = 3.02472; xsec_unc = 0.137093; return; }
    else if(stop_mass == 370) { xsec = 2.8077; xsec_unc = 0.138064; return; }
    else if(stop_mass == 375) { xsec = 2.61162; xsec_unc = 0.138477; return; }
    else if(stop_mass == 380) { xsec = 2.43031; xsec_unc = 0.136999; return; }
    else if(stop_mass == 385) { xsec = 2.26365; xsec_unc = 0.13728; return; }
    else if(stop_mass == 390) { xsec = 2.10786; xsec_unc = 0.13732; return; }
    else if(stop_mass == 395) { xsec = 1.9665; xsec_unc = 0.134737; return; }
    else if(stop_mass == 400) { xsec = 1.83537; xsec_unc = 0.136985; return; }
    else if(stop_mass == 405) { xsec = 1.70927; xsec_unc = 0.137114; return; }
    else if(stop_mass == 410) { xsec = 1.60378; xsec_unc = 0.135468; return; }
    else if(stop_mass == 415) { xsec = 1.49798; xsec_unc = 0.134453; return; }
    else if(stop_mass == 420) { xsec = 1.39688; xsec_unc = 0.136719; return; }
    else if(stop_mass == 425) { xsec = 1.31169; xsec_unc = 0.135013; return; }
    else if(stop_mass == 430) { xsec = 1.22589; xsec_unc = 0.133237; return; }
    else if(stop_mass == 435) { xsec = 1.14553; xsec_unc = 0.135478; return; }
    else if(stop_mass == 440) { xsec = 1.07484; xsec_unc = 0.137238; return; }
    else if(stop_mass == 445) { xsec = 1.01019; xsec_unc = 0.134187; return; }
    else if(stop_mass == 450) { xsec = 0.948333; xsec_unc = 0.134559; return; }
    else if(stop_mass == 455) { xsec = 0.890847; xsec_unc = 0.134587; return; }
    else if(stop_mass == 460) { xsec = 0.836762; xsec_unc = 0.134468; return; }
    else if(stop_mass == 465) { xsec = 0.787221; xsec_unc = 0.134149; return; }
    else if(stop_mass == 470) { xsec = 0.740549; xsec_unc = 0.134127; return; }
    else if(stop_mass == 475) { xsec = 0.697075; xsec_unc = 0.133926; return; }
    else if(stop_mass == 480) { xsec = 0.655954; xsec_unc = 0.134392; return; }
    else if(stop_mass == 485) { xsec = 0.618562; xsec_unc = 0.133705; return; }
    else if(stop_mass == 490) { xsec = 0.582467; xsec_unc = 0.133914; return; }
    else if(stop_mass == 495) { xsec = 0.549524; xsec_unc = 0.133691; return; }
    else if(stop_mass == 500) { xsec = 0.51848; xsec_unc = 0.133797; return; }
    else if(stop_mass == 505) { xsec = 0.489324; xsec_unc = 0.133608; return; }
    else if(stop_mass == 510) { xsec = 0.462439; xsec_unc = 0.133046; return; }
    else if(stop_mass == 515) { xsec = 0.436832; xsec_unc = 0.133703; return; }
    else if(stop_mass == 520) { xsec = 0.412828; xsec_unc = 0.13272; return; }
    else if(stop_mass == 525) { xsec = 0.390303; xsec_unc = 0.133443; return; }
    else if(stop_mass == 530) { xsec = 0.368755; xsec_unc = 0.133769; return; }
    else if(stop_mass == 535) { xsec = 0.348705; xsec_unc = 0.132706; return; }
    else if(stop_mass == 540) { xsec = 0.330157; xsec_unc = 0.132981; return; }
    else if(stop_mass == 545) { xsec = 0.312672; xsec_unc = 0.13277; return; }
    else if(stop_mass == 550) { xsec = 0.296128; xsec_unc = 0.132687; return; }
    else if(stop_mass == 555) { xsec = 0.280734; xsec_unc = 0.132363; return; }
    else if(stop_mass == 560) { xsec = 0.266138; xsec_unc = 0.13193; return; }
    else if(stop_mass == 565) { xsec = 0.251557; xsec_unc = 0.131731; return; }
    else if(stop_mass == 570) { xsec = 0.238537; xsec_unc = 0.133409; return; }
    else if(stop_mass == 575) { xsec = 0.226118; xsec_unc = 0.132741; return; }
    else if(stop_mass == 580) { xsec = 0.214557; xsec_unc = 0.131697; return; }
    else if(stop_mass == 585) { xsec = 0.203566; xsec_unc = 0.133257; return; }
    else if(stop_mass == 590) { xsec = 0.193079; xsec_unc = 0.132037; return; }
    else if(stop_mass == 595) { xsec = 0.183604; xsec_unc = 0.130973; return; }
    else if(stop_mass == 600) { xsec = 0.174599; xsec_unc = 0.132074; return; }
    else if(stop_mass == 605) { xsec = 0.166131; xsec_unc = 0.130154; return; }
    else if(stop_mass == 610) { xsec = 0.158242; xsec_unc = 0.13142; return; }
    else if(stop_mass == 615) { xsec = 0.150275; xsec_unc = 0.13285; return; }
    else if(stop_mass == 620) { xsec = 0.142787; xsec_unc = 0.130642; return; }
    else if(stop_mass == 625) { xsec = 0.136372; xsec_unc = 0.127962; return; }
    else if(stop_mass == 630) { xsec = 0.129886; xsec_unc = 0.132957; return; }
    else if(stop_mass == 635) { xsec = 0.123402; xsec_unc = 0.13016; return; }
    else if(stop_mass == 640) { xsec = 0.11795; xsec_unc = 0.127132; return; }
    else if(stop_mass == 645) { xsec = 0.112008; xsec_unc = 0.12808; return; }
    else if(stop_mass == 650) { xsec = 0.107045; xsec_unc = 0.129232; return; }
    else if(stop_mass == 655) { xsec = 0.102081; xsec_unc = 0.130012; return; }
    else if(stop_mass == 660) { xsec = 0.09725; xsec_unc = 0.129038; return; }
    else if(stop_mass == 665) { xsec = 0.0927515; xsec_unc = 0.129548; return; }
    else if(stop_mass == 670) { xsec = 0.0885084; xsec_unc = 0.130218; return; }
    else if(stop_mass == 675) { xsec = 0.0844877; xsec_unc = 0.130703; return; }
    else if(stop_mass == 680) { xsec = 0.0806192; xsec_unc = 0.131131; return; }
    else if(stop_mass == 685) { xsec = 0.0769099; xsec_unc = 0.131517; return; }
    else if(stop_mass == 690) { xsec = 0.0734901; xsec_unc = 0.132344; return; }
    else if(stop_mass == 695) { xsec = 0.0701805; xsec_unc = 0.132716; return; }
    else if(stop_mass == 700) { xsec = 0.0670476; xsec_unc = 0.133429; return; }
    else if(stop_mass == 705) { xsec = 0.0641426; xsec_unc = 0.13363; return; }
    else if(stop_mass == 710) { xsec = 0.0612942; xsec_unc = 0.133941; return; }
    else if(stop_mass == 715) { xsec = 0.0585678; xsec_unc = 0.134663; return; }
    else if(stop_mass == 720) { xsec = 0.0560753; xsec_unc = 0.134984; return; }
    else if(stop_mass == 725) { xsec = 0.0536438; xsec_unc = 0.135804; return; }
    else if(stop_mass == 730) { xsec = 0.0513219; xsec_unc = 0.135682; return; }
    else if(stop_mass == 735) { xsec = 0.0491001; xsec_unc = 0.136268; return; }
    else if(stop_mass == 740) { xsec = 0.0470801; xsec_unc = 0.136895; return; }
    else if(stop_mass == 745) { xsec = 0.045061; xsec_unc = 0.136816; return; }
    else if(stop_mass == 750) { xsec = 0.0431418; xsec_unc = 0.137455; return; }
    else if(stop_mass == 755) { xsec = 0.0413447; xsec_unc = 0.137833; return; }
    else if(stop_mass == 760) { xsec = 0.0396264; xsec_unc = 0.138518; return; }
    else if(stop_mass == 765) { xsec = 0.0379036; xsec_unc = 0.138537; return; }
    else if(stop_mass == 770) { xsec = 0.0363856; xsec_unc = 0.139334; return; }
    else if(stop_mass == 775) { xsec = 0.0348796; xsec_unc = 0.139597; return; }
    else if(stop_mass == 780) { xsec = 0.0334669; xsec_unc = 0.140267; return; }
    else if(stop_mass == 785) { xsec = 0.0320548; xsec_unc = 0.140406; return; }
    else if(stop_mass == 790) { xsec = 0.0307373; xsec_unc = 0.14115; return; }
    else if(stop_mass == 795) { xsec = 0.0295348; xsec_unc = 0.141397; return; }
    else if(stop_mass == 800) { xsec = 0.0283338; xsec_unc = 0.14171; return; }
    else if(stop_mass == 805) { xsec = 0.0272206; xsec_unc = 0.14241; return; }
    else if(stop_mass == 810) { xsec = 0.0261233; xsec_unc = 0.142891; return; }
    else if(stop_mass == 815) { xsec = 0.0251107; xsec_unc = 0.143632; return; }
    else if(stop_mass == 820) { xsec = 0.0241099; xsec_unc = 0.143805; return; }
    else if(stop_mass == 825) { xsec = 0.0230866; xsec_unc = 0.144428; return; }
    else if(stop_mass == 830) { xsec = 0.0221834; xsec_unc = 0.144791; return; }
    else if(stop_mass == 835) { xsec = 0.0213766; xsec_unc = 0.145511; return; }
    else if(stop_mass == 840) { xsec = 0.0204715; xsec_unc = 0.146131; return; }
    else if(stop_mass == 845) { xsec = 0.0197653; xsec_unc = 0.146602; return; }
    else if(stop_mass == 850) { xsec = 0.0189612; xsec_unc = 0.14702; return; }
    else if(stop_mass == 855) { xsec = 0.0182516; xsec_unc = 0.147648; return; }
    else if(stop_mass == 860) { xsec = 0.0175509; xsec_unc = 0.147944; return; }
    else if(stop_mass == 865) { xsec = 0.0168336; xsec_unc = 0.148528; return; }
    else if(stop_mass == 870) { xsec = 0.0162314; xsec_unc = 0.148772; return; }
    else if(stop_mass == 875) { xsec = 0.015625; xsec_unc = 0.149567; return; }
    else if(stop_mass == 880) { xsec = 0.0150143; xsec_unc = 0.150389; return; }
    else if(stop_mass == 885) { xsec = 0.0144112; xsec_unc = 0.150614; return; }
    else if(stop_mass == 890) { xsec = 0.0138979; xsec_unc = 0.151; return; }
    else if(stop_mass == 895) { xsec = 0.0133962; xsec_unc = 0.151325; return; }
    else if(stop_mass == 900) { xsec = 0.0128895; xsec_unc = 0.152026; return; }
    else if(stop_mass == 905) { xsec = 0.0123843; xsec_unc = 0.152968; return; }
    else if(stop_mass == 910) { xsec = 0.0119837; xsec_unc = 0.153089; return; }
    else if(stop_mass == 915) { xsec = 0.0114713; xsec_unc = 0.153678; return; }
    else if(stop_mass == 920) { xsec = 0.0110688; xsec_unc = 0.154082; return; }
    else if(stop_mass == 925) { xsec = 0.0106631; xsec_unc = 0.154806; return; }
    else if(stop_mass == 930) { xsec = 0.0102629; xsec_unc = 0.155313; return; }
    else if(stop_mass == 935) { xsec = 0.0098874; xsec_unc = 0.156066; return; }
    else if(stop_mass == 940) { xsec = 0.00952142; xsec_unc = 0.156055; return; }
    else if(stop_mass == 945) { xsec = 0.00916636; xsec_unc = 0.156849; return; }
    else if(stop_mass == 950) { xsec = 0.00883465; xsec_unc = 0.157177; return; }
    else if(stop_mass == 955) { xsec = 0.00851073; xsec_unc = 0.158094; return; }
    else if(stop_mass == 960) { xsec = 0.00820884; xsec_unc = 0.15844; return; }
    else if(stop_mass == 965) { xsec = 0.00791403; xsec_unc = 0.159216; return; }
    else if(stop_mass == 970) { xsec = 0.00763112; xsec_unc = 0.159742; return; }
    else if(stop_mass == 975) { xsec = 0.00735655; xsec_unc = 0.160548; return; }
    else if(stop_mass == 980) { xsec = 0.00710317; xsec_unc = 0.160626; return; }
    else if(stop_mass == 985) { xsec = 0.00684867; xsec_unc = 0.16144; return; }
    else if(stop_mass == 990) { xsec = 0.00660695; xsec_unc = 0.161813; return; }
    else if(stop_mass == 995) { xsec = 0.00637546; xsec_unc = 0.162158; return; }
    else if(stop_mass == 1000) { xsec = 0.00615134; xsec_unc = 0.162953; return; }
    else if(stop_mass == 1005) { xsec = 0.00593765; xsec_unc = 0.163716; return; }
    else if(stop_mass == 1010) { xsec = 0.00572452; xsec_unc = 0.163857; return; }
    else if(stop_mass == 1015) { xsec = 0.00553094; xsec_unc = 0.164628; return; }
    else if(stop_mass == 1020) { xsec = 0.00533968; xsec_unc = 0.164963; return; }
    else if(stop_mass == 1025) { xsec = 0.00514619; xsec_unc = 0.165762; return; }
    else if(stop_mass == 1030) { xsec = 0.00497235; xsec_unc = 0.165838; return; }
    else if(stop_mass == 1035) { xsec = 0.00479906; xsec_unc = 0.166646; return; }
    else if(stop_mass == 1040) { xsec = 0.00463806; xsec_unc = 0.166947; return; }
    else if(stop_mass == 1045) { xsec = 0.00447537; xsec_unc = 0.167071; return; }
    else if(stop_mass == 1050) { xsec = 0.00432261; xsec_unc = 0.167859; return; }
    else if(stop_mass == 1055) { xsec = 0.00417983; xsec_unc = 0.168637; return; }
    else if(stop_mass == 1060) { xsec = 0.00403886; xsec_unc = 0.168981; return; }
    else if(stop_mass == 1065) { xsec = 0.0038962; xsec_unc = 0.169794; return; }
    else if(stop_mass == 1070) { xsec = 0.00376343; xsec_unc = 0.169764; return; }
    else if(stop_mass == 1075) { xsec = 0.00364174; xsec_unc = 0.170634; return; }
    else if(stop_mass == 1080) { xsec = 0.00352093; xsec_unc = 0.170908; return; }
    else if(stop_mass == 1085) { xsec = 0.00339813; xsec_unc = 0.171929; return; }
    else if(stop_mass == 1090) { xsec = 0.00328695; xsec_unc = 0.172274; return; }
    else if(stop_mass == 1095) { xsec = 0.00317628; xsec_unc = 0.172617; return; }
    else if(stop_mass == 1100) { xsec = 0.00307413; xsec_unc = 0.173377; return; }
    else if(stop_mass == 1105) { xsec = 0.00297377; xsec_unc = 0.173822; return; }
    else if(stop_mass == 1110) { xsec = 0.00287148; xsec_unc = 0.174725; return; }
    else if(stop_mass == 1115) { xsec = 0.00278078; xsec_unc = 0.175091; return; }
    else if(stop_mass == 1120) { xsec = 0.00268873; xsec_unc = 0.175883; return; }
    else if(stop_mass == 1125) { xsec = 0.00260821; xsec_unc = 0.176126; return; }
    else if(stop_mass == 1130) { xsec = 0.00251529; xsec_unc = 0.176836; return; }
    else if(stop_mass == 1135) { xsec = 0.00243484; xsec_unc = 0.177128; return; }
    else if(stop_mass == 1140) { xsec = 0.00236295; xsec_unc = 0.177977; return; }
    else if(stop_mass == 1145) { xsec = 0.00228192; xsec_unc = 0.178507; return; }
    else if(stop_mass == 1150) { xsec = 0.00221047; xsec_unc = 0.179259; return; }
    else if(stop_mass == 1155) { xsec = 0.00213907; xsec_unc = 0.180255; return; }
    else if(stop_mass == 1160) { xsec = 0.00206845; xsec_unc = 0.180518; return; }
    else if(stop_mass == 1165) { xsec = 0.0020063; xsec_unc = 0.180954; return; }
    else if(stop_mass == 1170) { xsec = 0.00194569; xsec_unc = 0.181194; return; }
    else if(stop_mass == 1175) { xsec = 0.0018741; xsec_unc = 0.182145; return; }
    else if(stop_mass == 1180) { xsec = 0.00182266; xsec_unc = 0.183074; return; }
    else if(stop_mass == 1185) { xsec = 0.00176211; xsec_unc = 0.183375; return; }
    else if(stop_mass == 1190) { xsec = 0.00170006; xsec_unc = 0.184075; return; }
    else if(stop_mass == 1195) { xsec = 0.00164968; xsec_unc = 0.184438; return; }
    else if(stop_mass == 1200) { xsec = 0.00159844; xsec_unc = 0.185209; return; }
    else if(stop_mass == 1205) { xsec = 0.0015472; xsec_unc = 0.185977; return; }
    else if(stop_mass == 1210) { xsec = 0.00149657; xsec_unc = 0.186485; return; }
    else if(stop_mass == 1215) { xsec = 0.00145544; xsec_unc = 0.187347; return; }
    else if(stop_mass == 1220) { xsec = 0.00140288; xsec_unc = 0.188774; return; }
    else if(stop_mass == 1225) { xsec = 0.00136155; xsec_unc = 0.18989; return; }
    else if(stop_mass == 1230) { xsec = 0.00131271; xsec_unc = 0.188763; return; }
    else if(stop_mass == 1235) { xsec = 0.0012717; xsec_unc = 0.189588; return; }
    else if(stop_mass == 1240) { xsec = 0.00123066; xsec_unc = 0.19049; return; }
    else if(stop_mass == 1245) { xsec = 0.00119994; xsec_unc = 0.191442; return; }
    else if(stop_mass == 1250) { xsec = 0.0011583; xsec_unc = 0.193006; return; }
    else if(stop_mass == 1255) { xsec = 0.00112694; xsec_unc = 0.194441; return; }
    else if(stop_mass == 1260) { xsec = 0.00108716; xsec_unc = 0.194141; return; }
    else if(stop_mass == 1265) { xsec = 0.00105517; xsec_unc = 0.196361; return; }
    else if(stop_mass == 1270) { xsec = 0.00102241; xsec_unc = 0.196297; return; }
    else if(stop_mass == 1275) { xsec = 0.000991293; xsec_unc = 0.19762; return; }
    else if(stop_mass == 1280) { xsec = 0.000961012; xsec_unc = 0.197926; return; }
    else if(stop_mass == 1285) { xsec = 0.000932394; xsec_unc = 0.198682; return; }
    else if(stop_mass == 1290) { xsec = 0.000903404; xsec_unc = 0.199924; return; }
    else if(stop_mass == 1295) { xsec = 0.000876957; xsec_unc = 0.200777; return; }
    else if(stop_mass == 1300) { xsec = 0.000850345; xsec_unc = 0.201604; return; }
    else if(stop_mass == 1305) { xsec = 0.00082443; xsec_unc = 0.202883; return; }
    else if(stop_mass == 1310) { xsec = 0.00079983; xsec_unc = 0.20373; return; }
    else if(stop_mass == 1315) { xsec = 0.000775222; xsec_unc = 0.204622; return; }
    else if(stop_mass == 1320) { xsec = 0.000751372; xsec_unc = 0.205919; return; }
    else if(stop_mass == 1325) { xsec = 0.000728912; xsec_unc = 0.206884; return; }
    else if(stop_mass == 1330) { xsec = 0.000706867; xsec_unc = 0.207763; return; }
    else if(stop_mass == 1335) { xsec = 0.000685372; xsec_unc = 0.208587; return; }
    else if(stop_mass == 1340) { xsec = 0.000664649; xsec_unc = 0.209879; return; }
    else if(stop_mass == 1345) { xsec = 0.000644804; xsec_unc = 0.211487; return; }
    else if(stop_mass == 1350) { xsec = 0.000625155; xsec_unc = 0.212761; return; }
    else if(stop_mass == 1355) { xsec = 0.000606802; xsec_unc = 0.213529; return; }
    else if(stop_mass == 1360) { xsec = 0.000588512; xsec_unc = 0.214428; return; }
    else if(stop_mass == 1365) { xsec = 0.000570506; xsec_unc = 0.216584; return; }
    else if(stop_mass == 1370) { xsec = 0.000553379; xsec_unc = 0.216036; return; }
    else if(stop_mass == 1375) { xsec = 0.000536646; xsec_unc = 0.21775; return; }
    else if(stop_mass == 1380) { xsec = 0.000521404; xsec_unc = 0.218383; return; }
    else if(stop_mass == 1385) { xsec = 0.000505008; xsec_unc = 0.219675; return; }
    else if(stop_mass == 1390) { xsec = 0.000490353; xsec_unc = 0.221444; return; }
    else if(stop_mass == 1395) { xsec = 0.000476164; xsec_unc = 0.222016; return; }
    else if(stop_mass == 1400) { xsec = 0.000461944; xsec_unc = 0.222704; return; }
    else if(stop_mass == 1405) { xsec = 0.000448172; xsec_unc = 0.224911; return; }
    else if(stop_mass == 1410) { xsec = 0.000435082; xsec_unc = 0.225606; return; }
    else if(stop_mass == 1415) { xsec = 0.000422967; xsec_unc = 0.226095; return; }
    else if(stop_mass == 1420) { xsec = 0.000410381; xsec_unc = 0.22797; return; }
    else if(stop_mass == 1425) { xsec = 0.000398106; xsec_unc = 0.228949; return; }
    else if(stop_mass == 1430) { xsec = 0.000386792; xsec_unc = 0.231319; return; }
    else if(stop_mass == 1435) { xsec = 0.000375724; xsec_unc = 0.231724; return; }
    else if(stop_mass == 1440) { xsec = 0.000364616; xsec_unc = 0.232234; return; }
    else if(stop_mass == 1445) { xsec = 0.000353965; xsec_unc = 0.234637; return; }
    else if(stop_mass == 1450) { xsec = 0.000343923; xsec_unc = 0.234948; return; }
    else if(stop_mass == 1455) { xsec = 0.000333885; xsec_unc = 0.235468; return; }
    else if(stop_mass == 1460) { xsec = 0.000324344; xsec_unc = 0.23771; return; }
    else if(stop_mass == 1465) { xsec = 0.0003153; xsec_unc = 0.238004; return; }
    else if(stop_mass == 1470) { xsec = 0.00030583; xsec_unc = 0.240064; return; }
    else if(stop_mass == 1475) { xsec = 0.000296811; xsec_unc = 0.240314; return; }
    else if(stop_mass == 1480) { xsec = 0.000288149; xsec_unc = 0.239248; return; }
    else if(stop_mass == 1485) { xsec = 0.000279711; xsec_unc = 0.241257; return; }
    else if(stop_mass == 1490) { xsec = 0.000271724; xsec_unc = 0.241274; return; }
    else if(stop_mass == 1495) { xsec = 0.000264275; xsec_unc = 0.243545; return; }
    else if(stop_mass == 1500) { xsec = 0.000256248; xsec_unc = 0.24372; return; }
    else if(stop_mass == 1505) { xsec = 0.000248853; xsec_unc = 0.245827; return; }
    else if(stop_mass == 1510) { xsec = 0.000241844; xsec_unc = 0.246187; return; }
    else if(stop_mass == 1515) { xsec = 0.000234438; xsec_unc = 0.248442; return; }
    else if(stop_mass == 1520) { xsec = 0.000227374; xsec_unc = 0.248909; return; }
    else if(stop_mass == 1525) { xsec = 0.000221045; xsec_unc = 0.250895; return; }
    else if(stop_mass == 1530) { xsec = 0.000214431; xsec_unc = 0.248728; return; }
    else if(stop_mass == 1535) { xsec = 0.000208092; xsec_unc = 0.251043; return; }
    else if(stop_mass == 1540) { xsec = 0.000201748; xsec_unc = 0.253207; return; }
    else if(stop_mass == 1545) { xsec = 0.000196399; xsec_unc = 0.255641; return; }
    else if(stop_mass == 1550) { xsec = 0.000190474; xsec_unc = 0.255213; return; }
    else if(stop_mass == 1555) { xsec = 0.000185188; xsec_unc = 0.257329; return; }
    else if(stop_mass == 1560) { xsec = 0.000179263; xsec_unc = 0.256931; return; }
    else if(stop_mass == 1565) { xsec = 0.000174021; xsec_unc = 0.259111; return; }
    else if(stop_mass == 1570) { xsec = 0.000169176; xsec_unc = 0.258106; return; }
    else if(stop_mass == 1575) { xsec = 0.000163861; xsec_unc = 0.260597; return; }
    else if(stop_mass == 1580) { xsec = 0.000159583; xsec_unc = 0.262958; return; }
    else if(stop_mass == 1585) { xsec = 0.000154719; xsec_unc = 0.26195; return; }
    else if(stop_mass == 1590) { xsec = 0.000150506; xsec_unc = 0.264111; return; }
    else if(stop_mass == 1595) { xsec = 0.000145626; xsec_unc = 0.263077; return; }
    else if(stop_mass == 1600) { xsec = 0.000141382; xsec_unc = 0.265291; return; }
    else if(stop_mass == 1605) { xsec = 0.000137131; xsec_unc = 0.267424; return; }
    else if(stop_mass == 1610) { xsec = 0.000132187; xsec_unc = 0.26668; return; }
    else if(stop_mass == 1615) { xsec = 0.000127929; xsec_unc = 0.269117; return; }
    else if(stop_mass == 1620) { xsec = 0.000124086; xsec_unc = 0.267738; return; }
    else if(stop_mass == 1625) { xsec = 0.00011982; xsec_unc = 0.270483; return; }
    else if(stop_mass == 1630) { xsec = 0.000116042; xsec_unc = 0.268071; return; }
    else if(stop_mass == 1635) { xsec = 0.000112767; xsec_unc = 0.27127; return; }
    else if(stop_mass == 1640) { xsec = 0.000108936; xsec_unc = 0.269351; return; }
    else if(stop_mass == 1645) { xsec = 0.000105746; xsec_unc = 0.271783; return; }
    else if(stop_mass == 1650) { xsec = 0.000102693; xsec_unc = 0.27292; return; }
    else if(stop_mass == 1655) { xsec = 0.000100112; xsec_unc = 0.274445; return; }
    else if(stop_mass == 1660) { xsec = 9.75763e-05; xsec_unc = 0.275431; return; }
    else if(stop_mass == 1665) { xsec = 9.52062e-05; xsec_unc = 0.276946; return; }
    else if(stop_mass == 1670) { xsec = 9.29857e-05; xsec_unc = 0.277869; return; }
    else if(stop_mass == 1675) { xsec = 9.08285e-05; xsec_unc = 0.279347; return; }
    else if(stop_mass == 1680) { xsec = 8.87433e-05; xsec_unc = 0.281539; return; }
    else if(stop_mass == 1685) { xsec = 8.66618e-05; xsec_unc = 0.283509; return; }
    else if(stop_mass == 1690) { xsec = 8.46535e-05; xsec_unc = 0.284432; return; }
    else if(stop_mass == 1695) { xsec = 8.27102e-05; xsec_unc = 0.28591; return; }
    else if(stop_mass == 1700) { xsec = 8.07774e-05; xsec_unc = 0.287497; return; }
    else if(stop_mass == 1705) { xsec = 7.8666e-05; xsec_unc = 0.288194; return; }
    else if(stop_mass == 1710) { xsec = 7.6572e-05; xsec_unc = 0.290265; return; }
    else if(stop_mass == 1715) { xsec = 7.45994e-05; xsec_unc = 0.291193; return; }
    else if(stop_mass == 1720) { xsec = 7.25199e-05; xsec_unc = 0.293013; return; }
    else if(stop_mass == 1725) { xsec = 7.05189e-05; xsec_unc = 0.293697; return; }
    else if(stop_mass == 1730) { xsec = 6.85712e-05; xsec_unc = 0.294972; return; }
    else if(stop_mass == 1735) { xsec = 6.67296e-05; xsec_unc = 0.296167; return; }
    else if(stop_mass == 1740) { xsec = 6.49184e-05; xsec_unc = 0.297686; return; }
    else if(stop_mass == 1745) { xsec = 6.30949e-05; xsec_unc = 0.298524; return; }
    else if(stop_mass == 1750) { xsec = 6.13637e-05; xsec_unc = 0.299789; return; }
    else if(stop_mass == 1755) { xsec = 5.97301e-05; xsec_unc = 0.300928; return; }
    else if(stop_mass == 1760) { xsec = 5.80751e-05; xsec_unc = 0.302585; return; }
    else if(stop_mass == 1765) { xsec = 5.65479e-05; xsec_unc = 0.30366; return; }
    else if(stop_mass == 1770) { xsec = 5.49998e-05; xsec_unc = 0.305241; return; }
    else if(stop_mass == 1775) { xsec = 5.35686e-05; xsec_unc = 0.306718; return; }
    else if(stop_mass == 1780) { xsec = 5.20828e-05; xsec_unc = 0.306799; return; }
    else if(stop_mass == 1785) { xsec = 5.07079e-05; xsec_unc = 0.309201; return; }
    else if(stop_mass == 1790) { xsec = 4.93948e-05; xsec_unc = 0.310043; return; }
    else if(stop_mass == 1795) { xsec = 4.80635e-05; xsec_unc = 0.31138; return; }
    else if(stop_mass == 1800) { xsec = 4.67492e-05; xsec_unc = 0.312291; return; }
    else if(stop_mass == 1805) { xsec = 4.55055e-05; xsec_unc = 0.314321; return; }
    else if(stop_mass == 1810) { xsec = 4.42835e-05; xsec_unc = 0.315499; return; }
    else if(stop_mass == 1815) { xsec = 4.30744e-05; xsec_unc = 0.316302; return; }
    else if(stop_mass == 1820) { xsec = 4.19954e-05; xsec_unc = 0.317151; return; }
    else if(stop_mass == 1825) { xsec = 4.08527e-05; xsec_unc = 0.319048; return; }
    else if(stop_mass == 1830) { xsec = 3.97561e-05; xsec_unc = 0.319718; return; }
    else if(stop_mass == 1835) { xsec = 3.87041e-05; xsec_unc = 0.322028; return; }
    else if(stop_mass == 1840) { xsec = 3.76008e-05; xsec_unc = 0.32268; return; }
    else if(stop_mass == 1845) { xsec = 3.66914e-05; xsec_unc = 0.324529; return; }
    else if(stop_mass == 1850) { xsec = 3.56995e-05; xsec_unc = 0.325039; return; }
    else if(stop_mass == 1855) { xsec = 3.47689e-05; xsec_unc = 0.326767; return; }
    else if(stop_mass == 1860) { xsec = 3.38528e-05; xsec_unc = 0.328878; return; }
    else if(stop_mass == 1865) { xsec = 3.29644e-05; xsec_unc = 0.328975; return; }
    else if(stop_mass == 1870) { xsec = 3.20679e-05; xsec_unc = 0.329608; return; }
    else if(stop_mass == 1875) { xsec = 3.12583e-05; xsec_unc = 0.331541; return; }
    else if(stop_mass == 1880) { xsec = 3.04342e-05; xsec_unc = 0.333117; return; }
    else if(stop_mass == 1885) { xsec = 2.96516e-05; xsec_unc = 0.332866; return; }
    else if(stop_mass == 1890) { xsec = 2.88952e-05; xsec_unc = 0.336279; return; }
    else if(stop_mass == 1895) { xsec = 2.81145e-05; xsec_unc = 0.336845; return; }
    else if(stop_mass == 1900) { xsec = 2.73974e-05; xsec_unc = 0.338247; return; }
    else if(stop_mass == 1905) { xsec = 2.66796e-05; xsec_unc = 0.339708; return; }
    else if(stop_mass == 1910) { xsec = 2.59941e-05; xsec_unc = 0.339526; return; }
    else if(stop_mass == 1915) { xsec = 2.52784e-05; xsec_unc = 0.341137; return; }
    else if(stop_mass == 1920) { xsec = 2.46598e-05; xsec_unc = 0.342714; return; }
    else if(stop_mass == 1925) { xsec = 2.39932e-05; xsec_unc = 0.342328; return; }
    else if(stop_mass == 1930) { xsec = 2.33737e-05; xsec_unc = 0.34394; return; }
    else if(stop_mass == 1935) { xsec = 2.27623e-05; xsec_unc = 0.345138; return; }
    else if(stop_mass == 1940) { xsec = 2.21454e-05; xsec_unc = 0.346933; return; }
    else if(stop_mass == 1945) { xsec = 2.15924e-05; xsec_unc = 0.350815; return; }
    else if(stop_mass == 1950) { xsec = 2.10232e-05; xsec_unc = 0.349444; return; }
    else if(stop_mass == 1955) { xsec = 2.05211e-05; xsec_unc = 0.350155; return; }
    else if(stop_mass == 1960) { xsec = 1.98996e-05; xsec_unc = 0.352135; return; }
    else if(stop_mass == 1965) { xsec = 1.9408e-05; xsec_unc = 0.353328; return; }
    else if(stop_mass == 1970) { xsec = 1.88974e-05; xsec_unc = 0.354643; return; }
    else if(stop_mass == 1975) { xsec = 1.84612e-05; xsec_unc = 0.357904; return; }
    else if(stop_mass == 1980) { xsec = 1.79562e-05; xsec_unc = 0.358898; return; }
    else if(stop_mass == 1985) { xsec = 1.75673e-05; xsec_unc = 0.35989; return; }
    else if(stop_mass == 1990) { xsec = 1.70612e-05; xsec_unc = 0.360953; return; }
    else if(stop_mass == 1995) { xsec = 1.66228e-05; xsec_unc = 0.364709; return; }
    else if(stop_mass == 2000) { xsec = 1.62355e-05; xsec_unc = 0.365277; return; }
    else { xsec = 0.; xsec_unc = 0.; }
  }

}
