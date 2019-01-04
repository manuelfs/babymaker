// skim_ntuples.cxx: Skims reduced trees
// USAGE: ./plot/skim_ntuples.exe infolder outfolder [cuts=\"ht>500&&met>200\"] [njobs=-1] [ijob=-1]


#include "utilities.hh"

#include <ctime>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>     /* atoi */

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TDirectory.h"

using namespace std;
using std::cout;
using std::endl;

void onefile_skim(TString infiles, TString outfolder, TString cuts, TString tag);

int main(int argc, char *argv[]){
  time_t startTime;
  time(&startTime);

  if(argc < 3) {
    cout<<endl<<"Required at least 2 arguments: "
	<<"./plot/skim_ntuples.exe infolder outfolder [cuts=\"nleps==1&&st>500&&met>200&&njets>=6&&nbdm>=1&&mj14>250&&nveto==0\"] "
	<<"[njobs=0] [ijob=0] [filetag=\"\"]"<<endl<<endl;;
    return 1;
  }
  TString folder(argv[1]), outfolder(argv[2]), cuts="nleps==1&&st>500&&met>200&&njets>=6&&nbdm>=1&&mj14>250&&nveto==0";
  if(argc >= 4) cuts = argv[3]; 
  unsigned njobs(0), ijob(0);
  if(argc >= 6) {
    njobs = atoi(argv[4]);
    ijob  = atoi(argv[5]);
  }
  TString filetag="";
  if(argc >= 7) filetag = argv[6];

  TString tag = cuts; // Using tag to avoid file names too long for TFile  
  if(cuts=="standard") cuts="nleps>=1&&st>500&&met>100";
  if(cuts=="nl1met200ht500") cuts="nleps>=1&&ht>500&&met>200";
  if(cuts=="stdnj5") cuts="nleps>=1 && st>500 && met>100 && njets>=5";
  if(cuts=="abcd") cuts="nleps==1&&st>500&&met>200&&njets>=6&&nbdm>=1&&mj14>250&&nveto==0";
  if(cuts=="baseline") cuts="nleps==1&&st>500&&met>200&&njets>=6&&nbdm>=1";
  if(cuts=="sys_abcd") 
    cuts = "nleps==1&&max(st,Max$(sys_st))>500&&max(met,Max$(sys_met))>200&&max(njets,Max$(sys_njets))>=6&&max(nbdm,Max$(sys_nbdm))>=1&&max(mj14,Max$(sys_mj14))>250";
  if(cuts=="zcand")
    cuts = "nleps==2&&Max$(leps_pt)>40&&((elel_m>80&&elel_m<100)||(mumu_m>80&&mumu_m<100))";
  if(cuts=="dy_ht300")
    cuts = "nvleps==2&&nleps>=1&&Max$(leps_pt)>30&&((elelv_m>80&&elelv_m<100)||(mumuv_m>80&&mumuv_m<100))&&ht>300";
  if(cuts=="ttisr")
    cuts = "nleps==2&&Max$(leps_pt)>40&&nbdm==2";
  if(cuts=="wisr")
    cuts = "met>100&&Max$(leps_pt)>40&&nbdl==0";
  if(cuts=="wisrht200")
    cuts = "ht>200&&met>100&&Max$(leps_pt)>40&&nbdl==0";
  if(cuts=="qcd")
    cuts = "ht>1000&&met<50&&(nvmus+nvels)==0";
  if(cuts=="qcd_njet10")
     cuts = "ht>1000&&met<50&&(nvmus+nvels)==0&&njets>=10";
  if(cuts=="mm_std")
	 cuts="Sum$(mm_nleps>=1&&mm_st>500.&&mm_met>200.)>0";	     
  if(cuts=="mm_std_nj5mj250")
	 cuts="Sum$(mm_nleps>=1&&mm_st>500&&mm_met>200&&mm_njets>=5&&mm_mj14_lep>250)>0||Sum$(mm_nleps>=1&&mm_st>500&&mm_met>200&&mm_njets>=5&&mm_mj14_nolep>250)>0";
  if(cuts=="rpvfit")
    cuts = "max(st,Max$(sys_st))>1200&&max(njets,Max$(sys_njets))>=4&&max(nbdm,Max$(sys_nbdm))>=1&&(max(mj12,Max$(sys_mj12))>500)";
  if(cuts=="rpvregion")
    cuts = "((nleps==0&&ht>1500)||(nleps==1&&ht>1200))&&njets>=4&&mj12>=300&&(nbdm>=1||nbdm>=1)";
  if(cuts=="st1000")
    cuts = "max(st,Max$(sys_st))>1000";
  if(cuts=="llm60nj2")
    cuts = "(mumuv_m>60||elelv_m>60)&&njets>=2";
  if(cuts=="httrig") // Prescaled HT triggers for fake MET trigger efficiency measurement
    cuts = "pass&&nvleps==0&&(trig[11]||trig[12]||trig[47]||trig[48]||trig[49]||trig[50]||trig[51]||trig[52]||trig[53]||trig[54])";
  if(cuts=="httrig_w_leps") // Prescaled HT triggers
    cuts = "pass&&(trig[11]||trig[12]||trig[47]||trig[48]||trig[49]||trig[50]||trig[51]||trig[52]||trig[53]||trig[54])";
  
  //// Higgsino skims
  TString njcut = "njets>=4&&njets<=5";
  TString sys_njcut = "(njets==4||sys_njets[1]==4||sys_njets[2]==4||njets==5||sys_njets[1]==5||sys_njets[2]==5)";
  TString nbcut = "&&(nbt>=2||nbdt>=2)";
  TString sys_nbcut = "&&max(nbdt,Max$(sys_nbdt))>=2";
  TString zcand = "&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))>80&&(mumu_m*(mumu_m>0)+elel_m*(elel_m>0))<100";
  TString higtrim = "&&higd_drmax<2.2&&higd_dm<=40&&higd_am<=200";
  TString sys_higtrim = "&&min(higd_drmax,Min$(sys_higd_drmax))<2.2&&min(higd_dm,Min$(sys_higd_dm))<=40&&min(higd_am,Min$(sys_higd_am))<=200";
  if(cuts=="higqcd")   cuts = njcut      +"&& met>150 && nvleps==0";
  if(cuts=="higloose") cuts = njcut+nbcut+"&& met>150 && nvleps==0";
  if(cuts=="higtight") cuts = njcut+nbcut+"&& met>150 && nvleps==0 && ntks==0&&!low_dphi"+higtrim;
  if(cuts=="higsys")   cuts = sys_njcut+sys_nbcut+"&& max(met,Max$(sys_met))>150 && nvleps==0 && ntks==0&&!low_dphi"+sys_higtrim;
  if(cuts=="higlep1")  cuts = njcut+nbcut+"&& nleps==1 && Max$(leps_pt)>30";
  if(cuts=="higlep2")  cuts = njcut+zcand+"&& nleps==2 && Max$(leps_pt)>40";

  //// exploration skims
  if(cuts=="nl2nj3nbdm2zveto")
    cuts = "nleps>=2&&(elel_m<80||elel_m>100)&&(mumu_m<80||mumu_m>100)&&njets>=3&&nbdm>=2";

  if(cuts.Contains("mchi")){
    cuts.ReplaceAll("mchi","");
    cuts = "Sum$(mc_id==1000023&&mc_mass=="+cuts+")>0";
  }

  //// RA2/b skims
  TString pass="globalTightHalo2016Filter==1&&HBHENoiseFilter==1&&HBHEIsoNoiseFilter==1&&eeBadScFilter==1";
  pass += "&&EcalDeadCellTriggerPrimitiveFilter==1&&BadChargedCandidateFilter&&BadPFMuonFilter&&NVtx>0&&JetID";
  if(cuts=="ra2_qcd")     cuts = pass+"&&(@Electrons.size()+@Muons.size())==0&&NJets>=3";	     
  if(cuts=="ra2_ht300")   cuts = pass+"&&HT>300";	     
  if(cuts=="ra2_eht300")  cuts = pass+"&&Max$(Electrons.Pt()*(abs(Electrons.Eta())<2))>35&&HT>300";	     
  if(cuts=="ra2_zmht200") cuts = pass+"&&@ZCandidates.size()>=1&&MHT>200";	     
	     

  vector<TString> files = dirlist(folder, "*"+filetag+"*.root");
  unsigned nfiles(files.size()), ini(0), end(nfiles);
  if(njobs>0){
    if(ijob<1 || ijob>njobs){
      cout<<endl<<"You need to set the 5th argument between 1 and "<<njobs<<endl<<endl;
      return 1;
    }
    unsigned jobfiles = (nfiles+njobs-1)/njobs;
    unsigned nbigjobs = (nfiles+njobs-1)%njobs+1;
    if(ijob <= nbigjobs){
      ini = jobfiles*(ijob-1);
      end = ini + jobfiles;
    } else {
      ini = nbigjobs*jobfiles+(jobfiles-1)*(ijob-1-nbigjobs);
      end = ini + jobfiles-1;
    }
  }
  cout<<"Doing files "<<ini+1<<" to "<<end<<" out of "<<nfiles<<endl;
  for(unsigned file(ini); file < end; file++){
    cout<<file+1<<"/"<<nfiles<<": ";
    onefile_skim(folder+"/"+files[file], outfolder, cuts, tag);
  }

  time_t curTime;
  time(&curTime);
  int seconds = difftime(curTime,startTime);
  cout<<endl<<"Took "<< seconds << " seconds ("<<hoursMinSec(seconds)<<") to skim "<< end-ini<<" files."<<endl<<endl;
}

void onefile_skim(TString infiles, TString outfolder, TString cuts, TString tag){
  TString folder(infiles), outfile(infiles);
  folder.Remove(folder.Last('/')+1, folder.Length());

  // Finding outfile name
  outfile.Remove(0, outfile.Last('/')); outfile.ReplaceAll("*","");
  if(outfile.Contains(".root")) outfile.ReplaceAll(".root","_"+tag+".root");
  else outfile += ("_"+tag+".root");
  outfile.ReplaceAll(">=","GE"); outfile.ReplaceAll("<=","SE"); outfile.ReplaceAll("&","_");
  outfile.ReplaceAll(">","G"); outfile.ReplaceAll("<","S"); outfile.ReplaceAll("=","");
  outfile.ReplaceAll("(",""); outfile.ReplaceAll(")",""); outfile.ReplaceAll("+","");
  outfile.ReplaceAll("[",""); outfile.ReplaceAll("]",""); outfile.ReplaceAll("|","_");
  outfile.ReplaceAll("$",""); outfile.ReplaceAll(",","_"); outfile.ReplaceAll("!","NOT");
  outfile.ReplaceAll(" ",""); outfile.ReplaceAll("@",""); 
  outfile = outfolder+outfile;

  // Checking if output file exists
  TString outname(outfile);
  outname.ReplaceAll(outfolder, ""); outname.ReplaceAll("/", "");
  vector<TString> outfiles = dirlist(outfolder, outname);
  if(outfiles.size()>0) {
    cout<<"File "<<outfile<<" exists. Exiting"<<endl;
    return;
  }

  gSystem->mkdir(outfolder, kTRUE);
  TFile out_rootfile(outfile, "CREATE");
  if(out_rootfile.IsZombie() || !out_rootfile.IsOpen()) return;
  out_rootfile.cd();
  TChain tree("tree");
  int nfiles = tree.Add(infiles);
  TChain treeglobal("treeglobal");
  treeglobal.Add(infiles);

  //cout<<"Skimming the "<<nfiles<<" files in "<<infiles<<endl;
  long nentries(tree.GetEntries());
  TTree *ctree = tree.CopyTree(cuts);
  long centries=ctree->GetEntries();
  if(ctree) ctree->Write();
  else {
    cout<<"Could not find tree in "<<infiles<<endl;
    return;
  }

  //// treeglobal
  TFile infile(infiles);
  if(infile.Get("treeglobal") != 0){
    out_rootfile.cd();
    TTree *ctreeglobal = treeglobal.CopyTree("1");
    if(ctreeglobal)   ctreeglobal->Write();
    //else cout<<"Could not find treeglobal in "<<infiles<<endl;
  }
  out_rootfile.Close();
  cout<<"Written "<<centries<<" entries to "<<outfile<<" from "<<nfiles<<" files and "<<nentries<<" entries."<<endl;
}

