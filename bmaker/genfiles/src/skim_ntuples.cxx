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
	<<"./plot/skim_ntuples.exe infolder outfolder [cuts=\"nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj>250\"] "
	<<"[njobs=-1] [ijob=-1]"<<endl<<endl;;
    return 1;
  }
  TString folder(argv[1]), outfolder(argv[2]), cuts="nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj>250";
  if(argc >= 4) cuts = argv[3]; 
  unsigned njobs(0), ijob(0);
  if(argc >= 6) {
    njobs = atoi(argv[4]);
    ijob  = atoi(argv[5]);
  }
  TString tag = cuts; // Using tag to avoid file names too long for TFile
  if(cuts=="abcd") cuts="nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj>250";
  if(cuts=="sys_abcd") 
    cuts = "nleps==1&&max(ht,Max$(sys_ht))>500&&max(met,Max$(sys_met))>200&&max(njets,Max$(sys_njets))>=6&&max(nbm,Max$(sys_nbm))>=1&&max(mj,Max$(sys_mj))>250";
  if(cuts=="zisr")
    cuts = "nvleps==2&&nleps>=1&&Max$(leps_pt)>30&&((elelv_m>80&&elelv_m<100)||(mumuv_m>80&&mumuv_m<100))";
  if(cuts=="dy_ht300")
    cuts = "nvleps==2&&nleps>=1&&Max$(leps_pt)>30&&((elelv_m>80&&elelv_m<100)||(mumuv_m>80&&mumuv_m<100))&&ht>300";
  if(cuts=="ttisr")
    cuts = "nvleps==2&&nleps>=1&&Max$(leps_pt)>30&&njets>=2&&nbm==2";
  if(cuts=="ttdilep_ht300")
    cuts = "nels==1&&nmus==1&&Max$(leps_pt)>30&&ht>300&&met>100&&nbm>=1";
  if(cuts=="qcd")
    cuts = "ht>1000&&met<50&&(nvmus+nvels)==0";
  if(cuts=="qcd_njet10")
     cuts = "ht>1000&&met<50&&(nvmus+nvels)==0&&njets>=10";
	     

  vector<TString> files = dirlist(folder, "*.root");
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
    onefile_skim(folder+"/"+files[file], outfolder, cuts, tag);
  }

  time_t curTime;
  time(&curTime);
  cout<<endl<<"Took "<< difftime(curTime,startTime)<<" seconds to skim "<< end-ini<<" files."<<endl<<endl;
}

void onefile_skim(TString infiles, TString outfolder, TString cuts, TString tag){
  TString folder(infiles), outfile(infiles);
  folder.Remove(folder.Last('/')+1, folder.Length());

  // Finding outfile name
  outfile.Remove(0, outfile.Last('/')); outfile.ReplaceAll("*","");
  if(outfile.Contains(".root")) outfile.ReplaceAll(".root","_"+tag+".root");
  else outfile += ("_"+tag+".root");
  outfile.ReplaceAll(">=","ge"); outfile.ReplaceAll("<=","se"); outfile.ReplaceAll("&","_");
  outfile.ReplaceAll(">","g"); outfile.ReplaceAll("<","s"); outfile.ReplaceAll("=","");
  outfile.ReplaceAll("(",""); outfile.ReplaceAll(")",""); outfile.ReplaceAll("+","");
  outfile.ReplaceAll("[",""); outfile.ReplaceAll("]",""); outfile.ReplaceAll("|","_");
  outfile.ReplaceAll("$",""); outfile.ReplaceAll(",","_");
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
  TTree *ctreeglobal = treeglobal.CopyTree("1");
  if(ctree) ctree->Write();
  else cout<<"Could not find tree in "<<infiles<<endl;
  if(ctreeglobal)   ctreeglobal->Write();
  else cout<<"Could not find treeglobal in "<<infiles<<endl;
  out_rootfile.Close();
  cout<<"Written "<<outfile<<" from "<<nfiles<<" files and "<<nentries<<" entries."<<endl;
}

