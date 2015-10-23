// skim_ntuples.cxx: Skims reduced trees
// USAGE: ./plot/skim_ntuples.exe infolder outfolder [cuts=\"ht>500&&met>200\"] [njobs=-1] [ijob=-1]


#include "utilities.hh"

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

void onefile_skim(TString infiles, TString outfolder, TString cuts);

int main(int argc, char *argv[]){

  if(argc < 3) {
    cout<<endl<<"Required at least 2 arguments: "
	<<"./plot/skim_ntuples.exe infolder outfolder [cuts=\"ht>500&&met>200\"] "
	<<"[njobs=-1] [ijob=-1]"<<endl<<endl;;
    return 1;
  }
  TString folder(argv[1]), outfolder(argv[2]), cuts="ht>500&&met>200";
  if(argc >= 4) cuts = argv[3]; 
  unsigned njobs(0), ijob(0);
  if(argc >= 6) {
    njobs = atoi(argv[4]);
    ijob  = atoi(argv[5]);
  }
  vector<TString> files = dirlist(folder, ".root");
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
    onefile_skim(folder+"/"+files[file], outfolder, cuts);
  }
  return 0;
}

void onefile_skim(TString infiles, TString outfolder, TString cuts){
  TString folder(infiles), outfile(infiles);
  folder.Remove(folder.Last('/')+1, folder.Length());

  // Finding outfile name
  outfile.Remove(0, outfile.Last('/')); outfile.ReplaceAll("*","");
  if(outfile.Contains(".root")) outfile.ReplaceAll(".root","_"+cuts+".root");
  else outfile += ("_"+cuts+".root");
  outfile.ReplaceAll(">=","ge"); outfile.ReplaceAll("<=","se"); outfile.ReplaceAll("&","_");
  outfile.ReplaceAll(">","g"); outfile.ReplaceAll("<","s"); outfile.ReplaceAll("=","");
  outfile.ReplaceAll("(",""); outfile.ReplaceAll(")",""); outfile.ReplaceAll("+","");
  outfile.ReplaceAll("[",""); outfile.ReplaceAll("]",""); outfile.ReplaceAll("|","_");
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

