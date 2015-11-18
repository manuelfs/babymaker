// skim_scan.cxx: Separate signal scan by mass point
// USAGE: ./plot/skim_ntuples.exe infolder outfolder [cuts=\"ht>500&&met>200\"] [njobs=-1] [ijob=-1]


#include "utilities.hh"

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>     /* atoi */
#include "TH2.h"
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
	<<"./plot/skim_ntuples.exe infolder outfolder [cuts=\"\"] "
	<<"[njobs=-1] [ijob=-1]"<<endl<<endl;;
    return 1;
  }
  TString folder(argv[1]), outfolder(argv[2]), cuts="";
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
  TChain tree("tree");
  int nfiles = tree.Add(infiles);
  TChain treeglobal("treeglobal");
  treeglobal.Add(infiles);

  //Project tree into mass plane to find points for skim
  TH2F * mass_plane = new TH2F("mglu_vs_mlsp","mglu_vs_mlsp",3000,-0.5,2999.5,3000,-0.5,2999.5);
  tree.Project("mglu_vs_mlsp","mgluino:mlsp","","colz");
 
  //find relevant range to loop over
  int ini_x = mass_plane->FindFirstBinAbove(0,1);
  int last_x = mass_plane->FindLastBinAbove(0,1);
  int ini_y = mass_plane->FindFirstBinAbove(0,2);
  int last_y = mass_plane->FindLastBinAbove(0,2);

  //load all pairs as cuts into vector
  vector<TString> pair_cuts;
  for(int iy=ini_y; iy<=last_y; iy++){
    for(int ix=ini_x; ix<=last_x; ix++){
      if(mass_plane->GetBinContent(ix,iy) > 0){
	int mglui = static_cast<int>(mass_plane->GetYaxis()->GetBinCenter(iy));
	int mlsp = static_cast<int>(mass_plane->GetXaxis()->GetBinCenter(ix));
	if(cuts =="") pair_cuts.push_back(Form("mgluino==%i&&mlsp==%i",mglui,mlsp));
	else pair_cuts.push_back(Form("mgluino==%i&&mlsp==%i&&",mglui,mlsp)+cuts); 
      }
    }
  }
 

  TString folder(infiles);
  vector<TString > outfiles;
  folder.Remove(folder.Last('/')+1, folder.Length());
  gSystem->mkdir(outfolder, kTRUE);
  // Finding outfile names
  for(unsigned int ip = 0; ip<pair_cuts.size(); ip++){
   
    TString outfile=infiles;
    outfile.Remove(0, outfile.Last('/')); outfile.ReplaceAll("*","");
    if(outfile.Contains(".root")) outfile.ReplaceAll(".root","_"+pair_cuts.at(ip)+".root");
    else outfile += ("_"+pair_cuts.at(ip)+".root");
    outfile.ReplaceAll(">=","ge"); outfile.ReplaceAll("<=","se"); outfile.ReplaceAll("&&","_");
    outfile.ReplaceAll(">","g"); outfile.ReplaceAll("<","s"); outfile.ReplaceAll("=","");
    outfile.ReplaceAll("(",""); outfile.ReplaceAll(")",""); outfile.ReplaceAll("+","");
    outfile.ReplaceAll("[",""); outfile.ReplaceAll("]",""); outfile.ReplaceAll("|","_");
    outfile = outfolder+outfile;
  
    //cout<<"outfile is"<<outfile<<endl;
    // Checking if output file exists
    TString outname(outfile);
    outname.ReplaceAll(outfolder, ""); outname.ReplaceAll("/", "");
    vector<TString> existing_outfiles = dirlist(outfolder, outname);
    if(existing_outfiles.size()>0) {
      cout<<"File "<<outfile<<" exists. Skipping"<<endl;
      continue;
    }
      

    //  cout<<"creating outfile"<<endl;
    TFile out_rootfile(outfile, "CREATE");
    if(out_rootfile.IsZombie() || !out_rootfile.IsOpen()) return;
    out_rootfile.cd();


    //cout<<"Skimming the "<<nfiles<<" files in "<<infiles<<endl;
    long nentries(tree.GetEntries());
    TTree * ctree = tree.CopyTree(pair_cuts.at(ip));
    TTree * ctreeglobal = treeglobal.CopyTree("1");
    if(ctree) ctree->Write();
    else cout<<"Could not find tree in "<<infiles<<endl;
    if(ctreeglobal)  ctreeglobal->Write();
    else cout<<"Could not find treeglobal in "<<infiles<<endl;
    cout<<"Written "<<outfile<<" from "<<nfiles<<" files and "<<nentries<<" entries."<<endl;
    out_rootfile.Close();
    
  }
}

