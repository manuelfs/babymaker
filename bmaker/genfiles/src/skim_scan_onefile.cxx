// skim_scan_onefile.cxx: Separate signal scan by mass point


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

int main(int argc, char *argv[]){

  if(argc < 1) {
    cout<<endl<<"Requires 1 argument: "
	<<"./run/skim_scan_onefile.exe infile";
    return 1;
  }

  TString infiles(argv[1]);
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
	pair_cuts.push_back(Form("mgluino==%i&&mlsp==%i",mglui,mlsp));
      }
    }
  }
 

  /*TString folder(infiles);
  vector<TString > outfiles;
  folder.Remove(folder.Last('/')+1, folder.Length());*/
  // gSystem->mkdir(outfolder, kTRUE);
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
    outfile.ReplaceAll("/","");
    //outfile = outfolder+outfile;
  
    //cout<<"outfile is"<<outfile<<endl;
    // Checking if output file exists
    /*TString outname(outfile);
    //outname.ReplaceAll(outfolder, ""); outname.ReplaceAll("/", "");
    vector<TString> existing_outfiles = dirlist(outfolder, outname);
    if(existing_outfiles.size()>0) {
      cout<<"File "<<outfile<<" exists. Skipping"<<endl;
      continue;
      }*/
      
    //  cout<<"creating outfile"<<endl;
    TFile out_rootfile(outfile, "CREATE");
    if(out_rootfile.IsZombie() || !out_rootfile.IsOpen()) return 1;
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

