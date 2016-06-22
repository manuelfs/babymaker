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
#include <ctime>

using namespace std;

int main(int argc, char *argv[]){
  time_t startTime;
  time(&startTime);

  if(argc < 1) {
    cout<<endl<<"Requires 1 argument: "<<"./run/skim_scan_onefile.exe infile";
    return 1;
  }

  TString infiles(argv[1]);
  infiles += "/*.root";
  TString outpath = "";
  if (argc>1) outpath = argv[2];
  TChain tree("tree");
  tree.Add(infiles);
  TChain treeglobal("treeglobal");
  treeglobal.Add(infiles);  
  //long nentries(tree.GetEntries());

  TString outfolder = outpath;
  outfolder.Remove(outfolder.Last('/')+1, outfolder.Length());
  if(outfolder!="") gSystem->mkdir(outfolder, kTRUE);

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
  vector<TString> mass_tag;
  int Npoints=0;
  for(int iy=ini_y; iy<=last_y; iy++){
    for(int ix=ini_x; ix<=last_x; ix++){
      if(mass_plane->GetBinContent(ix,iy) > 0){
        int mglui = static_cast<int>(mass_plane->GetYaxis()->GetBinCenter(iy));
        int mlsp = static_cast<int>(mass_plane->GetXaxis()->GetBinCenter(ix));
        pair_cuts.push_back(Form("mgluino==%i&&mlsp==%i",mglui,mlsp));
        mass_tag.push_back(Form("mGluino-%i_mLSP-%i",mglui,mlsp));
        cout<<"Found mass point "<<mass_tag.back()<<endl;
	Npoints++;
      }
    }
  }

  time_t curTime;
  time(&curTime);
  // char time_c[100];
  //struct tm * timeinfo = localtime(&curTime);
  //strftime(time_c,100,"%Y-%m-%d %H:%M:%S",timeinfo);
  int seconds(floor(difftime(curTime,startTime)+0.5));
  cout<<endl<<"Took "<<seconds<<" seconds to fill chain, project and find the "<<Npoints<<" mass pairs"<<endl<<endl;
 
  // Finding outfile names
  for(unsigned int ip = 0; ip<pair_cuts.size(); ip++){
    time_t startloop;
    time(&startloop);
    TString outfile = "";
    if (outpath==""){
      outfile=infiles;
      //outfile.Remove(0, outfile.Last('/')); outfile.ReplaceAll("*","");
      if(outfile.Contains(".root")) outfile.ReplaceAll(".root","_"+pair_cuts.at(ip)+".root");
      else outfile += ("_"+pair_cuts.at(ip)+".root");
      outfile.ReplaceAll(">=","ge"); outfile.ReplaceAll("<=","se"); outfile.ReplaceAll("&&","_");
      outfile.ReplaceAll(">","g"); outfile.ReplaceAll("<","s"); outfile.ReplaceAll("=","");
      outfile.ReplaceAll("(",""); outfile.ReplaceAll(")",""); outfile.ReplaceAll("+","");
      outfile.ReplaceAll("[",""); outfile.ReplaceAll("]",""); outfile.ReplaceAll("|","_");
      outfile.ReplaceAll("/","");
      //outfile = outfolder+outfile;
    } else {
      outfile = outpath;
      outfile.ReplaceAll("MASS_TAG",mass_tag[ip]);
    }
      
    TFile out_rootfile(outfile.Data(), "RECREATE");
    if(out_rootfile.IsZombie() || !out_rootfile.IsOpen()) return 1;
    out_rootfile.cd();

    //cout<<"Skimming the "<<nfiles<<" files in "<<infiles<<endl;
    TTree * ctree = tree.CopyTree(pair_cuts.at(ip));
    TTree * ctreeglobal = treeglobal.CopyTree("1");
    if(ctree) ctree->Write();
    else cout<<"Could not find tree in "<<infiles<<endl;
    if(ctreeglobal)  ctreeglobal->Write();
    else cout<<"Could not find treeglobal in "<<infiles<<endl;
    
    time_t endloop;
    time(&endloop);
    int secs(floor(difftime(endloop,startloop)+0.5));
    cout<<"Written "<<outfile<<" with "<<ctree->GetEntries()<<" entries in "<<secs<<" seconds"<<endl;
    out_rootfile.Close();
  }
  time(&curTime);
  cout<<"Took "<< difftime(curTime,startTime)<<" seconds to skim "<< pair_cuts.size()<<" files."<<endl;
}

