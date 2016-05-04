// merge_ntuples.cxx: Merges ntuples tree and treeglobal

#include <ctime>
#include <fstream>
#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  time_t startTime, endtime;
  time(&startTime);

  if(argc < 4) {
    cout<<endl<<"Requires at least 3 arguments: ./plot/merge_ntuples.exe infolder outfolder ntuples_tag <out_tag>"<<endl<<endl;;
    return 1;
  }

  TString folder(argv[1]), outfolder(argv[2]), ntuples_tag(argv[3]), out_tag("");
  gSystem->mkdir(outfolder, kTRUE);
  if(argc>4) out_tag = argv[4];
  TString ntuples = folder+"/*"+ntuples_tag+"*";

  // Merging tree TTrees
  TChain chain("tree");
  int nfiles = chain.Add(ntuples);
  TString outname = outfolder+"/mergedbaby_"+ntuples_tag+"_"+out_tag+"_nfiles_"+to_string(nfiles)+".root";
  chain.Merge(outname);

  // Merging treeglobal TTrees
  TChain chaing("treeglobal");
  chaing.Add(ntuples);
  TTree *tglobal = chaing.CopyTree("1");
  tglobal->SetDirectory(0);
  TFile rootfile(outname, "UPDATE");
  rootfile.cd();
  tglobal->Write();
  rootfile.Close();

  time(&endtime); 
  cout<<"Took "<<difftime(endtime, startTime)<<" seconds to merge "<<ntuples<<" into "<<outname<<endl<<endl;  
 }
