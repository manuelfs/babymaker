// merge_ntuples.cxx: Merges ntuples tree and treeglobal

#include <ctime>
#include <fstream>
#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TString.h"

using namespace std;
using std::cout;
using std::endl;

int main(int argc, char *argv[]){
  time_t startTime, endtime;
  time(&startTime);

  if(argc < 3) {
    cout<<endl<<"Required at least 2 arguments: "
	<<"./plot/merge_ntuples.exe ntuples name"<<endl<<endl;;
    return 1;
  }

  TString ntuples(argv[1]), rootname(argv[2]);
  if(!rootname.Contains(".root")) rootname += ".root";

  // Merging tree TTrees
  TChain chain("tree");
  chain.Add(ntuples);
  chain.Merge(rootname);

  // Merging treeglobal TTrees
  TChain chaing("treeglobal");
  chaing.Add(ntuples);
  TTree *tglobal = chaing.CopyTree("1");
  tglobal->SetDirectory(0);
  TFile rootfile(rootname, "UPDATE");
  rootfile.cd();
  tglobal->Write();
  rootfile.Close();

  time(&endtime); 
  cout<<"Took "<<difftime(endtime, startTime)<<" seconds to merge "<<ntuples<<" into "<<rootname<<endl<<endl;  
 }
