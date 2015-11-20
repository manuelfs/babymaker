#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "utilities.hh"

using namespace std;

void change_branch_one(TString indir, TString name, TString outdir, float newweight);

int main(int argc, char *argv[]){

  if(argc<3){
    cout<<"Format: ./run/change_branch.exe <infolder> <outfolder> <weight>"<<endl;
    return 1;
  }

  TString folder(argv[1]), outfolder(argv[2]), string_weight(argv[3]);
  if(!folder.EndsWith("/")) folder.Append("/");
  if(!outfolder.EndsWith("/")) outfolder.Append("/");
  float weight = string_weight.Atof();

  vector<TString> files = dirlist(folder,".root");
  
  for(unsigned int i=0; i<files.size(); i++){
    change_branch_one(folder, files[i], outfolder, weight);
  }
}

void change_branch_one(TString indir, TString name, TString outdir, float newweight){
    
  //Set up old tree and branch
  TFile* oldfile = new TFile(indir+name);
  TTree* oldtree = static_cast<TTree*>(oldfile->Get("tree"));

  float weight_ = 0;
  oldtree->SetBranchAddress("weight", &weight_);

  //Set up new tree
  name.ReplaceAll(".root","_mod.root");
  TFile* newfile = new TFile(outdir+name,"recreate");
  TTree* newtree = oldtree->CloneTree(0);

  //Loop and fill events with new weights
  int nentries = oldtree->GetEntries();
  for(int i=0; i<nentries; i++){
    oldtree->GetEntry(i);
    weight_ = newweight;
    newtree->Fill();
  }
  //Save tree
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}
