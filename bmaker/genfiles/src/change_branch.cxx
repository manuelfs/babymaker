#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "utilities.hh"

using namespace std;

void change_branch_one(TString indir, TString name, TString outdir, TString var_type, TString var, TString var_val);

int main(int argc, char *argv[]){

  if(argc<5){
    cout<<"Format: ./run/change_branch.exe <infolder> <outfolder> <branch type> <branch name> <branch value>"<<endl;
    return 1;
  }

  TString folder(argv[1]), outfolder(argv[2]), var_type(argv[3]), var(argv[4]), var_val(argv[5]);
  if(!folder.EndsWith("/")) folder.Append("/");
  if(!outfolder.EndsWith("/")) outfolder.Append("/");

  vector<TString> files = dirlist(folder,".root");
  
  for(unsigned int i=0; i<files.size(); i++){
    change_branch_one(folder, files[i], outfolder, var_type, var, var_val);
  }
}

void change_branch_one(TString indir, TString name, TString outdir, TString var_type, TString var, TString var_val){
    
  //Set up old tree and branch
  TFile* oldfile = new TFile(indir+name);
  TTree* oldtree = static_cast<TTree*>(oldfile->Get("tree"));
  
  int new_var_int_;
  float new_var_flt_;
  double new_var_dbl_;
  bool new_var_bool_;

  if(var_type=="int"){           new_var_int_ = 0;     oldtree->SetBranchAddress(var, &new_var_int_); }
  else if(var_type=="float"){    new_var_flt_ = 0;   oldtree->SetBranchAddress(var, &new_var_flt_); }
  else if(var_type=="double"){   new_var_dbl_ = 0;  oldtree->SetBranchAddress(var, &new_var_dbl_); }
  else if(var_type=="bool"){     new_var_bool_ = 0;    oldtree->SetBranchAddress(var, &new_var_bool_); }
  else {cout<<"Branch type invalid: Use only \"int\", \"float\", \"double\", or \"bool\""<<endl; exit(0);}

  //Set up new tree
  name.ReplaceAll(".root","_mod.root");
  TFile* newfile = new TFile(outdir+name,"recreate");
  TTree* newtree = oldtree->CloneTree(0);

  //Loop and fill events with new weights
  int nentries = oldtree->GetEntries();
  for(int i=0; i<nentries; i++){
    oldtree->GetEntry(i);

    if(var_type=="int")             new_var_int_ = var_val.Atoi();
    else if(var_type=="float")      new_var_flt_ = var_val.Atof();
    else if(var_type=="double")     new_var_dbl_ = static_cast<double>(var_val.Atof());
    else if(var_type=="bool")       new_var_bool_ = var_val.Atoi();
    else {cout<<"Branch type invalid: Not setting variable correctly"<<endl; exit(0);}

    newtree->Fill();
  }
  //Save tree
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}
