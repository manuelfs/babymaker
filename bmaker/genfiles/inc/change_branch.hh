#ifndef H_CHANGE_BRANCH
#define H_CHANGE_BRANCH

#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "utilities.hh"

using namespace std;

void change_branch_one(TString indir, TString name, TString outdir, TString var_type, TString var, TString var_val){

  //Set multiply
  bool multiply=false;
  if(var_val.BeginsWith("*") || var_val.EndsWith("*")){
    var_val.Remove(TString::kBoth,'*');
    multiply=true;
  }

  //Handle bools
  if(var_val=="true")
    var_val=="1";
  else if(var_val=="false")
    var_val=="0";
  
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

    if(!multiply){
      if(var_type=="int")             new_var_int_ = var_val.Atoi();
      else if(var_type=="float")      new_var_flt_ = var_val.Atof();
      else if(var_type=="double")     new_var_dbl_ = static_cast<double>(var_val.Atof());
      else if(var_type=="bool")       new_var_bool_ = var_val.Atoi();
    }
    else {
      if(var_type=="int")           { new_var_int_ *= var_val.Atoi(); }
      else if(var_type=="float")    { new_var_flt_ *= var_val.Atof(); }
      else if(var_type=="double")   { new_var_dbl_ *= var_val.Atof(); }
      else if(var_type=="bool")     { new_var_bool_ *= var_val.Atof(); cout<<"Warning multiplying branch of type \"bool\""<<endl;}
    }
    newtree->Fill();
  }
  //Save tree
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}

#endif
