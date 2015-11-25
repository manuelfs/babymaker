#ifndef H_CHANGE_BRANCH
#define H_CHANGE_BRANCH

#include <deque>
#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include <vector>
#include "utilities.hh"

using namespace std;

void change_branch_one(TString indir, TString name, TString outdir, vector<TString> var_type, vector<TString> var, vector<TString> var_val){

  if(var_type.size()!=var.size() || var_type.size()!=var_val.size())
    { cout<<"Error: Branch vectors are not the same size"<<endl; exit(0); }

  const int nvar = var.size();

  //Setup
  vector<bool> multiply(nvar,false);
  for(int isetup=0; isetup<nvar; isetup++){
    //Set multiply
    if(var_val[isetup].BeginsWith("*") || var_val[isetup].EndsWith("*")){
      var_val[isetup].Remove(TString::kBoth,'*');
      multiply[isetup]=true;
    }
    
    //Handle bools
    if(var_val[isetup]=="true")
      var_val[isetup]="1";
    else if(var_val[isetup]=="false")
      var_val[isetup]="0";
  }

  //Set up old tree and branch
  TFile* oldfile = new TFile(indir+name);
  TTree* oldtree = static_cast<TTree*>(oldfile->Get("tree"));

  vector<int> new_var_int_(nvar,-999);
  vector<float> new_var_flt_(nvar,-999);
  vector<double> new_var_dbl_(nvar,-999);
  deque<bool> new_var_bool_(nvar, false); // vector<bool>is not a vector in C++... so can't pass bools by reference
  
  //Branches
  for(int ibch=0; ibch<nvar; ibch++){
    if(var_type[ibch]=="int")         {  new_var_int_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_int_[ibch]); }
    else if(var_type[ibch]=="float")  {  new_var_flt_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_flt_[ibch]); }
    else if(var_type[ibch]=="double") {  new_var_dbl_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_dbl_[ibch]); }
    else if(var_type[ibch]=="bool")   {  new_var_bool_[ibch] = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_bool_[ibch]);
      if(multiply[ibch])  cout<<"Warning: Multiplying branch of type \"bool\". Skipping branch."<<endl;}
    else {cout<<"Error: Branch type invalid. Use only \"int\", \"float\", \"double\", or \"bool\""<<endl; exit(0);}
  }

  //Set up new tree
  name.ReplaceAll(".root","_mod.root");
  TFile* newfile = new TFile(outdir+name,"recreate");
  TTree* newtree = oldtree->CloneTree(0);

  //Loop and fill events with new weights
  int nentries = oldtree->GetEntries();
  for(int i=0; i<nentries; i++){
    oldtree->GetEntry(i);
    
    //Set vars
    for(int iset=0; iset<nvar; iset++){
      if(!multiply[iset]){
	if(var_type[iset]=="int")             new_var_int_[iset]  = var_val[iset].Atoi();
	else if(var_type[iset]=="float")      new_var_flt_[iset]  = var_val[iset].Atof();
	else if(var_type[iset]=="double")     new_var_dbl_[iset]  = static_cast<double>(var_val[iset].Atof());
	else if(var_type[iset]=="bool")       new_var_bool_[iset] = var_val[iset].Atoi();
      }
      else {
	if(var_type[iset]=="int")             new_var_int_[iset]  *= var_val[iset].Atoi(); 
	else if(var_type[iset]=="float")      new_var_flt_[iset]  *= var_val[iset].Atof(); 
	else if(var_type[iset]=="double")     new_var_dbl_[iset]  *= var_val[iset].Atof(); 
      }
    }
    newtree->Fill();
  }
  //Save tree
  newtree->AutoSave();
  delete oldfile;
  delete newfile;
}

#endif
