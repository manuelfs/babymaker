//----------------------------------------------------------------------------
// utilities - Various functions used accross the code
//----------------------------------------------------------------------------


#include <cmath>
#include <deque>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>   // setw

#include "TCollection.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TList.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TRegexp.h"

#include "cross_sections.hh"
#include "utilities.hh"

using namespace std;

//// Changes the value of branches in TTrees
int change_branch_one(TString indir, TString name, TString outdir, vector<TString> var_type, vector<TString> var,  
                      vector<vector<TString> > var_val, int totentries){

  if(var_type.size()!=var.size() || var_type.size()!=var_val.size())
    { cout<<"[Change Branch One] ERROR: Branch vectors are not the same size"<<endl; exit(0); }

  const int nvar = var.size();

  //Setup
  vector<varType> ivar_type;
  vector<bool> isLep(nvar, false);
  vector<vector<bool> > multiply(nvar);
  for(int isetup=0; isetup<nvar; isetup++){
    size_t vlength(var_val[isetup].size());
    multiply[isetup] = vector<bool>(vlength, false);
    //Set multiply
    for(unsigned int idx=0; idx<vlength; idx++){
      if(var_val[isetup][idx].BeginsWith("*") || var_val[isetup][idx].EndsWith("*")){
        var_val[isetup][idx].Remove(TString::kBoth,'*');
        multiply[isetup][idx]=true;
      }
    }

    //Handle bools
    if(var_type[isetup]=="bool" || var_type[isetup]=="vbool"){
      for(unsigned int idx=0; idx<var_val[isetup].size(); idx++){
        if(var_val[isetup][idx]=="true")
          var_val[isetup][idx]="1";
        else if(var_val[isetup][idx]=="false")
          var_val[isetup][idx]="0";
      }
    }

    // Converting string types to int and bool for speed
    if(var[isetup].Contains("_lep"))  isLep[isetup] = true;

    if(var_type[isetup]=="int")          ivar_type.push_back(kInt);
    else if(var_type[isetup]=="float")   ivar_type.push_back(kFloat);
    else if(var_type[isetup]=="double")  ivar_type.push_back(kDouble);
    else if(var_type[isetup]=="bool")    ivar_type.push_back(kBool);
    else if(var_type[isetup]=="vint")    ivar_type.push_back(kvInt);
    else if(var_type[isetup]=="vfloat")  ivar_type.push_back(kvFloat);
    else if(var_type[isetup]=="vdouble") ivar_type.push_back(kvDouble);
    else if(var_type[isetup]=="vbool")   ivar_type.push_back(kvBool);
    else {
      cout<<"var_type "<<var_type[isetup]<<" not supported. Exiting"<<endl;
      exit(1);
    }
  } // Loop over variables

  //Set up old tree and branch
  TFile* oldfile = new TFile(indir+name);
  TTree* oldtree = static_cast<TTree*>(oldfile->Get("tree"));
  TTree* oldtreeglobal = static_cast<TTree*>(oldfile->Get("treeglobal"));
  
  vector<int> new_var_int_(nvar,-999);
  vector<float> new_var_flt_(nvar,-999);
  vector<double> new_var_dbl_(nvar,-999);
  deque<bool> new_var_bool_(nvar, false); // vector<bool>is not a vector in C++... so can't pass bools by reference
  vector<vector<int> * > new_var_vint_(nvar);
  vector<vector<float> * > new_var_vflt_(nvar);
  vector<vector<double> * > new_var_vdbl_(nvar);
  vector<vector<bool> * > new_var_vbool_(nvar);
  
  // Hack to protect total weight from NaN, and not include w_pu
  float w_lumi(1.), w_lumi_old(1.), w_lep(1.), w_fs_lep(1.), w_btag_old(1.), w_corr(1.);
  float w_isr_old=1., w_pu_old=1.;
  //Branches
  int nleps_=0, mgluino_(0);
  float eff_trig_(0);
  oldtree->SetBranchAddress("eff_trig",&eff_trig_);
  oldtree->SetBranchAddress("nleps",&nleps_);
  oldtree->SetBranchAddress("mgluino",&mgluino_);

  for(int ibch=0; ibch<nvar; ibch++){
    switch (ivar_type[ibch]){
    default:
    case kInt:
      new_var_int_[ibch]   = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_int_[ibch]);  break;
    case kFloat:
      new_var_flt_[ibch]   = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_flt_[ibch]);  break;
    case kDouble:
      new_var_dbl_[ibch]   = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_dbl_[ibch]);  break;
    case kBool:
      new_var_bool_[ibch]  = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_bool_[ibch]); break;
    case kvInt:
      new_var_vint_[ibch]  = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_vint_[ibch]); break;
    case kvFloat:
      new_var_vflt_[ibch]  = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_vflt_[ibch]); break;
    case kvDouble:
      new_var_vdbl_[ibch]  = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_vdbl_[ibch]); break;
    case kvBool:
      new_var_vbool_[ibch] = 0;  oldtree->SetBranchAddress(var[ibch], &new_var_vbool_[ibch]);break;
    }
  }

  //Set up new tree
  name.ReplaceAll(".root","_renorm.root");
  TFile* newfile = new TFile(outdir+name,"recreate");
  TTree* newtree = oldtree->CloneTree(0);
  double eff_jetid(1.); // If not FastSim, we apply as pass
  if(name.Contains("FSPremix")) eff_jetid = 0.99;

  //Loop and fill events with new weights
  int nentries = oldtree->GetEntries();
  time_t begtime, endtime;
  time(&begtime);
  for(int i=0; i<nentries; i++){
    //if(i==10) exit(0);
    oldtree->GetEntry(i);

    w_corr = 1;
    if(i%500000==0 || i==nentries-1) {
      time(&endtime);
      int seconds(difftime(endtime,begtime));
      cout<<"Doing entry "<<setw(10)<<addCommas(i+1)<<" of "<<addCommas(nentries)
	  <<"    Took "<<setw(6)<<seconds<<" seconds at "
	  <<setw(4)<<roundNumber(i,1,seconds*1000.)<<" kHz"<<endl;
    }

    // float minpdf(1e10);
    
    //Set vars
    for(int iset=0; iset<nvar; iset++){
      // // Hack to recompute sys_pdf[1] which had a 1e-3 cut
      // if(var[iset].Contains("w_pdf")){
      //   for(unsigned int isys=0;isys<new_var_vflt_[iset]->size();isys++)  
      //     if(new_var_vflt_[iset]->at(isys) < minpdf) minpdf = new_var_vflt_[iset]->at(isys);
      // }  
      
      //// Events with leptons already have the correct w_lep
      if(isLep[iset] && nleps_!=0) {
        // Hack to protect total weight from NaN, and not include w_pu
        if(var[iset].Contains("w_lep"))    w_lep    = noNaN(new_var_flt_[iset]);
        if(var[iset].Contains("w_fs_lep")) w_fs_lep = noNaN(new_var_flt_[iset]);
        continue; // For lepton scale factors    
      }
      
      //// Saving values of variables before renormalization
      if(var[iset].Contains("w_pu"))  {w_pu_old  = noNaN(new_var_flt_[iset]); w_corr *= var_val[iset][0].Atof(); }
      if(var[iset].Contains("w_isr"))  {w_isr_old  = noNaN(new_var_flt_[iset]); w_corr *= var_val[iset][0].Atof(); }
      //if(var[iset].Contains("w_toppt"))  {w_toppt_old  = noNaN(new_var_flt_[iset]); w_corr *= var_val[iset][0].Atof(); }
      if(var[iset].Contains("w_btag"))   {w_btag_old   = noNaN(new_var_flt_[iset]); w_corr *= var_val[iset][0].Atof(); }
      if(var[iset].Contains("w_lumi"))   {w_lumi_old  = noNaN(new_var_flt_[iset]);  }

      // // Hack for empty pdf branches
      // if(var[iset].Contains("w_pdf")){
      //   if(new_var_vflt_[iset]->size()==0){
      //     new_var_vflt_[iset]->resize(100,1);
      //     if (i==0) cout<<"\n[Change Branch One] WARNING: Empty branch of \"w_pdf\". Setting values to 1"<<endl;
      //     continue;
      //   }
      // }
      // else if(var[iset].Contains("sys_pdf")){
      //   if(new_var_vflt_[iset]->size()==0){
      //     new_var_vflt_[iset]->resize(2,1);
      //     if (i==0) cout<<"[Change Branch One] WARNING: Empty branch of \"sys_pdf\". Setting values to 1"<<endl;
      //     continue;
      //   }
      // }
      for(unsigned int vidx=0; vidx<var_val[iset].size(); vidx++){
        if(!multiply[iset][vidx]){
          switch (ivar_type[iset]){
          default:
          case kInt:     new_var_int_[iset]              =  var_val[iset][vidx].Atoi();  break;
          case kFloat:   new_var_flt_[iset]              =  var_val[iset][vidx].Atof();  break;
          case kDouble:  new_var_dbl_[iset]              =  static_cast<double>(var_val[iset][vidx].Atof());  break;
          case kBool:    new_var_bool_[iset]             =  var_val[iset][vidx].Atoi();  break;
          case kvInt:    new_var_vint_[iset]->at(vidx)   =  var_val[iset][vidx].Atoi();  break;
          case kvFloat:  new_var_vflt_[iset]->at(vidx)   =  var_val[iset][vidx].Atof();  break;
          case kvDouble: new_var_vdbl_[iset]->at(vidx)   =  static_cast<double>(var_val[iset][vidx].Atof());  break;
          case kvBool:   new_var_vbool_[iset]->at(vidx)  =  var_val[iset][vidx].Atoi();  break;
          }
        } else {
          switch (ivar_type[iset]){
          default:
          case kInt:     new_var_int_[iset]             *=  var_val[iset][vidx].Atoi();  break;
          case kFloat:   new_var_flt_[iset]             *=  var_val[iset][vidx].Atof();  break;
          case kDouble:  new_var_dbl_[iset]             *=  static_cast<double>(var_val[iset][vidx].Atof());  break;
          case kvInt:    new_var_vint_[iset]->at(vidx)  *=  var_val[iset][vidx].Atoi();  break;
          case kvFloat:  new_var_vflt_[iset]->at(vidx)  *=  var_val[iset][vidx].Atof();  break;
          case kvDouble: new_var_vdbl_[iset]->at(vidx)  *=  static_cast<double>(var_val[iset][vidx].Atof());  break;
          case kBool:
          case kvBool:
            cout<<"[Change Branch One] WARNING: You cannot multiply Booleans. Skipping branch"<<endl;
            break;
          }
        } // if multiply
        if(ivar_type[iset] == kFloat)  {
	  if(isnan(new_var_flt_[iset])) {
	    new_var_flt_[iset] = 1.;
	    if(i==0) cout<<endl<<"==== WARNING: Branch \""<<var[iset]<<"\" is NaN ====="<<endl<<endl;
	  } else new_var_flt_[iset] = noNaN(new_var_flt_[iset]);
	}
        if(ivar_type[iset] == kvFloat) new_var_vflt_[iset]->at(vidx) = noNaN(new_var_vflt_[iset]->at(vidx));
      } // Loop over elements in each variable
      // if(var[iset].Contains("sys_pdf")) new_var_vflt_[iset]->at(1) = minpdf*var_val[iset][1].Atof(); 

      // Hack to protect total weight from NaN, and not include w_pu
      if(var[iset].Contains("w_lep"))    {w_lep    = new_var_flt_[iset]; w_corr *= var_val[iset][0].Atof(); }
      if(var[iset].Contains("w_fs_lep")) {w_fs_lep = new_var_flt_[iset]; w_corr *= var_val[iset][0].Atof(); }
      if(var[iset].Contains("w_lumi"))   {
	w_lumi   = new_var_flt_[iset]; 
	if(w_lumi_old<0) {
	  w_lumi *= -1;
	  new_var_flt_[iset] *= -1;
	}
	w_corr *= w_lumi/w_lumi_old; // Not currently used
      }
    } // Loop over variables
    // Hack to protect total weight from NaN, and not include w_pu
    for(int iset=0; iset<nvar; iset++)
      if(var[iset].Contains("weight")){
        new_var_flt_[iset] = w_lep * w_fs_lep * w_btag_old * w_isr_old * w_pu_old * w_lumi * eff_jetid 
	  * var_val[iset][0].Atof();
      }
    newtree->Fill();
  }
  
  unsigned int new_nev_sample(0);
  float new_xsec(0.);
  float xsec(0.), exsec(0.);
  oldtreeglobal->SetBranchAddress("nev_sample", &new_nev_sample);
  if(mgluino_>0) oldtreeglobal->SetBranchAddress("xsec", &new_xsec);
  TTree* newtreeglobal = oldtreeglobal->CloneTree(0);
  long nentriesg = oldtreeglobal->GetEntries();
  for(int i=0; i<nentriesg; i++){
    oldtreeglobal->GetEntry(i);
    new_nev_sample = totentries;
    if(mgluino_>0) {
      if(name.Contains("T1") || name.Contains("T5")) xsec::signalCrossSection(mgluino_, xsec, exsec);
      else  xsec::stopCrossSection(mgluino_, xsec, exsec);
      new_xsec = xsec;
    }
    newtreeglobal->Fill();
  }
  //Save tree
  newtree->AutoSave();
  newtreeglobal->AutoSave();
  delete oldfile;
  delete newfile;
  return nentries;
}

// Returns list of directorites or files in folder
vector<TString> dirlist(const TString &folder,
                        const TString &inname,
                        const TString &tag){
  TRegexp regex_tag(tag,true), regex_inname(inname, true);
  TString pwd(gSystem->pwd());
  vector<TString> v_dirs;
  TSystemDirectory dir(folder, folder);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=static_cast<TSystemFile*>(next()))) {
      fname = file->GetName();
      if (inname=="dir") {
        if ((file->IsDirectory() && !fname.Contains(".") && fname.Contains(regex_tag))) v_dirs.push_back(fname);
      } else  if(fname.Contains(regex_inname)) v_dirs.push_back(fname);
    }
  } // if(files)
  gSystem->cd(pwd); // The TSystemDirectory object seems to change current folder
  return v_dirs;
}

int change_branch_one(TString indir, TString name, TString outdir, vector<TString> var_type, vector<TString> var, 
		      vector<TString> var_val, TString newname){

  if(var_type.size()!=var.size() || var_type.size()!=var_val.size())
    { cout<<"[Change Branch One] Error: Branch vectors are not the same size"<<endl; exit(0); }

  const int nvar = var.size();
  bool verbose = false;

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
  TFile* oldfile = new TFile(indir+"/"+name);
  TTree* oldtree = static_cast<TTree*>(oldfile->Get("tree"));
  TTree* oldtreeglobal = static_cast<TTree*>(oldfile->Get("treeglobal"));
  
  vector<int> new_var_int_(nvar,-999);
  vector<float> new_var_flt_(nvar,-999);
  vector<double> new_var_dbl_(nvar,-999);
  deque<bool> new_var_bool_(nvar, false); // vector<bool>is not a vector in C++... so can't pass bools by reference
  vector<vector<int> * > new_var_vint_(nvar);
  vector<vector<float> * > new_var_vflt_(nvar);
  vector<vector<double> * > new_var_vdbl_(nvar);
  vector<vector<bool> * > new_var_vbool_(nvar);

  //Branches
  int nleps_=0;
  oldtree->SetBranchAddress("nleps",&nleps_);

  for(int ibch=0; ibch<nvar; ibch++){
    if(var_type[ibch]=="int")          {  new_var_int_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_int_[ibch]);   }
    else if(var_type[ibch]=="float")   {  new_var_flt_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_flt_[ibch]);   }
    else if(var_type[ibch]=="double")  {  new_var_dbl_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_dbl_[ibch]);   }
    else if(var_type[ibch]=="bool")    {  new_var_bool_[ibch] = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_bool_[ibch]);
      if(multiply[ibch])  cout<<"[Change Branch One] Warning: Multiplying branch of type \"bool\". Skipping branch."<<endl;         }
    else if(var_type[ibch]=="vint")    {  new_var_vint_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_vint_[ibch]); }
    else if(var_type[ibch]=="vfloat")  {  new_var_vflt_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_vflt_[ibch]); }
    else if(var_type[ibch]=="vdouble") {  new_var_vdbl_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_vdbl_[ibch]); }
    else if(var_type[ibch]=="vbool")   {  new_var_vbool_[ibch]  = 0;    oldtree->SetBranchAddress(var[ibch], &new_var_vbool_[ibch]); 
      if(multiply[ibch])  cout<<"[Change Branch One] Warning: Multiplying branch of type \"vbool\". Skipping branch."<<endl;        }

    else {cout<<"[Change Branch One] Error: Branch type invalid. Use only \"int\", \"float\", \"double\", \"bool\", \"vint\", \"vfloat\", \"vdouble\", or\"vbool\""<<endl; exit(0);}
  }

  //Set up new tree
  if(newname=="empty"){
    newname = name;
    newname.ReplaceAll(".root","_mod.root");
  }
  TFile* newfile = new TFile(outdir+newname,"recreate");
  TTree* newtree = oldtree->CloneTree(0);
  TTree* newtreeglobal = oldtreeglobal->CloneTree();

  //Loop and fill events with new weights
  int nentries = oldtree->GetEntries();
  for(int i=0; i<nentries; i++){
    oldtree->GetEntry(i);

    if(verbose && ((i<100&&i%10==0) || (i<1000&&i%100==0) || (i<10000&&i%1000==0) || (i%10000==0))){
      if(isatty(1)){
        printf("\r[Change Branch One] Processsing File: %i / %i (%i%%)",i,nentries,static_cast<int>((i*100./nentries)));
        fflush(stdout);
        if((i<100&&i+10>=nentries) || (i<1000&&i+100>=nentries) || (i<10000&&i+1000>=nentries) || (i>=10000&&i+10000>=nentries)) printf("\n");
      }
    }
    
    //Set vars
    for(int iset=0; iset<nvar; iset++){
      if(var[iset].Contains("_lep")&&nleps_!=0) continue; // For lepton scale factors
      if(!multiply[iset]){
        if(var_type[iset]=="int")             new_var_int_[iset]  = var_val[iset].Atoi();
        else if(var_type[iset]=="float")      new_var_flt_[iset]  = var_val[iset].Atof();
        else if(var_type[iset]=="double")     new_var_dbl_[iset]  = static_cast<double>(var_val[iset].Atof());
        else if(var_type[iset]=="bool")       new_var_bool_[iset] = var_val[iset].Atoi();
        else if(var_type[iset]=="vint")           
          for(unsigned int vidx=0; vidx<new_var_vint_[iset]->size(); vidx++)
            new_var_vint_[iset]->at(vidx) = var_val[iset].Atoi();
        else if(var_type[iset]=="vfloat")         
          for(unsigned int vidx=0; vidx<new_var_vflt_[iset]->size(); vidx++)
            new_var_vflt_[iset]->at(vidx) = var_val[iset].Atof();
        else if(var_type[iset]=="vdouble")        
          for(unsigned int vidx=0; vidx<new_var_vdbl_[iset]->size(); vidx++)
            new_var_vdbl_[iset]->at(vidx) = static_cast<double>(var_val[iset].Atof());
        else if(var_type[iset]=="vbool")          
          for(unsigned int vidx=0; vidx<new_var_vbool_[iset]->size(); vidx++)
            new_var_vbool_[iset]->at(vidx) = var_val[iset].Atoi();
      }
      else {
        if(var_type[iset]=="int")             new_var_int_[iset]  *= var_val[iset].Atoi(); 
        else if(var_type[iset]=="float")      new_var_flt_[iset]  *= var_val[iset].Atof(); 
        else if(var_type[iset]=="double")     new_var_dbl_[iset]  *= var_val[iset].Atof(); 
        else if(var_type[iset]=="vint")           
          for(unsigned int vidx=0; vidx<new_var_vint_[iset]->size(); vidx++)
            new_var_vint_[iset]->at(vidx) *= var_val[iset].Atoi();
        else if(var_type[iset]=="vfloat")         
          for(unsigned int vidx=0; vidx<new_var_vflt_[iset]->size(); vidx++)
            new_var_vflt_[iset]->at(vidx) *= var_val[iset].Atof();
        else if(var_type[iset]=="vdouble")        
          for(unsigned int vidx=0; vidx<new_var_vdbl_[iset]->size(); vidx++)
            new_var_vdbl_[iset]->at(vidx) *= static_cast<double>(var_val[iset].Atof());
      }
    }
    newtree->Fill();
  }
  
  //Save tree
  newtree->AutoSave();
  newtreeglobal->AutoSave();
  delete oldfile;
  delete newfile;
  return nentries;
}

bool eigen2x2(float matrix[2][2], float &eig1, float &eig2){
  float root = pow(matrix[0][0],2) + pow(matrix[1][1],2)-2*matrix[0][0]*matrix[1][1]+4*matrix[0][1]*matrix[1][0];
  if(root<0) return false;

  eig1 = (matrix[0][0]+matrix[1][1]+sqrt(root))/2.;
  eig2 = (matrix[0][0]+matrix[1][1]-sqrt(root))/2.;
  return true;
}

bool id_big2small(const int_double& left, const int_double& right){
  return left.second > right.second;
}

bool dd_small2big(const double_double& left, const double_double& right){
  return left.first < right.first;
}

bool dd_big2small(const double_double& left, const double_double& right){
  return left.first > right.first;
}

long double DeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(fabs(phi2-phi1), 2.L*PI);
  return dphi>PI ? 2.L*PI-dphi : dphi;
}

long double SignedDeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(phi2-phi1, 2.L*PI);
  if(dphi>PI){
    return dphi-2.L*PI;
  }else if(dphi<-PI){
    return dphi+2.L*PI;
  }else{
    return dphi;
  }
}

float dR(float eta1, float eta2, float phi1, float phi2) {
  return AddInQuadrature(eta1-eta2, DeltaPhi(phi1,phi2));
}

TString roundNumber(double num, int decimals, double denom){
  if(denom==0) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }

  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  return result;
}

TString addCommas(double num){
  TString result(""); result += num;
  int posdot(result.First('.'));
  if(posdot==-1) posdot = result.Length();
  for(int ind(posdot-3); ind > 0; ind -= 3)
    result.Insert(ind, ",");
  return result;
}

long double AddInQuadrature(long double x, long double y){
  if(fabs(y)>fabs(x)){
    const long double temp = y;
    y=x;
    x=temp;
  }
  if(x==0.) return y;
  const long double rat=y/x;
  return fabs(x)*sqrt(1.0L+rat*rat);
}

long double GetMass(long double e, long double px, long double py, long double pz){
  px/=e; py/=e; pz/=e;
  return fabs(e)*sqrt(1.0L-px*px-py*py-pz*pz);
}

long double GetMT(long double m1, long double pt1, long double phi1,
                  long double m2, long double pt2, long double phi2){
  return sqrt(m1*m1+m2*m2+2.L*(sqrt((m1*m1+pt1*pt1)*(m2*m2+pt2*pt2))-pt1*pt2*cos(phi2-phi1)));
}

long double GetMT(long double pt1, long double phi1,
                  long double pt2, long double phi2){
  //Faster calculation in massless case
  return sqrt(2.L*pt1*pt2*(1.L-cos(phi2-phi1)));
}

bool Contains(const string& text, const string& pattern){
  return text.find(pattern) != string::npos;
}

vector<string> Tokenize(const string& input,
                        const string& tokens){
  char* ipt(new char[input.size()+1]);
  memcpy(ipt, input.data(), input.size());
  ipt[input.size()]=static_cast<char>(0);
  char* ptr(strtok(ipt, tokens.c_str()));
  vector<string> output(0);
  while(ptr!=NULL){
    output.push_back(ptr);
    ptr=strtok(NULL, tokens.c_str());
  }
  return output;
}

void get_count_and_uncertainty(TTree& tree,
                               const string& cut,
                               double& count,
                               double& uncertainty){
  const string hist_name("temp");
  TH1D temp(hist_name.c_str(), "", 1, -1.0, 1.0);
  tree.Project(hist_name.c_str(), "0.0", cut.c_str());
  count=temp.IntegralAndError(0,2,uncertainty);
}

void AddPoint(TGraph& graph, const double x, const double y){
  graph.SetPoint(graph.GetN(), x, y);
}

string execute(const string &cmd){
  FILE *pipe = popen(cmd.c_str(), "r");
  if(!pipe) throw runtime_error("Could not open pipe.");
  const size_t buffer_size = 128;
  char buffer[buffer_size];
  string result = "";
  while(!feof(pipe)){
    if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
  }

  pclose(pipe);
  return result;
}

string RemoveTrailingNewlines(string str){
  while(!str.empty() && str.at(str.length()-1) == '\n'){
    str.erase(str.length()-1);
  }
  return str;
}

vector<double> LinearSpacing(size_t npts, double low, double high){
  vector<double> pts(npts,low+0.5*(high-low));
  if(npts>1){
    double gap = (high-low)/(npts-1.0);
    for(size_t pt = 0; pt < npts; ++pt){
      pts.at(pt) = low+pt*gap;
    }
  }
  return pts;
}

TString hoursMinSec(long seconds){
  int minutes((seconds/60)%60), hours(seconds/3600);
  TString hhmmss("");
  if(hours<10) hhmmss += "0";
  hhmmss += hours; hhmmss += ":";
  if(minutes<10) hhmmss += "0";
  hhmmss += minutes; hhmmss += ":";
  if((seconds%60)<10) hhmmss += "0";
  hhmmss += seconds%60; 

  return hhmmss;
}

void mergeNtuples(vector<TString> ntuples, TString outname){

  // Merging tree TTrees
  TChain chain("tree");
  for(size_t ind=0; ind<ntuples.size(); ind++) chain.Add(ntuples[ind]);
  chain.Merge(outname);

  // Merging treeglobal TTrees
  TChain chaing("treeglobal");
  for(size_t ind=0; ind<ntuples.size(); ind++) chaing.Add(ntuples[ind]);
  TTree *tglobal = chaing.CopyTree("1");
  tglobal->SetDirectory(0);
  TFile rootfile(outname, "UPDATE");
  rootfile.cd();
  tglobal->Write();
  rootfile.Close();

}
