//----------------------------------------------------------------------------
// utilities - Various functions used accross the code
//----------------------------------------------------------------------------

#ifndef INT_ROOT
#include "utilities.hh"
#endif

#include <cmath>
#include <deque>
#include <iostream>
#include <string>
#include <stdexcept>

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

using namespace std;

// Returns list of directorites or files in folder
vector<TString> dirlist(const TString &folder,
                        const TString &inname,
                        const TString &tag){
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
        if ((file->IsDirectory() && !fname.Contains(".") && fname.Contains(tag))) v_dirs.push_back(fname);
      } else  if(fname.Contains(inname)) v_dirs.push_back(fname);
    }
  } // if(files)
  gSystem->cd(pwd); // The TSystemDirectory object seems to change current folder
  return v_dirs;
}

void change_branch_one(TString indir, TString name, TString outdir, vector<TString> var_type, vector<TString> var, vector<TString> var_val){

  if(var_type.size()!=var.size() || var_type.size()!=var_val.size())
    { cout<<"[Change Branch One] Error: Branch vectors are not the same size"<<endl; exit(0); }

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
  name.ReplaceAll(".root","_mod.root");
  TFile* newfile = new TFile(outdir+name,"recreate");
  TTree* newtree = oldtree->CloneTree(0);
  TTree* newtreeglobal = oldtreeglobal->CloneTree();

  //Loop and fill events with new weights
  int nentries = oldtree->GetEntries();
  for(int i=0; i<nentries; i++){
    oldtree->GetEntry(i);

    if((i<100&&i%10==0) || (i<1000&&i%100==0) || (i<10000&&i%1000==0) || (i%10000==0)){
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
}

void change_branch_one(TString indir, TString name, TString outdir, vector<TString> var_type, vector<TString> var,  vector<vector<TString> > var_val){

  if(var_type.size()!=var.size() || var_type.size()!=var_val.size())
    { cout<<"[Change Branch One] Error: Branch vectors are not the same size"<<endl; exit(0); }

  const int nvar = var.size();
  
  //Setup
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
  }

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

  //Branches
  int nleps_=0;
  oldtree->SetBranchAddress("nleps",&nleps_);

  for(int ibch=0; ibch<nvar; ibch++){
    if(var_type[ibch]=="int")          {  new_var_int_[ibch]   = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_int_[ibch]);    }
    else if(var_type[ibch]=="float")   {  new_var_flt_[ibch]   = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_flt_[ibch]);    }
    else if(var_type[ibch]=="double")  {  new_var_dbl_[ibch]   = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_dbl_[ibch]);    }
    else if(var_type[ibch]=="bool")    {  new_var_bool_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_bool_[ibch]);   }
    else if(var_type[ibch]=="vint")    {  new_var_vint_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_vint_[ibch]);   }
    else if(var_type[ibch]=="vfloat")  {  new_var_vflt_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_vflt_[ibch]);   }
    else if(var_type[ibch]=="vdouble") {  new_var_vdbl_[ibch]  = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_vdbl_[ibch]);   }
    else if(var_type[ibch]=="vbool")   {  new_var_vbool_[ibch] = 0;     oldtree->SetBranchAddress(var[ibch], &new_var_vbool_[ibch]);  }

    else {cout<<"[Change Branch One] Error: Branch type invalid. Use only \"int\", \"float\", \"double\", \"bool\", \"vint\", \"vfloat\", \"vdouble\", or\"vbool\""<<endl; exit(0);}
  }

  //Set up new tree
  name.ReplaceAll(".root","_mod.root");
  TFile* newfile = new TFile(outdir+name,"recreate");
  TTree* newtree = oldtree->CloneTree(0);
  TTree* newtreeglobal = oldtreeglobal->CloneTree();

  //Loop and fill events with new weights
  int nentries = oldtree->GetEntries();
  for(int i=0; i<nentries; i++){
    oldtree->GetEntry(i);

    if((i<100&&i%10==0) || (i<1000&&i%100==0) || (i<10000&&i%1000==0) || (i%10000==0) || (i+1==nentries)){
      if(isatty(1)){
        printf("\r[Change Branch One] Processsing File: %i / %i (%i%%)",i+1,nentries,static_cast<int>((i*100./nentries)));
        fflush(stdout);
        if((i<100&&i+10>=nentries) || (i<1000&&i+100>=nentries) || (i<10000&&i+1000>=nentries) || (i>=10000&&i+10000>=nentries)) printf("\n");
      }
    }
    
    //Set vars
    for(int iset=0; iset<nvar; iset++){
      if(var[iset].Contains("_lep")&&nleps_!=0) continue; // For lepton scale factors    
      for(unsigned int vidx=0; vidx<var_val[iset].size(); vidx++){
	if(!multiply[iset][vidx]){
	  if(var_type[iset]=="int")             new_var_int_[iset]  = var_val[iset][vidx].Atoi();
	  else if(var_type[iset]=="float")      new_var_flt_[iset]  = var_val[iset][vidx].Atof();
	  else if(var_type[iset]=="double")     new_var_dbl_[iset]  = static_cast<double>(var_val[iset][vidx].Atof());
	  else if(var_type[iset]=="bool")       new_var_bool_[iset] = var_val[iset][vidx].Atoi();
	  else if(var_type[iset]=="vint")       new_var_vint_[iset]->at(vidx) = var_val[iset][vidx].Atoi();
	  else if(var_type[iset]=="vfloat")     new_var_vflt_[iset]->at(vidx) = var_val[iset][vidx].Atof();
	  else if(var_type[iset]=="vdouble")    new_var_vdbl_[iset]->at(vidx) = static_cast<double>(var_val[iset][vidx].Atof());
	  else if(var_type[iset]=="vbool")      new_var_vbool_[iset]->at(vidx) = var_val[iset][vidx].Atoi();
	}
	else {
	  if(var_type[iset]=="int")             new_var_int_[iset]  *= var_val[iset][vidx].Atoi(); 
	  else if(var_type[iset]=="float")      new_var_flt_[iset]  *= var_val[iset][vidx].Atof(); 
	  else if(var_type[iset]=="double")     new_var_dbl_[iset]  *= var_val[iset][vidx].Atof(); 
	  else if(var_type[iset]=="vint")       new_var_vint_[iset]->at(vidx) *= var_val[iset][vidx].Atoi();
	  else if(var_type[iset]=="vfloat")     new_var_vflt_[iset]->at(vidx) *= var_val[iset][vidx].Atof();
	  else if(var_type[iset]=="vdouble")    new_var_vdbl_[iset]->at(vidx) *= static_cast<double>(var_val[iset][vidx].Atof());
	}
      }
    }
    newtree->Fill();
  }
  
  //Save tree
  newtree->AutoSave();
  newtreeglobal->AutoSave();
  delete oldfile;
  delete newfile;
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

