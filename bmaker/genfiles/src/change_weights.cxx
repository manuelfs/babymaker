#include <iostream>
#include "TChain.h"
#include "TRegexp.h"
#include <string>
#include <vector>

#include "change_branch.hh"
#include "utilities.hh"

using namespace std;

int main(int argc, char *argv[]){

  if(argc<1){
    cout<<"Format: ./run/change_weights.exe <infolder>=. <sample>=\"*.root\" <outfolder>=infolder"<<endl;
    return 1;
  }
  
  TString folder("."), sample("*.root"), outfolder(folder);
  if(argc>=2) folder=argv[1]; 
  if(argc>=3) sample=argv[2]; 
  if(argc>=4) outfolder=argv[3]; 
  if(!folder.EndsWith("/")) folder.Append("/");
  if(!outfolder.EndsWith("/")) outfolder.Append("/");

  TChain ch("tree"); ch.Add(folder+sample);
  TChain gl("treeglobal"); gl.Add(folder+sample);
  
  vector<TString> var_type, var, var_val;

  double xsec = gl.GetMaximum("xsec");
  int nentries = ch.GetEntries();
  int npos = ch.GetEntries("weight>0");
  int nneg = ch.GetEntries("weight<0");
  double weight = ch.GetMaximum("weight");

  int nentries_eff = npos-nneg;
  double weight_eff = xsec*1000/nentries_eff;

  double weight_corr = weight_eff/weight;

  cout<<"[Change Weights] Xsec: "<<xsec<<", Number of entries: "<<nentries<<", Weight: "<<weight<<endl;
  cout<<"[Change Weights] Xsec: "<<xsec<<", Effective number of entries: "<<nentries_eff<<", Effective weight: "<<weight_eff<<endl;

  var_type.push_back("float");
  var.push_back("weight");
  var_val.push_back("*"+to_string(weight_corr));

  vector<TString> files = dirlist(folder,".root");
  TRegexp regex(sample,true);
  
  for(unsigned int i=0; i<files.size(); i++){
    if(files[i].Contains(regex))
      change_branch_one(folder, files[i], outfolder, var_type, var, var_val);
  }
}
