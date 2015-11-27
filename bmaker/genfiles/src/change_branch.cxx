#include <iostream>
#include "TFile.h"
#include "TRegexp.h"
#include "TString.h"
#include "TTree.h"

#include "utilities.hh"

using namespace std;

int main(int argc, char *argv[]){

  if(argc<5){
    cout<<"Format: ./run/change_branch.exe <infolder> <outfolder> <sample> <branch_type1> <branch_name1> <branch_value1> ... <branch_valueN>"<<endl;
    cout<<"Accepts up to N branches folowing the format above."<<endl;
    cout<<"<sample> is an indentifier for files in the folder, e.g \"*.root\" or \"*TTJets*SingleLep*\""<<endl;
    cout<<"<branch_type> must be \"int\", \"float\", \"double\", \"bool\", \"vint\", \"vfloat\", \"vdouble\", or \"vbool\""<<endl;
    cout<<"Accepts multiplication of branches, i.e. <branch_value>=\"*3.1415\" or <branch_value>=\"3.1415*\""<<endl;
    return 1;
  }

  TString folder(argv[1]), outfolder(argv[2]), sample(argv[3]);
  if(!folder.EndsWith("/")) folder.Append("/");
  if(!outfolder.EndsWith("/")) outfolder.Append("/");

  vector<TString> var_type, var, var_val;
  for(int ivar=4; ivar<argc; ivar+=3){
    var_type.push_back(argv[ivar]);
    var.push_back(argv[ivar+1]);
    var_val.push_back(argv[ivar+2]);
  }

  vector<TString> files = dirlist(folder,".root");
  TRegexp regex(sample,true);
  
  for(unsigned int i=0; i<files.size(); i++){
    if(files[i].Contains(regex)){
      cout<<"[Change Branch] File #"<<i+1<<endl;
      change_branch_one(folder, files[i], outfolder, var_type, var, var_val);
    }
  }
}
