#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "utilities.hh"

using namespace std;

int main(int argc, char *argv[]){

  if(argc<5){
    cout<<"Format: ./run/change_branch.exe <infolder> <outfolder> <branch_type1> <branch_name1> <branch_value1> ... <branch_valueN>"<<endl;
    cout<<"Accepts up to N branches folowing the format above."<<endl;
    cout<<"Accepts multiplication of branches, i.e. <branch_value>=\"*3.1415\" or <branch_value>=\"3.1415*\""<<endl;
    cout<<"<branch_type> must be \"int\", \"float\", \"double\", or \"bool\""<<endl;
    return 1;
  }

  TString folder(argv[1]), outfolder(argv[2]);
  if(!folder.EndsWith("/")) folder.Append("/");
  if(!outfolder.EndsWith("/")) outfolder.Append("/");

  vector<TString> var_type, var, var_val;
  for(int ivar=3; ivar<argc; ivar+=3){
    var_type.push_back(argv[ivar]);
    var.push_back(argv[ivar+1]);
    var_val.push_back(argv[ivar+2]);
  }

  vector<TString> files = dirlist(folder,".root");
  
  for(unsigned int i=0; i<files.size(); i++){
    change_branch_one(folder, files[i], outfolder, var_type, var, var_val);
  }
}
