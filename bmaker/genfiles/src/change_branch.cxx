#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "utilities.hh"

#include "change_branch.hh"

using namespace std;

int main(int argc, char *argv[]){

  if(argc<5){
    cout<<"Format: ./run/change_branch.exe <infolder> <outfolder> <branch type> <branch name> <branch value>"<<endl;
    cout<<"<branch type> must be \"int\", \"float\", \"double\", or \"bool\""<<endl
    cout<<"Accepts multiplication of branches, i.e. <branch value>=\"*3.1415\" or <branch value>=\"3.1415*\""<<endl;
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
