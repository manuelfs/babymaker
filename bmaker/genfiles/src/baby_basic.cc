// baby_basic: basic version of baby_base to handle reduce tree ntuples
//File generated with generate_baby.exe

#include "baby_base.hh"

#include "baby_basic.hh"

#include <stdexcept>
#include <string>
#include <vector>

#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"

using namespace std;

baby_basic::baby_basic():
  baby_base(){
}

baby_basic::baby_basic(const string &filename):
  baby_base(filename){
}

void baby_basic::Fill(){
  baby_base::Fill();
  //Resetting variables
}

string baby_basic::Type() const{
  return "basic";
}

baby_basic::~baby_basic(){
}

void baby_basic::GetEntry(const long entry){
  baby_base::GetEntry(entry);

}

