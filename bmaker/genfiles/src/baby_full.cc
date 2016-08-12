// baby_full: full version of baby_base to handle reduce tree ntuples
//File generated with generate_baby.exe

#include "baby_base.hh"

#include "baby_full.hh"

#include <stdexcept>
#include <string>
#include <vector>

#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"

using namespace std;

baby_full::baby_full():
  baby_base(){
}

baby_full::baby_full(const string &filename):
  baby_base(filename){
}

void baby_full::Fill(){
  baby_base::Fill();
  //Resetting variables
}

string baby_full::BabyType() const{
  return "full";
}

baby_full::~baby_full(){
}

void baby_full::GetEntry(const long entry){
  baby_base::GetEntry(entry);

}

