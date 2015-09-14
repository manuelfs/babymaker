// baby_basic: basic version of baby to handle reduce tree ntuples
//File generated with generate_baby.exe

#include "babymaker/bmaker/interface/baby.hh"

#include "babymaker/bmaker/interface/baby_basic.hh"

#include <stdexcept>
#include <string>
#include <vector>

#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"

using namespace std;

baby_basic::baby_basic():
  baby(),
  njets_(0),
  b_njets_(tree_.Branch("njets", &njets_)),
  c_njets_(false){
}

baby_basic::baby_basic(const string &filename):
  baby(filename),
  njets_(0),
  b_njets_(NULL),
  c_njets_(false){
  chain_.SetBranchAddress("njets", &njets_, &b_njets_);
}

void baby_basic::Fill(){
  baby::Fill();
  //Resetting variables
  njets_ = static_cast<int >(bad_val_);
}

string baby_basic::Type() const{
  return "basic";
}

baby_basic::~baby_basic(){
}

void baby_basic::GetEntry(const long entry){
  baby::GetEntry(entry);

  c_njets_ = false;
}

int  const & baby_basic::njets() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_njets_ && b_njets_){
    b_njets_->GetEntry(entry_);
    c_njets_ = true;
  }
  return njets_;
}

int  & baby_basic::njets(){
  if(read_only_ && !c_njets_ && b_njets_){
    b_njets_->GetEntry(entry_);
    c_njets_ = true;
  }
  return njets_;
}

