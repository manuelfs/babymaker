// baby_basic_: basic_ version of baby to handle reduce tree ntuples
//File generated with generate_baby.exe

#include "babymaker/bmaker/interface/baby.hh"

#include "babymaker/bmaker/interface/baby_basic_.hh"

#include <stdexcept>
#include <string>
#include <vector>

#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"

using namespace std;

baby_basic_::baby_basic_():
  baby(),
  nevents_(0),
  b_nevents_(tree_.Branch("nevents", &nevents_)),
  c_nevents_(false){
}

baby_basic_::baby_basic_(const string &filename):
  baby(filename),
  nevents_(0),
  b_nevents_(NULL),
  c_nevents_(false){
  chain_.SetBranchAddress("nevents", &nevents_, &b_nevents_);
}

void baby_basic_::Fill(){
  baby::Fill();
  //Resetting variables
  nevents_ = static_cast<int >(bad_val_);
}

string baby_basic_::Type() const{
  return "basic_";
}

baby_basic_::~baby_basic_(){
}

void baby_basic_::GetEntry(const long entry){
  baby::GetEntry(entry);

  c_nevents_ = false;
}

int  const & baby_basic_::nevents() const{
  if(!read_only_){
    throw std::logic_error("Trying to write to const tree.");
  }
  if(!c_nevents_ && b_nevents_){
    b_nevents_->GetEntry(entry_);
    c_nevents_ = true;
  }
  return nevents_;
}

int  & baby_basic_::nevents(){
  if(read_only_ && !c_nevents_ && b_nevents_){
    b_nevents_->GetEntry(entry_);
    c_nevents_ = true;
  }
  return nevents_;
}

