// baby: base class to handle reduce tree ntuples
//File generated with generate_baby.exe

#include "babymaker/bmaker/interface/baby.hh"

#include <stdexcept>
#include <string>
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
using namespace std;

bool baby::VectorLoader::loaded_ = false;

baby::VectorLoader baby::vl_ = baby::VectorLoader();

baby::VectorLoader::VectorLoader(){
  if(!loaded_){
    gROOT->ProcessLine("#include <vector>");
    loaded_ = true;
  }
}

const double baby::bad_val_ = -999.;

baby::baby():
  chain_("junk", "junk"),
  tree_("tree", "tree"),
  entry_(0),
  read_only_(false){
}

baby::baby(const string &filename):
  chain_("tree","tree"),
  tree_("junk","junk"),
  entry_(0),
  read_only_(true){
  chain_.Add(filename.c_str());
}

void baby::Fill(){
  if(read_only_){
    throw std::logic_error("Trying to write to read-only tree");
  }else{
    tree_.Fill();
  }

  //Resetting variables
}

void baby::Write(){
  if(read_only_){
    throw std::logic_error("Trying to write to read-only tree.");
  }else{
    tree_.Write();
  }
}

string baby::Type() const{
  return "";
}

baby::~baby(){
}

int baby::Add(const std::string &filename){
  if(!read_only_){
    throw std::logic_error("Trying to add files to tree opened for writing.");
  }
  return chain_.Add(filename.c_str());
}

bool baby::PassString(TString cut){
 bool result = true;
 return result;
}

long baby::GetEntries() const{
  if(read_only_){
    return chain_.GetEntries();
  }else{
    return tree_.GetEntries();
  }
}

void baby::GetEntry(const long entry){
  if(!read_only_){
    throw std::logic_error("Trying to read from write-only tree.");
  }

  entry_ = chain_.LoadTree(entry);
}

int  const & baby::nevents() const{
  throw std::logic_error("nevents does not exist in this baby version.");
}

int  const & baby::njets() const{
  throw std::logic_error("njets does not exist in this baby version.");
}

int  & baby::nevents(){
  throw std::logic_error("nevents does not exist in this baby version.");
}

int  & baby::njets(){
  throw std::logic_error("njets does not exist in this baby version.");
}

#include "babymaker/bmaker/interface/baby_basic.hh"
#include "babymaker/bmaker/interface/baby_basic_.hh"
baby* NewTree(const std::type_info &type){

  if(type == typeid(baby)) return new baby;
  else if(type == typeid(baby_basic)) return static_cast<baby*>(new baby_basic);
  else if(type == typeid(baby_basic_)) return static_cast<baby*>(new baby_basic_);
  else return new baby;
}

