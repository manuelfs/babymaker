// baby_basic_: basic_ version of baby to handle reduce tree ntuples
// File generated with generate_baby.exe

#ifndef H_BABY_BASIC_
#define H_BABY_BASIC_

#include <vector>
#include <string>

#include "TTree.h"
#include "TChain.h"

#include "babymaker/bmaker/interface/baby.hh"

class baby_basic_ : public baby{
public:
  baby_basic_(); // Constructor to create tree
  baby_basic_(const std::string &filename); // Constructor to read tree

  virtual void GetEntry(const long entry);

  virtual void Fill();

  virtual std::string Type() const;

  virtual ~baby_basic_();

  virtual int  const & nevents() const;
  virtual int  & nevents();

private:
  int  nevents_;
  TBranch *b_nevents_;
  mutable bool c_nevents_;
};

#endif
