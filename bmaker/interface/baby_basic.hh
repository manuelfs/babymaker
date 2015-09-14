// baby_basic: basic version of baby to handle reduce tree ntuples
// File generated with generate_baby.exe

#ifndef H_BABY_BASIC
#define H_BABY_BASIC

#include <vector>
#include <string>

#include "TTree.h"
#include "TChain.h"

#include "babymaker/bmaker/interface/baby.hh"

class baby_basic : public baby{
public:
  baby_basic(); // Constructor to create tree
  baby_basic(const std::string &filename); // Constructor to read tree

  virtual void GetEntry(const long entry);

  virtual void Fill();

  virtual std::string Type() const;

  virtual ~baby_basic();

  virtual int  const & njets() const;
  virtual int  & njets();

private:
  int  njets_;
  TBranch *b_njets_;
  mutable bool c_njets_;
};

#endif
