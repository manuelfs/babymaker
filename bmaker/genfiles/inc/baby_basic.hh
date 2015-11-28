// baby_basic: basic version of baby_base to handle reduce tree ntuples
// File generated with generate_baby.exe

#ifndef H_BABY_BASIC
#define H_BABY_BASIC

#include <vector>
#include <string>

#include "TTree.h"
#include "TChain.h"

#include "baby_base.hh"

class baby_basic : public baby_base{
public:
  baby_basic(); // Constructor to create tree
  baby_basic(const std::string &filename); // Constructor to read tree

  virtual void GetEntry(const long entry);

  virtual void Fill();

  virtual std::string Type() const;

  virtual ~baby_basic();


private:
};

#endif
