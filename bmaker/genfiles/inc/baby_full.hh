// baby_full: full version of baby_base to handle reduce tree ntuples
// File generated with generate_baby.exe

#ifndef H_BABY_FULL
#define H_BABY_FULL

#include <vector>
#include <string>

#include "TTree.h"
#include "TChain.h"

#include "baby_base.hh"

class baby_full : public baby_base{
public:
  baby_full(); // Constructor to create tree
  baby_full(const std::string &filename); // Constructor to read tree

  virtual void GetEntry(const long entry);

  virtual void Fill();

  virtual std::string BabyType() const;

  virtual ~baby_full();


private:
};

#endif
