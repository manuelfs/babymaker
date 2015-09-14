// baby: base class to handle reduced tree ntuples
// File generated with generate_baby.exe

#ifndef H_BABY
#define H_BABY

#include <vector>
#include <string>

#include "TTree.h"
#include "TChain.h"
class baby{
public:
  baby(); // Constructor to create tree
  baby(const std::string &filename); // Constructor to read tree

  int Add(const std::string &filename);
  long GetEntries() const;
  virtual void GetEntry(const long entry);
  bool PassString(TString cut);

  virtual void Fill();
  void Write();

  virtual std::string Type() const;

  static const double bad_val_;

  virtual ~baby();


  __attribute__((noreturn)) virtual int  const & nevents() const;
  __attribute__((noreturn)) virtual int  & nevents();
  __attribute__((noreturn)) virtual int  const & njets() const;
  __attribute__((noreturn)) virtual int  & njets();

protected:
  TChain chain_;
  TTree tree_;
  long entry_;
  const bool read_only_;

private:
  class VectorLoader{
  public:
    VectorLoader();
  private:
    static bool loaded_;
  };

  static VectorLoader vl_;
};

baby* NewTree(const std::type_info &type);

#include "babymaker/bmaker/interface/baby_basic.hh"
#include "babymaker/bmaker/interface/baby_basic_.hh"

#endif
