// MC_TOOLS: Functions that deal with MC particles

#include <iostream>
#include "babymaker/bmaker/interface/mc_tools.hh"
#include "TRegexp.h"

using namespace std;
using namespace utilities;

// Checks if "id" is one of the immediate daughters of "mc"
bool mc_tools::hasDaughter(const reco::GenParticle &mc, size_t id){
  for(size_t idau(0); idau < mc.numberOfDaughters(); idau++) 
    if(abs(mc.daughter(idau)->pdgId())==id) return true;

  return false;
}


int mc_tools::numChargeDaughters(const reco::GenParticle &mc){
  int ncharged=0;
  for(size_t idau(0); idau < mc.numberOfDaughters(); idau++) 
    if(abs(mc.daughter(idau)->charge())!=0) ncharged++;

  return ncharged;
}


// Checks if "mc" is "id" and does not decay to itself
bool mc_tools::isLast(const reco::GenParticle &mc, size_t id){
  return (abs(mc.pdgId())==id && !hasDaughter(mc, id));
}
// Checks if "mc" comes from a W or a tau from a W
bool mc_tools::fromWOrWTau(const reco::GenParticle &mc){
  const reco::GenParticle *mcMom;
  int momId = abs(mom(mc, mcMom));
  if(momId == 24) return true;
  if(momId == 15 && abs(mom(*mcMom, mcMom)) == 24) return true;
  return false;
}


bool mc_tools::fromTau(const reco::GenParticle &mc){
  const reco::GenParticle *mcMom;
  int momId = abs(mom(mc, mcMom));
  if(momId == 15 && abs(mom(*mcMom, mcMom)) == 24) return true;
  return false;
}


// Returns the id and the pointer to the mother of "mc" not being itself
int mc_tools::mom(const reco::GenParticle &mc, const reco::GenParticle *&mcMom){
  mcMom = static_cast<const reco::GenParticle *>(mc.mother());
  if(mcMom){
    if(mcMom->pdgId() == mc.pdgId()) return mom(*mcMom, mcMom);
    else return mcMom->pdgId();
  } else return 0;
}

// Returns the index of the mother in the "indices" list
int mc_tools::getMomIndex(const reco::GenParticle &mc, vector<pair<int, const reco::GenParticle *> > indices){
  const reco::GenParticle * mcMom = static_cast<const reco::GenParticle *>(mc.mother());
  if(mcMom){
    for(size_t ind=0; ind<indices.size(); ind++)
      if(indices[ind].second == mcMom) return static_cast<int>(ind);
    return getMomIndex(*mcMom, indices);
  } else return -1;
}

// Checks if "mc" eventually decays to "id" (after it stops decays to itself)
bool mc_tools::decaysTo(const reco::GenParticle &mc, size_t id, const reco::GenParticle *&mcDau){
  bool idInDaughters(false);
  for(size_t idau(0); idau < mc.numberOfDaughters(); idau++) {
    size_t dauid(abs(mc.daughter(idau)->pdgId()));
    if(dauid == id) {
      idInDaughters = true;
      mcDau = static_cast<const reco::GenParticle *>(mc.daughter(idau));
    }
    if(dauid == abs(mc.pdgId())) return decaysTo(*static_cast<const reco::GenParticle*>((mc.daughter(idau))), id, mcDau);
  } // Loop over daughters
  return idInDaughters;
}
// Code for identifying gluon splitting based on 
// https://raw.githubusercontent.com/cms-btv-pog/RecoBTag-PerformanceMeasurements/7_4_X/plugins/BTagAnalyzer.cc
bool mc_tools::isFromGSP(const reco::Candidate *c) 
{
  // for example, in case a matched genParton does not exist
  if(!c) return false;

  bool isFromGSP = false;
  if( c->numberOfMothers() == 1 ) {
    const reco::Candidate* dau = c;
    const reco::Candidate* mom = c->mother();
    while( dau->numberOfMothers() == 1 && !( isHardProcess(mom->status()) && (abs(mom->pdgId())==4 || abs(mom->pdgId())==5) ) ) {
      if( abs(mom->pdgId())==21 && (abs(c->pdgId())==4 || abs(c->pdgId())==5) )
	{
	  isFromGSP = true;
	  break;
	}
      dau = mom;
      mom = dau->mother();
    }
  }

  return isFromGSP;
}

bool mc_tools::isHardProcess(const int status)
{
  // Assume Pythia8 status codes
  if( status>=21 && status<=29 ) return true;
  else return false;

}

void mc_tools::printParticle(const reco::GenParticle &mc){
  int momid(0);
  if(mc.mother()) momid = mc.mother()->pdgId();
  cout<<setw(8)<<mc.pdgId()<<",  mom "<<setw(8)<<momid<<",  status "<<setw(3)<<mc.status()
      <<",  (pt,eta,phi) = ("<<setw(7)<<roundNumber(mc.pt(),2)<<", "
      <<setw(5)<<roundNumber(mc.eta(),2)<<", "<<setw(5)<<roundNumber(mc.phi(),2)
      <<"), mass "<<setw(6)<<roundNumber(mc.mass(),1)
      <<",  nDau "<<setw(2)<<mc.numberOfDaughters()<<": ";
  for(size_t idau(0); idau < mc.numberOfDaughters(); idau++) {
    cout<<mc.daughter(idau)->pdgId();
    if(idau < mc.numberOfDaughters()-1) cout<<", ";
  }
  cout<<endl;

}

void mc_tools::getMassPoints(TString mpoints, int &mgluino, int &mlsp){

  TString ml(mpoints), mg(mpoints);
  if (mpoints.Contains("TChiHH")){
    mg.Remove(0, ml.Last('_')+1);
    mg.ReplaceAll("\n","");
    ml = "1";
    if(!mg.IsFloat()){
      cout<<"ERROR: Improper parsing of mass points"<<endl;
      cout<<"Mass of gluino parsed as \""<<mg<<"\""<<endl;
      exit(0);
    }
  } else {
    ml.Remove(0, ml.Last('_')+1);
    ml.ReplaceAll("\n","");
    mg.Remove(mg.Last('_'));
    mg.Remove(0, mg.Last('_')+1);
    
    if(!mg.IsFloat() || !ml.IsFloat()){
      cout<<"ERROR: Improper parsing of mass points"<<endl;
      cout<<"Mass of gluino parsed as \""<<mg<<"\""<<endl;
      cout<<"Mass of LSP parsed as \""<<ml<<"\""<<endl;
      exit(0);
    }
  }

  mgluino = mg.Atoi();
  mlsp = ml.Atoi();
}

mc_tools::mc_tools(){
}

mc_tools::~mc_tools(){ 
}

