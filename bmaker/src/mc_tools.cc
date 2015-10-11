// MC_TOOLS: Functions that deal with MC particles

#include <iostream>
#include "babymaker/bmaker/interface/mc_tools.hh"

using namespace std;
using namespace utilities;

// Checks if "id" is one of the immediate daughters of "mc"
bool mc_tools::hasDaughter(const reco::GenParticle &mc, size_t id){
  for(size_t idau(0); idau < mc.numberOfDaughters(); idau++) 
    if(abs(mc.daughter(idau)->pdgId())==id) return true;

  return false;
}

// Checks if "mc" is "id" and does not decay to itself
bool mc_tools::isLast(const reco::GenParticle &mc, size_t id){
  return (abs(mc.pdgId())==id && !hasDaughter(mc, id));
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

void mc_tools::printParticle(const reco::GenParticle &mc){
  int momid(0);
  if(mc.mother()) momid = mc.mother()->pdgId();
  cout<<setw(8)<<mc.pdgId()<<",  mom "<<setw(8)<<momid<<",  status "<<setw(3)<<mc.status()
      <<",  (pt,eta,phi) = ("<<setw(7)<<roundNumber(mc.pt(),2)<<", "
      <<setw(5)<<roundNumber(mc.eta(),2)<<", "<<setw(5)<<roundNumber(mc.phi(),2)
      <<"),  nDau "<<setw(2)<<mc.numberOfDaughters()<<": ";
  for(size_t idau(0); idau < mc.numberOfDaughters(); idau++) {
    cout<<mc.daughter(idau)->pdgId();
    if(idau < mc.numberOfDaughters()-1) cout<<", ";
  }
  cout<<endl;

}

mc_tools::mc_tools(){
}

mc_tools::~mc_tools(){ 
}

