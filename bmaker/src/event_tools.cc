//// PHYS_OBJECTS: Common physics objects definitions
//// Function names follow the first-lowercase, following words-uppercase. No underscores

// System include files
#include <algorithm>
#include <iostream>
#include <fstream>

// FW include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

// User include files
#include "babymaker/bmaker/interface/event_tools.hh"
#include "babymaker/bmaker/interface/utilities.hh"

using namespace std;
using namespace utilities;

bool event_tools::isInJSON(string type, int run, int lumiblock){
  //if(type=="golden") return inJSON(VRunLumi2015golden, run, lumiblock);
  if(type=="golden") return true; // Applying this in bmaker_*_cfg.py
  if(type=="nonblind") return inJSON(VRunLumi2015nonblind, run, lumiblock);

  return true;
}

bool event_tools::hasGoodPV(edm::Handle<reco::VertexCollection> vtx){
  bool one_good_pv(false);
  for(unsigned ipv(0); ipv < vtx->size(); ipv++){
    const double pv_rho(sqrt(pow(vtx->at(ipv).x(),2) + pow(vtx->at(ipv).y(),2)));
    if(vtx->at(ipv).ndof()>4 && fabs(vtx->at(ipv).z())<24. && pv_rho<2.0 && vtx->at(ipv).isFake()==0){
      one_good_pv = true;
      break;
    }
  } // Loop over vertices
  return one_good_pv;
}

bool event_tools::passBeamHalo(int run, int event){
  if(!doBeamHalo) return true;
  if(badBeamHaloEvents.find(run) == badBeamHaloEvents.end()) return true;
  if(badBeamHaloEvents[run].find(event) == badBeamHaloEvents[run].end()) return true;
  return false;
}

void event_tools::fillBeamHaloMap(string eventList){
  cout<<"BABYMAKER::event_tools: Reading CSC Beam Halo filter file "<<eventList.c_str()<<endl;
  ifstream file(eventList.c_str());
  TString run_s, event_s;
  int nlines(0), run, event;
  while(file){
    nlines++;
    file >> run_s;
    if(run_s.CountChar(':') < 2) break;
    event_s = run_s;
    run_s.Remove(run_s.Index(':'), run_s.Length());
    event_s.Remove(0, event_s.Index(':')+1); event_s.Remove(0, event_s.Index(':')+1);
    run = run_s.Atoi(); event = event_s.Atoi();
    if(badBeamHaloEvents.find(run) == badBeamHaloEvents.end()) badBeamHaloEvents[run] = set<int>(); // New run
    if(badBeamHaloEvents[run].find(event) == badBeamHaloEvents[run].end()) badBeamHaloEvents[run].insert(event);// New event
    //cout<<"Pushed run "<<run<<" and event "<<event<<endl;
  }

  for(map<int, set<int> >::const_iterator it = badBeamHaloEvents.begin(); it != badBeamHaloEvents.end(); ++it) {
    run = it->first;
    //cout << run << "  "<<endl;
  }



}

event_tools::event_tools(TString outname):
  VRunLumi2015nonblind(MakeVRunLumi("nonblind")){

  doBeamHalo = false;
  string command("printf ${CMSSW_BASE}/src/babymaker/data/csc_beamhalo_filter/eventlist_");
  string eventList = execute(command.c_str());
  vector<string> datasets = {"SingleElectron", "SingleMuon", "MET", "HTMHT", "JetHT", "DoubleMuon", "DoubleEG"};
  for(size_t ind(0); ind < datasets.size(); ind++){
    if(outname.Contains(datasets[ind].c_str())){
	doBeamHalo = true;
	eventList += datasets[ind];
	eventList += "_csc2015.txt";
	break;
      }
  } // Loop over datasets for which we have CSC Beam Halo event lists with bad events
  
  if(doBeamHalo) fillBeamHaloMap(eventList);
}
// const std::vector<std::vector<int> > VRunLumi2015golden(MakeVRunLumi("golden"));

event_tools::~event_tools(){
}
