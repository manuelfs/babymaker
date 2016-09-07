// Event level tools

#ifndef H_EVENT_TOOLS
#define H_EVENT_TOOLS

// System include files
#include <string>
#include <map>
#include <set>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h" // For VertexCollection


#include "TString.h"

#include "babymaker/bmaker/interface/in_json.hh"

class event_tools{
public:
  bool isInJSON(std::string type, int run, int lumiblock);
  bool hasGoodPV(edm::Handle<reco::VertexCollection> vtx);
  bool passBeamHalo(int run, int event);
  bool passFSMET(edm::Handle<pat::JetCollection> alljets, edm::Handle<edm::View <reco::GenJet> > genjets);
  void fillBeamHaloMap(std::string eventList);

  const std::vector<std::vector<int> > VRunLumi2015nonblind;
  const std::vector<std::vector<int> > VRunLumi2016json2p6;

  std::map<int, std::set<int> > badBeamHaloEvents;
  bool doBeamHalo;

  static int type(const std::string &name);

  event_tools(TString outname);
  ~event_tools();
};

#endif
