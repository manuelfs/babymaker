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
  if(type=="json2p6") return inJSON(VRunLumi2016json2p6, run, lumiblock);

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

//// Clean up the FastSim MET events with unmatched genJets
//// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSRecommendationsICHEP16#Cleaning_up_of_fastsim_jets_from
bool event_tools::passFSMET(edm::Handle<pat::JetCollection> alljets, edm::Handle<edm::View <reco::GenJet> > genjets){
  for (size_t ijet(0); ijet < alljets->size(); ijet++) {
    const pat::Jet &jet = (*alljets)[ijet];
    if(fabs(jet.eta()) >= 2.5 || jet.pt()<=20 || jet.chargedHadronEnergyFraction()>=0.1) continue;
    
    bool matched = false;
    for (size_t gjet(0); gjet < genjets->size(); gjet++) {
      const reco::GenJet &genjet = (*genjets)[gjet];
      double dr(deltaR(jet, genjet));
      if(dr < 0.3) {
	matched = true;
	break;
      }
    } // Loop over genjets
    if(!matched) {
      // cout<<endl<<endl<<"GENJETs"<<endl;
      // for (size_t gjet(0); gjet < genjets->size(); gjet++) {
      // 	const pat::Jet &genjet = (*genjets)[gjet];
      // 	cout<<"("<<genjet.pt()<<", "<<genjet.eta()<<", "<<genjet.phi()<<")"<<endl;
      // }
      // cout<<"Jet unmatched ("<<jet.pt()<<", "<<jet.eta()<<", "<<jet.phi()<<")"<<endl;
      return false;
    }
  } // Loop over jets
  
  return true;
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


int event_tools::type(const string &name){
  int sample = -999, category = -9, bin = -99;
  if(contains(name, "Run201")){ sample = 0;
    if(contains(name, "SingleElectron")){ category = 0;
    }else if(contains(name, "SingleMuon")){ category = 1;
    }else if(contains(name, "DoubleEG")){ category = 2;
    }else if(contains(name, "DoubleMuon")){ category = 3;
    }else if(contains(name, "MET")){ category = 4;
    }else if(contains(name, "HTMHT")){ category = 5;
    }else if(contains(name, "JetHT")){ category = 6;
    }
    auto pos = name.find("Run201")+7;
    if(pos < name.size()
       && isalpha(name.at(pos))
       && pos-1 < name.size()
       && isdigit(name.at(pos-1))){
      int run = toupper(name.at(pos))-65;
      int year = name.at(pos-1)-48;
      bin = (year-5)*26+run+1;
      //2015A=1, ..., 2015D=4, ..., 2016A=27, 2016B=28, 2016C=29, ...
      if(bin > 99) bin = 0;//Sanity check
    }
  }else if((contains(name, "TTJets") || contains(name, "TT_")) && !contains(name, "TTTT_")){ sample = 1;
    if(contains(name, "TTJets_Tune")){ category = 0; bin = 0;
    }else if(contains(name, "SingleLept")){ category = 1; bin = 0;
      if(contains(name, "genMET-150")) bin = 1;
    }else if(contains(name, "DiLept")){ category = 2; bin = 0;
      if(contains(name, "genMET-150")) bin = 1;
    }else if(contains(name, "TTJets_HT")){ category = 3;
      if(contains(name, "HT-600to800")){ bin = 0;
      }else if(contains(name, "HT-800to1200")){ bin = 1;
      }else if(contains(name, "HT-1200to2500")){ bin = 2;
      }else if(contains(name, "HT-2500toInf")){ bin = 3;
      }
    }else if(contains(name, "TT_")){ category = 4; bin = 0;
    }else if(contains(name, "TTJets_Mtt")){ category = 5; bin = 0;
    }
  }else if(contains(name, "WJets") && !contains(name, "TTWJets")){ sample = 2;
    if(contains(name, "WJetsToLNu_Tune")){ category = 0; bin = 0;
    }else if(contains(name, "WJetsToLNu_HT")){ category = 1;
      if(contains(name, "HT-70To100")){ bin = 0;
      }else if(contains(name, "HT-100To200")){ bin = 1;
      }else if(contains(name, "HT-200To400")){ bin = 2;
      }else if(contains(name, "HT-400To600")){ bin = 3;
      }else if(contains(name, "HT-600To800")){ bin = 4;
      }else if(contains(name, "HT-800To1200")){ bin = 5;
      }else if(contains(name, "HT-1200To2500")){ bin = 6;
      }else if(contains(name, "HT-2500ToInf")){ bin = 7;
      }else if(contains(name, "HT-600ToInf")){ bin = 10;
      }
    }else if(contains(name, "WJetsToQQ_HT")){ category = 2;
      if(contains(name, "HT-600ToInf")){ bin = 0;
      }
    }
  }else if(contains(name, "ST_")){ sample = 3;
    if(contains(name, "ST_s-channel")){ category = 0; bin = 0;
    }else if(contains(name, "ST_t-channel_top")){ category = 1; bin = 0;
    }else if(contains(name, "ST_t-channel_antitop")){ category = 2; bin = 0;
    }else if(contains(name, "ST_tW_top")){ category = 3; bin = 0;
    }else if(contains(name, "ST_tW_antitop")){ category = 4; bin = 0;
    }
  }else if(contains(name, "TTWJets")){ sample = 4;
    if(contains(name, "TTWJetsToLNu")){ category = 0; bin = 0;
    }else if(contains(name, "TTWJetsToQQ")){ category = 1; bin = 0;
    }
  }else if(contains(name, "TTZ")){ sample = 5;
    if(contains(name, "TTZToLLNuNu")){ category = 0; bin = 0;
    }else if(contains(name, "TTZToQQ")){ category = 1; bin = 0;
    }
  }else if(contains(name, "DYJetsToLL")){ sample = 6;
    if(contains(name, "DYJetsToLL_M-50_Tune")){ category = 0; bin = 0;
    }else if(contains(name, "DYJetsToLL_M-50_HT")){ category = 1;
      if(contains(name, "HT-70to100")){ bin = 0;
      }else if(contains(name, "HT-100to200")){ bin = 1;
      }else if(contains(name, "HT-200to400")){ bin = 2;
      }else if(contains(name, "HT-400to600")){ bin = 3;
      }else if(contains(name, "HT-600to800")){ bin = 4;
      }else if(contains(name, "HT-800to1200")){ bin = 5;
      }else if(contains(name, "HT-1200to2500")){ bin = 6;
      }else if(contains(name, "HT-2500toInf")){ bin = 7;
      }else if(contains(name, "HT-600toInf")){ bin = 10;
      }
    }
  }else if(contains(name, "QCD")){ sample = 7;
    if(contains(name, "QCD_HT")){ category = 0;
      if(contains(name, "HT50to100")){ bin = 0;
      }else if(contains(name, "HT100to200")){ bin = 1;
      }else if(contains(name, "HT200to300")){ bin = 2;
      }else if(contains(name, "HT300to500")){ bin = 3;
      }else if(contains(name, "HT500to700")){ bin = 4;
      }else if(contains(name, "HT700to1000")){ bin = 5;
      }else if(contains(name, "HT1000to1500")){ bin = 6;
      }else if(contains(name, "HT1500to2000")){ bin = 7;
      }else if(contains(name, "HT2000toInf")){ bin = 8;
      }
    }else if(contains(name, "QCD_Pt")){ category = 1;
      if(contains(name, "5to10")){ bin = 0;
      }else if(contains(name, "10to15")){ bin = 1;
      }else if(contains(name, "15to30")){ bin = 2;
      }else if(contains(name, "30to50")){ bin = 3;
      }else if(contains(name, "50to80")){ bin = 4;
      }else if(contains(name, "80to120")){ bin = 5;
      }else if(contains(name, "120to170")){ bin = 6;
      }else if(contains(name, "170to300")){ bin = 7;
      }else if(contains(name, "300to470")){ bin = 8;
      }else if(contains(name, "470to600")){ bin = 9;
      }else if(contains(name, "600to800")){ bin = 10;
      }else if(contains(name, "800to1000")){ bin = 11;
      }else if(contains(name, "1000to1400")){ bin = 12;
      }else if(contains(name, "1400to1800")){ bin = 13;
      }else if(contains(name, "1800to2400")){ bin = 14;
      }else if(contains(name, "2400to3200")){ bin = 15;
      }else if(contains(name, "3200toInf")){ bin = 16;
      }
    }
  }else if(contains(name, "ZJets")){ sample = 8;
    if(contains(name, "ZJetsToNuNu")){ category = 0;
      if(contains(name, "HT-70To100")){ bin = 0;
      }else if(contains(name, "HT-100To200")){ bin = 1;
      }else if(contains(name, "HT-200To400")){ bin = 2;
      }else if(contains(name, "HT-400To600")){ bin = 3;
      }else if(contains(name, "HT-600To800")){ bin = 4;
      }else if(contains(name, "HT-800To1200")){ bin = 5;
      }else if(contains(name, "HT-1200To2500")){ bin = 6;
      }else if(contains(name, "HT-2500ToInf")){ bin = 7;
      }else if(contains(name, "HT-600ToInf")){ bin = 10;
      }
    }else if(contains(name, "ZJetsToQQ")){ category = 1;
      if(contains(name, "HT600toInf")){ bin = 0;
      }
    }
  }else if(contains(name, "ttH")){ sample = 9;
    if(contains(name, "ttHJetTobb")){ category = 0; bin = 0;
    }
  }else if(contains(name, "TTGJets")){ sample = 10;
    if(contains(name, "TTGJets_Tune")){ category = 0; bin = 0;
    }
  }else if(contains(name, "TTTT")){ sample = 11;
    if(contains(name, "TTTT_Tune")){ category = 0; bin = 0;
    }
  }else if(contains(name, "WH_") && !contains(name,"TChiWH")){ sample = 12;
    if(contains(name, "WH_HToBB_WToLNu")){ category = 0; bin = 0;
    }
  }else if(contains(name, "ZH_")){ sample = 13;
    if(contains(name, "ZH_HToBB_ZToNuNu")){ category = 0; bin = 0;
    }
  }else if(contains(name, "WW") && !contains(name,"TChiHH")){ sample = 14;
    if(contains(name, "WWToLNuQQ")){ category = 0; bin = 0;
    }else if(contains(name, "WWTo2L2Nu")){ category = 1; bin = 0;
    }
  }else if(contains(name, "WZ") && !contains(name,"TChiHH")){ sample = 15;
    if(contains(name, "WZTo1L3Nu")){ category = 0; bin = 0;
    }else if(contains(name, "WZTo1L1Nu2Q")){ category = 1; bin = 0;
    }else if(contains(name, "WZTo2L2Q")){ category = 2; bin = 0;
    }else if(contains(name, "WZTo3LNu")){ category = 3; bin = 0;
    }
  }else if(contains(name, "ZZ") && !contains(name,"TChiHH")){ sample = 16;
    if(contains(name, "ZZ_Tune")){ category = 0; bin = 0;
    }
  }else if(contains(name, "T1tttt")){ sample = 100; category = 0; bin = 0;
  }else if(contains(name, "T2tt")){ sample = 101; category = 0; bin = 0;
  }else if(contains(name, "T1bbbb")){ sample = 102; category = 0; bin = 0;
  }else if(contains(name, "T2bb")){ sample = 103; category = 0; bin = 0;
  }else if(contains(name, "T1qqqq")){ sample = 104; category = 0; bin = 0;
  }else if(contains(name, "RPV")){ sample = 105; category = 0; bin = 0;
  }else if(contains(name, "TChiHH")){ sample = 106; category = 0; bin = 0;
  }else if(contains(name, "TChiWH")){ sample = 107; category = 0; bin = 0;
  }

  if(sample < 0 || category < 0 || bin < 0
     || category > 9 || bin > 99){
    DBG("Could not find type code for " << name << ": sample=" << sample << ", category=" << category << ", bin=" << bin);
    int code = -(1000*abs(sample)+100*abs(category)+abs(bin));
    if(code >= 0) code = -99999;
    return code;
  }else{
    int code = 1000*sample+100*category+bin;
    if(code < 0 || code > 107000){
      DBG("Type code out of range for " << name << ": sample=" << sample << ", category=" << category << ", bin=" << bin);
    }
    return code;
  }
}

event_tools::event_tools(TString outname){

  doBeamHalo = outname.Contains("Run2015");
  if(doBeamHalo) { 
    string command("printf ${CMSSW_BASE}/src/babymaker/data/csc_beamhalo_filter/csc2015_ee4sc_Jan13.txt");
    string eventList = execute(command.c_str());
    fillBeamHaloMap(eventList);
  }
}
// const std::vector<std::vector<int> > VRunLumi2015golden(MakeVRunLumi("golden"));

event_tools::~event_tools(){
}
