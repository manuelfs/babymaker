// combine_datasets: Finds all unique events in a list of datasets

#include <ctime>

#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <unistd.h>  // getopt
#include <iomanip>   // setw

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"

#include "utilities.hh"

using namespace std; 

int main(int argc, char *argv[]){
  time_t startTime, curTime;
  time(&startTime);

  TString file_datasets("txt/singlelep.txt"), infolder(""), outfolder("out/");
  int begrun(-1), endrun(-1);
  int c(0);
  while((c=getopt(argc, argv, "f:i:o:b:e:"))!=-1){
    switch(c){
    case 'i':
      infolder=optarg;
      break;
    case 'b':
      begrun=atoi(optarg);
      break;
    case 'e':
      endrun=atoi(optarg);
      break;
    case 'o':
      outfolder=optarg;
      break;
    case 'f':
      file_datasets=optarg;
      break;
    default:
      break;
    }
  }
  if(file_datasets=="" || infolder==""){
    cout<<endl<<"Specify input folder and datasets: "
	<<"./run/combine_datasets.exe -i <infolder> -o <outfolder=out> -f <file_datasets=txt/singlelep.txt> -b  <begrun=-1> -e <endrun=-1>"<<endl<<endl;
    return 1;
  }

  TString run_s="_runs"; run_s += begrun; 
  if(endrun>begrun){
    run_s += "-"; run_s += endrun;
  }
  if(begrun>0){
    if(endrun<begrun){
      cout<<"You set begrun to "<<begrun<<", and endrun to "<<endrun
	  <<", but endrun has to be >= to begrun. Exiting"<<endl<<endl;
      return 1;
    }
    cout<<"Combining "<<run_s<<" of ntuples in "<<infolder<<endl;
  }

  vector<TString> datasets;
  TString buffer, basename("Run2016");
  ifstream indata(file_datasets);
  while(indata){
    indata >> buffer;
    if(buffer!=""){
      datasets.push_back(buffer);
      basename += ("_"+buffer);
    }
  }
  if(begrun>0) basename += run_s;

  map<int, map<int, set<Long64_t> > > runs;
  Long64_t event;
  int run, lumiblock;

  for(unsigned idata(0); idata < datasets.size(); idata++){
    TChain chain("tree"), treeglobal("treeglobal");
    TString filename(infolder+"/*"+datasets[idata]+"*.root");
    int files = chain.Add(filename);
    if(files<1) {
      cout<<"No files found for "<<filename<<endl;
      continue;
    }
    treeglobal.Add(filename);
    gSystem->mkdir(outfolder, kTRUE);
    TString outname(outfolder+"/baby_");
    outname += idata;
    outname += "_"+basename;
    outname += ".root";
    TFile outfile(outname, "RECREATE");
    outfile.cd();

    TTree *outtree(chain.CloneTree(0));

    // TBranch *b_event = chain.Branch("event", &event);
    // TBranch *b_run = chain.Branch("run", &run);
    TBranch *b_event(nullptr), *b_lumiblock(nullptr), *b_run(nullptr);
    chain.SetBranchAddress("event", &event, &b_event);
    chain.SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
    chain.SetBranchAddress("run", &run, &b_run);

    long entries(chain.GetEntries()), tree_entry;

    cout<<endl<<"Doing "<<files<<" files in "<<filename<<" with "<<entries<<" entries"<<endl;
    time(&startTime);
    for(int entry(0); entry<entries; entry++){
      if(entry!=0 && entry%250000==0) {
	time(&curTime);
	int seconds(difftime(curTime,startTime));
	
	cout<<"Doing entry "<<setw(10)<<addCommas(entry)<<" of "<<addCommas(entries)
	    <<"    Took "<<setw(6)<<seconds<<" seconds at "
	    <<setw(4)<<roundNumber(entry,1,seconds*1000.)<<" kHz"<<endl;
      }
      
      // Load "run" first, and check if it's in the range we care about
      tree_entry = chain.LoadTree(entry);
      b_run->GetEntry(tree_entry);
      if(begrun>0 && (run<begrun || run>endrun)) continue;
      b_lumiblock->GetEntry(tree_entry);
      b_event->GetEntry(tree_entry);

      if(runs.find(run) == runs.end()) runs.emplace(run, map<int, set<Long64_t> >{}); // New run
      auto &lumiblocks = runs.at(run);
      if(lumiblocks.find(lumiblock) == lumiblocks.end()) lumiblocks.emplace(lumiblock, set<Long64_t>{}); // New lumiblock
      auto &events = lumiblocks.at(lumiblock);
      if(events.find(event) == events.end()){ // New event
	events.emplace(event);
	// You need to load all branches to copy them into outtree
	chain.GetEntry(entry);
	outtree->Fill();
      } 
    } // Loop over entries
    outtree->Write();
    treeglobal.CloneTree(-1,"fast");
    outfile.Write();
    outfile.Close();
    time(&curTime);
    cout<<"Took "<<difftime(curTime,startTime) <<" seconds to write "<<outname<<endl;

  } // Loop over datasets

  // for(auto it = events.cbegin(); it != events.cend(); ++it) {
  //   cout << it->first  <<", ";
  // } // Needs c++11

  if(false){
    TString txtname(outfolder+"/runs_"+basename+".txt");
    ofstream txtfile(txtname);
    int prevrun(0);
    for(map<int, map<int, set<Long64_t> > >::const_iterator it = runs.begin(); it != runs.end(); ++it) {
      run = it->first;
      if(run/1000 != prevrun){
	prevrun = run/1000;
	txtfile<<endl;
      }
      txtfile << run << "  ";
    }
    txtfile<<endl;
    txtfile.close();
    cout<<endl<<"Written run numbers in "<<txtname<<endl;
  }
  cout<<endl<<endl;

  return 0;
}
