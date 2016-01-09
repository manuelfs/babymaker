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
  int c(0);
  while((c=getopt(argc, argv, "f:i:o:"))!=-1){
    switch(c){
    case 'i':
      infolder=optarg;
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
	<<"./run/combine_datasets.exe -i <infolder> -o <outfolder=out> -f <file_datasets=txt/singlelep.txt>"<<endl<<endl;
    return 1;
  }

  vector<TString> datasets;
  TString buffer, basename("Run2015D");
  ifstream indata(file_datasets);
  while(indata){
    indata >> buffer;
    if(buffer!=""){
      datasets.push_back(buffer);
      basename += ("_"+buffer);
    }
  }

  map<int, set<Long64_t> > events;
  Long64_t event;
  int run;

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
    TString outname(outfolder+"/baby_"+basename+"_");
    outname += idata; outname += ".root";
    TFile outfile(outname, "RECREATE");
    outfile.cd();

    TTree *outtree(chain.CloneTree(0));

    chain.SetBranchAddress("event", &event);
    chain.SetBranchAddress("run", &run);

    long entries(chain.GetEntries());
    // entries = 100;

    cout<<endl<<"Doing "<<files<<" files in "<<filename<<" with "<<entries<<" entries"<<endl;
    for(int entry(0); entry<entries; entry++){
      chain.GetEntry(entry);
      if(entry!=0 && entry%250000==0) {
	time(&curTime);
	int seconds(difftime(curTime,startTime));
	
	cout<<"Doing entry "<<setw(10)<<addCommas(entry)<<" of "<<addCommas(entries)
	    <<"    Took "<<setw(6)<<seconds<<" seconds at "
	    <<setw(4)<<roundNumber(entry,1,seconds*1000.)<<" kHz"<<endl;
      }
      
      if(events.find(run) == events.end()) events[run] = set<Long64_t>(); // New run
      if(events[run].find(event) == events[run].end()){ // New event
	events[run].insert(event);
	outtree->Fill();
      } 
    } // Loop over entries
    outtree->Write();
    treeglobal.CloneTree(-1,"fast");
    outfile.Write();
    outfile.Close();
    time(&curTime);
    cout<<"Took "<<difftime(curTime,startTime) <<" seconds to write "<<outname<<endl;
    time(&startTime);

  } // Loop over datasets

  // for(auto it = events.cbegin(); it != events.cend(); ++it) {
  //   cout << it->first  <<", ";
  // } // Needs c++11

  TString txtname(outfolder+"/runs_"+basename+".txt");
  ofstream txtfile(txtname);
  int prevrun(0);
  for(map<int, set<Long64_t> >::const_iterator it = events.begin(); it != events.end(); ++it) {
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

  cout<<endl<<endl;

  return 0;
}
