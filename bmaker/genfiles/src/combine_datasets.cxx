// combine_datasets: Finds all unique events in a list of cfA files

#include <ctime>

#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <unistd.h>  // getopt

#include "TString.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

using namespace std; 

int main(int argc, char *argv[]){
  time_t startTime, curTime;
  time(&startTime);

  TString file_datasets("txt/datasamples/alldata.txt"), infolder("");
  int c(0);
  while((c=getopt(argc, argv, "f:i:"))!=-1){
    switch(c){
    case 'i':
      infolder=optarg;
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
	<<"./run/combine_datasets.exe -i <infolder> -f <file_datasets>"<<endl<<endl;
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

  map<int, set<int> > events;
  int event, run;

  for(unsigned idata(0); idata < datasets.size(); idata++){
    TChain chain("tree"), treeglobal("treeglobal");
    TString filename(infolder+"/*"+datasets[idata]+"*.root");
    int files = chain.Add(filename);
    if(files<1) {
      cout<<"No files found for "<<filename<<endl;
      continue;
    }
    treeglobal.Add(filename);
    TString outname("out/baby_"+basename+"_");
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
	cout<<"Doing entry "<<entry<<" of "<<entries<<endl;
      }
      
      if(events.find(run) == events.end()) events[run] = set<int>(); // New run
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

  cout<<endl<<endl;

  return 0;
}
