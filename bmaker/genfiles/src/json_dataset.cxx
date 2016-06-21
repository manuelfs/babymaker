// json_ntuple: Writes out a JSON file with the runs and lumis in ntuple

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

  TString infolder(""), outfolder("out/"), tag("");
  int c(0);
  while((c=getopt(argc, argv, "f:i:t:"))!=-1){
    switch(c){
    case 'i':
      infolder=optarg;
      break;
    case 't':
      tag=optarg;
      break;
    default:
      break;
    }
  }
  if(tag=="" || infolder==""){
    cout<<endl<<"Specify input folder and tag: "
	<<"./run/json_ntuple.exe -i <infolder> -o <outfolder=out> -t <tag>"<<endl<<endl;
    return 1;
  }
  gSystem->mkdir(outfolder, kTRUE);

  map<int, set<int> > lumis;
  int lumi;
  int run;

  TChain chain("tree");
  TString filename(infolder+"/*"+tag+"*.root");
  int files = chain.Add(filename);
  if(files<1) {
    cout<<endl<<"No files found for "<<filename<<". Exiting"<<endl<<endl;
    return 0;
  }

  TBranch *b_lumi(NULL), *b_run(NULL);
  chain.SetBranchAddress("lumiblock", &lumi, &b_lumi);
  chain.SetBranchAddress("run", &run, &b_run);

  long entries(chain.GetEntries()), tree_entry;

  time(&startTime);
  cout<<endl<<"Doing "<<files<<" files in "<<filename<<" with "<<entries<<" entries"<<endl;
  for(int entry(0); entry<entries; entry++){
    if(entry!=0 && entry%1000000==0) {
      time(&curTime);
      int seconds(difftime(curTime,startTime));
	
      cout<<"Doing entry "<<setw(10)<<addCommas(entry)<<" of "<<addCommas(entries)
	  <<"    Took "<<setw(6)<<seconds<<" seconds at "
	  <<setw(4)<<roundNumber(entry,1,seconds*1000.)<<" kHz"<<endl;
    }
      
    tree_entry = chain.LoadTree(entry);
    b_run->GetEntry(tree_entry);
    b_lumi->GetEntry(tree_entry);

    if(lumis.find(run) == lumis.end()) lumis[run] = set<int>(); // New run
    lumis[run].insert(lumi);
  } // Loop over entries

  TString txtname(outfolder+"/json_"+tag+".txt");
  ofstream txtfile(txtname);
  txtfile<<"{";
  for(map<int, set<int> >::const_iterator it = lumis.begin(); it != lumis.end(); ++it) {
    if(it != lumis.begin()) txtfile<<"],"<<endl;
    run = it->first;
    txtfile<<"\""<<run<<"\": [";
    int prev_lumi = -99;
    size_t ind = 0;
    for (set<int>::iterator itlumi = lumis[run].begin(); itlumi != lumis[run].end(); ++itlumi){
      lumi = *itlumi;
      if(itlumi == lumis[run].begin()) txtfile<<"["<<lumi<<", ";
      else {
	if(lumi > prev_lumi+1) {
	  txtfile<<prev_lumi<<"], ["<<lumi<<", ";
	}
      } // if not first lumi
      if(ind == lumis[run].size()-1) txtfile<<lumi<<"]";
      prev_lumi = lumi;
      ind++;
    } // Loop over lumis

  }
  txtfile<<"]}"<<endl;
  txtfile.close();
  time(&curTime);
  int seconds(difftime(curTime,startTime));
  cout<<endl<<"Written json in "<<txtname<<". Took "<<seconds<<" seconds"<<endl;
  cout<<endl<<endl;

  return 0;
}
