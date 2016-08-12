// find_w_isr: Finds the average w_isr in the TTJets HT bins for renormalization


#include "utilities.hh"

#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#include "TChain.h"
#include "TString.h"
#include "TH1D.h"

using namespace std;


int main(){
  time_t startTime;
  time(&startTime);

  TString hname="histo";
  TH1D histo(hname,"",1000,0,1.5); 
  TString folder="/net/cms26/cms26r0/babymaker/babies/2016_08_10/to_renormalize/manuelf/tt/";
  vector<TString> files({"TTJets_Tune", "TTJets_HT-600to800", "TTJets_HT-800to1200", "TTJets_HT-1200to2500",
	"TTJets_HT-2500toInf"});
  vector<TString> vars({"w_isr", "sys_isr[0]", "sys_isr[1]"});
  vector<double> average;

  vector<TChain*> chains;
  for(size_t ind=0; ind<files.size(); ind++){
    chains.push_back(new TChain("tree"));
    chains.back()->Add(folder+"/*"+files[ind]+"*.root");
    // cout<<endl<<"== "<<files[ind]<<": "<<addCommas(chains.back()->GetEntries())<<" entries in "
    // 	<<nfiles<<" files"<<endl;
    if(ind>0) cout<<"  if(sample.Contains(\""<<files[ind]<<"\")) {"<<endl;
    for(size_t ivar=0; ivar<vars.size(); ivar++){
      chains.back()->Project(hname, vars[ivar],"","goff");
      double mean = histo.GetMean();
      if(ind==0) average.push_back(mean);
      else {
	TString wanted = "wanted_"+vars[ivar];
	//wanted.ReplaceAll("[",""); wanted.ReplaceAll("]","");
	cout<<"    "<<wanted<<" = "<<roundNumber(mean,4,average[ivar])<<";"<<endl;
      }
    } // Loop over variables
    if(ind>0) cout<<"  }"<<endl;
  }



  time_t curTime;
  time(&curTime);
  int seconds = difftime(curTime,startTime);
  cout<<endl<<"Took "<< seconds << " seconds ("<<hoursMinSec(seconds)<<") to find average weights"<<endl<<endl;
}

