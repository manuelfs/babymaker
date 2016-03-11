// merge_glu_stop: Merges gluino and stop scans

#include <iostream>
#include <ctime>
#include <string>
#include <vector>

#include "TChain.h"
#include "TError.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "utilities.hh"

using namespace std;

int main(int argc, char *argv[]){
  time_t begtime, endtime;
  time(&begtime);

  if(argc<1){
    cout<<"Format: ./run/merge_glu_stop.exe glufile <outfolder> <stopfolder>"<<endl;
    return 1;
  }
  
  // Take command line arguments
  TString glufile, stopfolder("/net/cms2/cms2r0/babymaker/babies/2016_02_09/mc/T2tt/");
  TString outfolder("mergedGluStop/"); 
  if(argc>=2) glufile = argv[1]; 
  if(argc>=3) outfolder=argv[2]; 
  if(argc>=4) stopfolder=argv[3]; 
  if(!outfolder.EndsWith("/")) outfolder.Append("/");
  gSystem->mkdir(outfolder, kTRUE);

  TString outfile(glufile);
  outfile.Remove(0, outfile.Last('/')+1);
  outfile.ReplaceAll("T5tttt","T5tttt-Stop");
  outfile = outfolder+"/"+outfile;
  
  // Parsing masses
  TString mlsp_s(glufile);
  mlsp_s.Remove(0,mlsp_s.Index("LSP-")+4); 
  mlsp_s.Remove(mlsp_s.Index('_'), mlsp_s.Length());
  int mlsp(mlsp_s.Atoi());

  // Finding desired stop mass
  int mstop = mlsp + 175;
  TString stop_lsp("*mGluino-"+to_string(mstop)+"_mLSP-"+mlsp_s+"*");
  vector<TString> stopfiles = dirlist(stopfolder, stop_lsp);
  if(stopfiles.size()==0){
    stopfolder = "new_stops/";
    stopfiles = dirlist(stopfolder, stop_lsp);
  }
  if(stopfiles.size()==0 || stopfiles.size()>=2){
    cout<<"Found "<<stopfiles.size()<<" files matching "<<stop_lsp<<". Exiting"<<endl<<endl;
    exit(1);
  }

  vector<TString> ntuples = {glufile, stopfolder+stopfiles[0]};
  mergeNtuples(ntuples, outfile);
  cout<<endl<<"Merged T5tttt and T2tt ntuples at "<<outfile<<endl;

  time(&endtime); 
  int seconds = difftime(endtime, begtime);
  cout<<endl<<"Took "<<seconds<<" seconds ("<<hoursMinSec(seconds)<<") "<<endl<<endl;
}
