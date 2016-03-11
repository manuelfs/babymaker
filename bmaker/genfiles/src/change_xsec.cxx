// change_xsec: Changes cross section of ntuple

#include <iostream>
#include <ctime>
#include <string>
#include <vector>

#include "TChain.h"
#include "TError.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "utilities.hh"
#include "cross_sections.hh"

using namespace std;

int main(int argc, char *argv[]){
  time_t begtime, endtime;
  time(&begtime);

  if(argc<1){
    cout<<"Format: ./run/change_xsec.exe infile <outfolder>"<<endl;
    return 1;
  }
  
  // Take command line arguments
  TString infile("/net/cms2/cms2r0/babymaker/babies/2016_02_09/mc/T2tt/baby_SMS-T2tt_mGluino-550_mLSP-375_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9-v1_renorm.root"), outfolder("new_stops"); 
  if(argc>=2) infile = argv[1]; 
  if(argc>=3) outfolder=argv[2]; 
  if(!outfolder.EndsWith("/")) outfolder.Append("/");
  gSystem->mkdir(outfolder, kTRUE);

  TString folder(infile), file(infile);
  folder.Remove(folder.Last('/')+1, folder.Length());
  file.Remove(0, file.Last('/')+1);

  // Parsing masses
  TString mstop_s(infile);
  mstop_s.Remove(0,mstop_s.Index("ino-")+4); 
  mstop_s.Remove(mstop_s.Index('_'), mstop_s.Length());
  int mstop(mstop_s.Atoi());
  TString mlsp_s(infile);
  mlsp_s.Remove(0,mlsp_s.Index("LSP-")+4); 
  mlsp_s.Remove(mlsp_s.Index('_'), mlsp_s.Length());
  int mlsp(mlsp_s.Atoi());
  TString stop_lsp = "mGluino-"+mstop_s+"_mLSP-"+mlsp_s;


  float xsecOri, xsecNew, exsec;
  xsec::stopCrossSection(mstop, xsecOri, exsec);

  for(int mstopNew = mstop+25; mstopNew <= 1500; mstopNew += 25){
    xsec::stopCrossSection(mstopNew, xsecNew, exsec);
    TString factorXsec("*"); factorXsec += xsecNew/xsecOri;
    cout<<"mstop "<<mstop<<", xsecOri "<<xsecOri<<", new stop "<<mstopNew<<", xsecNew "<<xsecNew<<", factor "<<factorXsec<<endl;
    TString newname = file;
    TString stopNew_lsp = "mGluino-"+to_string(mstopNew)+"_mLSP-"+to_string(mlsp+mstopNew-mstop);
    newname.ReplaceAll(stop_lsp, stopNew_lsp);
    vector<TString> var_types({"float", "float"}), vars({"weight","w_lumi"}), var_vals({factorXsec, factorXsec});
    change_branch_one(folder, file, outfolder, var_types, vars, var_vals, newname);
    cout<<"Saved "<<outfolder<<newname<<endl;
  }

  time(&endtime); 
  int seconds = difftime(endtime, begtime);
  cout<<endl<<"Took "<<seconds<<" seconds ("<<hoursMinSec(seconds)<<") "<<endl<<endl;
}
