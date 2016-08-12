// change_weight: Produces new ntuples with the weights and systematics normalized to average 1

#include <iostream>
#include <ctime>
#include <string>
#include <vector>
#include <iomanip>   // setw

#include "TChain.h"
#include "TError.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "baby_full.hh"
#include "utilities.hh"
#include "cross_sections.hh"

using namespace std;

namespace{
  const int NSYSTS = 19;
}

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches       
  time_t begtime, endtime;
  time(&begtime);

  if(argc<2){
    cout<<"Format: ./run/change_weights.exe <infolder>=. <sample>=\"*.root\" <outfolder>=infolder"<<endl<<endl<<endl;
    return 1;
  }
  
  // Take command line arguments
  TString folder("."), sample("*.root"); 
  if(argc>=2) folder=argv[1]; 
  if(argc>=3) sample=argv[2]; 
  TString outfolder(folder);
  if(argc>=4) outfolder=argv[3]; 
  if(!folder.EndsWith("/")) folder.Append("/");
  if(!outfolder.EndsWith("/")) outfolder.Append("/");
  gSystem->mkdir(outfolder, kTRUE);

  float xsec(0.);
  if (!sample.Contains("SMS")){
    TChain chglob("treeglobal");
    chglob.Add((folder+sample).Data());
    size_t foundent = chglob.GetEntries();
    if (foundent==0) {
      cout<<"[Change Weights] ERROR: No entries found! Exit."<<endl;
      exit(1);
    } else {
      cout<<"[Change Weights] Found "<<foundent<<" entries."<<endl;
    }
    chglob.SetBranchAddress("xsec", &xsec);
    chglob.GetEntry(0); 
    cout<<"Found xsec = "<<xsec<<endl;
  }
  
  //Setup chains
  baby_full ch((folder+sample).Data());         
  ch.GetEntry(0); //Set branches to get size
 
  vector<TString> var_type, var; 
  vector<vector<TString> > var_val(NSYSTS);
  double nent_eff=0, nent=0, sum_weff_l0=0, sum_weff_l1=0, sum_btag=0, sum_pu=0, sum_toppt=0, sum_wisr=0;
  // vector<double> sum_wpdf(ch.w_pdf().size(),0);
  vector<double> sum_bctag(ch.sys_bctag().size()), sum_udsgtag(ch.sys_udsgtag().size());
  vector<double> sum_fs_bctag(ch.sys_fs_bctag().size()), sum_fs_udsgtag(ch.sys_fs_udsgtag().size());
  vector<double> sum_isr(ch.sys_isr().size());
  // vector<double> sum_spdf(2,0.);
  vector<double> sum_mur(ch.sys_mur().size()), sum_muf(ch.sys_muf().size()), sum_murf(ch.sys_murf().size());
  double nent_zlep=0, sum_wlep=0, sum_fs_wlep=0;
  vector<double> sum_slep(ch.sys_lep().size()), sum_fs_slep(ch.sys_fs_lep().size());
  // double sum_pdf_min(0.);
 
  int nentries = ch.GetEntries();
  //Loop over events and get sum of weights
  time(&begtime);
  for(int ientry=0; ientry<nentries; ientry++){
    ch.GetEntry(ientry); 

    if(ientry%500000==0 || ientry==nentries-1) {
      time(&endtime);
      int seconds(difftime(endtime,begtime));
      cout<<"Doing entry "<<setw(10)<<addCommas(ientry+1)<<" of "<<addCommas(nentries)
	  <<"    Took "<<setw(6)<<seconds<<" seconds at "
	  <<setw(4)<<roundNumber(ientry,1,seconds*1000.)<<" kHz"<<endl;
    }
    float lsign = ch.w_lumi()>0 ? 1:-1;
    nent_eff += lsign;
    nent++;

    // nent_eff ++;
    sum_btag  +=  noNaN(ch.w_btag());
    sum_pu    +=  noNaN(ch.w_pu());
    sum_toppt +=  noNaN(ch.w_toppt());
    sum_wisr  +=  noNaN(ch.w_isr());
    // for(unsigned int isys=0;isys<ch.w_pdf().size();isys++)           sum_wpdf[isys]        +=  noNaN(ch.w_pdf().at(isys));
    for(unsigned int isys=0;isys<ch.sys_bctag().size();isys++)       sum_bctag[isys]       +=  noNaN(ch.sys_bctag().at(isys));
    for(unsigned int isys=0;isys<ch.sys_udsgtag().size();isys++)     sum_udsgtag[isys]     +=  noNaN(ch.sys_udsgtag().at(isys));
    for(unsigned int isys=0;isys<ch.sys_fs_bctag().size();isys++)    sum_fs_bctag[isys]    +=  noNaN(ch.sys_fs_bctag().at(isys));
    for(unsigned int isys=0;isys<ch.sys_fs_udsgtag().size();isys++)  sum_fs_udsgtag[isys]  +=  noNaN(ch.sys_fs_udsgtag().at(isys));
    for(unsigned int isys=0;isys<ch.sys_isr().size();isys++)         sum_isr[isys]         +=  noNaN(ch.sys_isr().at(isys));
    for(unsigned int isys=0;isys<ch.sys_mur().size();isys++)         sum_mur[isys]         +=  noNaN(ch.sys_mur().at(isys));
    for(unsigned int isys=0;isys<ch.sys_muf().size();isys++)         sum_muf[isys]         +=  noNaN(ch.sys_muf().at(isys));
    for(unsigned int isys=0;isys<ch.sys_murf().size();isys++)        sum_murf[isys]        +=  noNaN(ch.sys_murf().at(isys));
    // if (ch.sys_pdf().size()==0) { 
    //   for(unsigned int isys=0;isys<2;isys++)  sum_spdf[isys]  +=  1; 
    // } else {
    //   for(unsigned int isys=0;isys<ch.sys_pdf().size();isys++)       sum_spdf[isys]        +=  noNaN(ch.sys_pdf().at(isys));
    // }

    // // Hack to recompute sys_pdf[1] which had a 1e-3 cut
    // float minpdf(1e10);
    // for(unsigned int isys=0;isys<ch.w_pdf().size();isys++)  
    //   if(ch.w_pdf()[isys] < minpdf) minpdf = ch.w_pdf()[isys];
    // if(ch.w_pdf().size()==0)
    //   sum_pdf_min += 1;
    // else
    //   sum_pdf_min += minpdf;
				      
    double weight = noNaN(ch.w_lep()) * noNaN(ch.w_fs_lep()) * noNaN(ch.w_btag()) *
      noNaN(ch.w_isr()) * noNaN(ch.w_pu());
    //Lepton weights
    if(ch.nleps()==0) {
      nent_zlep += 1;
      sum_weff_l0  += weight;
    } else{
      sum_weff_l1  += weight;
      sum_wlep     +=  ch.w_lep();
      sum_fs_wlep  +=  ch.w_fs_lep();
      for(unsigned int isys=0;isys<ch.sys_lep().size();isys++)        sum_slep[isys]        +=  ch.sys_lep().at(isys);
      for(unsigned int isys=0;isys<ch.sys_fs_lep().size();isys++)     sum_fs_slep[isys]     +=  ch.sys_fs_lep().at(isys);
    }
  }
  // // Hack to recompute sys_pdf[1] which had a 1e-3 cut
  // sum_spdf[1] = sum_pdf_min;

  const float luminosity = 1000.;
  if (sample.Contains("SMS")){
    float exsec(0.);
    if(sample.Contains("T1") || sample.Contains("T5")) xsec::signalCrossSection(ch.mgluino(), xsec, exsec);
    else  xsec::stopCrossSection(ch.mgluino(), xsec, exsec);
  }
  float w_lumi = xsec*luminosity / nent_eff;
  //float w_lumi_corr = w_lumi / fabs(ch.w_lumi());

  // Average w_toppt in bins of ht_isr_me found in inclusive TTJets (SingleLept and DiLept)
  double wanted_toppt(1.);
  if(sample.Contains("TTJets_HT-600to800"))   wanted_toppt = 0.8577;
  if(sample.Contains("TTJets_HT-800to1200"))  wanted_toppt = 0.8352;
  if(sample.Contains("TTJets_HT-1200to2500")) wanted_toppt = 0.8201;
  if(sample.Contains("TTJets_HT-2500toInf"))  wanted_toppt = 0.8198;

  // Average values found with find_w_isr.cxx
  double wanted_w_isr(1.);
  vector<double> wanted_sys_isr(2,1.);
  if(sample.Contains("TTJets_HT-600to800")) {
    wanted_w_isr = 0.7838;
    wanted_sys_isr[0] = 0.8965;
    wanted_sys_isr[1] = 0.6604;
  }
  if(sample.Contains("TTJets_HT-800to1200")) {
    wanted_w_isr = 0.7600;
    wanted_sys_isr[0] = 0.8851;
    wanted_sys_isr[1] = 0.6230;
  }
  if(sample.Contains("TTJets_HT-1200to2500")) {
    wanted_w_isr = 0.7365;
    wanted_sys_isr[0] = 0.8739;
    wanted_sys_isr[1] = 0.5861;
  }
  if(sample.Contains("TTJets_HT-2500toInf")) {
    wanted_w_isr = 0.7254;
    wanted_sys_isr[0] = 0.8686;
    wanted_sys_isr[1] = 0.5687;
  }

  time(&endtime); 
  int seconds = difftime(endtime, begtime);
  float hertz = nentries; hertz /= seconds;
  cout<<"[Change Weights] Completed in "<<seconds<<" seconds ("<<hoursMinSec(seconds)<<") for "<<nentries
      <<" events -> "<<roundNumber(hertz,1,1000)<<" kHz, "<<roundNumber(1000,2,hertz)<<" ms per event"<<endl<<endl;

  //Set type and var name
  var_type.push_back("float");      var.push_back("weight");
  var_type.push_back("float");      var.push_back("w_btag");            
  var_type.push_back("float");      var.push_back("w_pu");              
  var_type.push_back("float");      var.push_back("w_toppt");           
  var_type.push_back("float");      var.push_back("w_isr");            
  // var_type.push_back("vfloat");     var.push_back("w_pdf");             
  var_type.push_back("vfloat");     var.push_back("sys_bctag");         
  var_type.push_back("vfloat");     var.push_back("sys_udsgtag");       
  var_type.push_back("vfloat");     var.push_back("sys_fs_bctag");      
  var_type.push_back("vfloat");     var.push_back("sys_fs_udsgtag");    
  var_type.push_back("vfloat");     var.push_back("sys_isr");
  var_type.push_back("vfloat");     var.push_back("sys_pdf");
  var_type.push_back("vfloat");     var.push_back("sys_mur");
  var_type.push_back("vfloat");     var.push_back("sys_muf");
  var_type.push_back("vfloat");     var.push_back("sys_murf");
  var_type.push_back("float");      var.push_back("w_lep");
  var_type.push_back("float");      var.push_back("w_fs_lep");
  var_type.push_back("vfloat");     var.push_back("sys_lep");
  var_type.push_back("vfloat");     var.push_back("sys_fs_lep");
  var_type.push_back("float");      var.push_back("w_lumi");            

  // cout<<"sum_weff "<<sum_weff<<", wcorr "<<nent_eff*w_lumi_corr/sum_weff<<", nent "<<nent_eff<<
  //   ", w_lumi_corr "<<w_lumi_corr<<endl;
  //Calculate weights
  float w_corr_l0 = (nent-sum_wlep)/nent_zlep * (nent-sum_fs_wlep)/nent_zlep;
  if(nent_zlep==0) w_corr_l0 = 1.;
  var_val[0].push_back("*"+to_string(nent*wanted_w_isr/(sum_weff_l0*w_corr_l0 + sum_weff_l1)));
  var_val[1].push_back("*"+to_string(nent/sum_btag));
  var_val[2].push_back("*"+to_string(nent/sum_pu));
  var_val[3].push_back("*"+to_string(nent*wanted_toppt/sum_toppt));
  var_val[4].push_back("*"+to_string(nent*wanted_w_isr/sum_wisr));
  // for(unsigned int idx=0;idx<sum_wpdf.size();idx++)         var_val[4].push_back("*"+to_string(nent/sum_wpdf[idx]));
  for(unsigned int idx=0;idx<sum_bctag.size();idx++)        var_val[5].push_back("*"+to_string(nent/sum_bctag[idx]));
  for(unsigned int idx=0;idx<sum_udsgtag.size();idx++)      var_val[6].push_back("*"+to_string(nent/sum_udsgtag[idx]));
  for(unsigned int idx=0;idx<sum_fs_bctag.size();idx++)     var_val[7].push_back("*"+to_string(nent/sum_fs_bctag[idx]));
  for(unsigned int idx=0;idx<sum_fs_udsgtag.size();idx++)   var_val[8].push_back("*"+to_string(nent/sum_fs_udsgtag[idx]));
  for(unsigned int idx=0;idx<sum_isr.size();idx++)          var_val[9].push_back("*"+to_string(nent*wanted_sys_isr[idx]/sum_isr[idx]));
  // for(unsigned int idx=0;idx<sum_spdf.size();idx++)         var_val[10].push_back("*"+to_string(nent/sum_spdf[idx]));
  for(unsigned int idx=0;idx<sum_mur.size();idx++)          var_val[11].push_back("*"+to_string(nent/sum_mur[idx]));
  for(unsigned int idx=0;idx<sum_muf.size();idx++)          var_val[12].push_back("*"+to_string(nent/sum_muf[idx]));
  for(unsigned int idx=0;idx<sum_murf.size();idx++)         var_val[13].push_back("*"+to_string(nent/sum_murf[idx]));
  //Calculate lepton weights
  var_val[14].push_back("*"+to_string((nent-sum_wlep)/nent_zlep));
  var_val[15].push_back("*"+to_string((nent-sum_fs_wlep)/nent_zlep));
  for(unsigned int idx=0;idx<sum_slep.size();idx++)         
    var_val[16].push_back("*"+to_string((nent-sum_slep[idx])/nent_zlep));
  for(unsigned int idx=0;idx<sum_fs_slep.size();idx++)      
    var_val[17].push_back("*"+to_string((nent-sum_fs_slep[idx])/nent_zlep));
  var_val[18].push_back(roundNumber(w_lumi,12));


  int totentries(0);
  vector<TString> files = dirlist(folder,sample);
  for(unsigned int i=0; i<files.size(); i++){
    cout<<"[Change Weights] File "<<i+1<<"/"<<files.size()<<": "<<files[i]<<endl;
    totentries += change_branch_one(folder, files[i], outfolder, var_type, var, var_val, nentries);
  }

  time(&endtime); 
  seconds = difftime(endtime, begtime);
  hertz = totentries; hertz /= seconds;
  cout<<endl<<"Took "<<seconds<<" seconds ("<<hoursMinSec(seconds)<<") for "<<totentries
      <<" events -> "<<roundNumber(hertz,1,1000)<<" kHz, "<<roundNumber(1000,2,hertz)<<" ms per event"<<endl<<endl;
}
