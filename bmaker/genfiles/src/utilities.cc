//----------------------------------------------------------------------------
// utilities - Various functions used accross the code
//----------------------------------------------------------------------------

#ifndef INT_ROOT
#include "utilities.hh"
#endif

#include <cmath>

#include <iostream>
#include <string>
#include <stdexcept>

#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TList.h"
#include "TCollection.h"
#include "TH1D.h"
#include "TTree.h"
#include "TGraph.h"

using namespace std;

// Returns list of directorites or files in folder
vector<TString> dirlist(const TString &folder,
                        const TString &inname,
                        const TString &tag){
  TString pwd(gSystem->pwd());
  vector<TString> v_dirs;
  TSystemDirectory dir(folder, folder);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=static_cast<TSystemFile*>(next()))) {
      fname = file->GetName();
      if (inname=="dir") {
        if ((file->IsDirectory() && !fname.Contains(".") && fname.Contains(tag))) v_dirs.push_back(fname);
      } else  if(fname.Contains(inname)) v_dirs.push_back(fname);
    }
  } // if(files)
  gSystem->cd(pwd); // The TSystemDirectory object seems to change current folder
  return v_dirs;
}

bool eigen2x2(float matrix[2][2], float &eig1, float &eig2){
  float root = pow(matrix[0][0],2) + pow(matrix[1][1],2)-2*matrix[0][0]*matrix[1][1]+4*matrix[0][1]*matrix[1][0];
  if(root<0) return false;

  eig1 = (matrix[0][0]+matrix[1][1]+sqrt(root))/2.;
  eig2 = (matrix[0][0]+matrix[1][1]-sqrt(root))/2.;
  return true;
}

bool id_big2small(const int_double& left, const int_double& right){
  return left.second > right.second;
}

bool dd_small2big(const double_double& left, const double_double& right){
  return left.first < right.first;
}

bool dd_big2small(const double_double& left, const double_double& right){
  return left.first > right.first;
}

long double DeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(fabs(phi2-phi1), 2.L*PI);
  return dphi>PI ? 2.L*PI-dphi : dphi;
}

long double SignedDeltaPhi(long double phi1, long double phi2){
  long double dphi = fmod(phi2-phi1, 2.L*PI);
  if(dphi>PI){
    return dphi-2.L*PI;
  }else if(dphi<-PI){
    return dphi+2.L*PI;
  }else{
    return dphi;
  }
}

float dR(float eta1, float eta2, float phi1, float phi2) {
  return AddInQuadrature(eta1-eta2, DeltaPhi(phi1,phi2));
}

TString roundNumber(double num, int decimals, double denom){
  if(denom==0) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }

  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  return result;
}

TString addCommas(double num){
  TString result(""); result += num;
  int posdot(result.First('.'));
  if(posdot==-1) posdot = result.Length();
  for(int ind(posdot-3); ind > 0; ind -= 3)
    result.Insert(ind, ",");
  return result;
}

long double AddInQuadrature(long double x, long double y){
  if(fabs(y)>fabs(x)){
    const long double temp = y;
    y=x;
    x=temp;
  }
  if(x==0.) return y;
  const long double rat=y/x;
  return fabs(x)*sqrt(1.0L+rat*rat);
}

long double GetMass(long double e, long double px, long double py, long double pz){
  px/=e; py/=e; pz/=e;
  return fabs(e)*sqrt(1.0L-px*px-py*py-pz*pz);
}

long double GetMT(long double m1, long double pt1, long double phi1,
                  long double m2, long double pt2, long double phi2){
  return sqrt(m1*m1+m2*m2+2.L*(sqrt((m1*m1+pt1*pt1)*(m2*m2+pt2*pt2))-pt1*pt2*cos(phi2-phi1)));
}

long double GetMT(long double pt1, long double phi1,
                  long double pt2, long double phi2){
  //Faster calculation in massless case
  return sqrt(2.L*pt1*pt2*(1.L-cos(phi2-phi1)));
}

bool Contains(const string& text, const string& pattern){
  return text.find(pattern) != string::npos;
}

vector<string> Tokenize(const string& input,
                        const string& tokens){
  char* ipt(new char[input.size()+1]);
  memcpy(ipt, input.data(), input.size());
  ipt[input.size()]=static_cast<char>(0);
  char* ptr(strtok(ipt, tokens.c_str()));
  vector<string> output(0);
  while(ptr!=NULL){
    output.push_back(ptr);
    ptr=strtok(NULL, tokens.c_str());
  }
  return output;
}

void get_count_and_uncertainty(TTree& tree,
                               const string& cut,
                               double& count,
                               double& uncertainty){
  const string hist_name("temp");
  TH1D temp(hist_name.c_str(), "", 1, -1.0, 1.0);
  tree.Project(hist_name.c_str(), "0.0", cut.c_str());
  count=temp.IntegralAndError(0,2,uncertainty);
}

void AddPoint(TGraph& graph, const double x, const double y){
  graph.SetPoint(graph.GetN(), x, y);
}

string execute(const string &cmd){
  FILE *pipe = popen(cmd.c_str(), "r");
  if(!pipe) throw runtime_error("Could not open pipe.");
  const size_t buffer_size = 128;
  char buffer[buffer_size];
  string result = "";
  while(!feof(pipe)){
    if(fgets(buffer, buffer_size, pipe) != NULL) result += buffer;
  }

  pclose(pipe);
  return result;
}

string RemoveTrailingNewlines(string str){
  while(!str.empty() && str.at(str.length()-1) == '\n'){
    str.erase(str.length()-1);
  }
  return str;
}

vector<double> LinearSpacing(size_t npts, double low, double high){
  vector<double> pts(npts,low+0.5*(high-low));
  if(npts>1){
    double gap = (high-low)/(npts-1.0);
    for(size_t pt = 0; pt < npts; ++pt){
      pts.at(pt) = low+pt*gap;
    }
  }
  return pts;
}

