//----------------------------------------------------------------------------
// utilities - Various functions used accross the code
//----------------------------------------------------------------------------

#ifndef H_UTILITIES
#define H_UTILITIES

#include <cstddef>
#include <cstdio>
#include <iostream>
#include <cmath>

#include <string>
#include <vector>

#include "TString.h"
#include "TTree.h"
#include "TGraph.h"

typedef std::pair<int,double> int_double;
typedef std::pair<double,double> double_double;
const long double PI = acos(-1.L);
enum varType{kDouble, kvDouble, kFloat, kvFloat, kInt, kvInt, kBool, kvBool};

int change_branch_one(TString indir, TString name, TString outdir, std::vector<TString> var_type, std::vector<TString> var, 
		      std::vector<std::vector<TString> > var_val, int totentries);
int change_branch_one(TString indir, TString name, TString outdir, std::vector<TString> var_type, std::vector<TString> var, 
		       std::vector<TString> var_val);
bool eigen2x2(float matrix[2][2], float &eig1, float &eig2);
bool id_big2small(const int_double& left, const int_double& right);
bool dd_big2small(const double_double& left, const double_double& right);
bool dd_small2big(const double_double& left, const double_double& right);
long double DeltaPhi(long double phi1, long double phi2);
long double SignedDeltaPhi(long double phi1, long double phi2);
float dR(float eta1, float eta2, float phi1, float phi2);
TString roundNumber(double num, int decimals, double denom=1.);
TString addCommas(double num);
long double AddInQuadrature(long double x, long double y);
long double GetMass(long double e, long double px, long double py, long double pz);
long double GetMT(long double m1, long double pt1, long double phi1,
                  long double m2, long double pt2, long double phi2);
long double GetMT(long double pt1, long double phi1,
                  long double pt2, long double phi2);
bool Contains(const std::string& text, const std::string& pattern);

std::vector<std::string> Tokenize(const std::string& input,
                                  const std::string& tokens=" ");
void get_count_and_uncertainty(TTree& tree,
                               const std::string& cut,
                               double& count,
                               double& uncertainty);
void AddPoint(TGraph& graph, const double x, const double y);
TString hoursMinSec(long seconds);
std::vector<TString> dirlist(const TString &folder,
                             const TString &inname="dir",
                             const TString &tag="");

template<class T>
T noNaN(T val, T defval=1.){
  if(isnan(val)) {
    std::cout<<"Value is NaN. Returning "<<defval<<std::endl;
    return defval;
  } else return val;
}

template<class T>
short Sign(T val){
  return (T(0) < val) - (val < T(0));
}

std::string execute(const std::string &cmd);
std::string RemoveTrailingNewlines(std::string str);

std::vector<double> LinearSpacing(size_t npts, double low, double high);

#endif
