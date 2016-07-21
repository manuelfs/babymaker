#include <locale>
#include <iostream>
#include <iomanip>
#include <string>
#include <set>

#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "utilities.hh"

using namespace std;

struct Variable{
  TString name;
  long zip_size, tot_size;

  bool operator<(const Variable &var) const{
    // return (zip_size<var.zip_size)
    //   || (zip_size==var.zip_size && tot_size<var.tot_size)
    //   || (zip_size==var.zip_size && tot_size==var.tot_size && name<var.name);
    return (name>var.name);
  }
};

class comma_numpunct : public numpunct<char>{
protected:
  virtual char do_thousands_sep() const{
    return ',';
  }

  virtual string do_grouping() const{
    return "\03";
  }
};

int findVar(const Variable &var, vector<Variable> &vars){
  int index = -1;
  for(size_t ivar=0; ivar<vars.size(); ivar++){
    if(var.name == vars[ivar].name){
      index = ivar;
      break;
    }
  }

  return index;
}

bool varCompare(const Variable &var1, const Variable &var2) {
  return var1.zip_size > var2.zip_size;
}

int main(int argc, char *argv[]){
  cout.imbue(locale(locale(), new comma_numpunct));

  if(argc<2) {
    cout<<endl<<"Format is: ./run/disk_usage.exe file <treename=tree>"<<endl<<endl;
    return 1;
  }

  TString nametree="tree", filename = argv[1];
  if(argc>=3) nametree = argv[2];

  TChain chain(nametree);
  if(!chain.Add(filename) || !chain.GetListOfLeaves()) {
    cout<<endl<<"No tree found in "<<filename<<endl<<endl;
    return 1;
  }

  long zip_sum(0), tot_sum(0);
  Ssiz_t max_length(0);
    
  vector<Variable> vars, allvars;
  set<Variable>::iterator it;

  for(int i = 0; i < chain.GetListOfLeaves()->GetSize(); ++i){
    TBranch *b = static_cast<TLeaf*>(chain.GetListOfLeaves()->At(i))->GetBranch();
    if(!b) continue;
    Variable v;
    v.name = b->GetName();
    v.zip_size = b->GetZipBytes();
    v.tot_size = b->GetTotBytes();

    allvars.push_back(v);
    zip_sum += v.zip_size;
    tot_sum += v.tot_size;
    if(v.name.Length() > max_length) max_length = v.name.Length();

    if(v.name.First('_')>=0) v.name.Remove(v.name.First('_'), v.name.Length()-v.name.First('_'));
    int index = findVar(v, vars);
    if(index == -1) vars.push_back(v);
    else {
      vars[index].zip_size += v.zip_size;
      vars[index].tot_size += v.tot_size;
    }
  } // Loop over branches

  int wbytes=16, wother=11;
  TString sep = "      ";
  sort(vars.begin(), vars.end(), varCompare);
  sort(allvars.begin(), allvars.end(), varCompare);
  cout << "Finding sizes of branches in "<<endl<<filename << endl;
  cout <<endl<< setw(max_length) << "Branch name" << ' '
       << setw(wbytes) << "Bytes" << ' '
       << setw(wother) << "Frac. [%]" << ' '
       << setw(wother) << "Cumulative"  << sep
       << setw(max_length) << "Branch group name" << ' '
       << setw(wbytes) << "Bytes" << ' '
       << setw(wother) << "Frac. [%]" << ' '
       << setw(wother) << "Cumulative" << endl;
  for(int ind=0; ind<(max_length+wbytes+2*wother+3); ind++) cout << "=";
  cout << sep;
  for(int ind=0; ind<(max_length+wbytes+2*wother+3); ind++) cout << "=";
  cout << endl 
       << setw(max_length) << "Total" << ' '
       << setw(wbytes) << zip_sum << ' '
       << setw(wother) << "100.00" << ' '
       << setw(wother) << "-" << sep
       << setw(max_length) << "Total" << ' '
       << setw(wbytes) << zip_sum << ' '
       << setw(wother) << "100.00" << ' '
       << setw(wother) << "-" << endl;
  long running_total(0), tot2=0;
  for(size_t ivar=0; ivar<allvars.size(); ivar++){
    running_total += allvars[ivar].zip_size;
    double this_frac = (100.0*allvars[ivar].zip_size)/zip_sum;
    double tot_frac = (100.0*running_total)/zip_sum;
    cout << setw(max_length) << allvars[ivar].name << " "
	 << setw(wbytes) << allvars[ivar].zip_size << " "
	 << setw(wother) << roundNumber(this_frac,2) << " "
	 << setw(wother) << roundNumber(tot_frac,2) << sep;
    if(ivar<vars.size()){
      tot2 += vars[ivar].zip_size;
      this_frac = (100.0*vars[ivar].zip_size)/zip_sum;
      tot_frac = (100.0*tot2)/zip_sum;
      cout << setw(max_length) << vars[ivar].name << " "
	   << setw(wbytes) << vars[ivar].zip_size << " "
	   << setw(wother) << roundNumber(this_frac,2) << " "
	   << setw(wother) << roundNumber(tot_frac,2);

    }
    cout<<endl;
  } // Loop over all branches
  cout<<endl;
}
