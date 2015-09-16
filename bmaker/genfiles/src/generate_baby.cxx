#include "generate_baby.hh"

#include <cstring>

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <set>

#include <unistd.h>

using namespace std;

string ToCaps(string str){
  for(string::iterator it = str.begin();
      it != str.end();
      ++it){
    *it = toupper(*it);
  }
  return str;
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

vector<string> Tokenize(const string& input,
                        const string& tokens = " "){
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

string FixName(string name){
  //Variable names can have alphanumeric characters and underscores only
  string allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_";

  //Remove illegal characters
  size_t pos = name.find_first_not_of(allowed);
  while(pos < name.size()){
    name.at(pos) = '_';
    pos = name.find_first_not_of(allowed);
  }

  //Replace double underscore with single underscore
  pos = name.rfind("__");
  while(pos < name.size()){
    name.replace(pos, 2, "_");
    pos = name.rfind("__");
  }

  //Remove leading and trailing spaces
  pos = 0;
  for(pos = 0; pos < name.size(); ++pos){
    if(name.at(pos) != ' ') break;
  }
  size_t endpos = name.size();
  for(endpos = name.size(); endpos != 0; --endpos){
    if(name.at(endpos-1) != ' ') break;
  }

  return name.substr(pos, endpos-pos);
}

set<Variable> GetVariables(const string &file_name){
  string allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890_";
  set<Variable> vars;

  ifstream infile(("../../variables/"+file_name).c_str());
  string line;
  while(getline(infile, line)){
    size_t start = line.find_first_not_of(" ");
    if(start >= line.size() || line.at(start) == '#' || line.at(start) == '/') continue;

    //Replace double space with single space
    size_t pos = line.rfind("  ");
    while(pos < line.size()){
      line.replace(pos, 2, " ");
      pos = line.rfind("  ");
    }
    size_t end = line.find_last_of(allowed)+1;
    size_t split = line.rfind(' ', end)+1;

    vars.insert(Variable(line.substr(start, split-start),
                         line.substr(split, end-split)));
  }
  infile.close();

  return vars;
}

int main(){
  vector<string> file_names = Tokenize(execute("ls ../../variables/ 2> /dev/null"), "\n");

  vector<pair<string, set<Variable> > > sep_vars(file_names.size());
  vector<string> fixed_names(file_names.size());
  for(size_t ifile = 0; ifile < file_names.size(); ++ifile){
    fixed_names.at(ifile) = FixName(file_names.at(ifile));
    sep_vars.at(ifile).first = fixed_names.at(ifile);
    sep_vars.at(ifile).second = GetVariables(file_names.at(ifile));
  }

  set<Variable> all_vars;
  for(size_t ifile = 0; ifile < sep_vars.size(); ++ifile){
    for(set<Variable>::const_iterator var = sep_vars.at(ifile).second.begin();
        var != sep_vars.at(ifile).second.end();
        ++var){
      all_vars.insert(*var);
    }
  }

  set<Variable> com_vars;
  if(sep_vars.size()){
    for(set<Variable>::const_iterator var = sep_vars.at(0).second.begin();
        var != sep_vars.at(0).second.end();
        ++var){
      bool found_in_all = true;
      for(size_t ifile = 1; found_in_all && ifile < sep_vars.size(); ++ifile){
        if(sep_vars.at(ifile).second.find(*var) == sep_vars.at(ifile).second.end()){
          found_in_all = false;
        }
      }
      if(found_in_all){
        com_vars.insert(*var);
      }
    }
    for(set<Variable>::const_iterator var = com_vars.begin();
        var != com_vars.end();
        ++var){
      for(size_t ifile = 0; ifile < sep_vars.size(); ++ifile){
        sep_vars.at(ifile).second.erase(*var);
      }
    }
  }

  WriteBaseHeader(all_vars, com_vars, fixed_names);
  WriteBaseSource(all_vars, com_vars, fixed_names);

  for(size_t ifile = 0; ifile < sep_vars.size(); ++ifile){
    WriteSepHeader(sep_vars.at(ifile));
    WriteSepSource(sep_vars.at(ifile));
  }

}

bool Contains(const string &text, const string &pattern){
  return text.find(pattern) != string::npos;
}

void WriteBaseHeader(const set<Variable> &all_vars,
                     const set<Variable> &com_vars,
                     const vector<string> &names){

  ofstream file("../interface/baby_base.hh");

  file << "// baby_base: base class to handle reduced tree ntuples\n";
  file << "// File generated with generate_baby.exe\n\n";

  file << "#ifndef H_BABY_BASE\n";
  file << "#define H_BABY_BASE\n\n";

  file << "#include <vector>\n";
  file << "#include <string>\n\n";

  file << "#include \"TTree.h\"\n";
  file << "#include \"TChain.h\"\n";
  //  file << "#include \"TTreeFormula.h\"\n\n";

  file << "class baby_base{\n";
  file << "public:\n";
  file << "  baby_base(); // Constructor to create tree\n";
  file << "  baby_base(const std::string &filename); // Constructor to read tree\n\n";

  file << "  int Add(const std::string &filename);\n";

  file << "  long GetEntries() const;\n";
  file << "  virtual void GetEntry(const long entry);\n";
  file << "  bool PassString(TString cut);\n\n";

  file << "  virtual void Fill();\n";
  file << "  void Write();\n\n";

  file << "  virtual std::string Type() const;\n\n";

  file << "  static const double bad_val_;\n\n";

  file << "  virtual ~baby_base();\n\n";

  for(set<Variable>::const_iterator var = com_vars.begin();
      var != com_vars.end();
      ++var){
    file << "  " << var->type_ << " const & " << var->name_ << "() const;\n";
    file << "  " << var->type_ << " & " << var->name_ << "();\n";
  }
  file << '\n';

  for(set<Variable>::const_iterator var = all_vars.begin();
      var != all_vars.end();
      ++var){
    if(com_vars.find(*var) != com_vars.end()) continue;
    file << "  __attribute__((noreturn)) virtual "
         << var->type_ << " const & " << var->name_ << "() const;\n";
    file << "  __attribute__((noreturn)) virtual "
         << var->type_ << " & " << var->name_ << "();\n";
  }
  file << "  TChain chain_;\n";
  file << "  TTree tree_;\n";
  file << "  long entry_;\n";
  file << '\n';

  file << "protected:\n";
  file << "  const bool read_only_;\n\n";

  file << "private:\n";
  file << "  class VectorLoader{\n";
  file << "  public:\n";
  file << "    VectorLoader();\n";
  file << "  private:\n";
  file << "    static bool loaded_;\n";
  file << "  };\n\n";

  file << "  static VectorLoader vl_;\n";
  for(set<Variable>::const_iterator var = com_vars.begin();
      var != com_vars.end();
      ++var){
    file << "  " << var->type_ << ' ' << var->name_ << "_;\n";
    if(Contains(var->type_, "vector")){
      file << "  " << var->type_ << " *p_" << var->name_ << "_;\n";
    }
    file << "  TBranch *b_" << var->name_ << "_;\n";
    file << "  mutable bool c_" << var->name_ << "_;\n";
  }
  file << "};\n\n";

  file << "baby_base* NewTree(const std::type_info &type);\n\n";

  for(size_t i = 0; i < names.size(); ++i){
    file << "#include \"babymaker/bmaker/interface/baby_" << names.at(i) << ".hh\"\n";
  }
  file <<'\n';

  file << "#endif" << endl;

  file.close();
}

void WriteBaseSource(const set<Variable> &all_vars,
                     const set<Variable> &com_vars,
                     const vector<string> &names){
  ofstream file("../src/baby_base.cc");

  file << "// baby_base: base class to handle reduce tree ntuples\n";
  file << "//File generated with generate_baby.exe\n\n";

  file << "#include \"babymaker/bmaker/interface/baby_base.hh\"\n\n";

  file << "#include <stdexcept>\n";
  file << "#include <string>\n";
  file << "#include <iostream>\n";
  file << "#include <vector>\n\n";

  file << "#include \"TROOT.h\"\n";
  file << "#include \"TTree.h\"\n";
  file << "#include \"TBranch.h\"\n";
  file << "#include \"TChain.h\"\n";
  //  file << "#include \"TTreeFormula.h\"\n\n";

  file << "using namespace std;\n\n";

  file << "bool baby_base::VectorLoader::loaded_ = false;\n\n";

  file << "baby_base::VectorLoader baby_base::vl_ = baby_base::VectorLoader();\n\n";

  file << "baby_base::VectorLoader::VectorLoader(){\n";
  file << "  if(!loaded_){\n";
  file << "    gROOT->ProcessLine(\"#include <vector>\");\n";
  file << "    loaded_ = true;\n";
  file << "  }\n";
  file << "}\n\n";

  file << "const double baby_base::bad_val_ = -999.;\n\n";

  file << "baby_base::baby_base():\n";
  file << "  chain_(\"junk\", \"junk\"),\n";
  file << "  tree_(\"tree\", \"tree\"),\n";
  file << "  entry_(0),\n";
  if(com_vars.size()){
    const set<Variable>::const_iterator com_end_2 = --com_vars.end();
    file << "  read_only_(false),\n";
    for(set<Variable>::const_iterator var = com_vars.begin();
        var != com_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &p_" << var->name_ << "_)),\n";
      }else{
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &" << var->name_ << "_)),\n";
      }
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << com_end_2->name_ << "_(0),\n";
    file << "  b_" << com_end_2->name_ << "_(tree_.Branch(\"" << com_end_2->name_ << "\", &" << com_end_2->name_ << "_)),\n";
    file << "  c_" << com_end_2->name_ << "_(false){\n";
  }else{
    file << "  read_only_(false){\n";
  }
  file << "}\n\n";

  file << "baby_base::baby_base(const string &filename):\n";
  file << "  chain_(\"tree\",\"tree\"),\n";
  file << "  tree_(\"junk\",\"junk\"),\n";
  file << "  entry_(0),\n";
  if(com_vars.size()){
    const set<Variable>::const_iterator com_end_2 = --com_vars.end();
    file << "  read_only_(true),\n";
    for(set<Variable>::const_iterator var = com_vars.begin();
        var != com_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
      }
      file << "  b_" << var->name_ << "_(NULL),\n";
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << com_end_2->name_ << "_(0),\n";
    if(Contains(com_end_2->type_, "vector")){
      file << "  p_" << com_end_2->name_ << "_(&" << com_end_2->name_ << "_),\n";
    }
    file << "  b_" << com_end_2->name_ << "_(NULL),\n";
    file << "  c_" << com_end_2->name_ << "_(false){\n";
  }else{
    file << "  read_only_(true){\n";
  }
  file << "  chain_.Add(filename.c_str());\n";
  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &p_" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }else{
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }
  }
  file << "}\n\n";

  file << "void baby_base::Fill(){\n";
  file << "  if(read_only_){\n";
  file << "    throw std::logic_error(\"Trying to write to read-only tree\");\n";
  file << "  }else{\n";
  file << "    tree_.Fill();\n";
  file << "  }\n\n";

  file << "  //Resetting variables\n";
  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  " << var->name_ << "_.clear();\n";
    }else if(Contains(var->type_, "tring")){
      file << "  " << var->name_ << "_ = \"\";\n";
    }else{
      file << "  " << var->name_ << "_ = static_cast<" << var->type_ << ">(bad_val_);\n";
    }
  }
  file << "}\n\n";

  file << "void baby_base::Write(){\n";
  file << "  if(read_only_){\n";
  file << "    throw std::logic_error(\"Trying to write to read-only tree.\");\n";
  file << "  }else{\n";
  file << "    tree_.Write();\n";
  file << "  }\n";
  file << "}\n\n";

  file << "string baby_base::Type() const{\n";
  file << "  return \"\";\n";
  file << "}\n\n";

  file << "baby_base::~baby_base(){\n";
  file << "}\n\n";

  file << "int baby_base::Add(const std::string &filename){\n";
  file << "  if(!read_only_){\n";
  file << "    throw std::logic_error(\"Trying to add files to tree opened for writing.\");\n";
  file << "  }\n";
  file << "  return chain_.Add(filename.c_str());\n";
  file << "}\n\n";


  file << "bool baby_base::PassString(TString cut){\n";
  //  file << " TTreeFormula f(\"formula\",cut, &chain_);\n";
  //  file << " bool result = f.EvalInstance(0);\n";
  file << " bool result = true;\n";
  file << " return result;\n";
  file << "}\n\n";
  
  file << "long baby_base::GetEntries() const{\n";
  file << "  if(read_only_){\n";
  file << "    return chain_.GetEntries();\n";
  file << "  }else{\n";
  file << "    return tree_.GetEntries();\n";
  file << "  }\n";
  file << "}\n\n";

  file << "void baby_base::GetEntry(const long entry){\n";
  file << "  if(!read_only_){\n";
  file << "    throw std::logic_error(\"Trying to read from write-only tree.\");\n";
  file << "  }\n\n";

  for(set<Variable>::const_iterator var = com_vars.begin(); var!= com_vars.end(); ++var){
    file << "  c_" << var->name_ << "_ = false;\n";
  }
  file << "  entry_ = chain_.LoadTree(entry);\n";
  file << "}\n\n";

  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    file << var->type_ << " const & baby_base::" << var->name_ << "() const{\n";
    file << "  if(!read_only_){\n";
    file << "    throw std::logic_error(\"Trying to write to const tree.\");\n";
    file << "  }\n";
    file << "  if(!c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = com_vars.begin(); var != com_vars.end(); ++var){
    file << var->type_ << " & baby_base::" << var->name_ << "(){\n";
    file << "  if(read_only_ && !c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = all_vars.begin(); var != all_vars.end(); ++var){
    if(com_vars.find(*var) != com_vars.end()) continue;
    file << var->type_ << " const & baby_base::" << var->name_ << "() const{\n";
    file << "  throw std::logic_error(\"" << var->name_
         << " does not exist in this baby_base version.\");\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = all_vars.begin(); var != all_vars.end(); ++var){
    if(com_vars.find(*var) != com_vars.end()) continue;
    file << var->type_ << " & baby_base::" << var->name_ << "(){\n";
    file << "  throw std::logic_error(\"" << var->name_
         << " does not exist in this baby_base version.\");\n";
    file << "}\n\n";
  }

  for(size_t i = 0; i < names.size(); ++i){
    file << "#include \"babymaker/bmaker/interface/baby_" << names.at(i) << ".hh\"\n";
  }
  file << "baby_base* NewTree(const std::type_info &type){\n\n";
  file << "  if(type == typeid(baby_base)) return new baby_base;\n";
  for(size_t i = 0; i < names.size(); ++i){
    file << "  else if(type == typeid(baby_" << names.at(i) << ")) return static_cast<baby_base*>(new baby_" << names.at(i) << ");\n";
  }
  file << "  else return new baby_base;\n";
  file << "}\n\n";

  file.close();
}

void WriteSepHeader(const pair<string, set<Variable> > &sep_vars){
  string name = sep_vars.first;
  string NAME = ToCaps(name);
  set<Variable> vars = sep_vars.second;
  ofstream file(("../interface/baby_"+name+".hh").c_str());

  file << "// baby_" << name << ": " << name << " version of baby_base to handle reduce tree ntuples\n";
  file << "// File generated with generate_baby.exe\n\n";

  file << "#ifndef H_BABY_" << NAME << "\n";
  file << "#define H_BABY_" << NAME << "\n\n";

  file << "#include <vector>\n";
  file << "#include <string>\n\n";

  file << "#include \"TTree.h\"\n";
  file << "#include \"TChain.h\"\n\n";

  file << "#include \"babymaker/bmaker/interface/baby_base.hh\"\n\n";

  file << "class baby_" << name << " : public baby_base{\n";
  file << "public:\n";
  file << "  baby_" << name << "(); // Constructor to create tree\n";
  file << "  baby_" << name << "(const std::string &filename); // Constructor to read tree\n\n";

  file << "  virtual void GetEntry(const long entry);\n\n";

  file << "  virtual void Fill();\n\n";

  file << "  virtual std::string Type() const;\n\n";

  file << "  virtual ~baby_" << name << "();\n\n";

  for(set<Variable>::const_iterator var = vars.begin();
      var != vars.end();
      ++var){
    file << "  virtual " << var->type_ << " const & " << var->name_ << "() const;\n";
    file << "  virtual " << var->type_ << " & " << var->name_ << "();\n";
  }
  file << '\n';

  file << "private:\n";
  for(set<Variable>::const_iterator var = vars.begin();
      var != vars.end();
      ++var){
    file << "  " << var->type_ << ' ' << var->name_ << "_;\n";
    if(Contains(var->type_, "vector")){
      file << "  " << var->type_ << " *p_" << var->name_ << "_;\n";
    }
    file << "  TBranch *b_" << var->name_ << "_;\n";
    file << "  mutable bool c_" << var->name_ << "_;\n";
  }
  file << "};\n\n";

  file << "#endif" << endl;

  file.close();
}

void WriteSepSource(const pair<string, set<Variable> > &sep_vars){
  string name = sep_vars.first;
  string NAME = ToCaps(name);
  set<Variable> vars = sep_vars.second;
  ofstream file(("../src/baby_"+name+".cc").c_str());

  file << "// baby_" << name << ": " << name << " version of baby_base to handle reduce tree ntuples\n";
  file << "//File generated with generate_baby.exe\n\n";

  file << "#include \"babymaker/bmaker/interface/baby_base.hh\"\n\n";
  file << "#include \"babymaker/bmaker/interface/baby_" << name << ".hh\"\n\n";

  file << "#include <stdexcept>\n";
  file << "#include <string>\n";
  file << "#include <vector>\n\n";

  file << "#include \"TTree.h\"\n";
  file << "#include \"TBranch.h\"\n";
  file << "#include \"TChain.h\"\n\n";

  file << "using namespace std;\n\n";

  file << "baby_" << name << "::baby_" << name << "():\n";
  if(vars.size()){
    const set<Variable>::const_iterator vars_end_2 = --vars.end();
    file << "  baby_base(),\n";
    for(set<Variable>::const_iterator var = vars.begin();
        var != vars_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &p_" << var->name_ << "_)),\n";
      }else{
        file << "  b_" << var->name_ << "_(tree_.Branch(\"" << var->name_ << "\", &" << var->name_ << "_)),\n";
      }
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << vars_end_2->name_ << "_(0),\n";
    if(Contains(vars_end_2->type_, "vector")){
      file << "  p_" << vars_end_2->name_ << "_(&" << vars_end_2->name_ << "_),\n";
      file << "  b_" << vars_end_2->name_ << "_(tree_.Branch(\"" << vars_end_2->name_ << "\", &p_" << vars_end_2->name_ << "_)),\n";
    }else{
      file << "  b_" << vars_end_2->name_ << "_(tree_.Branch(\"" << vars_end_2->name_ << "\", &" << vars_end_2->name_ << "_)),\n";
    }
    file << "  c_" << vars_end_2->name_ << "_(false){\n";
  }else{
    file << "  baby_base(){\n";
  }
  file << "}\n\n";

  file << "baby_" << name << "::baby_" << name << "(const string &filename):\n";
  if(vars.size()){
    const set<Variable>::const_iterator vars_end_2 = --vars.end();
    file << "  baby_base(filename),\n";
    for(set<Variable>::const_iterator var = vars.begin();
        var != vars_end_2;
        ++var){
      file << "  " << var->name_ << "_(0),\n";
      if(Contains(var->type_, "vector")){
        file << "  p_" << var->name_ << "_(&" << var->name_ << "_),\n";
      }
      file << "  b_" << var->name_ << "_(NULL),\n";
      file << "  c_" << var->name_ << "_(false),\n";
    }
    file << "  " << vars_end_2->name_ << "_(0),\n";
    if(Contains(vars_end_2->type_, "vector")){
      file << "  p_" << vars_end_2->name_ << "_(&" << vars_end_2->name_ << "_),\n";
    }
    file << "  b_" << vars_end_2->name_ << "_(NULL),\n";
    file << "  c_" << vars_end_2->name_ << "_(false){\n";
  }else{
    file << "  baby_base(filename){\n";
  }
  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &p_" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }else{
      file << "  chain_.SetBranchAddress(\"" << var->name_ << "\", &" << var->name_ << "_, &b_" << var->name_ << "_);\n";
    }
  }
  file << "}\n\n";

  file << "void baby_" << name << "::Fill(){\n";
  file << "  baby_base::Fill();\n";

  file << "  //Resetting variables\n";
  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    if(Contains(var->type_, "vector")){
      file << "  " << var->name_ << "_.clear();\n";
    }else if(Contains(var->type_, "tring")){
      file << "  " << var->name_ << "_ = \"\";\n";
    }else{
      file << "  " << var->name_ << "_ = static_cast<" << var->type_ << ">(bad_val_);\n";
    }
  }
  file << "}\n\n";

  file << "string baby_" << name << "::Type() const{\n";
  file << "  return \"" << name << "\";\n";
  file << "}\n\n";

  file << "baby_" << name << "::~baby_" << name << "(){\n";
  file << "}\n\n";

  file << "void baby_" << name << "::GetEntry(const long entry){\n";
  file << "  baby_base::GetEntry(entry);\n\n";

  for(set<Variable>::const_iterator var = vars.begin(); var!= vars.end(); ++var){
    file << "  c_" << var->name_ << "_ = false;\n";
  }
  file << "}\n\n";

  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    file << var->type_ << " const & baby_" << name << "::" << var->name_ << "() const{\n";
    file << "  if(!read_only_){\n";
    file << "    throw std::logic_error(\"Trying to write to const tree.\");\n";
    file << "  }\n";
    file << "  if(!c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  for(set<Variable>::const_iterator var = vars.begin(); var != vars.end(); ++var){
    file << var->type_ << " & baby_" << name << "::" << var->name_ << "(){\n";
    file << "  if(read_only_ && !c_" << var->name_ << "_ && b_" << var->name_ <<"_){\n";
    file << "    b_" << var->name_ << "_->GetEntry(entry_);\n";
    file << "    c_" << var->name_ << "_ = true;\n";
    file << "  }\n";
    file << "  return " << var->name_ << "_;\n";
    file << "}\n\n";
  }

  file.close();
}

