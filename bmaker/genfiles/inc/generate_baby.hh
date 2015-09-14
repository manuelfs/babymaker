#ifndef H_GENERATE_BABY
#define H_GENERATE_BABY

#include <vector>
#include <set>
#include <string>

class Variable{
public:
  Variable():
    type_(""),
    name_(""){
  }

  Variable(const std::string &type,
           const std::string &name):
    type_(type),
    name_(name){
  }

  bool operator<(const Variable& var) const{
    return type_<var.type_ || (type_==var.type_ && name_<var.name_);
  }

  std::string type_, name_;
};

bool Contains(const std::string &text, const std::string &pattern);

void WriteBaseHeader(const std::set<Variable> &all_vars,
                     const std::set<Variable> &com_vars,
                     const std::vector<std::string> &names);
void WriteBaseSource(const std::set<Variable> &all_vars,
                     const std::set<Variable> &com_vars,
                     const std::vector<std::string> &names);
void WriteSepHeader(const std::pair<std::string, std::set<Variable> > &vars);
void WriteSepSource(const std::pair<std::string, std::set<Variable> > &vars);

#endif
