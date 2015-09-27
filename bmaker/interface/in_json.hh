#ifndef H_IN_JSON
#define H_IN_JSON

#include <vector>
#include <string>

std::vector<std::vector<int> > MakeVRunLumi(std::string input);
bool inJSON(std::vector<std::vector<int> > VRunLumi, int Run, int LS);
void CheckVRunLumi(std::vector<std::vector<int> > VRunLumi);
void CheckVRunLumi2(std::vector<std::vector<int> > VVRunLumi);

#endif //H_IN_JSON
