
#include <vector>
#include <set>
#include "ProblemGraph.h"
#include "PricingGraph.h"
using namespace std;

bool is_simple_path(const vector<int> & path);

void check_a_path(const vector<int> & path_e_array, ProblemGraph * probgraph,  const int s, const int t);

void check_qt(ProblemGraph * problemgraph, Quaternion & quat);

void negative_check(PricingGraph * pricinggraph);