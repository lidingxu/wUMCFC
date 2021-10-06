#define SCIP_DEBUG
#include <vector>
#include <set>
#include "ProblemGraph.h"
#include "PricingGraph.h"
using namespace std;

bool is_simple_path(const vector<int> & path){
    set<int> my_set;
    for(int i: path){
        auto it = my_set.find(i);
        if(it != my_set.end()){
            my_set.insert(i);
        }
        else{
            if(!my_set.empty()){
                return false;
            }    
        }
    }
    return true;
} 

void check_a_path(const vector<int> & path_e_array, ProblemGraph * problemgraph, const int s, const int t){
    vector<int> path_v_array;
    for(int i = 0; i < path_e_array.size() - 1; i++){
        assert(problemgraph->e_array[path_e_array[i]].e_head == problemgraph->e_array[path_e_array[i + 1]].e_tail);
    }
	// get the edge array and compute the path cost
	for (int j = 0; j < path_e_array.size(); j++) {
			Edge_Prob & e = problemgraph->e_array[path_e_array[j]];
			path_v_array.push_back(e.e_tail);
	}
	path_v_array.push_back(problemgraph->e_array[path_e_array.back()].e_head);
    assert(path_v_array[0] == s);
    assert(path_v_array.back() == t);
    assert(is_simple_path(path_v_array));
}

void check_qt(ProblemGraph * problemgraph, Quaternion & quat){
    assert(problemgraph->e_array[quat.e12_ind].e_tail == quat.v1_ind);
    assert(problemgraph->e_array[quat.e12_ind].e_head == quat.v2_ind);
    assert(problemgraph->e_array[quat.e23_ind].e_tail == quat.v2_ind);
    assert(problemgraph->e_array[quat.e23_ind].e_head == quat.v3_ind);
    assert(problemgraph->e_array[quat.e32_ind].e_tail == quat.v3_ind);
    assert(problemgraph->e_array[quat.e32_ind].e_head == quat.v2_ind);
    assert(problemgraph->e_array[quat.e21_ind].e_tail == quat.v2_ind);
    assert(problemgraph->e_array[quat.e21_ind].e_head == quat.v1_ind);
}

void negative_check(PricingGraph * pricinggraph){
    for(auto e: pricinggraph->e_array){
        for(auto e_: *pricinggraph->getOutEdges(e.e_head)){
            assert(e.e_cost + pricinggraph->e_array[e_].e_cost >= 0);
        }
    }
}