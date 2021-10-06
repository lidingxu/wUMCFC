/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   PricingGraph.cpp
 * @brief  PricingGraph class for NTRT data
 * @author Liding XU
 */

 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <iostream>
#include <utility>
#include <list>
#include <vector>
#include <algorithm>
#include <deque>
#include <queue>

#include "objscip/objscip.h"
#include "Graph.h"
#include "ProblemGraph.h"
#include "PricingGraph.h"
#include "check.h"

using namespace std;

/** construct a pricing graph by problem graph */
PricingGraph::PricingGraph(
	Conflict * conflict_, /**< pointer to the problem graph */
	ProblemGraph * probgraph_ /**< pointer to the problem graph */
) {
	// construct and reserve the memomory
	conflict = conflict_;
	probgraph = probgraph_;
	vertex_map = vector<vector<int>>(probgraph_->nv);
	clique_map = vector<vector<int>>(conflict_->nc);
	nv = 2* probgraph_->ne;
	v_array.reserve(nv);
	v_to_in_edges = vector<vector<int>>(nv);
	v_to_out_edges = vector<vector<int>>(nv);
	int v_ind = 0;
	for (int i = 0; i < probgraph_->ne; i++) {
		// construct edge and its auxiliary
		addVertex(v_ind, false, i);
		v_ind++;
		addVertex(v_ind, true, i);
		v_ind++;
	}
	int e_ind = 0;
	for (int i = 0; i < probgraph_->ne; i++) {
		addEdge(e_ind, 2 * i, 2 * i + 1, true, i, i,  probgraph->e_array[i].e_cost);
		e_ind++;
		int v = probgraph_->e_array[i].e_head;
		for(int out: probgraph_->v_to_out_edges[v]){
			addEdge(e_ind, 2 * i + 1, 2 * out, false, i, out);
			e_ind++;
		}
	}
	ne = e_ind;
};

/** add an edge with index and cost */
void PricingGraph::addEdge(
	const int e_ind, /**< edge index */
	const int tail,  /**< vertex index of edge tail */
	const int head, /**< vertex index of edge head */
	const bool is_auxiliary, /**< auxiliary indicator: (i,j) -> (\bar{i}, \bar{j})*/
	const int tail_e_ind, /**< edge index of tail*/
	const int head_e_ind, /**< edge index of head*/
	const SCIP_Real e_bs_cost_ /**< edge cost */
) {
	e_array.push_back(Edge_Pricing(e_ind, tail, head, is_auxiliary, tail_e_ind, head_e_ind, e_bs_cost_));
	edge_map[make_pair(tail, head)] = e_ind;
	// update out edges
	v_to_out_edges[tail].push_back(e_ind);
	// update in edges
	v_to_in_edges[head].push_back(e_ind);
	if (is_auxiliary) {
		vertex_map[probgraph->e_array[head_e_ind].e_head].push_back(e_ind);
		for (int c_ind : *(conflict->getCINDS(head_e_ind))) {
			clique_map[c_ind].push_back(e_ind);
		}
	}
	edge_pair_map[make_pair(tail_e_ind, head_e_ind)] = e_ind;
}

/** (Experimental)  Bellman-Ford algorithm to find shortest path from s to t*/
SCIP_Real PricingGraph::BellmanFord(
	SCIP * scip, /**< scip pointer */
	int s, /**< source vertex index of problem graph*/
	int t,  /**< target vertex index of problem graph */
	vector <int>& path /**< shortest path edge array from s to t in problem graph*/
) {
	SCIP_Real shortest_len = SCIP_DEFAULT_INFINITY;
	vector<int> shortest_preds;
	shortest_preds.reserve(nv);
	int shortest_s_e_ind = -1;
	int shortest_t_e_ind = -1;
	vector<SCIP_Real> dists(nv);
	vector<int> preds(nv);
	for (int s_e_ind : *probgraph->getOutEdges(s)) {
		s_e_ind = 2 * s_e_ind;
		for (int t_e_ind : *probgraph->getInEdges(t)) {
			t_e_ind = t_e_ind * 2 + 1;
			fill(dists.begin(), dists.end(), SCIP_DEFAULT_INFINITY);
			fill(preds.begin(), preds.end(), -1);
			dists[s_e_ind] = 0;
			for (int i = 0; i < nv - 1; i++) {
				for (int e_ind = 0; e_ind < ne; e_ind++) {
					if (e_array[e_ind].is_forbidden) {
						continue;
					}
					int tail = e_array[e_ind].e_tail;
					int head = e_array[e_ind].e_head;
					SCIP_Real cost = getCost(e_ind);
					if (SCIPisLT(scip, dists[tail] + cost, dists[head])) {
						dists[head] = dists[tail] + cost;
						preds[head] = tail;
					}
				}
			}
			if( preds[t_e_ind] == -1){
				continue;
			}
			#ifdef SCIP_DEBUG
			SCIP_Bool negative_cycle = FALSE;
			for(auto e: e_array){
				if (e.is_forbidden) {
						continue;
				}
				int tail = e.e_tail;
				int head = e.e_head;
				if(dists[tail] == SCIP_DEFAULT_INFINITY){
					continue;
				}
				if(SCIPisLT(scip, dists[tail] + e.e_cost , dists[head])){
					negative_cycle = TRUE;
				}
			}
			if(negative_cycle){
				//SCIPdebugMessage("------!\n");
				SCIPerrorMessage("negative cycle deteced bias!\n");
				#ifdef CHECK_PAIR
				int count = 0;
				int all = 0;
				for(auto e: e_array){
        			for(auto e_: *getOutEdges(e.e_head)){
						SCIP_Real val = e.e_cost + e_array[e_].e_cost;
						if(val < 0){
							//SCIPdebugMessage("val:%lf e1:%s e2:%s\n",val, e.is_auxiliary ? "aux" : "non", e_array[e_].is_auxiliary ? "aux":"non");
           					//assert(e.e_cost + e_array[e_].e_cost >= 0);
							count++;
						}
						all++;
        			}
    			}
				SCIPdebugMessage("ratio neg:%lf\n", (count + 0.0)  / all);
				#endif
			}
			#endif
			if (dists[t_e_ind] < shortest_len) {
				shortest_len = dists[t_e_ind];
				shortest_preds = preds;
				shortest_s_e_ind = s_e_ind;
				shortest_t_e_ind = t_e_ind;
			}
		}
	}
	if(shortest_t_e_ind == -1){
		path.clear();
		return SCIP_DEFAULT_INFINITY;
	}
	vector<int> rev_pricing_path;
	vector<SCIP_Bool> visited(nv, FALSE); 
	rev_pricing_path.reserve(nv);
	int v = shortest_t_e_ind;
	while (v != shortest_s_e_ind) {
		// no path exists
		if(visited[v]){
			SCIPdebugMessage("revisited!!!!!!!!\n");
		}
		assert(!visited[v]);
		visited[v] = TRUE;
		rev_pricing_path.push_back(v);
		v = shortest_preds[v];
		assert(v >= 0);
		assert(v < nv);
	}
	rev_pricing_path.push_back(shortest_s_e_ind);
	path.clear();
	for (int i = rev_pricing_path.size() - 1; i >= 0; i -= 2) {
		int prob_e_ind = v_array[rev_pricing_path[i]].prob_e_ind;
		path.push_back(prob_e_ind);
	}
	return shortest_len;
}


/** (Experimental)  Shortest Paths Faster Algorithm (SPFA) to find shortest path from s to t*/
SCIP_Real PricingGraph::SPFA_ADJ(
	SCIP * scip, /**< scip pointer */
	int s, /**< source vertex index of problem graph*/
	int t,  /**< target vertex index of problem graph */
	vector <int>& path /**< shortest path edge array from s to t in problem graph*/
) {
	SCIP_Real shortest_len = SCIP_DEFAULT_INFINITY;
	vector<int> shortest_preds;
	shortest_preds.reserve(nv);
	int shortest_s_e_ind = -1;
	int shortest_t_e_ind = -1;
	vector<SCIP_Real> dists(nv);
	vector<int> preds(nv);
	vector<SCIP_Bool> inque(nv);
	// construct adjcacent list
	
	vector<list<int>> adjacent_list(nv, list<int>());
	for(int v = 0; v < nv; v++){
		for(int e_ind: *getOutEdges(v)){
			if (!e_array[e_ind].is_forbidden) {
				adjacent_list[v].push_back(e_ind);
			}
		}
	}

	for (int s_e_ind : *probgraph->getOutEdges(s)) {
		// loop through all outgoing edges from s
		s_e_ind = 2 * s_e_ind;
		// initialize distance array to infinity
		fill(dists.begin(), dists.end(), SCIP_DEFAULT_INFINITY);
		// initialize predecessor array to -1 (null)
		fill(preds.begin(), preds.end(), -1);
		// initialize queue indicator array to FALSE 
		fill(inque.begin(), inque.end(), FALSE);
		dists[s_e_ind] = 0;
		deque<int> v_deque;
		// push source to queue
		v_deque.push_back(s_e_ind);
		inque[s_e_ind] = TRUE;
		while(!v_deque.empty()){
			// pop queue
			int v = v_deque.front();
			v_deque.pop_front();
			inque[v] = FALSE;
			for(int e_ind: adjacent_list[v]){
				int head = e_array[e_ind].e_head;
				SCIP_Real tmp_dist = dists[v] + getCost(e_ind);
				// relaxation
				if (SCIPisLT(scip, tmp_dist, dists[head])) {
					dists[head] = tmp_dist;
					preds[head] = v;
					if(!inque[head]){
						v_deque.push_back(head);
						inque[head] = TRUE;
					}
				}
			}
		}
		for (int t_e_ind : *probgraph->getInEdges(t)) {
			t_e_ind = t_e_ind * 2 + 1;
			if( preds[t_e_ind] == -1){
				continue;
			}
			if (dists[t_e_ind] < shortest_len) {
				shortest_len = dists[t_e_ind];
				shortest_preds = preds;
				shortest_s_e_ind = s_e_ind;
				shortest_t_e_ind = t_e_ind;
			}
		}
	}
	if(shortest_t_e_ind == -1){
		path.clear();
		return SCIP_DEFAULT_INFINITY;
	}
	vector<int> rev_pricing_path;
	vector<SCIP_Bool> visited(nv, FALSE); 
	rev_pricing_path.reserve(nv);
	int v = shortest_t_e_ind;
	while (v != shortest_s_e_ind) {
		// no path exists
		//assert(!visited[v]);
		visited[v] = TRUE;
		rev_pricing_path.push_back(v);
		v = shortest_preds[v];
		//assert(v >= 0);
		//assert(v < nv);
	}
	rev_pricing_path.push_back(shortest_s_e_ind);
	path.clear();
	for (int i = rev_pricing_path.size() - 1; i >= 0; i -= 2) {
		int prob_e_ind = v_array[rev_pricing_path[i]].prob_e_ind;
		path.push_back(prob_e_ind);
	}
	return shortest_len;
}


/** (Experimental)  Shortest Paths Faster Algorithm (SPFA) to find shortest path from s to t*/
SCIP_Real PricingGraph::SPFA(
	SCIP * scip, /**< scip pointer */
	int s, /**< source vertex index of problem graph*/
	int t,  /**< target vertex index of problem graph */
	vector <int>& path /**< shortest path edge array from s to t in problem graph*/
) {
	SCIP_Real shortest_len = SCIP_DEFAULT_INFINITY;
	vector<int> shortest_preds;
	shortest_preds.reserve(nv);
	int shortest_s_e_ind = -1;
	int shortest_t_e_ind = -1;
	vector<SCIP_Real> dists(nv);
	vector<int> preds(nv);
	vector<SCIP_Bool> inque(nv);
	// construct adjcacent list
	

	for (int s_e_ind : *probgraph->getOutEdges(s)) {
		// loop through all outgoing edges from s
		s_e_ind = 2 * s_e_ind;
		// initialize distance array to infinity
		fill(dists.begin(), dists.end(), SCIP_DEFAULT_INFINITY);
		// initialize predecessor array to -1 (null)
		fill(preds.begin(), preds.end(), -1);
		// initialize queue indicator array to FALSE 
		fill(inque.begin(), inque.end(), FALSE);
		dists[s_e_ind] = 0;
		deque<int> v_deque;
		// push source to queue
		v_deque.push_back(s_e_ind);
		inque[s_e_ind] = TRUE;
		while(!v_deque.empty()){
			// pop queue
			int v = v_deque.front();
			v_deque.pop_front();
			inque[v] = FALSE;
			//for(int e_ind: adjacent_list[v]){

			for(int e_ind: *getOutEdges(v)){
				if (e_array[e_ind].is_forbidden) {
						continue;
				}
				int head = e_array[e_ind].e_head;
				SCIP_Real tmp_dist = dists[v] + getCost(e_ind);
				// relaxation
				if (SCIPisLT(scip, tmp_dist, dists[head])) {
					dists[head] = tmp_dist;
					preds[head] = v;
					if(!inque[head]){
						v_deque.push_back(head);
						inque[head] = TRUE;
					}
				}
			}
		}
		for (int t_e_ind : *probgraph->getInEdges(t)) {
			t_e_ind = t_e_ind * 2 + 1;
			if( preds[t_e_ind] == -1){
				continue;
			}
			if (dists[t_e_ind] < shortest_len) {
				shortest_len = dists[t_e_ind];
				shortest_preds = preds;
				shortest_s_e_ind = s_e_ind;
				shortest_t_e_ind = t_e_ind;
			}
		}
	}
	if(shortest_t_e_ind == -1){
		path.clear();
		return SCIP_DEFAULT_INFINITY;
	}
	vector<int> rev_pricing_path;
	vector<SCIP_Bool> visited(nv, FALSE); 
	rev_pricing_path.reserve(nv);
	int v = shortest_t_e_ind;
	while (v != shortest_s_e_ind) {
		// no path exists
		//assert(!visited[v]);
		visited[v] = TRUE;
		rev_pricing_path.push_back(v);
		v = shortest_preds[v];
		//assert(v >= 0);
		//assert(v < nv);
	}
	rev_pricing_path.push_back(shortest_s_e_ind);
	path.clear();
	for (int i = rev_pricing_path.size() - 1; i >= 0; i -= 2) {
		int prob_e_ind = v_array[rev_pricing_path[i]].prob_e_ind;
		path.push_back(prob_e_ind);
	}
	return shortest_len;
}


/** (Experimental) Shortest Paths Faster Algorithm (SPFA) to find shortest path from s to t*/
SCIP_Real PricingGraph::SPFA_PSEUDO(
	SCIP * scip, /**< scip pointer */
	int s, /**< source vertex index of problem graph*/
	int t,  /**< target vertex index of problem graph */
	vector <int>& path /**< shortest path edge array from s to t in problem graph*/
) {
	SCIP_Real shortest_len = SCIP_DEFAULT_INFINITY;
	int shortest_s_e_ind = -1;
	const int pseudo_s = -2;
	int shortest_t_e_ind = -1;
	vector<SCIP_Real> dists(nv, SCIP_DEFAULT_INFINITY);
	vector<int> preds(nv, -1);
	vector<SCIP_Bool> inque(nv, FALSE);

	deque<int> v_deque;
	for (int s_e_ind : *probgraph->getOutEdges(s)) {
		s_e_ind *=  2;
		dists[s_e_ind] = 0;
		preds[s_e_ind] = pseudo_s;
		// push source to queue
		v_deque.push_back(s_e_ind);
		inque[s_e_ind] = TRUE;
	}
	int count = 0;
	while(!v_deque.empty()){
		// pop queue
		int v = v_deque.front();
		v_deque.pop_front();
		inque[v] = FALSE;
		//for(int e_ind: adjacent_list[v]){
		for(int e_ind: *getOutEdges(v)){
			if (e_array[e_ind].is_forbidden) {
					continue;
			}
			int head = e_array[e_ind].e_head;
			SCIP_Real tmp_dist = dists[v] + getCost(e_ind) + 1e-8;
			// relaxation
			if (tmp_dist - dists[head] < -1e-9) {
				dists[head] = tmp_dist;
				preds[head] = v;
				if(!inque[head]){
					v_deque.push_back(head);
					inque[head] = TRUE;
				}
			}
		}
		#ifdef SCIP_DEBUG
			count++;
			if(count >4000){
				//SCIPdebugMessage("%d\n", count);
			}
		#endif
	}
	for (int t_e_ind : *probgraph->getInEdges(t)) {
		t_e_ind = t_e_ind * 2 + 1;
		if( preds[t_e_ind] == -1){
			continue;
		}
		if (dists[t_e_ind] < shortest_len) {
			shortest_len = dists[t_e_ind];
			shortest_t_e_ind = t_e_ind;
		}
	}
	if(shortest_t_e_ind == -1){
		path.clear();
		return SCIP_DEFAULT_INFINITY;
	}
	vector<int> rev_pricing_path;
	//vector<SCIP_Bool> visited(nv, FALSE); 
	rev_pricing_path.reserve(nv);
	int v = shortest_t_e_ind;
	while (v != pseudo_s) {
		// no path exists
		//assert(!visited[v]);
		//visited[v] = TRUE;
		rev_pricing_path.push_back(v);
		v = preds[v];
		//assert(v >= 0);
		//assert(v < nv);
	}
	path.clear();
	for (int i = rev_pricing_path.size() - 1; i >= 0; i -= 2) {
		int prob_e_ind = v_array[rev_pricing_path[i]].prob_e_ind;
		path.push_back(prob_e_ind);
	}
	//printf("%ld \n", path.size());
	return shortest_len;
}

class Cmp
{
  vector<SCIP_Real> * dists;
public:
  Cmp(vector<SCIP_Real> * dists_)
    {dists = dists_;}
  bool operator() (const int& a, const int& b) const
  {
    return (*dists)[a] > (*dists)[b];
  }
};


/** Two Way Dijkstra Algorithm to find shortest path from s to t*/
SCIP_Real PricingGraph::TWDijkstra(
	SCIP * scip, /**< scip pointer */
	int s, /**< source vertex index of problem graph*/
	int t,  /**< target vertex index of problem graph */
	vector <int>& path /**< shortest path edge array from s to t in problem graph*/
) {
	SCIP_Real shortest_len = SCIP_DEFAULT_INFINITY;
	const int pseudo_s = -2;


	vector<SCIP_Real> dists(nv, SCIP_DEFAULT_INFINITY);
	vector<int> preds(nv, -1);

	Cmp cmp(&dists);
	priority_queue<int,vector<int>,Cmp> pq (cmp);

	// push s_e_ind to pq
	for (int s_e_ind : *probgraph->getOutEdges(s)) {
		int p_v_ind = 2 * s_e_ind;
		int p_e_ind = getEdgeByEnds(p_v_ind, p_v_ind + 1)->e_ind;
		if (e_array[p_e_ind].is_forbidden) {
				continue;
		}
		#ifdef SCIP_DEBUG
		if(!e_array[p_e_ind].is_auxiliary){
			SCIPdebugMessage("not auxiliary\n");
		}
		#endif
		dists[p_v_ind + 1] = getCost(p_e_ind);
		preds[p_v_ind + 1] =  pseudo_s;
		// push source to queue
		pq.push(p_v_ind + 1);
	}

	int find_p_v = nv;
	while(!pq.empty()){
		int p_v_ind = pq.top();
		pq.pop();
		if(v_array[p_v_ind].is_auxiliary && probgraph->e_array[v_array[p_v_ind].prob_e_ind].e_head == t){
			shortest_len = dists[p_v_ind];
			find_p_v = p_v_ind;
			break;
		}
		for(int e_ind: *getOutEdges(p_v_ind)){
			if (e_array[e_ind].is_forbidden) {
					continue;
			}
			SCIP_Real tmp_dist = dists[p_v_ind] + getCost(e_ind) + 1e-6;
			int head = e_array[e_ind].e_head;
			#ifdef SCIP_DEBUG
				if(v_array[head].is_auxiliary){
					SCIPdebugMessage("v auxiliary\n");
				}
			#endif
			for(int e_ind_: *getOutEdges(head)){
				if (e_array[e_ind_].is_forbidden) {
					continue;
				}
				SCIP_Real tmp_dist_ = tmp_dist + getCost(e_ind_);
				int head_ = e_array[e_ind_].e_head;
				//(SCIPisNegative(scip, getCost(e_ind) + getCost(e_ind_))){
				//	SCIP_Real cost = getCost(e_ind) + getCost(e_ind_);
				//}
				//printf("%lf %lf %lf\n",  getCost(e_ind) + getCost(e_ind_), getCost(e_ind),  getCost(e_ind_));
				if( tmp_dist_ - dists[head_] < -1e-8){
					dists[head_] = tmp_dist_;
					preds[head_] = p_v_ind;
					pq.push(head_);
				}
			}
		}
	}

	if(find_p_v == nv){
		path.clear();
		return SCIP_DEFAULT_INFINITY;
	}
	vector<int> rev_pricing_path;
	rev_pricing_path.reserve(nv);

	int v = find_p_v;
	while (v != pseudo_s) {
		// no path exists
		//assert(!vis[v]);
		//vis[v] = TRUE;
		rev_pricing_path.push_back(v_array[v].prob_e_ind);
		v = preds[v];
		//assert(v != -1);
		//assert(v < nv);
	}
	path.clear();
	for (int i = rev_pricing_path.size() - 1; i >= 0; i--) {
		path.push_back(rev_pricing_path[i]);
	}
	//printf("%d %d %d %d %lf\n", probgraph->e_array[path[0]].e_tail,probgraph->e_array[path.back()].e_head, s, t, shortest_len);
	return shortest_len;
}