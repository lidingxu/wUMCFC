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

/**@file   PricingGraph.h
 * @brief  PricingGraph class for NTRT data
 * @author Liding XU
 */

 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRICINGGRAPH_H__
#define __PRICINGGRAPH_H__

#include <iostream>
#include <utility>
#include <map>
#include <vector>
#include <algorithm>

#include "objscip/objscip.h"
#include "Graph.h"
#include "ProblemGraph.h"

using namespace std;


/** Vertex Pricing class */
class Vertex_Pricing
{
public:
	/** constructs the vertex */
	Vertex_Pricing(
		const int v_ind_, /**< vertex index */
		const bool is_auxiliary_, /**< auxiliary indicator \bar{e}*/
		const int prob_e_ind_ /**< corresponding edge index in the problem */
	) : v_ind(v_ind_), is_auxiliary(is_auxiliary_), prob_e_ind(prob_e_ind_) {};

	int v_ind; /**< vertex index */
	bool is_auxiliary; /**< auxiliary indicator \bar{e}*/
	int prob_e_ind; /**< corresponding edge index in the problem */
};


/** Edge Pricing class */
class Edge_Pricing
{
public:
	/** construct the edge with index only */
	Edge_Pricing(
		const int ind_, /**< edge index */
		const int tail_,  /**< vertex index of edge tail */
		const int head_, /**< vertex index of edge head */
		const bool is_auxiliary_, /**< auxiliary indicator: (i,j) -> (\bar{i}, \bar{j})*/
		const int tail_e_ind_, /**< edge index of tail*/
		const int head_e_ind_, /**< edge index of head*/
		const SCIP_Real e_bs_cost_ = 0 /**< edge cost */
	) : e_ind(ind_), e_tail(tail_), e_head(head_), e_cost(0), is_auxiliary(is_auxiliary_),
		tail_e_ind(tail_e_ind_), head_e_ind(head_e_ind_), e_bs_cost(e_bs_cost_) {};

	int e_ind; /**< edge index */
	int e_tail; /**< vertex index of edge tail */
	int e_head; /**< vertex index of edge head */
	SCIP_Real e_cost; /**< edge cost */
	const SCIP_Real e_bs_cost; /**< base cost*/
	SCIP_Bool is_auxiliary; /**< auxiliary indicator: (i,j) -> (\bar{i}, \bar{j})*/
	SCIP_Bool is_forbidden; /**< is edge forbidden by branch decisions */
	int tail_e_ind; /**< edge index of tail*/
	int head_e_ind; /**< edge index of tail*/
};


/** PricingGraph class */
class PricingGraph : public Graph<Vertex_Pricing, Edge_Pricing>
{
public:
	/** construct a pricing graph by problem graph */
	PricingGraph(
		Conflict * conflict_, /**< pointer to the problem graph */
		ProblemGraph * probgraph_ /**< pointer to the problem graph */
		);

	/** add a vertex with index and cost*/
	inline void addVertex(
		const int v_ind_, /**< vertex index*/
		const bool is_auxiliary_, /**< auxiliary indicator \bar{e}*/
		const int prob_e_ind_ /**< corresponding edge index in the problem */
	) {
		v_array.push_back(Vertex_Pricing(v_ind_, is_auxiliary_, prob_e_ind_));
	};

	/** add an edge with index and cost */
	void addEdge(
		const int e_ind, /**< edge index */
		const int tail,  /**< vertex index of edge tail */
		const int head, /**< vertex index of edge head */
		const bool is_auxiliary, /**< auxiliary indicator: (i,j) -> (\bar{i}, \bar{j})*/
		const int tail_e_ind, /**< edge index of tail*/
		const int head_e_ind, /**< edge index of tail*/
		const SCIP_Real e_bs_cost_ = 0 /**< edge cost */
	);

	/** reset the edge to non forbbiden */
	void resetNonForbidden() {
		for (int i = 0; i < ne; i++) {
			e_array[i].is_forbidden = FALSE;
		}
	};

	/** reset the edge to non forbbiden */
	void setForbidden(
		vector<int>& forb_e_array /**< forbidden edge of problem graph array */
	) {
		for (int i = 0; i < forb_e_array.size(); i++) {
			int e_ind = getEIndbyEPair(forb_e_array[i], forb_e_array[i]);
			assert(e_ind != -1);
			if (e_ind != -1) {
				e_array[e_ind].is_forbidden = TRUE;
			}
		}
	}

	/** Bellman-Ford algorithm to find shortest path from s to t*/
	virtual SCIP_Real BellmanFord(
		SCIP * scip, /**< scip pointer */
		int s, /**< source vertex index of problem graph*/
		int t,  /**< target vertex index of problem graph */
		vector <int>& path /**< shortest path array from s to t in pricing graph*/
	);

	/** Shortest Paths Faster Algorithm (SPFA) to find shortest path from s to t*/
	SCIP_Real SPFA(
		SCIP * scip, /**< scip pointer */
		int s, /**< source vertex index of problem graph*/
		int t,  /**< target vertex index of problem graph */
		vector <int>& path /**< shortest path array from s to t in pricing graph*/
	);

	/** Shortest Paths Faster Algorithm (SPFA) to find shortest path from s to t*/
	SCIP_Real SPFA_ADJ(
		SCIP * scip, /**< scip pointer */
		int s, /**< source vertex index of problem graph*/
		int t,  /**< target vertex index of problem graph */
		vector <int>& path /**< shortest path array from s to t in pricing graph*/
	);


	/** Shortest Paths Faster Algorithm (SPFA) to find shortest path from s to t*/
	SCIP_Real SPFA_PSEUDO(
		SCIP * scip, /**< scip pointer */
		int s, /**< source vertex index of problem graph*/
		int t,  /**< target vertex index of problem graph */
		vector <int>& path /**< shortest path edge array from s to t in problem graph*/
	);

	/** Two Way Dijkstra Algorithm to find shortest path from s to t*/
	SCIP_Real TWDijkstra(
		SCIP * scip, /**< scip pointer */
		int s, /**< source vertex index of problem graph*/
		int t,  /**< target vertex index of problem graph */
		vector <int>& path /**< shortest path edge array from s to t in problem graph*/
	); 

	/** reset the edge cost to the basic cost for reduced cost case */
	void resetReducedCost() {
		for (int i = 0; i < ne; i++) {
			e_array[i].e_cost = e_array[i].e_bs_cost;
		}
	}

	/** reset the edge cost to the basic cost for reduced cost case */
	void resetFarkas() {
		for (int i = 0; i < ne; i++) {
			e_array[i].e_cost = 0;
		}
	}

	/**< return the user-defined cost on edge e*/
	inline SCIP_Real getCost(
		const int e_ind  /**< vertex index*/
	) {
		return e_array[e_ind].e_cost;
	};

	/**< return the pricing edge index by the problem edge pair, -1 does not exist */
	int getEIndbyEPair(
		const int e_start, /**< starting edge*/
		const int e_end /**< end edge */
	) {
		auto it = edge_pair_map.find(make_pair(e_start, e_end));
		if (it != edge_pair_map.end()) {
			return it->second;
		}
		else {
			return -1;
		}
	}
	
	/**< return edge index array by the problem vertex */
	vector<int> * getEIndsofV(
		const int v_ind /**< vertex index */
	) {
		return &vertex_map[v_ind];
	}

	/**< return edge index array by the problem clique */
	vector<int> * getEIndsofClq(
		const int v_ind /**< vertex index */
	) {
		return &clique_map[v_ind];
	}

	ProblemGraph * probgraph;
	Conflict* conflict;
	map<pair<int, int>, int> edge_pair_map; /** map the edge pair in problem graph to the edge in the price edge*/
	vector<vector<int>> vertex_map; /** map the vertex in problemg graph to pricing edges of which head edge's head is the vertex*/
	vector<vector<int>> clique_map; /** map the clique to pricing edges of which its ends is in the clique*/
};

#endif