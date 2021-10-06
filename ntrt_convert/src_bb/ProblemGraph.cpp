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

/**@file   ProblemGraph.cpp
 * @brief  ProblemGraph class for NTRT data
 * @author Liding XU
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <list>
#include <algorithm>

#include "ProblemGraph.h"


/** Vertex class */
Vertex_Prob::Vertex_Prob(
	const int v_ind_, /**< vertex index */
	const SCIP_Real v_cost_ /**< vertex cost */
) : v_ind(v_ind_), v_cost(v_cost_) {
	v_cons_coeff = 0;
};


/** Edge class */
Edge_Prob::Edge_Prob(
	const int e_ind_, /**< edge index */
	const int e_tail_,  /**< vertex index of edge tail */
	const int e_head_, /**< vertex index of edge head */
	const SCIP_Real e_cost_, /**< edge cost */
	const SCIP_Real e_capacity_ /**< edge capacity */
) : e_ind(e_ind_), e_tail(e_tail_), e_head(e_head_), e_cost(e_cost_), e_capacity(e_capacity_) {};


/** Opportunity Quaternion*/
Quaternion::Quaternion(
	const int v1_ind_,
	const int v2_ind_,
	const int v3_ind_,
	const int e12_ind_,
	const int e23_ind_,
	const int e32_ind_,
	const int e21_ind_,
	const SCIP_Real rate_, /**< coding rate */
	const SCIP_Real redcost_ /**< reduced cost */
) :v1_ind(v1_ind_), v2_ind(v2_ind_), v3_ind(v3_ind_), e12_ind(e12_ind_),
e23_ind(e23_ind_), e32_ind(e32_ind_), e21_ind(e21_ind_), rate(rate_), redcost(redcost_) {};

/** ProblemGraph class */

/** construct a graph with nv_ vertcies and ne_ edges where vertices and edges are added larter on */
ProblemGraph::ProblemGraph(
	const int nv_, /**< the number of vertices */
	const int ne_  /**< the number of edges */
) : Graph(nv_, ne_) {};

/** copy constructor */
ProblemGraph::ProblemGraph(
	const ProblemGraph & probgraph_ /**< the graph to copy */
) : Graph(probgraph_) {
	quat_array = probgraph_.quat_array;
	nq = probgraph_.nq;
};

/** add a vertex with index and cost*/
void ProblemGraph::addVertex(
	const int v_ind_, /**< vertex index*/
	const SCIP_Real v_cost_ /**< vertex cost*/
) {
	v_array.push_back(Vertex_Prob(v_ind_, v_cost_));
};

/**< return the user-defined cost on edge e*/
inline SCIP_Real ProblemGraph::getCost(
	const int e_ind_  /**< vertex index*/
) {
	return flow_val / e_array[e_ind_].e_capacity;
};

/** add an edge with index and cost */
void ProblemGraph::addEdge(
	const int e_ind, /**< edge index*/
	const int tail,  /**< vertex index of edge tail */
	const int head, /**< vertex index of edge head */
	const SCIP_Real e_cost, /**< edge cost */
	const SCIP_Real e_capacity/**< edge capacity */
) {
	e_array.push_back(Edge_Prob(e_ind, tail, head, e_cost, e_capacity));
	// update edge map
	edge_map[make_pair(tail, head)] = e_ind;
	// update out edges
	v_to_out_edges[tail].push_back(e_ind);
	// update in edges
	v_to_in_edges[head].push_back(e_ind);

	// add to vertex coefficient
 	v_array[head].v_cons_coeff += 1 / (2 * e_capacity);

	// (v1,v2) (v2,v3), (v3,v2), (v2,v1)->, v2 = tail, v1 = head
	int v2_ind = tail;
	int v1_ind = head;
	int e21_ind = e_ind;
	for (auto e23_ind : v_to_out_edges[v2_ind]) {
		// get v3
		Edge_Prob* e23 = getEdgeByIND(e23_ind);
		int v3_ind = e23->e_head;
		if (v3_ind == v1_ind) {
			continue;
		}
		// check whether (v3, v2) and (v1, v2)
		Edge_Prob* e32 = getEdgeByEnds(v3_ind, v2_ind);
		Edge_Prob* e12 = getEdgeByEnds(v1_ind, v2_ind);
		// check whether quaternion exists
		if (e32 != NULL && e12 != NULL) {
			int e32_ind = e32->e_ind;
			int e12_ind = e12->e_ind;
			SCIP_Real rate = MIN(1 / e_capacity, 1 / e23->e_capacity);
			SCIP_Real red_cost = MIN(e_cost, e23->e_cost);
			quat_array.push_back(Quaternion(v1_ind, v2_ind, v3_ind, e12_ind,
				e23_ind, e32_ind, e21_ind, rate, red_cost));
		}
	}

	// (v1,v2)-> (v2,v3), (v3,v2), (v2,v1), v1 = tail, v2 = head
	v1_ind = tail;
	v2_ind = head;
	int e12_ind = e_ind;
	for (auto e23_ind : v_to_out_edges[v2_ind]) {
		// get v3
		Edge_Prob* e23 = getEdgeByIND(e23_ind);
		int v3_ind = e23->e_head;
		if (v3_ind == v1_ind) {
			continue;
		}
		// check whether (v3, v2) and (v1, v2)
		Edge_Prob* e32 = getEdgeByEnds(v3_ind, v2_ind);
		Edge_Prob* e21 = getEdgeByEnds(v2_ind, v1_ind);
		// check whether quaternion exists
		if (e32 != NULL && e21 != NULL) {
			int e32_ind = e32->e_ind;
			int e21_ind = e21->e_ind;
			SCIP_Real rate = MIN(1 / e21->e_capacity, 1 / e23->e_capacity);
			SCIP_Real red_cost = MIN(e21->e_cost, e23->e_cost);
			quat_array.push_back(Quaternion(v1_ind, v2_ind, v3_ind, e12_ind,
				e23_ind, e32_ind, e21_ind, rate, red_cost));
		}
	}
	nq = quat_array.size();
};

/** ConflictClique class */
ConflictClique::ConflictClique(
	const int c_ind_, /**< clique index */
	const int ne_  /**< the number of edges */
) : c_ind(c_ind_), ne(ne_) {
	/* reserve the memory */
	e_ind_array.reserve(ne_);
};

/** add an edge by index */
void ConflictClique::addEdge(
	const int e_ind_ /**< edge index */
) {
	e_ind_array.push_back(e_ind_);
};



/** Conflict class */

/**< construct a coonflict management class */
Conflict::Conflict(
	const ProblemGraph * probgraph, /**< the pointer to the problem graph */
	const int nc_  /**< the number ofconflict cliques */
) : nc(nc_) {
	c_array.reserve(nc_);
	edge_to_c_inds = vector<vector<int>>(probgraph->ne);
};

/** copy constructor */
Conflict::Conflict(
	const Conflict & conflict_ /**< the conflict to copy*/
) {
	probgraph = conflict_.probgraph;
	nc = conflict_.nc;
	c_array = conflict_.c_array;
	edge_to_c_inds = conflict_.edge_to_c_inds;
};

/** add an clique by index and its size*/
void Conflict::addClique(
	const int c_ind_, /**< clique index */
	const int ne_,  /**< the number of edges */
	const vector<int> & e_array /**< edges in the clique */
) {
	c_array.push_back(ConflictClique(c_ind_, ne_));
	int c_ind = c_array.size() - 1;
	for (int e_ind : e_array) {
		c_array[c_ind].addEdge(e_ind);
		edge_to_c_inds[e_ind].push_back(c_ind);
	}
};

/** return the clique indices of an edge */
vector<int> *  Conflict::getCINDS(
	const int e_ind_ /**< edge index */
) {
	return &edge_to_c_inds[e_ind_];
};
