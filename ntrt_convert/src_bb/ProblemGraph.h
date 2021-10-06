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

/**@file   ProblemGraph.h
 * @brief  ProblemGraph class for NTRT data
 * @author Liding XU
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PROBLEMGRAPH_H__
#define __PROBLEMGRAPH_H__

#include <utility>
#include <vector>
#include <algorithm>

#include "Graph.h"

using namespace std;

/** Vertex Problem class */
class Vertex_Prob
{
public:
	/** constructs the vertex */
	Vertex_Prob(
		const int v_ind_, /**< vertex index */
		const SCIP_Real v_cost_ /**< vertex cost */
		);

	int v_ind; /**< vertex index */
	SCIP_Real v_cost; /**< vertex cost */
	SCIP_Real v_cons_coeff; /**< vertex coeffs*/
};


/** Edge Problem class */
class Edge_Prob
{
public:
	/** construct the edge with index only */
	Edge_Prob(
		const int e_ind_, /**< edge index */
		const int e_tail_,  /**< vertex index of edge tail */
		const int e_head_, /**< vertex index of edge head */
		const SCIP_Real e_cost_, /**< edge cost */
		const SCIP_Real e_capacity_ /**< edge capacity */
		);

	int e_ind; /**< edge index */
	int e_tail; /**< vertex index of edge tail */
	int e_head; /**< vertex index of edge head */
	SCIP_Real e_cost; /**< edge cost */
	SCIP_Real e_capacity; /**< edge capacity */

};

/** Opportunity Quaternion*/
class Quaternion {
public:
	Quaternion(
		const int v1_ind_,
		const int v2_ind_,
		const int v3_ind_,
		const int e12_ind_,
		const int e23_ind_,
		const int e32_ind_,
		const int e21_ind_,
		const SCIP_Real rate_, /**< coding rate */
		const SCIP_Real redcost_ /**< reduced cost */
	);
	int v1_ind, v2_ind, v3_ind; /**< opportunity quaternion vertices: {v1, v2 , v3}*/
	int e12_ind, e23_ind, e32_ind, e21_ind; /**< opportunity quaternion edges: (v1,v2) ,(v2,v3), (v3,v2), (v2,v1)*/
	SCIP_Real rate; /**< coding rate */
	SCIP_Real redcost; /**< reduced cost */
};

/** ProblemGraph class */
class ProblemGraph : public Graph<Vertex_Prob, Edge_Prob>
{
public: 
	/** construct a graph with nv_ vertcies and ne_ edges where vertices and edges are added larter on */
	ProblemGraph(
		const int nv_, /**< the number of vertices */
		const int ne_  /**< the number of edges */
	);

	/** copy constructor */
	ProblemGraph(
		const ProblemGraph & probgraph_ /**< the graph to copy */
	);

	/** add a vertex with index and cost*/
	void addVertex(
		const int v_ind_, /**< vertex index*/
		const SCIP_Real v_cost_ /**< vertex cost*/
		);

	/** add an edge with index and cost */
	void addEdge(
		const int e_ind, /**< edge index*/
		const int tail,  /**< vertex index of edge tail */
		const int head, /**< vertex index of edge head */
		const SCIP_Real e_cost, /**< edge cost */
		const SCIP_Real e_capacity/**< edge capacity */
	);

	/**< return the user-defined cost on edge e*/
	SCIP_Real getCost(
		const int e_ind_  /**< vertex index*/
	);
	
	int nq; /**< the number of opportunity quaternions*/
	vector<Quaternion> quat_array;  /**< array stores opportunity quaternions of the graph */
	SCIP_Real flow_val; /**< current flow value*/
};



/** ConflictClique class */
class ConflictClique {
public:
	/** construct a conflict clique with ne_ edges where edges are added larter on */
	ConflictClique(
		const int c_ind_, /**< clique index */
		const int ne_  /**< the number of edges */
	);


	/** add an edge by index */
	void addEdge(
		const int e_ind_ /**< edge index */
	);

	int c_ind; /**< clique index */
	int ne; /**< the number of edges in the conflict clique */
	vector<int> e_ind_array; /**< array stores the indices of edges in the conflict clique */

};


/** Conflict class */
class Conflict {
public:
	/**< construct a coonflict management class */
	Conflict(
		const ProblemGraph * probgraph, /**< the pointer to the problem graph */
		const int nc_  /**< the number ofconflict cliques */
	);

	/** copy constructor */
	Conflict(
		const Conflict & conflict_ /**< the conflict to copy*/
	);


	/** add an clique by index and its size*/
	void addClique(
		const int c_ind_, /**< clique index */
		const int ne_,  /**< the number of edges */
		const vector<int> & e_array /**< edges in the clique */
	);

	/** return the clique indices of an edge */
	vector<int> * getCINDS(
		const int e_ind_ /**< edge index */
	);	

	ProblemGraph * probgraph; /**< pointer to the problem_graph */
	int nc; /**< the number of conflict cliques, the index starts from 0 to nc-1 */
	vector<ConflictClique> c_array; /**< array stores the indices of conflict cliques */
	vector<vector<int>> edge_to_c_inds; /**< vector stores the index array of belonging conflict cliques for each edge */
};
#endif