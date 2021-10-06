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

/**@file   Graph.h
 * @brief  Graph class 
 * @author Liding XU
 */

 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <iostream>
#include <utility>
#include <map>
#include <vector>

#include "objscip/objscip.h"
#include "scip/scip_numerics.h"

using namespace std;

/** Graph class */
template<typename Vertex_T, typename Edge_T> class Graph
{
public:
	Graph(){
		nv = 0;
		ne = 0;
	};
	/** construct a graph with nv_ vertcies and ne_ edges where edges is added later on */
	Graph(
		const int nv_, /**< the number of vertices */
		const int ne_  /**< the number of edges*/
	) : nv(nv_), ne(ne_) {
		/* reserve the memory */
		v_array.reserve(nv_);
		e_array.reserve(ne_);
		v_to_in_edges = vector<vector<int>>(nv_);
		v_to_out_edges = vector<vector<int>>(nv_);
	};

	/** copy constructor */
	Graph(
		const Graph & graph_ /**< the graph to copy */
	) {
		nv = graph_.nv;
		ne = graph_.ne;
		v_array = graph_.v_array;
		e_array = graph_.e_array;
		v_to_in_edges = graph_.v_to_in_edges;
		v_to_out_edges = graph_.v_to_out_edges;
		edge_map = graph_.edge_map;
	};

	/** return the pointer to the edge by index e_ind_ */
	inline Edge_T * getEdgeByIND(
		const int e_ind_  /**< vertex index*/
	) {
		return &e_array[e_ind_];
	};

	/** return the pointer to the edge by index e_tail_, e_head_, if no such edge return NULL */
	inline Edge_T * getEdgeByEnds(
		const int e_tail, /**< vertex index of edge tail */
		const int e_head /**< bvertex index of edge head */
	) {
		auto it = edge_map.find(make_pair(e_tail, e_head));
		if (it == edge_map.end()) {
			return NULL;
		}
		else {
			return &e_array[it->second];
		}
	};

	/** return out edges of vertex v_ind_*/
	inline vector<int> * getOutEdges(
		const int v_ind_ /**< vertex index*/
	) {
		return &v_to_out_edges[v_ind_];
	};

	/** return in edges of vertex v_ind_*/
	inline vector<int> * getInEdges(
		const int v_ind_ /**< vertex index*/
	) {
		return &v_to_in_edges[v_ind_];
	};

	/**< return the user-defined cost on edge e*/
	inline SCIP_Real getCost(
		const int e_ind_  /**< vertex index*/
	) {
		return e_array[e_ind_].e_cost;
	};


	/** Bellman-Ford algorithm to find shortest path from s to t*/
	SCIP_Real BellmanFord(
		int s, /**< source vertex index*/
		int t,  /**< target vertex index*/
		vector<int>& path /**< shortest path array from s to t*/
	) {
		vector<SCIP_Real> dists(nv, SCIP_DEFAULT_INFINITY);
		vector<int> preds(nv, -1);

		dists[s] = 0;

		for (int i = 0; i < nv - 1; i++) {
			for (int e_ind = 0; e_ind < ne; e_ind++) {
				int tail = e_array[e_ind].e_tail;
				int head = e_array[e_ind].e_head;
				SCIP_Real cost = getCost(e_ind);
				if (dists[tail] + cost < dists[head]) {
					dists[head] = dists[tail] + cost;
					preds[head] = tail;
				}
			}
		}
		path.clear();
		vector<int> rev_path;
		int v = t;
		while (v != s) {
			rev_path.push_back(v);
			v = preds[v];
		}
		rev_path.push_back(s);
		for(int i = rev_path.size() - 1; i >= 0; i--){
			path.push_back(rev_path[i]);
		}
		return dists[t];
	};

	int nv; /**< the number of vertices, the index starts from 0 to nv-1 */
	int ne; /**< the number of edges, the index starts from 0 to ne-1 */
	vector<Vertex_T> v_array; /**< array stores vertices of the graph */
	vector<Edge_T> e_array; /**< array stores edges of the graph */
	map<pair<int, int>, int> edge_map; /**< map to edge index by its tail and head vertex indice */
	vector<vector<int>> v_to_in_edges; /**< vector store the index array of in edges for each vertex */
	vector<vector<int>> v_to_out_edges; /**< vector store the index array of in edges for each vertex */
};

#endif