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

/**@file   ProbDataNTRT.h
 * @brief  C++ problem data for TSP
 * @author Liding XU
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __NTRTPROBDATA_H__
#define __NTRTPROBDATA_H__

#include "objscip/objscip.h"
#include "scip/cons_linear.h"
#include "ProblemGraph.h"
#include <map>
#include <list>
#include <utility>

using namespace scip;

/** Variable classes */
class Path {
public:
	/** default constructor */
	Path(
		const int p_class_, /**< path class index */
		const SCIP_Real p_cost_, /**< path cost */
		const vector<int>& v_array_, /**< vertices in the path */
		const vector<int>& e_array_ /**< edges in the path */
	) : p_class(p_class_), p_cost(p_cost_), v_array(v_array_), e_array(e_array_) {};


	int p_class; /**< path class index */
	SCIP_Real p_cost; /**< path cost */
	vector<int> v_array; /**< vertices in the path */
	vector<int> e_array; /**< edges in the path */
};

class PathVar:public Path {
public:
	/** default constructor */
	PathVar(
		const int p_class_, /**< path class index */
		const SCIP_Real p_cost_, /**< path cost */
		const vector<int>& v_array_, /**< vertices in the path */
		const vector<int>& e_array_, /**< edges in the path */
		SCIP_VAR * p_var_ /**< pointer to the path var */
	): Path(p_class_, p_cost_, v_array_, e_array_), p_var(p_var_) {};

	/** constructor by the base Path class */
	PathVar(
		const Path & path, /**< base path reference */
		SCIP_VAR * p_var_ /**< pointer to the path var */
	) : Path::Path(path), p_var(p_var_) {};
	
	SCIP_VAR * p_var; /**< pointer to the path var */
};

class VertexVar{
public:
	/** default constructor */
	VertexVar(
		const int v_ind_, /**< vertex index */
		SCIP_VAR * v_var_ /**< pointer to the vertex var */
	) : v_ind(v_ind_), v_var(v_var_) {};



	SCIP_VAR * v_var; /**< pointer to the vertex var */
	int v_ind; /**< vertex index */
};

class CodingVar {
public:
	/** default constructor */
	CodingVar(
		const int v1_ind_,
		const int v2_ind_,
		const int v3_ind_,
		const int e12_ind_,
		const int e23_ind_,
		const int e32_ind_,
		const int e21_ind_,
		const SCIP_Real rate_,  /**< communication rate of coding channles */
		const SCIP_Real redcost_,  /**< reduced communication cost of broadcasting */
		SCIP_VAR * cod_var_ /**< pointer to the coding var */
	) :v1_ind(v1_ind_), v2_ind(v2_ind_), v3_ind(v3_ind_), e12_ind(e12_ind_),
		e23_ind(e23_ind_), e32_ind(e32_ind_), e21_ind(e21_ind_), rate(rate_), redcost(redcost_), cod_var(cod_var_) {};

	SCIP_VAR * cod_var; /**< pointer to the coding var */
	SCIP_Real rate;  /**< communication rate of coding channles */
	SCIP_Real redcost;  /**< reduced communication cost of broadcasting */
	int v1_ind; /**< coding vertex v1 index */
	int v2_ind; /**< coding vertex v2 index */
	int v3_ind; /**< coding vertex v3 index */
	int e12_ind, e23_ind, e32_ind, e21_ind;  /**< opportunity quaternion edges: (v1,v2) ,(v2,v3), (v3,v2), (v2,v1)*/
};

/** coding opportunity constraint */
class COCons {
public:
	/** default constructor */
	COCons(
		SCIP_CONS* co_cons_,  /**< pointer to hold the created constraint */
		const int start_e_, /**< start coding edge */
		const int end_e_ /**< end coding edge */
	): co_cons(co_cons_), start_e(start_e_), end_e(end_e_){};

	SCIP_CONS* co_cons;  /**< pointer to hold the created constraint */
	int start_e; /**< start coding edge */
	int end_e; /**< end coding edge */
};

/** vertex constraint */
class VCons {
public:
	/** default constructor */
	VCons(
		SCIP_CONS* v_cons_,  /**< pointer to hold the created constraint */
		const int v_ind_ /**< vertex index */
	) : v_cons(v_cons_), v_ind(v_ind_) {};

	SCIP_CONS* v_cons;  /**< pointer to hold the created constraint */
	int v_ind; /**< vertex index */
};

/** unsplittable constraint */
class UPCons {
public:
	/** default constructor */
	UPCons(
		SCIP_CONS* up_cons_,  /**< pointer to hold the created constraint */
		const int f_ind_, /**< flow index */
		const int s_, /**< source vertex */
		const int t_ /**< target vertex */
	) : up_cons(up_cons_), f_ind(f_ind_), s(s_), t(t_) {};

	SCIP_CONS* up_cons;  /**< pointer to hold the created constraint */
	int f_ind; /**< flow index */
	int s; /**< source vertex */
	int t; /**< target vertex */
};

/** clique constraint */
class ClqCons {
public:
	/** default constructor */
	ClqCons(
		SCIP_CONS* clq_cons_,  /**< pointer to hold the created constraint */
		const int c_ind_ /**< clique index */
	) : clq_cons(clq_cons_), c_ind(c_ind_) {};

	SCIP_CONS* clq_cons;  /**< pointer to hold the created constraint */
	int c_ind; /**< clique index */
};

/** Variable data which is attached to path variables.
 */
struct SCIP_VarData
{
   int id; // -1 path var, otherwise vertex id
   void* var;
};


extern SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
);

extern SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int id, /**< -1 path var, otherwise vertex id */
   void* var
);

/** SCIP user problem data for NTRT */
class ProbDataNTRT : public ObjProbData
{
public:

   /** default constructor */
   ProbDataNTRT(
		ProblemGraph * problemgraph_,   /**< problem graph data */
		Conflict * conflict_,  /**< conflict clique data */
		const int nf_,  /**< the number of flows*/
		const vector<SCIP_Real> & flows_, /**< demand flows*/
		const vector<pair<int, int>> & st_pairs_ /**< sourcer and target pairs*/
      )
      : problemgraph(problemgraph_), conflict(conflict_), nf(nf_), 
	   flows(flows_), st_pairs(st_pairs_){}

   /** copy constructor */
   ProbDataNTRT(
	   ProblemGraph * problemgraph_,   /**< problem graph data */
	   Conflict * conflict_,  /**< conflict clique data */
	   const int nf_,  /**< the number of flows*/
	   const vector<SCIP_Real> & flows_, /**< demand flows*/
	   const vector<pair<int, int>> & st_pairs_, /**< sourcer and target pairs*/
	   map<pair<int, int>, int> & co_cons_map_, /**< map edge pairs to its coding opportunity contraint index*/
	   vector<int> & v_cons_map_, /**< map vertex to its vertex contraint index*/
		 map<int, int> & cod_cons_to_codvar_ /** map the coding constraints to the coding variable involved */
   ) : problemgraph(problemgraph_), conflict(conflict_), nf(nf_), flows(flows_), 
	   st_pairs(st_pairs_), co_cons_map(co_cons_map_), v_cons_map(v_cons_map_), cod_cons_to_codvar(cod_cons_to_codvar_) {}

   /**< destructor */
   ~ProbDataNTRT();

   /**< utility functions*/

   /** release scip reference in probelme data*/
   SCIP_RETCODE releaseAll(
	   SCIP*              scip                /**< SCIP data structure */
   );

   /** create initial columns */
   SCIP_RETCODE createInitialColumns(
	   SCIP*                 scip               /**< SCIP data structure */
   );

   /** add coding opportunity constraint */
   SCIP_RETCODE addCOConsVar(
	   SCIP * 	scip, /**< SCIP data structure */
	   const int v1_ind_,
	   const int v2_ind_,
	   const int v3_ind_,
	   const int e12_ind_,
	   const int e23_ind_,
	   const int e32_ind_,
	   const int e21_ind_,
	   const SCIP_Real rate_,
	   const SCIP_Real cost_
   );

   /** add vertex variable and its data in scip and problem data*/
   SCIP_RETCODE addVertexConsVar(
	   SCIP * 	scip, /**< SCIP data structure */
	   const int v_ind_, /**< vertex index */
	   const SCIP_Real v_cost_ /**< path cost */
   );

   /** add unsplittable constraint */
   SCIP_RETCODE addUPCons(
	   SCIP * 	scip, /**< SCIP data structure */
	   const int f_ind_, /**< flow index */
	   const int s, /**< source index */
	   const int t /**< target index */
   );

   /** add clique constraints with all coding varibales */
   SCIP_RETCODE addClqCons(
	   SCIP * 	scip, /**< SCIP data structure */
	   const int c_ind_ /**< clique index */
   );

   /** add path variable and its data in scip and problem data*/
   SCIP_RETCODE addPathVar(
	   SCIP * 	scip, /**< SCIP data structure */
	   const int p_class_, /**< path class index */
	   const SCIP_Real p_cost_, /**< path cost */
	   const vector<int>& v_array_, /**< vertices in the path */
	   const vector<int>& e_array_, /**< edges in the path */
	   const SCIP_Bool is_pricing /**< indicate pricing variable or not*/
   );
      
   /** check the vertex is a source or a target*/
   bool isSourceTarget(
	   const int v_ind_ /**< vertex index */
   );
   
   /** return the number of path variables */
   int getNumPathVars();

   /** check the flow is fixed */
   bool isFlowFixed(
	   int f /**< flow index */
   );

   /** destructor of user problem data to free original user data (called when original problem is freed)
    *
    *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_delorig(
      SCIP*              scip                /**< SCIP data structure */
      );
   
   /** destructor of user problem data to free transformed user data (called when transformed problem is freed)
    *
    *  If the "*deleteobject" flag in the scip_trans() method was set to TRUE, this method is not needed,
    *  because all the work to delete the user problem data can be done in the destructor of the user problem
    *  data object. If the "*deleteobject" flag was set to FALSE, and the user problem data object stays alive
    *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
    *  longer needed.
    */
   virtual SCIP_RETCODE scip_deltrans(
      SCIP*              scip                /**< SCIP data structure */
      );

   /** creates user data of transformed problem by transforming the original user problem data
    *  (called after problem was transformed)
    *
    *  The user has two possibilities to implement this method:
    *   1. Return the pointer to the original problem data object (this) as pointer to the transformed problem data
    *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
    *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
    *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
    *      solving process is terminated.
    *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
    *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
    *      destructor of the object if the transformed problem data is no longer needed.
    */
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP data structure */
      ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
      SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
      );


    int nf; /**< the number of flows*/
	ProblemGraph * problemgraph;   /**< problem graph data */
	Conflict * conflict;  /**< conflict clique data */
	vector<SCIP_Real> flows; /**< demand flows */
	vector<pair<int, int>> st_pairs; /**< sourcer and target pairs */

	vector<list<PathVar>> p_vars; /**< path variables for each flow */
	vector<VertexVar> v_vars; /**< vertex variables */
	vector<CodingVar> cod_vars; /**< coding variables */

	vector<COCons> co_conss; /**< coding opportunity constraints */
	vector<VCons> v_conss; /**< vertex constraints */
	vector<UPCons> up_conss; /**< unsplittable constraints */
	vector<ClqCons> clq_conss; /**< clique constraints */

	map<pair<int, int>, int> co_cons_map; /**< map edge pairs to its coding opportunity contraint index*/
	map<int, int> cod_cons_to_codvar; /** map the coding constraints to the coding variable involved */
	vector<int> v_cons_map; /**< map vertex to its vertex contraint index*/
	double t_add = 0;
	double t_var = 0;
	int path_add= 0;
	int former = 0;
	vector<int> cts;
};/*lint !e1712*/


#endif
