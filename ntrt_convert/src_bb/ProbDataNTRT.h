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

class EdgeVar{
public:
	/* default constructor */
	EdgeVar(
		const int e_ind_, /**< edge index >*/
		const int f_ind_, /**< flow index>*/
		SCIP_VAR * e_var_ /**< pointer to the edge var >*/	
	) : e_ind(e_ind_), f_ind(f_ind_), e_var(e_var_){};

	int e_ind; /**< edge index >*/
	int f_ind; /**< flow index>*/
	SCIP_VAR * e_var; /**< pointer to the edge var >*/	
};

class OpportunityVar{
public:
	/* default constructor */
	OpportunityVar(
		const int start_ind_, /**< edge e1 index >*/
		const int end_ind_, /**< edge e2 index >*/
		const int f_ind_, /**< flow index>*/
		SCIP_VAR * op_var_ /**< pointer to the edge var >*/	
	) : start_ind(start_ind_), end_ind(end_ind_), f_ind(f_ind_), op_var(op_var_){};

	int start_ind; /**< edge e1 index >*/
	int end_ind; /**< edge e2 index >*/
	int f_ind; /**< flow index>*/
	SCIP_VAR * op_var; /**< pointer to the edge var >*/	
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
/** flow conserve constraint */
class FCCons{
public:
	/** default constructor */
	FCCons(
		SCIP_CONS * fc_cons_, /**< pointer to the flow conserve constraint */
		const int flow_id_, /**< flow index id */ 
		const int vertex_id_ /**< vertex id */
	): fc_cons(fc_cons_), flow_id(flow_id_), vertex_id(vertex_id_){};

	SCIP_CONS * fc_cons; /**< pointer to the flow conserve constraint */
	int flow_id; /**< flow index id */ 
	int vertex_id; /**< vertex id */
};

/** directed opportunity constraint */
class DOCons{
public:
	/** default constructor */
	DOCons(
		SCIP_CONS * do_cons_, /**< pointer to directed opportunity constraint */
		const int flow_id_, /**< flow index id */ 
		const int e_ind_ /**< end coding edge */
	): do_cons(do_cons_), flow_id(flow_id_), e_ind(e_ind_){};

	SCIP_CONS * do_cons; /**< pointer to directed opportunity constraint */
	int flow_id; /**< flow index id */ 
	int e_ind; /**< end coding edge */
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

   /** add edge variables and flow conserve constraint */
   SCIP_RETCODE addEdgeVarFCCons(
	   SCIP * scip,
	   const int s,
	   const int t,
	   const int f_id
   );

	/** add Opportunity variable */
    SCIP_RETCODE addOpVar(
		SCIP * scip, 
   		const int e1_ind, 
		const int e2_ind,
		const int f_id);

	/** add directed opportunity constraint */
    SCIP_RETCODE addDOCons(
		SCIP * scip, 
   		const int f_id, 
		const int e_ind);

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

   /** add clique constraints with all coding varibales */
   SCIP_RETCODE addClqCons(
	   SCIP * 	scip, /**< SCIP data structure */
	   const int c_ind_ /**< clique index */
   );
      
   /** check the vertex is a source or a target*/
   bool isSourceTarget(
	   const int v_ind_ /**< vertex index */
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
	vector<pair<int, int>> st_pairs; /**< source and target pairs */

	vector<VertexVar> v_vars; /**< vertex variables */
	vector<vector<EdgeVar>> e_vars; /**< edge variables */
	vector<vector<OpportunityVar>> op_vars; /**< opportunity variables */
	vector<CodingVar> cod_vars; /**< coding variables */

	vector<vector<FCCons>> fc_conss; /**< flow conserve constraints */
	vector<vector<DOCons>> do_conss; /**< directed opportunity constraints */
	vector<COCons> co_conss; /**< coding opportunity constraints */
	vector<VCons> v_conss; /**< vertex constraints */
	vector<ClqCons> clq_conss; /**< clique constraints */

	map<pair<int, int>, int> co_cons_map; /**< map edge pairs to its coding opportunity contraint index*/
	map<int, int> cod_cons_to_codvar; /** map the coding constraints to the coding variable involved */
	vector<int> v_cons_map; /**< map vertex to its vertex contraint index*/
};/*lint !e1712*/


#endif
