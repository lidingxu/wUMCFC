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

/**@file   ProbDataNTRT.cpp
 * @brief  C++ problem data for NTRT
 * @author Liding XU
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <algorithm>
#include "objscip/objscip.h"
#include "scip/struct_cons.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"
#include "scip/var.h"
#include "ProbDataNTRT.h"
#include "check.h"


using namespace scip;
using namespace std;



/** ProbDataNTRT destructor */
ProbDataNTRT::~ProbDataNTRT()
{
	// delete graph 
	if (problemgraph != NULL) {
		delete problemgraph;
	}
	// delete conflict
	if (conflict != NULL) {
		delete conflict;
	}
}

/** release scip reference in probelme data*/
SCIP_RETCODE ProbDataNTRT::releaseAll(
	SCIP*              scip                /**< SCIP data structure */
) {
	// release edge varibles
	for(int i = 0; i < nf; i++){
		for(int j = 0; j < problemgraph->ne; j++){
			SCIP_CALL(SCIPreleaseVar(scip, &e_vars[i][j].e_var));
		}
	}

	// release vertex variables
	for (int i = 0; i < v_vars.size(); i++) {
		SCIP_CALL(SCIPreleaseVar(scip, &v_vars[i].v_var));
	}

    // release opportunity variables
	for(int i = 0; i < nf; i++){
		for(int j = 0; j < op_vars[i].size(); j++){
			SCIP_CALL(SCIPreleaseVar(scip, &op_vars[i][j].op_var));
		}
	}

	// release coding variables
	for (int i = 0; i < cod_vars.size(); i++) {
		SCIP_CALL(SCIPreleaseVar(scip, &cod_vars[i].cod_var));
	}

	// release flow conserve constraints
	for (int i = 0; i < nf; i++) {
		for(int j = 0; j < fc_conss[i].size(); j++){
			SCIP_CALL(SCIPreleaseCons(scip, &fc_conss[i][j].fc_cons));
		}
	}

	// release directed opportunity constraints
	for(int i = 0; i < nf; i ++){
		for(int j = 0; j < do_conss[i].size(); j++){
			SCIP_CALL(SCIPreleaseCons(scip, &do_conss[i][j].do_cons));
		}
	}

	// release coding oppportunity constraints
	for (int i = 0; i < co_conss.size(); i++) {
		SCIP_CALL(SCIPreleaseCons(scip, &co_conss[i].co_cons));
	}

	// release vertex constraints
	for (int i = 0; i < v_conss.size(); i++) {
		SCIP_CALL(SCIPreleaseCons(scip, &v_conss[i].v_cons));
	}

	// release clique constraints
	for (int i = 0; i < clq_conss.size(); i++) {
		SCIP_CALL(SCIPreleaseCons(scip, &clq_conss[i].clq_cons));
	}
	return SCIP_OKAY;
}
/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbDataNTRT::scip_delorig(
   SCIP*              scip                /**< SCIP data structure */
   )
{
	SCIP_CALL(releaseAll(scip));
	return SCIP_OKAY;
}

/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbDataNTRT::scip_deltrans(
   SCIP*              scip                /**< SCIP data structure */
   )
{ 
   SCIP_CALL(releaseAll(scip));
   return SCIP_OKAY;
}


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
SCIP_RETCODE ProbDataNTRT::scip_trans(
   SCIP*              scip,               /**< SCIP data structure */
   ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
   SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
   )
{  /*lint --e{715}*/
   assert( objprobdata != NULL );
   assert( deleteobject != NULL );

   assert( problemgraph != NULL );
   assert( conflict != NULL);

   // create and cpature transformed path varibles 
   SCIPdebugMessage("start transform \n");
   ProblemGraph * transgraph = new ProblemGraph(*problemgraph);
   Conflict * transconflict = new Conflict(*conflict);


   // allocate memory for target prob data

	ProbDataNTRT * transprobdata = new ProbDataNTRT(transgraph, transconflict, nf,
	   flows, st_pairs, co_cons_map, v_cons_map, cod_cons_to_codvar);
 
  // transform and capture transformed flow conservative constraints
  for(int i = 0; i < fc_conss.size(); i++){
	  transprobdata->fc_conss.push_back(vector<FCCons>());
	  for(int j = 0; j < fc_conss[i].size(); j++){
		  transprobdata->fc_conss[i].push_back(fc_conss[i][j]);
		  SCIP_CALL(SCIPtransformCons(scip, fc_conss[i][j].fc_cons, &transprobdata->fc_conss[i][j].fc_cons));
	  }
  }

  // transfrom and capture transformed directed opportunity constraints
  for(int i = 0; i < do_conss.size(); i++){
	  transprobdata->do_conss.push_back(vector<DOCons>());
	  for(int j = 0; j < do_conss[i].size(); j++){
		  transprobdata->do_conss[i].push_back(do_conss[i][j]);
		  SCIP_CALL(SCIPtransformCons(scip, do_conss[i][j].do_cons, &transprobdata->do_conss[i][j].do_cons));
	  }
  }

  // transform and capture transformed coding oppportunity constraints
   for (int i = 0; i < co_conss.size(); i++) {
	   transprobdata->co_conss.push_back(co_conss[i]);
	   SCIP_CALL(SCIPtransformCons(scip, co_conss[i].co_cons, &transprobdata->co_conss[i].co_cons));
   }
   // transform and capture transformed vertex constraints
   for (int i = 0; i < v_conss.size(); i++) {
	   transprobdata->v_conss.push_back(v_conss[i]);
	   SCIP_CALL(SCIPtransformCons(scip, v_conss[i].v_cons, &transprobdata->v_conss[i].v_cons));
   }
   // transform and cpature transformed clique constraints
   for (int i = 0; i < clq_conss.size(); i++) {
	   transprobdata->clq_conss.push_back(clq_conss[i]);
	   SCIP_CALL(SCIPtransformCons(scip, clq_conss[i].clq_cons, &transprobdata->clq_conss[i].clq_cons));
   }

   // create and capture transformed edge variables
   for(int i = 0 ; i < e_vars.size(); i++){
	   transprobdata->e_vars.push_back(vector<EdgeVar>());
	   for(int j = 0; j < e_vars[i].size(); j++){
		   transprobdata->e_vars[i].push_back(e_vars[i][j]);
		   SCIP_CALL(SCIPtransformVar(scip, e_vars[i][j].e_var, &transprobdata->e_vars[i][j].e_var));
	   }
   }

   // create and capture transformed opportunity variables
   for(int i = 0; i < op_vars.size(); i++){
	   transprobdata->op_vars.push_back(vector<OpportunityVar>());
	   for(int j = 0; j < op_vars[i].size(); j++){
		   transprobdata->op_vars[i].push_back(op_vars[i][j]);
		   SCIP_CALL(SCIPtransformVar(scip, op_vars[i][j].op_var, &transprobdata->op_vars[i][j].op_var));
	   }
   }

   // create and cpature transformed vertex variables
   for (int i = 0; i < v_vars.size(); i++) {
	   transprobdata->v_vars.push_back(v_vars[i]);
	   SCIP_CALL(SCIPtransformVar(scip, v_vars[i].v_var, &(transprobdata->v_vars[i].v_var)));
   }
   // create and cpature transformed coding variables
   for (int i = 0; i < cod_vars.size(); i++) {
	   transprobdata->cod_vars.push_back(cod_vars[i]);
	   SCIP_CALL(SCIPtransformVar(scip, cod_vars[i].cod_var, &(transprobdata->cod_vars[i].cod_var)));
   }

   SCIPdebugMessage("end transform \n");
   assert( transprobdata != NULL );
   *objprobdata = transprobdata;           
   
   *deleteobject = TRUE;

   return SCIP_OKAY;
}      

/** create initial columns */
SCIP_RETCODE ProbDataNTRT::createInitialColumns(
	SCIP*                 scip               /**< SCIP data structure */
) {
	// create edge variables
	e_vars.reserve(nf);
	fc_conss.reserve(nf);
	//printf("added..");
	for(int i = 0; i < nf; i++){
		e_vars.push_back(vector<EdgeVar>());
		fc_conss.push_back(vector<FCCons>());
		//SCIPdebugMessage("added");
		SCIP_CALL(addEdgeVarFCCons(scip, st_pairs[i].first, st_pairs[i].second, i));
	}

	// create opportunity variables
	op_vars.reserve(nf);
	for(int i = 0; i < nf ; i++){
		op_vars.push_back(vector<OpportunityVar>());
		for(int e1_ind = 0;  e1_ind < problemgraph->ne; e1_ind++){
			int head = problemgraph->e_array[e1_ind].e_head;
			if(head == st_pairs[i].first || head == st_pairs[i].second){
				continue;
			}
			for(int e2_ind: *problemgraph->getOutEdges(head)){
				SCIP_CALL(addOpVar(scip, e1_ind, e2_ind, i));
			}
		}
	} 

	// create directed opportunity constraint
	do_conss.reserve(nf);
	for(int i = 0; i < nf; i++){
		do_conss.push_back(vector<DOCons>());
		for(int e_ind = 0; e_ind < problemgraph->ne; e_ind++){
			SCIP_CALL(addDOCons(scip, i, e_ind));
		}
	}

	// create opportunity constraints and coding variables
	co_conss.reserve(problemgraph->nq);
	for (int i = 0; i < problemgraph->nq; i++) {
		//check_qt(problemgraph, problemgraph->quat_array[i]);
		SCIP_CALL(addCOConsVar(scip, problemgraph->quat_array[i].v1_ind,
			problemgraph->quat_array[i].v2_ind,
			problemgraph->quat_array[i].v3_ind,
			problemgraph->quat_array[i].e12_ind,
			problemgraph->quat_array[i].e23_ind,
			problemgraph->quat_array[i].e32_ind,
			problemgraph->quat_array[i].e21_ind,
			problemgraph->quat_array[i].rate,
			problemgraph->quat_array[i].redcost
		));
	}
	
	SCIPdebugMessage("--opportunity constraints and varibales created!\n");

	// create clique constraint with all coding varibales
	clq_conss.reserve(conflict->nc);
	for (int i = 0; i < conflict->nc; i++) {
		SCIP_CALL(addClqCons(scip, i));
	}

	SCIPdebugMessage("--clique constraint created!\n");
	   
	// create vertex variables and constraints
	vector<Vertex_Prob> & v_array = problemgraph->v_array;
	v_conss.reserve(v_array.size());
	v_cons_map = vector<int>(v_array.size(), -1);
	for (int i = 0; i < v_array.size(); i++) {
		if (!isSourceTarget(i)) {
			//SCIP_CALL(addVertexConsVar(scip, i, v_array[i].v_cost));
		}
	}

	return SCIP_OKAY;
}

/** utility functions*/

/** add edge variables and flow conserve constraint */
SCIP_RETCODE ProbDataNTRT::addEdgeVarFCCons(
	   SCIP * scip,
	   const int s,
	   const int t,
	   const int f_id
){

	for(int i= 0; i < problemgraph->ne ; i++){
		// create edge variable 
		SCIP_VAR* e_var;
		char e_name[250];
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&e_var, /**<	pointer to edge object*/
			e_name, /**<	name of variable, or NULL for automatic name creation*/
			0, /**< lower bound of variable*/
			1, /**< upper bound of variable*/
			flows[f_id] * problemgraph->e_array[i].e_cost, /**< objective function value*/
			SCIP_VARTYPE_BINARY, /**<	type of variable*/
			TRUE, /**<	should var's column be present in the initial root LP?*/
			FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));
		/* add variable to SCIP and to the problem data */
		SCIP_CALL(SCIPaddVar(scip, e_var));
		SCIP_CALL(SCIPcaptureVar(scip, e_var));
		e_vars[f_id].push_back(EdgeVar(i, f_id, e_var));
		SCIP_CALL(SCIPreleaseVar(scip, &e_var));
		//SCIPdebugMessage("e_var added");
	}
	
	for(int i = 0; i < problemgraph->nv ; i++){
		// create flow conserve constraint
		if(i == s || i == t){
			SCIP_CONS* in_con;
			SCIP_CONS* out_con;
			vector<SCIP_VAR *> in_e_var_ptrs;
			vector<SCIP_Real> in_coeffs;
			vector<SCIP_VAR *> out_e_var_ptrs;
			vector<SCIP_Real> out_coeffs;
			SCIP_Real ins, outs;
			for(int e_ind: *problemgraph->getInEdges(i)){
				in_e_var_ptrs.push_back(e_vars[f_id][e_ind].e_var);
				in_coeffs.push_back(1);
			}
			for(int e_ind: *problemgraph->getOutEdges(i)){
				out_e_var_ptrs.push_back(e_vars[f_id][e_ind].e_var);
				out_coeffs.push_back(1);
			}
			if(i == s){
				ins = 0;
				outs = 1;
			}
			else if(i == t){
				ins = 1;
				outs = 0;
			}
			SCIP_CALL(SCIPcreateConsLinear(
				scip, /**< SCIP data structure */
				&in_con, /**< 	pointer to hold the created constraint */
				"flow conserve constraint", /**< name of constraint */
				in_e_var_ptrs.size(), /**< number of nonzeros in the constraint */
				in_e_var_ptrs.data(), /**< array with variables of constraint entries */
				in_coeffs.data(), /**< array with coefficients of constraint entries */
				ins,					/* lhs */
				ins,                    /* rhs */
				true,                   /* initial */
				false,                  /* separate */
				true,                   /* enforce */
				true,                   /* check */
				true,                   /* propagate */
				false,                  /* local */
				false,                   /* modifiable */
				false,                  /* dynamic */
				false,                  /* removable */
				false));               /* stickingatnode */
			SCIP_CALL(SCIPaddCons(scip, in_con));
			SCIP_CALL(SCIPcaptureCons(scip, in_con));
			fc_conss[f_id].push_back(FCCons(in_con, f_id, i));
			SCIP_CALL(SCIPreleaseCons(scip, &in_con));
			SCIP_CALL(SCIPcreateConsLinear(
				scip, /**< SCIP data structure */
				&out_con, /**< 	pointer to hold the created constraint */
				"flow conserve constraint", /**< name of constraint */
				out_e_var_ptrs.size(), /**< number of nonzeros in the constraint */
				out_e_var_ptrs.data(), /**< array with variables of constraint entries */
				out_coeffs.data(), /**< array with coefficients of constraint entries */
				outs,					/* lhs */
				outs,                    /* rhs */
				true,                   /* initial */
				false,                  /* separate */
				true,                   /* enforce */
				true,                   /* check */
				true,                   /* propagate */
				false,                  /* local */
				false,                   /* modifiable */
				false,                  /* dynamic */
				false,                  /* removable */
				false));               /* stickingatnode */
			SCIP_CALL(SCIPaddCons(scip, out_con));
			SCIP_CALL(SCIPcaptureCons(scip, out_con));
			fc_conss[f_id].push_back(FCCons(out_con, f_id, i));
			SCIP_CALL(SCIPreleaseCons(scip, &out_con));
		}
		else{
			SCIP_CONS* con;
			SCIP_CONS* con_1;
			SCIP_CONS* con_2;
			vector<SCIP_VAR *> e_var_ptrs, e_var_ptrs_1, e_var_ptrs_2;
			vector<SCIP_Real> coeffs, coeffs_1, coeffs_2;
			for(int e_ind: *problemgraph->getInEdges(i)){
				e_var_ptrs.push_back(e_vars[f_id][e_ind].e_var);
				coeffs.push_back(1);
				e_var_ptrs_1.push_back(e_vars[f_id][e_ind].e_var);
				coeffs_1.push_back(1);
			}
			for(int e_ind: *problemgraph->getOutEdges(i)){
				e_var_ptrs.push_back(e_vars[f_id][e_ind].e_var);
				coeffs.push_back(-1);
				e_var_ptrs_2.push_back(e_vars[f_id][e_ind].e_var);
				coeffs_2.push_back(1);
			}
			SCIP_CALL(SCIPcreateConsLinear(
				scip, /**< SCIP data structure */
				&con, /**< 	pointer to hold the created constraint */
				"flow conserve constraint", /**< name of constraint */
				e_var_ptrs.size(), /**< number of nonzeros in the constraint */
				e_var_ptrs.data(), /**< array with variables of constraint entries */
				coeffs.data(), /**< array with coefficients of constraint entries */
				0.0,					/* lhs */
				0.0,                    /* rhs */
				true,                   /* initial */
				false,                  /* separate */
				true,                   /* enforce */
				true,                   /* check */
				true,                   /* propagate */
				false,                  /* local */
				false,                   /* modifiable */
				false,                  /* dynamic */
				false,                  /* removable */
				false));               /* stickingatnode */
			SCIP_CALL(SCIPaddCons(scip, con));
			SCIP_CALL(SCIPcaptureCons(scip, con));
			fc_conss[f_id].push_back(FCCons(con, f_id, i));
			SCIP_CALL(SCIPreleaseCons(scip, &con));

			SCIP_CALL(SCIPcreateConsLinear(
				scip, /**< SCIP data structure */
				&con_1, /**< 	pointer to hold the created constraint */
				"flow conserve constraint", /**< name of constraint */
				e_var_ptrs_1.size(), /**< number of nonzeros in the constraint */
				e_var_ptrs_1.data(), /**< array with variables of constraint entries */
				coeffs_1.data(), /**< array with coefficients of constraint entries */
				0.0,					/* lhs */
				1.0,                    /* rhs */
				true,                   /* initial */
				false,                  /* separate */
				true,                   /* enforce */
				true,                   /* check */
				true,                   /* propagate */
				false,                  /* local */
				false,                   /* modifiable */
				false,                  /* dynamic */
				false,                  /* removable */
				false));               /* stickingatnode */
			SCIP_CALL(SCIPaddCons(scip, con_1));
			SCIP_CALL(SCIPcaptureCons(scip, con_1));
			fc_conss[f_id].push_back(FCCons(con_1, f_id, i));
			SCIP_CALL(SCIPreleaseCons(scip, &con_1));
		}
	}
	//printf("fccons added");
	return SCIP_OKAY;
}

/** add opportunity variable */
SCIP_RETCODE ProbDataNTRT::addOpVar(
	SCIP * scip, 
   	const int e1_ind, 
	const int e2_ind,
	const int f_id)
{
	SCIP_VAR* op_var;
	SCIP_CALL(SCIPcreateVar(
		scip, /**<	SCIP data structure*/
		&op_var, /**<	pointer to edge object*/
		"op_var", /**<	name of variable, or NULL for automatic name creation*/
		0, /**< lower bound of variable*/
		SCIPinfinity(scip), /**< upper bound of variable*/
		0, /**<	objective function value*/
		SCIP_VARTYPE_CONTINUOUS, /**<	type of variable*/
		TRUE, /**<	should var's column be present in the initial root LP?*/
		FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
		NULL, NULL, NULL, NULL, NULL
		));
	/* add variable to SCIP and to the problem data */
	SCIP_CALL(SCIPaddVar(scip, op_var));
	SCIP_CALL(SCIPcaptureVar(scip, op_var));
	op_vars[f_id].push_back(OpportunityVar(e1_ind, e2_ind, f_id, op_var));
	SCIP_CALL(SCIPreleaseVar(scip, &op_var));
	return SCIP_OKAY;
}


/** add directed opportunity constraint */
SCIP_RETCODE ProbDataNTRT::addDOCons(
	SCIP * scip, 
	const int f_id, 
	const int e_ind)
{
	SCIP_CONS* con1;
	vector<SCIP_VAR *> e_var_ptrs;
	vector<SCIP_Real> coeffs;
	e_var_ptrs.push_back(e_vars[f_id][e_ind].e_var);
	coeffs.push_back(flows[f_id]);
	for(int i = 0; i < op_vars[f_id].size(); i++){
		if(op_vars[f_id][i].end_ind == e_ind){
			e_var_ptrs.push_back(op_vars[f_id][i].op_var);
			coeffs.push_back(-1);
		}
	}
	// not source or target vertex
	if(coeffs.size() > 1){
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con1, /**< 	pointer to hold the created constraint */
		"directed opportunity constraint", /**< name of constraint */
		e_var_ptrs.size(), /**< number of nonzeros in the constraint */
		e_var_ptrs.data(), /**< array with variables of constraint entries */
		coeffs.data(), /**< array with coefficients of constraint entries */
		0.0,					/* lhs */
		0.0,                    /* rhs */
		true,                   /* initial */
		false,                  /* separate */
		true,                   /* enforce */
		true,                   /* check */
		true,                   /* propagate */
		false,                  /* local */
		false,                   /* modifiable */
		false,                  /* dynamic */
		false,                  /* removable */
		false));               /* stickingatnode */
	SCIP_CALL(SCIPaddCons(scip, con1));
	SCIP_CALL(SCIPcaptureCons(scip, con1));
	do_conss[f_id].push_back(DOCons(con1, f_id, e_ind));
	SCIP_CALL(SCIPreleaseCons(scip, &con1));
	}
	

	SCIP_CONS* con2;
	e_var_ptrs.clear();
	coeffs.clear();
	e_var_ptrs.push_back(e_vars[f_id][e_ind].e_var);
	coeffs.push_back(flows[f_id]);
	for(int i = 0; i < op_vars[f_id].size(); i++){
		if(op_vars[f_id][i].start_ind == e_ind){
			e_var_ptrs.push_back(op_vars[f_id][i].op_var);
			coeffs.push_back(-1);
		}
	}
	if(coeffs.size() > 1){
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con2, /**< 	pointer to hold the created constraint */
		"directed opportunity constraint", /**< name of constraint */
		e_var_ptrs.size(), /**< number of nonzeros in the constraint */
		e_var_ptrs.data(), /**< array with variables of constraint entries */
		coeffs.data(), /**< array with coefficients of constraint entries */
		0.0,					/* lhs */
		0.0,                    /* rhs */
		true,                   /* initial */
		false,                  /* separate */
		true,                   /* enforce */
		true,                   /* check */
		true,                   /* propagate */
		false,                  /* local */
		false,                   /* modifiable */
		false,                  /* dynamic */
		false,                  /* removable */
		false));               /* stickingatnode */
	SCIP_CALL(SCIPaddCons(scip, con2));
	SCIP_CALL(SCIPcaptureCons(scip, con2));
	do_conss[f_id].push_back(DOCons(con2, f_id, e_ind));
	SCIP_CALL(SCIPreleaseCons(scip, &con2));
	}
	return SCIP_OKAY;
}

/** add coding opportunity constraint and variable */
SCIP_RETCODE ProbDataNTRT::addCOConsVar(
	SCIP * 	scip, /**< SCIP data structure */
	const int v1_ind_,
	const int v2_ind_,
	const int v3_ind_,
	const int e12_ind_,
	const int e23_ind_,
	const int e32_ind_,
	const int e21_ind_,
	const SCIP_Real rate_,
	const SCIP_Real redcost_
) {
	// create coding variable 
	SCIP_VAR* cod_var;
	SCIP_CALL(SCIPcreateVar(
		scip, /**<	SCIP data structure*/
		&cod_var, /**< 	pointer to variable object*/
		"cod_var", /**<	name of variable, or NULL for automatic name creation*/
		-SCIPinfinity(scip), /**< lower bound of variable*/
		SCIPinfinity(scip), /**< upper bound of variable*/
		-redcost_, /**<	objective function value*/
		SCIP_VARTYPE_CONTINUOUS, /**<	type of variable*/
		TRUE, /**<	should var's column be present in the initial root LP?*/
		FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
		NULL, NULL, NULL, NULL, NULL
	));
	/* add variable to SCIP and to the problem data */
	SCIP_CALL(SCIPaddVar(scip, cod_var));
	//SCIP_CALL(SCIPchgVarLbLazy(scip, cod_var, 0.0));
	// capture the variable
	SCIP_CALL(SCIPcaptureVar(scip, cod_var));
	cod_vars.push_back(CodingVar(v1_ind_, v2_ind_, v3_ind_,
		e12_ind_, e23_ind_, e32_ind_, e21_ind_, rate_, redcost_, cod_var));
	/* release variable for the reader. */

	// create coding opportunity constraint for edge pair e12, e23
	SCIP_CONS* con123;
	vector<SCIP_VAR *> e_var_ptrs;
	vector<SCIP_Real> coeffs;
	e_var_ptrs.push_back(cod_var);
	coeffs.push_back(1);
	for(int f_id = 0; f_id < nf; f_id++){
		for(int i = 0; i < op_vars[f_id].size(); i++){
			if(op_vars[f_id][i].start_ind == e12_ind_ && op_vars[f_id][i].end_ind == e23_ind_){
				e_var_ptrs.push_back(op_vars[f_id][i].op_var);
				coeffs.push_back(-1);
				break;
			}
		}
	}
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con123, /**< 	pointer to hold the created constraint */
		"coding constraint", /**< 	name of constraint */ 
		e_var_ptrs.size(), /**< number of nonzeros in the constraint */
		e_var_ptrs.data(), /**< array with variables of constraint entries */
		coeffs.data(), /**< array with coefficients of constraint entries */
		-SCIPinfinity(scip),    /* lhs */
		0.0,                    /* rhs */
		true,                   /* initial */
		false,                  /* separate */
		true,                   /* enforce */
		true,                   /* check */
		true,                   /* propagate */
		false,                  /* local */
		false,                   /* modifiable */
		false,                  /* dynamic */
		false,                  /* removable */
		false));               /* stickingatnode */
	SCIP_CALL(SCIPaddCons(scip, con123));
	SCIP_CALL(SCIPcaptureCons(scip, con123));
	co_conss.push_back(COCons(con123, e12_ind_, e23_ind_));
	co_cons_map[make_pair(e12_ind_, e23_ind_)] = co_conss.size() - 1;
	cod_cons_to_codvar[co_conss.size() - 1] = cod_vars.size() - 1;
	SCIP_CALL(SCIPreleaseCons(scip, &con123));

	// create coding opportunity constraint for edge pair e32, e21
	SCIP_CONS* con321;
	e_var_ptrs.clear();
	coeffs.clear();
	e_var_ptrs.push_back(cod_var);
	coeffs.push_back(1);
	for(int f_id = 0; f_id < nf; f_id++){
		for(int i = 0; i < op_vars[f_id].size(); i++){
			if(op_vars[f_id][i].start_ind == e32_ind_ && op_vars[f_id][i].end_ind == e21_ind_){
				e_var_ptrs.push_back(op_vars[f_id][i].op_var);
				coeffs.push_back(-1);
				break;
			}
		}
	}
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con321, /**< 	pointer to hold the created constraint */
		"coding constraint", /**< 	name of constraint */
		e_var_ptrs.size(), /**< number of nonzeros in the constraint */
		e_var_ptrs.data(), /**< array with variables of constraint entries */
		coeffs.data(), /**< array with coefficients of constraint entries */
		-SCIPinfinity(scip),    /* lhs */
		0.0,                    /* rhs */
		true,                   /* initial */
		false,                  /* separate */
		true,                   /* enforce */
		true,                   /* check */
		true,                   /* propagate */
		false,                  /* local */
		false,                   /* modifiable */
		false,                  /* dynamic */
		false,                  /* removable */
		false));               /* stickingatnode */
	SCIP_CALL(SCIPaddCons(scip, con321));
	SCIP_CALL(SCIPcaptureCons(scip, con321));
	co_conss.push_back(COCons(con321, e32_ind_, e21_ind_));
	co_cons_map[make_pair(e32_ind_, e21_ind_)] = co_conss.size() - 1;
	cod_cons_to_codvar[co_conss.size() - 1] = cod_vars.size() - 1;
	SCIP_CALL(SCIPreleaseCons(scip, &con321));

	SCIP_CALL(SCIPreleaseVar(scip, &cod_var));
	return SCIP_OKAY;
}

/** add vertex constraint and variable */
SCIP_RETCODE ProbDataNTRT::addVertexConsVar(
	SCIP * 	scip, /**< SCIP data structure */
	const int v_ind_, /**< vertex index */
	const SCIP_Real v_cost_ /**< vertex cost */
) {
	SCIP_VAR* v_var;
	SCIP_CALL(SCIPcreateVar(
		scip, /**< SCIP data structure*/
		&v_var, /**< pointer to variable object*/
		"vertex_var", /**< name of variable, or NULL for automatic name creation*/
		0.0, /**< lower bound of variable*/
		1.0, /**< upper bound of variable*/
		v_cost_, /**< objective function value*/
		SCIP_VARTYPE_BINARY, /**< type of variable*/
		TRUE, /**< should var's column be present in the initial root LP?*/
		FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
		NULL, NULL, NULL, NULL, NULL
	));
	/* add variable to SCIP and to the problem data */
	SCIP_CALL(SCIPaddVar(scip, v_var));
	// capture the variable
	SCIP_CALL(SCIPcaptureVar(scip, v_var));
	v_vars.push_back(VertexVar(v_ind_, v_var));
	
	// create vertex constraint
	SCIP_CONS* con;
	vector<SCIP_VAR *> e_var_ptrs;
	vector<SCIP_Real> coeffs;
	e_var_ptrs.push_back(v_var);
	coeffs.push_back(-1);
	for(int f_id = 0; f_id < nf; f_id++){
		for(int i = 0; i < problemgraph->ne; i++){ 
			if(problemgraph->e_array[i].e_head == v_ind_){
				e_var_ptrs.push_back(e_vars[f_id][i].e_var);
				//coeffs.push_back(flows[f_id]/(2*problemgraph->e_array[i].e_capacity));
				coeffs.push_back( 0.99 / nf);
			}
		}
	}
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con, /**< 	pointer to hold the created constraint */
		"vertex constraint", /**< 	name of constraint */
		e_var_ptrs.size(), /**< number of nonzeros in the constraint */
		e_var_ptrs.data(), /**< array with variables of constraint entries */
		coeffs.data(), /**< array with coefficients of constraint entries */
		-SCIPinfinity(scip),    /* lhs */
		0.0,                    /* rhs */
		true,                   /* initial */
		false,                  /* separate */
		true,                   /* enforce */
		true,                   /* check */
		true,                   /* propagate */
		false,                  /* local */
		false,                   /* modifiable */
		false,                  /* dynamic */
		false,                  /* removable */
		false));               /* stickingatnode */
	SCIP_CALL(SCIPaddCons(scip, con));
	SCIP_CALL(SCIPcaptureCons(scip, con));
	v_conss.push_back(VCons(con, v_ind_));
	v_cons_map[v_ind_] = v_conss.size() - 1;
	SCIP_CALL(SCIPreleaseCons(scip, &con));

	/* release variable */
	SCIP_CALL(SCIPreleaseVar(scip, &v_var));

	return SCIP_OKAY;
}

/** add clique constraints with all coding varibales */
SCIP_RETCODE ProbDataNTRT::addClqCons(
	SCIP * 	scip, /**< SCIP data structure */
	const int c_ind_ /**< clique index */
) {
	// create clique constraint
	vector<int> & e_ind_array = conflict->c_array[c_ind_].e_ind_array;
	// sort the edge indices belonging to the clique
	sort(e_ind_array.begin(), e_ind_array.end());
	vector<SCIP_VAR *> e_var_ptrs;
	vector<SCIP_Real> coeffs;
	for (int i = 0; i < cod_vars.size(); i++) {
		// check coding varibale is in the clique
		if (binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e12_ind) &&
			binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e23_ind) &&
			binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e32_ind) &&
			binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e21_ind)) {
			// add
			e_var_ptrs.push_back(cod_vars[i].cod_var);
			coeffs.push_back(-cod_vars[i].rate);
		}
	}
	SCIP_CONS* con;
	for(int i: conflict->c_array[c_ind_].e_ind_array){
		for(int f_id = 0; f_id < nf; f_id++){
			e_var_ptrs.push_back(e_vars[f_id][i].e_var);
			coeffs.push_back( flows[f_id] / problemgraph->e_array[i].e_capacity);
		}
	}
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con, /**< 	pointer to hold the created constraint */
		"clique constraint", /**< name of constraint */
		e_var_ptrs.size(), /**< number of nonzeros in the constraint */
		e_var_ptrs.data(), /**< array with variables of constraint entries */
		coeffs.data(), /**< array with coefficients of constraint entries */
		-SCIPinfinity(scip),    /* lhs */
		1.0,                    /* rhs */
		true,                   /* initial */
		false,                  /* separate */
		true,                   /* enforce */
		true,                   /* check */
		true,                   /* propagate */
		false,                  /* local */
		true,                   /* modifiable */
		false,                  /* dynamic */
		false,                  /* removable */
		false));               /* stickingatnode */
	SCIP_CALL(SCIPaddCons(scip, con));
	SCIP_CALL(SCIPcaptureCons(scip, con));
	clq_conss.push_back(ClqCons(con, c_ind_));
	SCIP_CALL(SCIPreleaseCons(scip, &con));
	return SCIP_OKAY;
}

/** check the vertex is a source or a target*/
bool ProbDataNTRT::isSourceTarget(
	const int v_ind_ /**< vertex index */
	) {
	for (auto pair : st_pairs) {
		if (pair.first == v_ind_ || pair.second == v_ind_) {
			return TRUE;
		}
	}
	return FALSE;
}