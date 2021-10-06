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


SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
   )
{
   SCIPfreeBlockMemory(scip, vardata);
   return SCIP_OKAY;
};

SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int id, // -1 path var, otherwise vertex id
   void* var
)
{
   SCIP_CALL( SCIPallocBlockMemory(scip, vardata));
   (*vardata)->id = id;
   (*vardata)->var = var;
   return SCIP_OKAY;
};

SCIP_DECL_VARTRANS(vardataTrans){
	assert(sourcedata != NULL);
	assert(sourcevar != NULL);
	SCIP_CALL(vardataCreate(scip, targetdata, sourcedata->id, sourcedata->var));
	//SCIPdebugMessage("yes");
	return SCIP_OKAY;
};

SCIP_DECL_VARDELORIG(vardataDelOrig)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
};

SCIP_DECL_VARDELTRANS(vardataDelTrans)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
};




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
	// release path varibles
	for(int i = 0; i < nf; i++){
		for (auto it = p_vars[i].begin(); it != p_vars[i].end(); it++) {
			SCIP_CALL(SCIPreleaseVar(scip, &it->p_var));
		}
	}

	// release vertex variables
	for (int i = 0; i < v_vars.size(); i++) {
		SCIP_CALL(SCIPreleaseVar(scip, &v_vars[i].v_var));
	}

	// release coding variables
	for (int i = 0; i < cod_vars.size(); i++) {
		SCIP_CALL(SCIPreleaseVar(scip, &cod_vars[i].cod_var));
	}

	// release coding oppportunity constraints
	for (int i = 0; i < co_conss.size(); i++) {
		SCIP_CALL(SCIPreleaseCons(scip, &co_conss[i].co_cons));
	}

	// release vertex constraints
	for (int i = 0; i < v_conss.size(); i++) {
		SCIP_CALL(SCIPreleaseCons(scip, &v_conss[i].v_cons));
	}

	// release unsplittable constraints
	for (int i = 0; i < up_conss.size(); i++) {
		SCIP_CALL(SCIPreleaseCons(scip, &up_conss[i].up_cons));
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
   SCIPdebugMessage("start transform !!!!!!!!!!!\n");
   ProblemGraph * transgraph = new ProblemGraph(*problemgraph);
   Conflict * transconflict = new Conflict(*conflict);

   // allocate memory for target prob data



	ProbDataNTRT * transprobdata = new ProbDataNTRT(transgraph, transconflict, nf,
	   flows, st_pairs, co_cons_map, v_cons_map, cod_cons_to_codvar);

	transprobdata->cts = cts;

      // transform and cpature transformed coding oppportunity constraints
   for (int i = 0; i < co_conss.size(); i++) {
	   transprobdata->co_conss.push_back(COCons(co_conss[i]));
	   SCIP_CALL(SCIPtransformCons(scip, co_conss[i].co_cons, &transprobdata->co_conss[i].co_cons));
   }
   // transform and cpature transformed vertex constraints
   for (int i = 0; i < v_conss.size(); i++) {
	   transprobdata->v_conss.push_back(VCons(v_conss[i]));
	   SCIP_CALL(SCIPtransformCons(scip, v_conss[i].v_cons, &transprobdata->v_conss[i].v_cons));
   }
   // transform and cpature transformed unsplittable constraints
   for (int i = 0; i < up_conss.size(); i++) {
	   transprobdata->up_conss.push_back(UPCons(up_conss[i]));
	   SCIP_CALL(SCIPtransformCons(scip, up_conss[i].up_cons, &transprobdata->up_conss[i].up_cons));
   }

   // transform and cpature transformed clique constraints
   for (int i = 0; i < clq_conss.size(); i++) {
	   transprobdata->clq_conss.push_back(ClqCons(clq_conss[i]));
	   SCIP_CALL(SCIPtransformCons(scip, clq_conss[i].clq_cons, &transprobdata->clq_conss[i].clq_cons));
   }

   // create and cpature transformed path varibles 
   SCIPdebugMessage("start transform \n");
   transprobdata->p_vars = vector<list<PathVar>>(p_vars.size());
   for (int i = 0; i < p_vars.size(); i++) {
	   for (auto it = p_vars[i].begin(); it != p_vars[i].end(); it++) {
		   transprobdata->p_vars[i].push_back(PathVar(*it));
		   SCIP_CALL(SCIPtransformVar(scip, it->p_var, &(transprobdata->p_vars[i].back().p_var)));
	   }
   }
   // create and cpature transformed vertex variables
   for (int i = 0; i < v_vars.size(); i++) {
	   transprobdata->v_vars.push_back(VertexVar(v_vars[i]));
	   SCIP_CALL(SCIPtransformVar(scip, v_vars[i].v_var, &(transprobdata->v_vars[i].v_var)));
   }
   // create and cpature transformed coding variables
   for (int i = 0; i < cod_vars.size(); i++) {
	   transprobdata->cod_vars.push_back(CodingVar(cod_vars[i]));
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
	// use network coding or not?
	SCIP_Bool use_coding;
	SCIP_CALL(SCIPgetBoolParam(scip,  "ntrt/netcoding", &use_coding));

	// create opportunity constraints and variables
	cts = vector<int>(conflict->nc, 0);
	co_conss.reserve(problemgraph->nq);
	for (int i = 0; i < problemgraph->nq; i++) {
		//check_qt(problemgraph, problemgraph->quat_array[i]);
		// if use coding, add coding opportunity  constraint and variables
		if(use_coding){
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
	}
	
	SCIPdebugMessage("--%ld opportunity constraints and varibales created!\n", co_conss.size());

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
	

	SCIPdebugMessage("--vertex varibales and constraints created!\n");

	// create unsplittable constraints
	up_conss.reserve(nf);
	for (int i = 0; i < nf; i++) {
		SCIP_CALL(addUPCons(scip, i, st_pairs[i].first, st_pairs[i].second));
	}
	
	SCIPdebugMessage("--unsplittable constraints creaed!\n");
   
	// create initial path variables
	p_vars = vector<list<PathVar>>(nf);
	for (int i = 0; i < nf; i++) {
		problemgraph->flow_val = flows[i];
		vector<int> v_array;
		// run Bell-Ford algorithm to find the path vertex array
		problemgraph->BellmanFord(st_pairs[i].first, st_pairs[i].second, v_array);
		vector<int> e_array;
		// get the edge array and compute the path cost
		SCIP_Real p_cost = 0;
		for (int j = 1; j < v_array.size(); j++) {
			int tail = v_array[j - 1];
			int head = v_array[j];
			Edge_Prob * e_ptr = problemgraph->getEdgeByEnds(tail, head);
			assert(e_ptr != NULL);
			e_array.push_back(e_ptr->e_ind);
			p_cost += e_ptr->e_cost;
		}
		// add the path variable
		assert(is_simple_path(v_array));
		SCIPdebugMessage("-flow %d\n", i);
		SCIP_CALL(addPathVar(scip, i, p_cost, v_array, e_array, FALSE));
	}
	SCIPdebugMessage("--initial columns created!\n");
	return SCIP_OKAY;
}

/** utility functions*/

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
	SCIP_Real coeff_ = 1;
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con123, /**< 	pointer to hold the created constraint */
		"coding constraint", /**< 	name of constraint */ 
		1, /**< number of nonzeros in the constraint */
		&cod_var, /**< 	array with variables of constraint entries */
		&coeff_, /**< array with coefficients of constraint entries */
		-SCIPinfinity(scip),    /* lhs */
		0.0,                    /* rhs */
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
	SCIP_CALL(SCIPaddCons(scip, con123));
	SCIP_CALL(SCIPcaptureCons(scip, con123));
	co_conss.push_back(COCons(con123, e12_ind_, e23_ind_));
	co_cons_map[make_pair(e12_ind_, e23_ind_)] = co_conss.size() - 1;
	cod_cons_to_codvar[co_conss.size() - 1] = cod_vars.size() - 1;
	SCIP_CALL(SCIPreleaseCons(scip, &con123));

	// create coding opportunity constraint for edge pair e32, e21
	SCIP_CONS* con321;
	SCIP_Real coeff = 1;
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con321, /**< 	pointer to hold the created constraint */
		"coding constraint", /**< 	name of constraint */
		1, /**< number of nonzeros in the constraint */
		&cod_var, /**< 	array with variables of constraint entries */
		&coeff, /**< array with coefficients of constraint entries */
		-SCIPinfinity(scip),    /* lhs */
		0.0,                    /* rhs */
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
	
	SCIP_VARDATA * vardata;
	SCIP_CALL(vardataCreate(scip, &vardata, v_vars.size() - 1, NULL));
	/* set callback functions */
  	SCIPvarSetData(v_var, vardata);
	SCIPvarSetDelorigData(v_var, vardataDelOrig);
	SCIPvarSetTransData(v_var, vardataTrans);
    SCIPvarSetDeltransData(v_var, vardataDelTrans);

	// create vertex constraint
	SCIP_CONS* con;
	SCIP_Real coeff = -1;
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con, /**< 	pointer to hold the created constraint */
		"vertex constraint", /**< 	name of constraint */
		1, /**< number of nonzeros in the constraint */
		&v_var, /**< 	array with variables of constraint entries */
		&coeff, /**< array with coefficients of constraint entries */
		-SCIPinfinity(scip),    /* lhs */
		0.0,                    /* rhs */
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
	v_conss.push_back(VCons(con, v_ind_));
	v_cons_map[v_ind_] = v_conss.size() - 1;
	SCIP_CALL(SCIPreleaseCons(scip, &con));

	/* release variable */
	SCIP_CALL(SCIPreleaseVar(scip, &v_var));

	return SCIP_OKAY;
}

/** add unsplittable constraint */
SCIP_RETCODE ProbDataNTRT::addUPCons(
	SCIP * 	scip, /**< SCIP data structure */
	const int f_ind_, /**< flow index */
	const int s, /**< source index */
	const int t /**< target index */
) {
	// create vertex constraint
	SCIP_CONS* con;
	SCIP_CALL(SCIPcreateConsSetpart(
		scip, /**< SCIP data structure */
		&con, /**< 	pointer to hold the created constraint */
		"up constraint", /**< 	name of constraint */
		0, /**< number of nonzeros in the constraint */
		NULL, /**< 	array with variables of constraint entries */
		TRUE,                   /* initial */
		FALSE,                  /* separate */
		TRUE,                   /* enforce */
		TRUE,                   /* check */
		TRUE,                   /* propagate */
		FALSE,                  /* local */
		TRUE,                   /* modifiable */
		FALSE,                  /* dynamic */
		FALSE,                  /* removable */
		FALSE));               /* stickingatnode */

	SCIP_CALL(SCIPaddCons(scip, con));
	SCIP_CALL(SCIPcaptureCons(scip, con));
	up_conss.push_back(UPCons(con, f_ind_, s, t));
	SCIP_CALL(SCIPreleaseCons(scip, &con));
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
	vector<SCIP_VAR *> codvar_array;
	vector<SCIP_Real> coef_array;
	for (int i = 0; i < cod_vars.size(); i++) {
		// check coding varibale is in the clique
		if (binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e12_ind) &&
			binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e23_ind) &&
			binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e32_ind) &&
			binary_search(e_ind_array.begin(), e_ind_array.end(), cod_vars[i].e21_ind)) {
			// add
			codvar_array.push_back(cod_vars[i].cod_var);
			coef_array.push_back(-cod_vars[i].rate);
		}
	}
	SCIP_CONS* con;
	SCIP_CALL(SCIPcreateConsLinear(
		scip, /**< SCIP data structure */
		&con, /**< 	pointer to hold the created constraint */
		"clique constraint", /**< name of constraint */
		codvar_array.size(), /**< number of nonzeros in the constraint */
		codvar_array.data(), /**< array with variables of constraint entries */
		coef_array.data(), /**< array with coefficients of constraint entries */
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


/** add path variable and its data in scip and problem data */
SCIP_RETCODE ProbDataNTRT::addPathVar(
	SCIP * 	scip, /**< SCIP data structure */
	const int p_class_, /**< path class index */
	const SCIP_Real p_cost_, /**< path cost */
	const vector<int>& v_array_, /**< vertices in the path */
	const vector<int>& e_array_, /**< edges in the path */
	const SCIP_Bool is_pricing /**< indicate pricing variable or not*/
) {
	SCIP_VAR* p_var;
	char p_name[20] = "path_var";
	//static int ct = 0;
	//SCIPdebugMessage("ct:%d\n", ct++);
	/* add variable to SCIP and to the problem data */
	//assert(p_cost_ *  flows[p_class_] > 0);
	if (is_pricing) {
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&p_var, /**< 	pointer to variable object*/
			p_name, /**< name of variable, or NULL for automatic name creation*/
			0.0, /**<	lower bound of variable*/
			1.0, /**< 	upper bound of variable */
			p_cost_ *  flows[p_class_], /**<	objective function value */
			SCIP_VARTYPE_BINARY, /**< type of variable */
			FALSE, /**<	should var's column be present in the initial root LP?*/
			TRUE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));
		SCIP_CALL(SCIPaddPricedVar(scip, p_var, 1.0));
	}
	else {
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&p_var, /**< 	pointer to variable object*/
			p_name, /**< name of variable, or NULL for automatic name creation*/
			0.0, /**<	lower bound of variable*/
			1.0, /**< 	upper bound of variable */
			p_cost_ *  flows[p_class_], /**<	objective function value */
			SCIP_VARTYPE_BINARY, /**< type of variable */
			TRUE, /**<	should var's column be present in the initial root LP?*/
			TRUE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));
		SCIP_CALL(SCIPaddVar(scip, p_var));
	}
	SCIP_CALL(SCIPchgVarUbLazy(scip, p_var, 1.0));
	// capture the variable
	SCIP_CALL(SCIPcaptureVar(scip, p_var));
	p_vars[p_class_].push_back(PathVar(p_class_, p_cost_, v_array_, e_array_, p_var));

	SCIP_VARDATA * vardata;
	SCIP_CALL(vardataCreate(scip, &vardata, -1,  (void *) &p_vars[p_class_].back() ) );
	/* set callback functions */
  	SCIPvarSetData(p_var, vardata);
	if(!is_pricing){
		SCIPvarSetDelorigData(p_var, vardataDelOrig);
		SCIPvarSetTransData(p_var, vardataTrans);
    	SCIPvarSetDeltransData(p_var, vardataDelTrans);
	}
	else{
		SCIPvarSetDeltransData(p_var, vardataDelTrans);
	}

	//SCIPdebugMessage("path var added!\n");
	// add the path variable in the coding opportunity constraint
	for (int i = 0; i < e_array_.size() - 1; i++) {
		int start_ind = e_array_[i];
		int end_ind = e_array_[i + 1];
		auto it = co_cons_map.find(make_pair(start_ind, end_ind));
		if (it != co_cons_map.end()) {
			SCIP_CALL(SCIPaddCoefLinear(scip, co_conss[it->second].co_cons, p_var , -flows[p_class_]));
		}
	}
	//SCIPdebugMessage("path var added in cod constraint!\n");

	path_add++;
	for (int c_ind = 0; c_ind < conflict->nc; c_ind++) {
		SCIP_Real coef = 0;
		for (int e_ind: e_array_) {
			if(conflict->c_array[c_ind].e_set.find(e_ind) != conflict->c_array[c_ind].e_set.end()){
				coef += flows[p_class_] / problemgraph->e_array[e_ind].e_capacity;
			}
		}
		if(coef > 0){
			SCIP_CALL(SCIPaddCoefLinear(scip, clq_conss[c_ind].clq_cons, p_var, coef));
			cts[c_ind]++;			
		}
	}
	//SCIPdebugMessage("path var added in clqiue constraint!\n");

	// add the path variable in the unsplittable constraint
	SCIP_CALL(SCIPaddCoefSetppc(scip, up_conss[p_class_].up_cons, p_var));

	//SCIPdebugMessage("path var added in up constraint!\n");

	/* release variable for the reader. */
	SCIP_CALL(SCIPreleaseVar(scip, &p_var));
	//SCIPdebugMessage("path var released!\n");

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

/** return the number of path variables */
int ProbDataNTRT::getNumPathVars() {
	int num = 0;
	for (int i = 0; i < nf; i++) {
		num += p_vars[i].size();
	}	
	return num;
}

/** check the flow is fixed */
bool ProbDataNTRT::isFlowFixed(
	int f /**< flow index */
) {
	bool fixed_one = FALSE;
	SCIP_VAR * one_var;
	// check wether path varibale in the flow is fixed
	for (auto it = p_vars[f].begin(); it != p_vars[f].end(); it++) {
		if (SCIPvarGetLbLocal(it->p_var) > 0.5) {
			one_var = it->p_var;
			fixed_one = TRUE;
			break;
		}
	}
	#ifdef SCIP_DEBUG_
	if(fixed_one){
		for(auto it = p_vars[f].begin(); it != p_vars[f].end(); it++){
			if(it->p_var != one_var){
				assert(SCIPvarGetUbLocal(it->p_var) < 0.5);
			}
		}
	}
	#endif
	return fixed_one;
}