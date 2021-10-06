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

/**@file   ConshdlrNTRT.cpp
 * @brief  C++ constraint handler for edge forbidden
 * @author Liding XU
 */

 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "objscip/objscip.h"
#include "ProbDataNTRT.h"
#include "ConshdlrNTRT.h"
#include "scip/cons_linear.h"
#include <list>

using namespace scip;
using namespace std;

/** create constraint data */
static
SCIP_RETCODE consdataCreate(
	SCIP*                 scip,               /**< SCIP data structure */
	SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
	const vector<int> &  e_inds_,            /**< forbidden edge index */
	const int 			 f_id_,                /**< flow id */
	SCIP_NODE*            node                /**< the node in the B&B-tree at which the cons is sticking */
)
{
	assert(scip != NULL);
	assert(consdata != NULL);

	int n_e = e_inds_.size();

	int * e_inds = new int[n_e];
	for(int i = 0; i < n_e; i++){
		e_inds[i] = e_inds_[i];
	}	
	SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

	(*consdata)->e_inds = e_inds;
	(*consdata)->n_e = n_e;
	(*consdata)->f_id = f_id_;
	(*consdata)->npropagatedvars = 0;
	(*consdata)->npropagations = 0;
	(*consdata)->propagated = FALSE;
	(*consdata)->node = node;

	return SCIP_OKAY;
}

/** fixes a variable to zero if the corresponding packings are not valid for this constraint/node (due to branching) */
static
SCIP_RETCODE checkVariable(
	SCIP*                 scip,               /**< SCIP data structure */
	PathVar*              pathvar,                /**< path variable to check  */
	int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
	SCIP_Bool*            cutoff,              /**< pointer to store if a cutoff was detected */
	const int    e_ind /**< foribdden edge index */
)
{
	SCIP_Bool fixed;
	SCIP_Bool infeasible;

	assert(scip != NULL);
	assert(pathvar != NULL);
	assert(nfixedvars != NULL);
	assert(cutoff != NULL);

	/* if variables is locally fixed to zero continue */
	if (SCIPvarGetUbLocal(pathvar->p_var) < 0.5)
		return SCIP_OKAY;

	/* if the path does not contain this edge continue */
	bool contain = FALSE;
	for (int ind : pathvar->e_array) {
		if (ind == e_ind) {
			contain = TRUE;
			break;
		}
	}
	if (!contain) {
		return SCIP_OKAY;
	}

	SCIP_CALL(SCIPfixVar(scip, pathvar->p_var, 0.0, &infeasible, &fixed));
	if (infeasible)
	{
		assert(SCIPvarGetLbLocal(pathvar->p_var) > 0.5);
		(*cutoff) = TRUE;
	}
	else
	{
		assert(fixed);
		(*nfixedvars)++;
	}

	return SCIP_OKAY;
}

/** fixes variables to zero if the corresponding packings are not valid for this sonstraint/node (due to branching) */
static
SCIP_RETCODE consdataFixVariables(
	SCIP*                 scip,               /**< SCIP data structure */
	list<PathVar>&  p_vars,           /**< generated variables */
	const int e_ind,    /**< forbidden edges*/
	SCIP_RESULT*          result              /**< pointer to store the result of the fixing */
)
{
	int nfixedvars;
	int v;
	SCIP_Bool cutoff;

	nfixedvars = 0;
	cutoff = FALSE;

	for (auto it = p_vars.begin(); it != p_vars.end() && !cutoff ; it++) {
		SCIP_CALL(checkVariable(scip, &(*it), &nfixedvars, &cutoff, e_ind));
	}

	if (cutoff)
		*result = SCIP_CUTOFF;
	else if (nfixedvars > 0)
		*result = SCIP_REDUCEDDOM;

	return SCIP_OKAY;
}

/** frees specific constraint data
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 */
SCIP_DECL_CONSDELETE(ConshdlrNTRT::scip_delete)
{
   assert(consdata != NULL);
   assert(*consdata != NULL);
   delete[] (*consdata)->e_inds;
   SCIPfreeBlockMemory(scip, consdata);
	return SCIP_OKAY;
}/*lint !e715*/

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrNTRT::scip_trans) {
	//SCIPdebugMessage("trans conshdlr\n");
	SCIP_CONSDATA* sourcedata;
	SCIP_CONSDATA* targetdata;

	assert(sourcecons != NULL);
	assert(targetcons != NULL);
	assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
	sourcedata = SCIPconsGetData(sourcecons);
	assert(sourcedata != NULL);

	/* create constraint data for target constraint */
	vector<int> e_inds(sourcedata->e_inds, sourcedata->e_inds+sourcedata->n_e);
	SCIP_CALL(consdataCreate(scip, &targetdata, e_inds, sourcedata->f_id, sourcedata->node));

	/* create target constraint */
	SCIP_CALL(SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
		SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
		SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
		SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
		SCIPconsIsStickingAtNode(sourcecons)));

	return SCIP_OKAY;
}


/** domain propagation method of constraint handler
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched, but did not find any domain reductions
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
 */
SCIP_DECL_CONSPROP(ConshdlrNTRT::scip_prop) {
	assert(scip != NULL);
	assert(result != NULL);

	ProbDataNTRT * probdata = NULL;
	probdata = dynamic_cast<ProbDataNTRT *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);
	
	*result = SCIP_DIDNOTFIND;

	int numvars = probdata->getNumPathVars();
	for (int c_ind = 0; c_ind < nconss; c_ind++)
	{
		SCIP_CONSDATA* consdata = NULL;
		consdata = SCIPconsGetData(conss[c_ind]);
		assert(consdata != NULL);
		if (!consdata->propagated)
		{
			if(consdata->f_id == -1){
				for(int i = 0; i < probdata->nf; i++){
					for(int j = 0; j < consdata->n_e; j++){
						SCIP_CALL(consdataFixVariables(scip, probdata->p_vars[i], consdata->e_inds[j], result));
					}
				}
			}
			else{
				for(int j = 0; j < consdata->n_e; j++){
					SCIP_CALL(consdataFixVariables(scip, probdata->p_vars[consdata->f_id], consdata->e_inds[j], result));
				}
			}
			consdata->npropagations++;

			if (*result != SCIP_CUTOFF)
			{
				consdata->propagated = TRUE;
				consdata->npropagatedvars = numvars;
			}
			else {
				break;
			}
		}	
	}
	return SCIP_OKAY;
}

/** constraint activation notification method of constraint handler
 *
 *  @see SCIP_DECL_CONSACTIVE(x) in @ref type_cons.h
 */
SCIP_DECL_CONSACTIVE(ConshdlrNTRT::scip_active) {
	SCIP_CONSDATA* consdata = NULL;

	//SCIPdebugMessage("active debug\n");
	assert(scip != NULL);
	assert(cons != NULL);

	ProbDataNTRT * probdata = NULL;
	probdata = dynamic_cast<ProbDataNTRT *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);

	consdata = SCIPconsGetData(cons);
	assert(consdata != NULL);

	int num_path_vars = probdata->getNumPathVars();
	assert(consdata->npropagatedvars <= num_path_vars);

	if (consdata->npropagatedvars != num_path_vars)
	{
		consdata->propagated = FALSE;
		SCIP_CALL(SCIPrepropagateNode(scip, consdata->node));
	}

	return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler
 *
 *  @see SCIP_DECL_CONSDEACTIVE(x) in @ref type_cons.h
 */
SCIP_DECL_CONSDEACTIVE(ConshdlrNTRT::scip_deactive) {
	SCIP_CONSDATA* consdata;

	assert(scip != NULL);
	assert(cons != NULL);

	consdata = SCIPconsGetData(cons);
	assert(consdata != NULL);
	assert(consdata->propagated || SCIPgetNChildren(scip) == 0);

	ProbDataNTRT * probdata = NULL;
	probdata = dynamic_cast<ProbDataNTRT *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);

	/* set the number of propagated variables to current number of variables in SCIP */
	consdata->npropagatedvars = probdata->getNumPathVars();
	return SCIP_OKAY;
}

/** feasibility check method of constraint handler for primal solutions */
SCIP_DECL_CONSCHECK(ConshdlrNTRT::scip_check) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(result != NULL);

	/* do nothing */
	*result = SCIP_FEASIBLE;

	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
SCIP_DECL_CONSENFOLP(ConshdlrNTRT::scip_enfolp) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(result != NULL);

	/* do nothing */
	*result = SCIP_FEASIBLE;

	return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
SCIP_DECL_CONSENFOPS(ConshdlrNTRT::scip_enfops) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(result != NULL);

	/* do nothing */
	*result = SCIP_FEASIBLE;

	return SCIP_OKAY;
}
/** variable rounding lock method of constraint handler */
SCIP_DECL_CONSLOCK(ConshdlrNTRT::scip_lock) {
	assert(scip != NULL);
	assert(conshdlr != NULL);
	assert(cons != NULL);

	SCIPdebugMessage("Locking method for store graph constraint: <%s>.\n", SCIPconsGetName(cons));

	return SCIP_OKAY;
}


/** creates and captures an edge forbidden constraint */
SCIP_RETCODE SCIPcreateConsForbiddenEdge(
	SCIP*                 scip,               /**< SCIP data structure */
	SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
	const char*           name,               /**< name of constraint */
	const vector<int> &   e_inds,            /**< forbidden edge index */
	const int 			  f_id,                /**< flow id */
	SCIP_NODE*            node,               /**< the node in the B&B-tree at which the cons is sticking */
	SCIP_Bool             local               /**< is constraint only valid locally? */
)
{
	SCIP_CONSHDLR* conshdlr = NULL;
	SCIP_CONSDATA* consdata = NULL;

	/* find the forbidden edge constraint handler */
	conshdlr = SCIPfindConshdlr(scip, "edge forbidden");
	if (conshdlr == NULL)
	{
		SCIPerrorMessage("subtour constraint handler not found\n");
		return SCIP_PLUGINNOTFOUND;
	}

	/* create constraint data */
	SCIP_CALL(consdataCreate(scip, &consdata, e_inds, f_id, node));

	/* create constraint */
	SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
		local, FALSE, FALSE, FALSE, TRUE));

	//SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
	//	local, modifiable, dynamic, removable, FALSE));

	return SCIP_OKAY;
}

/** return the forbidden edge arrray */
vector<int> getForbiddenEdges(
	SCIP*                 scip,               /**< SCIP data structure */
	const int f_id  /**< flow id */
) {
	SCIP_CONSHDLR* conshdlr = NULL;
	SCIP_CONS** conss = NULL;
	SCIP_CONS* cons = NULL;
	int nconss;
	/* find the forbidden edge constraint handler */
	conshdlr = SCIPfindConshdlr(scip, "edge forbidden");
	assert(scip != NULL);
	assert(conshdlr != NULL);

	/* collect all branching decision constraints */
	conss = SCIPconshdlrGetConss(conshdlr);
	nconss = SCIPconshdlrGetNConss(conshdlr);
	vector<int> forb_e_array;
	for (int i = 0; i < nconss; i++) {
		if(!SCIPconsIsActive(conss[i])){
			continue;
		}
		SCIP_CONSDATA* consdata = NULL;
		consdata = SCIPconsGetData(conss[i]);
		assert(consdata != NULL);
		if(consdata->f_id == f_id || consdata->f_id == -1){
			for(int j = 0 ; j < consdata->n_e; j++){
				forb_e_array.push_back(consdata->e_inds[j]);
			}
		}
	}
	return forb_e_array;
}
