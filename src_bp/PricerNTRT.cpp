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


/**@file PricerNTRT.h
 * @brief NTRT pricer plugin
 * @author Liding XU
 */
#define SCIP_DEBUG
#define CHECK_DUAL2
#include "PricerNTRT.h"
#include "ProbDataNTRT.h"
#include "ProblemGraph.h"
#include "ConshdlrNTRT.h"
#include <time.h>
#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include "check.h"
#include "scip/clock.h"
#include "scip/type_set.h"
#include "objscip/objscip.h"
#include "scip/struct_cons.h"
#include "scip/cons_linear.h"
#include "scip/var.h"
#include "scip/struct_scip.h"
#include "scip/cons_setppc.h"

#define getDual(isfarkas, scip, cons) (isfarkas ? -SCIPgetDualfarkasLinear(scip, cons): -SCIPgetDualsolLinear(scip, cons))
#define getDualSetppc(isfarkas, scip, cons) (isfarkas ? -SCIPgetDualfarkasSetppc(scip, cons): -SCIPgetDualsolSetppc(scip, cons))

using namespace std;
using namespace scip;




/** initialization method of variable pricer (called after problem was transformed)
 *
 *  Because SCIP transformes the original problem in preprocessing, we need to get the references to
 *  the variables and constraints in the transformed problem from the references in the original
 *  problem.
 */
SCIP_DECL_PRICERINIT(PricerNTRT::scip_init)
{
	ProbDataNTRT * probdata = NULL;
	probdata = dynamic_cast<ProbDataNTRT*>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);
	pricinggraph = NULL;
	pricinggraph = new PricingGraph(probdata->conflict, probdata->problemgraph);
	return SCIP_OKAY;
} /*lint !e715*/

/** perform pricing
 *
 *  @todo compute shortest length restricted tour w.r.t. duals
 */
SCIP_RETCODE PricerNTRT::pricing(
	SCIP*                 scip,               /**< SCIP data structure */
	bool                  isfarkas            /**< whether we perform Farkas pricing */
) const{
	ProbDataNTRT * probdata = NULL;
	probdata = dynamic_cast<ProbDataNTRT *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);
	// branch decisions forbid edges
	ProblemGraph * problemgraph = NULL;
	problemgraph = probdata->problemgraph;
	assert(problemgraph != NULL);
	Conflict * conflict = NULL;
	conflict = probdata->conflict;
	assert(conflict != NULL);

	if(isfarkas){
		pricinggraph->resetFarkas();
	}
	else{
		pricinggraph->resetReducedCost();
	}
		for (int i = 0; i < probdata->co_conss.size(); i++) {
			SCIP_Real dual =  getDual(isfarkas, scip, probdata->co_conss[i].co_cons);
			//assert(SCIPisDualfeasGE(scip, dual, 0));
			int e_ind = pricinggraph->getEIndbyEPair(probdata->co_conss[i].start_e, probdata->co_conss[i].end_e);
			assert(e_ind != -1);
			pricinggraph->e_array[e_ind].e_cost -= dual;
		}
		// set clique constraints' duals to pricing graph cost
		for (int i = 0; i < probdata->clq_conss.size(); i++) {
			SCIP_Real dual =  getDual(isfarkas, scip, probdata->clq_conss[i].clq_cons);
			//assert(SCIPisDualfeasGE(scip, dual, 0));
			vector<int> * e_inds = pricinggraph->getEIndsofClq(probdata->clq_conss[i].c_ind);
			for (int e_ind : *e_inds) {
				int head_e_ind = pricinggraph->e_array[e_ind].head_e_ind;
				pricinggraph->e_array[e_ind].e_cost += dual / problemgraph->e_array[head_e_ind].e_capacity;
			}
		}
		// set vertex constraints' duals to pricing graph cost
		for (int i = 0; i < probdata->v_conss.size(); i++) {
			SCIP_Real dual =  getDual(isfarkas, scip, probdata->v_conss[i].v_cons);
			//assert(SCIPisDualfeasGE(scip, dual, 0));
			int v_ind = probdata->v_conss[i].v_ind;
			vector<int> * e_inds = pricinggraph->getEIndsofV(probdata->v_conss[i].v_ind);
			for (int e_ind : *e_inds) {
				int head_e_ind = pricinggraph->e_array[e_ind].head_e_ind;
				pricinggraph->e_array[e_ind].e_cost += dual / (2 * problemgraph->e_array[head_e_ind].e_capacity);
			}
		}
		bool isfind = false;
		// get the dual of unsplitable constraints and compute the shortest path
		for (int f = 0; f < probdata->nf; f++) {
			// if flow is fixed, continue
			if (probdata->isFlowFixed(f)) {
				continue;
			}
			// reset the forbidden edges
			vector<int> forb_e_array = getForbiddenEdges(scip, f);
			pricinggraph->resetNonForbidden();
			pricinggraph->setForbidden(forb_e_array);
			// get the dual
			SCIP_Real dual = getDualSetppc(isfarkas, scip, probdata->up_conss[f].up_cons);
			vector<int> path_e_array;
			vector<int> path_e_array_;
			path_e_array.reserve(problemgraph->nv);
			//SCIPdebugMessage("two way dijkstra\n");
			SCIP_Real len = pricinggraph->TWDijkstra(scip, probdata->st_pairs[f].first, probdata->st_pairs[f].second, path_e_array);
			//SCIP_Real len = pricinggraph->SPFA_PSEUDO(scip, probdata->st_pairs[f].first, probdata->st_pairs[f].second, path_e_array);
			if(path_e_array.empty()){
				continue;
			}
			if (SCIPisDualfeasNegative(scip, dual + len * probdata->flows[f])) {
				SCIP_Bool is_same = FALSE;
				vector<int> path_v_array;
				path_v_array.reserve(problemgraph->nv);
				// get the edge array and compute the path cost
				SCIP_Real p_cost = 0;
				for (int j = 0; j < path_e_array.size(); j++) {
					Edge_Prob & e = problemgraph->e_array[path_e_array[j]];
					p_cost += e.e_cost;
					path_v_array.push_back(e.e_tail);
				}
				path_v_array.push_back(problemgraph->e_array[path_e_array.back()].e_head);
				// add the path variable
				SCIP_CALL(probdata->addPathVar(scip, f, p_cost, path_v_array, path_e_array, TRUE));
				//SCIPdebugMessage("var added!\n");
				isfind = true;
			}
		}
		return SCIP_OKAY;
}



/** Pricing of additional variables if LP is feasible.
 *
 *  - get the values of the dual variables you need
 *  - construct the reduced-cost arc lengths from these values
 *  - find the shortest admissible path with respect to these lengths
 *  - if this tour has negative reduced cost, add it to the LP
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : at least one improving variable was found, or it is ensured that no such variable exists
 *  - SCIP_DIDNOTRUN  : the pricing process was aborted by the pricer, there is no guarantee that the current LP solution is optimal
 */
SCIP_DECL_PRICERREDCOST(PricerNTRT::scip_redcost)
{

	/* set result pointer, see above */
	(*result) = SCIP_DIDNOTRUN;

	pricing(scip, FALSE);
	(*result) = SCIP_SUCCESS;
	return SCIP_OKAY;
} /*lint !e715*/


/** Pricing of additional variables if LP is infeasible.
 *
 *  - get the values of the dual Farks multipliers you need
 *  - construct the reduced-cost arc lengths from these values
 *  - find the shortest admissible path with respect to these lengths
 *  - if this tour has negative reduced cost, add it to the LP
 */
SCIP_DECL_PRICERFARKAS(PricerNTRT::scip_farkas)
{
	pricing(scip, TRUE);
	return SCIP_OKAY;
} /*lint !e715*/





  