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

/**@file   BranchNTRT.cpp
 * @brief  event handler for variable deletion in NTRT
 * @author Liding XU
 */

 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include "objscip/objscip.h"
#include "ProblemGraph.h"
#include "BranchNTRT.h"
#include "ProbDataNTRT.h"
#include "ConshdlrNTRT.h"
#include <algorithm>
#include <unordered_set>
#include <utility>


bool myfunction (const pair<int, SCIP_Real> & i , const pair<int, SCIP_Real> & j) { return (i.second >j.second); }



void find_branch_edges(SCIP * scip,vector<int> & e1_inds, vector<int> & e2_inds, vector<PathVar *>& p_lst, vector<SCIP_Real> & sols, int f_id, int s, int t, ProblemGraph * probgraph){
	int trys = s;
	SCIP_Bool all_include = TRUE;
	int iter = 0;
	while(all_include){
		// test whether all paths include the trys node
		for(int i =0; i < p_lst.size(); i++){
			if(p_lst[i]->v_array[iter] != trys){
				all_include = FALSE;
				break;
			}
		}
		// if all include updated trys node
		if(all_include){
			iter++;
			s = trys;
			trys = p_lst[0]->v_array[iter];
		}
	}
	iter--;
	e1_inds.clear();
	e2_inds.clear();
	SCIP_Real max_flow = 0;
	int e1_ind = -1, e2_ind = -1;
	int v = s;
	vector<int> out_es;
	vector<int> fb_es = getForbiddenEdges(scip, f_id);
	sort(fb_es.begin(), fb_es.end());
	//printf("%ld", fb_es.size());
	for(int e_ind: *probgraph->getOutEdges(v)){
		if(!binary_search(fb_es.begin(), fb_es.end(), e_ind)){
			out_es.push_back(e_ind);
		}
	}
	out_es = *probgraph->getOutEdges(v);
	//assert(out_es.size() >= 2);
	vector<SCIP_Real> e_fracs(out_es.size(), 0);
	for(int j = 0; j < out_es.size(); j++){	
		for(int i = 0; i < p_lst.size(); i++){
			if(p_lst[i]->e_array[iter] == out_es[j]){
				e_fracs[j] += sols[i];
				break;
			}
		}
	}
	for(int i = 0; i < out_es.size(); i++){
		if(e_fracs[i] > max_flow){
			max_flow = e_fracs[i];
			e1_ind = out_es[i];
		}
	}
	max_flow = 0;
	for(int i = 0; i < out_es.size(); i++){
		if(e_fracs[i] > max_flow && out_es[i] != e1_ind){
			max_flow = e_fracs[i];
			e2_ind = out_es[i];
		}
	}
	for(int e_ind: out_es){
		if(e_ind == e1_ind || e_ind == e2_ind){
			continue;
		}
		else if(e1_inds.size() < e2_inds.size()){
			e1_inds.push_back(e_ind);
		}
		else{
			e2_inds.push_back(e_ind);
		}
	}
	//assert(e1_ind != -1);
	//assert(e2_ind != -1);
	e1_inds.push_back(e1_ind);
	e2_inds.push_back(e2_ind);
}


/** branching execution method for fractional LP solutions */
SCIP_DECL_BRANCHEXECLP(BranchNTRT::scip_execlp) {
	//SCIPdebugMessage("start branch...\n-----------------\n-------------\n----------------\n");
	SCIP_VAR** lpcands;
	SCIP_Real* lpcandsfrac;
	SCIP_Real* lpcandssol;
	int nlpcands;
	assert(scip != NULL);
	ProbDataNTRT * probdata = NULL;
	*result = SCIP_DIDNOTRUN;
	// get and convert the problem data
	probdata = dynamic_cast<ProbDataNTRT *>(SCIPgetObjProbData(scip));
	assert(probdata != NULL);
	//SCIPdebugMessage("add_t:%lf look_t:%lf\n", probdata->t_var, probdata->t_add);
	/* get fractional LP candidates */
	SCIP_CALL(SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, NULL, &nlpcands, NULL));
	//SCIPdebugMessage("\n");
	assert(nlpcands > 0);
	//SCIPdebugMessage("branch cands gotten...the number of Branch candidate %d %lf\n", nlpcands, SCIPgetDualbound(scip));
	// collect fractional path variables
	vector<VertexVar *> frac_vvars;
	vector<SCIP_Real> v_fracs;
	//vector<int> num_fracs(probdata->nf, 0);
	vector<vector<PathVar *>> frac_pvars(probdata->nf);
	vector<vector<SCIP_Real>> p_sols(probdata->nf);
	vector<SCIP_Bool> frac_flows(probdata->nf, FALSE);
	SCIP_Bool exist_p_frac = FALSE;	
	for (int i = 0; i < nlpcands; i++) {
		SCIP_VAR * frac_var = lpcands[i];
		SCIP_Real sol  = lpcandssol[i]; 
		SCIP_VARDATA * vardata = SCIPvarGetData(frac_var);
		if(vardata->id == -1){
			PathVar * pvar = (PathVar *) vardata->var;
			int f_id = pvar->p_class;
			frac_pvars[f_id].push_back(pvar);
			p_sols[f_id].push_back(sol);
			//num_fracs[f_id] += 1;
			frac_flows[f_id] = TRUE;
			exist_p_frac = TRUE;
		}
		else if(!exist_p_frac){
			frac_vvars.push_back(&probdata->v_vars[vardata->id]);
			assert(frac_var  == probdata->v_vars[vardata->id].v_var);
			v_fracs.push_back(sol);
			//printf("v\n");
		}
	}
	ProblemGraph * probgraph = probdata->problemgraph;
	// the branch rules only deals with path variables
	if(!exist_p_frac){
	    assert(v_fracs.size() > 0);
		//SCIPdebugMessage("----FRACTIONAL VERTEX LEFT!!!\n");
		int max_id = 0;
		SCIP_Real max_frac = 0;
		//printf("%ld\n", v_fracs.size());
		for(int i = 0; i < v_fracs.size(); i++){
			if(max_frac < v_fracs[i]){
				max_frac = v_fracs[i];
				max_id = i;
			}
		}
		int v_ind = frac_vvars[max_id]->v_ind;
		vector<int> outes = *probgraph->getOutEdges(v_ind);
		SCIP_NODE* downchild;
		SCIP_NODE* upchild;
		SCIP_CONS* conse;
		SCIP_CALL(SCIPcreateChild(scip, &downchild, 0.0, SCIPgetLocalTransEstimate(scip)));
		SCIP_CALL(SCIPcreateChild(scip, &upchild, 0.0, SCIPgetLocalTransEstimate(scip)));
		SCIP_CALL(SCIPchgVarUbNode(scip, downchild, frac_vvars[max_id]->v_var, 0.0));
		SCIP_CALL(SCIPchgVarLbNode(scip, upchild, frac_vvars[max_id]->v_var, 1.0));
		/* create corresponding constraints */
		SCIP_CALL(SCIPcreateConsForbiddenEdge(scip, &conse, "edge by path var", outes, -1, downchild, TRUE));

		/* add constraints to nodes */
		SCIP_CALL(SCIPaddConsNode(scip, downchild, conse, NULL));

		/* release constraints */
		SCIP_CALL(SCIPreleaseCons(scip, &conse));	
	}
	else{
		//SCIPdebugMessage("----FRACTIONAL PATH EXISTS!!!\n");
		// find the flow with the max number of paths
		int max_id = 0;
		SCIP_Real max_flow = 0;
		//int max_num = 0;
		for(int i = 0; i < probdata->nf; i++){
			if(frac_flows[i] && probdata->flows[i] > max_flow){
				max_flow = probdata->flows[i];
				max_id = i;
			}
		}
		// find the source and target of the flow
		int s = probdata->st_pairs[max_id].first;
		int t = probdata->st_pairs[max_id].second;
		
		// the partition set for branching
		vector<int> e1_inds;
		vector<int> e2_inds;

		find_branch_edges(scip, e1_inds, e2_inds, frac_pvars[max_id], p_sols[max_id], max_id, s, t, probdata->problemgraph);

		SCIP_NODE* childe1;
		SCIP_NODE* childe2;
		SCIP_CONS* conse1;
		SCIP_CONS* conse2;

		/* create the branch-and-bound tree child nodes of the current node */
		SCIP_CALL(SCIPcreateChild(scip, &childe1, 0.0, SCIPgetLocalTransEstimate(scip)));
		SCIP_CALL(SCIPcreateChild(scip, &childe2, 0.0, SCIPgetLocalTransEstimate(scip)));

		/* create corresponding constraints */
		SCIP_CALL(SCIPcreateConsForbiddenEdge(scip, &conse1, "edge by path var", e1_inds, max_id, childe1, TRUE));
		SCIP_CALL(SCIPcreateConsForbiddenEdge(scip, &conse2, "edge by path var", e2_inds, max_id, childe2, TRUE));

		/* add constraints to nodes */
		SCIP_CALL(SCIPaddConsNode(scip, childe1, conse1, NULL));
		SCIP_CALL(SCIPaddConsNode(scip, childe2, conse2, NULL));

		/* release constraints */
		SCIP_CALL(SCIPreleaseCons(scip, &conse1));
		SCIP_CALL(SCIPreleaseCons(scip, &conse2));		
	}

	*result = SCIP_BRANCHED;

	//SCIPdebugMessage("----FRACTIONAL PATH BRANCHED!!!\n");
	return SCIP_OKAY;
}


