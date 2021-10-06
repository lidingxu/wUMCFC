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
#include "ProbDataNTRT.h"
#include "BranchNTRT.h"

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
	ProblemGraph * probgraph = probdata->problemgraph;
	/* get fractional LP candidates */
	SCIP_CALL(SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, NULL, &nlpcands, NULL));
	//SCIPdebugMessage("\n");
	assert(nlpcands > 0);
	//SCIPdebugMessage("branch cands gotten...the number of Branch candidate %d\n", nlpcands);
	// collect fractional path variables
	vector<VertexVar *> frac_vvars;
	vector<SCIP_Real> v_fracs;
	vector<vector<SCIP_Bool>> frac_evars;
	vector<vector<SCIP_Real>> frac_val_evars;
	//vector<vector<SCIP_Real>> frac_sols;
	vector<SCIP_Bool> frac_flows(probdata->nf, FALSE);
	for(int i = 0; i < probdata->nf; i++){
		frac_evars.push_back(vector<SCIP_Bool>(probgraph->ne, FALSE));
		frac_val_evars.push_back(vector<SCIP_Real>(probgraph->ne, 0));
		//frac_sols.push_back(vector<SCIP_Real>(probgraph->ne, 0));
	}
	SCIP_Bool exist_p_frac = FALSE;
	for(int j = 0; j < probdata->nf; j++){
		for (int k = 0; k < probgraph->ne; k++) {
			SCIP_Real sol = SCIPgetVarSol(scip, probdata->e_vars[j][k].e_var);
			if(SCIPisPositive(scip, sol) && SCIPisLT(scip, sol, 1)){
				frac_evars[j][k] = TRUE;
				frac_val_evars[j][k] = sol;
				frac_flows[j] = TRUE;
				exist_p_frac = TRUE;
			}
		}
	}
	for(int j = 0; j < probdata->v_vars.size(); j++){
		SCIP_Real sol = SCIPgetVarSol(scip, probdata->v_vars[j].v_var);
		if(SCIPisPositive(scip, sol) && SCIPisLT(scip, sol, 1)){
			frac_vvars.push_back(&probdata->v_vars[j]);
			v_fracs.push_back(sol);
		}
	}
	// the branch rules only deals with path variables
	if(!exist_p_frac){
		//SCIPdebugMessage("----FRACTIONAL VERTEX LEFT!!!\n");
		int max_id = 0;
		SCIP_Real max_frac = 0;
		for(int i = 0; i < v_fracs.size(); i++){
			if(max_frac < v_fracs[i]){
				max_frac = v_fracs[i];
				max_id = i;
			}
		}
		SCIP_NODE* downchild;
		SCIP_NODE* upchild;

		SCIP_CALL(SCIPbranchVar(scip, frac_vvars[max_id]->v_var, &downchild, NULL, &upchild)); 		
	}
	else{
		//SCIPdebugMessage("----FRACTIONAL PATH EXISTS!!!\n");
		// find the flow with the max number of paths
		int max_id = -1;
		SCIP_Real max_flow = -1;
		for(int i = 0; i < probdata->nf; i++){
			if(frac_flows[i] && probdata->flows[i] > max_flow){
				max_flow = probdata->flows[i];
				max_id = i;
			}
		}
		assert(max_id != -1);
		int v = probdata->st_pairs[max_id].first;
		vector<int> e1_inds, e2_inds;
		int ct = 0;
		while(true){
			int fracnum = 0;
			for(int e_ind: *probgraph->getOutEdges(v)){
				//SCIPdebugMessage("%d, %d,  %s:%lf, %lf\n", max_id, e_ind, frac_evars[max_id][e_ind]?"y":"n", frac_val_evars[max_id][e_ind], SCIPgetVarSol(scip, probdata->e_vars[max_id][e_ind].e_var));
				if(frac_evars[max_id][e_ind]){
					fracnum++;
				}
			}
			if(fracnum >= 2){
				vector<int> out_es;
				for(int e_ind: *probgraph->getOutEdges(v)){
					if(SCIPvarGetUbLocal(probdata->e_vars[max_id][e_ind].e_var) > 0.5 ){
						out_es.push_back(e_ind);
					}
				}
				SCIP_Real max_flow = 0;
				int e1_ind = -1, e2_ind = -1;
				for(int e_ind: out_es){
					if(frac_evars[max_id][e_ind] && frac_val_evars[max_id][e_ind] > max_flow){
						max_flow = frac_val_evars[max_id][e_ind];
						e1_ind = e_ind;
					}
				}
				//SCIPdebugMessage("%d, %d, fe: %lf\n",max_id, e1_ind, max_flow);
				max_flow = 0;
				for(int e_ind: out_es){
					if(frac_evars[max_id][e_ind] && frac_val_evars[max_id][e_ind] > max_flow && e_ind != e1_ind){
						max_flow = frac_val_evars[max_id][e_ind];
						e2_ind = e_ind;
					}
				}
				//SCIPdebugMessage("%d, %d, se: %lf\n", max_id, e2_ind, max_flow);
				int toe1s = (out_es.size() - 2) / 2;
				for(int e_ind: out_es){
					if(e_ind == e1_ind || e_ind == e2_ind){
						continue;
					}	
					else if(e1_inds.size() < toe1s){
						e1_inds.push_back(e_ind);
					}		
					else{
						e2_inds.push_back(e_ind);
					}
				}
				assert(e1_ind != e2_ind);
				e1_inds.push_back(e1_ind);
				e2_inds.push_back(e2_ind);
				break;
			}
			else{
				SCIP_Bool find = FALSE;
				for(int e_ind: *probgraph->getOutEdges(v)){
					if(SCIPisEQ(scip, SCIPgetVarSol(scip, probdata->e_vars[max_id][e_ind].e_var), 1)){
						v = probgraph->e_array[e_ind].e_head;
						find = TRUE;
						break;
					}
				}
				if(!find){
					for (int k = 0; k < probgraph->ne; k++) {
						SCIP_Real sol = SCIPgetVarSol(scip, probdata->e_vars[max_id][k].e_var);
						SCIPdebugMessage("(%lf,%lf, %s),", sol, frac_val_evars[max_id][k], frac_evars[max_id][k]? "TRUE":"FALSE");
					}
				}
				if(!find){
					int s =1;	
				}
				assert(find);
			}
		}

		SCIP_NODE* childe1;
		SCIP_NODE* childe2;
		/* create the branch-and-bound tree child nodes of the current node */
		SCIP_CALL(SCIPcreateChild(scip, &childe1, 0.0, SCIPgetLocalTransEstimate(scip)));
		SCIP_CALL(SCIPcreateChild(scip, &childe2, 0.0, SCIPgetLocalTransEstimate(scip)));


		for(int e_ind: e1_inds){
			SCIP_CALL( SCIPchgVarUbNode(scip, childe1, probdata->e_vars[max_id][e_ind].e_var, 0) );
		}	

		for(int e_ind: e2_inds){
			SCIP_CALL( SCIPchgVarUbNode(scip, childe2, probdata->e_vars[max_id][e_ind].e_var, 0) );
		}	
	}

	*result = SCIP_BRANCHED;

	//SCIPdebugMessage("----FRACTIONAL PATH BRANCHED!!!\n");
	return SCIP_OKAY;
}



