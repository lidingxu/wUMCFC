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

/**@file   ConshdlrNTRT.h
 * @brief  C++ constraint handler for edge forbidden
 * @author Liding XU
 */

 /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONSHDLRNTRT_H__
#define __CONSHDLRNTRT_H__

#include "objscip/objscip.h"
#include "ProbDataNTRT.h"

using namespace scip;

struct SCIP_ConsData
{
	int*  e_inds;  /**< ponters to forbidden edges*/
	int n_e; /**< the number of forbidden edges */
	int 			 f_id;                /**< flow id */
	SCIP_NODE* node; /**< the node in the B&B-tree at which the cons is sticking */
	bool propagated; /**< is constraint already propagated? */
	int npropagatedvars; /**< number of variables that existed, the last time, the related node was
		*   propagated, used to determine whether the constraint should be
		*   repropagated*/
	int npropagations;  /**< stores the number propagations runs of this constraint */
};


/** C++ constraint handler for TSP subtour elimination constraints */
class ConshdlrNTRT : public ObjConshdlr
{
public:
	/** default constructor */
	ConshdlrNTRT(
		SCIP* scip /**< SCIP data structure */
	)
		: ObjConshdlr(scip, /**< SCIP data structure */
			"edge forbidden", /**< 	name of constraint handler */
			"stores the local edge forbbiden branch decisions", /**< description of constraint handler*/
			0, /**< priority of the constraint handler for separation */
			0, /**< . priority of the constraint handler for constraint enforcing */
			9999999, /**< . priority of the constraint handler for checking infeasibility (and propagation) */
			-1,  /**< frequency for separating cuts; zero means to separate only in the root node */
			1,  /**< . frequency for propagating domains; zero means only preprocessing propagation */
			1, /**< . frequency for using all instead of only the useful constraints in separation,
                *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
			0,  /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
			FALSE, /**< should separation method be delayed, if other separators found cuts? */
			FALSE,   /**< should propagation method be delayed, if other propagators found reductions? */
			TRUE,  /**< should the constraint handler be skipped, if no constraints are available? */
			SCIP_PROPTIMING_BEFORELP, /**< positions in the node solving loop where propagation method of constraint handlers should be executed */
			SCIP_PRESOLTIMING_FAST  /**< timing mask of the constraint handler's presolving method */)
	{
	}

	/** feasibility check method of constraint handler for primal solutions */
	virtual SCIP_DECL_CONSCHECK(scip_check);

	/** constraint enforcing method of constraint handler for LP solutions */
	virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

	/** constraint enforcing method of constraint handler for pseudo solutions */
	virtual SCIP_DECL_CONSENFOPS(scip_enfops);

	/** variable rounding lock method of constraint handler */
	virtual SCIP_DECL_CONSLOCK(scip_lock);

	/** frees specific constraint data
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 */
	virtual SCIP_DECL_CONSDELETE(scip_delete);

	/** transforms constraint data into data belonging to the transformed problem */
	virtual SCIP_DECL_CONSTRANS(scip_trans);


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
	virtual SCIP_DECL_CONSPROP(scip_prop);

	/** constraint activation notification method of constraint handler
	 *
	 *  @see SCIP_DECL_CONSACTIVE(x) in @ref type_cons.h
	 */
	virtual SCIP_DECL_CONSACTIVE(scip_active);

	/** constraint deactivation notification method of constraint handler
	 *
	 *  @see SCIP_DECL_CONSDEACTIVE(x) in @ref type_cons.h
	 */
	virtual SCIP_DECL_CONSDEACTIVE(scip_deactive);
	
}; /*lint !e1712*/

/** creates and captures an edge forbidden constraint */
SCIP_RETCODE SCIPcreateConsForbiddenEdge(
	SCIP*                 scip,               /**< SCIP data structure */
	SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
	const char*           name,               /**< name of constraint */
	const vector<int> & e_inds,           /**< forbidden edge index */
	const int 			  f_id,                /**< flow id */
	SCIP_NODE*            node,               /**< the node in the B&B-tree at which the cons is sticking */
	SCIP_Bool             local               /**< is constraint only valid locally? */
);

/** return the forbidden edge arrray */
vector<int> getForbiddenEdges(
	SCIP*                 scip,              /**< SCIP data structure */
	int f_id /**< flow id*/
);

#endif
