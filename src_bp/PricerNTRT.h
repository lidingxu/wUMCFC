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

#ifndef __PRICERNTRT_H__
#define __PRICERNTRT_H__

#include "objscip/objscip.h"
#include "scip/pub_var.h"
#include "ProbDataNTRT.h"
#include "scip/scip_numerics.h"
#include "PricingGraph.h"

#include <vector>
#include <list>

using namespace std;
using namespace scip;


/** pricer class */
class PricerNTRT : public ObjPricer
{
public:

	/** Constructs the pricer object with the data needed */
	PricerNTRT(
		SCIP*                               scip,       /**< SCIP pointer */
		const char* VRP_PRICER_NAME /** < Pricer name */
	) : ObjPricer(scip, VRP_PRICER_NAME, "Finds path with negative reduced cost.", 0, TRUE) 
	{}

	/** Destructs the pricer object. */
	virtual ~PricerNTRT() {
		if (pricinggraph != NULL) {
			delete pricinggraph;
		}
	};

	/** initialization method of variable pricer (called after problem was transformed) */
	virtual SCIP_DECL_PRICERINIT(scip_init);

	/** reduced cost pricing method of variable pricer for feasible LPs */
	virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

	/** farkas pricing method of variable pricer for infeasible LPs */
	virtual SCIP_DECL_PRICERFARKAS(scip_farkas);

	/** perform pricing */
	SCIP_RETCODE pricing(
		SCIP*              scip,               /**< SCIP data structure */
		bool               isfarkas            /**< whether we perform Farkas pricing */
	) const;

	PricingGraph * pricinggraph;

	
};

#endif
