/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   table_xyz.c
 * @ingroup OTHER_CFILES
 * @brief  xyz statistics table
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "table_xyz.h"
#include "objscip/objscip.h"
#include "ConshdlrNTRT.h"
#include "ProbDataNTRT.h"
#include "ConshdlrNTRT.h"
#include "ProblemGraph.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace std;

#define TABLE_NAME              "xyz"
#define TABLE_DESC              "xyz statistics table"
#define TABLE_POSITION          -1000        /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE    SCIP_STAGE_SOLVED    /**< output of the statistics table is only printed from this stage onwards */




/*
 * Data structures
 */

/* TODO: fill in the necessary statistics table data */

/** statistics table data */
struct SCIP_TableData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of statistics table
 */

/* TODO: Implement all necessary statistics table methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for statistics table plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_TABLECOPY(tableCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableCopyXyz NULL
#endif


/** destructor of statistics table to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_TABLEFREE(tableFreeXyz)
{  
   SCIP_TABLEDATA* tabledata;
   tabledata = SCIPtableGetData(table);
   assert(tabledata != NULL);
   SCIPfreeMemory(scip, &tabledata);
   SCIPtableSetData(table, NULL);
   return SCIP_OKAY;
}
#else
#define tableFreeXyz NULL
#endif


/** initialization method of statistics table (called after problem was transformed) */
#if 0
static
SCIP_DECL_TABLEINIT(tableInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableInitXyz NULL
#endif


/** deinitialization method of statistics table (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_TABLEEXIT(tableExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableExitXyz NULL
#endif


/** solving process initialization method of statistics table (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_TABLEINITSOL(tableInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz statistics table not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define tableInitsolXyz NULL
#endif


/** solving process deinitialization method of statistics table (called before branch and bound process data is freed) */
static
SCIP_DECL_TABLEEXITSOL(tableExitsolXyz)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputXyz)
{  /*lint --e{715}*/
   assert(scip != NULL);
   ProbDataNTRT* probdata = NULL;
	probdata = dynamic_cast<ProbDataNTRT*>(SCIPgetObjProbData(scip));
   SCIP_SOL * sol =  SCIPgetBestSol(scip);
   SCIP_Real dec_cost = 0;
   for(auto cod_var: probdata->cod_vars){
      dec_cost += SCIPgetSolVal(scip, sol, cod_var.cod_var)* cod_var.redcost;
   }
   SCIPinfoMessage(scip, file, "decreased cost: %lf\n", dec_cost);
   return SCIP_OKAY;
}





/*
 * statistics table specific interface methods
 */

/** creates the xyz statistics table and includes it in SCIP */
SCIP_RETCODE SCIPincludeTableXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;
   

   /* create xyz statistics table data */
   tabledata = NULL;
   /* TODO: (optional) create statistics table specific data here */

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME, TABLE_DESC, TRUE,
         tableCopyXyz, tableFreeXyz, tableInitXyz, tableExitXyz,
         tableInitsolXyz, tableExitsolXyz, tableOutputXyz,
         tabledata, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   /* add xyz statistics table parameters */
   /* TODO: (optional) add statistics table specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
