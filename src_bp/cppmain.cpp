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
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   NetworkRouting/src/cmain.c
 * @brief  Main file for NetworkRouting pricing example
 * @author Liding XU
 *
 *  This the file contains the \ref main() main function of the projects. This includes all the default plugins of
 *  \SCIP and the ones which belong to that projects. After that is starts the interactive shell of \SCIP or processes
 *  the shell arguments if given.
 */

 /* include SCIP components */

#include <iostream>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "table_xyz.h"

#include "ReaderNTRT.h"
#include "ConshdlrNTRT.h"
#include "PricerNTRT.h"
#include "BranchNTRT.h"

using namespace scip;
using namespace std;

/** creates and runs a SCIP instance with default and NTRT plugins */
static
SCIP_RETCODE runSCIP(
	int                        argc,          /**< number of arguments from the shell */
	char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
)
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include reader for networkrouting instances */
   SCIP_CALL( SCIPincludeObjReader(scip, new ReaderNTRT(scip), TRUE) );

   /* include networkrouting constraint handler */
   SCIP_CALL( SCIPincludeObjConshdlr(scip, new ConshdlrNTRT(scip), TRUE) );

   /* include edge forbidden branch handler */
   SCIP_CALL(SCIPincludeObjBranchrule(scip, new BranchNTRT(scip), TRUE));

   SCIP_CALL(SCIPsetRealParam(scip, "limits/gap", 1e-4));
   SCIP_CALL(SCIPsetRealParam(scip, "limits/absgap", 1e-6));
   SCIP_CALL(SCIPsetIntParam(scip, "timing/clocktype", 1));

   /* add user-defined statistics table*/
   SCIP_CALL( SCIPincludeTableXyz(scip));
   /* include networkrouting pricer  */
   static const char* NTRT_PRICER_NAME = "NTRT_Pricer";
   SCIP_CALL( SCIPincludeObjPricer(scip, new PricerNTRT(scip, NTRT_PRICER_NAME), TRUE) );

   /* for column generation instances, disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );

   /* set the parameter of network coding, default: FALSE */
   SCIP_CALL( SCIPaddBoolParam(scip, "ntrt/netcoding","use network coding",  NULL, FALSE, TRUE,  NULL, NULL) );

   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /**********************************
    * Process command line arguments *
    **********************************/
   SCIPdebugMessage("All plugins added\n");
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/
   /*
	for(auto cod_var: probdata->cod_vars){
		SCIP_Real cost = max(probdata->problemgraph->getCost(cod_var.e21_ind), probdata->problemgraph->getCost(cod_var.e23_ind));
		cod_var.cod_var
	}*/
   printf("finished");
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
