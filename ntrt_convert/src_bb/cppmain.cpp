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
#define SCIP_DEBUG

#include <iostream>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "BranchNTRT.h"
#include "ReaderNTRT.h"

using namespace std;
using namespace scip;
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

   /* include edge forbidden branch handler */
   //SCIP_CALL(SCIPincludeObjBranchrule(scip, new BranchNTRT(scip), TRUE));

   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIP_CALL(SCIPsetRealParam(scip, "limits/time", 3600));
   //SCIP_CALL(SCIPsetIntParam(scip, "constraints/setppc/propfreq", -1));
   //SCIP_CALL(SCIPsetIntParam(scip, "constraints/linear/propfreq", -1));
   //SCIP_CALL(SCIPsetIntParam(scip, "constraints/components", -1));
   //SCIP_CALL(SCIPsetIntParam(scip, "constraints/linear/sepafreq", -1));
   //SCIP_CALL(SCIPsetIntParam(scip, "constraints/linear/tightenboundsfreq", -1));
   //SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/rangedrowpropagation", FALSE));
   //SCIP_CALL(SCIPsetBoolParam(scip, "constraints/linear/rangedrowartcons", FALSE));
   //SCIP_CALL(SCIPsetBoolParam(scip, "conflict/enable", FALSE));

   /**********************************
    * Process command line arguments *
    **********************************/
   SCIPdebugMessage("All plugins added\n");
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/

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
   else{
      return  0;
   }
}