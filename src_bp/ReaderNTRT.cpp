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

/**@file   ReaderNtrt.cpp
 * @brief  NetworkRouting problem reader
 * @author Liding XU
 *
 * This file implements the reader/parser used to read the networkrouting input data. For more details see \ref NETWORKROUTING_READER.
 *
 * @page NETWORKROUTING_READER Parsing the input format and creating the problem
 *
 * In the <code>data</code> directory you find a few data files which contain each one networkrouting problem. They have
 * the following structure. In the first line you find four integer numbers. The first one gives you the number of 
 * vertexes, the second the number of edges, the third the number of flows conflict cliques, and the fourth the number of
 * conflict cliques respectively. The first part of lines each contains an integer and an float which specify the
 * cost of the vertex. The second part of lines each contains two integers and two floats which specify
 * the edge from the first vertex to the second one, the weight of this edge and the capacity of this edge. The third
 * part of lines, for every two lines, one line contains the size of a conflict clique edges and the other contains a list of vertice 
 * in the clique. The fourth part of lines each contains two integer and one float
 * which specify the demand from the source of first integer to the target of the second integer. 
 *
 * For parsing that data, we implemented a reader plugin for \SCIP. A reader has several callback methods and at least
 * one interface methods (the one including the reader into \SCIP). For our purpose we only implemented the \ref
 * READERREAD "READERREAD" callback and the interface method which adds the reader plugin to \SCIP.

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <utility>

#include "objscip/objscip.h"

#include "ReaderNTRT.h"
#include "ProblemGraph.h"
#include "ProbDataNTRT.h"

using namespace scip;
using namespace std;

/** destructor of file reader to free user data (called when SCIP is exiting) */
SCIP_DECL_READERFREE(ReaderNTRT::scip_free)
{
	return SCIP_OKAY;
} /*lint !e715*/


/** problem writing method of reader; NOTE: if the parameter "genericnames" is TRUE, then
 *  SCIP already set all variable and constraint names to generic names; therefore, this
 *  method should always use SCIPvarGetName() and SCIPconsGetName();
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the writing to the file stream, it should return
 *  with RETCODE SCIP_WRITEERROR.
 */
SCIP_DECL_READERWRITE(ReaderNTRT::scip_write)
{
	*result = SCIP_DIDNOTRUN;

	return SCIP_OKAY;
} /*lint !e715*/

/** problem reading method of reader
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
 */
SCIP_DECL_READERREAD(ReaderNTRT::scip_read) {
	*result = SCIP_DIDNOTRUN;

	ProblemGraph * problemgraph = NULL;
	Conflict * conflict = NULL;

    SCIPdebugMessage("Start read!\n");
	// open the file
	ifstream filedata(filename);
	if (!filedata) {
		return SCIP_READERROR;
	}
	filedata.clear();

	// read the number of vertice, edges, cliques and demands
	int nv, ne, nc, nf;
	filedata >> nv >> ne >> nf >> nc;

	SCIPdebugMessage("Start read! nv:%d, ne:%d, nc:%d, nf:%d \n", nv, ne, nc, nf);
	// create the graph
	problemgraph = new ProblemGraph(nv, ne);
	assert(problemgraph != NULL);
	// read the vertex cost
	int received_nv = 0;
	while (!filedata.eof()) {
		SCIP_Real cost;
		filedata >> cost;
		problemgraph->addVertex(received_nv, cost);
		received_nv++;
		if (received_nv == nv) {
			break;
		}
	}
	assert(received_nv == nv);

	SCIPdebugMessage("--vertex read completed!\n");

	// read the edge ends, cost and capacity
	int received_ne = 0;
	while (!filedata.eof()) {
		SCIP_Real cost, capacity;
		int tail, head;
		filedata >> tail >> head >> cost >> capacity;
		problemgraph->addEdge(received_ne, tail, head, cost, capacity);
		received_ne++;
		if (received_ne == ne) {
			break;
		}
	}
	assert(received_ne == ne);

	SCIPdebugMessage("--edge read completed!\n");

	// create the conflict
	conflict = new Conflict(problemgraph, nc);
	// read the clique
	int received_nc = 0;
	while (!filedata.eof()) {
		int clique_size;
		filedata >> clique_size;
		int received_ne = 0;
		vector<int> e_array;
		while (!filedata.eof()) {
			int e_ind;
			filedata >> e_ind;
			e_array.push_back(e_ind);
			received_ne++;
			if (received_ne == clique_size) {
				break;
			}
		}
		assert(received_ne == clique_size);
		// add clique
		conflict->addClique(received_nc, clique_size, e_array);
		received_nc++;
		if (received_nc == nc) {
			break;
		}
	}
	assert(received_nc == nc);

    SCIPdebugMessage("--clique read completed!\n");

	// create the flows and source and target pairs of demands
	vector<SCIP_Real> flows = vector<SCIP_Real>(nf);
	vector<pair<int, int>> st_pairs = vector<pair<int, int>>(nf);
	// read the flow
	int received_nf = 0;
	while (!filedata.eof()) {
		filedata >> st_pairs[received_nf].first >> st_pairs[received_nf].second >> flows[received_nf];
		received_nf++;
		if (received_nf == nf) {
			break;
		}
	}
	assert(received_nf == nf);

	SCIPdebugMessage("--flow read completed!\n");

	// create the problem's data structure
	ProbDataNTRT * problemdata = NULL;
	problemdata = new ProbDataNTRT(problemgraph, conflict, nf, flows, st_pairs);
	assert(problemdata != NULL);
	SCIPdebugMessage("--problem data completed!\n");
	// deletedobject
	SCIP_CALL(SCIPcreateObjProb(scip, filename, problemdata, TRUE));

	SCIPdebugMessage("objprob created and creating inital solutions!\n");

	// create initial columns
	SCIP_CALL(problemdata->createInitialColumns(scip));
 
	static const char* NTRT_PRICER_NAME = "NTRT_Pricer";
    SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, NTRT_PRICER_NAME)) );
	
	*result = SCIP_SUCCESS;

	SCIPdebugMessage("--reader read completed!\n");
	return SCIP_OKAY;
}

/**@name Interface methods
 *
 * @{
 */


/**@} */
