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

/**@file   reader_sub.cpp
 * @brief  D-optimal design problem reader
 * @author Liding XU
 * This file implements the reader/parser used to read the cbp input data. For more details see \ref NETWORKROUTING_READER.
 *
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

#include "probdata.h"
#include "reader_sub.h"

using namespace scip;
using namespace std;

/** destructor of file reader to free user data (called when SCIP is exiting) */
SCIP_DECL_READERFREE(ReaderSubmodular::scip_free)
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
SCIP_DECL_READERWRITE(ReaderSubmodular::scip_write)
{
	*result = SCIP_DIDNOTRUN;
	return SCIP_OKAY;
} /*lint !e715*/

/**@} */

/** problem reading method of reader
 *
 *  possible return values for *result:
 *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
 *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
 *
 *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
 */
SCIP_DECL_READERREAD(ReaderSubmodular::scip_read) {
	*result = SCIP_DIDNOTRUN;

   	SCIPdebugMessage("Start read!\n");
	ifstream filedata(filename);
	if (!filedata) {
		return SCIP_READERROR;
	}
	filedata.clear();

	// read parameters
    // SCIP_Real epsilon;
    SCIP_Real epsilon;
	int numvars, dim, card;
	filedata >> numvars >> dim >> card >> epsilon;


	vector<vector<SCIP_Real>> A(dim, vector<SCIP_Real>(numvars, 0));

	// open the file


	// read matrix A
	for(int i = 0; i < numvars; i++){
		for(int j = 0; j < dim; j++){
			filedata >> A[j][i];
			//printf("%f ", A[i,j]);
		}
		//printf("\n");
	}


	epsilon = sqrt(epsilon);
	SCIPdebugMessage("numvars:%d dim:%d card:%d\n", numvars, dim, card, epsilon);
	// create the problem's data structure
	
	ProbData * problemdata = NULL;
	problemdata = new ProbData(numvars, dim, A, card, epsilon);
	assert(problemdata != NULL);
	SCIPdebugMessage("--problem data completed!\n");
	SCIP_CALL(SCIPcreateObjProb(scip, filename, problemdata, FALSE));


	SCIP_CALL(problemdata->createInitial(scip));
   
   	*result = SCIP_SUCCESS;

	SCIPdebugMessage("--reader read completed!\n");
	return SCIP_OKAY;
}


/**@} */
