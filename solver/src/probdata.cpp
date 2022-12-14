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

/**@file   probdata.cpp
 * @brief  Problem data for D-optimal design problem
 * @author Liding XU
 * /
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <string.h>
#include <utility>  
#include <math.h> 

#include "probdata.h"
#include "objscip/objscip.h"
#include "scip/struct_cons.h"
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#include "scip/expr_sum.h"
#include "scip/expr_var.h"
#include "scip/expr_pow.h"
#include "scip/expr_log.h"
#include "scip/expr_product.h"
#include "scip/var.h"
#include "scip/scip.h"
#include <stdlib.h>

using namespace scip;
using namespace std;



/** ProbData destructor */
ProbData::~ProbData()
{

}


/** release scip reference in probelme data*/
SCIP_RETCODE ProbData::releaseAll(
	SCIP*              scip                /**< SCIP data structure */
) {
	// release
	for (int i = 0; i < bin_vars.size(); i++) {
		SCIP_CALL(SCIPreleaseVar(scip, &bin_vars[i]));
	}


	for(int i = 0; i < numvars; i++){
		for(int j = 0; j < dim; j++){
			SCIP_CALL(SCIPreleaseVar(scip, &Z[i][j]));
			SCIP_CALL(SCIPreleaseVar(scip, &t[i][j]));
		}
	}

	for(int j = 0; j < dim; j++){
		SCIP_CALL(SCIPreleaseVar(scip, &t[numvars][j]));
	}

	for(int i = 0; i < dim; i++){
		for(int j = 0; j < dim; j++){
			SCIP_CALL(SCIPreleaseVar(scip, &J[i][j]));
			SCIP_CALL(SCIPreleaseVar(scip, &epsZ2[i][j]));
			SCIP_CALL(SCIPreleaseVar(scip, &epsZ[i][j]));
		}
	}

	SCIP_CALL(SCIPreleaseVar(scip, &obj_var));

	for (int i = 0; i < conss.size(); i++){
		SCIP_CALL(SCIPreleaseCons(scip, &conss[i]));
	}

	//SCIPdebugMessage("freed %d %d %d\n ", int(sizevar), numvar, numcons);
	return SCIP_OKAY;
}


/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbData::scip_delorig(
   SCIP*              scip                /**< SCIP data structure */
   )
{
	SCIP_CALL(releaseAll(scip));
	return SCIP_OKAY;
}







/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed)
 *
 *  The user has two possibilities to implement this method:
 *   1. Return the pointer to the original problem data object (this) as pointer to the transformed problem data
 *      object. The user may modify some internal attributes, but he has to make sure, that these modifications are
 *      reversed in the scip_deltrans() method, such that the original problem data is restored. In this case,
 *      he should set *deleteobject to FALSE, because the problem data must not be destructed by SCIP after the
 *      solving process is terminated.
 *   2. Call the copy constructor of the problem data object and return the created copy as transformed problem
 *      data object. In this case, he probably wants to set *deleteobject to TRUE, thus letting SCIP call the
 *      destructor of the object if the transformed problem data is no longer needed.
 */
SCIP_RETCODE ProbData::scip_trans(
   SCIP*              scip,               /**< SCIP data structure */
   ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
   SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
   )
{  /*lint --e{715}*/
   	assert( objprobdata != NULL );
   	assert( deleteobject != NULL );

   	// create and cpature transformed path varibles 
   	SCIPdebugMessage("start transform !!!!!!!!!!!\n");

    // allocate memory for target prob data
	ProbData * transprobdata = new ProbData(numvars, dim, A, card, epsilon);
	transprobdata->fullvalue = fullvalue;
	transprobdata->emptyvalue = emptyvalue;
	transprobdata->card = card;
	transprobdata->knapweights = knapweights;
	transprobdata->is_nature  = is_nature;
	transprobdata->gradient_cut = gradient_cut;

	SCIP_VAR* var;
	for (int i = 0; i < bin_vars.size(); i++) {
		SCIP_CALL(SCIPtransformVar(scip, bin_vars[i], &var));
		transprobdata->bin_vars.push_back(var);
	}


	for(int i = 0; i < numvars; i++){
		transprobdata->Z.push_back(vector<SCIP_VAR*>());
		transprobdata->t.push_back(vector<SCIP_VAR*>());
		for(int j = 0; j < dim; j++){		
			SCIP_CALL(SCIPtransformVar(scip, Z[i][j], &var));
			transprobdata->Z[i].push_back(var);
			SCIP_CALL(SCIPtransformVar(scip, t[i][j], &var));
			transprobdata->t[i].push_back(var);
		}
	}

	transprobdata->t.push_back(vector<SCIP_VAR*>());
	for(int j = 0; j < dim; j++){		
		SCIP_CALL(SCIPtransformVar(scip, t[numvars][j], &var));
		transprobdata->t[numvars].push_back(var);
	}


	for(int i = 0; i < dim; i++){
		transprobdata->J.push_back(vector<SCIP_VAR*>());
		transprobdata->epsZ2.push_back(vector<SCIP_VAR*>());
		transprobdata->epsZ.push_back(vector<SCIP_VAR*>());
		for(int j = 0; j < dim; j++){
			SCIP_CALL(SCIPtransformVar(scip, J[i][j], &var));
			transprobdata->J[i].push_back(var);
			SCIP_CALL(SCIPtransformVar(scip, epsZ[i][j], &var));
			transprobdata->epsZ[i].push_back(var);
			SCIP_CALL(SCIPtransformVar(scip, epsZ2[i][j], &var));
			transprobdata->epsZ2[i].push_back(var);
		}
	}



	SCIP_CALL(SCIPtransformVar(scip, obj_var, &var));
	transprobdata->obj_var = var;

	//SCIPdebugMessage("transform  4\n");

	for (int i = 0; i < conss.size(); i++){
		//SCIPprintCons( scip, conss[i], NULL ); 	
		//SCIPdebugMessage("%d/%d\n", i, conss.size() );
		SCIP_CONS * cons;
		SCIP_CALL(SCIPtransformCons(scip, conss[i], &cons ));
		transprobdata->conss.push_back(cons);
	}

	

   SCIPdebugMessage("end transform \n");
   assert( transprobdata != NULL );
   *objprobdata = transprobdata;           
   
   *deleteobject = FALSE;

   return SCIP_OKAY;
}      

/** destructor of user problem data to free original user data (called when original problem is freed)
 *
 *  If the "deleteobject" flag in the SCIPcreateObjProb() method was set to TRUE, this method is not needed,
 *  because all the work to delete the user problem data can be done in the destructor of the user problem
 *  data object. If the "deleteobject" flag was set to FALSE, and the user problem data object stays alive
 *  after the SCIP problem is freed, this method should delete all the problem specific data that is no
 *  longer needed.
 */
SCIP_RETCODE ProbData::scip_deltrans(
   SCIP*              scip                /**< SCIP data structure */
   )
{ 
   SCIP_CALL(releaseAll(scip));
   return SCIP_OKAY;
}



/** create variables and initial  constraints */
SCIP_RETCODE ProbData::createInitial(
	SCIP*                 scip               /**< SCIP data structure */
) {   

	// add binary variables
	
	for(int i = 0; i < numvars; i++){
		SCIP_VAR* bin_var;
		string tmp = "b"+ std::to_string(i);
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&bin_var, /**< 	pointer to variable object*/
			tmp.c_str(), /**< name of variable, or NULL for automatic name creation*/
			0, /**<	lower bound of variable*/
			1, /**< 	upper bound of variable */
			0, /**<	objective function value */
			SCIP_VARTYPE_BINARY, /**< type of variable */
			TRUE, /**<	should var's column be present in the initial root LP?*/
			FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));
		SCIP_CALL(SCIPaddVar(scip, bin_var));
		SCIP_CALL(SCIPcaptureVar(scip, bin_var));
		bin_vars.push_back(bin_var);
		SCIP_CALL(SCIPreleaseVar(scip, &bin_var));
	}

	// add objective variable (in minimization form, max to min - )
	SCIP_VAR * obj_var_;
	SCIP_CALL(SCIPcreateVar(
		scip, /**<	SCIP data structure*/
		&obj_var_, /**< 	pointer to variable object*/
		"obj_var", /**< name of variable, or NULL for automatic name creation*/
		-SCIPinfinity(scip), /**<	lower bound of variable*/
		SCIPinfinity(scip), /**< 	upper bound of variable */
		-1, /**<	objective function value */
		SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
		TRUE, /**<	should var's column be present in the initial root LP?*/
		FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
		NULL, NULL, NULL, NULL, NULL
	));
	SCIP_CALL(SCIPaddVar(scip, obj_var_));
	obj_var = obj_var_;
	SCIP_CALL(SCIPcaptureVar(scip, obj_var_));	
	SCIP_CALL(SCIPreleaseVar(scip, &obj_var_));



	emptyvalue = 2 * log(epsilon);

	// build MISOCP formulation
	SCIP_CONS * cons;
	// create variables
	for(int i = 0; i < numvars; i++){
		// Z
		Z.push_back(vector<SCIP_VAR *> ());
		for(int j = 0; j < dim; j++){
			SCIP_VAR * zij;
			string str = "z"+ std::to_string(i) + "_" + std::to_string(j);
			SCIP_CALL(SCIPcreateVar(
				scip, /**<	SCIP data structure*/
				&zij, /**< 	pointer to variable object*/
				str.c_str(), /**< name of variable, or NULL for automatic name creation*/
				-SCIPinfinity(scip), /**<	lower bound of variable*/
				SCIPinfinity(scip), /**< 	upper bound of variable */
				0, /**<	objective function value */
				SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
				TRUE, /**<	should var's column be present in the initial root LP?*/
				FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
				NULL, NULL, NULL, NULL, NULL
			));		
			SCIP_CALL(SCIPaddVar(scip, zij));
			SCIP_CALL(SCIPcaptureVar(scip, zij));
			Z[i].push_back(zij);
			SCIP_CALL(SCIPreleaseVar(scip, &zij));	
		}

		// t
		t.push_back(vector<SCIP_VAR *> ());
		for(int j = 0; j < dim; j++){
			SCIP_VAR * tij;
			string str = "t"+ std::to_string(i) + "_" + std::to_string(j);
			SCIP_CALL(SCIPcreateVar(
				scip, /**<	SCIP data structure*/
				&tij, /**< 	pointer to variable object*/
				str.c_str(), /**< name of variable, or NULL for automatic name creation*/
				0, /**<	lower bound of variable*/
				SCIPinfinity(scip), /**< 	upper bound of variable */
				0, /**<	objective function value */
				SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
				TRUE, /**<	should var's column be present in the initial root LP?*/
				FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
				NULL, NULL, NULL, NULL, NULL
			));			
			SCIP_CALL(SCIPaddVar(scip, tij));
			SCIP_CALL(SCIPcaptureVar(scip, tij));
			t[i].push_back(tij);
			SCIP_CALL(SCIPreleaseVar(scip, &tij));
		}	

	}

	// epsZ, t
	for(int j1 = 0; j1 < dim; j1++){
		epsZ.push_back(vector<SCIP_VAR *> ());
		for(int j2 = 0; j2 < dim; j2++){
			SCIP_VAR * z;
			string str = "epsz"+ std::to_string(j1) + "_" + std::to_string(j2);
			SCIP_CALL(SCIPcreateVar(
				scip, /**<	SCIP data structure*/
				&z, /**< 	pointer to variable object*/
				str.c_str(), /**< name of variable, or NULL for automatic name creation*/
				-SCIPinfinity(scip), /**<	lower bound of variable*/
				SCIPinfinity(scip), /**< 	upper bound of variable */
				0, /**<	objective function value */
				SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
				TRUE, /**<	should var's column be present in the initial root LP?*/
				FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
				NULL, NULL, NULL, NULL, NULL
			));		
			SCIP_CALL(SCIPaddVar(scip, z));
			SCIP_CALL(SCIPcaptureVar(scip, z));
			epsZ[j1].push_back(z);
			SCIP_CALL(SCIPreleaseVar(scip, &z));
		}
	}

	t.push_back(vector<SCIP_VAR *> ());
	for(int j = 0; j < dim; j++){
		SCIP_VAR * tij;
		string str = "t"+ std::to_string(numvars) + "_" + std::to_string(j);
		SCIP_CALL(SCIPcreateVar(
			scip, /**<	SCIP data structure*/
			&tij, /**< 	pointer to variable object*/
			str.c_str(), /**< name of variable, or NULL for automatic name creation*/
			0, /**<	lower bound of variable*/
			SCIPinfinity(scip), /**< 	upper bound of variable */
			0, /**<	objective function value */
			SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
			TRUE, /**<	should var's column be present in the initial root LP?*/
			FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
			NULL, NULL, NULL, NULL, NULL
		));			
		SCIP_CALL(SCIPaddVar(scip, tij));
		SCIP_CALL(SCIPcaptureVar(scip, tij));
		t[numvars].push_back(tij);
		SCIP_CALL(SCIPreleaseVar(scip, &tij));
	}

	// J
	for(int i = 0; i < dim; i++){
		J.push_back(vector<SCIP_VAR *> ());
		for(int j = 0; j < dim; j++){
			SCIP_VAR * Jij;
			string str = "J"+ std::to_string(i) + "_" +std::to_string(j);
			SCIP_CALL(SCIPcreateVar(
				scip, /**<	SCIP data structure*/
				&Jij, /**< 	pointer to variable object*/
				str.c_str(), /**< name of variable, or NULL for automatic name creation*/
				j > i ? 0 : (j == i ? 0 : -SCIPinfinity(scip)), /**< lower bound of variable*/
				j > i ? 0 : SCIPinfinity(scip), /**< 	upper bound of variable */
				0, /**<	objective function value */
				SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
				TRUE, /**<	should var's column be present in the initial root LP?*/
				FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
				NULL, NULL, NULL, NULL, NULL
			));			
			SCIP_CALL(SCIPaddVar(scip, Jij));
			SCIP_CALL(SCIPcaptureVar(scip, Jij));
			J[i].push_back(Jij);
			SCIP_CALL(SCIPreleaseVar(scip, &Jij));
		}
	}

	
	// build constraint
	// \sum Ai Zi = J, J lower-trigangular
	for(int j1 = 0; j1 < dim; j1++){
		for(int j2 = j1; j2 < dim; j2++){
			vector<SCIP_VAR*> vars(numvars + 2);
			vector<SCIP_Real> weights(numvars + 2);
			for(int  i = 0; i < numvars; i++){
				SCIP_Real Aij1 = A[j1][i];
				SCIP_VAR* Zij2 =  Z[i][j2];
			    weights[i] = Aij1;
				vars[i] = Zij2;
			}	
			vars[numvars] = epsZ[j1][j2];
			weights[numvars] = epsilon;
			vars[numvars + 1] = J[j1][j2];
			weights[numvars + 1] = -1;
			string str = "A" + std::to_string(j1) + "Z" + std::to_string(j2) + "=J";
			SCIP_CALL(SCIPcreateConsLinear(
				scip,               /**< SCIP data structure */
				&cons,        /**< pointer to hold the created constraint */
				str.c_str(),             /**< name of constraint */
				numvars + 2,            /**< number of variables in the constraint */
				vars.data(),    /**< array with variables of constraint entries */
				weights.data(),
				0,
				0,             
				TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
															*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
				TRUE,                /**< should the constraint be separated during LP processing?
															*   Usually set to TRUE. */
				TRUE,               /**< should the constraint be enforced during node processing?
															*   TRUE for model constraints, FALSE for additional, redundant constraints. */
				TRUE,               /**< should the constraint be checked for feasibility?
															*   TRUE for model constraints, FALSE for additional, redundant constraints. */
				TRUE,               /**< should the constraint be propagated during node processing?
															*   Usually set to TRUE. */
				FALSE,
				FALSE,              /**< is constraint only valid locally?
															*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
				FALSE,              /**< is constraint subject to aging?
															*   Usually set to FALSE. Set to TRUE for own cuts which
															*   are separated as constraints. */
				FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
															*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
				FALSE               /**< should the constraint always be kept at the node where it was added, even
															*   if it may be moved to a more global node?
															*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
			));
			SCIP_CALL(SCIPaddCons(scip, cons));
			SCIP_CALL(SCIPcaptureCons(scip, cons));
			conss.push_back(cons);
			SCIP_CALL(SCIPreleaseCons(scip, &cons));	
		}
	}

	// \sum_{i} tij \le Jjj
	for(int j = 0; j < dim; j++){
		vector<SCIP_VAR *> vars(numvars + 2, NULL);
		vector<SCIP_Real> weights(numvars + 2, 1);
		for(int i = 0; i < numvars + 1; i++){
			vars[i] = t[i][j];
		}
		// test: weights[numvars] = 0;
		vars[numvars + 1] = J[j][j];
		weights[numvars + 1] = -1;
		string str = "sumt" + std::to_string(j) + "<=J" + std::to_string(j);
		SCIP_CALL(SCIPcreateConsLinear(
			scip,               /**< SCIP data structure */
			&cons,        /**< pointer to hold the created constraint */
			str.c_str(),             /**< name of constraint */
			numvars + 2,            /**< number of variables in the constraint */
			vars.data(),    /**< array with variables of constraint entries */
			weights.data(),
			-SCIPinfinity(scip),
			0,             
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
		));
		SCIP_CALL(SCIPaddCons(scip, cons));
		SCIP_CALL(SCIPcaptureCons(scip, cons));
		conss.push_back(cons);
		SCIP_CALL(SCIPreleaseCons(scip, &cons));
	}

	// zij^2 \le tij wij
	for(int i = 0; i < numvars; i++){
		for(int j = 0; j < dim; j++){
			string str = "soc" + std::to_string(i) + std::to_string(j);
			SCIP_VAR * vars1[2] = {Z[i][j], t[i][j]};
			SCIP_VAR * vars2[2] = {Z[i][j], bin_vars[i]};
			SCIP_Real coefs[2] = {1, -1}; 
			SCIP_CALL(SCIPcreateConsQuadraticNonlinear(
				scip,               	/**< SCIP data structure */
				&cons,       			/**< pointer to hold the created constraint */
				str.c_str(),            /**< name of constraint */
				0,   					/**< 	number of linear terms  */
				NULL,					/**<  array with variables in linear part */
				NULL,					/**< array with coefficients of variables in linear part  */
				2,						/**< number of quadratic terms */
				vars1, 	/**< array with first variables in quadratic terms */
				vars2, 	/**< array with second variables in quadratic terms */
				coefs,	/** array with coefficients of quadratic terms  */
				-SCIPinfinity(scip),
				0, 
				TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
															*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
				TRUE,                /**< should the constraint be separated during LP processing?
															*   Usually set to TRUE. */
				TRUE,               /**< should the constraint be enforced during node processing?
															*   TRUE for model constraints, FALSE for additional, redundant constraints. */
				TRUE,               /**< should the constraint be checked for feasibility?
															*   TRUE for model constraints, FALSE for additional, redundant constraints. */
				TRUE,               /**< should the constraint be propagated during node processing?
															*   Usually set to TRUE. */
				FALSE,              /**< is constraint only valid locally?
															*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
				FALSE,              /**< is constraint subject to aging?
															*   Usually set to FALSE. Set to TRUE for own cuts which
															*   are separated as constraints. */
				FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
															*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
				FALSE               /**< should the constraint always be kept at the node where it was added, even
															*   if it may be moved to a more global node?
															*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
			));	
			SCIP_CALL(SCIPaddCons(scip, cons));
			SCIP_CALL(SCIPcaptureCons(scip, cons));
			conss.push_back(cons);
			SCIP_CALL(SCIPreleaseCons(scip, &cons));

			// this is needed for SCIP 8.0.1
			// linearize +/-  2Z[i][j] <= t[i][j] + bin_vars[i] 
			SCIP_VAR * vars3[3] = {Z[i][j], t[i][j], bin_vars[i]};
			SCIP_Real wts3[3]  = {2, -1, -1};
			for(int k = 0; k < 2; k++){
				wts3[0] = k ? 2 : -2;
				SCIP_CALL(SCIPcreateConsLinear(
					scip,               /**< SCIP data structure */
					&cons,        /**< pointer to hold the created constraint */
					"linear1",             /**< name of constraint */
					3,            /**< number of variables in the constraint */
					vars3,    /**< array with variables of constraint entries */
					wts3,
					-SCIPinfinity(scip),
					0,             
					TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
																*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
					TRUE,                /**< should the constraint be separated during LP processing?
																*   Usually set to TRUE. */
					FALSE,               /**< should the constraint be enforced during node processing?
																*   TRUE for model constraints, FALSE for additional, redundant constraints. */
					FALSE,               /**< should the constraint be checked for feasibility?
																*   TRUE for model constraints, FALSE for additional, redundant constraints. */
					FALSE,               /**< should the constraint be propagated during node processing?
																*   Usually set to TRUE. */
					FALSE,              /**< is constraint only valid locally?
																*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
					FALSE,					/**< is constraint modifiable (subject to column generation)? Usually set to FALSE. In column generation applications, set to TRUE if pricing adds coefficients to this constraint. */
					FALSE,              /**< is constraint subject to aging?
																*   Usually set to FALSE. Set to TRUE for own cuts which
																*   are separated as constraints. */
					TRUE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
																*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
					FALSE               /**< should the constraint always be kept at the node where it was added, even
																*   if it may be moved to a more global node?
																*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
				));
				SCIP_CALL(SCIPaddCons(scip, cons));
				SCIP_CALL(SCIPcaptureCons(scip, cons));
				conss.push_back(cons);
				SCIP_CALL(SCIPreleaseCons(scip, &cons));
			}
		}
	}

	// soc for epsZ
	for(int j1 = 0; j1 < dim; j1++)
	{
		epsZ2.push_back(vector<SCIP_VAR *> ());
		for(int j2 = 0; j2 < dim; j2++){
			// epsZ[j1][j2]
			SCIP_VAR * z;
			string str = "eps2z"+ std::to_string(j1) + "_" + std::to_string(j2);
			SCIP_CALL(SCIPcreateVar(
				scip, /**<	SCIP data structure*/
				&z, /**< 	pointer to variable object*/
				str.c_str(), /**< name of variable, or NULL for automatic name creation*/
				-SCIPinfinity(scip), /**<	lower bound of variable*/
				SCIPinfinity(scip), /**< 	upper bound of variable */
				0, /**<	objective function value */
				SCIP_VARTYPE_CONTINUOUS, /**< type of variable */
				TRUE, /**<	should var's column be present in the initial root LP?*/
				FALSE, /**<	is var's column removable from the LP (due to aging or cleanup)?*/
				NULL, NULL, NULL, NULL, NULL
			));		
			SCIP_CALL(SCIPaddVar(scip, z));
			SCIP_CALL(SCIPcaptureVar(scip, z));
			epsZ2[j1].push_back(z);
			SCIP_CALL(SCIPreleaseVar(scip, &z));
		}
	}

	for(int j2 = 0; j2 < dim; j2++)
	{
		for(int j1 = 0; j1 < dim; j1++)
		{
			// epsZ[j1][j2]^2 <= epsZ2[j1][j2]
			SCIP_VAR* linearvars[1] = {epsZ2[j1][j2]};
			SCIP_Real linearcoefs[1] = {-1};
			SCIP_VAR* quadvars[1] = {epsZ[j1][j2]};
			SCIP_Real quadcoefs[1] = {1};
			string str = "epsZ^2<=epsZ2 "+ std::to_string(j1) + "_" +std::to_string(j2) ;
			SCIP_CALL(SCIPcreateConsQuadraticNonlinear(
				scip,               	/**< SCIP data structure */
				&cons,       			/**< pointer to hold the created constraint */
				str.c_str(),            /**< name of constraint */
				1,   					/**< 	number of linear terms  */
				linearvars,					/**<  array with variables in linear part */
				linearcoefs,					/**< array with coefficients of variables in linear part  */
				1,						/**< number of quadratic terms */
				quadvars, 	/**< array with first variables in quadratic terms */
				quadvars, 	/**< array with second variables in quadratic terms */
				quadcoefs,	/** array with coefficients of quadratic terms  */
				-SCIPinfinity(scip),
				0, 
				TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
															*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
				TRUE,                /**< should the constraint be separated during LP processing?
															*   Usually set to TRUE. */
				TRUE,               /**< should the constraint be enforced during node processing?
															*   TRUE for model constraints, FALSE for additional, redundant constraints. */
				TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
				TRUE,               /**< should the constraint be propagated during node processing?
															*   Usually set to TRUE. */
				FALSE,              /**< is constraint only valid locally?
															*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
				FALSE,              /**< is constraint subject to aging?
															*   Usually set to FALSE. Set to TRUE for own cuts which
															*   are separated as constraints. */
				FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
															*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
				FALSE               /**< should the constraint always be kept at the node where it was added, even
															*   if it may be moved to a more global node?
															*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
			));	
			SCIP_CALL(SCIPaddCons(scip, cons));
			SCIP_CALL(SCIPcaptureCons(scip, cons));
			conss.push_back(cons);
			SCIP_CALL(SCIPreleaseCons(scip, &cons));	

			// linearize epsZ^2 <= epsZ2
			str += "linearize";
			SCIP_VAR * vars3[2] = {epsZ[j1][j2], epsZ2[j1][j2]};
			SCIP_Real wts3[2]  = {2, 1};
			for(int k = 0; k < 2; k++){
				wts3[0] = k ? 2 : -2;
				SCIP_CALL(SCIPcreateConsLinear(
					scip,               /**< SCIP data structure */
					&cons,        /**< pointer to hold the created constraint */
					str.c_str(),             /**< name of constraint */
					2,            /**< number of variables in the constraint */
					vars3,    /**< array with variables of constraint entries */
					wts3,
					-1,
					SCIPinfinity(scip),             
					TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
																*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
					TRUE,                /**< should the constraint be separated during LP processing?
																*   Usually set to TRUE. */
					FALSE,               /**< should the constraint be enforced during node processing?
																*   TRUE for model constraints, FALSE for additional, redundant constraints. */
					FALSE,               /**< should the constraint be checked for feasibility?
																*   TRUE for model constraints, FALSE for additional, redundant constraints. */
					FALSE,               /**< should the constraint be propagated during node processing?
																*   Usually set to TRUE. */
					FALSE,              /**< is constraint only valid locally?
																*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
					FALSE,					/**< is constraint modifiable (subject to column generation)? Usually set to FALSE. In column generation applications, set to TRUE if pricing adds coefficients to this constraint. */
					FALSE,              /**< is constraint subject to aging?
																*   Usually set to FALSE. Set to TRUE for own cuts which
																*   are separated as constraints. */
					TRUE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
																*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
					FALSE               /**< should the constraint always be kept at the node where it was added, even
																*   if it may be moved to a more global node?
																*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
				));
				SCIP_CALL(SCIPaddCons(scip, cons));
				SCIP_CALL(SCIPcaptureCons(scip, cons));
				conss.push_back(cons);
				SCIP_CALL(SCIPreleaseCons(scip, &cons));
			}		
		}

		// soc 
		string str = "soc" + std::to_string(numvars) + std::to_string(j2);
		vector<SCIP_VAR *> vars(dim + 1);
		vector<SCIP_Real> coefs(dim + 1, 1);
		for(int j1 = 0; j1 < dim; j1++){
			vars[j1] = epsZ2[j1][j2];
		}
		vars[dim] = t[numvars][j2];
		coefs[dim] = -1;
		SCIP_CALL(SCIPcreateConsLinear(
			scip,               /**< SCIP data structure */
			&cons,        /**< pointer to hold the created constraint */
			str.c_str(),             /**< name of constraint */
			dim + 1,            /**< number of variables in the constraint */
			vars.data(),    /**< array with variables of constraint entries */
			coefs.data(),
			-SCIPinfinity(scip),
			0,             
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
		));
		SCIP_CALL(SCIPaddCons(scip, cons));
		SCIP_CALL(SCIPcaptureCons(scip, cons));
		conss.push_back(cons);
		SCIP_CALL(SCIPreleaseCons(scip, &cons));
	}

	// objective function form
	bool logdet_form = false;
	if(logdet_form)
	{
		vector<SCIP_EXPR *> children(dim);
		for(int  j = 0; j < dim; j++){
			SCIP_EXPR* varexpr;
			SCIP_EXPR* logexpr;
			SCIP_CALL(SCIPcreateExprVar(scip, &varexpr, J[j][j], NULL, NULL) );
			//SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, numvars, exprs.data(), coefs.data(), 0.0, NULL, NULL) );
			SCIP_CALL(SCIPcreateExprLog(scip, &logexpr, varexpr, NULL, NULL)); 
			children[j] = logexpr;
			SCIP_CALL(SCIPreleaseExpr(scip, &varexpr));
		}

		SCIP_EXPR * objexpr;
		SCIP_CALL(SCIPcreateExprVar(scip, &objexpr, obj_var, NULL, NULL) );
		children.push_back(objexpr);
		
		//SCIPdebugMessage("%d %f\n", dim, 1.0 / dim);
		vector<SCIP_Real> coeffs(dim + 1, 1.0 / dim);
		coeffs[dim] = -1;

	
		SCIP_EXPR * topexpr;
		SCIP_CALL( SCIPcreateExprSum(scip, &topexpr, dim + 1, children.data(), coeffs.data(), 0.0, NULL, NULL) );


		SCIP_CALL(SCIPcreateConsNonlinear(
			scip,               /**< SCIP data structure */
			&cons,               /**< pointer to hold the created constraint */
			"obj_cons",               /**< name of constraint */
			topexpr,               /**< expression of constraint (must not be NULL) */
			emptyvalue,                /**< left hand side of constraint */
			SCIPinfinity(scip),                /**< right hand side of constraint */
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
			));


		SCIP_CALL(SCIPaddCons(scip, cons));
		SCIP_CALL(SCIPcaptureCons(scip, cons));
		conss.push_back(cons);
		SCIP_CALL(SCIPreleaseCons(scip, &cons));
	
		SCIP_CALL(SCIPreleaseExpr(scip, &topexpr));

		for(int j = 0; j < dim + 1; j++){
			SCIP_CALL(SCIPreleaseExpr(scip, &children[j]));
		}	


    
		vector<SCIP_VAR*> vars(dim + 1, NULL);
		vector<SCIP_Real>  wts(dim + 1, 0);
		for(int i = 0; i < dim; i++){
			vars[i] = J[i][i];
			wts[i] = 1.0 / dim;
		}
		vars[dim] = obj_var;
		wts[dim] = -1;

		// manually linearize
		SCIP_CALL(SCIPcreateConsLinear(
			scip,               /**< SCIP data structure */
			&cons,        /**< pointer to hold the created constraint */
			"linearize_obj",             /**< name of constraint */
			dim+1,            /**< number of variables in the constraint */
			vars.data(),    /**< array with variables of constraint entries */
			wts.data(),
			1 + emptyvalue,
			SCIPinfinity(scip),           
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
		));
		// this may not needed for SCIP 8.0.1
		//SCIP_CALL(SCIPaddCons(scip, cons));
		//SCIP_CALL(SCIPcaptureCons(scip, cons));
		//conss.push_back(cons);
		//SCIP_CALL(SCIPreleaseCons(scip, &cons));		
	}
	else
	{
		vector<SCIP_EXPR *> children(dim);
		for(int  j = 0; j < dim; j++){
			SCIP_EXPR* varexpr;
			SCIP_EXPR* powexpr;
			SCIP_CALL(SCIPcreateExprVar(scip, &varexpr, J[j][j], NULL, NULL) );
			//SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, numvars, exprs.data(), coefs.data(), 0.0, NULL, NULL) );
			SCIP_CALL(SCIPcreateExprPow(scip, &powexpr, varexpr, 1.0 / dim, NULL, NULL)); 
			children[j] = powexpr;
			SCIP_CALL(SCIPreleaseExpr(scip, &varexpr));
		}

		SCIP_EXPR * prodexpr;
		SCIPcreateExprProduct(scip,  &prodexpr, dim, children.data(), 1, NULL, NULL);
		for(int j = 0; j < dim; j++){
			SCIP_CALL(SCIPreleaseExpr(scip, &children[j]));
		}

		SCIP_EXPR * objexpr;
		SCIP_CALL(SCIPcreateExprVar(scip, &objexpr, obj_var, NULL, NULL) );


		SCIP_EXPR * exprs[2] = {prodexpr, objexpr};
		SCIP_Real  coeffs[2] = {1., -1.};
		SCIP_EXPR * topexpr;
		SCIP_CALL( SCIPcreateExprSum(scip, &topexpr, 2, exprs, coeffs, 0.0, NULL, NULL) );


		SCIP_CALL(SCIPcreateConsNonlinear(
			scip,               /**< SCIP data structure */
			&cons,               /**< pointer to hold the created constraint */
			"obj_cons",               /**< name of constraint */
			topexpr,               /**< expression of constraint (must not be NULL) */
			0,                /**< left hand side of constraint */
			SCIPinfinity(scip),                /**< right hand side of constraint */
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
			));
		
		SCIP_CALL(SCIPreleaseExpr(scip, &objexpr));
		SCIP_CALL(SCIPreleaseExpr(scip, &prodexpr));
		SCIP_CALL(SCIPreleaseExpr(scip, &topexpr));

		SCIP_CALL(SCIPaddCons(scip, cons));
		SCIP_CALL(SCIPcaptureCons(scip, cons));
		conss.push_back(cons);
		SCIP_CALL(SCIPreleaseCons(scip, &cons));
		
		vector<SCIP_VAR*> vars(dim + 1, NULL);
		vector<SCIP_Real>  wts(dim + 1,0);
		for(int i = 0; i < dim; i++){
			vars[i] = J[i][i];
			wts[i] = 1.0 / dim;
		}
		vars[dim] = obj_var;
		wts[dim] = -1;

		// manually linearize

		SCIP_CALL(SCIPcreateConsLinear(
			scip,               /**< SCIP data structure */
			&cons,        /**< pointer to hold the created constraint */
			"linear_obj",             /**< name of constraint */
			dim+1,            /**< number of variables in the constraint */
			vars.data(),    /**< array with variables of constraint entries */
			wts.data(),
			0,
			SCIPinfinity(scip),           
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
		));
		SCIP_CALL(SCIPaddCons(scip, cons));
		SCIP_CALL(SCIPcaptureCons(scip, cons));
		conss.push_back(cons);
		SCIP_CALL(SCIPreleaseCons(scip, &cons));
	}


	// add cardinality constraint
	if(has_knapcons){
		//knapweights = vector<SCIP_Real>(numvars, 1);
		SCIP_CALL(SCIPcreateConsLinear(
			scip,               /**< SCIP data structure */
			&cons,        /**< pointer to hold the created constraint */
			"knapsack",             /**< name of constraint */
			numvars,            /**< number of variables in the constraint */
			bin_vars.data(),    /**< array with variables of constraint entries */
			knapweights.data(),
			0,
			card,             
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
		));
		SCIP_CALL(SCIPaddCons(scip, cons));
		SCIP_CALL(SCIPcaptureCons(scip, cons));
		conss.push_back(cons);
		SCIP_CALL(SCIPreleaseCons(scip, &cons));
	}
	else if(has_cardcons){
		SCIPdebugMessage("%f\n", card);
		knapweights = vector<SCIP_Real>(numvars, 1);
		SCIP_CALL(SCIPcreateConsLinear(
			scip,               /**< SCIP data structure */
			&cons,        /**< pointer to hold the created constraint */
			"card",             /**< name of constraint */
			numvars,            /**< number of variables in the constraint */
			bin_vars.data(),    /**< array with variables of constraint entries */
			knapweights.data(),
			card,
			card,             
			TRUE,               /**< should the LP relaxation of constraint be in the initial LP?
														*   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
			TRUE,                /**< should the constraint be separated during LP processing?
														*   Usually set to TRUE. */
			TRUE,               /**< should the constraint be enforced during node processing?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be checked for feasibility?
														*   TRUE for model constraints, FALSE for additional, redundant constraints. */
			TRUE,               /**< should the constraint be propagated during node processing?
														*   Usually set to TRUE. */
			FALSE,
			FALSE,              /**< is constraint only valid locally?
														*   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
			FALSE,              /**< is constraint subject to aging?
														*   Usually set to FALSE. Set to TRUE for own cuts which
														*   are separated as constraints. */
			FALSE,              /**< should the relaxation be removed from the LP due to aging or cleanup?
														*   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
			FALSE               /**< should the constraint always be kept at the node where it was added, even
														*   if it may be moved to a more global node?
														*   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
		));
		SCIP_CALL(SCIPaddCons(scip, cons));
		SCIP_CALL(SCIPcaptureCons(scip, cons));
		conss.push_back(cons);
		SCIP_CALL(SCIPreleaseCons(scip, &cons));		
	}

	return SCIP_OKAY;
}

/**@} */
