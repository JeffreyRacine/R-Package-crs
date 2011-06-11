/*
 * 		Copyright (C) 2010 Zhenghua Nie. All Rights Reserved.
 * 		This code is published under the Common Public License.
 *
 * 		File:   snomadr.cpp
 * 		Author: Zhenghua Nie
 * 		Date:   Mon 16 May 2011
 *
 * 		We use ipoptr developed by Jelmer Ypma as the prototype of this package.
 * 		Some code is copied and edited from ipoptr. 
 * 		Please reference the license of ipoptr.
 *
 * 		This file defines the main function snomadRSolve that provides
 * 		an interface to NOMAD from R.
 * 		The function converts and R object containing objective function,
 * 		constraints, etc. into a Parameters , solves the problem,
 * 		and returns the result.
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License,  or
 *    (at your option) any later version.
 *      
 *    This program is distributed WITHOUT ANY WARRANTY. See the
 *    GNU General Public License for more details.
 *           
 *    If you do not have a copy of the GNU General Public License,  
 *    write to the Free Software Foundation, Inc., 
 *    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
#include "snomadr.h"

using namespace std;

static SEXP thefun,theenv;
SEXP showArgs1(SEXP largs);

bool file_exists(const char*path)
{
		struct stat sbuf;
		return stat(path,&sbuf)==0;
}

/* The following code is modified from  R_CheckUserInterrupt */
 void R_Interrupt(void)  
 { 
// R_CheckStack();		 
 /* This is the point where GUI systems need to do enough event 
 processing to determine whether there is a user interrupt event 
 pending. Need to be careful not to do too much event 
 processing though: if event handlers written in R are allowed 
 to run at this point then we end up with concurrent R 
 evaluations and that can cause problems until we have proper 
 concurrency support. LT */ 
//#if ( defined(HAVE_AQUA) || defined(Win32) ) 
// R_ProcessEvents(); 
// #else 
 onintr(); 
// #endif /* Win32 */ 
 } 

//
// Extracts element with name 'str' from R object 'list'
// and returns that element.
//
SEXP getListElement(SEXP list, std::string str)
{
		SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
		int i;
		for (i = 0; i < length(list); i++)
				if(str.compare(CHAR(STRING_ELT(names, i))) == 0) {
						elmt = VECTOR_ELT(list, i);
						break;
				}
		return elmt;
}


//
// Set options specified in R object opts in IpoptApplcation app.
// If we set the option to approximate the Hessian, then IpoptRNLP
// needs to know that, so it can return false when calling eval_h,
// instead of trying to find an R function that evaluates the Hessian.
//
// SEXP opts is an R list containing three sub-lists, with names
// integer, string and numeric. These sub-lists contain the actual
// options and there values that were specified by the user.
// Passing the options in this way makes it easier to call different
// SetValue function in IpoptApplication of the different types.
//
void setApplicationOptions(NOMAD::Parameters & p, SEXP opts ) {

		FILE *fp;
		const char opt_file[]="tmp.opt";

		fp=fopen(opt_file, "w");

		// extract the sub-lists with options of the different types into separate lists
		SEXP opts_integer = getListElement(opts, "integer");
		SEXP opts_numeric = getListElement(opts, "numeric");
		SEXP opts_string = getListElement(opts, "string");

		//	showArgs1(opts_integer);
		//	showArgs1(opts_numeric);
		//		showArgs1(opts_string);

		// loop over the integer options and set them
		SEXP opts_integer_names;
		opts_integer_names = getAttrib(opts_integer, R_NamesSymbol);
		for (int list_cnt=0;list_cnt<length( opts_integer );list_cnt++) {

				SEXP opt_value;
				PROTECT(opt_value = AS_INTEGER(VECTOR_ELT(opts_integer, list_cnt)));

				fprintf(fp, "%s\t%d\n", CHAR(STRING_ELT(opts_integer_names,  list_cnt)),  INTEGER(opt_value)[0]);

				UNPROTECT(1);	
		}

		// loop over the numeric options and set them
		SEXP opts_numeric_names;
		opts_numeric_names = getAttrib(opts_numeric, R_NamesSymbol);
		for (int list_cnt=0;list_cnt<length( opts_numeric );list_cnt++) {

				SEXP opt_value;
				PROTECT(opt_value = VECTOR_ELT(opts_numeric, list_cnt));

				fprintf(fp, "%s\t%.10g\n", CHAR(STRING_ELT(opts_numeric_names,  list_cnt)),  REAL(opt_value)[0]);

				UNPROTECT(1);	

		}


		// loop over the string options and set them
		SEXP opts_string_names;
		opts_string_names = getAttrib(opts_string, R_NamesSymbol);
		for (int list_cnt=0;list_cnt<length( opts_string );list_cnt++) {

				// opt_value will contain the first (should be the only one) element of the list
				SEXP opt_value;
				PROTECT(opt_value = STRING_ELT(VECTOR_ELT(opts_string, list_cnt),0));

				fprintf(fp, "%s\t%s\n", CHAR(STRING_ELT(opts_string_names,  list_cnt)),  CHAR(opt_value));


				UNPROTECT(1);	
		}

		fclose(fp);

		p.read(opt_file);

		remove(opt_file);

		return;
}


SEXP showArgs1(SEXP largs)
{
		int i, nargs = LENGTH(largs);
		Rcomplex cpl;
		SEXP el, names = getAttrib(largs, R_NamesSymbol);
		const char *name;

		for(i = 0; i < nargs; i++) {
				el = VECTOR_ELT(largs, i);
				name = isNull(names) ? "" : CHAR(STRING_ELT(names, i));
				switch(TYPEOF(el)) {
						case REALSXP:
								Rprintf("[%d] '%s' %f\n", i+1, name, REAL(el)[0]);
								break;
						case LGLSXP:
						case INTSXP:
								Rprintf("[%d] '%s' %d\n", i+1, name, INTEGER(el)[0]);
								break;
						case CPLXSXP:
								cpl = COMPLEX(el)[0];
								Rprintf("[%d] '%s' %f + %fi\n", i+1, name, cpl.r, cpl.i);
								break;
						case STRSXP:
								Rprintf("[%d] '%s' %s\n", i+1, name,
												CHAR(STRING_ELT(el, 0)));
								break;
						default:
								Rprintf("[%d] '%s' R type\n", i+1, name);
				}
		}
		return(R_NilValue);
}
//

double *eval_f( int m, int n, double *x)
{
		R_CheckUserInterrupt();

		SEXP rargs,Rcall,result;

		// Allocate memory for a vector of reals.
		// This vector will contain the elements of x,
		// x is the argument to the R function R_eval_f
		PROTECT(rargs = allocVector(REALSXP,n));
		for (int i=0;i<n;i++) {
				REAL(rargs)[i] = x[i];
		}

		// evaluate R function R_eval_f with the control x as an argument
		PROTECT(Rcall = lang2(thefun,rargs));
		PROTECT(result = eval(Rcall,theenv));

		// recode the return value from SEXP to Number
		double *ret_value;
		ret_value = (double *)malloc(sizeof(double)*m);
		for(int i=0;i<m;i++)
				ret_value[i]=REAL(result)[i];

		UNPROTECT(3);

		return(ret_value);
}


class RMy_Evaluator : public NOMAD::Evaluator {
		public:
				RMy_Evaluator  ( const NOMAD::Parameters & p ) :
						NOMAD::Evaluator ( p ) {};

				~RMy_Evaluator ( void ) {};

				bool eval_x ( NOMAD::Eval_Point   & x          ,
								const NOMAD::Double & h_max      ,
								bool                & count_eval   ) const {

						R_CheckUserInterrupt();

						double *x0;
						int n = x.get_n();
						int m = x.get_m();
						x0 = (double *)malloc(sizeof(double)*n);
						for(int i=0;i<n;i++)
						{
								x0[i]=x[i].value();
						}


						double *ret_value = eval_f(m, n, x0);
						for(int i=0;i<m;i++)
								x.set_bb_output  ( i , ret_value[i]  ); // objective value

						count_eval = true; // count a black-box evaluation

						free(x0);
						free(ret_value);

						return true;       // the evaluation succeeded
				}
};

SEXP print_solution(double obj_value, double *x, int n, int iter, int nmulti, NOMAD::stop_type status)
{
		int num_return_elements = 5;
		SEXP R_result_list;

		// R_result_list is a member object, which has been protected in the constructor
		// and will be unprotected in the destructor.
		PROTECT(R_result_list = allocVector(VECSXP, num_return_elements));

		// attach names to the return list
		SEXP names;
		PROTECT(names = allocVector(STRSXP, num_return_elements));

		SET_STRING_ELT(names, 0, mkChar("status"));
		SET_STRING_ELT(names, 1, mkChar("message"));
		SET_STRING_ELT(names, 2, mkChar("iterations"));
		SET_STRING_ELT(names, 3, mkChar("objective"));
		SET_STRING_ELT(names, 4, mkChar("solution"));
		setAttrib(R_result_list, R_NamesSymbol, names);

		// convert status to an R object
		SEXP R_status;


		PROTECT(R_status = allocVector(INTSXP,1));
		INTEGER(R_status)[0] = (int) status;


		SEXP R_status_message;
		PROTECT(R_status_message = allocVector(STRSXP, 1));
		if(nmulti < 1 ) {  /*0: snomadRSolve,  >0: smultnomadRSolve*/
				SET_STRING_ELT(R_status_message, 0, mkChar(stop_message[status]));
		}
		else
		{
				char *mes;
				mes = (char *)malloc(sizeof(char)*250);
				sprintf(mes, "Multiple mads runs - [%d]", nmulti);
				SET_STRING_ELT(R_status_message, 0, mkChar(mes));
				free(mes);
		}

		// convert value of objective function to an R object
		SEXP R_objective;
		PROTECT(R_objective = allocVector(REALSXP,1));
		REAL(R_objective)[0] = obj_value;

		// convert the value of the controls to an R object
		SEXP R_solution;
		PROTECT(R_solution = allocVector(REALSXP,n));
		for (int i=0;i<n;i++) {
				REAL(R_solution)[i] = x[i];
		}

		// convert number of iterations to an R object and add to the result_list
		SEXP R_num_iterations;
		PROTECT(R_num_iterations = allocVector(INTSXP,1));
		INTEGER(R_num_iterations)[0] = iter;


		// add elements to the list
		SET_VECTOR_ELT(R_result_list, 0, R_status);
		SET_VECTOR_ELT(R_result_list, 1, R_status_message);
		SET_VECTOR_ELT(R_result_list, 2, R_num_iterations);
		SET_VECTOR_ELT(R_result_list, 3, R_objective);
		SET_VECTOR_ELT(R_result_list, 4, R_solution);

		UNPROTECT(num_return_elements+2);

		return(R_result_list);
}
// we want this function to be available in R, so we put extern around it.
extern "C" {

		SEXP snomadRSolve( SEXP args )
		{
				R_CheckUserInterrupt();

				SEXP solution;

				//		showArgs1(args);

				PROTECT(solution=args);

				NOMAD::Display out(std::cout);
				out.precision(NOMAD::DISPLAY_PRECISION_STD);

				theenv = getListElement(args, "snomadr.environment");
				thefun = getListElement(args, "eval.f");

				try{
						R_CheckUserInterrupt();

						NOMAD::Parameters param(out);
						SEXP sN = getListElement(args,"n");
						int N = INTEGER(sN)[0];

						param.set_DIMENSION(N);

						// Default options
						//				param.set_DIRECTION_TYPE(NOMAD::GPS_2N_RAND);
						//				param.set_OPPORTUNISTIC_EVAL(true);
						//				param.set_MIN_POLL_SIZE(0.001);
						//				param.set_MIN_MESH_SIZE(0.001);
						//  			param.set_INITIAL_MESH_SIZE(0.01);
						param.set_MAX_BB_EVAL(10000);

						vector<NOMAD::bb_input_type> bbin(N);
						NOMAD::Point x0(N);
						NOMAD::Point lb(N);
						NOMAD::Point ub(N);

						SEXP sbbin = getListElement(args, "bbin");
						SEXP slb = getListElement(args, "lower.bounds");
						SEXP sub = getListElement(args, "upper.bounds");
						SEXP sx0 = getListElement(args, "x0");
						int  print_output = INTEGER(getListElement(args, "print.output"))[0];

						for(int i= 0;i<N;i++)
						{
								R_CheckUserInterrupt();

								switch(INTEGER(sbbin)[i])
								{
										case 0:
										default:
												bbin[i]= NOMAD::CONTINUOUS; 
												x0[i]= REAL(sx0)[i];
												lb[i]= REAL(slb)[i];
												ub[i]= REAL(sub)[i];
												break;
										case 1:
												bbin[i]= NOMAD::INTEGER; 
												x0[i]= (int)(REAL(sx0)[i]);
												lb[i]= (int)(REAL(slb)[i]);
												ub[i]= (int)(REAL(sub)[i]);
												break;
										case 2:
												bbin[i]= NOMAD::CATEGORICAL ; 
												x0[i]= (int)(REAL(sx0)[i]);
												lb[i]= (int)(REAL(slb)[i]);
												ub[i]= (int)(REAL(sub)[i]);
												break;
										case 3:
												bbin[i]= NOMAD::BINARY;
												x0[i]= (int)(REAL(sx0)[i]);
												lb[i]= (int)(REAL(slb)[i]);
												ub[i]= (int)(REAL(sub)[i]);
												break;

								}
						}
						param.set_BB_INPUT_TYPE(bbin);
						param.set_X0(x0);
						param.set_LOWER_BOUND(lb);
						param.set_UPPER_BOUND(ub);

						SEXP sbbout = getListElement(args, "bbout");

						vector<NOMAD::bb_output_type> bbot(LENGTH(sbbout));

						for(int i=0;i<LENGTH(sbbout);i++){
								bbot[i] = NOMAD::bb_output_type(INTEGER(sbbout)[i]);
						}

						param.set_BB_OUTPUT_TYPE(bbot);

						param.set_DISPLAY_STATS("bbe ( sol ) obj");

						/* set other options in R  */
						SEXP opts;
						opts = getListElement(args, "options");
						setApplicationOptions(param, opts );

						/* set other options in nomad.opt */
						if(file_exists("nomad.opt"))
								param.read("nomad.opt");

						param.check();
						RMy_Evaluator ev(param);

						NOMAD::Mads mads(param,&ev);
						NOMAD::stop_type status = mads.run();

						if ( status == NOMAD::CTRL_C ) R_Interrupt();  //)exit(-1);

						/* output the results of statistics, we may need to learn more about Stats */
						NOMAD::Stats &stats = mads.get_stats();

						if(print_output > 0)
						{
								printf("\niterations: %d\n", stats.get_iterations());
								printf("time:       %d\n", stats.get_real_time());
						}
						/* output the results */
						const NOMAD::Eval_Point *cur_x;
						cur_x= mads.get_best_feasible();

						NOMAD::Double obj_value = cur_x->get_f();

						NOMAD::Point best_x(N);
						best_x= *cur_x;

						double *sol_x;
						sol_x = (double *)malloc(sizeof(double)*N);

						for(int i=0;i<N;i++)
						{
								sol_x[i]= best_x[i].value();
						}

						solution = print_solution(obj_value.value(), sol_x, N, mads.get_cache().size(), 0,  status);

						free(sol_x);

				}
				catch(std::exception &e){
						std::cerr<<"\nNOMAD has been interrupted ("<<e.what()<<")\n\n";
				}

				NOMAD::Slave::stop_slaves(out);
				NOMAD::end();

				UNPROTECT(1);

				return(solution);
		}

} // extern "C"

/* for multiple restart */

void LH_values_for_var_i ( int     ind ,
				int     p   ,
				NOMAD::Point & x, const NOMAD::Point &lb, const NOMAD::Point  & ub, 
				const vector<NOMAD::bb_input_type> &bbin) {

		NOMAD::Random_Pickup rp(p);
		int    i;
		NOMAD::Double v;
		double UB = ub[ind].value();
		double LB = lb[ind].value();

		for ( i = 0 ; i < p ; ++i ) {
				double w = (UB - LB)/p;
				v = LB + ( i + rand()/NOMAD::D_INT_MAX ) * w;
				if( bbin[ind]!= NOMAD::CONTINUOUS) 
				{
						x[rp.pickup()] =(int) v.value();
				}
				else{
						x[rp.pickup()] = v;
				}
		}
}

/*----------------------------------------*/
/*  LH search used to generate x0 points  */
/*----------------------------------------*/
void LH_x0 ( int n , int p , vector<NOMAD::Point *> & x0_pts, const NOMAD::Point &lb, const NOMAD::Point & ub, const vector<NOMAD::bb_input_type>  &bbin) {

		//
		// pts contains n points of dimension p: each of these points contains
		// p different values for each variable:
		NOMAD::Point ** pts = new NOMAD::Point * [n] , * x;

		int pm1 = p - 1 , i;

		// creation of p search points:
		for ( int k = 0 ; k < p ; ++k ) {
				R_CheckUserInterrupt();

				x = new NOMAD::Point (n);

				for ( i = 0 ; i < n ; ++i ) {

						if ( k==0 ) {
								pts[i] = new NOMAD::Point(p);

								LH_values_for_var_i ( i , p , *pts[i], lb, ub, bbin );

						}

						(*x)[i] = (*pts[i])[k];

						if ( k == pm1 )
								delete pts[i];
				}

				x0_pts.push_back ( x );

		}

		delete [] pts;
}




// we want this function to be available in R, so we put extern around it.
extern "C" {

		SEXP smultinomadRSolve( SEXP args )
		{
				R_CheckUserInterrupt();
				SEXP solution;

				//	showArgs1(args);

				PROTECT(solution=args);

				NOMAD::Display out(std::cout);
				out.precision(NOMAD::DISPLAY_PRECISION_STD);

				theenv = getListElement(args, "snomadr.environment");
				thefun = getListElement(args, "eval.f");

				srand(unsigned(time(NULL)));
				rand();

				try{
						R_CheckUserInterrupt();

						vector<NOMAD::Point*> x0_pts;

						NOMAD::Parameters param(out);
						SEXP sN = getListElement(args,"n");
						int N = INTEGER(sN)[0];
						int i;

						SEXP snmulti = getListElement(args, "nmulti");
						int nmulti= INTEGER(snmulti)[0];

						param.set_DIMENSION(N);

						// Default options
						//				param.set_DIRECTION_TYPE(NOMAD::GPS_2N_RAND);
						//				param.set_OPPORTUNISTIC_EVAL(true);
						//				param.set_MIN_POLL_SIZE(0.001);
						//				param.set_MIN_MESH_SIZE(0.001);
						//	  		param.set_INITIAL_MESH_SIZE(0.01);
						param.set_MAX_BB_EVAL(100);


						vector<NOMAD::bb_input_type> bbin(N);
						NOMAD::Point lb(N);
						NOMAD::Point ub(N);

						SEXP sbbin = getListElement(args, "bbin");
						SEXP slb = getListElement(args, "lower.bounds");
						SEXP sub = getListElement(args, "upper.bounds");
						int  print_output = INTEGER(getListElement(args, "print.output"))[0];

						for( i= 0;i<N;i++)
						{
								R_CheckUserInterrupt();

								switch(INTEGER(sbbin)[i])
								{
										case 0:
										default:
												bbin[i]= NOMAD::CONTINUOUS; 
												lb[i]= REAL(slb)[i];
												ub[i]= REAL(sub)[i];
												break;
										case 1:
												bbin[i]= NOMAD::INTEGER; 
												lb[i]= (int)(REAL(slb)[i]);
												ub[i]= (int)(REAL(sub)[i]);
												break;
										case 2:
												bbin[i]= NOMAD::CATEGORICAL ; 
												lb[i]= (int)(REAL(slb)[i]);
												ub[i]= (int)(REAL(sub)[i]);
												break;
										case 3:
												bbin[i]= NOMAD::BINARY;
												lb[i]= (int)(REAL(slb)[i]);
												ub[i]= (int)(REAL(sub)[i]);
												break;

								}
						}
						param.set_BB_INPUT_TYPE(bbin);
						param.set_LOWER_BOUND(lb);
						param.set_UPPER_BOUND(ub);

						/* initial points */
						LH_x0(N,nmulti,x0_pts, lb, ub, bbin);

						SEXP sx0 = getListElement(args, "x0");  //we may use this as  best_x

						if(LENGTH(sx0) == N)    /* x0 as one of the  initial points */
						{
								for( i=0;i<N;i++) {
										(*x0_pts[0])[i]= REAL(sx0)[i];
								}
						}
						else {

								// read best_x.txt:
								ifstream fin ( "best_x.txt");

								if ( !fin.fail() )
										for ( i = 0 ; i < N ; ++i )
												fin >> (*x0_pts[0])[i];

								fin.close();
						}

						SEXP sbbout = getListElement(args, "bbout");

						vector<NOMAD::bb_output_type> bbot(LENGTH(sbbout));

						for( i=0;i<LENGTH(sbbout);i++){
								bbot[i] = NOMAD::bb_output_type(INTEGER(sbbout)[i]);
						}

						param.set_BB_OUTPUT_TYPE(bbot);

						param.set_DISPLAY_DEGREE (0);
						param.set_DISPLAY_STATS("bbe ( sol ) obj");

						/* set other options in R  */
						SEXP opts;
						opts = getListElement(args, "options");
						setApplicationOptions(param, opts );

						/* set other options in nomad.opt */
						if(file_exists("nomad.opt"))
								param.read("nomad.opt");


						param.set_X0 ( *x0_pts[0] );

						param.check();
						// display all starting points:
						if(print_output > 0 ){
								out << endl;
								for ( int j = 0 ; j < nmulti ; ++j )
										out << "starting point # " << j << ": ( " << *x0_pts[j] << " )" << endl;
								out << endl;
						}


						RMy_Evaluator ev(param);

						const NOMAD::Eval_Point * cur_x;
						NOMAD::Point              best_x (N);
						NOMAD::Double             best_f = NOMAD::INF , worst_f = 0.0 , avg_f = 0.0;

						// MADS runs:
						// ----------
						int bbe = 0;
						i = 0;
						NOMAD::stop_type status ;

						while ( true ) {

								R_CheckUserInterrupt();

								// algorithm creation:
								NOMAD::Mads mads ( param , &ev );
								status = mads.run();

								if ( status == NOMAD::CTRL_C ) R_Interrupt();  //)exit(-1);


								bbe += mads.get_cache().size();

								// displays and remember the best point:
								if( print_output > 0 ) 
										out << "run #" << setw(2) << i << ": ";
								cur_x = mads.get_best_feasible();
								if ( cur_x ) {

										if(print_output > 0 ) 
												out << "f=" << cur_x->get_f() << endl;

										if ( cur_x->get_f() < best_f ) {

												best_f = cur_x->get_f();
												best_x = *cur_x;
										}

										if ( cur_x->get_f() > worst_f )
												worst_f = cur_x->get_f();

										avg_f += cur_x->get_f();

								}
								else{
										if(print_output > 0 ) 
												out << "NULL" << endl;
								}

								if ( ++i == nmulti )
										break;

								// preparation of next run:
								mads.reset();
								param.reset_X0();
								param.set_X0 ( *x0_pts[i] );
								param.check();
						}


						// display the solution:
						if(print_output > 0 ) {
								out << endl << "bb eval : " << bbe << endl
										<< "best    : " << best_f;
								out << endl
										<< "worst   : " << worst_f << endl
										<< "solution: ";
								out << "x = ( ";
								best_x.display ( out , " " , -1 , -1 );
								out << " ) f(x) = " << best_f.value();
								out << endl << endl;
						}

						ofstream fout ( "best_x.txt" );
						fout << setprecision(32);
						best_x.display ( fout , " " , -1 , -1 );
						fout.close();

						// delete x0 points:
						for ( i = 0 ; i < nmulti ; ++i )
								delete x0_pts[i];



						double obj_value = best_f.value();

						double *sol_x;
						sol_x = (double *)malloc(sizeof(double)*N);

						for(int i=0;i<N;i++)
						{
								sol_x[i]= best_x[i].value();
						}

						solution = print_solution(obj_value, sol_x, N, bbe, nmulti, status);

						free(sol_x);

				}
				catch(std::exception &e){
						std::cerr<<"\nNOMAD has been interrupted ("<<e.what()<<")\n\n";
				}

				NOMAD::Slave::stop_slaves(out);
				NOMAD::end();

				UNPROTECT(1);

				return(solution);
		}

} // extern "C"
