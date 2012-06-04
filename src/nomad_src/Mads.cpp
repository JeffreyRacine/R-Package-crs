/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.5.1        */
/*                                                                                     */
/*  Copyright (C) 2001-2012  Mark Abramson        - the Boeing Company, Seattle        */
/*                           Charles Audet        - Ecole Polytechnique, Montreal      */
/*                           Gilles Couture       - Ecole Polytechnique, Montreal      */
/*                           John Dennis          - Rice University, Houston           */
/*                           Sebastien Le Digabel - Ecole Polytechnique, Montreal      */
/*                           Christophe Tribes    - Ecole Polytechnique, Montreal      */
/*                                                                                     */
/*  funded in part by AFOSR and Exxon Mobil                                            */
/*                                                                                     */
/*  Author: Sebastien Le Digabel                                                       */
/*                                                                                     */
/*  Contact information:                                                               */
/*    Ecole Polytechnique de Montreal - GERAD                                          */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada                  */
/*    e-mail: nomad@gerad.ca                                                           */
/*    phone : 1-514-340-6053 #6928                                                     */
/*    fax   : 1-514-340-5665                                                           */
/*                                                                                     */
/*  This program is free software: you can redistribute it and/or modify it under the  */
/*  terms of the GNU Lesser General Public License as published by the Free Software   */
/*  Foundation, either version 3 of the License, or (at your option) any later         */
/*  version.                                                                           */
/*                                                                                     */
/*  This program is distributed in the hope that it will be useful, but WITHOUT ANY    */
/*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    */
/*  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   */
/*                                                                                     */
/*  You should have received a copy of the GNU Lesser General Public License along     */
/*  with this program. If not, see <http://www.gnu.org/licenses/>.                     */
/*                                                                                     */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad               */
/*-------------------------------------------------------------------------------------*/
/**
  \file   Mads.cpp
  \brief  MADS algorithm (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-20
  \see    Mads.hpp
*/
#include "Mads.hpp"
#include "RNG.hpp"
using namespace std;  //zhenghua
/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
bool NOMAD::Mads::_force_quit          = false;
bool NOMAD::Mads::_flag_check_bimads   = true;
bool NOMAD::Mads::_flag_reset_mesh     = true;
bool NOMAD::Mads::_flag_reset_barriers = true;
bool NOMAD::Mads::_flag_p1_active      = false;

/*---------------------------------------------------------*/
/*       force quit (static, called by pressing ctrl-c)    */
/*---------------------------------------------------------*/
void NOMAD::Mads::force_quit ( int signalValue )
{
  NOMAD::Mads::_force_quit = true;
  NOMAD::Evaluator_Control::force_quit();
  NOMAD::Evaluator::force_quit();

#ifdef USE_TGP
  NOMAD::TGP_Output_Model::force_quit();
#endif
}

/*---------------------------------------------------------*/
/*                access to the flags (static)             */
/*---------------------------------------------------------*/
void NOMAD::Mads::get_flags ( bool & flag_check_bimads   ,
			      bool & flag_reset_mesh     ,			    
			      bool & flag_reset_barriers ,
			      bool & flag_p1_active        )
{
  flag_check_bimads   = _flag_check_bimads;
  flag_reset_mesh     = _flag_reset_mesh;
  flag_reset_barriers = _flag_reset_barriers;
  flag_p1_active      = _flag_p1_active;
}

/*---------------------------------------------------------*/
/*                      initializations                    */
/*---------------------------------------------------------*/
/*  . only to be invoked by constructors                   */
/*  . private method                                       */
/*---------------------------------------------------------*/
void NOMAD::Mads::init ( void )
{
  NOMAD::Mads::_force_quit = false;

  if ( !NOMAD::Slave::is_master() )
    return;

  // Mads::force_quit() will be called if ctrl-c is pressed:
  signal ( SIGINT  , NOMAD::Mads::force_quit );
#ifndef WINDOWS
  signal ( SIGPIPE , NOMAD::Mads::force_quit );  // (ctrl-c during a "| more")
#endif
#ifdef USE_MPI
  signal ( SIGTERM , NOMAD::Mads::force_quit );
#endif

  // random number generator seed initialization:
	bool valid_seed=NOMAD::RNG::set_seed(_p.get_seed());
	if ( !valid_seed )
		throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,"seed for random number generator not initialized properly!" );

  // model searches initialization:
    if ( _p.has_model_search() ) {
#ifdef USE_TGP
    if ( _p.get_model_search(1) == NOMAD::TGP_MODEL )
      _model_search1 = new TGP_Model_Search ( _p );
#endif
    if ( _p.get_model_search(1) == NOMAD::QUADRATIC_MODEL )
      _model_search1 = new Quad_Model_Search ( _p );
#ifdef USE_TGP
    if ( _p.get_model_search(2) == NOMAD::TGP_MODEL )
      _model_search2 = new TGP_Model_Search ( _p );
#endif
    if ( _p.get_model_search(2) == NOMAD::QUADRATIC_MODEL )
      _model_search2 = new Quad_Model_Search ( _p );
  }

#ifdef USE_TGP
  _ev_control.set_last_TGP_model ( NULL );
#endif

  // VNS search initialization:
  if ( _p.get_VNS_search() )
    _VNS_search = new VNS_Search ( _p );

  // cache search initialization:
  if ( _p.get_cache_search() )
    _cache_search = new Cache_Search ( _p );
}

/*---------------------------------------------------------*/
/*                       destructor                        */
/*---------------------------------------------------------*/
NOMAD::Mads::~Mads ( void )
{
  delete _pareto_front;
  delete _model_search1;
  delete _model_search2;
  delete _VNS_search;
  delete _cache_search;
  delete _L_curve;

  if ( _extended_poll && !_user_ext_poll)
    delete _extended_poll;
}

/*---------------------------------------------------------*/
/*                         reset                           */
/*---------------------------------------------------------*/
/*  default values for parameters: keep_barriers = false   */
/*                                 keep_stats    = false   */
/*---------------------------------------------------------*/
void NOMAD::Mads::reset ( bool keep_barriers , bool keep_stats )
{
  // evaluator control:
#ifdef USE_TGP
  _ev_control.set_last_TGP_model ( NULL );
#endif

  // user search:
  _user_search = NULL;

  // model search #1:
  if ( _p.get_model_search(1) != NOMAD::NO_MODEL ) {
    if ( _model_search1 )
      _model_search1->reset();
    else {
      if ( _p.get_model_search(1) == NOMAD::TGP_MODEL ) {
#ifdef USE_TGP
	_model_search1 = new TGP_Model_Search  ( _p ) ;
#endif
      }
      else
	_model_search1 = new Quad_Model_Search ( _p );
    }
  }
  else {
    delete _model_search1;
    _model_search1 = NULL;
  }

  // model search #2:
  if ( _p.get_model_search(2) != NOMAD::NO_MODEL ) {
    if ( _model_search2 )
      _model_search2->reset();
    else {
      if ( _p.get_model_search(2) == NOMAD::TGP_MODEL ) {
#ifdef USE_TGP
	_model_search2 = new TGP_Model_Search  ( _p ) ;
#endif
      }
      else
	_model_search2 = new Quad_Model_Search ( _p );
    }
  }
  else {
    delete _model_search2;
    _model_search2 = NULL;
  }

  // VNS search:
  if ( _p.get_VNS_search() ) {
    if ( _VNS_search )
      _VNS_search->reset();
    else
      _VNS_search = new VNS_Search ( _p );
  }
  else {
    delete _VNS_search;
    _VNS_search = NULL;
  }

  // cache search:
  if ( _p.get_cache_search() ) {
    if ( _cache_search )
      _cache_search->reset();
    else
      _cache_search = new Cache_Search ( _p );
  }
  else {
    delete _cache_search;
    _cache_search = NULL;
  }

  // barriers:
  _flag_reset_barriers = !keep_barriers;
  if ( _flag_reset_barriers ) {
    _true_barrier.reset();
    _sgte_barrier.reset();
  }

  // first success:
  _first_succ = NULL;

  // stats:
  if ( !keep_stats )
    _stats.reset();

  // mesh:
  NOMAD::Mesh::init ( _p.get_mesh_update_basis().value() ,
		      _p.get_mesh_coarsening_exponent () ,
		      _p.get_mesh_refining_exponent   () ,
		      _p.get_initial_mesh_index       ()   );
}

/*---------------------------------------------------------*/
/*            algorithm execution (single-objective)       */
/*---------------------------------------------------------*/
NOMAD::stop_type NOMAD::Mads::run ( void )
{
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_gen_dd();
  NOMAD::stop_type  stop_reason = NOMAD::UNKNOWN_STOP_REASON;

#ifdef USE_MPI

  // init the slaves:
  bool stop_slaves_here = false;

  if ( NOMAD::Slave::is_master() ) {
    if ( !NOMAD::Slave::are_running() ) {
      NOMAD::Slave::init_slaves ( out );
      stop_slaves_here = true;
    }
  }
  else {
    NOMAD::Slave s ( _p , _ev_control.get_evaluator() );
    s.run();
    return stop_reason;
  }

#endif

  try {

    // check an extended poll if there are categorical
    // variables and disable extended poll otherwise:
    if ( _p.get_signature()->has_categorical() ) {

      if ( _user_ext_poll && !_extended_poll )
	throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
	"categorical variables: user extended poll object is NULL" );

      if ( _p.get_extended_poll_enabled() && !_user_ext_poll ) {

	if ( _extended_poll )
	  delete _extended_poll;
	_extended_poll = new NOMAD::Extended_Poll ( _p );

	std::string error_str;	
	if ( !_extended_poll->set_neighbors_exe ( error_str ) )
	  throw NOMAD::Exception ( "Mads.cpp" , __LINE__ , error_str );
      }
    }
    else if ( _extended_poll ) {
      if ( !_user_ext_poll )
	delete _extended_poll;
      _extended_poll = NULL;
    }
    
    // check if Mads::run() has been called for multi-objective:
    if ( NOMAD::Mads::_flag_check_bimads && _p.get_nb_obj() > 1 )
      throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
	    "Mads::run() called for multi-objective instead of Mads::multi_run()" );

#ifndef R_VERSION
    if ( display_degree == NOMAD::NORMAL_DISPLAY ||  display_degree == NOMAD::FULL_DISPLAY )
       out << std::endl << NOMAD::open_block ( "MADS run" );
       
    if ( display_degree == NOMAD::NORMAL_DISPLAY ) {
      _ev_control.display_stats ( true                   ,
    				  out                    ,
    				  _p.get_display_stats() ,
    				  NULL                   ,
    				  false                  ,
    				  NULL                     );
      out << std::endl << std::endl;
    }
#endif

    // barriers init:
    if ( _flag_reset_barriers ) {
      _true_barrier.reset();
      _sgte_barrier.reset();
    }
    
    // first success reset:
    const NOMAD::Point * old_fs = _first_succ;
    _first_succ = NULL;

    // evaluator control init:
    _ev_control.reset();

    // reset the extended poll:
    if ( _extended_poll && _p.get_extended_poll_enabled() )
      _extended_poll->reset();

    // mesh init/reset:
    if ( _flag_reset_mesh )
      NOMAD::Mesh::init ( _p.get_mesh_update_basis().value() ,
			  _p.get_mesh_coarsening_exponent () ,
			  _p.get_mesh_refining_exponent   () ,
			  _p.get_initial_mesh_index       ()   );
    
    NOMAD::success_type       success , last_success;
    int                       nb_search_pts;
    bool                      count_search;
    bool                      stop           = false;
    const NOMAD::Eval_Point * new_feas_inc   = NULL;
    const NOMAD::Eval_Point * new_infeas_inc = NULL;
    
    stop_reason = NOMAD::NO_STOP;

    // x0 eval:
    eval_x0 ( stop , stop_reason );

    // phase one: if no feasible starting point:
    bool phase_one_done = false;
    if ( stop                          &&
	 stop_reason == NOMAD::X0_FAIL &&
	 _p.has_EB_constraints()       &&
	 ( _stats.get_eval() > 0 ||
	   ( _p.get_opt_only_sgte() && _stats.get_sgte_eval() > 0 ) ) ) {
      
      phase_one_done = true;
      Phase_One_Search p1s ( _p );
      p1s.search ( *this          ,
		   nb_search_pts  ,
		   stop           ,
		   stop_reason    ,
		   success        ,
		   count_search   ,
		   new_feas_inc   ,
		   new_infeas_inc   );
    }
    
    // initial Latin-Hypercube (LH) search:
    if ( !stop && !phase_one_done && _p.get_LH_search_p0() > 0 ) {  
      
      LH_Search lh ( _p , true , _flag_p1_active );
      int       nb_search_pts;
      
      lh.search ( *this          ,
		  nb_search_pts  ,
		  stop           ,
		  stop_reason    ,
		  success        ,
		  count_search   ,
		  new_feas_inc   ,
		  new_infeas_inc   );
      
      if ( success == NOMAD::FULL_SUCCESS )
	_stats.add_LH_success();
      if ( count_search )
	_stats.add_nb_LH_searches();
      _stats.add_LH_pts ( nb_search_pts );
    }
    
    // no iterations allowed:
    if ( !stop && _p.get_max_iterations() == 0 ) {
      stop        = true;
      stop_reason = NOMAD::MAX_ITER_REACHED;
    }
    
    // L_curve initialization:
    delete _L_curve;
    _L_curve = NULL;
    const NOMAD::Double L_curve_target = _p.get_L_curve_target();
    if ( L_curve_target.is_defined() ) {
      _L_curve = new NOMAD::L_Curve ( L_curve_target );
      const NOMAD::Eval_Point * best_feasible = get_best_feasible();
      if ( best_feasible )
	_L_curve->insert ( _stats.get_bb_eval() , best_feasible->get_f() );
    }

    int max_cfi = _p.get_max_consecutive_failed_iterations();
    int nb_cfi  = 0;

    success = last_success = NOMAD::UNSUCCESSFUL;

    // MADS iterations:
    while ( !stop ) {

      iteration ( stop           ,
		  stop_reason    ,
		  success        ,
		  new_feas_inc   ,
		  new_infeas_inc   );

      if ( success == NOMAD::UNSUCCESSFUL && last_success == NOMAD::UNSUCCESSFUL )
	++nb_cfi;
      else
	nb_cfi = (success == NOMAD::UNSUCCESSFUL) ? 1 : 0;

      last_success = success;

      // check the consecutive number of failed iterations:
      if ( max_cfi > 0 && nb_cfi > max_cfi ) {
	stop        = true;
	stop_reason = NOMAD::MAX_CONS_FAILED_ITER;
      }

    }

    // parallel version:
#ifdef USE_MPI
    
    // asynchronous mode: wait for the evaluations in progress:
    if ( _p.get_asynchronous() ) {
      std::list<const NOMAD::Eval_Point *> evaluated_pts;
      _ev_control.wait_for_evaluations ( NOMAD::ASYNCHRONOUS ,
					 _true_barrier       ,
					 _sgte_barrier       ,
					 _pareto_front       ,
					 stop                ,
					 stop_reason         ,
					 success             ,
					 evaluated_pts         );
    }
    
    // update stats:
    _stats.set_MPI_data_size ( NOMAD::Slave::get_data_sent() +
			       NOMAD::Slave::get_data_rcvd()   );
    
#endif
    
    // final cache save (overwrite=true):
    _ev_control.save_caches ( true );
    
    // final displays:
    const NOMAD::Eval_Point * bf = get_best_feasible();
    bool write_stats = bf &&
      ( bf->get_tag()        != _ev_control.get_last_stats_tag() ||
	_stats.get_bb_eval() != _ev_control.get_last_stats_bbe()    );
    
    const std::string & stats_file_name = _p.get_stats_file_name();

    if ( !stats_file_name.empty() ) {
      if ( write_stats && !_p.get_display_all_eval() ) {
	_ev_control.stats_file ( stats_file_name , bf , true , NULL );
      }
      if ( !bf ) {
	std::ofstream fout ( (_p.get_problem_dir() + stats_file_name).c_str() );
	if ( fout.fail() ) {
      	  if ( out.get_gen_dd() != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
	    out << std::endl
		<< "Warning (" << "Mads.cpp" << ", " << __LINE__
		<< "): could not save information in stats file \'"
		<< stats_file_name << "\'" << std::endl << std::endl;
	}
	else
	  fout << "no feasible solution has been found in "
	       << _stats.get_bb_eval() << " evaluations"
	       << std::endl;
	fout.close();
      }
    }
    
    if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY ) {   
      
      // final stats:
      if ( display_degree == NOMAD::NORMAL_DISPLAY && bf && write_stats &&
	   !_p.get_display_all_eval() )
	_ev_control.display_stats ( false                  ,
				    out                    ,
				    _p.get_display_stats() ,
				    bf                     ,
				    true                   ,
				    NULL                     );
#ifndef R_VERSION
      std::ostringstream msg;
      msg << "end of run (" << stop_reason << ")";
      out << std::endl << NOMAD::close_block ( msg.str() );
#endif
    }
    
    // mono-objective final displays:
    if ( _p.get_nb_obj() == 1 ) {
      
      if ( display_degree == NOMAD::FULL_DISPLAY )
	out << std::endl << NOMAD::open_block ( "NOMAD final display" );

#ifndef R_VERSION
      display();
#endif
      
      if ( display_degree == NOMAD::FULL_DISPLAY )
	out.close_block();
    }
    
    // reset the first success to its previous value:
    _first_succ = old_fs;

  } // end of the try block
  
  catch ( std::exception & e ) {
    
#ifdef USE_MPI
    if ( NOMAD::Slave::are_running() )
      NOMAD::Slave::stop_slaves ( out );
#endif

    throw NOMAD::Exception ( "Mads.cpp" , __LINE__ , e.what() );
  }
  
  // stop the slaves:
#ifdef USE_MPI
  if ( NOMAD::Slave::are_running() && stop_slaves_here )
    NOMAD::Slave::stop_slaves ( out );
#endif

  return stop_reason;
}

/*----------------------------------------------------------------------*/
/*      launch a single optimization for multi-objective optimization   */
/*----------------------------------------------------------------------*/
/*  . the display_degree is given as a parameter since it corresponds   */
/*    to the original iterative display degree before all degrees have  */
/*    been set to zero                                                  */
/*  . private method                                                    */
/*----------------------------------------------------------------------*/
void NOMAD::Mads::multi_launch_single_opt
( NOMAD::dd_type               display_degree ,
  int                          mads_runs      ,
  int                          overall_bbe    ,
  NOMAD::Multi_Obj_Evaluator & ev             , 
  int                        & stagnation_cnt ,
  NOMAD::Stats               & multi_stats    ,
  bool                       & stop           ,
  NOMAD::stop_type           & stop_reason      )
{
  // max number of bb evaluations for one MADS run:
  int max_bbe = _p.get_max_bb_eval();

  // size of the Pareto front before the MADS run:
  int tmp = _pareto_front->size();

  // current MADS run:
  int cur_mads_run = multi_stats.get_mads_runs();

  // displays:
  const NOMAD::Display & out = _p.out();

  if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY) {
    out << "MADS run " << std::setw(2) << cur_mads_run + 1;
    if ( mads_runs > 0 )
      out << "/" << mads_runs;
    out << " ...";
  }

  // run single-objective MADS (it also updates the Pareto front):
  NOMAD::Mads::set_flag_check_bimads ( false );
  NOMAD::stop_type single_run_stop_reason = run();
  NOMAD::Mads::set_flag_check_bimads ( true );

  if ( single_run_stop_reason == NOMAD::CTRL_C              ||
       single_run_stop_reason == NOMAD::ERROR               ||
       single_run_stop_reason == NOMAD::UNKNOWN_STOP_REASON ||
       single_run_stop_reason == NOMAD::X0_FAIL             ||
       single_run_stop_reason == NOMAD::F_TARGET_REACHED    ||
       single_run_stop_reason == NOMAD::P1_FAIL                ) {
    stop        = true;
    stop_reason = single_run_stop_reason;
  }

  // update MULTI-MADS stats from MADS stats:
  multi_stats.update ( _stats , false ); // for_search = false
  multi_stats.add_mads_run();

  int nb_new_pts = _pareto_front->size() - tmp;
  int global_bbe = multi_stats.get_bb_eval();

  // displays:
  if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY) {

    // display basic stats on the terminated run:
    out << "... OK [bb eval="    << std::setw(3) << _stats.get_bb_eval()
	<< "] [overall bb eval=" << std::setw(5) << global_bbe
	<< "] [# dominant pts="  << std::setw(4) << _pareto_front->size()
	<< "] [# new pts="       << std::setw(4) << nb_new_pts << "]";
	
    // display f1, f2, and f:
    const NOMAD::Eval_Point * bf = get_best_feasible();
    if ( bf ) {

      const NOMAD::Point & bbo = bf->get_bb_outputs();

      out << " [f1=" << bbo[ev.get_i1()]
	  << " f2=" << bbo[ev.get_i2()];
      if ( display_degree == NOMAD::FULL_DISPLAY )
	out << " f="  << bf->get_f();
      out << "]";
    }
    out << std::endl;
  }

  if ( _stats.get_bb_eval() == 0 && nb_new_pts == 0 )
    ++stagnation_cnt;
  else
    stagnation_cnt = 0;

  // stop ?
  if ( !stop ) {

    // test the number of MADS runs:
    if ( mads_runs > 0 ) {
      if ( multi_stats.get_mads_runs() >= mads_runs ) {
	stop        = true;
	stop_reason = NOMAD::MULTI_NB_MADS_RUNS_REACHED;
      }
    }

    // test if no new Pareto point has been generated for 50*n MADS runs:
    else {
      if ( stagnation_cnt > 50 * _p.get_nb_free_variables() ) {
	stop        = true;
	stop_reason = NOMAD::MULTI_STAGNATION;
      }
    }
  }

  if ( overall_bbe >= 0 && global_bbe >= overall_bbe ) {
    stop        = true;
    stop_reason = NOMAD::MULTI_MAX_BB_REACHED;
  }
  
  bool user_calls_enabled = _p.get_user_calls_enabled();

  if ( !stop ) {

    // ell is the mesh index on which the last run terminated:
    // int ell = NOMAD::Mesh::get_mesh_index();

    // reset MADS:
    reset();

    // this strategy deciding the initial mesh size
    // was used with versions < 3.4
    // if ( cur_mads_run > 1 )
    //  _p.set_INITIAL_MESH_INDEX ( (ell > 5) ? 5 : ell );

    // modify MAX_BB_EVAL for single runs (in order to have
    // less than overall_bbe blackbox evaluations):
    if ( overall_bbe >= 0 && global_bbe + max_bbe > overall_bbe )
      _p.set_MAX_BB_EVAL ( overall_bbe - global_bbe );
  }
  
  // set the number of MADS runs for the general Stats object:
  _stats.set_mads_runs ( multi_stats.get_mads_runs() );

  // call the user-defined function Multi_Obj_Evaluator::update_mads_run():
  if ( user_calls_enabled )
    ev.update_mads_run ( _stats         ,
			 _ev_control    ,
			 _true_barrier  ,
			 _sgte_barrier  ,
			 *_pareto_front   );
}

/*--------------------------------------------------------------------------*/
/*  compute and set the minimal poll size for multi-objective optimization  */
/*  (private)                                                               */
/*--------------------------------------------------------------------------*/
void NOMAD::Mads::multi_set_min_poll_size ( const NOMAD::Point & lb        ,
					    const NOMAD::Point & ub        ,
					    const NOMAD::Point & delta_p_0 ,
					    NOMAD::Double        delta_j     )
{
  delta_j /= sqrt ( NOMAD::Mesh::get_mesh_update_basis() );

  int          n = delta_p_0.size();
  NOMAD::Point delta_p_min (n);

  for ( int i = 0 ; i < n ; ++i ) {

    // set a relative value:
    if ( lb[i].is_defined() && ub[i].is_defined() )
      delta_p_min[i] = delta_j * ( ub[i] - lb[i] );

    // set an absolute value:
    else
      delta_p_min[i] = delta_j;

    // compare to delta_p_0:
    if ( delta_p_min[i] > delta_p_0[i] )
      delta_p_min[i] = delta_p_0[i];
  }

  _p.set_MIN_POLL_SIZE ( delta_p_min );
}

/*---------------------------------------------------------*/
/*                 algorithm execution (multi)             */
/*---------------------------------------------------------*/
NOMAD::stop_type NOMAD::Mads::multi_run ( void )
{
  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_gen_dd();
  NOMAD::stop_type  stop_reason = NOMAD::UNKNOWN_STOP_REASON;

  // init the slaves:

#ifdef USE_MPI

  bool stop_slaves_here = false;

  if ( NOMAD::Slave::is_master() ) {
    if ( !NOMAD::Slave::are_running() ) {
      NOMAD::Slave::init_slaves ( out );
      stop_slaves_here = true;
    }
  }
  else {
    NOMAD::Slave s ( _p , _ev_control.get_evaluator() );
    s.run();
    return stop_reason;
  }

#endif

  try {

    // objective indexes:
    NOMAD::Multi_Obj_Evaluator::set_obj_indexes ( _p.get_index_obj() );

    // bounds:
    const NOMAD::Point & lb = _p.get_lb();
    const NOMAD::Point & ub = _p.get_ub();

    // MULTI-MADS stopping criteria:
    int  mads_runs      = _p.get_multi_nb_mads_runs();    // max number of MADS runs
    int  overall_bbe    = _p.get_multi_overall_bb_eval(); // max number of total bb eval.
    bool use_delta_crit = _p.get_multi_use_delta_crit();  // use the delta term. crit.
    int  stagnation_cnt = 0;

    if ( mads_runs == 0 )
      throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
	    "Mads::multi_run(): parameter MULTI_NB_MADS_RUNS is not positive" );

    if ( _p.get_nb_obj() != 2 )
      throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
      "Mads::multi_run(): NOMAD current version handles a maximum of two objectives" );

    // remember cache save period:
    int old_csp = _p.get_cache_save_period();

    // remember L_CURVE_TARGET:
    NOMAD::Double old_lct = _p.get_L_curve_target();
    
    // remember solution file:
    std::string old_sol_file = _p.get_solution_file();
    
    // remember the original LH search parameters:
    int lh_p0 = _p.get_LH_search_p0();
    int lh_pi = _p.get_LH_search_pi();
    
    // remember the original minimal poll size:
    const NOMAD::Point original_min_poll_size = _p.get_min_poll_size();
    
    // remember display degrees:
    NOMAD::dd_type iter_dd = out.get_iter_dd();
    std::string    old_dd;
    out.get_display_degree ( old_dd );

    // save list of starting points:
    std::string x0_cache_file = _p.get_x0_cache_file();
    std::vector<NOMAD::Point *> x0s;
    {
      const std::vector<NOMAD::Point *> & x0s_tmp = _p.get_x0s();
      size_t nx0 = x0s_tmp.size() , k;
      for ( k = 0 ; k < nx0 ; ++k )
	x0s.push_back ( new Point ( *x0s_tmp[k] ) );
    }

    // compute delta_p_0:
    NOMAD::Point delta_p_0 ( _p.get_dimension() );
    _p.get_signature()->get_mesh().get_delta_p ( delta_p_0                   ,
						 _p.get_initial_mesh_index()   );
    
    if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
      out << std::endl << NOMAD::open_block ( "multi-MADS run" ) << std::endl;
    
    bool stop   = false;
    stop_reason = NOMAD::NO_STOP;
    
    // MULTI-MADS stats:
    NOMAD::Stats multi_stats ( _stats );
    
    // access to the evaluator (downcast to a Multi_Obj_Evaluator):
    NOMAD::Multi_Obj_Evaluator * ev =
      static_cast<NOMAD::Multi_Obj_Evaluator*> ( _ev_control.get_evaluator() );
    if ( !ev->is_multi_obj() )
      throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
      "Mads::multi_run(): associated Evaluator object is not a Multi_Obj_Evaluator" );

    // parameters modifications:
    // -------------------------
    
    // STATS_FILE:
    const std::string            old_stats_file_name = _p.get_stats_file_name();
    const std::list<std::string> old_stats_file      = _p.get_stats_file();
    _p.reset_stats_file();
    
    // MAX_BB_EVAL:
    int max_bbe = _p.get_max_bb_eval();
    if ( overall_bbe >= 0 && ( max_bbe < 0 || overall_bbe < max_bbe ) )
      _p.set_MAX_BB_EVAL ( overall_bbe );
    
    // disable display:
    _p.set_DISPLAY_DEGREE ( NOMAD::NO_DISPLAY );
    
    // disable solution file:
    _p.set_SOLUTION_FILE ( "" );
    
    // disable CACHE_SAVE_PERIOD:
    _p.set_CACHE_SAVE_PERIOD ( -1 );
    
    // disable L_CURVE_TARGET:
    _p.set_L_CURVE_TARGET ( NOMAD::Double() );
    
    // LH_SEARCH and MAX_BB_EVAL adjustment:
    if ( lh_p0 > 0 ) {
      _p.set_LH_SEARCH ( lh_p0 , 0 );
      if ( max_bbe >= 0 ) {
	int bbe = max_bbe + lh_p0;
	if ( overall_bbe >= 0 && bbe > overall_bbe )
	  bbe = overall_bbe;
	_p.set_MAX_BB_EVAL ( bbe );
      }
    }
    
    // parameters validation:
    _p.check ( true ,    // remove_history_file  = true
	       true ,    // remove_solution_file = true
	       true   ); // remove_stats_file    = true

    // Pareto front initialization:
    delete _pareto_front;
    _pareto_front = new NOMAD::Pareto_Front;
    
    // initial optimizations ( minimize f1(x) or f2(x) ):
    // --------------------------------------------------
    const NOMAD::Eval_Point * best_f2;
    int i;

    for ( i = 0 ; i < 2 ; ++i ) {
      
      if ( stop )
	break;

      // minimize f2:
      if ( i == 1 ) {
	
	// new starting point:
	best_f2 = _pareto_front->get_best_f2();
	if ( best_f2 ) {  
	  _p.set_EXTERN_SIGNATURE ( best_f2->get_signature() );
	  _p.reset_X0();
	  _p.set_X0 ( *best_f2 );
	}
	
	// LH_SEARCH:
	if ( lh_pi > 0 )
	  _p.set_LH_SEARCH ( lh_pi , 0 );
	else if ( lh_p0 > 0 )
	  _p.set_LH_SEARCH ( 0 , 0 );
	
	// MAX_BB_EVAL:
	if ( max_bbe >= 0 ) {
	  int bbe = max_bbe + ( (lh_pi > 0 ) ? lh_pi : 0 );
	  if ( overall_bbe >= 0 ) {
	    if ( bbe > overall_bbe )
	      bbe = overall_bbe;
	    int global_bbe = multi_stats.get_bb_eval();
	    if ( global_bbe + bbe > overall_bbe )
	      bbe = overall_bbe - global_bbe;
	  }
	  _p.set_MAX_BB_EVAL ( bbe );
	}
	
	if ( _p.to_be_checked() )
	  _p.check ( false ,    // remove_history_file  = false
		     true  ,    // remove_solution_file = true
		     true    ); // remove_stats_file    = true
      }
      
      // set weights/reference:
      ev->set_weights ( 1-i , i );
      ev->set_ref     ( NULL    );
      
      // launch the single optimization:
      multi_launch_single_opt ( iter_dd        ,
				mads_runs      ,
				overall_bbe    ,
				*ev            ,
				stagnation_cnt ,
				multi_stats    ,
				stop           ,
				stop_reason      );
    }
    
    const NOMAD::Point        * ref;
    const NOMAD::Pareto_Point * xj;
    NOMAD::Double               delta_j;
    
    // the LH search is disabled:
    _p.set_LH_SEARCH ( 0 , 0 );
    
    // MAX_BB_EVAL reset:
    if ( max_bbe > 0 && ( lh_p0 > 0 || lh_pi > 0 ) ) {
      int bbe = max_bbe;
      if ( overall_bbe >= 0 ) {
	if ( bbe > overall_bbe )
	  bbe = overall_bbe;
	
	int global_bbe = multi_stats.get_bb_eval();
	if ( global_bbe + bbe > overall_bbe )
	  bbe = overall_bbe - global_bbe;
      }
      _p.set_MAX_BB_EVAL ( bbe );
    }

    // MULTI-MADS main loop:
    // ---------------------
    const NOMAD::Eval_Point * x0_tmp;

    while ( !stop ) {
      
      // get the reference point from the Pareto front:
      ref = _pareto_front->get_ref ( xj , delta_j );
      
      if ( !xj ) {
	stop        = true;
	stop_reason = NOMAD::MULTI_NO_PARETO_PTS;
	break;
      }
      
      // use delta as stopping criterion:
      if ( use_delta_crit ) {
	if ( delta_j.is_defined() && delta_j > 0.0 )
	  multi_set_min_poll_size ( lb , ub , delta_p_0 , delta_j );
	else
	  _p.set_MIN_POLL_SIZE ( original_min_poll_size );
      }
      
      // new starting point:
      x0_tmp = xj->get_element();
      _p.set_EXTERN_SIGNATURE ( x0_tmp->get_signature() );
      _p.reset_X0();
      _p.set_X0 ( *x0_tmp );

      _p.check ( false ,    // remove_history_file  = false
		 true  ,    // remove_solution_file = true
		 true    ); // remove_stats_file    = true

      // a reference point has been found: optimization
      // with reference-based function:
      if ( ref ) {
	
	// set reference:
	ev->set_ref ( ref );
	
	// launch the single optimization:
	multi_launch_single_opt ( iter_dd        ,
				  mads_runs      ,
				  overall_bbe    ,
				  *ev            ,
				  stagnation_cnt ,
				  multi_stats    ,
				  stop           ,
				  stop_reason      );
	
	delete ref;
	ev->set_ref ( NULL );
      }
      
      // no reference available: two optimizations ( f1(x) and f2(x) ):
      else {

	// for the stagnation check:
	const NOMAD::Eval_Point * pp_before;
	int  stagnation_cnt_before , overall_bbe_before;
	bool check_1 = false;

	// loop on f1 and f2:
	for ( i = 0 ; i < 2 ; ++i ) {
	  
	  if ( stop )
	    break;
	  
	  // minimize f2:
	  if ( i == 1 ) {
	
	    // new starting point:
	    best_f2 = _pareto_front->get_best_f2();
	    if ( best_f2 ) {  
	      _p.set_EXTERN_SIGNATURE ( best_f2->get_signature() );
	      _p.reset_X0();
	      _p.set_X0 ( *best_f2 );
	    }
	    else
	      _p.set_X0 ( *x0_tmp );
	    	
	    _p.check ( false ,    // remove_history_file  = false
		       true  ,    // remove_solution_file = true
		       true    ); // remove_stats_file    = true
	  }

	  // set weights/reference:
	  ev->set_weights ( 1-i , i );
	  ev->set_ref     ( NULL    );

	  stagnation_cnt_before = stagnation_cnt;
	  overall_bbe_before    = overall_bbe;
	  pp_before = ( _pareto_front->size() == 1 ) ?
	    _pareto_front->begin() : NULL;

	  // launch the single optimization:
	  multi_launch_single_opt ( iter_dd        ,
				    mads_runs      ,
				    overall_bbe    ,
				    *ev            ,
				    stagnation_cnt ,
				    multi_stats    ,
				    stop           ,
				    stop_reason      );

	  // stagnation check:
 	  if ( stagnation_cnt > stagnation_cnt_before &&
	       overall_bbe == overall_bbe_before      &&
	       _pareto_front->size() == 1             &&
	       _pareto_front->begin() == pp_before       ) {
	    
	    if ( i == 0 )
	      check_1 = true;
	    else if ( check_1 ) {
	      stop        = true;
	      stop_reason = NOMAD::MULTI_STAGNATION;
	    }
	  }
	}
      }

    } // end of MULTI-MADS main loop
      // ---------------------------

    // parameters re-initialization and final displays:
    {
      _p.reset_X0();
      size_t nx0 = x0s.size();
      if ( nx0 > 0 ) {
	for ( size_t k = 0 ; k < nx0 ; ++k ) {
	  _p.set_X0 ( *x0s[k] );
	  delete x0s[k];
	}
      }
      else if ( !x0_cache_file.empty() )
	_p.set_X0 ( x0_cache_file );
    }
    _p.set_MAX_BB_EVAL       ( max_bbe );
    _p.set_DISPLAY_DEGREE    ( old_dd  );
    _p.set_STATS_FILE        ( old_stats_file_name , old_stats_file );
    _p.set_SOLUTION_FILE     ( old_sol_file );
    _p.set_LH_SEARCH         ( lh_p0 , lh_pi );
    _p.set_CACHE_SAVE_PERIOD ( old_csp );
    _p.set_L_CURVE_TARGET    ( old_lct );
    
    if ( use_delta_crit )
      _p.set_MIN_POLL_SIZE ( original_min_poll_size );
    
    _p.check ( false ,    // remove_history_file  = false
	       true  ,    // remove_solution_file = true
	       true    ); // remove_stats_file    = true

    // reset MADS stats from MULTI-MADS stats:
    _stats = multi_stats;
    
    // final cache save (overwrite=true):
    _ev_control.save_caches ( true );
    
#ifndef R_VERSION
    if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY) {
      std::ostringstream msg;
      msg << "end of run (" << stop_reason << ")";
      out << std::endl << NOMAD::close_block ( msg.str() ) << std::endl;
    }
#endif

    // multi-objective final displays:
    if ( _p.get_nb_obj() > 1 ) {

      if ( display_degree == NOMAD::FULL_DISPLAY )
	out.open_block ( "NOMAD final display" );

      display();

      if ( display_degree == NOMAD::FULL_DISPLAY )
	out.close_block();
    }

  } // end of the try block
  
  catch ( std::exception & e ) {
    
#ifdef USE_MPI
    if ( NOMAD::Slave::are_running() )
      NOMAD::Slave::stop_slaves ( out );
#endif

    throw NOMAD::Exception ( "Mads.cpp" , __LINE__ , e.what() );
  }
  
  // stop the slaves:
#ifdef USE_MPI
  if ( NOMAD::Slave::are_running() && stop_slaves_here )
    NOMAD::Slave::stop_slaves ( out );
#endif

  return stop_reason;
}

/*---------------------------------------------------------*/
/*                 one MADS iteration (private)            */
/*---------------------------------------------------------*/
void NOMAD::Mads::iteration ( bool                     & stop           ,
			      NOMAD::stop_type         & stop_reason    ,
			      NOMAD::success_type      & success        ,
			      const NOMAD::Eval_Point *& new_feas_inc   ,
			      const NOMAD::Eval_Point *& new_infeas_inc   )
{

  bool forbid_poll_size_stop = false;

  // force quit (by pressing ctrl-c):
  if ( !stop && NOMAD::Mads::_force_quit ) {
    stop        = true;
    stop_reason = NOMAD::CTRL_C;
    return;
  }

  // poll center selection:
  ( ( _p.get_opt_only_sgte() ) ?
    _sgte_barrier : _true_barrier ).select_poll_center ( success );

  // displays:
  const NOMAD::Display & out = _p.out();
  if ( out.get_iter_dd() == NOMAD::FULL_DISPLAY )
    out << std::endl
	<< NOMAD::open_block ( "MADS iteration "
			       + NOMAD::itos ( _stats.get_iterations() ) )
	<< std::endl;
  display_iteration_begin();

  // SEARCH:
  // -------
  search ( stop , stop_reason , success , new_feas_inc , new_infeas_inc );

  // POLL:
  // -----
  if ( success != NOMAD::FULL_SUCCESS )
    poll ( stop                  ,
	   stop_reason           ,
	   success               ,
	   new_feas_inc          ,
	   new_infeas_inc        ,
	   forbid_poll_size_stop   );

  // UPDATES:
  // --------
  int old_ell = NOMAD::Mesh::get_mesh_index();

  if ( !stop ) {

    // mesh update:
    NOMAD::Mesh::update ( success );

    // check the min mesh/poll sizes stopping criteria:
    _p.get_signature()->get_mesh().check_min_mesh_sizes ( _p.get_max_mesh_index()      ,
							  NOMAD::Mesh::get_mesh_index(),
							  stop                         ,
							  stop_reason );

    // if the Delta_k^p stopping criterion is met with integer variables,
    // the last set of directions must have a minimal coordinate of 1;
    // otherwise the stopping criterion is disabled at this iteration:
    if ( forbid_poll_size_stop && stop && stop_reason == NOMAD::DELTA_P_MIN_REACHED ) {
      stop = false;
      stop_reason = NOMAD::NO_STOP;
    }

    // display:
    if ( _p.out().get_iter_dd() == NOMAD::FULL_DISPLAY )
      _p.out() << std::endl << NOMAD::open_block ( "Mesh update" )
	       << "previous mesh index: " << old_ell << std::endl
	       << "new mesh index     : " << NOMAD::Mesh::get_mesh_index() << std::endl
	       << NOMAD::close_block() << std::endl;

    // periodic cache save (overwrite=false):
    if ( _p.get_cache_save_period() > 0 &&
	 _stats.get_iterations()%_p.get_cache_save_period() ==
	 _p.get_cache_save_period() - 1 )
      _ev_control.save_caches ( false );
  }

  // number of iterations:
  _stats.add_iteration();
  if ( !stop                       &&
       _p.get_max_iterations() > 0 &&
       _stats.get_iterations() >= _p.get_max_iterations() ) {
    stop        = true;
    stop_reason = NOMAD::MAX_ITER_REACHED;
  }

  // max cache memory:
  if ( !stop                           &&
       _p.get_max_cache_memory() > 0.0 &&
       _ev_control.get_cache().size_of() > 1048576*_p.get_max_cache_memory() ) {
    stop        = true;
    stop_reason = NOMAD::MAX_CACHE_MEMORY_REACHED;
  }

  // L_CURVE_TARGET stopping criterion:
  if ( _L_curve && !stop ) {
    int bbe = _stats.get_bb_eval();
    if ( success == NOMAD::FULL_SUCCESS ) {
      if ( new_feas_inc )
	_L_curve->insert ( bbe , new_feas_inc->get_f() );
    }
    else if ( success == NOMAD::UNSUCCESSFUL && _L_curve->check_stop ( bbe ) ) {
      stop = true;
      stop_reason = NOMAD::L_CURVE_TARGET_REACHED;
    }
  }

  // call the user-defined function Evaluator::update_iteration():
  if ( _p.get_user_calls_enabled() ) {
    bool stop_before = stop;
    _ev_control.get_evaluator()->update_iteration ( success        ,
						    _stats         ,
						    _ev_control    ,
						    _true_barrier  ,
						    _sgte_barrier  ,
						    *_pareto_front ,
						    stop             );
    if ( !stop_before && stop )
      stop_reason = NOMAD::USER_STOPPED;
  }

  // if the algorithms stops, we set the mesh index to the value
  // it had before the mesh update:
  if ( stop )
    Mesh::set_mesh_index ( old_ell );

  // displays at the end of an iteration:
  display_iteration_end ( stop           ,
			  stop_reason    ,
			  success        ,
			  new_feas_inc   ,
			  new_infeas_inc   );
  // displays:
  if ( out.get_iter_dd() == NOMAD::FULL_DISPLAY )
    out << std::endl
	<< NOMAD::close_block ( "end of iteration "
				+ NOMAD::itos ( _stats.get_iterations()-1 ) );
}

/*---------------------------------------------------------*/
/*                       the poll (private)                */
/*---------------------------------------------------------*/
void NOMAD::Mads::poll ( bool                     & stop                  ,
			 NOMAD::stop_type         & stop_reason           ,
			 NOMAD::success_type      & success               ,
			 const NOMAD::Eval_Point *& new_feas_inc          ,
			 const NOMAD::Eval_Point *& new_infeas_inc        ,
			 bool                     & forbid_poll_size_stop   )
{
  forbid_poll_size_stop = false;

  if ( stop )
    return;

  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_poll_dd();

  success            = NOMAD::UNSUCCESSFUL;
  new_feas_inc       = NULL;
  new_infeas_inc     = NULL;

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl << NOMAD::open_block ( "MADS poll" ) << std::endl;

  NOMAD::Eval_Point                         * pt;
  const NOMAD::Eval_Point                   * x;
  const NOMAD::Direction                    * dir;
  int                                         i , k;
  NOMAD::poll_center_type                     pc_type;
  int                                         n;
  int                                         m          = _p.get_bb_nb_outputs();
  const std::vector<NOMAD::bb_input_type>   & bbit       = _p.get_bb_input_type();
  int                                         offset     = 0;
  int                                         mesh_index = NOMAD::Mesh::get_mesh_index();
  const NOMAD::Barrier                      & barrier    = get_active_barrier();
  std::vector<NOMAD::Signature *>             signatures;
  NOMAD::Signature                          * cur_signature;
  std::list<NOMAD::Direction>::const_iterator it , end;

  // poll centers:
  const NOMAD::Eval_Point * poll_centers[2] , * poll_center;
  poll_centers[0] = barrier.get_poll_center();
  poll_centers[1] = (_p.use_sec_poll_center()) ?
    barrier.get_sec_poll_center() : NULL;

  if ( !poll_centers[0] && !poll_centers[1] )
    throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
			     "Mads::poll(): could not get a poll center" );

  // loop on the two poll centers:
  // -----------------------------
  NOMAD::poll_type i_pc = NOMAD::PRIMARY;
  poll_center           = poll_centers[NOMAD::PRIMARY];

  while ( true ) {

    if ( poll_center ) {

      if ( display_degree == NOMAD::FULL_DISPLAY ) {
	if ( i_pc == NOMAD::SECONDARY )
	  out << "secondary ";
	out << "poll center: ( ";
	poll_center->Point::display ( out, " ", 2, NOMAD::Point::get_display_limit() );
	out << " )" << std::endl;
      }

      // get the poll center's signature:
      cur_signature = poll_center->get_signature();
      signatures.push_back ( cur_signature );
	
      if ( !cur_signature )
	throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
	      "Mads::poll(): could not get the poll center's signature" );
      n = cur_signature->get_n();
      if ( n != poll_center->size() )
	throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
	      "Mads::poll(): the poll center has an incompatible signature" );
      
      // get directions from the signature:
      std::list<NOMAD::Direction> dirs;

      cur_signature->get_directions ( dirs                          ,
				      i_pc                          ,
				      *poll_center                  ,
				      _first_succ                   ,     // can be NULL
				      NOMAD::Mesh::get_mesh_index()   );

      if ( !stop && dirs.empty() ) {
	if ( display_degree == NOMAD::FULL_DISPLAY )
	  out << "Mads::poll(): could not get directions: stop"
	      << std::endl << NOMAD::close_block() << std::endl;
	stop        = true;
	stop_reason = NOMAD::MESH_PREC_REACHED;
	return;
      }

      // poll center type:
      pc_type = ( poll_center->is_feasible ( _p.get_h_min() ) ) ?
	NOMAD::FEASIBLE : NOMAD::INFEASIBLE;

      end = dirs.end();

      int nb_dirs = static_cast<int>(dirs.size());

      // displays:
      if ( display_degree == NOMAD::FULL_DISPLAY ) {
	out << std::endl
	    << NOMAD::open_block ( "poll directions (include mesh size parameter)" );
	k = offset;
	for ( it = dirs.begin() ; it != end ; ++it ) {
	  out << "dir ";
	  out.display_int_w ( k++ , nb_dirs );
	  out << " : " << *it << std::endl;
	}
	out.close_block();
      }

      NOMAD::Random_Pickup rp ( nb_dirs );

      // creation of the poll trial points:
      k = 0;
      for ( it = dirs.begin() ; it != end ; ++it ) {
	
	it->set_index ( offset + k );
	
	dir = &(*it); 
		  
	
	pt = new NOMAD::Eval_Point ( n , m );
		  
	// pt = poll_center + dir: with a particular case for binary variables
	// equal to 1 with dir=1: the variables are set to 0 (1+1=0 in binary):
	for ( i = 0 ; i < n ; ++i )
	{
	  (*pt)[i] =
	    ( bbit[i]==NOMAD::BINARY && (*dir)[i]==1.0 && (*poll_center)[i]==1.0 ) ?
	    0.0 : (*pt)[i] = (*poll_center)[i] + (*dir)[i];
	}
		  
	// we check that the new poll trial point is different than the poll center
	// (this happens when the mesh size becomes too small):
	if ( !stop && pt->Point::operator == ( *poll_center ) ) {
	  delete pt;
	  if ( display_degree == NOMAD::FULL_DISPLAY )
	    out << "Mads::poll(): could not generate poll trial points: stop"
		<< std::endl << NOMAD::close_block() << std::endl;
	  stop        = true;
	  stop_reason = NOMAD::MESH_PREC_REACHED;
	  return;
	}
	  
	pt->set_signature        ( cur_signature );
	pt->set_direction        ( dir           );
	pt->set_mesh_index       ( &mesh_index   );
	pt->set_poll_center_type ( pc_type       );
	  
	// random direction?
	if ( NOMAD::dir_is_random ( dir->get_type() ) )
	  pt->set_rand_eval_priority ( rp.pickup() );

	_ev_control.add_eval_point ( pt                      ,
				     display_degree          ,
				     _p.get_snap_to_bounds() ,
				     NOMAD::Double()         ,
				     NOMAD::Double()         ,
				     NOMAD::Double()         ,
				     NOMAD::Double()           );
	++k;
      }
      
      offset = nb_dirs;
	
      // display the re-ordered list of poll trial points:
      if ( display_degree == NOMAD::FULL_DISPLAY ) {
	const std::set<NOMAD::Priority_Eval_Point> & poll_pts
	  = _ev_control.get_eval_lop();
	out << std::endl
	    << NOMAD::open_block ( "re-ordered list of "
				   + NOMAD::itos ( poll_pts.size() )
				   + " poll trial points" );
	std::set<NOMAD::Priority_Eval_Point>::const_iterator
	  end2 = poll_pts.end() , it2;
	for ( it2 = poll_pts.begin() ; it2 != end2 ; ++it2 ) {
	  x =  it2->get_point();
	  x->display_tag ( out );
	  out << " : ( ";
	  x->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
	  out << " )";
	  if ( x->get_direction() )
	    out << " (dir " << x->get_direction()->get_index() << ")";
	  out << std::endl;
	}
	out.close_block();
      }
    }
   
    // loop increment:
    if ( i_pc == NOMAD::PRIMARY ) {
      i_pc        = NOMAD::SECONDARY;
      poll_center = poll_centers[NOMAD::SECONDARY];
    }
    else
      break;
  }

  _stats.add_poll_pts ( _ev_control.get_nb_eval_points() );

  // the directions are checked to satisfy a minimum
  // poll size with integer variables:
  check_directions ( forbid_poll_size_stop );

  // eval_list_of_points (poll):
  // ---------------------------
  new_feas_inc   = NULL;
  new_infeas_inc = NULL;

  _ev_control.eval_list_of_points ( NOMAD::POLL    ,
				    _true_barrier  ,
				    _sgte_barrier  ,
				    _pareto_front  ,
				    stop           ,
				    stop_reason    ,
				    new_feas_inc   ,
				    new_infeas_inc ,
				    success          );
  
  // extended poll for categorical variables:
  // ----------------------------------------
  if ( !stop                          &&
       _extended_poll                 &&
       success != NOMAD::FULL_SUCCESS &&
       _p.get_extended_poll_enabled()    ) {

    // display:
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl << NOMAD::open_block ( "MADS extended poll" ) << std::endl;

#ifdef USE_MPI
    // asynchronous mode: wait for the evaluations in progress:
    if ( _p.get_asynchronous() ) {
      std::list<const NOMAD::Eval_Point *> evaluated_pts;
      _ev_control.wait_for_evaluations ( NOMAD::ASYNCHRONOUS ,
					 _true_barrier       ,
					 _sgte_barrier       ,
					 _pareto_front       ,
					 stop                ,
					 stop_reason         ,
					 success             ,
					 evaluated_pts         );
    }
#endif

    // reset the extended poll object:
    _extended_poll->poll_reset();

    // call the user defined method changing the categorical variables
    // (this creates the list of extended poll points):
    _extended_poll->construct_extended_points ( *barrier.get_poll_center() );

    // add the signatures in use to the list of poll signatures:
    {
      const std::set<NOMAD::Signature_Element> &
	tmp = _extended_poll->get_poll_signatures();
      std::set<NOMAD::Signature_Element>::const_iterator it , end = tmp.end();
      for ( it = tmp.begin() ; it != end ; ++it )
	signatures.push_back ( it->get_signature() );
    }

    // execute the extended poll:
    int nb_ext_poll_pts;
    _extended_poll->run ( *this           ,
			  nb_ext_poll_pts ,
			  stop            ,
			  stop_reason     ,
			  success         ,
			  new_feas_inc    ,
			  new_infeas_inc    );

    // stats updates:
    _stats.add_ext_poll_pts ( nb_ext_poll_pts );
    if ( success == NOMAD::FULL_SUCCESS )
      _stats.add_ext_poll_succ();
    _stats.add_nb_ext_polls();

    // display:
    if ( display_degree == NOMAD::FULL_DISPLAY )
      out << std::endl << NOMAD::close_block ( "end of extended poll" ) << std::endl;
  }
  
  // stats updates:
  if ( success == NOMAD::FULL_SUCCESS )
    _stats.add_poll_success();
  _stats.add_nb_poll_searches();

  // success directions (feasible and infeasible):
  update_success_directions ( new_feas_inc   , signatures , true  );
  update_success_directions ( new_infeas_inc , signatures , false );

  // check the PEB constraints: if we have a new best infeasible
  // incumbent from another infeasible incumbent
  // ( active_barrier.check_PEB_constraints() ):
  if ( _p.get_barrier_type() == NOMAD::PEB_P &&
       new_infeas_inc                        &&
       new_infeas_inc->get_poll_center_type() == NOMAD::INFEASIBLE )
    ( ( _p.get_opt_only_sgte() ) ?
      _sgte_barrier : _true_barrier ).check_PEB_constraints
      ( *new_infeas_inc , display_degree==NOMAD::FULL_DISPLAY );
  
  // final display:
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << NOMAD::close_block ( "end of poll" );
}

/*----------------------------------------------------------------*/
/*             check the directions after the poll step           */
/*----------------------------------------------------------------*/
/*  ensures that the last set of poll directions is small enough  */
/*  with integer variables                                        */
/*----------------------------------------------------------------*/
void NOMAD::Mads::check_directions ( bool & forbid_poll_size_stop )
{
  if ( !_p.get_min_poll_size_defined() ) {

    NOMAD::Double        v , min;
    const NOMAD::Point * dir;
    int                  i , n;

    const NOMAD::Signature * signature;

    const std::set<NOMAD::Priority_Eval_Point> & poll_pts
      = _ev_control.get_eval_lop();
    std::set<NOMAD::Priority_Eval_Point>::const_iterator
      end = poll_pts.end() , it;
    for ( it = poll_pts.begin() ; it != end ; ++it ) {
	
      signature = it->get_point()->get_signature();

      if ( signature ) {
      
	dir = it->get_point()->get_direction();

	if ( dir ) {
	
	  n = dir->size();

	  if ( n == signature->get_n() ) {
	    
	    const std::vector<NOMAD::bb_input_type> & bbit
	      = signature->get_input_types();   
	    
	    for ( i = 0 ; i < n ; ++i ) {
	      if ( bbit[i] == NOMAD::INTEGER ) {
		v = (*dir)[i].abs();
		if ( v.is_defined() && v > 0.0 && ( !min.is_defined() || v < min ) )
		  min = v;
	      }
	    }
	  }
	}
      }
    }
    
    if ( min.is_defined() && min > 1.0 )
      forbid_poll_size_stop = true;
  }  
}

/*---------------------------------------------------------*/
/*    update of the success directions, after the poll     */
/*    (private)                                            */
/*---------------------------------------------------------*/
void NOMAD::Mads::update_success_directions
( const NOMAD::Eval_Point         * new_inc    ,
  std::vector<NOMAD::Signature *> & signatures ,
  bool                              feasible     ) const
{
  if ( new_inc && new_inc->get_direction() ) {

    const NOMAD::Direction * dir       = new_inc->get_direction();
    NOMAD::Signature       * signature = new_inc->get_signature();
  
    if ( !signature )
      throw NOMAD::Exception ( "Mads.cpp" , __LINE__ ,
	    "Mads::update_success_directions(): new incumbent has no signature" );

    if ( feasible )
      new_inc->get_signature()->set_feas_success_dir ( *dir );
    else
      new_inc->get_signature()->set_infeas_success_dir ( *dir );

  }
}

/*---------------------------------------------------------*/
/*                      the search (private)               */
/*---------------------------------------------------------*/
void NOMAD::Mads::search ( bool                     & stop           ,
			   NOMAD::stop_type         & stop_reason    ,
			   NOMAD::success_type      & success        ,
			   const NOMAD::Eval_Point *& new_feas_inc   ,
			   const NOMAD::Eval_Point *& new_infeas_inc   )
{
  int                    nb_search_pts;
  bool                   count_search;
  int                    mads_iteration  = _stats.get_iterations();
  const NOMAD::Display & out             = _p.out();
  NOMAD::dd_type         display_degree  = out.get_search_dd();
  NOMAD::success_type    last_it_success = success;
  success                                = NOMAD::UNSUCCESSFUL;

  // first display:
  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl << NOMAD::open_block ( "MADS search" );

  // 1. speculative search:
  if ( _p.get_speculative_search() ) {

    if ( new_feas_inc || new_infeas_inc ) {

      Speculative_Search ss ( _p );

      ss.search ( *this          ,
		  nb_search_pts  ,
		  stop           ,
		  stop_reason    ,
		  success        ,
		  count_search   ,
		  new_feas_inc   ,
		  new_infeas_inc   );

      if ( success == NOMAD::FULL_SUCCESS )
	_stats.add_spec_success();
      if ( count_search )
	_stats.add_nb_spec_searches();
      _stats.add_spec_pts ( nb_search_pts );
    }
  }

  // 2. user search:
  if ( success != NOMAD::FULL_SUCCESS && _user_search ) {

    // initial user search display:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      std::ostringstream oss;
      oss << NOMAD::USER_SEARCH;
      out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    }

    // the search:
    _user_search->search ( *this          ,
			   nb_search_pts  ,
			   stop           ,
			   stop_reason    ,
			   success        ,
			   count_search   ,
			   new_feas_inc   ,
			   new_infeas_inc   );

    // update stats:
    if ( success == NOMAD::FULL_SUCCESS )
      _stats.add_usr_srch_success();
    if ( count_search )
      _stats.add_nb_usr_searches();
    _stats.add_usr_srch_pts ( nb_search_pts );

    // final user search display:
    if ( display_degree == NOMAD::FULL_DISPLAY ) {
      std::ostringstream oss;
      oss << "end of " << NOMAD::USER_SEARCH << " (" << success << ")";
      out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
    }
  }

  // 3. cache search:
  if ( success != NOMAD::FULL_SUCCESS && _p.get_cache_search() ) {

    // the search:
    _cache_search->search ( *this          ,
			    nb_search_pts  ,
			    stop           ,
			    stop_reason    ,
			    success        ,
			    count_search   ,
			    new_feas_inc   ,
			    new_infeas_inc   );

    // update stats:
    if ( success == NOMAD::FULL_SUCCESS )
      _stats.add_CS_success();
    if ( count_search )
      _stats.add_nb_cache_searches();
    _stats.add_CS_pts ( nb_search_pts );
  }

  // 4. Model Searches (stats are updated inside the searches):
  if ( success != NOMAD::FULL_SUCCESS && _p.has_model_search() ) {

#ifdef USE_MPI
    // asynchronous mode: wait for the evaluations in progress:
    if ( _p.get_asynchronous() ) {
      std::list<const NOMAD::Eval_Point *> evaluated_pts;
      _ev_control.wait_for_evaluations ( NOMAD::ASYNCHRONOUS ,
					 _true_barrier       ,
					 _sgte_barrier       ,
					 _pareto_front       ,
					 stop                ,
					 stop_reason         ,
					 success             ,
					 evaluated_pts         );
    }
#endif

    // model search #1:
    _model_search1->search ( *this          ,
			     nb_search_pts  ,
			     stop           ,
			     stop_reason    ,
			     success        ,
			     count_search   ,
			     new_feas_inc   ,
			     new_infeas_inc   );

    // save the TGP model for the model ordering:
    if ( _p.get_model_search(1) == NOMAD::TGP_MODEL ) {
#ifdef USE_TGP
      _ev_control.set_last_TGP_model
	( static_cast<NOMAD::TGP_Model_Search *>(_model_search1)->get_model() );
#endif
    }
    // model search #2:
    if ( success != NOMAD::FULL_SUCCESS && _model_search2 ) {

#ifdef USE_MPI
      // asynchronous mode: wait for the evaluations in progress:
      if ( _p.get_asynchronous() ) {
	std::list<const NOMAD::Eval_Point *> evaluated_pts;
	_ev_control.wait_for_evaluations ( NOMAD::ASYNCHRONOUS ,
					   _true_barrier       ,
					   _sgte_barrier       ,
					   _pareto_front       ,
					   stop                ,
					   stop_reason         ,
					   success             ,
					   evaluated_pts         );
      }
#endif
      _model_search2->search ( *this          ,
			       nb_search_pts  ,
			       stop           ,
			       stop_reason    ,
			       success        ,
			       count_search   ,
			       new_feas_inc   ,
			       new_infeas_inc   );

      // save the TGP model for the model ordering:
      if ( _p.get_model_search(2) == NOMAD::TGP_MODEL ) {
#ifdef USE_TGP
	_ev_control.set_last_TGP_model
	  ( static_cast<NOMAD::TGP_Model_Search *>(_model_search2)->get_model() );
#endif
      }
    }
  }

  // 5. VNS search:
  if ( _p.get_VNS_search()                    &&
       success         != NOMAD::FULL_SUCCESS &&
       last_it_success == NOMAD::UNSUCCESSFUL &&
       NOMAD::Mesh::get_mesh_index() >= 0     &&
       _stats.get_iterations() > 0               ) {

    // check the VNS_trigger criterion:
    int bbe = _stats.get_bb_eval();
    if ( bbe==0 ||
	 _stats.get_VNS_bb_eval() / static_cast<float>(bbe) < _p.get_VNS_trigger() ) {

#ifdef USE_MPI
      // asynchronous mode: wait for the evaluations in progress:
      if ( _p.get_asynchronous() ) {
	std::list<const NOMAD::Eval_Point *> evaluated_pts;
	_ev_control.wait_for_evaluations ( NOMAD::ASYNCHRONOUS ,
					   _true_barrier       ,
					   _sgte_barrier       ,
					   _pareto_front       ,
					   stop                ,
					   stop_reason         ,
					   success             ,
					   evaluated_pts         );
      }
#endif

      _VNS_search->search ( *this          ,
			    nb_search_pts  ,
			    stop           ,
			    stop_reason    ,
			    success        ,
			    count_search   ,
			    new_feas_inc   ,
			    new_infeas_inc   );

      if ( success == NOMAD::FULL_SUCCESS )
	_stats.add_VNS_success();
      if ( count_search )
	_stats.add_nb_VNS_searches();
      _stats.add_VNS_pts ( nb_search_pts );
    }
  }

  // 6. Latin-Hypercube (LH) search:
  if ( success != NOMAD::FULL_SUCCESS && _p.get_LH_search_pi() > 0 ) {

    // for the first iteration: do not perform the
    // search if there was an initial LH search:
    if ( mads_iteration > 0 || _p.get_LH_search_p0() <= 0 ) {

      LH_Search lh ( _p , false , _flag_p1_active );

      lh.search ( *this          ,
		  nb_search_pts  ,
		  stop           ,
		  stop_reason    ,
		  success        ,
		  count_search   ,
		  new_feas_inc   ,
		  new_infeas_inc   );
      if ( success == NOMAD::FULL_SUCCESS )
	_stats.add_LH_success();
      if ( count_search )
	_stats.add_nb_LH_searches();
      _stats.add_LH_pts ( nb_search_pts );
    }
  }

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << NOMAD::close_block ( "end of search" );
}

/*---------------------------------------------------------*/
/*                       x0 eval (private)                 */
/*---------------------------------------------------------*/
void NOMAD::Mads::eval_x0 ( bool             & stop        ,
			    NOMAD::stop_type & stop_reason   )
{
  const std::vector<NOMAD::Point *> & x0s           = _p.get_x0s();
  const std::string                 & x0_cache_file = _p.get_x0_cache_file();
  if ( x0s.empty() && x0_cache_file.empty() )
    return;

  const NOMAD::Display    & out = _p.out();
  NOMAD::dd_type display_degree = out.get_gen_dd();

  if ( display_degree == NOMAD::FULL_DISPLAY )
    out << std::endl << NOMAD::open_block ( "starting point evaluation" );

  NOMAD::Eval_Point * pt;
  size_t              k;
  int                 m = _p.get_bb_nb_outputs();
  int                 n = _p.get_dimension();
  std::ostringstream  err;

  // x0s from vector Parameters::_x0s:
  // ---------------------------------
  size_t x0s_size = x0s.size();
  for ( k = 0 ; k < x0s_size ; ++k ) {

    // the current starting point has to be in dimension n:
    if ( x0s[k]->size() != n ) {
      err << "starting point ( " << *x0s[k] << " ) is not of dimension " << n;
      throw NOMAD::Exception ( "Mads.cpp" , __LINE__ , err.str() );
    }

    // creation of the Eval_Point:
    pt = new NOMAD::Eval_Point;
    pt->set           ( *x0s[k] , m        );
    pt->set_signature ( _p.get_signature() );

    _ev_control.add_eval_point ( pt              ,
				 display_degree  ,
				 false           ,
				 NOMAD::Double() ,
				 NOMAD::Double() ,
				 NOMAD::Double() ,
				 NOMAD::Double()    );
  }

  // x0 from a cache file:
  // ---------------------
  if ( !x0_cache_file.empty() ) {

    NOMAD::Cache            & cache = _ev_control.get_cache();
    const NOMAD::Eval_Point * x;

    // another cache file (this file won't be modified):
    if ( x0_cache_file != _p.get_cache_file() ) {

      NOMAD::Cache x0_cache ( out , ( _p.get_opt_only_sgte() ) ? NOMAD::SGTE  :
			                                         NOMAD::TRUTH   );
      std::string  file_name = _p.get_problem_dir() + x0_cache_file;

      if ( !x0_cache.load ( file_name , NULL , display_degree==NOMAD::FULL_DISPLAY ) ) {
 	err << "could not load (or create) the cache file " << file_name;
 	throw NOMAD::Exception ( "Mads.cpp" , __LINE__ , err.str() );
      }

      // we copy all the temporary cache points
      // into the list of points to be evaluated:
      x = x0_cache.begin();
      while ( x ) {

	pt = new NOMAD::Eval_Point;
	pt->set ( *x , m );

	if ( x->get_signature() )
	  pt->set_signature ( x->get_signature() );
	else if ( x->size() == n )
	  pt->set_signature ( _p.get_signature() );

	if ( pt->get_signature() )
	  _ev_control.add_eval_point ( pt              ,
				       display_degree  ,
				       false           ,
				       NOMAD::Double() ,
				       NOMAD::Double() ,
				       NOMAD::Double() ,
				       NOMAD::Double()   );
	else {
	  if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
	    out << std::endl << "Warning (" << "Mads.cpp" << ", " << __LINE__
		<< "): could not use the starting point " << *pt
		<< " (no signature)" << std::endl << std::endl;
	  delete pt;
	}

	x = x0_cache.next();
      }

      // insertion of this temporary cache in the algorithm's cache:
      cache.insert ( x0_cache );
    }

    // x0 cache file and the algorithm's cache file are the same:
    else {

      x = cache.begin();
      while ( x ) {
	pt = &NOMAD::Cache::get_modifiable_point ( *x );

	if ( x->get_signature() )
	  pt->set_signature ( x->get_signature() );
	else if ( x->size() == n )
	  pt->set_signature ( _p.get_signature() );

	if ( pt->get_signature() )
	  _ev_control.add_eval_point ( pt              ,
				       display_degree  ,
				       false           ,
				       NOMAD::Double() ,
				       NOMAD::Double() ,
				       NOMAD::Double() ,
				       NOMAD::Double()    );
	else {
	  if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
	    out << std::endl << "Warning (" << "Mads.cpp" << ", " << __LINE__
		<< "): could not use the starting point " << *pt
		<< "(no signature)" << std::endl;
	}
	x = cache.next();
      }
    }
  }

  // display of all starting points:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {

    const std::set<NOMAD::Priority_Eval_Point> & pts = _ev_control.get_eval_lop();

    // one starting point:
    if ( pts.size() == 1 ) {
      out << std::endl << "x0 eval point: ( ";
      pts.begin()->get_point()->Point::display ( out                               ,
						 " "                               ,
						 2                                 ,
						 NOMAD::Point::get_display_limit()   );
      out << " )" << std::endl;
    }

    // several starting points:
    else
      _ev_control.display_eval_lop ( NOMAD::X0_EVAL );
  }

  NOMAD::success_type       success;   
  const NOMAD::Eval_Point * new_feas_inc   = NULL;
  const NOMAD::Eval_Point * new_infeas_inc = NULL;

  // eval_list_of_points (x0):
  // -------------------------
  _ev_control.eval_list_of_points ( NOMAD::X0_EVAL ,
				    _true_barrier  ,
				    _sgte_barrier  ,
				    _pareto_front  ,
				    stop           ,
				    stop_reason    ,
				    new_feas_inc   ,
				    new_infeas_inc ,
				    success          );
  if ( !stop &&
       ( success == NOMAD::UNSUCCESSFUL      ||
	 (!new_feas_inc && !new_infeas_inc ) ||
	 ( _p.get_barrier_type() == NOMAD::EB &&
	   !get_active_barrier().get_best_feasible() ) ) ) {
    stop        = true;
    stop_reason = NOMAD::X0_FAIL;
  }

  // displays:
  display_iteration_end ( stop           ,
			  stop_reason    ,
			  success        ,
			  new_feas_inc   ,
			  new_infeas_inc   );

  // stop the algorithm if no iterations are allowed:
  if ( !stop && _p.get_max_iterations() == 0 ) {
    stop        = true;
    stop_reason = NOMAD::MAX_ITER_REACHED;
  }

  if ( display_degree == NOMAD::FULL_DISPLAY ) 
    out << std::endl << NOMAD::close_block ( "end of starting point evaluation" );
}

/*---------------------------------------------------------*/
/*                  display the Pareto front               */
/*---------------------------------------------------------*/
void NOMAD::Mads::display_pareto_front ( void ) const
{   
  if ( !_pareto_front )
    return;

  const std::string    & stats_file_name = _p.get_stats_file_name();
  const NOMAD::Display & out             = _p.out();
  NOMAD::dd_type         display_degree  = out.get_gen_dd();

  // loop on the Pareto points:
  if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
    out << std::endl << NOMAD::open_block ( "Pareto front" ) << std::endl;

  const NOMAD::Eval_Point * cur = _pareto_front->begin();
  while ( cur ) {

    if ( cur->is_eval_ok() && cur->is_feasible ( _p.get_h_min() ) ) {

      const std::list<int>           & index_obj = _p.get_index_obj();
      std::list<int>::const_iterator   it , end  = index_obj.end();
      const NOMAD::Point             & bbo       = cur->get_bb_outputs();
      int                              i         = 0;
      NOMAD::Point multi_obj ( static_cast<int>(index_obj.size()) );

      for ( it = index_obj.begin() ; it != end ; ++it )
	multi_obj[i++] = bbo[*it];

      if ( !stats_file_name.empty() )
	_ev_control.stats_file ( stats_file_name , cur , true , &multi_obj );

      if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY && !_p.get_display_stats().empty() )
	_ev_control.display_stats ( false                  ,
				    out                    ,
				    _p.get_display_stats() ,
				    cur                    ,
				    true                   ,
				    &multi_obj               );
    } 
    cur = _pareto_front->next();
  }

  if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
    out << NOMAD::close_block();

  // other stats:
  if ( display_degree == NOMAD::FULL_DISPLAY ) {

    out << std::endl << "number of pts : " << _pareto_front->size() << std::endl;

    NOMAD::Double delta , surf;
    _pareto_front->get_delta_surf ( delta , surf  ,
				    _p.get_multi_f_bounds() ); // f1_min, f1_max,
    // f2_min, f2_max
    out << "delta_j       : " << delta << std::endl
	<< "surf          : ";
    if ( surf.is_defined() )
      out << 100*surf << "%" << std::endl;
    else
      out << NOMAD::Double()
	  << " (define valid MULTI_F_BOUNDS values to access this output)"
	  << std::endl;
  }
  else if ( display_degree != NOMAD::NO_DISPLAY && display_degree != NOMAD::MINIMAL_DISPLAY)
    out << std::endl << "number of Pareto points: " << _pareto_front->size()
	<< std::endl;
}

/*---------------------------------------------------------*/
/*                        MADS display                     */
/*---------------------------------------------------------*/
void NOMAD::Mads::display ( const NOMAD::Display & out ) const
{
  NOMAD::dd_type display_degree = out.get_gen_dd();

  if ( !NOMAD::Slave::is_master() )
    return;

  // 0. no display:
  // --------------
  if ( display_degree == NOMAD::NO_DISPLAY || display_degree == NOMAD::MINIMAL_DISPLAY) {

    // there may be a pareto front to write as a stats file:
    if ( _pareto_front           &&
	 !_pareto_front->empty() &&
	 !_p.get_stats_file_name().empty() )
      display_pareto_front();
    
    return;
  }

  // incumbents:
  const NOMAD::Eval_Point * bf = get_best_feasible();
  const NOMAD::Eval_Point * bi = get_best_infeasible();

  // 1. detailed display:
  // --------------------
  if ( display_degree == NOMAD::FULL_DISPLAY ) {

    // cache:
    out << std::endl
	<< NOMAD::open_block ( "cache" )
	<< ( _p.get_opt_only_sgte() ?
	     _ev_control.get_sgte_cache() : _ev_control.get_cache() )
	<< NOMAD::close_block();

    // constraints:
    if ( _p.has_constraints() )
      out << std::endl
	  << NOMAD::open_block ( "constraints handling") << std::endl
	  << get_active_barrier()
	  << NOMAD::close_block();

    // Pareto front:
    if ( _pareto_front ) {
      if ( _pareto_front->empty() )
	out << std::endl << "Pareto front empty" << std::endl;
      else
	display_pareto_front();
    }

    // stats:
    out << std::endl
	<< NOMAD::open_block ( "stats" )
	<< _stats
	<< NOMAD::close_block();

    // model stats:
#ifdef DEBUG
    display_model_stats ( out );
#endif

    // miscellaneous:
    if ( !_pareto_front ) {
      out << std::endl
	  << NOMAD::open_block ( "miscellaneous" )
	  << "mesh index      : min="
	  << NOMAD::Mesh::get_min_mesh_index() << ", max="
	  << NOMAD::Mesh::get_max_mesh_index() << ", last="
	  << NOMAD::Mesh::get_mesh_index() << std::endl
	  << "best feas. sol  : ";

      if ( bf ) {
	out << "( ";
	bf->Point::display ( out , " " , -1 , -1 );
	out << " ) h=" << bf->get_h()
	    << " f="   << bf->get_f() << std::endl;
      }
      else
	out << "none" << std::endl;
      out << "best infeas. sol: ";
      if ( bi ) {
	out << "( ";
	bi->Point::display ( out , " " , -1 , -1 );
	out << " ) h=" << bi->get_h()
	    << " f="   << bi->get_f() << std::endl;
      }
      else
	out << "none" << std::endl;

      out.close_block();
    }
  }

  // 2. normal display:
  // ------------------
  else {

    // blackbox evaluations:
    out << std::endl
	<< "blackbox evaluations    : " << _stats.get_bb_eval() << std::endl;
      
    // output stats:
    if ( _stats.get_stat_sum().is_defined() )
      out << "stat sum                : " << _stats.get_stat_sum() << std::endl;
    if ( _stats.get_stat_avg().is_defined() )
      out << "stat avg                : " << _stats.get_stat_avg() << std::endl;
  
    // Pareto front (multi-objective optimization):
    if ( _pareto_front ) {
      out << "number of MADS runs     : " << _stats.get_mads_runs() << std::endl;
      if ( _pareto_front->empty() )
	out << "Pareto front            : empty" << std::endl;
      else
	display_pareto_front();
    }

    // single-objective optimization (display of best solutions):
    else {

      if ( !bf && !bi )
	out << "no solution" << std::endl;
      else {
	
	if ( bi ) {
	  out << "best infeasible solution: ( ";
	  bi->Point::display ( out , " " , -1 , -1 );
	  out << " ) h=" << bi->get_h()
	      << " f="  << bi->get_f() << std::endl;
	}

	out << "best feasible solution  : ";

	if ( bf ) {
	  out << "( ";
	  bf->Point::display ( out , " " , -1 , -1 );
	  out << " ) h=" << bf->get_h()
	      << " f="  << bf->get_f() << std::endl;
	}
	else
	  out << "no feasible solution has been found" << std::endl;
      }
    }
  }
}

/*---------------------------------------------------------*/
/*                    display model stats                  */
/*---------------------------------------------------------*/
void NOMAD::Mads::display_model_stats ( const NOMAD::Display & out ) const
{
  if ( _model_search1 )
    out << std::endl << NOMAD::open_block ( "model search #1 stats" )
	<< *_model_search1 << NOMAD::close_block();
  if ( _model_search2 )
    out << std::endl << NOMAD::open_block ( "model search #2 stats" )
	<< *_model_search2 << NOMAD::close_block();
  if ( _p.get_model_eval_sort() != NOMAD::NO_MODEL ) {
    out << std::endl << NOMAD::open_block ( "model ordering stats" );
    _ev_control.display_model_ordering_stats ( out );
    out << NOMAD::close_block();
  }
}

/*---------------------------------------------------------*/
/*  display mesh and poll sizes for a given signature and  */
/*  the current mesh index (private)                       */
/*---------------------------------------------------------*/
void NOMAD::Mads::display_deltas ( const NOMAD::Signature & s ) const
{
  NOMAD::Point delta_m , delta_p;
  s.get_mesh().get_delta_m ( delta_m , NOMAD::Mesh::get_mesh_index() );
  s.get_mesh().get_delta_p ( delta_p , NOMAD::Mesh::get_mesh_index() );
  _p.out() << "mesh size            : ( " << delta_m << " )" << std::endl
	   << "poll size            : ( " << delta_p << " )"
	   << std::endl;
}

/*-------------------------------------------------------*/
/*  displays at the beginning of an iteration (private)  */
/*-------------------------------------------------------*/
void NOMAD::Mads::display_iteration_begin ( void ) const
{
  const NOMAD::Display & out = _p.out();
  if ( out.get_iter_dd() != NOMAD::FULL_DISPLAY )
    return;

  // incumbents:
  const NOMAD::Eval_Point * bf = get_best_feasible();
  const NOMAD::Eval_Point * bi = get_best_infeasible();
  const NOMAD::Signature  * s1 = NULL;

  out << "blackbox evaluations : " << _stats.get_bb_eval() << std::endl;
#ifdef USE_MPI
  if ( _p.get_asynchronous() )
    out << "eval. in progress    : " << _ev_control.get_nb_eval_in_progress()
	<< std::endl;
#endif
  out << "best feas. solution  : ";
  if ( bf ) {
    out << "( ";
    bf->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
    out << " ) h=" << bf->get_h()
	<< " f="   << bf->get_f()
	<< std::endl;
  }
  else
    out << "none" << std::endl;
  out << "best infeas. solution: ";
  if ( bi ) {
    out << "( ";
    bi->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
    out << " ) h=" << bi->get_h()
	<< " f="   << bi->get_f()
	<< std::endl;
  }
  else
    out << "none" << std::endl;

  out << "poll center          : ";
  const NOMAD::Eval_Point * poll_center = get_active_barrier().get_poll_center();
  if ( poll_center ) {
    out << "( ";
    poll_center->Point::display ( out , " " , 2 , NOMAD::Point::get_display_limit() );
    out << " )" << std::endl;

    s1 = poll_center->get_signature();
    if (s1)
      display_deltas ( *s1 );
  }
  else
    out << "none" << std::endl;
  
  const NOMAD::Eval_Point * sec_poll_center
    = get_active_barrier().get_sec_poll_center();
  
  if ( sec_poll_center ) {
    out << "sec. poll center     : ( ";
    sec_poll_center->Point::display ( out                               ,
				      " "                               ,
				      2                                 ,
				      NOMAD::Point::get_display_limit()   );
    out << " )" << std::endl;
    const NOMAD::Signature * s2 = sec_poll_center->get_signature();
    if ( s2 && (!s1 || s1 != s2) )
      display_deltas ( *s2 );
  }

  out << "mesh index           : " << NOMAD::Mesh::get_mesh_index() << std::endl
      << "h_max                : "
      << get_active_barrier().get_h_max() << std::endl;
}

/*---------------------------------------------------------*/
/*       displays at the end of an iteration (private)     */
/*---------------------------------------------------------*/
void NOMAD::Mads::display_iteration_end
( bool                      stop           ,
  NOMAD::stop_type          stop_reason    ,
  NOMAD::success_type       success        ,
  const NOMAD::Eval_Point * new_feas_inc   ,
  const NOMAD::Eval_Point * new_infeas_inc   ) const
{
  const NOMAD::Display & out = _p.out();

  if ( out.get_iter_dd() != NOMAD::FULL_DISPLAY )
    return;

  out << std::endl
      << "terminate MADS       : ";
  out.display_yes_or_no ( stop );
  out << std::endl;
  if ( stop ) {
    out << "termination cause    : " << stop_reason;
    if ( stop_reason==NOMAD::X0_FAIL &&
	 !_flag_p1_active            &&
	 _p.has_EB_constraints()        )
      out << " (phase one will be performed)";
    out << std::endl;
  }
  out << "iteration status     : " << success << std::endl;
  out << "new feas. incumbent  : ";
  if ( new_feas_inc )
    out << *new_feas_inc;
  else
    out << "none" << std::endl;
  out << "new infeas. incumbent: ";
  if ( new_infeas_inc )
    out << *new_infeas_inc;
  else
    out << "none" << std::endl;
}
