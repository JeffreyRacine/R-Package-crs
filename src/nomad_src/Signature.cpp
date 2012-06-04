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
  \file   Signature.cpp
  \brief  Evaluation point signature (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Signature.hpp
*/
#include "Signature.hpp"
using namespace std;  //zhenghua
/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
#ifdef MEMORY_DEBUG
int NOMAD::Signature::_cardinality     = 0;
int NOMAD::Signature::_max_cardinality = 0;
#endif

/*--------------------------------------------------*/
/*                    constructor 1                 */
/*--------------------------------------------------*/
NOMAD::Signature::Signature
( int                                                     n                  ,
  const std::vector<NOMAD::bb_input_type>               & input_types        ,
  const NOMAD::Point                                    & initial_mesh_size  ,
  const NOMAD::Point                                    & min_mesh_size      ,
  const NOMAD::Point                                    & min_poll_size      ,
  const NOMAD::Point                                    & lb                 ,
  const NOMAD::Point                                    & ub                 ,
  const NOMAD::Point                                    & scaling            ,
  const NOMAD::Point                                    & fixed_variables    ,
  const std::vector<bool>                               & periodic_variables ,
  const std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> & var_groups )
  :  _mesh ( NULL  ) ,
     _std  ( false )
{
  init ( n                  ,
	 input_types        ,
	 initial_mesh_size  ,
	 min_mesh_size      ,
	 min_poll_size      ,
	 lb                 ,
	 ub                 ,
	 scaling            ,
	 fixed_variables    ,
	 periodic_variables ,
	 var_groups           );

#ifdef MEMORY_DEBUG
  ++NOMAD::Signature::_cardinality;
  if ( NOMAD::Signature::_cardinality > NOMAD::Signature::_max_cardinality )
    ++NOMAD::Signature::_max_cardinality;
#endif
}

/*--------------------------------------------------*/
/*                    constructor 2                 */
/*--------------------------------------------------*/
NOMAD::Signature::Signature
( int                                       n                  ,
  const std::vector<NOMAD::bb_input_type> & input_types        ,
  const NOMAD::Point                      & initial_mesh_size  ,
  const NOMAD::Point                      & lb                 ,
  const NOMAD::Point                      & ub                 ,
  const std::set<NOMAD::direction_type>   & direction_types    ,
  const std::set<NOMAD::direction_type>   & sec_poll_dir_types ,
  int                                       halton_seed        ,
  const NOMAD::Display                    & out                  )
  : _mesh ( NULL  ) ,
    _std  ( false )
{
  if ( static_cast<int> ( input_types.size() ) != n )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::Signature(): bad argument: input_types" );

  // automatic creation of groups of variables:
  // ------------------------------------------
  std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> var_groups;

  {
    std::set<int> vi_cbi;  // list of cont./bin./int. variables
    std::set<int> vi_cat;  // list of categorical variables

    for ( int i = 0 ; i < n ; ++i )
      if (  input_types[i] != NOMAD::CATEGORICAL )
	vi_cbi.insert(i);
      else
	vi_cat.insert(i);
    
    // creation of a group for cont./bin./int. variables:
    if ( !vi_cbi.empty() )
      var_groups.insert ( new NOMAD::Variable_Group ( vi_cbi             ,
						      direction_types    ,
						      sec_poll_dir_types ,
						      halton_seed++      ,
						      out                  ) );
    
    // creation of a group for categorical variables:
    if ( !vi_cat.empty() )
      var_groups.insert ( new NOMAD::Variable_Group ( vi_cat            ,
						      direction_types   ,
						      sec_poll_dir_types,
						      -1                , // no Halton seed
						      out                 ) );
  }

  // init:
  // -----
  init ( n                   ,
	 input_types         ,
	 initial_mesh_size   ,
	 NOMAD::Point()      ,
	 NOMAD::Point()      ,
	 lb                  ,
	 ub                  ,
	 NOMAD::Point()      ,
	 NOMAD::Point()      ,
	 std::vector<bool>() ,
	 var_groups            );

  // delete the temporary groups of variables:
  std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp>::const_iterator
    it , end = var_groups.end();

  for ( it = var_groups.begin() ; it != end ; ++it )
    delete *it;

#ifdef MEMORY_DEBUG
  ++NOMAD::Signature::_cardinality;
  if ( NOMAD::Signature::_cardinality > NOMAD::Signature::_max_cardinality )
    ++NOMAD::Signature::_max_cardinality;
#endif
}

/*--------------------------------------------------*/
/*                  copy constructor                */
/*--------------------------------------------------*/
NOMAD::Signature::Signature ( const NOMAD::Signature & s )
  :  _lb                 ( s._lb                     ) ,
     _ub                 ( s._ub                     ) ,
     _scaling            ( s._scaling                ) ,
     _fixed_variables    ( s._fixed_variables        ) ,
     _input_types        ( s._input_types            ) ,
     _all_continuous     ( s._all_continuous         ) ,
     _has_categorical    ( s._has_categorical        ) ,
     _periodic_variables ( s._periodic_variables     ) ,
     _mesh               ( new NOMAD::Mesh(*s._mesh) ) ,
     _std                ( false                     ) ,
     _feas_success_dir   ( s._feas_success_dir       ) ,
     _infeas_success_dir ( s._infeas_success_dir     )
{
  std::list<NOMAD::Variable_Group *>::const_iterator it , end = s._var_groups.end();
  for ( it = s._var_groups.begin() ; it != end ; ++it )
    _var_groups.push_back ( new NOMAD::Variable_Group(**it) ); 

#ifdef MEMORY_DEBUG
  ++NOMAD::Signature::_cardinality;
  if ( NOMAD::Signature::_cardinality > NOMAD::Signature::_max_cardinality )
    ++NOMAD::Signature::_max_cardinality;
#endif  
}

/*--------------------------------------------------*/
/*                    destructor                    */
/*--------------------------------------------------*/
NOMAD::Signature::~Signature ( void )
{
  clear(); 
#ifdef MEMORY_DEBUG
  --NOMAD::Signature::_cardinality;
#endif
}

/*--------------------------------------------------------------------------------*/
/*  clear (private method called by destructor and after an exception is thrown)  */
/*--------------------------------------------------------------------------------*/
void NOMAD::Signature::clear ( void )
{
  delete _mesh;
  _mesh            = NULL;
  _all_continuous  = true;
  _has_categorical = false;
  _std             = false;
  reset_var_groups();
  _feas_success_dir.clear();
  _infeas_success_dir.clear();
  _lb.clear();
  _ub.clear();
  _scaling.clear();
  _fixed_variables.clear();
  _input_types.clear();
  _periodic_variables.clear();
}

/*--------------------------------------------------*/
/*            feasible success direction            */
/*--------------------------------------------------*/
void NOMAD::Signature::set_feas_success_dir ( const NOMAD::Direction & d )
{
  if ( d.size() != static_cast<int>(_input_types.size()) )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::set_feas_success_dir(): bad direction" );
  _feas_success_dir = d;
}

/*--------------------------------------------------*/
/*            infeasible success direction          */
/*--------------------------------------------------*/
void NOMAD::Signature::set_infeas_success_dir ( const NOMAD::Direction & d )
{
  if ( d.size() != static_cast<int>(_input_types.size()) )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::set_infeas_success_dir(): bad direction" );
  _infeas_success_dir = d;
}

/*--------------------------------------------------*/
/*               initializations (private)          */
/*--------------------------------------------------*/
void NOMAD::Signature::init
( int                                       n                  ,
  const std::vector<NOMAD::bb_input_type> & input_types        ,
  const NOMAD::Point                      & initial_mesh_size  ,
  const NOMAD::Point                      & min_mesh_size      ,
  const NOMAD::Point                      & min_poll_size      ,
  const NOMAD::Point                      & lb                 ,
  const NOMAD::Point                      & ub                 ,
  const NOMAD::Point                      & scaling            ,
  const NOMAD::Point                      & fixed_variables    ,
  const std::vector<bool>                 & periodic_variables ,
  const std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> & var_groups )
{
  // reset success directions:
  _feas_success_dir.clear();
  _infeas_success_dir.clear();

  // check the dimension (n):
  if ( n <= 0 )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
				   "NOMAD::Signature::init(): bad argument: n" );

  // bounds:
  if ( lb.empty() )
    _lb.reset ( n );
  else if ( lb.size() == n )
    _lb = lb;
  else
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
				   "NOMAD::Signature::init(): bad argument: lb" );
  if ( ub.empty() )
    _ub.reset ( n );
  else if ( ub.size() == n )
    _ub = ub;
  else
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
				   "NOMAD::Signature::init(): bad argument: ub" );

  // scaling:
  if ( scaling.empty() )
    _scaling.reset ( n );
  else if ( scaling.size() == n )
    _scaling = scaling;
  else
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
				   "NOMAD::Signature::init(): bad argument: scaling" );
  int i;
  for ( i = 0 ; i < n ; ++i )
    if ( _scaling[i].is_defined() && _scaling[i] == 0.0 )
      throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	    "NOMAD::Signature::init(): bad argument: scaling (zero value)" );

  // fixed variables:
  if ( fixed_variables.empty() )
    _fixed_variables.reset ( n );
  else if ( fixed_variables.size() == n )
    _fixed_variables = fixed_variables;
  else
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::init(): bad argument: fixed_variables" );

  // periodic variables:
  _periodic_variables = periodic_variables;
  if ( !_periodic_variables.empty() ) {
    
    if ( static_cast<int>(periodic_variables.size()) != n )
      throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	    "NOMAD::Signature::init(): bad argument: periodic_variables" );

    for ( i = 0 ; i < n ; ++i )
      if ( _periodic_variables[i] && ( !_lb[i].is_defined() || !_ub[i].is_defined() ) )
	throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	      "NOMAD::Signature::init(): incompatible periodic variables" );
  }

  // input types:
  if ( static_cast<int>(input_types.size()) != n )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::init(): bad argument: input_types" );
  _input_types     = input_types;
  _all_continuous  = true;
  _has_categorical = false;

  for ( i = 0 ; i < n ; ++i ) {
    if ( _input_types[i] == NOMAD::CATEGORICAL ) {
      _has_categorical = true;
      _all_continuous  = false;
      break;
    }
    if ( _input_types[i] != NOMAD::CONTINUOUS ) {
      _all_continuous = false;
      if ( _has_categorical )
	break;
    }
  }

  // variable groups:
  reset_var_groups();

  std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp>::const_iterator
    end = var_groups.end() , it;
  for ( it = var_groups.begin() ; it != end ; ++it ) {

    if ( !(*it)->check ( _fixed_variables , input_types , NULL ) )
      throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	    "NOMAD::Signature::init(): incompatible variable group" );

    _var_groups.push_back ( new NOMAD::Variable_Group (**it) );
  }

  // mesh:
  if ( initial_mesh_size.size() != n                             ||
       (min_mesh_size.is_defined() && min_mesh_size.size() != n) ||
       (min_poll_size.is_defined() && min_poll_size.size() != n)    )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::init(): mesh arguments with different sizes" );
  delete _mesh;
  
  try {
    _mesh = new NOMAD::Mesh ( initial_mesh_size , min_mesh_size , min_poll_size );
  }
  catch ( NOMAD::Exception & e ) {
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ ,
					      *this , e.what() );
  }

}

/*--------------------------------------------------*/
/*                       reset                      */
/*--------------------------------------------------*/
void NOMAD::Signature::reset
( int                                       n                        ,
  const std::vector<NOMAD::bb_input_type> & input_types              ,
  const NOMAD::Point                      & initial_mesh_size        ,
  const NOMAD::Point                      & min_mesh_size            ,
  const NOMAD::Point                      & min_poll_size            ,
  const NOMAD::Point                      & lb                       ,
  const NOMAD::Point                      & ub                       ,
  const NOMAD::Point                      & scaling                  ,
  const NOMAD::Point                      & fixed_variables          ,
  const std::vector<bool>                 & periodic_variables       ,
  const std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> & var_groups   )
{
  reset_var_groups();
  init ( n                  ,
	 input_types        ,
	 initial_mesh_size  ,
	 min_mesh_size      ,
	 min_poll_size      ,
	 lb                 ,
	 ub                 ,
	 scaling            ,
	 fixed_variables    ,
	 periodic_variables ,
	 var_groups           );
}

/*--------------------------------------------------*/
/*                reset_var_groups (private)        */
/*--------------------------------------------------*/
void NOMAD::Signature::reset_var_groups ( void )
{
  std::list<NOMAD::Variable_Group *>::const_iterator
    end = _var_groups.end() , it;
  for ( it = _var_groups.begin() ; it != end ; ++it )
    delete *it;
  _var_groups.clear();
}

/*----------------------------------------------------------*/
/*           check the compatibility of a point             */
/*   . we only check the number of variables                */
/*   . other stuff (fixed variables, binary variables,...)  */
/*     will be checked by Eval_Point::check()               */
/*----------------------------------------------------------*/
bool NOMAD::Signature::is_compatible ( const NOMAD::Point & x ) const
{
  if ( get_n() != x.size() || !_mesh || _var_groups.empty() )
    return false;
  return true;
}

/*-----------------------------------------------------*/
/*     compute the directions (they include delta_m)   */
/*-----------------------------------------------------*/
void NOMAD::Signature::get_directions ( std::list<NOMAD::Direction> & dirs          ,
					NOMAD::poll_type              poll          ,
					const NOMAD::Point          & poll_center   ,
					const NOMAD::Point          * first_success ,  // can be NULL
					int                           mesh_index      ) const
{
  NOMAD::Direction                          * pd;
  int                                         i;
  std::list<NOMAD::Direction>::const_iterator it_dir , end_dir;
  std::set<int>::const_iterator               it_vi  , end_vi;

  // get delta_m (mesh size parameter):
  int          n = get_n();
  NOMAD::Point delta_m (n);
	NOMAD::Point delta_p(n);	
  _mesh->get_delta_m ( delta_m , mesh_index );
  _mesh->get_delta_p ( delta_p , mesh_index );	

  // loop on variable groups:
  std::list<NOMAD::Variable_Group*>::const_iterator end_vg = _var_groups.end() , it_vg;
  for ( it_vg = _var_groups.begin() ; it_vg != end_vg ; ++it_vg ) {

    const std::set<int> & var_indexes = (*it_vg)->get_var_indexes();

    // get the directions for the current group of variables:
    std::list<NOMAD::Direction> dirs_nc;
    (*it_vg)->get_directions ( dirs_nc             ,
			       poll                ,
			       poll_center         ,
			       first_success       ,
			       mesh_index          ,
			       _feas_success_dir   ,
			       _infeas_success_dir   );

    // scale with delta_m and resize the directions to size n;
    // also round integer and binary variables:
    end_dir = dirs_nc.end();
    for ( it_dir = dirs_nc.begin() ; it_dir != end_dir ; ++it_dir ) {
      dirs.push_back ( NOMAD::Direction ( n , 0.0 , it_dir->get_type() ) );
      pd = &(*(--dirs.end()));

      end_vi = var_indexes.end();
      i      = 0;
      for ( it_vi = var_indexes.begin() ; it_vi != end_vi ; ++it_vi ) {

	(*pd)[*it_vi] = delta_m[*it_vi] * (*it_dir)[i++];

	// integer variables:		    
	  if ( _input_types[*it_vi] == NOMAD::INTEGER ) {
		  if ( (*pd)[*it_vi] >= delta_p[*it_vi]/3.0 )
			  (*pd)[*it_vi] = ceil  ( (*pd)[*it_vi].value() );
		  else if ( (*pd)[*it_vi] <= -delta_p[*it_vi]/3.0 )
			  (*pd)[*it_vi] = floor ( (*pd)[*it_vi].value() );
		  else	
		  {
			  double x=(*pd)[*it_vi].value();
			  (*pd)[*it_vi] = (x>0)? floor(x+0.5): ceil(x-0.5);
		  }
		}		  
		  
	// binary variables:
	else if ( _input_types[*it_vi] == NOMAD::BINARY ) {
	  if ( (*pd)[*it_vi] != 0.0 )
	    (*pd)[*it_vi] = 1.0;
	}

	// categorical variables: set direction=0:
	else if ( _input_types[*it_vi] == NOMAD::CATEGORICAL )
	  (*pd)[*it_vi] = 0.0;
      }
    }
  }
}

/*----------------------------------------------------------------*/
/*  get just one direction for a given mesh (used by VNS search)  */
/*----------------------------------------------------------------*/
void NOMAD::Signature::get_one_direction ( NOMAD::Direction & dir          ,
					   int                mesh_index   ,
					   int                halton_index   ) const
{
	int                           i;
	std::set<int>::const_iterator it_vi  , end_vi;
	
	// get delta_m (mesh size parameter):
	int          n = get_n();
	NOMAD::Point delta_m (n);
	NOMAD::Point delta_p(n);	
	_mesh->get_delta_m ( delta_m , mesh_index );
	_mesh->get_delta_p ( delta_p , mesh_index );	
	
	dir.reset    ( n , 0.0 );
	dir.set_type ( NOMAD::UNDEFINED_DIRECTION );
	
	// loop on variable groups:
	std::list<NOMAD::Variable_Group*>::const_iterator end_vg = _var_groups.end() , it_vg;
	for ( it_vg = _var_groups.begin() ; it_vg != end_vg ; ++it_vg ) {
		
		const std::set<int> & var_indexes = (*it_vg)->get_var_indexes();
		
		// get the direction for the current group of variables:
		NOMAD::Direction dir_nc;
		
		if ( (*it_vg)->get_one_direction ( dir_nc , mesh_index , halton_index++ ) ) {
			
			// scale with delta_m and round integer and binary variables:
			end_vi = var_indexes.end();
			i      = 0;
			for ( it_vi = var_indexes.begin() ; it_vi != end_vi ; ++it_vi ) {
				
				dir[*it_vi] = delta_m[*it_vi] * dir_nc[i++];
				
				// integer variables:		  	  
				if ( _input_types[*it_vi] == NOMAD::INTEGER ) {
					if ( dir[*it_vi] >= delta_p[*it_vi]/3.0 )
						dir[*it_vi] = ceil  ( dir[*it_vi].value() );
					else if ( dir [*it_vi] <= -delta_p[*it_vi]/3.0 )
						dir[*it_vi] = floor ( dir[*it_vi].value() );
					else
					{
						double x=dir[*it_vi].value();
						dir[*it_vi] = (x>0)? floor(x+0.5): ceil(x-0.5);
					}
				}
				
				// binary variables:
				else if ( _input_types[*it_vi] == NOMAD::BINARY ) {
					if ( dir[*it_vi] != 0.0 )
						dir[*it_vi] = 1.0;
				}
				
				// categorical variables: set direction=0:
				else if ( _input_types[*it_vi] == NOMAD::CATEGORICAL )
					dir[*it_vi] = 0.0;
			}
		}
	}
}

/*----------------------------------*/
/*             scaling              */
/*    (done before an evaluation)   */
/*----------------------------------*/
void NOMAD::Signature::scale ( NOMAD::Point & x )
{
  int n = get_n();
  if ( n != x.size() )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::scale(x): x.size() != signature.size()" );
  NOMAD::Double sci;
  for ( int i = 0 ; i < n ; ++i ) {
    sci = _scaling[i];
    if ( sci.is_defined() )
      x[i] *= sci;
  }
}

/*----------------------------------*/
/*             unscaling            */
/*    (done after an evaluation)    */
/*----------------------------------*/
void NOMAD::Signature::unscale ( NOMAD::Point & x )
{
  int n = get_n();
  if ( n != x.size() )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::unscale(x): x.size() != signature.size()" );
  NOMAD::Double sci;
  for ( int i = 0 ; i < n ; ++i ) {
    sci = _scaling[i];
    if ( sci.is_defined() )
      x[i] /= sci;
  }
}

/*-----------------------------------------------------------------------*/
/*                             snap to bounds                            */
/*-----------------------------------------------------------------------*/
/*  . returns true if x has been modified                                */
/*  . supposes that treat_periodic_variables() has already been invoked  */
/*    (if periodic variables have been treated, then bounds are          */
/*     satisfied and there is no need to snap anymore                    */
/*-----------------------------------------------------------------------*/
bool NOMAD::Signature::snap_to_bounds ( NOMAD::Point & x , NOMAD::Direction * direction )
{
  int n = get_n();
  if ( n != x.size() )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
	  "NOMAD::Signature::snap_to_bounds(x): x.size() != signature.size()" );

  bool modified        = false;
  bool no_periodic_var = _periodic_variables.empty();

  for ( int i = 0 ; i < n ; ++i )
    
    if ( no_periodic_var || !_periodic_variables[i] ) {

      const NOMAD::Double & ubi = _ub[i];
      NOMAD::Double       & xi  = x[i];
      if ( ubi.is_defined() && xi > ubi ) {
	if (direction)
	  (*direction)[i] += ubi - xi;
	xi       = ubi;
	modified = true;
      }
      const NOMAD::Double & lbi = _lb[i];
      if ( lbi.is_defined() && xi < lbi ) {
	if (direction)
	  (*direction)[i] += (lbi - xi);
	xi       = lbi;
	modified = true;
      }
    }
  
  return modified;
}

/*--------------------------------------------------*/
/*             treat the periodic variables         */
/*         returns true if x has been modified      */
/*--------------------------------------------------*/
bool NOMAD::Signature::treat_periodic_variables ( NOMAD::Point            & x       ,
						  const NOMAD::Direction *  old_dir ,
						  NOMAD::Direction       *& new_dir   )
{
  if ( _periodic_variables.empty() )
    return false;

  int n = get_n();
  if ( n != x.size() )
    throw NOMAD::Signature::Signature_Error ( "Signature.cpp" , __LINE__ , *this ,
    "NOMAD::Signature::treat_periodic_variables(x): x.size() != signature.size()" );

  new_dir = ( old_dir ) ? new NOMAD::Direction (*old_dir) : NULL;

  bool modified = false;

  for ( int i = 0 ; i < n ; ++i ) {

    NOMAD::bb_input_type bbit = _input_types[i];
    
    if ( _periodic_variables[i] && !_fixed_variables[i].is_defined() &&
	 ( bbit == NOMAD::CONTINUOUS || bbit == NOMAD::INTEGER ) ) {

      const NOMAD::Double & ubi = _ub[i];
      const NOMAD::Double & lbi = _lb[i];
      NOMAD::Double       & xi  = x[i];
  
      bool chk = false;

      NOMAD::Double new_x = xi;
      while ( new_x > ubi ) {
	new_x += lbi - ubi;
	chk    = true;
      }

      if ( !chk ) {
	while ( new_x < lbi ) {
	  new_x += ubi - lbi;
	  chk    = true;
	}
      }

      if ( chk ) {

	if ( bbit == NOMAD::INTEGER )
	  new_x.round();
	
	if (new_dir)
	  (*new_dir)[i] += new_x - xi;

	x[i]     = new_x;
	modified = true;
      }
    }
  }
  
  return modified;
}

/*--------------------------------------------------*/
/*                      display                     */
/*--------------------------------------------------*/
void NOMAD::Signature::display ( const NOMAD::Display & out ) const
{
  if ( _std )
    out << "(standard signature)" << std::endl;

  // dimension:
  out << "n                 : "  << get_n() << std::endl;

  // bounds:
  out << "lb                : ";
  if ( _lb.is_defined() )
    out << "( " << _lb << ")";
  else
    out << "none";
  out << std::endl;
  out << "ub                : ";
  if ( _ub.is_defined() )
    out << "( " << _ub << ")";
  else
    out << "none";
  out << std::endl;

  // scaling:
  out << "scaling           : ";
  if ( _scaling.is_defined() )
    out << "( " << _scaling << ")";
  else
    out << "none";
  out << std::endl;

  // fixed variables:
  out << "fixed variables   : ";
  if ( _fixed_variables.is_defined() )
    out << "( " << _fixed_variables << ")";
  else
    out << "none";
  out << std::endl;
  
  // input types:
  out << "input types       : (" << _input_types << " )" << std::endl;

  // periodic variables:
  out << "periodic variables: ";
  if ( _periodic_variables.empty() )
    out << "none" << std::endl;
  else {
    size_t pvs = _periodic_variables.size();
    out << "{";
    for ( size_t k = 0 ; k < pvs ; ++k )
      out << _periodic_variables[k] << " ";
    out << "}" << std::endl;
  }
  
  // success directions:
  out << "feas. succ. dir.  : ";
  if ( _feas_success_dir.is_defined() )
    out << _feas_success_dir << std::endl;
  else
    out << "none";
  out << std::endl
      << "infeas. succ. dir.: ";
  if ( _infeas_success_dir.is_defined() )
    out << _infeas_success_dir;
  else
    out << "none";
  out << std::endl;
  
  // variable groups:
  out << NOMAD::open_block ( "variable groups" );
  int i = 0;
  std::list<NOMAD::Variable_Group *>::const_iterator end = _var_groups.end() , it;
  for ( it = _var_groups.begin() ; it != end ; ++it ) {
    out << NOMAD::open_block ( "group #" + NOMAD::itos(i++) )
	<< **it << NOMAD::close_block();
  }
  out.close_block();

  // mesh:
  out << NOMAD::open_block ( "mesh" ) << *_mesh
      << NOMAD::close_block();
}

/*-------------------------------------------*/
/*           comparison operator '<'         */
/*-------------------------------------------*/
/*  (success directions are not considered)  */
/*-------------------------------------------*/
bool NOMAD::Signature::operator < ( const NOMAD::Signature & s ) const
{
  if ( this == &s )
    return false;

  // standard signature: not checked: standard and non-standard signatures can be ==
  // ( this is tested in Parameters::check() )

  // dimension:
  // ----------
  int  n = _lb.size();
  int sn = s._lb.size();
  
  if ( n < sn )
    return true;
  if ( sn < n )
    return false;

  // variable groups:
  // ----------------
  int nvg1 = _var_groups.size();
  int nvg2 = s._var_groups.size();
  if ( nvg1 != nvg2 )
    return (nvg1 < nvg2);

  std::list<NOMAD::Variable_Group*>::const_iterator
    it1 =   _var_groups.begin() ,
    it2 = s._var_groups.begin() ,
    end =   _var_groups.end();
  
  while ( it1 != end ) {
    if ( **it1 < **it2 )
      return true;
    if ( **it2 < **it1 )
      return false;
    ++it1;
    ++it2;
  }

  // first check on the periodic variables:
  // --------------------------------------
  bool p1_empty = _periodic_variables.empty();
  bool p2_empty = s._periodic_variables.empty();
  if ( p1_empty != p2_empty )
    return p1_empty; 

  // first check on the mesh:
  // ------------------------
  const NOMAD::Point & delta_m_0     = _mesh->get_initial_mesh_size();
  const NOMAD::Point & delta_m_min   = _mesh->get_min_mesh_size();
  const NOMAD::Point & delta_p_min   = _mesh->get_min_poll_size();
  bool                 chkm          = delta_m_min.is_defined();
  bool                 chkp          = delta_p_min.is_defined();

  const NOMAD::Point & s_delta_m_0   = s._mesh->get_initial_mesh_size();
  const NOMAD::Point & s_delta_m_min = s._mesh->get_min_mesh_size();
  const NOMAD::Point & s_delta_p_min = s._mesh->get_min_poll_size();
  bool                 s_chkm        = s_delta_m_min.is_defined();
  bool                 s_chkp        = s_delta_p_min.is_defined();

  if ( _mesh != s._mesh ) {
    if ( chkm != s_chkm )
      return !chkm;
    if ( chkp != s_chkp )
      return !chkp;
  }

  /*---------------------------*/
  /*  loop on all coordinates  */
  /*---------------------------*/
  for ( int i = 0 ; i < n ; ++i ) {
    
    // input types:
    // ------------
    if ( _input_types[i] < s._input_types[i] )
      return true;
    if ( s._input_types[i] < _input_types[i] )
      return false;
    
    // bounds:
    // -------
    if ( _lb[i].comp_with_undef ( s._lb[i] ) )
      return true;
    if ( s._lb[i].comp_with_undef ( _lb[i] ) )
      return false;
    if ( _ub[i].comp_with_undef ( s._ub[i] ) )
      return true;
    if ( s._ub[i].comp_with_undef ( _ub[i] ) )
      return false;

    // scaling:
    // --------
    if ( _scaling[i].comp_with_undef ( s._scaling[i] ) )
      return true;
    if ( s._scaling[i].comp_with_undef ( _scaling[i] ) )
      return false;

    // fixed variables:
    // ----------------
    if ( _fixed_variables[i].comp_with_undef ( s._fixed_variables[i] ) )
      return true;
    if ( s._fixed_variables[i].comp_with_undef ( _fixed_variables[i] ) )
      return false;

    // periodic variables:
    // -------------------
    if ( !p1_empty && _periodic_variables[i] != s._periodic_variables[i] )
      return _periodic_variables[i];

    // mesh:
    // -----
    if ( _mesh != s._mesh ) {
      if ( delta_m_0[i].comp_with_undef ( s_delta_m_0[i] ) )
	return true;
      if ( s_delta_m_0[i].comp_with_undef ( delta_m_0[i] ) )
	return false;
      if ( chkm ) {
	if ( delta_m_min[i].comp_with_undef ( s_delta_m_min[i] ) )
	  return true;
	if ( s_delta_m_min[i].comp_with_undef ( delta_m_min[i] ) )
	  return false;
      }
      if ( chkp ) {
	if ( delta_p_min[i].comp_with_undef ( s_delta_p_min[i] ) )
	  return true;
	if ( s_delta_p_min[i].comp_with_undef ( delta_p_min[i] ) )
	  return false;
      }
    }
  }

  // both signatures are equal:
  return false;
}
