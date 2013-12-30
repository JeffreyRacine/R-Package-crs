/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.6.2        */
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
  \file   Mesh.cpp
  \brief  Class for the MADS mesh (implementation)
  \author Sebastien Le Digabel
  \date   2010-04-06
  \see    Mesh.hpp
*/
#include "Mesh.hpp"
using namespace std; // zhenghua

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
double NOMAD::Mesh::_mesh_update_basis        = 4.0;
int    NOMAD::Mesh::_mesh_coarsening_exponent =  1;
int    NOMAD::Mesh::_mesh_refining_exponent   = -1;
int    NOMAD::Mesh::_initial_mesh_index       =  0;
int    NOMAD::Mesh::_mesh_index               =  0;
int    NOMAD::Mesh::_min_mesh_index           =  0;
int    NOMAD::Mesh::_max_mesh_index           =  0;
int    NOMAD::Mesh::_max_halton_index         = -1;

/*-----------------------------------------------------------*/
/*                          constructor                      */
/*-----------------------------------------------------------*/
NOMAD::Mesh::Mesh ( const NOMAD::Point & delta_m_0   ,
		    const NOMAD::Point & delta_m_min ,
		    const NOMAD::Point & delta_p_min   )
 : _initial_mesh_size ( delta_m_0   ) ,
   _min_mesh_size     ( delta_m_min ) ,
   _min_poll_size     ( delta_p_min )
{
  bool chkp  = delta_p_min.is_defined();
  bool chkm  = delta_m_min.is_defined();
  int  k , n = delta_m_0.size();

  if ( !delta_m_0.is_complete() )
    throw NOMAD::Exception (  "Mesh.cpp" , __LINE__ ,
	  "NOMAD::Mesh::Mesh(): delta_m_0 has undefined values" );

  if ( chkp && delta_p_min.size() != n )
    throw NOMAD::Exception ( "Mesh.cpp" , __LINE__ ,
                 "NOMAD::Mesh::Mesh(): delta_m_0 and delta_p_min have different sizes" );

  if ( chkm && delta_m_min.size() != n )
    throw NOMAD::Exception ( "Mesh.cpp" , __LINE__ ,
                 "NOMAD::Mesh::Mesh(): delta_m_0 and delta_m_min have different sizes" );
  
  // we check that delta_p_min <= delta_p_0 and that delta_m_min <= delta_m_0:
  if ( chkp || chkm ) {

    NOMAD::Point * dp0 = NULL , * dm0 = NULL;

    if ( chkp ) {
      dp0 = new NOMAD::Point ( n );
      get_delta_p ( *dp0 , NOMAD::Mesh::_initial_mesh_index );
    }

    if ( chkm ) {
      dm0 = new NOMAD::Point ( n );
      get_delta_m ( *dm0 , NOMAD::Mesh::_initial_mesh_index );
    }

    std::string error;

    for ( k = 0 ; k < n ; ++k ) {
  
      if ( dp0                            &&
	   _min_poll_size[k].is_defined() &&
	   (*dp0)[k] < _min_poll_size[k]     ) {
 	error = "NOMAD::Mesh::Mesh(): delta_p_0 < delta_p_min";
	break;
       }

      if ( dm0                            &&
	   _min_mesh_size[k].is_defined() &&
	   (*dm0)[k] < _min_mesh_size[k]     ) {
	error = "NOMAD::Mesh::Mesh(): delta_m_0 < delta_m_min";
	break;
      }

    }

    delete dp0;
    delete dm0;

    if ( !error.empty() )
      throw NOMAD::Exception ( "Mesh.cpp" , __LINE__ , error );
  }
}

/*-----------------------------------------------------------*/
/*          Mesh static members initialization (static)      */
/*-----------------------------------------------------------*/
void NOMAD::Mesh::init ( double mesh_update_basis        ,
			 int    mesh_coarsening_exponent ,
			 int    mesh_refining_exponent   ,
			 int    initial_mesh_index         )
{
	NOMAD::Mesh::_mesh_update_basis        = mesh_update_basis;
	NOMAD::Mesh::_mesh_coarsening_exponent = mesh_coarsening_exponent;
	NOMAD::Mesh::_mesh_refining_exponent   = mesh_refining_exponent;
	NOMAD::Mesh::_mesh_index           =
	NOMAD::Mesh::_initial_mesh_index =
    NOMAD::Mesh::_min_mesh_index     =
    NOMAD::Mesh::_max_mesh_index     = initial_mesh_index;
	NOMAD::Mesh::_max_halton_index     = -1;
}

/*-----------------------------------------------------------*/
/*                    update the mesh (static)               */
/*-----------------------------------------------------------*/
void NOMAD::Mesh::update ( NOMAD::success_type success , int & mesh_index )
{
	// defaults:
	//  full success: lk = lk - 1
	//  failure     : lk = lk + 1
	
	if ( success == NOMAD::FULL_SUCCESS ) {
		mesh_index -= NOMAD::Mesh::_mesh_coarsening_exponent;
		if ( mesh_index < -NOMAD::L_LIMITS )
			mesh_index = -NOMAD::L_LIMITS;
	}	
	else if ( success == NOMAD::UNSUCCESSFUL )
		mesh_index -= NOMAD::Mesh::_mesh_refining_exponent;
	
	if ( mesh_index > _max_mesh_index )
		_max_mesh_index = mesh_index;
	
	
	if ( mesh_index < _min_mesh_index )
		_min_mesh_index = mesh_index;
}

/*-----------------------------------------------------------*/
/*              manually set the mesh index (static)         */
/*-----------------------------------------------------------*/
void NOMAD::Mesh::set_mesh_index ( int mesh_index )
{
  NOMAD::Mesh::_mesh_index = mesh_index;
  if ( mesh_index > _max_mesh_index )
    _max_mesh_index = mesh_index;
  if ( mesh_index < _min_mesh_index )
    _min_mesh_index = mesh_index;
}

/*-----------------------------------------------------------*/
/*                           display                         */
/*-----------------------------------------------------------*/
void NOMAD::Mesh::display ( const NOMAD::Display & out ) const
{
  out << "n                       : " << get_n()                   << std::endl
      << "mesh update basis       : " << _mesh_update_basis        << std::endl
      << "mesh coarsening exponent: " << _mesh_coarsening_exponent << std::endl
      << "mesh refining exponent  : " << _mesh_refining_exponent   << std::endl
      << "initial mesh index      : " << _initial_mesh_index       << std::endl
      << "initial mesh size       : "
      << "(" << _initial_mesh_size << " )" << std::endl;
  out << "minimal mesh size       : ";
  if ( _min_mesh_size.is_defined() )
    out << "(" << _min_mesh_size     << " )" << std::endl;
  else
    out << "none";
  out << std::endl
      << "minimal poll size       : ";
  if ( _min_poll_size.is_defined() )
    out << "(" << _min_poll_size     << " )" << std::endl;
  else
    out << "none";
  out << std::endl;
}

/*----------------------------------------------------------*/
/*  check the stopping conditions on the minimal poll size  */
/*  and on the minimal mesh size                            */
/*----------------------------------------------------------*/
void NOMAD::Mesh::check_min_mesh_sizes ( int                max_mesh_index ,
					 int                mesh_index     ,
					 bool             & stop           ,
					 NOMAD::stop_type & stop_reason      ) const
{
  if ( stop )
    return;

  // 1. mesh index tests:
  if ( abs ( mesh_index ) > NOMAD::L_LIMITS )
  {
    stop        = true;
    stop_reason = NOMAD::L_LIMITS_REACHED;
  }
  if ( max_mesh_index != NOMAD::UNDEFINED_L && mesh_index > max_mesh_index ) 
  {
    stop        = true;
    stop_reason = NOMAD::L_MAX_REACHED;
  }

  // 2. delta_k^p (poll size) tests:
  if ( check_min_poll_size_criterion ( mesh_index ) ) 
  {
    stop        = true;
    stop_reason = NOMAD::DELTA_P_MIN_REACHED;
  }

  // 3. delta_k^m (mesh size) tests:
  if ( check_min_mesh_size_criterion ( mesh_index ) ) 
  {
    stop        = true;
    stop_reason = NOMAD::DELTA_M_MIN_REACHED;
  }
}

/*-----------------------------------------------------------*/
/*              check the minimal poll size (private)        */
/*-----------------------------------------------------------*/
bool NOMAD::Mesh::check_min_poll_size_criterion ( int mesh_index ) const
{
  if ( !_min_poll_size.is_defined() )
    return false;
  NOMAD::Point delta_p;
  return get_delta_p ( delta_p , mesh_index );
}

/*-----------------------------------------------------------*/
/*              check the minimal mesh size (private)        */
/*-----------------------------------------------------------*/
bool NOMAD::Mesh::check_min_mesh_size_criterion ( int mesh_index ) const
{
  if ( !_min_mesh_size.is_defined() )
    return false;
  NOMAD::Point delta_m;
  return get_delta_m ( delta_m , mesh_index );
}

/*----------------------------------------------------------------*/
/*  get delta_m (mesh size parameter)                             */
/*       Delta^m_k = Delta^m_0 \tau^{ell_0^+ - ell_k^+}           */
/*----------------------------------------------------------------*/
bool NOMAD::Mesh::get_delta_m ( NOMAD::Point & delta_m , int mesh_index ) const
{
	int n = get_n();
	delta_m.reset ( n );
	
	// power_of_tau = tau^{ max{0,l0} - max{0,lk} }:
	NOMAD::Double power_of_tau
    = pow ( _mesh_update_basis ,
		   ( (_initial_mesh_index > 0) ? _initial_mesh_index : 0) -
		   ( (mesh_index          > 0) ? mesh_index          : 0)   );
	
	bool stop    = false;
	bool mms_def = _min_mesh_size.is_defined();
	
	// Delta^m_k = power_of_tau * Delta^m_0:
	for ( int i = 0 ; i < n ; ++i )
	{
		delta_m[i] = _initial_mesh_size[i] * power_of_tau;
		
		if ( mms_def && !stop && _min_mesh_size[i].is_defined() && delta_m[i] < _min_mesh_size[i] )
			stop = true;
		
	}
	
	return stop;
}

/*----------------------------------------------------------------*/
/*  get delta_m (mesh size parameter)                             */
/*       Delta^m_k = Delta^m_0 \tau^{ell_0^+ - ell_k^+}           */
/*----------------------------------------------------------------*/
/*  static version                                                */
/*----------------------------------------------------------------*/
void NOMAD::Mesh::get_delta_m ( NOMAD::Point       & delta_m            ,
				const NOMAD::Point & initial_mesh_size  ,
				double               mesh_update_basis  ,
				int                  initial_mesh_index ,
				int                  mesh_index           )
{
	
	int n = initial_mesh_size.size();
	delta_m.reset ( n );
	
	// power_of_tau = tau^{ max{0,l0} - max{0,lk} }:
	NOMAD::Double power_of_tau
    = pow ( mesh_update_basis ,
		   ( (initial_mesh_index > 0) ? initial_mesh_index : 0) -
		   ( (mesh_index         > 0) ? mesh_index         : 0)   );
	
	// Delta^m_k = power_of_tau * Delta^m_0:
	for ( int i = 0 ; i < n ; ++i )
		delta_m[i] = initial_mesh_size[i] * power_of_tau;
}

/*-------------------------------------------------------------------*/
/*  get delta_p (poll size parameter)                                */
/*       Delta^p_k = Delta^m_k \tau^{ |ell_k|/2 }                    */
/*                 = Delta^m_0 \tau^{ell_0^+ - ell_k^+ + |ell_k|/2}  */
/*-------------------------------------------------------------------*/
/*  the function also returns true if all values are < delta_p_min   */
/*-------------------------------------------------------------------*/
bool NOMAD::Mesh::get_delta_p ( NOMAD::Point & delta_p , int mesh_index ) const
{
	int n = get_n();
	delta_p.reset ( n );
	
	// power_of_tau = tau^{ max{0,l0} - max{0,lk} + |lk|/2}:
	NOMAD::Double power_of_tau
    = pow ( _mesh_update_basis , abs(mesh_index) / 2.0             +
		   ( (_initial_mesh_index > 0) ? _initial_mesh_index : 0) -
		   ( (mesh_index          > 0) ? mesh_index          : 0)   );
	
	bool stop    = true;
	
	bool mps_def = _min_poll_size.is_defined();
	
	// Delta^p_k = power_of_tau * Delta^m_0:
	for ( int i = 0 ; i < n ; ++i ) 
	{
		delta_p[i] = _initial_mesh_size[i] * power_of_tau;
		if ( !mps_def || !_min_poll_size[i].is_defined() || delta_p[i] >= _min_poll_size[i] )
			stop = false;
	}
	
	return stop;
}

/*-------------------------------------------------------------------*/
/*  get delta_p (poll size parameter)                                */
/*       Delta^p_k = Delta^m_k \tau^{ |ell_k|/2 }                    */
/*                 = Delta^m_0 \tau^{ell_0^+ - ell_k^+ + |ell_k|/2}  */
/*-------------------------------------------------------------------*/
/*  static version                                                   */
/*-------------------------------------------------------------------*/
void NOMAD::Mesh::get_delta_p ( NOMAD::Point       & delta_p            ,
				const NOMAD::Point & initial_mesh_size  ,
				double               mesh_update_basis  ,
				int                  initial_mesh_index ,
				int                  mesh_index           ) {

  int n = initial_mesh_size.size();
  delta_p.reset ( n );

  // power_of_tau = tau^{ max{0,l0} - max{0,lk} + |lk|/2}:
  NOMAD::Double power_of_tau
    = pow ( mesh_update_basis , abs(mesh_index) / 2.0             +
	    ( (initial_mesh_index > 0) ? initial_mesh_index : 0) -
	    ( (mesh_index         > 0) ? mesh_index         : 0)   );

  // Delta^p_k = power_of_tau * Delta^m_0:
  for ( int i = 0 ; i < n ; ++i )
    delta_p[i] = initial_mesh_size[i] * power_of_tau;
}
