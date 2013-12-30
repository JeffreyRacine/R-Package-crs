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
  \file   Mesh.hpp
  \brief  Class for the MADS mesh (headers)
  \author Sebastien Le Digabel
  \date   2010-04-06
  \see    Mesh.cpp
*/
#ifndef __MESH__
#define __MESH__

#include "Point.hpp"

namespace NOMAD {

  /// Class for the MADS mesh.
  /**
     - The mesh in NOMAD is defined with the basic directions and the
       mesh size parameter Delta^m_k.
     - The mesh size parameter is defined with the mesh index (the integer ell_k).
     - Use the static function \c NOMAD::Mesh::get_mesh_index() to access the value
       of the mesh index ell_k.
     - The poll size parameter Delta^p_k is not used to define the mesh but
       to define the poll trial points.
     - At each MADS iteration the mesh is updated with
       Delta^m_k+1 = tau^w+ Delta^m_k and w+ >= 0 (dominating iteration)
       or with
       Delta^m_k+1 = tau^w- Delta^m_k and w- < 0 (iteration failure).
       The mesh is not changed after improving iterations.
     - Mesh and poll size parameters are stored as NOMAD::Point objects
       (one value for each variable) in order to ensure automatic scaling.
     - See the MADS papers for more details on the mesh.
  */
  class Mesh {

  private:

    static double _mesh_update_basis;        ///< Mesh update basis       (tau)
    static int    _mesh_coarsening_exponent; ///< Mesh coarsening exponent (w+)
    static int    _mesh_refining_exponent;   ///< Mesh refining exponent   (w-)
    static int    _initial_mesh_index;       ///< Initial mesh index    (ell_0)
    static int    _mesh_index;               ///< Mesh index            (ell_k)
    
    static int    _min_mesh_index;           ///< Minimal value reached by the mesh index
    static int    _max_mesh_index;           ///< Maximal value reached by the mesh index

    static int    _max_halton_index;         ///< Maximal Halton index for Ortho-MADS
    
    NOMAD::Point  _initial_mesh_size;        ///< Delta^m_0 (abs values; one delta/var)
    
    /// Equivalent of \c _min_poll_size for the mesh size parameter Delta^m_k.
    NOMAD::Point  _min_mesh_size;

    /// Delta^p_min (abs. values, one delta/var.); Can be undefined or incomplete.
    NOMAD::Point  _min_poll_size;

    /*--------------------------------------------------------------*/

    /// Private affectation operator.
    /**
       \param m The right-hand side object -- \b IN.
    */
    const Mesh & operator = ( const Mesh & m );

    /// Check the minimal poll size criterion.
    bool check_min_poll_size_criterion ( int mesh_index = Mesh::_mesh_index ) const;

    /// Check the minimal mesh size criterion.
    bool check_min_mesh_size_criterion ( int mesh_index = Mesh::_mesh_index ) const;

    /*--------------------------------------------------------------*/

  public:
    
    /// Constructor.
    /**
       \param delta_m_0    Initial mesh size Delta^m_0                      -- \b IN.
       \param delta_p_min  Minimal poll size Delta^p_min (may be undefined) -- \b IN.
       \param delta_m_min  Minimal mesh size Delta^m_min (may be undefined) -- \b IN.
    */
    Mesh ( const NOMAD::Point & delta_m_0   ,
	   const NOMAD::Point & delta_p_min ,
	   const NOMAD::Point & delta_m_min   );
    
    /// Copy constructor.
    /**
       \param m The copied object -- \b IN.
    */
    Mesh ( const Mesh & m )
      : _initial_mesh_size ( m._initial_mesh_size ) ,
	_min_mesh_size     ( m._min_mesh_size     ) ,
	_min_poll_size     ( m._min_poll_size     )   {}
    
    /// Destructor.
    virtual ~Mesh ( void ) {}
    
    /// Static members initialization.
    /**
       \param mesh_update_basis         Mesh update basis       (tau) -- \b IN.
       \param mesh_coarsening_exponent  Mesh coarsening exponent (w+) -- \b IN.
       \param mesh_refining_exponent    Mesh refining exponent   (w-) -- \b IN.
       \param initial_mesh_index        Initial mesh index ell_0      -- \b IN.
    */
    static void init ( double mesh_update_basis        ,
		       int    mesh_coarsening_exponent ,
		       int    mesh_refining_exponent   ,
		       int    initial_mesh_index         );
    
    /// Access to the initial mesh size.
    /**
       \return A NOMAD::Point for the initial mesh size.
    */
    const NOMAD::Point & get_initial_mesh_size ( void ) const
    {
      return _initial_mesh_size;
    }

    /// Access to the minimal mesh size.
    /**
       \return A NOMAD::Point for the minimal mesh size.
    */
    const NOMAD::Point & get_min_mesh_size ( void ) const { return _min_mesh_size; }

    /// Access to the minimal poll size.
    /**
       \return A NOMAD::Point for the minimal poll size.
    */
    const NOMAD::Point & get_min_poll_size ( void ) const { return _min_poll_size; }
    
    /// Access to the mesh index.
    /**
       \return An integer with the mesh index.
    */
    static int get_mesh_index ( void  ) { return Mesh::_mesh_index; }

    /// Manually set the mesh index.
    /**
       \param mesh_index The mesh index ell_k -- \b IN.
    */
    static void set_mesh_index ( int mesh_index );
    
    /// Access to the maximal mesh index.
    /**
       This corresponds to the maximal value of ell_k reached so far.
       \return An integer with the maximal mesh index.
    */
    static int get_max_mesh_index ( void ) { return Mesh::_max_mesh_index; }


    /// Manually set the maximum mesh index.
    /**
       This corresponds to the maximum value of ell_k reached so far.
       \param h The maximum mesh index -- \b IN.
	   */
	  static void set_max_mesh_index ( int h ) { _max_mesh_index = h; }
	  
    /// Access to the minimal mesh index.
    /**
       This corresponds to the minimal value of ell_k reached so far.
       \return An integer with the minimal mesh index.
    */
    static int get_min_mesh_index ( void ) { return Mesh::_min_mesh_index; }
    
	  /// Manually set the minimum mesh index.
	  /**
       This corresponds to the minimum value of ell_k reached so far.
       \param h The minimum mesh index -- \b IN.
    */
    static void set_min_mesh_index ( int h ) { _min_mesh_index = h; }

    /// Test if mesh index<max_mesh_index so far.
    /**
       \return True if mesh index lower than maximal mesh index False otherwise.
    */	
    // Use 1- Apply dynamic reduction of poll dir only if true
    // use 2- Ortho n+1. Apply selection n directions out of 2n based on target dir only if true. Otherwise the selection alternates +/-.  
    static bool mesh_index_is_not_max(void){return Mesh::_mesh_index<Mesh::_max_mesh_index; /*true*/ }  
	      
    /// Access to the mesh update basis tau.
    /**
       \return A double with the mesh update basis tau.
    */
    static double get_mesh_update_basis ( void ) { return Mesh::_mesh_update_basis; }
    
    /// Access to the maximal Halton index.
    /**
       This corresponds to the maximal value of the Halton index reached so far.
       \return An integer with the maximal Halton index.
    */
    static int get_max_halton_index ( void  ) { return _max_halton_index; }

    /// Manually set the maximal Halton index.
    /**
       This corresponds to the maximal value of the Halton index reached so far.
       \param h The maximal Halton index -- \b IN.
    */
    static void set_max_halton_index ( int h ) { _max_halton_index = h; }
    
    /// Update the mesh (#1).
    /**
       \param success    Type of success of the iteration           -- \b IN.
       \param mesh_index The mesh index before and after the update -- \b IN/OUT.
    */
    static void update ( NOMAD::success_type success , int & mesh_index );

    /// Update the mesh (#2).
    /**
       \param success Type of success of the iteration -- \b IN.
    */
    static void update ( NOMAD::success_type success )
    {
      Mesh::update ( success , Mesh::_mesh_index );
    }

    /// Access to the number of variables.
    /**
       This number of variables corresponds to the one used for the construction
       of the current NOMAD::Mesh object.
       \return An integer with the number of variables.
    */
    int get_n ( void ) const { return _initial_mesh_size.size(); }

    /// Access to the mesh size parameter Delta^m_k.
    /**
       - The mesh size parameter is computed with
         Delta^m_k = Delta^m_0 tau^{ell_0^+ - ell_k^+}.
       - It is a NOMAD::Point of size \c nc the number of free variables.
       \param delta_m    The mesh size parameter Delta^m_k -- \b OUT.
       \param mesh_index A mesh index ell_k -- \b IN
       -- \b optional (default = the current mesh index ell_k.)
       \return A boolean equal to \c true if all values are
               strictly inferior than the associated minimal
	       mesh size Delta^m_min
	       (stopping criterion MIN_MESH_SIZE).
    */
    virtual bool get_delta_m ( NOMAD::Point & delta_m                        ,
		       int            mesh_index = Mesh::_mesh_index   ) const;

    /// Access to the mesh size parameter Delta^m_k.
    /**
       - The mesh size parameter is computed with
         Delta^m_k = Delta^m_0 tau^{ell_0^+ - ell_k^+}.
       - It is a NOMAD::Point of size \c nc the number of free variables.
       - This method is static and allows to compute mesh sizes
         without a NOMAD::Mesh object.
       \param delta_m            Mesh size parameter Delta^m_k (size \c nc) -- \b OUT.
       \param initial_mesh_size  Initial mesh size Delta^m_0 (size \c nc)   -- \b IN.
       \param mesh_update_basis  Mesh update basis tau                      -- \b IN.
       \param initial_mesh_index Initial mesh index ell_0                   -- \b IN.
       \param mesh_index         Mesh index ell_k                           -- \b IN.
    */
    static void get_delta_m ( NOMAD::Point       & delta_m            ,
			      const NOMAD::Point & initial_mesh_size  ,
			      double               mesh_update_basis  ,
			      int                  initial_mesh_index ,
			      int                  mesh_index           );

    /// Access to the poll size parameter Delta^p_k.
    /**
       - The poll size parameter is computed with
         Delta^p_k = Delta^m_k tau^{ |ell_k|/2 }
         = Delta^m_0 tau^{ell_0^+ - ell_k^+ + |ell_k|/2}
       - It is a NOMAD::Point of size \c nc the number of free variables.
       \param delta_p    The poll size parameter Delta^p_k -- \b OUT.
       \param mesh_index A mesh index -- \b IN
       -- \b optional (default = the current mesh index ell_k.)
       \return A boolean equal to \c true if all values are
               strictly inferior than the associated minimal
	       poll size Delta^p_min
	       (stopping criterion MIN_POLL_SIZE).
    */
    virtual bool get_delta_p ( NOMAD::Point & delta_p                        ,
		       int            mesh_index = Mesh::_mesh_index   ) const;
    
    /// Access to the poll size parameter Delta^p_k.
    /**
       - The poll size parameter is computed with
         Delta^p_k = Delta^m_k tau^{ |ell_k|/2 }
         = Delta^m_0 tau^{ell_0^+ - ell_k^+ + |ell_k|/2}
       - It is a NOMAD::Point of size \c nc the number of free variables.
       - This method is static and allows to compute poll sizes
         without a NOMAD::Mesh object.
       \param delta_p            Poll size parameter Delta^p_k (size \c nc) -- \b OUT.
       \param initial_mesh_size  Initial mesh size Delta^m_0 (size \c nc)   -- \b IN.
       \param mesh_update_basis  Mesh update basis tau                      -- \b IN.
       \param initial_mesh_index Initial mesh index ell_0                   -- \b IN.
       \param mesh_index         Mesh index ell_k                           -- \b IN.
    */
    static void get_delta_p ( NOMAD::Point       & delta_p            ,
			      const NOMAD::Point & initial_mesh_size  ,
			      double               mesh_update_basis  ,
			      int                  initial_mesh_index ,
			      int                  mesh_index           );

    /// Check the stopping conditions on the minimal poll and mesh sizes.
    /**
       \param max_mesh_index Maximal mesh index ell_max -- \b IN.
       \param mesh_index     Current mesh index ell_k   -- \b IN.
       \param stop           Stop flag                  -- \b IN/OUT.
       \param stop_reason    Stop reason                -- \b OUT.
    */
    void check_min_mesh_sizes ( int                max_mesh_index ,
				int                mesh_index     ,
				bool             & stop           ,
				NOMAD::stop_type & stop_reason      ) const;
    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;
  };

  /// Display a NOMAD::Mesh object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param m   The NOMAD::Mesh object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
					      const NOMAD::Mesh    & m     )
  {
    m.display ( out );
    return out;
  }
}

#endif
