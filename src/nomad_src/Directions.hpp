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
  \file   Directions.hpp
  \brief  Set of polling directions (headers)
  \author Sebastien Le Digabel
  \date   2010-04-14
  \see    Directions.cpp
*/
#ifndef __DIRECTIONS__
#define __DIRECTIONS__

#include "Random_Pickup.hpp"
#include "Direction.hpp"
#include "Mesh.hpp"

namespace NOMAD {

  /// Set of polling directions.
  class Directions {

    /*--------------------------------------------------------------*/
  private:

    int _nc;  ///< Number of non-fixed variables.

    /// Types of the poll directions.
    std::set<NOMAD::direction_type> _direction_types;

    /// Types of the secondary poll directions.
    std::set<NOMAD::direction_type> _sec_poll_dir_types;

    bool _is_binary;      ///< Flag equal to \c true if all variables are binary.
    bool _is_categorical; ///< Flag equal to \c true if all variables are categorical.

    bool _is_orthomads;  //< Flag equal to \c true if Ortho-MADS directions are used.

    // LT-MADS attributes:
    NOMAD::Direction * _bl [1+2*NOMAD::L_LIMITS];  ///< Directions b(ell) (LT-MADS).
    int             _hat_i [1+2*NOMAD::L_LIMITS];  ///< b(ell) indexes (LT-MADS).
    bool   _lt_initialized;  ///< Flag equal to \c true if LT-MADS has been initialized.

    // Ortho-MADS members:
    int      * _primes;           ///< \c nc first prime numbers.
    int        _halton_seed;      ///< Halton seed \c t_0.
    int        _ortho_np1_cnt1;   ///< First counter for Ortho-MADS n+1 directions.
    int        _ortho_np1_cnt2;   ///< Second counter for Ortho-MADS n+1 directions.
    static int _max_halton_seed;  ///< Max Halton seed for all NOMAD::Directions objects.
    
    const NOMAD::Display & _out;  ///< Display.

    /*--------------------------------------------------------------*/

    /// Affectation operator.
    /**
       \param d The right-hand side object -- \b IN.
    */
    Directions & operator = ( const Directions & d );

    /// Ortho-MADS initializations.
    /**
       \param halton_seed Halton seed -- \b IN.
    */
    void ortho_mads_init ( int halton_seed );

    /// LT-MADS initializations.
    void lt_mads_init ( void );

    /// Compute the squared norm of \c normalized(2u_t-e) for Ortho-MADS.
    /**
       \param x     \b IN.
       \param norm  \b IN.
       \param b     \b IN.
       \param new_b \b OUT.
       \return The squared norm.
    */
    NOMAD::Double eval_ortho_norm ( const NOMAD::Double & x      ,
				    const NOMAD::Double & norm   ,
				    const NOMAD::Point  & b      ,
				    NOMAD::Point        & new_b    ) const;

    /// Householder transformation.
    /**
       - Generate \c _nc directions from the Halton direction.
       - Computes also \c H[i+nc] \c = \c -H[i] (completion to 2n directions).
       \param halton_dir     The Halton direction                -- \b IN.
       \param complete_to_2n Complete or not to \c 2n directions -- \b IN.
       \param H              The \c _nc directions               -- \b OUT.
    */
    void householder ( const NOMAD::Direction & halton_dir     ,
		       bool                     complete_to_2n ,
		       NOMAD::Direction      ** H                ) const;

    /// Get the expression of an integer \c t in inverse base \c p.
    /**
       \param t \b IN.
       \param p \b IN.
       \return The expression of \c t in inverse base \c p.
    */
    static NOMAD::Double get_phi ( int t , int p );

    /// Compute the Halton direction \c q(t,ell).
    /**
       \param halton_index Halton index \c t          -- \b IN.
       \param mesh_index   Mesh index \c ell          -- \b IN.
       \param alpha_t_l    Ortho-MADS \c alpha(t,ell) -- \b OUT.
       \param halton_dir   Halton direction           -- \b OUT.
       \return A boolean equal to \c true if the computation went well.
    */
    bool compute_halton_dir ( int                halton_index ,
			      int                mesh_index   ,
			      NOMAD::Double    & alpha_t_l    ,
			      NOMAD::Direction & halton_dir     ) const;

    /// Compute the target direction for Ortho-MADS n+1.
    /**
       \param feas_success_dir   Feasible success direction             -- \b IN.
       \param infeas_success_dir Infeasible success direction           -- \b IN.
       \param first_success      First success of the run (can be NULL) -- \b IN.
       \param poll_center        Poll center                            -- \b IN.
       \param target_dir         The target direction                   -- \b OUT.
    */
    void compute_ortho_target_dir ( const NOMAD::Direction & feas_success_dir   ,
				    const NOMAD::Direction & infeas_success_dir ,
				    const NOMAD::Point     * first_success      ,
				    const NOMAD::Point     & poll_center        ,
				    NOMAD::Point           & target_dir           );

    /// Access to the LT-MADS \c b(ell) direction.
    /**
       \param mesh_index Mesh index \c ell  -- \b IN.
       \param dtype Direction type          -- \b IN.
       \param hat_i LT-MADS \c hat{i} index -- \b IN/OUT.
       \return The LT-MADS \c b(ell) direction.
    */
    const NOMAD::Direction * get_bl ( int                     mesh_index ,
				      NOMAD::direction_type   dtype      ,
				      int                   & hat_i        );
    
    /// Create a new LT-MADS direction.
    /**
       If \c hat_i \c == \c -1, a new \c b(ell) direction is created
       and \c hat_i is set.
       \param mesh_index Mesh index \c ell           -- \b IN.
       \param dtype      Direction type              -- \b IN.
       \param diag_i     Diagonal index              -- \b IN.
       \param hat_i      LT-MADS \c hat{i} index     -- \b IN/OUT.
       \param dir        LT-MADS \c b(ell) direction -- \b OUT.
    */
    void create_lt_direction ( int                     mesh_index ,
			       NOMAD::direction_type   dtype      ,
			       int                     diag_i     ,
			       int                   & hat_i      ,
			       NOMAD::Direction     *& dir          );

    /// Permute the coordinates of a direction.
    /**
       \param dir                The direction      -- \b IN/OUT.
       \param permutation_vector Permutation vector -- \b IN.
    */
    void permute_coords ( NOMAD::Direction & dir                ,
			  const int        * permutation_vector   ) const;
    
    /// Compute binary directions.
    /**
       Only if all groups of variables are binary.
       \param dirs Set of directions -- \b OUT.
    */
    void compute_binary_directions ( std::list<NOMAD::Direction> & dirs ) const;
    
    /*--------------------------------------------------------------*/

  public:

    /// Constructor.
    /**
       The Halton seed will be automatically computed later
       if NOMAD::Parameters::_halton_seed==-1.
       \param nc                 Number of non-fixed variables          -- \b IN.
       \param direction_types    Types of the poll directions           -- \b IN.
       \param sec_poll_dir_types Types of the secondary poll directions -- \b IN.
       \param halton_seed        Halton seed \c t_0                     -- \b IN.
       \param out                The display                            -- \b IN.
    */
    Directions ( int                                     nc                 ,
		 const std::set<NOMAD::direction_type> & direction_types    ,
		 const std::set<NOMAD::direction_type> & sec_poll_dir_types ,
		 int                                     halton_seed        ,
		 const NOMAD::Display                  & out                  );

    /// Copy constructor.
    /**
       \param d The copied object -- \b IN.
    */
    Directions ( const Directions & d )
      : _nc                 ( d._nc                 ) ,
	_direction_types    ( d._direction_types    ) ,
	_sec_poll_dir_types ( d._sec_poll_dir_types ) ,
	_is_binary          ( d._is_binary          ) ,
	_is_categorical     ( d._is_categorical     ) ,
	_is_orthomads       ( d._is_orthomads       ) ,
	_lt_initialized     ( false                 ) ,
	_primes             ( NULL                  ) ,
	_halton_seed        ( d._halton_seed        ) ,
	_ortho_np1_cnt1     ( d._ortho_np1_cnt1     ) ,
	_ortho_np1_cnt2     ( d._ortho_np1_cnt2     ) ,
	_out                ( d._out                )   {}
    
    /// Destructor.
    virtual ~Directions ( void );
    
    /// Compute the directions for a given mesh.
    /**
       \param dirs               Set of directions                                          -- \b OUT.
       \param poll               Type of poll (primary or secondary)                        -- \b IN.
       \param poll_center        Poll center                                                -- \b IN.
       \param first_success      First success of the run (for Ortho-MADS n+1; can be NULL) -- \b IN.    
       \param mesh_index         Mesh index \c ell                                          -- \b IN.
       \param halton_index       Halton index \c t (set to \c -1 for a default value)       -- \b IN.
       \param feas_success_dir   Feasible success direction                                 -- \b IN.
       \param infeas_success_dir Infeasible success direction                               -- \b IN.
    */
    void compute ( std::list<NOMAD::Direction> & dirs               ,
		   NOMAD::poll_type              poll               ,
		   const NOMAD::Point          & poll_center        ,
		   const NOMAD::Point          * first_success      ,
		   int                           mesh_index         ,
		   int                           halton_index       ,
		   const NOMAD::Direction      & feas_success_dir   ,
		   const NOMAD::Direction      & infeas_success_dir   );
    
    /// Compute just one direction for a given mesh (used by VNS search).
    /**
       \param dir          One direction                       -- \b OUT.
       \param mesh_index   Mesh index \c ell                   -- \b IN.
       \param halton_index Halton index \c t (must be \c >0)   -- \b IN.
       \return A boolean equal to \c true if the computation went well.
    */
    bool compute_one_direction ( NOMAD::Direction & dir          ,
				 int                mesh_index   ,
				 int                halton_index   );
    
    /// Method for computing a Halton seed.
    /**
       The Halton seed for \c n variables is computed as the \c n th prime number.
       \param n Number of variables -- \b IN.
       \return The Halton seed \c t_0.
    */
    static int compute_halton_seed ( int n );

    /// Access to the max value of the Halton index \c t.
    /**
       \return The max value of the Halton index \c t.
    */
    static int get_max_halton_seed ( void ) { return Directions::_max_halton_seed; }

    /// Access to the Halton seed \c t_0.
    /**
       \return The Halton seed \c t_0.
    */
    int get_halton_seed ( void ) const { return _halton_seed; }

    /// Check if Ortho-MADS directions are used.
    /**
       \return A boolean equal to \c true if Ortho-MADS directions are used.
    */
    bool is_orthomads ( void ) const { return _is_orthomads; }

    /// Check if all variables are categorical.
    /**
       \return A boolean equal to \c true if all variables are categorical.
    */
    bool is_categorical ( void ) const { return _is_categorical; }

    /// Access to the poll direction types.
    /**
       \return Poll direction types.
    */
    const std::set<NOMAD::direction_type> & get_direction_types ( void ) const
    {
      return _direction_types;
    } 

    /// Access to the secondary poll direction types.
    /**
       \return Secondary poll direction types.
    */
    const std::set<NOMAD::direction_type> & get_sec_poll_dir_types ( void ) const
    {
      return _sec_poll_dir_types;
    }
    
    /// Reset directions for binary variables.
    void set_binary ( void );
    
    /// Reset directions for categorical variables.
    void set_categorical ( void );
    
    /// Comparison operator.
    /**
       \param  d The right-hand side object -- \b IN.
       \return A boolean equal to \c true if \c *this \c < \c d .
    */
    bool operator < ( const Directions & d ) const;

    /// Access to the NOMAD::Display member \c _out.
    /**
       \return The NOMAD::Display member \c _out.
    */
    const NOMAD::Display & out ( void ) const { return _out; }

    /// Display.
    /**
       \param out The NOMAD::Display object -- \b IN.
    */
    void display ( const NOMAD::Display & out ) const;

    /// Display.
    /**
       Uses the \c this->_out member as NOMAD::Display object.
    */
    void display ( void ) const { display ( _out ); }
  };
  
  /// Display a NOMAD::Directions object.
  /**
     \param out The NOMAD::Display object -- \b IN.
     \param d   The NOMAD::Directions object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
  */
  inline const NOMAD::Display & operator << ( const NOMAD::Display    & out ,
					      const NOMAD::Directions & d     )
  {
    d.display ( out );
    return out;
  }
}

#endif
