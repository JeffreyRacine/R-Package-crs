/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonsmooth Optimization by Mesh Adaptive Direct search - version 3.5        */
/*                                                                                     */
/*  Copyright (C) 2001-2010  Mark Abramson        - the Boeing Company, Seattle        */
/*                           Charles Audet        - Ecole Polytechnique, Montreal      */
/*                           Gilles Couture       - Ecole Polytechnique, Montreal      */
/*                           John Dennis          - Rice University, Houston           */
/*                           Sebastien Le Digabel - Ecole Polytechnique, Montreal      */
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
  \file   Phase_One_Evaluator.hpp
  \brief  NOMAD::Evaluator subclass for the phase one (headers)
  \author Sebastien Le Digabel
  \date   2010-04-22
  \see    Phase_One_Evaluator.cpp
*/
#ifndef __PHASE_ONE_EVALUATOR__
#define __PHASE_ONE_EVALUATOR__

#include "Evaluator.hpp"

namespace NOMAD {

  /// NOMAD::Evaluator subclass for the phase one.
  class Phase_One_Evaluator : public NOMAD::Evaluator {

  private:

    NOMAD::Evaluator & _basic_ev; ///< The original evaluator.

  public:

    /// Constructor.
    /**
       \param p  Parameters         -- \b IN.
       \param ev Original evaluator -- \b IN.
    */
    Phase_One_Evaluator ( const NOMAD::Parameters & p , NOMAD::Evaluator & ev )
      : NOMAD::Evaluator ( p  ) ,
	_basic_ev        ( ev )   {}

    /// Destructor.
    virtual ~Phase_One_Evaluator ( void ) {}

    /// User updates after a success.
    /**
       This virtual method is called every time a new (full) success is made.
       \param s Stats -- \b IN.
       \param x Successful point -- \b IN.
    */
    virtual void update_success ( const NOMAD::Stats & s , const NOMAD::Eval_Point & x )
    {
      _basic_ev.update_success ( s , x );
    }

    /// Evaluate the blackboxes at a given trial point.
    /**
       \param x The trial point -- \b IN/OUT.
       \param h_max      Maximal feasibility value \c h_max -- \b IN.
       \param count_eval Flag indicating if the evaluation has to be counted
                         or not -- \b OUT.
       \return A boolean equal to \c false if the evaluation failed.
     */
    virtual bool eval_x ( NOMAD::Eval_Point   & x          ,
			  const NOMAD::Double & h_max      ,
			  bool                & count_eval   ) const
    {
      return _basic_ev.eval_x ( x , h_max , count_eval );
    }

    /// User preprocessing of points before evaluations.
    /**
       This method is called before the evaluation of a list of points.
       \param pts List of points to preprocess -- \b IN/OUT.
    */
    virtual void list_of_points_preprocessing
    ( std::set<NOMAD::Priority_Eval_Point> & pts ) const
    {
      _basic_ev.list_of_points_preprocessing ( pts );
    }
    
    /// Objective computation.
    /**
       - Compute \c f(x) from the blackbox outputs of a point.
       - Special objective for MADS phase one.
       \param x The trial point -- \b IN/OUT.
    */
    virtual void compute_f ( NOMAD::Eval_Point & x ) const;
  };
}

#endif
