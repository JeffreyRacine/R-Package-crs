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
  \file   Quad_Model_Sorted_Point.hpp
  \brief  Interpolation point with distance to model center (headers)
  \author Sebastien Le Digabel
  \date   2010-11-15
  \see    Quad_Model_Sorted_Point.cpp
*/
#ifndef __QUAD_MODEL_SORTED_POINT__
#define __QUAD_MODEL_SORTED_POINT__

#include "Eval_Point.hpp"

namespace NOMAD {

  /// Class used to order interpolation points.
  class Quad_Model_Sorted_Point {

  private:

    NOMAD::Eval_Point * _x;    ///< The interpolation point.
    NOMAD::Double       _dist; ///< Distance to model center.

  public:

    /// Constructor.
    /**
       \param x      Interpolaton point -- \b IN.
       \param center Model center       -- \b IN.
    */
    Quad_Model_Sorted_Point ( NOMAD::Eval_Point * x , const NOMAD::Point & center );

    /// Copy constructor.
    /**
       \param x The copied object -- \b IN.
    */
    Quad_Model_Sorted_Point ( const Quad_Model_Sorted_Point & x )
      : _x ( x._x ) , _dist ( x._dist ) {}
  
    /// Affectation operator.
    /**
       \param x The right-hand side object -- \b IN.
       \return \c *this as the result of the affectation.
    */
    Quad_Model_Sorted_Point & operator = ( const Quad_Model_Sorted_Point & x );

    /// Comparison operator.
    /**
       \param x The right-hand side object -- \b IN.
       \return \c true if the current interpolation point is closer to the center.
    */
    bool operator < ( const Quad_Model_Sorted_Point & x ) const;
  
    /// Access to the interpolation point.
    /**
       \return The interpolation point.
    */
    NOMAD::Eval_Point * get_point ( void ) const { return _x; }
  };
}

#endif
