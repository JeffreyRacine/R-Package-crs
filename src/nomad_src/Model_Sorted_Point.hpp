/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search -                */
/*                                                                                 */
/*  NOMAD - version 3.9.1 has been created by                                      */
/*                 Charles Audet               - Ecole Polytechnique de Montreal   */
/*                 Sebastien Le Digabel        - Ecole Polytechnique de Montreal   */
/*                 Viviane Rochon Montplaisir - Ecole Polytechnique de Montreal   */
/*                 Christophe Tribes           - Ecole Polytechnique de Montreal   */
/*                                                                                 */
/*  The copyright of NOMAD - version 3.9.1 is owned by                             */
/*                 Sebastien Le Digabel        - Ecole Polytechnique de Montreal   */
/*                 Viviane Rochon Montplaisir - Ecole Polytechnique de Montreal   */
/*                 Christophe Tribes           - Ecole Polytechnique de Montreal   */
/*                                                                                 */
/*  NOMAD v3 has been funded by AFOSR and Exxon Mobil.                             */
/*                                                                                 */
/*  NOMAD v3 is a new version of NOMAD v1 and v2. NOMAD v1 and v2 were created     */
/*  and developed by Mark Abramson, Charles Audet, Gilles Couture, and John E.     */
/*  Dennis Jr., and were funded by AFOSR and Exxon Mobil.                          */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Ecole Polytechnique de Montreal - GERAD                                      */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/

/**
 \file   Model_Sorted_Point.hpp
 \brief  Interpolation point with distance to model center (headers)
 \author Sebastien Le Digabel
 \date   2010-11-15
 \see    Model_Sorted_Point.cpp
 */
#ifndef __MODEL_SORTED_POINT__
#define __MODEL_SORTED_POINT__

#include "Point.hpp"

namespace NOMAD {
    
    /// Class used to order interpolation points.
    class Model_Sorted_Point {
        
    private:
        
        NOMAD::Point  * _x;    ///< The point.
        NOMAD::Double   _dist; ///< Distance to center.
        
    public:
        
        /// Constructor 1/2.
        /**
         \param x      Interpolaton point -- \b IN.
         \param center Model center       -- \b IN.
         */
        Model_Sorted_Point ( NOMAD::Point * x , const NOMAD::Point & center );
        
        /// Constructor 2/2.
        /**
         \param x      Interpolaton point          -- \b IN.
         \param dist   Custom distance with center -- \b IN.
         */
        Model_Sorted_Point ( NOMAD::Point * x , const NOMAD::Double & dist )
        : _x ( x ) , _dist ( dist ) {}
        
        /// Copy constructor.
        /**
         \param x The copied object -- \b IN.
         */
        Model_Sorted_Point ( const Model_Sorted_Point & x )
        : _x ( x._x ) , _dist ( x._dist ) {}
        
        /// Affectation operator.
        /**
         \param x The right-hand side object -- \b IN.
         \return \c *this as the result of the affectation.
         */
        Model_Sorted_Point & operator = ( const Model_Sorted_Point & x );
        
        /// Destructor.
        virtual ~Model_Sorted_Point ( void ) {}
        
        /// Comparison operator.
        /**
         \param x The right-hand side object -- \b IN.
         \return \c true if the current interpolation point is closer to the center.
         */
        bool operator < ( const Model_Sorted_Point & x ) const;
        
        /// Access to the interpolation point.
        /**
         \return The interpolation point.
         */
        NOMAD::Point * get_point ( void ) const { return _x; }
        
        /// Access to the distance.
        /**
         \return The distance
         */
        const NOMAD::Double & get_dist ( void ) const { return _dist; }
    };
}

#endif
