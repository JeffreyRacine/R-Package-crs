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
 \file   L_Curve.hpp
 \brief  L_CURVE_TARGET stopping criterion (headers)
 \author Sebastien Le Digabel
 \date   2010-04-09
 \see    L_Curve.cpp
 */
#ifndef __L_CURVE__
#define __L_CURVE__

#include "Double.hpp"
#include "Uncopyable.hpp"

namespace NOMAD {
    
    /// Class implementing the L_CURVE_TARGET stopping criterion.
    class L_Curve : private NOMAD::Uncopyable {
        
    private:
        
        NOMAD::Double              _target;  ///< L_CURVE_TARGET parameter value.
        std::vector<NOMAD::Double> _f;       ///< List of objective values.
        std::vector<int          > _bbe;     ///< List of numbers of evaluations.
        
    public:
        
        /// Constructor.
        /**
         \param target L_CURVE_TARGET parameter value -- \b IN.
         */
        L_Curve ( const NOMAD::Double & target ) : _target ( target ) {}
        
        /// Destructor.
        virtual ~L_Curve ( void ) {}
        
        /// Insertion of a pair \c bbe/f in the lists \c _f and \c _bbe.
        /**
         \param bbe A new number of evaluations -- \b IN.
         \param f   A new objective value       -- \b IN.
         */
        void insert ( int bbe , const NOMAD::Double & f );
        
        /// Check the L_CURVE_TARGET stopping criterion.
        /**
         \param bbe An integer indicating a number of blackbox evaluations
         -- \b IN.
         \return A boolean equal to \c true if the method detects that
         the target will not be reached after bbe evaluations.
         */
        bool check_stop ( int bbe ) const;
    };
}
#endif
