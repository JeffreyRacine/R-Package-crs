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
  \file   Speculative_Search.hpp
  \brief  Speculative search (headers)
  \author Sebastien Le Digabel
  \date   2010-04-12
  \see    Speculative.cpp
*/
#ifndef __SPECULATIVE_SEARCH__
#define __SPECULATIVE_SEARCH__

#include "Search.hpp"
#include "Mads.hpp"

namespace NOMAD {

  /// Speculative search.
  /**
     The speculative search consists in looking further away along
     the successful direction after an improvement.
  */
  class Speculative_Search : public NOMAD::Search , private NOMAD::Uncopyable {

  public:

    /// Constructor.
    /**
       \param p Parameters -- \b IN.
    */
    Speculative_Search ( NOMAD::Parameters & p )
      : NOMAD::Search ( p , NOMAD::SPEC_SEARCH ) {}
    
    /// Destructor.
    virtual ~Speculative_Search ( void ) {}

    /// The speculative search.
    /**
       - x_k = x_{k-1} + Delta^m_{k-1} d \n
         s_k = x_{k-1} +        Delta^m_{k-1} d  if ell_{k-1} > 0 \n
                    or    tau * Delta^m_{k-1} d  otherwise.
       - The directions that we use already contain Delta^m:
         direction = Delta^m_{k-1} d
       - The following test
         \code
          if ( new_feas_inc || new_infeas_inc )
         \endcode
         is equal to \c true and has already been made in \c Mads.cpp.

      \param mads           NOMAD::Mads object invoking this search -- \b IN/OUT.
      \param nb_search_pts  Number of generated search points       -- \b OUT.
      \param stop           Stop flag                               -- \b IN/OUT.
      \param stop_reason    Stop reason                             -- \b OUT.
      \param success        Type of success                         -- \b OUT.
      \param count_search   Count or not the search                 -- \b OUT.
      \param new_feas_inc   New feasible incumbent                  -- \b IN/OUT.
      \param new_infeas_inc New infeasible incumbent                -- \b IN/OUT.
    */
    virtual void search ( NOMAD::Mads              & mads           ,
                          int                      & nb_search_pts  ,
                          bool                     & stop           ,
                          NOMAD::stop_type         & stop_reason    ,
                          NOMAD::success_type      & success        ,
                          bool                     & count_search   ,
                          const NOMAD::Eval_Point *& new_feas_inc   ,
                          const NOMAD::Eval_Point *& new_infeas_inc   );
  };
}

#endif
