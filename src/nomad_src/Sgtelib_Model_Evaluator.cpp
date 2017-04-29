/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct search - version 3.8.0      */
/*                                                                                     */
/*                                                                                     */
/*  NOMAD - version 3.8.0 has been created by                                          */
/*                 Charles Audet        - Ecole Polytechnique de Montreal              */
/*                 Sebastien Le Digabel - Ecole Polytechnique de Montreal              */
/*                 Christophe Tribes    - Ecole Polytechnique de Montreal              */
/*                                                                                     */
/*  The copyright of NOMAD - version 3.8.0 is owned by                                 */
/*                 Sebastien Le Digabel - Ecole Polytechnique de Montreal              */
/*                 Christophe Tribes    - Ecole Polytechnique de Montreal              */
/*                                                                                     */
/*  NOMAD v3 has been funded by AFOSR and Exxon Mobil.                                 */
/*                                                                                     */
/*  NOMAD v3 is a new version of NOMAD v1 and v2. NOMAD v1 and v2 were created and     */
/*  developed by Mark Abramson, Charles Audet, Gilles Couture and John E. Dennis Jr.,  */
/*  and were funded by AFOSR and Exxon Mobil.                                          */
/*                                                                                     */
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
 \file   Sgtelib_Model_Evaluator.cpp
 \brief  Interface between nomad evaluator and Sgtelib_Model_Manager.
 \author Bastien Talgorn
 \date   2013-04-25
 \see    Sgtelib_Model_Manager.cpp
 */

#include "Sgtelib_Model_Evaluator.hpp"


/*------------------------------------------------------------------------*/
/*                evaluate the sgtelib_model model at a given trial point           */
/*------------------------------------------------------------------------*/
bool NOMAD::Sgtelib_Model_Evaluator::eval_x ( NOMAD::Eval_Point   & x          ,
                                             const NOMAD::Double & h_max      ,
                                             bool                & count_eval   ) const
{
    return _sgtelib_model_manager->eval_x( &x         ,
                                          h_max      ,
                                          count_eval );
}



// #endif

