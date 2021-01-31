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
 \file   NelderMead_Simplex_Eval_Point.cpp
 \brief  Evaluation point with a NelderMead order priority (implementation)
 \author Christophe Tribes
 \date   2017-04-18
 \see    NelderMead_Simplex_Eval_Point.hpp
 */
#include "NelderMead_Simplex_Eval_Point.hpp"


NOMAD::Double NOMAD::NelderMead_Simplex_Eval_Point::_h_min = 0;

/*------------------------------------------------*/
/*                comparison operator             */
/*------------------------------------------------*/
/*  . x1.is_better_than(x2) returns true if       */
/*    x1=Best(x1,x2)                              */
/*------------------------------------------------*/
bool NOMAD::NelderMead_Simplex_Eval_Point::is_better_than ( const NOMAD::Set_Element<NOMAD::Eval_Point> & x ) const
{
    
    if ( this == &x )
        return false;
    
    const NOMAD::Eval_Point * x1 = get_element();
    const NOMAD::Eval_Point * x2 = x.get_element();
    
    // criterion 1: true f and h values:
    // -----------
    if ( dominates(*x1,*x2) )
        return true;
    
    if ( dominates(*x2,*x1) )
        return false;
    
    NOMAD::Double hx1 = x1->get_h();
    NOMAD::Double hx2 = x2->get_h();
    
    // criterion 2: take the point with the best h value (case no infeasible EB constraints):
    // ------------
    if ( hx1.is_defined() && hx2.is_defined() )
    {
        if ( hx1.value() < hx2.value() )
            return true;
        if ( hx1.value() > hx2.value() )
            return false;
    }
    
    // criterion 3a and 3b:
    //       - case x1 feasible or infeasible and x2 has infeasible EB constraints
    //       - case x2 feasible of infeasible and x1 has infeasible EB constraints
    // -------------------
    if ( hx1.is_defined() && ! hx2.is_defined() )
        return true;
    if ( hx2.is_defined() && ! hx1.is_defined() )
        return false;
    
    
    
    // criterion 4: compare the tags:
    // -------------
    return x1->get_tag() < x2->get_tag();
    
}


/*-----------------------------------------------------------------------*/
/*                comparison operator                                    */
/*-----------------------------------------------------------------------*/
/*  returns true if *this dominates < the right-hand side object x       */
/*-----------------------------------------------------------------------*/
bool NOMAD::NelderMead_Simplex_Eval_Point::dominates ( const NOMAD::Eval_Point & x1, const NOMAD::Eval_Point & x2 )
{
    
    NOMAD::Double fx2 = x2.get_f();
    NOMAD::Double fx1 = x1.get_f();
    
    NOMAD::Double hx2 = x2.get_h();
    NOMAD::Double hx1 = x1.get_h();
    
    if ( fx1.is_defined() && fx2.is_defined() )
    {
        
        // Compare when EB constraints are present ( case x1 PB non feasible and x2 EB non feasible: 0 < h(x1) < +Inf and h(x2) = +Inf => hx1 < hx2)
        if ( hx1.is_defined() && ! hx2.is_defined() &&  hx1.value() > 0 && fx1.value()  <= fx2.value() )
            return true;
        
        // Compare when EB constraints are present ( case x1 EB non feasible and x2 PB non feasible: 0 < h(x2) < +Inf and h(x1) = +Inf => hx2 < hx1 )
        if ( hx2.is_defined() && ! hx1.is_defined() )
            return false;
        
        if ( hx1.is_defined() && hx2.is_defined() )
        {
            
            // x1 is feasible:
            if ( hx1.value() <= _h_min.value() )
            {
                // both points are feasible:
                if ( hx2.value() <= _h_min.value()  )
                {
                    if ( fx1.value() < fx2.value() )
                        return true;
                    if ( fx2.value() < fx1.value() )
                        return false;
                }
                
                // x1 feasible and x2 infeasible:
                else
                    return false;
            }
            
            // x1 is infeasible:
            else
            {
                // x2 is feasible:
                if ( hx2.value() <= _h_min.value()  )
                    return false;
                
                // both points are infeasible:
                if ( ( hx1.value() < hx2.value() && fx1.value()  < fx2.value() ) ||
                    ( hx1.value() == hx2.value() && fx1.value()  < fx2.value() ) ||
                    ( hx1.value()  < hx2.value() && fx1.value() == fx2.value() )    )
                    return true;
                
                return false;
                
            }
        }
        
        // we only have f values:
        else
        {
            if ( fx1.value() < fx2.value() )
                return true;
            if ( fx2.value() < fx1.value() )
                return false;
        }
    }
    else
        throw NOMAD::Exception ( "NelderMead_Simplex_Eval_Point.cpp" , __LINE__ ,
                            "NelderMead_Simplex_Eval_Point::dominates(): could not compare points. Objective function value not defined." );
    return false;
}