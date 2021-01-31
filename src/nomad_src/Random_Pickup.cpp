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
 \file   Random_Pickup.cpp
 \brief  Class for randomly pick up integers (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-07
 \see    Random_Pickup.hpp
 */
#include "Random_Pickup.hpp"

/*---------------------------------------------------------*/
/*                         constructor                     */
/*---------------------------------------------------------*/
NOMAD::Random_Pickup::Random_Pickup ( int n )
: _n0   ( n          ) ,
_n    ( n          ) ,
_elts ( new int[n] )
{
    for ( int i = 0 ; i < n ; ++i )
        _elts[i] = i;
}

/*---------------------------------------------------------*/
/*                           reset                         */
/*---------------------------------------------------------*/
void NOMAD::Random_Pickup::reset ( void )
{
    _n = _n0;
    for ( int i = 0 ; i < _n ; ++i )
        _elts[i] = i;
}

/*---------------------------------------------------------*/
/*                randomly pick up an element              */
/*---------------------------------------------------------*/
int NOMAD::Random_Pickup::pickup ( void )
{
    if ( _n == 0 )
        return 0;
    int ind = NOMAD::RNG::rand()%_n;
    int tmp = _elts[ind];
    if ( ind < _n - 1 )
    {
        _elts[ind ] = _elts[_n-1];
        _elts[_n-1] = tmp;
    }
    --_n;
    
    return tmp;
}

/*---------------------------------------------------------*/
/*                   cancel the last pick up               */
/*---------------------------------------------------------*/
void NOMAD::Random_Pickup::cancel_last_pickup ( void )
{
    if ( _n < _n0 )
        ++_n;
}
