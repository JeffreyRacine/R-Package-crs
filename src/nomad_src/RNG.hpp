/*-------------------------------------------------------------------------------------*/
/*  NOMAD - Nonsmooth Optimization by Mesh Adaptive Direct search - version 3.6.2        */
/*                                                                                     */
/*  Copyright (C) 2001-2010  Mark Abramson        - the Boeing Company, Seattle        */
/*                           Charles Audet        - Ecole Polytechnique, Montreal      */
/*                           Gilles Couture       - Ecole Polytechnique, Montreal      */
/*                           John Dennis          - Rice University, Houston           */
/*                           Sebastien Le Digabel - Ecole Polytechnique, Montreal      */
/*                           Christophe Tribes    - Ecole Polytechnique, Montreal      */
/*                                                                                     */
/*  funded in part by AFOSR and Exxon Mobil                                            */
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
 \file   RNG.hpp
 \brief  Custom class for random number generator
 \author Christophe Tribes and Sebastien Le Digabel 
 \date   2011-09-28
 \see    RNG.cpp
 */

#ifndef __RNG__
#define __RNG__

#include "defines.hpp"
#include "Exception.hpp"

using namespace std;

namespace NOMAD {

	
	/// Class for random number generator 
	/**
		This class is used to set a seed for the random number generator and
	    get a random integer or a random double between two values.
	 */
	class RNG {
		
	public:

		/// Set seed
		/*
		 /param s The seed -- \b IN.
		 /return A boolean if the seed is acceptable, that is in [0,UINT32_MAX].
		 */
		static bool set_seed(int s);
		
		/// Get a random integer as uint32
		/** This function serves to obtain a random number \c 
		 /return An integer in the interval [0,UINT32_MAX].
		 */
		static uint32_t rand();
		
		/// Get a random number as double
		/*
			/param a Lower bound  -- \b IN.
			/param b Upper bound  -- \b IN.
		    /return A double in the interval [a,b].
		 */		
		static double rand(double a, double b){return a+((b-a)*NOMAD::RNG::rand())/UINT32_MAX;}
		
	private:
		static uint32_t x,y,z;  ///< Default parameter value for the random seed generator.
	};
}

#endif
