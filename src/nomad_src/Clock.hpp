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
 \file   Clock.hpp
 \brief  Clock class (headers)
 \author Sebastien Le Digabel
 \date   2010-04-02
 \see    Clock.cpp
 */
#ifndef __CLOCK__
#define __CLOCK__

#include <ctime>

#include "nomad_version.hpp"

#ifdef _MSC_VER
#pragma warning(disable:4275)
#pragma warning(disable:4251)
#ifdef DLL_EXPORTS
#define DLL_API __declspec(dllexport) 
#else
#define DLL_API __declspec(dllimport) 
#endif
#else
#define DLL_API
#endif


namespace NOMAD {
    
    /// Clock class.
    /**
     Time measurement.\n\n
     \b Example:
     \code
     Clock c;
     
     // some instructions here
     
     std::cout << "elapsed real time = " << c.get_real_time() << std::endl;
     std::cout << "elapsed CPU time  = " << c.get_CPU_time()  << std::endl;
     \endcode
     */
    class DLL_API Clock {
        
    private:
        
        time_t              _real_t0;          ///< Wall clock time measurement.
        clock_t             CPU_CPU_t0;           ///< CPU time measurement.  //zhenghua, for compiling on Solaris.
        static const double D_D_CLOCKS_PER_SEC; ///< System constant for CPU time measurement.       //zhenghua
        
    public:
        
        /// Constructor.
        Clock ( void ) : CPU_CPU_t0 ( clock() ) { time (&_real_t0); }
        
        /// Copy constructor.
        /**
         \param c The copied object -- \b IN.
         */
        Clock ( const Clock & c ) : _real_t0 ( c._real_t0 ) , CPU_CPU_t0 ( c.CPU_CPU_t0 ) {}
        
        /// Affectation operator.
        /**
         \param  c The right-hand side object -- \b IN.
         \return \c *this as the result of the affectation.
         */
        Clock & operator = ( const Clock & c )
        {
            _real_t0 = c._real_t0;
            CPU_CPU_t0  = c.CPU_CPU_t0;
            return *this;
        }
        
        /// Destructor.
        virtual ~Clock ( void ) {}
        
        /// Reset the clock.
        void reset ( void )
        {
            time ( &_real_t0 );
            CPU_CPU_t0 = clock();
        }
        
        /// Get wall clock time.
        /**
         \return The wall clock time.
         */
        int get_real_time ( void ) const;
        
        /// Get the CPU time.
        /**
         \return The CPU time.
         */
        double get_CPU_time ( void ) const
        {
            return ( clock() - CPU_CPU_t0 ) / D_D_CLOCKS_PER_SEC;
        }
    };
}

#endif
