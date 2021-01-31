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
 \file   Exception.hpp
 \brief  Custom class for exceptions (headers)
 \author Sebastien Le Digabel
 \date   2010-03-29
 \see    Exception.cpp
 */
#ifndef __NOMAD_EXCEPTION__
#define __NOMAD_EXCEPTION__

#include <sstream>

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
    
    /// Custom class for exceptions.
    /**
     NOMAD uses this type of exceptions.
     It indicates the file and line number at which a throw is made.
     
     \b Example
     
     \code
     throw NOMAD::Exception ( __FILE__ , __LINE__ , "an error message" );
     \endcode
     */
    class DLL_API Exception : public std::exception {
        
    private:
        
        mutable std::string _what;  ///< Error message.
        std::string         _file;  ///< File where the exception is thrown.
        int                 _line;  ///< Line number at which the exception is thrown.

    protected:
        std::string         _type_msg;  ///< Basic exception message indicating the type of exception
        
    public:
        
        /// Constructor.
        /**
         \param file A string corresponding to the file where the
         exception is thrown -- \b IN
         \param line An integer corresponding to the line number
         at which the exception is thrown -- \b IN.
         \param msg  A string corresponding to the error message -- \b IN.
         */
        Exception ( const std::string & file , int line , const std::string & msg )
        : _what ( msg  ) ,
        _file ( file ) ,
        _line ( line )  {
            _type_msg = "NOMAD::Exception thrown";
        }
        
        /// Destructor.
        virtual ~Exception ( void ) throw() {}
        
        /// Access to the error message.
        /**
         \return A string with the error message.
         */
        const char * what(void) const throw()
        {
            std::ostringstream oss;
            if (!_file.empty() || _line > 0)
                oss << _type_msg << " (" << _file << ", " << _line << ")";
            if (!_what.empty())
                oss << ": " << _what;
            _what = oss.str();
            return _what.c_str();
        }
    };
}

#endif
