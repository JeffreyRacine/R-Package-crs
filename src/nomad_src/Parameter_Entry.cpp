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
 \file   Parameter_Entry.cpp
 \brief  Parameter entry (implementation)
 \author Sebastien Le Digabel
 \date   2010-04-05
 \see    Parameter_Entry.hpp
 */
#include "Parameter_Entry.hpp"

/*-----------------------------------*/
/*  . constructor (from a string)    */
/*  . ignores all entries after '#'  */
/*-----------------------------------*/
NOMAD::Parameter_Entry::Parameter_Entry ( const std::string & entry           ,
                                         bool                remove_comments   )
: _ok                 ( false ) ,
_unique               ( true  ) ,
_next                 ( NULL  ) ,
_param_file           ( "" ),
_line                 ( 0 ) ,
_has_been_interpreted ( false )
{
    int                i , idst;
    std::string        s;
    
    std::istringstream in ( entry );
    in >> _name;
    
    if (_name.size()==0)
        return;
    
    if ( remove_comments && _name[0] == '#' )
        _name.clear();
    else
    {
        
        
        NOMAD::toupper ( _name );
        
        bool stats_file_name_read = false;
        
        while ( true )
        {
            in >> s;
            
            if ( in.fail() )
                break;
            
            // comment:
            if ( remove_comments && s[0]=='#' )
                break;
    
            // string with begining/ending quotes:
            bool had_quotes = false;
            if ( s[0] == '\"' || s[0] == '\'' )
            {
                
                had_quotes = true;
                char quote = s[0];
                
                s.erase ( s.begin() );
                
                if ( s[s.size()-1] == quote )
                    s.resize ( s.size() - 1 );
                
                else
                {
                    
                    std::string ss;
                    getline ( in , ss , quote );
                    
                    if ( in.fail() || !in.good())
                    {
                        _ok = false;
                        return;
                    }
                    
                    s = s + ss;
                }
            }
            
            // DISPLAY_STATS or STATS_FILE:
            if ( _name == "STATS_FILE" && !stats_file_name_read )
                stats_file_name_read = true;
            
            else if ( _name == "DISPLAY_STATS" || _name == "STATS_FILE" )
            {
                
                if ( had_quotes )
                {
                    _ok = false;
                    return;
                }
                
                std::string keyword , ss = s;
                bool        interpreted  = false;
                NOMAD::toupper ( ss );
                
                NOMAD::display_stats_type dst = NOMAD::DS_OBJ;
                while ( dst < NOMAD::DS_UNDEFINED )
                {
                    
                    keyword = NOMAD::Display::get_display_stats_keyword (dst);
                    
                    i = static_cast<int> ( ss.rfind ( keyword , ss.size()-1 ) );
                    
                    if ( i >= 0 )
                    {
                        std::string s1 = s.substr ( 0 , i );
                        std::string s2 = s.substr ( i+keyword.size() );
                        if ( !s1.empty() )
                            _values.push_back ( s1 );
                        _values.push_back ( keyword );
                        if ( !s2.empty() )
                            _values.push_back ( s2 );
                        _values.push_back ( std::string() );
                        interpreted = true;
                        break;
                    }
                    
                    // loop increment:
                    idst = dst;
                    ++idst;
                    dst = static_cast<NOMAD::display_stats_type> ( idst );
                }
                
                if ( !interpreted ) 
                {
                    _values.push_back ( s             );
                    _values.push_back ( std::string() );
                }
                
                continue;
            }
            
            // vector:
            if ( s.size() > 1 && ( s[0] == '[' || s[0] == '(' ) ) 
            {
                _values.push_back ( s[0]=='[' ? "[" : "(" );
                s.erase(s.begin());
            }
            int  sm1 = static_cast<int>(s.size()) - 1;
            char c   = s[sm1];
            if ( s.size() > 1 && ( c == ']' || c == ')' ) ) 
            {
                s.resize(sm1);
                _values.push_back (s);
                _values.push_back ( c==']' ? "]" : ")" );
                continue;
            }
            
            // other values:
            _values.push_back ( s );
        }
        
        if ( !_values.empty() )
            _ok = true;
    }
}

/*------------------------------*/
/*             display          */
/*------------------------------*/
void NOMAD::Parameter_Entry::display ( const NOMAD::Display & out ) const
{
    if ( _ok ) 
    {
        out << _name << ": ";
        std::list<std::string>::const_iterator end = _values.end();
        for ( std::list<std::string>::const_iterator it = _values.begin() ; it != end ; ++it )
            out << "[" << *it << "] ";
    }
}
