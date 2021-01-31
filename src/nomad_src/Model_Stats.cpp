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
 \file   Model_Stats.cpp
 \brief  Model stats (implementation)
 \author Sebastien Le Digabel
 \date   2010-09-24
 \see    Model_Stats.hpp
 */
#include "Model_Stats.hpp"

/*---------------------------------------------------------*/
/*                     affectation operator                */
/*---------------------------------------------------------*/
NOMAD::Model_Stats & NOMAD::Model_Stats::operator = ( const NOMAD::Model_Stats & s )
{
    _nb_truth            = s._nb_truth;
    _nb_sgte             = s._nb_sgte;
    _nb_MFN              = s._nb_MFN;
    _nb_WP_regression    = s._nb_WP_regression;
    _nb_regression       = s._nb_regression;
    _not_enough_pts      = s._not_enough_pts;
    _nb_Y_sets           = s._nb_Y_sets;
    _sum_nY              = s._sum_nY;
    _min_nY              = s._min_nY;
    _max_nY              = s._max_nY;
    _construction_error  = s._construction_error;
    _construction_time   = s._construction_time;
    _optimization_time   = s._optimization_time;
    _bad_cond            = s._bad_cond;
    MS_MS_nb_searches      = s.MS_MS_nb_searches;
    MS_MS_pts              = s.MS_MS_pts;
    MS_MS_success          = s.MS_MS_success;
    MS_MS_bb_eval          = s.MS_MS_bb_eval;
    MS_MS_sgte_eval        = s.MS_MS_sgte_eval;
    MS_MS_cache_hits       = s.MS_MS_cache_hits;
    MS_MS_rejected         = s.MS_MS_rejected;
    MS_MS_max_search_pts   = s.MS_MS_max_search_pts;
    MS_MS_nb_opt           = s.MS_MS_nb_opt;
    MS_MS_opt_error        = s.MS_MS_nb_opt;
    MS_MS_avg_model_eval   = s.MS_MS_avg_model_eval;
    MS_MS_max_model_eval   = s.MS_MS_max_model_eval;
    MS_MS_max_bbe          = s.MS_MS_max_bbe;
    ES_ES_nb_inside_radius = s.ES_ES_nb_inside_radius;
    ES_ES_nb_pts           = s.ES_ES_nb_pts;
    
    return *this;
}

/*---------------------------------------------------------*/
/*                      reset the stats                    */
/*---------------------------------------------------------*/
void NOMAD::Model_Stats::reset ( void )
{
    _nb_truth              =
    _nb_sgte             =
    _nb_MFN              =
    _nb_WP_regression    =
    _nb_regression       =
    _not_enough_pts      =
    _construction_error  =
    _bad_cond            =
    MS_MS_nb_searches      =
    MS_MS_pts              =
    MS_MS_success          =
    MS_MS_bb_eval          =
    MS_MS_sgte_eval        =
    MS_MS_cache_hits       =
    MS_MS_rejected         =
    MS_MS_max_search_pts   =
    MS_MS_nb_opt           =
    MS_MS_opt_error        =
    MS_MS_avg_model_eval   =
    MS_MS_max_model_eval   =
    MS_MS_max_bbe          =
    ES_ES_nb_inside_radius =
    ES_ES_nb_pts           =
    _nb_Y_sets           = 0;
    
    _min_nY = INT_MAX;
    _sum_nY = 0.0;
    
    _max_nY                = -1;
    _construction_time     =
    _optimization_time   = 0.0;
}

/*-------------------------------------------------------*/
/*  update stats from another NOMAD::Model_Stats object  */
/*-------------------------------------------------------*/
void NOMAD::Model_Stats::update ( const NOMAD::Model_Stats & s )
{
    _nb_truth            += s._nb_truth;
    _nb_sgte             += s._nb_sgte;
    _nb_MFN              += s._nb_MFN;
    _nb_WP_regression    += s._nb_WP_regression;
    _nb_regression       += s._nb_regression;
    _not_enough_pts      += s._not_enough_pts;
    _construction_error  += s._construction_error;
    _construction_time   += s._construction_time;
    _optimization_time   += s._optimization_time;
    _bad_cond            += s._bad_cond;
    MS_MS_nb_searches      += s.MS_MS_nb_searches;
    MS_MS_pts              += s.MS_MS_pts;
    MS_MS_success          += s.MS_MS_success;
    MS_MS_bb_eval          += s.MS_MS_bb_eval;
    MS_MS_sgte_eval        += s.MS_MS_sgte_eval;
    MS_MS_cache_hits       += s.MS_MS_cache_hits;
    MS_MS_rejected         += s.MS_MS_rejected;
    MS_MS_opt_error        += s.MS_MS_opt_error;
    MS_MS_max_bbe          += s.MS_MS_max_bbe;
    ES_ES_nb_inside_radius += s.ES_ES_nb_inside_radius;
    ES_ES_nb_pts           += s.ES_ES_nb_pts;
    _nb_Y_sets           += s._nb_Y_sets;
    _sum_nY              += s._sum_nY;
    
    _min_nY = ( _min_nY < s._min_nY ) ?
    _min_nY : s._min_nY;
    
    _max_nY = ( _max_nY > s._max_nY ) ?
    _max_nY : s._max_nY;
    
    MS_MS_max_model_eval = ( MS_MS_max_model_eval > s.MS_MS_max_model_eval ) ?
    MS_MS_max_model_eval : s.MS_MS_max_model_eval;
    
    MS_MS_max_search_pts = ( MS_MS_max_search_pts > s.MS_MS_max_search_pts ) ?
    MS_MS_max_search_pts : s.MS_MS_max_search_pts;
    
    if ( MS_MS_nb_opt + s.MS_MS_nb_opt == 0 )
        MS_MS_avg_model_eval = 0;
    else
        MS_MS_avg_model_eval = (   MS_MS_avg_model_eval*  MS_MS_nb_opt +
                              s.MS_MS_avg_model_eval*s.MS_MS_nb_opt   ) /
        ( MS_MS_nb_opt + s.MS_MS_nb_opt           );
    
    MS_MS_nb_opt += s.MS_MS_nb_opt;
}

/*----------------------------------------------------*/
/*                  update stats on |Y|               */
/*----------------------------------------------------*/
void NOMAD::Model_Stats::update_nY ( int nY )
{
    ++_nb_Y_sets;
    _sum_nY += nY;
    if ( nY > _max_nY )
        _max_nY = nY;
    if ( nY < _min_nY )
        _min_nY = nY;
}

/*----------------------------------------------------*/
/*               update _MS_max_search_pts            */
/*----------------------------------------------------*/
void NOMAD::Model_Stats::update_MS_max_search_pts ( int sp )
{
    if ( sp > MS_MS_max_search_pts )
        MS_MS_max_search_pts = sp;
}

/*-----------------------------------------------------------------------*/
/*  update stats _MS_nb_opt, _MS_max_model_eval, and _MS_avg_model_eval  */
/*-----------------------------------------------------------------------*/
void NOMAD::Model_Stats::update_MS_model_opt ( int eval )
{
    if ( eval > MS_MS_max_model_eval )
        MS_MS_max_model_eval = eval;
    ++MS_MS_nb_opt;
    MS_MS_avg_model_eval = ( (MS_MS_avg_model_eval*(MS_MS_nb_opt-1)) + eval ) / MS_MS_nb_opt;
}

/*---------------------------------------------------------*/
/*                          display                        */
/*---------------------------------------------------------*/
void NOMAD::Model_Stats::display ( const NOMAD::Display & out ) const
{
    out << "number of models built       : "   << get_nb_models()     << std::endl;
    if ( _nb_sgte > 0 )
        out << "number of truth models       : " << _nb_truth           << std::endl
        << "number of surrogate models   : " << _nb_sgte            << std::endl;
    out << "number of MFN interpolations : "   << _nb_MFN             << std::endl
    << "number of WP regressions     : "   << _nb_WP_regression   << std::endl
    << "number of quadr. regressions : "   << _nb_regression      << std::endl
    << "number of construction errors: "   << _construction_error << std::endl
    << "number of bad cond numbers   : "   << _bad_cond           << std::endl
    << "number of too small Y sets   : "   << _not_enough_pts     << std::endl
    << "min Y size                   : ";
    if ( _min_nY != INT_MAX )
        out << _min_nY;
    else
        out << "-";
    out << std::endl
    << "max Y size                   : ";
    if ( _max_nY != -1 )
        out << _max_nY;
    else
        out << "-";
    out << std::endl
    << "avg Y size                   : ";
    if ( get_avg_nY() != 0.0 )
        out << get_avg_nY();
    else
        out << "-";
    out << std::endl
    << "construction CPU time (s)    : "   << _construction_time  << std::endl;
    if ( MS_MS_nb_searches > 0 )
    {
        out << NOMAD::open_block ( "model searches" )
        << "number of searches                 : "   << MS_MS_nb_searches << std::endl
        << "number of search successes         : "   << MS_MS_success     << std::endl
        << "number of search points            : "   << MS_MS_pts         << std::endl
        << "number of blackbox evaluations     : "   << MS_MS_bb_eval     << std::endl;
        if ( MS_MS_sgte_eval > 0 )
            out << "number of sgte evaluations         : " << MS_MS_sgte_eval       << std::endl;
        out << "number of cache hits               : "   << MS_MS_cache_hits      << std::endl
        << "number of rejected candidates      : "   << MS_MS_rejected        << std::endl
        << "max number of trial points         : "   << MS_MS_max_search_pts  << std::endl
        << "number of optimizations            : "   << MS_MS_nb_opt          << std::endl
        << "number of optimization errors      : "   << MS_MS_opt_error       << std::endl
        << "number of max_bbe stops            : "   << MS_MS_max_bbe         << std::endl
        << "max number of model evaluations    : "   << MS_MS_max_model_eval  << std::endl
        << "average number of model evaluations: "   << MS_MS_avg_model_eval  << std::endl
        << "optimization CPU time (s)          : "   << _optimization_time  << std::endl
        << NOMAD::close_block();
    }
    
    if ( ES_ES_nb_pts > 0 ) 
    {
        out << NOMAD::open_block ( "model ordering" )
        << "number of points considered   : " << ES_ES_nb_pts << std::endl
        << "number of points inside radius: " << ES_ES_nb_inside_radius
        << " (";
        NOMAD::Double(ES_ES_nb_inside_radius*100.0/ES_ES_nb_pts).display ( out , "%.0f" );
        out    << "%)" << std::endl
        << NOMAD::close_block();
    }
}
