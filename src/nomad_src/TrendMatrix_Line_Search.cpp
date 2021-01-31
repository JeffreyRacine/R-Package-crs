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
 \file   TrendMatrix_Line_Search.cpp
 \brief  Trend matrix basic line search (implementation)
 \author Christophe Tribes
 \date   2017-05-3
 \see    TrendMatrix_Line_Search.hpp
 */
#include "TrendMatrix_Line_Search.hpp"

//bool option_random_in_cone = true;
bool option_random_in_cone = false;

/*----------------------------------------------------------------*/
/*                            the search                          */
/*----------------------------------------------------------------*/
/* Search parameters:                                             */
/* ------------------                                             */
/*  . Construct direction based on trend matrix around best_feas  */
/*    and best_infeas or just around best_feas;                   */
/*    default=around them both, not modifiable                    */
/*                                                                */
void NOMAD::TrendMatrix_Line_Search::search ( NOMAD::Mads              & mads           ,
                                             int                      & nb_search_pts  ,
                                             bool                     & stop           ,
                                             NOMAD::stop_type         & stop_reason    ,
                                             NOMAD::success_type      & success        ,
                                             bool                     & count_search   ,
                                             const NOMAD::Eval_Point *& new_feas_inc   ,
                                             const NOMAD::Eval_Point *& new_infeas_inc   )
{
    
    // TODO manages constant variables and variables groups
    // TODO make sure it is ok with binary, categorical, integer variables
    
    new_feas_inc  = new_infeas_inc = NULL;
    nb_search_pts = 0;
    success       = NOMAD::UNSUCCESSFUL;
    count_search  = false;
    
    _search_stats.reset();
    
    const NOMAD::Display    & out = _p.out();
    NOMAD::dd_type display_degree = out.get_search_dd();
    
    if ( stop )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "TrendMatrix_Line_Search::search(): not performed (stop flag is active)"
            << std::endl;
        return;
    }
    
    // black-box output types:
    const std::list<int>                     & index_obj_list = _p.get_index_obj();
    
    if ( index_obj_list.empty() )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "TrendMatrix_Line_Search::search(): not performed with no objective function"
            << std::endl;
        return;
    }
    
    // stats:
    NOMAD::Stats & stats = mads.get_stats();
    
    // initial displays:
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << NOMAD::TRENDMATRIX_LINE_SEARCH << " #"
        << stats.get_nb_TM_searches()+1;
        _out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    }

    // surrogate or truth model evaluations:
    NOMAD::eval_type ev_type =
    ( _p.get_opt_only_sgte() ) ? NOMAD::SGTE : NOMAD::TRUTH;
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << "Using " << ev_type << std::endl;
    }
    
    // active barrier:
    const NOMAD::Barrier & barrier =
    ( ev_type == NOMAD::SGTE ) ? mads.get_sgte_barrier() : mads.get_true_barrier();
    
    // current incumbents: xk[0]=best_feas and xk[1]=best_infeas:
    const NOMAD::Eval_Point * xk[2];
    xk[0] = barrier.get_best_feasible  ();
    xk[1] = barrier.get_best_infeasible();
    
    if ( !xk[0] && !xk[1] )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "end of " << NOMAD::TRENDMATRIX_LINE_SEARCH << " (no incumbent)";
            out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        return;
    }
    
    // from this point the search is counted:
    count_search = true;
    _search_stats.add_nb_TM_searches();
    
    NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
    
    // main loop on the two incumbents (feasible and infeasible):
    // ---------
    for ( int i_inc = 0 ; i_inc < 2 ; ++i_inc )
    {
        
        if ( xk[i_inc] )
        {
            
            // display the model center:
            if ( display_degree == NOMAD::FULL_DISPLAY )
            {
                out << std::endl << "Base point";
                if ( xk[0] && xk[1] )
                    out << " (" << i_inc+1 << "/2)";
                out << ": " << *xk[i_inc] << std::endl;
            }
            
            // get and check the signature:
            NOMAD::Signature * signature = xk[i_inc]->get_signature();
            if ( !signature )
            {
                if ( display_degree == NOMAD::FULL_DISPLAY )
                {
                    std::ostringstream oss;
                    oss << "end of " << NOMAD::MODEL_SEARCH << " (no signature)";
                    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
                }
                stats.update ( _search_stats , true);
                return;
            }
            
            // current mesh indices:
            NOMAD::Point mesh_indices = signature->get_mesh()->get_mesh_indices();
            
            
            int n = signature->get_n();
            if ( n != xk[i_inc]->size() )
            {
                if ( display_degree == NOMAD::FULL_DISPLAY )
                {
                    std::ostringstream oss;
                    oss << "end of " << NOMAD::TRENDMATRIX_LINE_SEARCH << " (incompatible signature)";
                    out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
                }
                stats.update   ( _search_stats , true );
                return;
            }
            
            // Set and get the trend direction
            ev_control.set_single_trend_direction ( *xk[i_inc] );
            const NOMAD::Point & trend_dir = ev_control.get_single_trend_direction();
            
            if ( trend_dir.is_defined() )
            {
                
                // compute the line search radius
                NOMAD::Point Delta , delta;
                signature->get_mesh()->get_Delta ( Delta );
                signature->get_mesh()->get_delta ( delta );
                
                
                NOMAD::Point line_search_radius = Delta;
                // line_search_radius *= _p.get_model_quad_radius_factor();
                
                if ( display_degree == NOMAD::FULL_DISPLAY )
                    out << "mesh indices         : ("   << mesh_indices  << ")"    << std::endl
                    << "mesh size parameter  : ( " << delta << " )" << std::endl
                    << "poll size parameter  : ( " << Delta << " )" << std::endl
                    << "line search radius   : ( " << line_search_radius << " )" << std::endl
                    << "trend direction      : ( " << trend_dir
                    << " )" << std::endl;
                
                NOMAD::Point search_dir;
                if ( option_random_in_cone )
                {
                    // Compute a random direction on a quasi N-Orthant
                    // Based on random direction on a N-sphere
                    //  see http://en.wikipedia.org/wiki/N-sphere
                    const int n = _p.get_dimension();
                    NOMAD::Point random_dir_in_cone( n );
                    for ( int i = 0 ; i < n ; ++i )
                    {
                        if ( trend_dir[i] == 0  )
                            random_dir_in_cone[i] = 0.0;
                        else
                        {
                            do random_dir_in_cone[i]=NOMAD::RNG::normal_rand(0,1);
                            while ( random_dir_in_cone[i] * trend_dir[i] < 0  );
                        }
                    }
                    if ( random_dir_in_cone.norm() > 0 )
                    {
                        random_dir_in_cone *= 1.0 /random_dir_in_cone.norm();
                        search_dir = -random_dir_in_cone ;
                    }
                    else
                       search_dir = -trend_dir;
                }
                else
                    search_dir = -trend_dir;
                    
                
                NOMAD::Point x_test = *xk[i_inc]+ search_dir*line_search_radius;
                
                
                // project to mesh (+round for integers), and create trial points:
                // ----------------------------------------------------------
                if ( x_test.is_defined() )
                    create_trial_point ( ev_control,
                                        x_test ,
                                        *xk[i_inc] ,
                                        *signature  ,
                                        mesh_indices  ,
                                        delta ,
                                        display_degree ,
                                        out              );
                
                
            }
        }
    } // end of main loop
        
    // evaluate the trial points:
    // --------------------------
    int bbe        = stats.get_bb_eval();
    // int sgte_eval  = stats.get_sgte_eval ();
    //int cache_hits = stats.get_cache_hits();
    
    new_feas_inc = new_infeas_inc = NULL;
    
    ev_control.disable_model_eval_sort();
    
    ev_control.eval_list_of_points ( _type                   ,
                                    mads.get_true_barrier() ,
                                    mads.get_sgte_barrier() ,
                                    mads.get_pareto_front() ,
                                    stop                    ,
                                    stop_reason             ,
                                    new_feas_inc            ,
                                    new_infeas_inc          ,
                                    success                   );
    
    ev_control.enable_model_eval_sort();
    
    // update stats:
    _search_stats.add_TM_bb_eval    ( stats.get_bb_eval   () - bbe        );
    
    if ( success == NOMAD::FULL_SUCCESS )
        _search_stats.add_TM_success();
    
    _search_stats.add_TM_pts ( nb_search_pts );
    
    // Update stats
    stats.update   ( _search_stats, true );
    
    // final display:
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "Trend matrix line search stats " ;
        _out << NOMAD::open_block ( oss.str() ) ;
        _out << " #pts: " << _search_stats.get_TM_pts() << std::endl ;

        _out << " #bb_eval: " << _search_stats.get_TM_bb_eval() << std::endl ;
        _out << " #success: " << _search_stats.get_TM_success() << std::endl;
        _out << NOMAD::close_block ( ) << std::endl;
        
        oss.str("");
        oss.clear();
        oss << "end of " << NOMAD::TRENDMATRIX_LINE_SEARCH << " (" << success << ")" << std::endl ;
        _out <<  NOMAD::close_block ( oss.str() ) << std::endl;
    }
}

/*---------------------------------------------------------------*/
/*        project to mesh and create a trial point (private)     */
/*---------------------------------------------------------------*/
void NOMAD::TrendMatrix_Line_Search::create_trial_point ( NOMAD::Evaluator_Control & ev_control,
                                                         NOMAD::Point                x,
                                                         const NOMAD::Eval_Point   & poll_center,
                                                         NOMAD::Signature          & signature ,
                                                         const NOMAD::Point        & mesh_indices ,
                                                         const NOMAD::Point        & delta        ,
                                                         NOMAD::dd_type              display_degree ,
                                                         const NOMAD::Display     & out              )
{
    
    if ( display_degree == NOMAD::FULL_DISPLAY )
    {
        out << "candidate (before projection)";
        out << ": ( " << x << " )" << std::endl;
    }
    
    
    // model search point:
    int n = x.size();
    
    // projection to mesh:
    x.project_to_mesh ( poll_center , delta , _p.get_lb() , _p.get_ub() );
    if ( display_degree == NOMAD::FULL_DISPLAY )
        out << "candidate (after projection) : ( "
        << x << " )" << std::endl;
    
    
    // Round for integer and binary variables:
    bool has_integer=false;
    bool has_binary=false;
    for (int i=0;i<n;i++)
    {
        if ( _p.get_bb_input_type()[i] == NOMAD::INTEGER )
        {
            has_integer=true;
            if ( x[i] >= 0.0 )
                x[i] = x[i].NOMAD::Double::ceil();
            else
                x[i] = x[i].NOMAD::Double::floor();
        }
        // binary variables:
        else if ( _p.get_bb_input_type()[i] == NOMAD::BINARY )
        {
            has_binary=true;
            if ( x[i] < 0.5 )
                x[i] = 0.0;
            else
                x[i] = 1.0;
            
        }
    }
    if ( has_integer && display_degree == NOMAD::FULL_DISPLAY )
        out << "candidate (after rounding integer) : ( "
        << x << " )" << std::endl;
    
    if ( has_binary && display_degree == NOMAD::FULL_DISPLAY )
        out << "candidate (after rounding binary) : ( "
        << x << " )" << std::endl;
    
    
    // compare x and center:
    if ( x == poll_center )
    {
        if ( display_degree == NOMAD::FULL_DISPLAY )
            out << "candidate rejected (candidate==model center)" << std::endl;
        return;
    }
    
    
    
    NOMAD::Eval_Point * tk = new NOMAD::Eval_Point;
    
    // if the search is optimistic, a direction is computed (this
    // will be used in case of success in the speculative search):
    if ( _p.get_model_search_optimistic() )
    {
        NOMAD::Direction dir ( n , 0.0 , NOMAD::TRENDMATRIX_SEARCH_DIR );
        dir.Point::operator = ( x - poll_center );
        tk->set_direction  ( &dir );
    }
    
    tk->set ( n , _p.get_bb_nb_outputs() );
    tk->set_signature  ( &signature  );
    tk->Point::operator = ( x );
    
    
    // we check that the candidate does not correspond to another candidate:
    const std::set<NOMAD::Priority_Eval_Point> & eval_lop
    = ev_control.get_eval_lop();
    
    std::set<NOMAD::Priority_Eval_Point>::const_iterator it , end = eval_lop.end();
    
    bool accept_point = true;
    for ( it = eval_lop.begin() ; it != end ; ++it )
        if ( it->get_point()->NOMAD::Point::operator == ( *tk ) )
        {
            accept_point = false;
            break;
        }
    
    if ( accept_point )
        // add the new point to the list of search trial points:
        ev_control.add_eval_point ( tk                     ,
                                   display_degree          ,
                                   _p.get_snap_to_bounds() ,
                                   NOMAD::Double()         ,
                                   NOMAD::Double()         ,
                                   NOMAD::Double()         ,
                                   NOMAD::Double()         );
    
}
