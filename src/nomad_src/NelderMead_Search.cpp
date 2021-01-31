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
 \file   NelderMead_Search.cpp
 \brief  NelderMead search (implementation)
 \author Christophe Tribes
 \date   2017-03-31
 \see    NelderMead_Search.hpp
 */
#include "NelderMead_Search.hpp"


void NOMAD::NelderMead_Search::search (NOMAD::Mads              & mads           ,
                                       int                      & nb_search_pts  ,
                                       bool                     & stop           ,
                                       NOMAD::stop_type         & stop_reason    ,
                                       NOMAD::success_type      & success        ,
                                       bool                     & count_search   ,
                                       const NOMAD::Eval_Point *& new_feas_inc   ,
                                       const NOMAD::Eval_Point *& new_infeas_inc   )
{
    new_feas_inc  = new_infeas_inc = NULL;
    const NOMAD::Eval_Point * step_new_feas_inc = NULL;
    const NOMAD::Eval_Point * step_new_infeas_inc = NULL;
    
    nb_search_pts = 0;
    int nb_success = 0;
    success       = NOMAD::UNSUCCESSFUL;
    NOMAD::success_type step_success = NOMAD::UNSUCCESSFUL;
    count_search  = false;
    
    _search_stats.reset();
    
    if ( stop )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "NelderMead_Search::search(): not performed (stop flag is active)"
            << std::endl;
        return;
    }
    
    // black-box input types:
    const std::list<int> & index_obj_list = _p.get_index_obj();
    
    // active cache:
    const NOMAD::Cache & cache = mads.get_cache();
    
    // evaluator for the problem (used to update f and h from a cache)
    const NOMAD::Evaluator * ev = mads.get_evaluator_control().get_evaluator();
    
    // NM search cannot be done on sgte for the moment
    if ( cache.get_eval_type() != NOMAD::TRUTH )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "NelderMead_Search::search(): not performed on surrogates"
            << std::endl;
        return;
    }
    
    if ( index_obj_list.empty() )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "NelderMead_Search::search(): not performed without objective function"
            << std::endl;
        return;
    }
    
    if ( index_obj_list.size() >= 2 )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "NelderMead_Search::search(): not performed on multi-objective problem"
            << std::endl;
        
        return;
    }
    
    // stats:
    NOMAD::Stats & stats = mads.get_stats();
    
    // initial displays:
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << NOMAD::NM_SEARCH << " #"
        << stats.get_nb_NM_SEARCHES ()+1;
        _out << std::endl << NOMAD::open_block ( oss.str() ) << std::endl;
    }
    
    // surrogate or truth model evaluations:
    NOMAD::eval_type ev_type =
    ( _p.get_opt_only_sgte() ) ? NOMAD::SGTE : NOMAD::TRUTH;
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << "Search using " << ev_type << std::endl;
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
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "end of " << NOMAD::NM_SEARCH << " (no incumbent)";
            _out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        return;
    }
    
    // from this point the search is counted:
    count_search = true;
    _search_stats.add_nb_NM_searches();
    
    // display the number of cache points:
    if ( _display_degree == NOMAD::FULL_DISPLAY )
        _out << std::endl << "number of points in cache: "
        << cache.size() << std::endl;
    
    NOMAD::Evaluator_Control & ev_control = mads.get_evaluator_control();
    
    bool NM_step_ok =true;
    bool create_trial_ok =true;
    bool eval_point_ok = true;
    
    // main loop on the two incumbents (feasible and infeasible)
    // A single valid incumbent is considered for NM
    // (primary or secondary but most likely primary)
    int i_inc;
    // ---------
    for (  i_inc = 0 ; i_inc < 2 ; ++i_inc )
    {
        
        if ( ! xk[i_inc] ||
            xk[i_inc]->size() != _n ||
            ! xk[i_inc]->is_complete() )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
            {
                if ( i_inc == 0 )
                    _out << "No primary poll center. Let us try with secondary poll center." << std::endl;
                if ( i_inc == 1 )
                {
                    _out << "No secondary poll center." << endl;
                    std::ostringstream oss;
                    oss << "end of " << NOMAD::NM_SEARCH <<" (no poll center)"; ;
                    return;
                }
            }
            if ( i_inc == 1 )
                return;
        }
        else
            break;
    }
    
    
    
    // display the model center:
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << std::endl << ((i_inc == 0 )? "primary":"secondary") << " poll center";
        _out << ": " << *xk[i_inc] << std::endl;
    }
    
    // get and check the signature:
    NOMAD::Signature * signature = xk[i_inc]->get_signature();
    if ( !signature )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "end of " << NOMAD::NM_SEARCH << " (no signature)";
            _out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        stats.update( _search_stats , true);
        return;
    }
    
    // Initialize NM step
    // -----------------
    stop = false;
    _step = INITIAL;
    
        
    // Pointer to the new point proposed by NM (except SHRINK phase)
    const NOMAD::Eval_Point * xr = NULL ;
    const NOMAD::Eval_Point * xe = NULL ;
    const NOMAD::Eval_Point * xic = NULL ;
    const NOMAD::Eval_Point * xoc = NULL ;
    
    // Pointer to a generic point of NM (will point to *xr, or *xe or *xic or *xoc )
    const NOMAD::Eval_Point * xt = NULL ;
    
    bool nm_stop = false;
    NM_stop_type nm_stop_reason =UNDEFINED_STEP;
    
    // Loop NelderMead steps
    // ---------------------
    while ( !nm_stop && ! stop )
    {
        
        
        // Perform a NM step:
        // ------------------
        NM_step_ok = NM_step ( cache         ,
                              *ev            ,
                              xk[i_inc]      ,
                              nm_stop        ,
                              nm_stop_reason );
        
        
        // count NM step even if NM_step not ok or if point not submitted
        if ( !stop && ! nm_stop )
            _search_stats.add_NM_step();
            
        
        
        // Create and evaluate a single trial point or multiple trial points (SHRINK step only)
        // ------------------------------------------------------------------------------------
        if ( !stop && !nm_stop && NM_step_ok && _nm_submitted_points.size() > 0 )
        {
            
            std::list<NOMAD::Eval_Point *>::iterator it_list = _nm_submitted_points.begin() ;
            
            for ( ; it_list != _nm_submitted_points.end() ; ++it_list )
            {
                // Create trial point and project to mesh (+round for integers):
                create_trial_ok = create_trial_point ( ev_control    ,
                                                      *it_list       ,
                                                      *xk[i_inc]      );
                
                
                if ( ! create_trial_ok )
                    break;
            }
            
            // Points submitted by NM step may not become trial points --> they are not inserted in the cache and will not be deleted at the end of Nomad --> manage delete befor continuing
            if ( ! create_trial_ok )
            {
                for ( it_list= _nm_submitted_points.begin() ; it_list != _nm_submitted_points.end() ; ++it_list )
                    delete *it_list;
                
                _nm_submitted_points.clear();
            }
            
            if ( _step == SHRINK && ! create_trial_ok )
            {
                nm_stop = true;
                nm_stop_reason = SHRINK_FAILED;
            }
            
            
            int bbe        = stats.get_bb_eval();
            // int sgte_eval  = stats.get_sgte_eval (); // Not yet
            int cache_hits = stats.get_cache_hits();
            if ( create_trial_ok )
            {
                step_new_feas_inc = step_new_infeas_inc = NULL;
                step_success = NOMAD::UNSUCCESSFUL ;
                
                
                ev_control.eval_list_of_points ( _type                   ,
                                                mads.get_true_barrier() ,
                                                mads.get_sgte_barrier() ,
                                                mads.get_pareto_front() ,
                                                stop                    ,
                                                stop_reason             ,
                                                step_new_feas_inc       ,
                                                step_new_infeas_inc     ,
                                                step_success            ,
                                                & _nm_evaluated_points     );
                
                // Clear the list of submitted nm points
                _nm_submitted_points.clear();
                
                // Check evaluation ok
                if ( _nm_evaluated_points.size() != 0 )
                {
                    xt= *(_nm_evaluated_points.begin());
                    
                    if ( ! xt->is_in_cache() || ! xt->is_eval_ok() )
                        eval_point_ok = false;
                    else
                        eval_point_ok = true ;
                }
                else
                    eval_point_ok  = false;
                
                // Update the success status and the incumbent
                if ( step_success > NOMAD::UNSUCCESSFUL )
                {
                    nb_success++;
                    if ( _p.get_NM_search_opportunistic() )
                    {
                        nm_stop = true;
                        nm_stop_reason = OPPORTUNISTIC_STOP;
                        
                    }
                    if ( step_success >= success )
                    {
                        success = step_success;
                        new_feas_inc = step_new_feas_inc;
                        new_infeas_inc = step_new_infeas_inc ;
                    }
                }
                
                if ( success == NOMAD::FULL_SUCCESS )
                    _search_stats.add_NM_success();
                
                if ( eval_point_ok )
                    nb_search_pts+=(int)_nm_evaluated_points.size();
            
                // update stats:
                _search_stats.add_NM_bb_eval    ( stats.get_bb_eval   () - bbe        );
                // _search_stats.add_NM_sgte_eval  ( stats.get_sgte_eval () - sgte_eval  );  // NM+sgte not yet available
                _search_stats.add_cache_hits ( stats.get_cache_hits() - cache_hits );
                _search_stats.add_NM_pts ( nb_search_pts );
            
                if ( _p.get_NM_search_max_trial_pts() > 0 && nb_search_pts > _p.get_NM_search_max_trial_pts() )
                {
                    nm_stop = true;
                    nm_stop_reason = MAX_SEARCH_POINTS_REACHED;
                    if ( _display_degree == NOMAD::FULL_DISPLAY )
                    {
                        _out << "the maximum number of NM trial points ("
                        << _p.get_NM_search_max_trial_pts() << ") has been reached";
                        _out << std::endl;
                    }
                    break;
                }
                int factor = _p.get_NM_search_max_trial_pts_nfactor();
                if ( factor > 0 && nb_search_pts > factor * _p.get_nb_free_variables() )
                {
                    nm_stop = true;
                    nm_stop_reason = MAX_SEARCH_POINTS_REACHED;
                    if ( _display_degree == NOMAD::FULL_DISPLAY )
                    {
                        _out << "the maximum number of NM trial points ("
                        << factor * _p.get_nb_free_variables() << ") has been reached";
                        _out << std::endl;
                    }
                    break;
                }
                
                if ( _p.get_NM_search_min_simplex_vol().value() > 0 && _simplex_von > 0 && _simplex_von < _p.get_NM_search_min_simplex_vol().value() )
                {
                    nm_stop = true;
                    nm_stop_reason = MIN_SIMPLEX_VOL_REACHED;
                    if ( _display_degree == NOMAD::FULL_DISPLAY )
                    {
                        _out << "the minimum simplex volume ("
                        << _p.get_NM_search_min_simplex_vol() << ") has been reached";
                        _out << std::endl;
                    }
                    break;
                }
            }
        }
        
        // Select next NM step
        // -------------------
        if ( !stop && ! nm_stop )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
            {
                _out << "Completed NM " ;
                switch ( _step )
                {
                    case INITIAL:
                        _search_stats.add_NM_initial_step();
                        _out << "initial step" <<endl;
                        break;
                    case REFLECT:
                        _search_stats.add_NM_reflect_step();
                        _out << "reflect step" <<endl;
                        break;
                    case EXPAND:
                        _search_stats.add_NM_expand_step();
                        _out << "expand step" <<endl;
                        break;
                    case OUTSIDE_CONTRACTION:
                        _search_stats.add_NM_outside_contraction_step();
                        _out << "outside contraction step"<<endl;
                        break;
                    case INSIDE_CONTRACTION:
                        _search_stats.add_NM_inside_contraction_step();
                        _out << "inside contraction step"<<endl;
                        break;
                    case SHRINK:
                        _search_stats.add_NM_shrink_step();
                        _out << "shrink step"<<endl;
                        break;
                    default:
                        _out << "unkown step"<<endl;
                        break;
                }
            }
            
            switch ( _step )
            {
                case INITIAL:
                    if ( NM_step_ok )
                        _step = REFLECT;
                    else
                    {
                        nm_stop = true;
                        nm_stop_reason = STEP_FAILED;
                    }
                    break;
                case REFLECT:
                    if ( NM_step_ok && create_trial_ok && eval_point_ok )
                    {
                        xr = xt;
                        if ( ! xr || ! xr->is_defined() )
                            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                    "NelderMead_Search::search(): could not get the reflected point" );
                        
                        if ( point_dominates_Y0( *xr ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The reflect point xr: "
                                << *xr << " dominates Y0. Next, perform Expand.";
                                _out << std::endl << std::endl;
                            }
                            _step = EXPAND;
                        }
                        else if ( Yn_dominates_point( *xr ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The reflect point xr: "
                                << *xr << " is dominated by Yn. Next, perform Inside Contraction.";
                                _out << std::endl << std::endl;
                            }
                            _step = INSIDE_CONTRACTION;
                        }
                        else if ( point_dominates_pts_in_Y( *xr , 2 ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The reflect point xr: " << std::endl
                                << *xr << " dominates at least 2 points of Y." << std::endl << std::endl ;
                            }
                            if ( insert_in_Y ( *xr ) )
                            {
                                if ( _display_degree == NOMAD::FULL_DISPLAY )
                                {
                                    _out << "Insertion in Y is successful. Next perform Reflect." << std::endl << std::endl;
                                }
                                _step = REFLECT;
                            }
                            else
                            {
                                if ( _display_degree == NOMAD::FULL_DISPLAY )
                                {
                                    _out << " Cannot insert xr in Y. Stop." << std::endl;
                                }
                                
                                // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
                                nm_stop = true;
                                nm_stop_reason = INSERTION_FAILED;
                            }
                        }
                        else if ( point_dominates_pts_in_Y( *xr , 1 ) || point_dominates_pts_in_Y( *xr , 0 ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The reflect point xr:" << std::endl
                                << *xr << " dominates 1 or 0 point of Y. Next, perform Outside Contraction.";
                                _out << std::endl;
                            }
                            _step = OUTSIDE_CONTRACTION;
                        }
                        else
                        {
                            _step = ( _perform_shrink ) ? SHRINK:COMPLETE ;
                            nm_stop = ( _perform_shrink ) ? false:true ;
                            if ( nm_stop )
                                nm_stop_reason = COMPLETED ;
                            else
                                throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                        "NelderMead_Search::search(): point dominates less than 0 points in Y!" );
                            
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The reflect point is xr:" << std::endl
                                << *xr << " Next, perform Shrink.";
                                _out << std::endl;
                            }
                            
                        }
                    }
                    else
                    {
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "Cannot create a proper reflect point xr. Next, perform Inside Contraction.";
                            _out << std::endl;
                        }
                        _step = INSIDE_CONTRACTION;
                    }
                    break;
                case EXPAND:
                    xe = xt;
                    if ( eval_point_ok && (! xe || ! xe->is_defined() ) )
                        throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                "NelderMead_Search::search(): could not get the expanded point" );
                    if ( eval_point_ok && point_dominates_Y0( *xe ) )
                    {
                        if ( ! xr || ! xr->is_defined() )
                            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                    "NelderMead_Search::search(): could not get the reflected point" );
                        
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "The expand point xe:" << std::endl << *xe << " dominates Y0. ";
                        }
                        
                        if ( insert_in_Y_best( *xr , *xe ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Insert in Y the best of xr and xe. Next perform Reflect. " << std::endl << std::endl;
                            }
                            _step = REFLECT;
                        }
                        else
                        {
                            
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The insertion in Y of the best of xr and xe dit not maintain a proper Y. Stop NM." << std::endl;
                            }
                            
                            // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
                            nm_stop = true;
                            nm_stop_reason = INSERTION_FAILED;
                        }
                    }
                    else
                    {
                        if ( ! xr || ! xr->is_defined() )
                            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                    "NelderMead_Search::search(): could not get the reflected point" );
                        
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "The expand point xe:" << std::endl << *xe << " do not dominates Y0 or cannot be evaluated. " << std::endl;
                        }
                        
                        if ( insert_in_Y( *xr ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Insert xr:" << std::endl << *xr << " in Y. Next perform Reflect. " << std::endl << std::endl ;
                            }
                            _step = REFLECT;
                        }
                        else
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Cannot insert xr:" << std::endl << *xr << " in Y. Stop." << std::endl;
                            }
                            
                            // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
                            nm_stop = true;
                            nm_stop_reason = INSERTION_FAILED;
                        }
                    }
                    break;
                case INSIDE_CONTRACTION:
                    xic = xt;
                    if ( ! create_trial_ok || ! eval_point_ok
                        || Yn_dominates_point( *xic ) )
                    {
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "Cannot create or evaluate xic or Yn dominates xic. " ;
                        }
                        
                        _step = ( _perform_shrink ) ? SHRINK:COMPLETE ;
                        nm_stop = ( _perform_shrink ) ? false:true ;
                        if ( nm_stop )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Stop NM." << std::endl << std::endl;
                            }
                            nm_stop_reason = COMPLETED ;
                        }
                        else if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "Next perform Shrink." << std::endl << std::endl;
                        }
                    }
                    else
                    {
                        if ( ! xic || ! xic->is_defined() )
                            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                    "NelderMead_Search::search(): could not get the interior contraction point" );
                        
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "The inside contraction point xic:" << std::endl
                            << *xic << " is not dominated by Yn." << std::endl << std::endl;
                        }
                        
                        if ( insert_in_Y( *xic ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Next perform Reflect." << std::endl << std::endl;
                            }
                            _step = REFLECT;
                        }
                        else
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The insert of xic in Y dit not maintain a proper Y. Stop NM." << std::endl << std::endl;
                            }
                            
                            // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
                            nm_stop = true;
                            nm_stop_reason = INSERTION_FAILED;
                        }
                    }
                    break;
                case OUTSIDE_CONTRACTION:
                    if ( NM_step_ok && create_trial_ok && eval_point_ok )
                    {
                        xoc = xt;
                        if ( ! xoc || ! xoc->is_defined() )
                            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                    "NelderMead_Search::search(): could not get the outside contraction point" );
                        if ( ! xr || ! xr->is_defined() )
                            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                    "NelderMead_Search::search(): could not get the reflected point" );
                        
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "The outside contraction point is xoc:" << std::endl << *xoc ;
                        }
                        
                        if ( insert_in_Y_best( *xr , *xoc ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The insertion of the best of xr and xoc in Y is valid. Next perform Reflect." << std::endl << std::endl;
                            }
                            _step = REFLECT;
                        }
                        else
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "The insertion of the best of xr and xoc in Y did not maintain a proper Y. Stop NM." << std::endl << std::endl ;
                            }
                            // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
                            nm_stop = true;
                            nm_stop_reason = INSERTION_FAILED;
                        }
                    }
                    else
                    {
                        if ( ! xr || ! xr->is_defined() )
                            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                                    "NelderMead_Search::search(): could not get the reflected point" );
                        
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                        {
                            _out << "Cannot create or evaluate xic.";
                        }
                        
                        if ( insert_in_Y( *xr ) )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Reflect point xr is successfully inserted in Y. Next perform Reflect." << std::endl << std::endl ;
                            }
                            _step = REFLECT;
                        }
                        else
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Insertion of reflect point did not maintain a proper Y. Stop NM." << std::endl;
                            }
                            
                            // The insertion did not maintain a proper Y (insufficient rank or failed insertion)
                            nm_stop = true;
                            nm_stop_reason = INSERTION_FAILED;
                        }
                    }
                    break;
                case SHRINK:
                    if ( _perform_shrink )
                    {
                        create_initial_sets_from_new_points (cache , nm_stop , nm_stop_reason );
                        if ( !nm_stop )
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Successfull creation of a new Shrinked simplex. Next perform Reflect." << std::endl << std::endl;
                            }
                            _step = REFLECT;
                        }
                        
                        else
                        {
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                            {
                                _out << "Cannot create a new Shrinked simplex. Stop NM." << std::endl;
                            }
                            nm_stop = true;
                            nm_stop_reason = SHRINK_FAILED;
                        }
                    }
                    else
                    {
                        nm_stop = true;
                        nm_stop_reason  = SHRINK_REQUEST;
                    }
                    break;
                default:
                    nm_stop = true;
                    nm_stop_reason = UNDEFINED_STEP;
                    break;
            }
        }
        
        // Empty the list of evaluated points at each NM step
        _nm_evaluated_points.clear();
        
        if ( nm_stop && _display_degree == FULL_DISPLAY )
        {
            if ( nm_stop_reason == INSERTION_FAILED )
                _out << "Stop because insertion of point in Y failed." << std::endl;
            else if ( nm_stop_reason == TOO_SMALL_SIMPLEX )
                _out << "Stop because simplex Y is too small." << std::endl;
            else if ( nm_stop_reason == SIMPLEX_RANK_INSUFFICIENT )
                _out << "Stop because rank of Dz[(y1-y0) (y2-y0)...(yn-y0)] is too small." << std::endl;
            else if ( nm_stop_reason == SHRINK_FAILED )
                _out << "Stop because insertion of shrink of Y failed." << std::endl;
            else if ( nm_stop_reason == COMPLETED )
                _out << "Stop because NM is completed." << std::endl;
            
        }
        if ( nm_stop_reason == STEP_FAILED )
        {
            throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                    "NelderMead_Search::search(): NM step failed");
        }
        
    }
    
    // Empty the set of points for the next NM
    _nm_Y.clear();
    
    // Update stats
    stats.update   ( _search_stats, true );
    
    // check if NM worked:
    if ( ! NM_step_ok )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "end of " << NOMAD::MODEL_SEARCH
            << " ( method error)";
            _out << std::endl << NOMAD::close_block ( oss.str() ) << std::endl;
        }
        return;
    }
    
    // final display:
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "Nelder-Mead stats " ;
        _out << NOMAD::open_block ( oss.str() ) ;
        _out << " #bb_eval: " << _search_stats.get_NM_bb_eval() << std::endl ;
        _out << " #success: " << _search_stats.get_NM_success() << std::endl;
        _out << " #expand step: " << _search_stats.get_NM_expand_step() << std::endl;
        _out << " #reflect_step: " << _search_stats.get_NM_reflect_step() << std::endl;
        _out << " #inside_contraction_step: " << _search_stats.get_NM_inside_contraction_step() << std::endl;
        _out << " #outside_contraction_step: " << _search_stats.get_NM_outside_contraction_step() << std::endl;
        _out << " #shrink_step: " << _search_stats.get_NM_shrink_step() << std::endl;
        _out << NOMAD::close_block ( ) << std::endl;

        oss.str("");
        oss.clear();
        oss << "end of " << NOMAD::NM_SEARCH << " (" << success << ")" << std::endl ;
        _out <<  NOMAD::close_block ( oss.str() ) << std::endl;
    }
}


/*--------------------------------------------------------*/
/*  Create Nelder-Mead simplex from NM evaluated points   */
/*--------------------------------------------------------*/
void NOMAD::NelderMead_Search::create_initial_sets_from_new_points (const NOMAD::Cache           & cache,
                                                                    bool                         & stop,
                                                                    NM_stop_type                 & stop_reason )
{
    
    
    // Clear the list of NM points
    _nm_Y.clear();
    

    size_t min_Y_size = _n_free+1;
    
    // Browse the new points.
    // Construct the new simplex from cache points
    std::list<const NOMAD::Eval_Point *>::iterator it;
    for ( it= _nm_evaluated_points.begin() ; it != _nm_evaluated_points.end() ; it++ )
    {
        const NOMAD::Eval_Point * cur = *it ;
        
        if ( cur->get_eval_status() == NOMAD::EVAL_OK &&
            cur->get_n          () == _n )
        {
            
            const NOMAD::Point & bbo_cur = cur->get_bb_outputs();
            
            if ( bbo_cur.is_complete() )
            {
                NOMAD::NelderMead_Simplex_Eval_Point Y ( cur );
                std::pair<std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator,bool> ret = _nm_Y.insert ( Y );
                
                if ( ! ret.second )
                {
                    stop_reason = INITIAL_FAILED;
                    stop = true;
                    if ( _display_degree == NOMAD::FULL_DISPLAY )
                        _out << "Stop NM because cannot insert a point in Y." << std::endl;
                    break;
                }
            }
            else
            {
                stop_reason = INITIAL_FAILED;
                stop = true;
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                    _out << "Stop NM because cannot insert a point in Y." << std::endl;
                break;
            }
        }
    }
    
    _nm_submitted_points.clear();
    
    if ( stop )
        return;
    
    
    // Update simplex characteristics (volumes and diameter)
    update_Y_characteristics();
    
    // Create Y0 and Yn
    make_list_Y0(stop,stop_reason);
    if ( stop )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "Cannot make Y0." << std::endl;
        return;
    }
    
    make_list_Yn(stop,stop_reason);
    if ( stop )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "Cannot make Yn." << std::endl;
        return;
    }
    
    display_Y_info();
    
    // not enough points (this is not counted as an error):
    if ( _nm_Y.size() < min_Y_size )
    {
        stop_reason = INITIAL_FAILED;
        stop = true;
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "Stop NM because not enough points in Y." << std::endl;
        return;
    }
}



/*----------------------------------------------------------------------------------*/
/*  Create initial sets of points for Nelder-Mead within a radius of current best   */
/*----------------------------------------------------------------------------------*/
void NOMAD::NelderMead_Search::create_initial_sets_from_cache (const NOMAD::Cache           & cache,
                                                               const NOMAD::Evaluator       & ev,
                                                               const NOMAD::Eval_Point      * xk,
                                                               bool                         & stop ,
                                                               NM_stop_type                 & stop_reason )
{
    // Clear the list of NM points
    _nm_Y.clear();
    
    
    size_t min_Y_size = _n_free + 1;
    
    // get and check the signature:
    NOMAD::Signature * signature = xk->get_signature();
    
    // current mesh indices:
    NOMAD::Point mesh_indices = signature->get_mesh()->get_mesh_indices();
    
    // compute the include radius: points in Y must be at
    // a max distance of ms_radius_factor times Delta^k:
    NOMAD::Point Delta = signature->get_mesh()->get_Delta ( );
    NOMAD::Point delta = signature->get_mesh()->get_delta ( );
    
    NOMAD::Point include_rectangle = Delta;
    
    include_rectangle *= _p.get_NM_search_include_factor();
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
        _out << "mesh indices           : ("   << mesh_indices  << ")"    << std::endl
        << "mesh size parameter    : ( " << delta << " )" << std::endl
        << "poll size parameter    : ( " << Delta << " )" << std::endl
        << "include rectangle size : ( " << include_rectangle
        << " )" << std::endl;
    
    
    const NOMAD::Point & bbo_xk  = xk->get_bb_outputs();
    int m = bbo_xk.size();
    
    // The set of points initially included
    std::set<NOMAD::NelderMead_Simplex_Eval_Point> T;
    
    // Set h_min
    NOMAD::NelderMead_Simplex_Eval_Point::set_h_min ( _p.get_h_min() );

    
#ifdef DEBUG
    {
        std::ostringstream oss;
        oss << "Insertion of potential points to include in initial Y: ";
        _out << open_block( oss.str() ) ;
    }
#endif
    
    // browse the cache:
    const NOMAD::Eval_Point * cur = cache.begin();
    while ( cur )
    {
        
        if ( cur->get_eval_status() == NOMAD::EVAL_OK &&
            cur->get_n          () == _n &&
            signature->is_compatible (*xk)             )
        {
            
            const NOMAD::Point & bbo_cur = cur->get_bb_outputs();
            
            if ( bbo_cur.is_complete() && check_outputs(bbo_cur, m ) )
            {
                // the center point has been found and put in list
                if ( *xk == *cur )
                {
                    NOMAD::NelderMead_Simplex_Eval_Point Y ( cur );
                    T.insert ( Y );
#ifdef DEBUG
                    _out << *(Y.get_element()) << endl;
#endif
                }
                // other points must be within the include region:
                else
                {
                    bool include = true;
                    for (int i = 0 ; i < _n ; i++ )
                    {
                        NOMAD::Double val = (*cur)[i] - (*xk)[i];
                        if ( val.abs() > include_rectangle[i]  )
                        {
                            include = false;
                            break;
                        }
                    }
                    if ( include )
                    {
                        
                        // Make sure to evaluate f or h for points in cache (important if cache is loaded from file)
                        ev.compute_f( NOMAD::Cache::get_modifiable_point ( *cur ) );
                        ev.compute_h( NOMAD::Cache::get_modifiable_point ( *cur ) );
                        
                        NOMAD::NelderMead_Simplex_Eval_Point Y ( cur );
                        std::pair<std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator,bool> ret = T.insert ( Y );
                        
#ifdef DEBUG
                        _out << *(Y.get_element()) << endl;
#endif
                        
                        
                        if ( ! ret.second )
                        {
                            stop_reason = INITIAL_FAILED;
                            stop = true;
                            if ( _display_degree == NOMAD::FULL_DISPLAY )
                                _out << "Cannot insert a point in Y." << std::endl;
                            return;
                        }
                    }
                }
            }
        }
        cur = cache.next();
    }
        
#ifdef DEBUG
    _out << close_block() << endl;
#endif
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << "Number of potential points to include in initial Y: ";
        _out << T.size() << std::endl;
    }
    
    // not enough points (this is not counted as an error):
    if ( T.size() < min_Y_size )
    {
        stop_reason = INITIAL_FAILED;
        stop = true;
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "Stop NM because not enough points in Y." << std::endl;
        return;
    }
    
    // Add points in simplex to obtain dim = n+1 and simplex affinely independant
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itT = T.begin();
    
#ifdef DEBUG
    {
        std::ostringstream oss;
        oss <<"Proceed to simplex creation";
        _out << open_block( oss.str() ) ;
    }
#endif
    
    // This strategy is extremely costly when evaluation number grows
    // Improvements required!
    if ( _p.get_NM_search_init_Y_iter() )
    {
        int k=0;
        
#ifdef DEBUG
        _out << "Use iterative method to select non-dominated points" << endl ;
#endif
        
        // Repeatedly
        //  - Parse T in search of non-dominated points --> S
        //  - Non-dominated points are removed from T and added in Y
        while ( T.size() != 0 && _nm_Y.size() < min_Y_size && !stop  )
        {
            // The set of non-dominated points at each loop
            std::vector<const NOMAD::Eval_Point *> S;
            
            for ( itT = T.begin(); itT !=T.end(); ++itT )
            {
                const NOMAD::Eval_Point * y = (*itT).get_element();
                
                // Test if one of the point in T dominates the candidate for insertion
                std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itT2 = T.begin() ;
                int flag =0;
                
                while( itT2 != T.end() )
                {
                    if ( itT2 == itT )
                    {
                        ++itT2;
                        continue;
                    }
                    
                    const NOMAD::Eval_Point * x = (*itT2).get_element();
                    
                    // Case x dominates y --> y cannot be included in Y
                    if ( NOMAD::NelderMead_Simplex_Eval_Point::dominates(*x,*y) )
                    {
                        flag = 1;
                        break;
                    }
                    
                    ++itT2;
                }
                
                // No point dominates y --> y is included in Y if rank is improved
                if ( flag == 0 )
                    S.push_back( y );
            }
            
            // Remove non-dominated points (in S) from T and add them in Y if rank is increased.
            std::vector<const NOMAD::Eval_Point *>::iterator itS;
            for ( itS=S.begin() ; itS != S.end() ; itS++ )
            {
#ifdef DEBUG
                _out << "k=" << k << ": Point zk:" << endl << *(*itS) ;
#endif
                // Add point in Y. Remove it if rank of Y not increased
                std::pair<std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator,bool> ret = _nm_Y.insert ( NOMAD::NelderMead_Simplex_Eval_Point( *itS) );
                
                if ( ! ret.second )
                {
                    stop_reason = INITIAL_FAILED;
                    stop = true;
                    if ( _display_degree == NOMAD::FULL_DISPLAY )
                        _out << "Stop NM because cannot insert a point in Y." << std::endl;
                    break;
                }
                
                int rank = 0;
                if ( _nm_Y.size() > 1 )
                {
                    rank = get_rank_DZ();
                    
                    if ( rank <= 0 )
                    {
                        stop_reason = INITIAL_FAILED;
                        stop = true;
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                            _out << "Cannot get rank of DZ=[(y1-y0 (y2-y0) ...(yk-y0)]." << std::endl;
                        break;
                        
                    }
                }
                
                
                // Erase last point or not
                if ( rank != k )
                {
#ifdef DEBUG
                    _out << " ---> zk NOT KEPT in Y " << endl;
#endif
                    _nm_Y.erase( NOMAD::NelderMead_Simplex_Eval_Point(*itS) );
                }
                else
                {
#ifdef DEBUG
                    _out << " ---> zk KEPT in Y " << endl;
#endif
                    k++;
                    
                }
                
                // Remove point from T
                T.erase( NOMAD::NelderMead_Simplex_Eval_Point( *itS) );
                
                // Stop if size is ok
                if ( _nm_Y.size() == min_Y_size )
                    break;
                
            }
            S.clear();
        }
        
#ifdef DEBUG
        _out << close_block() << endl;
#endif
        
        // not enough points or insufficient rank of simplex (this is not counted as an error):
        // Erase point
        if ( _nm_Y.size() < min_Y_size || get_rank_DZ() != _n_free )
        {
            stop_reason = INITIAL_FAILED;
            stop = true;
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Stop NM because not enough points in Y or rank of Y < n." << std::endl;
            return;
        }
        
        
        // Update simplex characteristics (volumes and diameter)
        update_Y_characteristics();
        
        make_list_Y0(stop,stop_reason);
        if ( stop )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Stop NM because cannot make Y0." << std::endl;
            return;
        }
        
        make_list_Yn(stop,stop_reason);
        if ( stop )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Stop NM because cannot make Yn." << std::endl;
            return;
        }
        
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss2;
            oss2 <<"After reduction";
            _out << open_block( oss2.str() ) ;
            display_Y_info();
            _out << close_block() <<endl;
        }
        return;
        
    }
    else
    {
        // First point is always added
        _nm_Y.insert ( *itT );
        
#ifdef DEBUG
        _out << "k=0: Point z0:" << endl << (*(*itT).get_point()) ;
        _out << " ---> z0 KEPT in Y " << endl;
#endif
        
        int count_feasible = 0;
        
        if ( _p.has_constraints()
            && (*(*itT).get_point()).get_h().is_defined()
            && (*(*itT).get_point()).is_feasible( _p.get_h_min() ) )
            count_feasible = 1 ;
        
        itT++;
        int k=1;
        while ( itT != T.end() && ! stop )
        {
            
            if ( ! _p.get_NM_search_init_Y_best_von() && _nm_Y.size() == min_Y_size )
                break;
            
#ifdef DEBUG
            _out << "k=" << k << ": Point zk:" << endl << (*(*itT).get_point()) ;
#endif
                        
            std::pair<std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator,bool> ret = _nm_Y.insert ( *itT );
            
            if ( ! ret.second )
            {
                stop_reason = INITIAL_FAILED;
                stop = true;
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                    _out << "Stop NM because cannot insert a point in Y." << std::endl;
                break;
            }
            
            int rank = get_rank_DZ();
            
            if ( rank <= 0 )
            {
                stop_reason = INITIAL_FAILED;
                stop = true;
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                    _out << "Cannot get rank of DZ=[(y1-y0 (y2-y0) ...(yk-y0)]." << std::endl;
                break;
                
            }
            
            // Erase last point or not
            if ( rank != k )
            {
                bool keep_zk_point = false;
                // Try to select the init Y having the best normalized volume
                // Method: keep last zk that should be erased and try to erase a substitute in Y instead that would give a smaller diam(Y) and maintain rank
                if ( _p.get_NM_search_init_Y_best_von() && rank > 2 )
                {
                    
#ifdef DEBUG
                    _out <<"Last zk should be removed because rank=n or rank is not increased."<< endl;
                    _out <<"Try to keep last zk anyway. Instead try to remove most distant points in Y to minimize the diameter." << endl;
#endif
                    
                    // Remove last zk to obtain a reference for the simplex diameter
                    _nm_Y.erase( *itT );
                    
                    // Update simplex diameter
                    update_Y_diameter();
                    
#ifdef DEBUG
                    if ( get_rank_DZ() == _n_free )
                    {
                        update_Y_characteristics();
                        display_Y_info();
                    }
                    _out <<"Before: diam(Y without zk)="<< _simplex_diam <<endl;
#endif
                    
                    const NOMAD::Eval_Point * original_simplex_diam_pt1;
                    const NOMAD::Eval_Point * original_simplex_diam_pt2;
                    
                    // First point on which the diameter is calculated
                    if ( _simplex_diam_pt1 && _simplex_diam_pt1->get_element()->is_defined() )
                    {
                        original_simplex_diam_pt1 = _simplex_diam_pt1->get_element() ;
                    }
                    else
                        continue;
                    
                    // Second point on which the diameter is calculated
                    if ( _simplex_diam_pt2 && _simplex_diam_pt2->get_element()->is_defined() )
                    {
                        
                        original_simplex_diam_pt2 = _simplex_diam_pt2->get_element() ;
                    }
                    else
                        continue;
                    
                    // The simplex diameter without zk
                    double simplex_diam_step0 = _simplex_diam;
                    
                    
                    // Put back zk
                    _nm_Y.insert ( NOMAD::NelderMead_Simplex_Eval_Point( *itT) );
                    
                    
                    bool delete_diam_pt1 = false;
                    bool delete_diam_pt2 = false;
                    
                    double simplex_diam_step1 = _simplex_diam;
                    double simplex_diam_step2 = _simplex_diam;
                    
                    
                    // Try to remove the first point
                    _nm_Y.erase ( NOMAD::NelderMead_Simplex_Eval_Point ( original_simplex_diam_pt1 )  );
                    
#ifdef DEBUG
                    _out << "Try to remove point " << *original_simplex_diam_pt1;
#endif
                    
                    
                    int rank = get_rank_DZ();
                    
                    if ( rank <= 0 )
                    {
                        stop_reason = INITIAL_FAILED;
                        stop = true;
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                            _out << "Minimize diameter 1: cannot get rank of DZ=[(y1-y0 (y2-y0) ...(yk-y0)]." << std::endl;
                        break;
                        
                    }
                    
                    update_Y_diameter();
                    
                    simplex_diam_step1 = _simplex_diam;
                    
#ifdef DEBUG
                    _out << "diam(Y')= " << simplex_diam_step1 <<endl;
#endif
                    
                    // Need to substract 1 to k because _nm_Y contains k points
                    if ( rank == k-1 && simplex_diam_step1 < simplex_diam_step0 )
                        delete_diam_pt1 = true;
                    
                    // Put the pt1 back into Y --> delete will be done later
                    _nm_Y.insert( NOMAD::NelderMead_Simplex_Eval_Point ( original_simplex_diam_pt1 ) );
                    
                    
                    // Try to remove the second point
                    _nm_Y.erase ( NOMAD::NelderMead_Simplex_Eval_Point ( original_simplex_diam_pt2 ) );
                    
#ifdef DEBUG
                    _out << "Try to remove point " << *original_simplex_diam_pt2;
#endif
                    
                    rank = get_rank_DZ();
                    
                    if ( rank <= 0 )
                    {
                        stop_reason = INITIAL_FAILED;
                        stop = true;
                        if ( _display_degree == NOMAD::FULL_DISPLAY )
                            _out << "Minimize diameter 2: cannot get rank of DZ=[(y1-y0 (y2-y0) ...(yk-y0)]." << std::endl;
                        break;
                        
                    }
                    
                    update_Y_diameter();
                    
                    simplex_diam_step2 = _simplex_diam;
                    
#ifdef DEBUG
                    _out << "diam(Y')= " << simplex_diam_step2 <<endl;
#endif
                    
                    // Need to substract 1 to k because _nm_Y contains k points
                    if ( rank == k-1 && simplex_diam_step2 < simplex_diam_step0 )
                        delete_diam_pt2 = true;
                    
                    // Put the pt2 back into Y --> delete will be done later
                    _nm_Y.insert( NOMAD::NelderMead_Simplex_Eval_Point ( original_simplex_diam_pt2 ) );
                    
                    
                    if ( delete_diam_pt1 && delete_diam_pt2 )
                    {
                        if ( simplex_diam_step2 < simplex_diam_step1 )
                        {
                            delete_diam_pt2 = true;
                            delete_diam_pt1 = false;
                        }
                        else
                        {
                            delete_diam_pt2 = false;
                            delete_diam_pt1 = true;
                            
                        }
                    }
                    
                    // Perform delete
                    if ( delete_diam_pt1 )
                    {
                        keep_zk_point = true;
#ifdef DEBUG
                        _out << "Point " << *original_simplex_diam_pt1;
                        _out << " ---> REMOVED from Y to decrease diameter" << endl;
#endif
                        _nm_Y.erase( NOMAD::NelderMead_Simplex_Eval_Point ( original_simplex_diam_pt1 ) );
                        
                        
                    }
                    if ( delete_diam_pt2 )
                    {
                        keep_zk_point = true;
                        
#ifdef DEBUG
                        _out << "Point " << *original_simplex_diam_pt2;
                        _out << " ---> REMOVED from Y to decrease diameter" << endl;
#endif
                        
                        _nm_Y.erase( NOMAD::NelderMead_Simplex_Eval_Point ( original_simplex_diam_pt2 ) );
                    }
                    
                }
                
                if ( ! keep_zk_point )
                {
#ifdef DEBUG
                    _out << " ---> zk NOT KEPT in Y " << endl;
#endif
                    _nm_Y.erase( *itT );
                }
                else
                {
#ifdef DEBUG
                    _out << "Point " << *(*itT).get_element();
                    _out << " ---> zk KEPT in Y " << endl;
#endif
                    // Increase feasible counter if a feasible point is kept (case with constraints)
                    if ( _p.has_constraints()
                        && (*(*itT).get_point()).get_h().is_defined()
                        && (*(*itT).get_point()).is_feasible( _p.get_h_min() ) )
                        count_feasible++;
                }
                itT++;
            }
            else
            {
#ifdef DEBUG
                _out << " ---> zk KEPT in Y " << endl;
#endif
                if ( (*(*itT).get_point()).is_feasible( _p.get_h_min() ) )
                    count_feasible++;
                k++;
                itT++;
            }
            
        }
        
#ifdef DEBUG
        _out << close_block() << endl;
#endif
        
        // not enough points or insufficient rank of simplex (this is not counted as an error):
        // Erase point
        if ( _nm_Y.size() < min_Y_size || get_rank_DZ() != _n_free )
        {
            stop_reason = INITIAL_FAILED;
            stop = true;
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Stop NM because not enough points in Y or rank of Y < n_free." << std::endl;
            return;
        }
        
        
        // Update simplex characteristics (volumes and diameter)
        update_Y_characteristics();
        
        make_list_Y0(stop,stop_reason);
        if ( stop )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Stop NM because cannot make Y0." << std::endl;
            return;
        }
        
        make_list_Yn(stop,stop_reason);
        if ( stop )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Stop NM because cannot make Yn." << std::endl;
            return;
        }
        
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss2;
            oss2 <<"After reduction";
            _out << open_block( oss2.str() ) ;
            display_Y_info();
            _out << close_block() <<endl;
        }
    }
    
}


/*---------------------------------------------------------------*/
/*        project to mesh and create a trial point (private)     */
/*---------------------------------------------------------------*/
bool NOMAD::NelderMead_Search::create_trial_point ( NOMAD::Evaluator_Control &  ev_control,
                                                   NOMAD::Eval_Point         *& x    ,
                                                   const NOMAD::Eval_Point   &  center         )
{
    
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << "candidate";
        if ( _proj_to_mesh )
            _out << " (before projection)";
        _out << ": " << std::endl << *x << std::endl;
    }
    
    // projection to mesh:
    if ( _proj_to_mesh )
    {
        x->project_to_mesh ( center , _p.get_signature()->get_mesh()->get_delta() , _p.get_lb() , _p.get_ub() );
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "candidate (after projection): " << std::endl << *x << std::endl;
    }
    
    
    // Round for integer and binary variables:
    bool has_integer=false;
    bool has_binary=false;
    for (int i=0 ; i < _n ; i++ )
    {
        if ( _p.get_bb_input_type()[i] == NOMAD::INTEGER )
        {
            has_integer=true;
            if ( (*x)[i] >= 0.0 )
                (*x)[i] = (*x)[i].NOMAD::Double::ceil();
            else
                (*x)[i] = (*x)[i].NOMAD::Double::floor();
        }
        // binary variables:
        else if ( _p.get_bb_input_type()[i] == NOMAD::BINARY )
        {
            has_binary=true;
            if ( (*x)[i] < 0.5 )
                (*x)[i] = 0.0;
            else
                (*x)[i] = 1.0;
            
        }
    }
    if ( has_integer && _display_degree == NOMAD::FULL_DISPLAY )
        _out << "candidate (after rounding integer) : ( "
        << *x << " )" << std::endl;
    
    if ( has_binary && _display_degree == NOMAD::FULL_DISPLAY )
        _out << "candidate (after rounding binary) : ( "
        << *x << " )" << std::endl;
    
#ifdef DEBUG
    _out << "Comparison of " << *x << " with " << _nm_Y.size() << " points in Y :" <<  endl;
    _out << "Points in cache: " << endl;
    ev_control.get_cache().display( _out );
#endif
    
    // compare x and all points in simplex
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it ;
    for ( it =_nm_Y.begin(); it != _nm_Y.end(); ++it )
    {
        if ( *x == *(*it).get_point() )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "candidate rejected (candidate already in simplex)" << std::endl;
            return false;
        }
    }

    // compare x and all points in cache
    // Warning: Do not use cache.find(x) --> problem with COOPMADS (point is tagged for evaluation
    // in watched_pts and the find will go in an infinite loop)
    const NOMAD::Cache & cache = (x->get_eval_type() == NOMAD::TRUTH ) ? ev_control.get_cache() : ev_control.get_sgte_cache();
    const NOMAD::Eval_Point * cur = cache.begin() ;
    while ( cur )
    {
        if ( *cur == *x )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "candidate rejected (candidate in cache)" << std::endl;
            return false;
        }
        cur=cache.next();
    }
    
    // WARNING: The direction is not set because it will impact the mesh anisotropy. The anisotropy is better set after polling.
    x->set ( static_cast<int>(_n) , _p.get_bb_nb_outputs() );
    x->set_signature  ( _p.get_signature()  );
    
    // add the new point for evaluation:
    ev_control.add_eval_point ( x                      ,
                               _display_degree         ,
                               _p.get_snap_to_bounds() ,
                               NOMAD::Double()         ,
                               NOMAD::Double()         ,
                               NOMAD::Double()         ,
                               NOMAD::Double()          );
    
    
    return true;
    
    
}

bool NOMAD::NelderMead_Search::NM_step (const NOMAD::Cache                   & cache         ,
                                        const NOMAD::Evaluator               & ev            ,
                                        const NOMAD::Eval_Point              * xk            ,
                                        bool                                 & stop          ,
                                        NM_stop_type                         & stop_reason      )
{
    
    switch ( _step )
    {
        case INITIAL:
            create_initial_sets_from_cache(cache, ev , xk,stop,stop_reason);
            break;
        case REFLECT:
            create_reflect_point(stop,stop_reason,1.0); // without expansion
            break;
        case EXPAND:
            create_reflect_point(stop,stop_reason,_delta_e);  // with expansion
            break;
        case OUTSIDE_CONTRACTION:
            create_reflect_point(stop,stop_reason,_delta_oc); // with outside contraction
            break;
        case INSIDE_CONTRACTION:
            create_reflect_point(stop,stop_reason,_delta_ic); // with inside contraction
            break;
        case SHRINK:
            if ( _perform_shrink )
                create_trial_shrink_points(stop,stop_reason);
            else
            {
                stop = true;
                stop_reason = COMPLETED ;
                return false;
            }
            
            break;
        default:
            stop = true;
            stop_reason = UNDEFINED_STEP;
            return false;
            break;
    }
    
    
    return true;
}

void NOMAD::NelderMead_Search::create_reflect_point (bool                & stop        ,
                                                     NM_stop_type        & stop_reason ,
                                                     const NOMAD::Double & delta )
{
    if ( delta <= -1.0 )
        throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                "NelderMead_Search::create_reflect_point(): delta cannot be <=-1" );
    
    // Clear the NM new points list
    _nm_submitted_points.clear();
    
    size_t min_Y_size = _n_free + 1;
    if ( _nm_Y.size() < min_Y_size )
    {
        stop = true;
        stop_reason = TOO_SMALL_SIMPLEX;
        return;
    }
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "Nelder Mead " ;
        if ( delta ==1 )
            oss << "reflection";
        else if ( delta > 1 )
            oss << "expansion";
        else if ( delta < 1 && delta >= 0 )
            oss << "outside contraction";
        else if ( delta < 0 )
            oss << "inside contraction";
        oss<< " ( delta=" << delta << " ) with " << _nm_Y.size() << " points: " ;
        _out << std::endl << open_block( oss.str() ) ;
    }
    
    // Determine the centroid of Y
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it;
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itBeforeEnd = _nm_Y.end();
    itBeforeEnd--;
    int i=0;
    NOMAD::Point yc(_n,0);
    
    for ( it = _nm_Y.begin() ; it != itBeforeEnd ; it++,i++)
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            _out << "y" << i << ":" << *(*it).get_point() << endl;
        }
        yc = yc + *(*it).get_point();
    }
    yc*=1.0/_n_free; // Average over free variables
    const NOMAD::Point & yn = *(*it).get_point();
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << "yn:" << *(*it).get_point() << endl;
        _out << "yc: x = (" << yc <<")"<< endl << endl;
    }
    
    NOMAD::Point d = yc - yn;
    d *= delta;
    
    // Creation of the reflected point
    NOMAD::Eval_Point * xr = new NOMAD::Eval_Point();
    xr->NOMAD::Point::operator= (yc + d);
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
       _out << "x(NM): " << *xr << endl;
    }
    
    if ( xr->NOMAD::Point::operator== ( yn ) )
    {
        delete xr;
        stop = true;
        if ( delta ==_delta_e )
            stop_reason = EXPANSION_FAILED;
        else if ( delta == 1 )
            stop_reason = REFLECT_FAILED;
        else if ( delta == _delta_ic )
            stop_reason = INSIDE_CONTRACTION_FAILED;
        else if ( delta == _delta_oc )
            stop_reason = OUTSIDE_CONTRACTION_FAILED;
        else
            stop_reason = UNDEFINED_STEP;
        
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            _out << "Reflected point too close to yn." << std::endl;
        }
        
    }
    else
        _nm_submitted_points.push_back( xr );
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << close_block() << endl;
    }
}

void NOMAD::NelderMead_Search::create_trial_shrink_points ( bool         & stop     ,
                                                           NM_stop_type & stop_reason  )
{
    if ( _gamma <= 0.0 || _gamma > 1)
        throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                "NelderMead_Search::create_trial_shrink_points(): valid values 1 >= gamma > 0" );
    
    size_t min_Y_size = _n_free + 1;
    if ( _nm_Y.size() < min_Y_size )
    {
        stop = true;
        stop_reason = TOO_SMALL_SIMPLEX;
        return;
    }
    
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "Nelder Mead shrink ( gamma=" << _gamma << " ) with " << _nm_Y.size() << " points: " ;
        _out << open_block( oss.str() ) ;
    }
    
    // Clear the NM new points list
    _nm_submitted_points.clear();
    
    // Shrink simplex Y
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it = _nm_Y.begin();
    const NOMAD::Eval_Point & y0 = *(*it).get_point();
    int i=0;
    for ( it = _nm_Y.begin() ; it !=_nm_Y.end(); ++it, ++i )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "y" << i << ":" << *(*it).get_point() << endl;
        
        const NOMAD::Point & yi = *(*it).get_point();
        
        // Creation of the shrink point
        
        NOMAD::Point y(y0);
        NOMAD::Point d(yi-y0);
        d*=_gamma;
        y = y + d;
        
        if ( i > 0 && y == yi )
        {
            stop = true;
            stop_reason = SHRINK_FAILED;
            
            if ( _display_degree == NOMAD::FULL_DISPLAY )
            {
                _out << "Shrink point too close to yi." << std::endl;
            }
            
            // Delete all the submitted points so far
            std::list<NOMAD::Eval_Point *>::iterator itS = _nm_submitted_points.begin();
            while ( itS != _nm_submitted_points.end() )
                delete *itS;
            
            _nm_submitted_points.clear();
            return;
        }
        NOMAD::Eval_Point * xs = new NOMAD::Eval_Point();
        xs->NOMAD::Point::operator= (y);
        _nm_submitted_points.push_back( xs );
        
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "--> xs" << i << ":" << *xs << endl;
        
    }
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
        _out << close_block() << endl;
    
    _nm_Y.clear();
}


bool NOMAD::NelderMead_Search::point_dominates_Y0( const NOMAD::Eval_Point &xt )
{
    if ( _nm_Y0.size()==0 )
        throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                "NelderMead_Search::point_dominates_Y0(): _nm_Y0 is empty" );
    
    // The first point of Y is used for comparison
    if ( _p.get_NM_search_use_only_Y() )
    {
        
        NOMAD::NelderMead_Simplex_Eval_Point Y_xt ( &xt );
        
        std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itY = _nm_Y.begin();
        
        if ( Y_xt.is_better_than( (*itY) ) )
            return true;
        else
            return false;
    }
    
    
    
    std::list<const NOMAD::Eval_Point *>::iterator itY0=_nm_Y0.begin();
    
    while ( itY0 != _nm_Y0.end() )
    {
        
        // xt < y ?
        if ( NOMAD::NelderMead_Simplex_Eval_Point::dominates( xt , *(*itY0) ) )
            return true;
        
        itY0++;
        
    }
    return false;
}

bool NOMAD::NelderMead_Search::point_dominates_pts_in_Y( const NOMAD::Eval_Point &xt , size_t nb_points_to_dominate )
{
    
    size_t n_dominates = 0;
    
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itY = _nm_Y.begin();
    
    while ( itY != _nm_Y.end() && n_dominates < nb_points_to_dominate)
    {
        const NOMAD::Eval_Point & y = *(*itY).get_point();
        
        // xt < y ?
        if ( NOMAD::NelderMead_Simplex_Eval_Point::dominates( xt , y ) )
            n_dominates++;
        
        itY++;
        
    }
    return ( n_dominates == nb_points_to_dominate );
}



bool NOMAD::NelderMead_Search::Yn_dominates_point( const NOMAD::Eval_Point &xt )
{
    
    if ( _nm_Yn.size()==0 )
        throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                "NelderMead_Search::Yn_dominates_point(): _nm_Yn is empty" );
    
    // The worst point of Y is used for comparison
    if ( _p.get_NM_search_use_only_Y() )
    {
        NOMAD::NelderMead_Simplex_Eval_Point Y_xt ( &xt );
        
        std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itY = _nm_Y.end();
        itY--;
        
        if ( (*itY).is_better_than(Y_xt) )
            return true;
        else
            return false;
    }
    
    // Without constraints, Yn contains a single point
    std::list<const NOMAD::Eval_Point *>::iterator itYn=_nm_Yn.begin();
    
    int flag = 0;
    
    while ( itYn != _nm_Yn.end() )
    {
        // y < xt ?
        if ( NOMAD::NelderMead_Simplex_Eval_Point::dominates( *(*itYn) , xt ) )
        {
            // One point of Yn dominates x
            flag = 1;
            break;
        }
        
        itYn++;
    }
    
    if ( flag == 1 )
        return true;
    
    // no point of Yn dominates xt --> check if h(yn) < h(xt) --> Yn dominates
    if (  _p.has_constraints() )
    {
        itYn=_nm_Yn.end();
        itYn--;  // Consider yn
        
        // Case with EB constraints and a point from Yn is infeasible
        if ( ! (*itYn)->get_h().is_defined() )
            return false;
        
        // Test also case where xt has no value for h (case with EB constraint)
        if ( ! xt.get_h().is_defined() || (*itYn)->get_h() < xt.get_h() )
            return true;
    }
    return false;
    
}

bool NOMAD::NelderMead_Search::insert_in_Y_best( const NOMAD::Eval_Point &x1 , const NOMAD::Eval_Point &x2 )
{
    
    NOMAD::NelderMead_Simplex_Eval_Point Y_x1 ( &x1 );
    NOMAD::NelderMead_Simplex_Eval_Point Y_x2 ( &x2 );
    
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itYn = _nm_Y.end();
    itYn--;
    
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        std::ostringstream oss;
        oss << "Insertion/Deletion of points in Y: " ;
        _out << std::endl << open_block( oss.str() ) ;
    }
    
    
    std::pair<std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator,bool> ret_x1 = _nm_Y.insert ( Y_x1 );
    
    // Insertion of x1 ok
    if ( ret_x1.second )
    {
        
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss2;
            oss2 << "Insertion in Y: " ;
            _out << open_block( oss2.str() ) ;
            _out << *(*(ret_x1.first)).get_point();
            _out << close_block();
        }
        
        std::pair<std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator,bool> ret_x2;
        
        bool same_x1_x2=false;
        
        // Try to insert x2 only if x1!=x2
        if ( x1 != x2 )
            ret_x2 = _nm_Y.insert ( Y_x2 );
        else
        {
            ret_x2.second = true;
            ret_x2.first = ret_x1.first;
            same_x1_x2 = true;
        }
        
        // Insertion of x2 ok
        if ( ret_x2.second )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
            {
                std::ostringstream oss2;
                oss2 << "Insertion in Y: " ;
                _out << open_block( oss2.str() ) ;
                _out << *(*(ret_x2.first)).get_point();
                _out << close_block();
            }
            
            
            // if x1 and x2 are after yn --> remove both x1 and x2 ==> insertion failed
            if ( std::distance(_nm_Y.begin(), ret_x1.first ) > std::distance(_nm_Y.begin(), itYn ) && std::distance(_nm_Y.begin(), ret_x2.first ) > std::distance(_nm_Y.begin(), itYn ) )
            {
                _nm_Y.erase( ret_x1.first );
                if ( ! same_x1_x2 )
                    _nm_Y.erase( ret_x2.first );
                
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                {
                    _out << "Cannot perform insertion because both points are dominated by yn" << endl;
                    _out << close_block();
                }
                return false;
            }
            
            
            // if x1 after x2 --> keep x2 remove x1
            if ( std::distance(_nm_Y.begin(), ret_x1.first ) > std::distance(_nm_Y.begin(), ret_x2.first ) )
            {
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                {
                    std::ostringstream oss2;
                    oss2 << "Delete from Y: " ;
                    _out << open_block( oss2.str() ) ;
                    _out << *(*(ret_x1.first)).get_point();
                    _out << close_block();
                }
                _nm_Y.erase( ret_x1.first );
            }
            else // if x2 after x1 --> keep x1 remove x2
            {
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                {
                    std::ostringstream oss2;
                    oss2 << "Insertion in Y: " ;
                    _out << open_block( oss2.str() ) ;
                    _out << *(*(ret_x1.first)).get_point();
                    _out << close_block();
                }
                if ( !same_x1_x2 )
                {
                    if ( _display_degree == NOMAD::FULL_DISPLAY )
                    {
                        std::ostringstream oss2;
                        oss2 << "Delete from Y: " ;
                        _out << open_block( oss2.str() ) ;
                        _out << *(*(ret_x2.first)).get_point();
                        _out << close_block();
                    }
                    _nm_Y.erase( ret_x2.first );
                }
            }
            
            // Remove yn
            std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it = _nm_Y.end();
            --it;
            if ( _display_degree == NOMAD::FULL_DISPLAY )
            {
                std::ostringstream oss2;
                oss2 << "Delete yn from Y: " ;
                _out << open_block( oss2.str() ) ;
                _out << *(*(it)).get_point();
                _out << close_block();
            }
            _nm_Y.erase ( it );
            
            // Update simplex characteristics (volumes and diameter)
            update_Y_characteristics();
            
            // Update lists Y0 and Yn
            bool stop;
            NM_stop_type stop_type;
            make_list_Y0(stop,stop_type);
            if ( stop )
            {
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                {
                    _out << "Cannot make Y0." << std::endl;
                    _out << close_block() << std::endl;
                }
                return false;
            }
            make_list_Yn(stop,stop_type);
            if ( stop )
            {
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                {
                    _out << "Cannot make list Yn." << std::endl;
                    _out << close_block() << std::endl;
                }
                return false;
            }
            
            
            if ( _display_degree == NOMAD::FULL_DISPLAY )
            {
                std::ostringstream oss2;
                oss2 <<"After insertion best";
                _out << open_block( oss2.str() ) ;
                display_Y_info();
                _out << close_block() <<endl;
            }
            
            if ( get_rank_DZ() != static_cast<int>(_nm_Y.size()) -1 )
            {
                if ( _display_degree == NOMAD::FULL_DISPLAY )
                {
                    _out << "Rank of DZ=[(y1-y0) (y2-y0) ... (yn-y0)] != n. Y is not a valid simplex. " << std::endl;
                    _out << close_block() << std::endl;
                }
                return false;
            }
            
        }
        else
        {
            _nm_Y.erase(ret_x2.first);
            
            // Update simplex characteristics (volumes and diameter)
            update_Y_characteristics();
            
            if ( _display_degree == NOMAD::FULL_DISPLAY )
            {
                std::ostringstream oss2;
                oss2 <<"After insertion cancelled";
                _out << open_block( oss2.str() ) ;
                display_Y_info();
                _out << close_block() <<endl;
            }
            
            return false;
            
        }
    }
    else
    {
        // Remove point from Y
        _nm_Y.erase(ret_x1.first);
        
        // Update simplex characteristics (volumes and diameter)
        update_Y_characteristics();
        
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss2;
            oss2 <<"After insertion cancelled";
            _out << open_block( oss2.str() ) ;
            display_Y_info();
            _out << close_block() <<endl;
        }
        
        return false;
    }
    if ( _display_degree == NOMAD::FULL_DISPLAY )
    {
        _out << close_block() <<endl;
    }
    
    return true;
}

bool NOMAD::NelderMead_Search::insert_in_Y( const NOMAD::Eval_Point &x )
{
    
    
    NOMAD::NelderMead_Simplex_Eval_Point Y_x ( &x );
    std::pair<std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator,bool> ret_x = _nm_Y.insert ( Y_x );
    
    // Insertion of xe ok
    if ( ret_x.second )
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "Insertion in NM simplex: " ;
            _out << open_block( oss.str() ) ;
            _out << *(*(ret_x.first)).get_point();
            _out << close_block();
        }
        
        // Remove yn
        std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it = _nm_Y.end();
        --it;
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss;
            oss << "Delete from NM simplex: " ;
            _out << open_block( oss.str() ) ;
            _out << *(*(it)).get_point();
            _out << close_block();
        }
        
        // The simplex is unchanged ==> insertion failed
        if ( it == ret_x.first )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Inserted point is last ==> insertion not successful, simplex unchanged." <<endl;
            
            _nm_Y.erase ( it );
            return false;
            
        }
        
        _nm_Y.erase ( it );
        
        // Update the simplex characteristics (volumes and diameter)
        update_Y_characteristics();
        
        
        // Update lists Y0 and Yn
        bool stop;
        NM_stop_type stop_type;
        make_list_Y0(stop,stop_type);
        if ( stop )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Cannot create list Y0 from Y" <<endl;
            return false;
        }
        make_list_Yn(stop,stop_type);
        if ( stop )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Cannot create list Yn from Y" <<endl;
            return false;
        }
        
        
        if ( _display_degree == NOMAD::FULL_DISPLAY )
        {
            std::ostringstream oss2;
            oss2 <<"After insertion";
            _out << open_block( oss2.str() ) ;
            display_Y_info();
            _out << close_block() <<endl;
        }
        
        if ( get_rank_DZ() != static_cast<int>(_nm_Y.size()) -1 )
        {
            if ( _display_degree == NOMAD::FULL_DISPLAY )
                _out << "Rank of DZ=[(y1-y0) (y2-y0) ... (yn-y0)] != n. Y is not a valid simplex. " << endl;
            return false;
        }
        
        
    }
    else
    {
        if ( _display_degree == NOMAD::FULL_DISPLAY )
            _out << "Cannot insert point in Y." <<endl;
        return false;
    }
    return true;
}

/*----------------------------------------------------------------*/
/*  Make list of undominated points from Nelder-Mead simplex Y0   */
/*----------------------------------------------------------------*/
void NOMAD::NelderMead_Search::make_list_Y0 (bool          & stop         ,
                                             NM_stop_type  & stop_reason  )
{
    
    stop = false;
    
    _nm_Y0.clear();
    
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itx;
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator ity=_nm_Y.begin();
    
    size_t max_size_Y0 = _nm_Y.size();
    if ( _p.get_NM_search_use_short_Y0() )
        max_size_Y0 = 2;
    
    _nm_Y0.push_back( (*ity).get_point() );
    ity++;
    while( ity != _nm_Y.end() && _nm_Y0.size() < max_size_Y0 )
    {
        int flag =0;
        const NOMAD::Eval_Point * y = (*ity).get_element();
        itx = _nm_Y.begin();
        
        while ( itx != _nm_Y.end() )
        {
            const NOMAD::Eval_Point * x = (*itx).get_element();
            
            // Case x dominates y --> y cannot be included in Y0
            if ( NOMAD::NelderMead_Simplex_Eval_Point::dominates(*x,*y) )
            {
                flag = 1;
                break;
            }
            
            ++itx;
            
        }
        
        // No point dominates y --> y is included in Y0
        if ( flag == 0 )
            _nm_Y0.push_back( y );
        
        ity++;
    }
}

/*--------------------------------------------------------------*/
/*  Make list of dominated points from Nelder-Mead simplex Yn   */
/*--------------------------------------------------------------*/
void NOMAD::NelderMead_Search::make_list_Yn (bool          & stop         ,
                                             NM_stop_type  & stop_reason  )
{
    
    stop = false;
    _nm_Yn.clear();
    
    
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itx;
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator ity=_nm_Y.begin();
    
    while( ity != _nm_Y.end() )
    {
        int flag =0;
        const NOMAD::Eval_Point * y = (*ity).get_element();
        itx=_nm_Y.begin();
        while ( itx != _nm_Y.end() )
        {
            
            const NOMAD::Eval_Point * x = (*itx).get_element();
            
            // Case y dominates x --> y cannot be included in Yn
            if ( NOMAD::NelderMead_Simplex_Eval_Point::dominates( *y, *x ) )
            {
                flag = 1;
                break;
            }
            
            ++itx;
            
        }
        
        // y dominates no points --> y is included in Y0
        if ( flag == 0 )
            _nm_Yn.push_back( y );
        
        ity++;
    }
    
}


/*---------------------------------------------------------*/
/*  check evaluation point outputs before the integration  */
/*  into NM set (private)                                  */
/*---------------------------------------------------------*/
bool NOMAD::NelderMead_Search::check_outputs ( const NOMAD::Point & bbo , int m ) const
{
    
    if ( bbo.size() != m )
        return false;
    
    for ( int i = 0 ; i < m ; ++i )
        if ( !bbo[i].is_defined() )
            return false;
    
    return true;
}


/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
void NOMAD::NelderMead_Search::update_Y_characteristics()
{
    
    // Update Y diameter
    // -----------------
    update_Y_diameter();
    
    
    // Update Y volumes
    // ----------------
    _simplex_von = -1;
    _simplex_vol = -1;
    
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it1 = _nm_Y.begin();
    
    if ( _nm_Y.size() != (size_t)_n_free + 1 )
    {
#ifdef DEBUG
        _out << "Cannot get the volume of simplex Y when its dimension is not n_free+1" << endl;
#endif
        return;
    }
    
    const NOMAD::Point * y0 = (*it1).get_point();
    if ( (size_t)y0->size() != (size_t)_n )
    {
#ifdef DEBUG
        _out << "Cannot get the volume of simplex Y when dimension of an element is not n" << endl;
#endif
        return;
    }
    
    NOMAD::Point Delta (static_cast<int>(_n),1.0);
    
    if ( _p.get_NM_search_scaled_DZ() )
    {
        const NOMAD::Signature * s = (*it1).get_point()->get_signature();
        
        // A point may not have a signature (case when points come from a cache file -> we suppose file cache points are initial points )
        if ( s != NULL )
            s->get_mesh()->get_Delta( Delta ) ;
        else
        {
            Delta = _p.get_initial_poll_size();
        }
    }    
    
    // Update volume
    //---------------
    it1 = _nm_Y.begin();
    
    // V : vector of yk-y0 ( for determinant, need square array (n_free x n_free) --> exclude fixed variables )
    double ** V = new double *[_n_free];
    for (size_t i = 0 ; i < (size_t)_n_free ; i++ )
        V[i]=new double [_n_free];
    
    int j = 0;
    it1++;
    while ( it1 != _nm_Y.end() )
    {
        int k = 0;
        for ( size_t i = 0; i < (size_t)_n ; i++ )
        {
            if ( ! _fixed_variables[i].is_defined() )
            {
                if ( k == _n_free )
                    throw NOMAD::Exception ( "NelderMead_Search.cpp" , __LINE__ ,
                                            "NelderMead_Search::update_Y_characteristics(): index > n_free" );
                
                V[j][k] =  ((*(*it1).get_point())[i].value() - (*y0)[i].value()) / Delta[i].value() ;
                k++;
            }
        }
        it1++;
        j++;
        
    }
    
    // Get determinant of DZ
    double det,nfact=1;
    
    bool success = NOMAD::get_determinant(V, det , _n_free);
    
    for ( size_t i = 0; i < (size_t)_n_free ; i++ )
        delete [] V[i];
    delete [] V;
    
    if ( success )
    {
#ifdef DEBUG
        _out << "The determinant of the matrix: det( [(y1-y0) (y2-y0) ... (ynf-y0)] ) = " <<  det << endl;
#endif
        
        // !n
        for ( size_t i=1 ; i < (size_t)_n_free+1 ; i++)
            nfact*=i;
        
        _simplex_vol = fabs( det )/nfact;  // Use fact(n_free) for volume
        
        if ( _simplex_diam > 0 )
            _simplex_von = _simplex_vol / pow(_simplex_diam,_n_free) ;  // Use n_free for the normalized volume
#ifdef DEBUG
        else
            _out << "Cannot get the normalized volume of simplex Y because simplex diameter <=0 " << endl;
#endif
        
        return;
    }
#ifdef DEBUG
    else
        _out << "Cannot get the volume of simplex Y because determinant failed." << endl;
#endif
    
    return ;
}


/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
void NOMAD::NelderMead_Search::update_Y_diameter()
{
    
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it1 = _nm_Y.begin();
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator it2;
    
    NOMAD::Point Delta (static_cast<int>(_n),1.0);
    
    if ( _p.get_NM_search_scaled_DZ() )
    {
        const NOMAD::Signature * s = (*it1).get_point()->get_signature();
        
        // A point may not have a signature (case when points come from a cache file -> we suppose file cache points are initial points )
        if ( s != NULL )
            s->get_mesh()->get_Delta( Delta ) ;
        else
        {
            Delta = _p.get_initial_poll_size();
        }
    }
    
    _simplex_diam = 0;
    for ( it1 = _nm_Y.begin() ; it1 != _nm_Y.end() ; it1++  )
    {
        it2 = it1;
        it2++;
        while ( it2 != _nm_Y.end() )
        {
            const NOMAD::Point DZ = ( *(*it1).get_point() - *(*it2).get_point() ) / Delta ;
            double length_DZ = pow(DZ.squared_norm().value(),0.5);
            
            if ( length_DZ > _simplex_diam )
            {
                _simplex_diam = length_DZ;
                _simplex_diam_pt1 = &(*it1);
                _simplex_diam_pt2 = &(*it2);
            }
            it2++;
        }
    }
    
    return ;
}


/*----------------------------------------------------------------------*/
/* Get the rank of the matrix DZ = [(y1-y0) (y2-y0) ... (yk-y0)]]       */
/*----------------------------------------------------------------------*/
int NOMAD::NelderMead_Search::get_rank_DZ() const
{
    // The dimension of DZ (k) is related to Y
    size_t k = _nm_Y.size() - 1 ;
    
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itY = _nm_Y.begin();
    
    const NOMAD::Point * y0 = (*itY).get_point();
    
    NOMAD::Point Delta (static_cast<int>(_n),1.0);
    
    if ( _p.get_NM_search_scaled_DZ() )
    {
        const NOMAD::Signature * s = (*itY).get_point()->get_signature();
        
        // A point may not have a signature (case when points come from a cache file -> we suppose file cache points are initial points )
        if ( s != NULL )
            s->get_mesh()->get_Delta( Delta ) ;
        else
        {
            Delta = _p.get_initial_poll_size();
        }
        
#ifdef DEBUG
        _out << "DZ is scaled with Delta: "<<Delta <<endl;
#endif
        
    }
    
    // DZ : vector of yk-y0 (multidimensional array)
    double ** DZ = new double *[k];
    for (size_t i = 0 ; i < k ; ++i )
        DZ[i]=new double [_n];
    
#ifdef DEBUG
    _out << "The rank of DZ=[";
#endif
    
    // Create DZ
    itY++;
    size_t j=0;
    while ( j < k )
    {
#ifdef DEBUG
        _out << " (" ;
#endif
        for ( size_t i = 0; i < (size_t)_n ; i++ )
        {
            DZ[j][i] = ((*(*itY).get_point())[i].value()- (*y0)[i].value()  ) / Delta[i].value();
#ifdef DEBUG
            _out << DZ[j][i] << " ";
#endif
        }
        j++;
        itY++;
#ifdef DEBUG
        _out << ")" ;
#endif
        
    }
    
    // Get the rank
    int rank= NOMAD::get_rank(DZ , k , _n , _p.get_NM_search_rank_eps().value() );
    
#ifdef DEBUG
    _out << " ] equals " << rank <<endl;
#endif
    
    for (size_t i=0 ; i < k ; ++i)
        delete [] DZ[i];;
    delete [] DZ;
    
    return rank;
}



/*---------------------------------------------------------*/
/*---------------------------------------------------------*/
void NOMAD::NelderMead_Search::display_Y_info()const
{
    
    _out << "Number of points in the simplex Y: " << _nm_Y.size() <<endl;
    
    if ( _simplex_vol > 0 )
        _out << "The volume of the simplex: " << _simplex_vol <<  endl;
    else
        _out << "WARNING: Evaluation of the simplex volume failed." << endl;
    
    if ( _simplex_diam > 0 )
        _out << "The diameter of the simplex: " << _simplex_diam <<  endl;
    else
        _out << "WARNING: Evaluation of the simplex diameter failed." << endl;
    
    if ( _simplex_von > 0 )
        _out << "The normalized volume of the simplex: " << _simplex_von <<  endl;
    else
        _out << "WARNING: Evaluation of the simplex diameter failed." << endl;
    
    
#ifdef DEBUG
    std::set<NOMAD::NelderMead_Simplex_Eval_Point>::iterator itY = _nm_Y.begin();
    _out << "The simplex Y contains: " << std::endl;
    for ( itY =_nm_Y.begin(); itY != _nm_Y.end(); ++itY)
        _out << *((*itY).get_point()) ;
    _out << endl;
    get_rank_DZ();
    
#endif
    
    if ( _p.has_constraints() )
    {
        _out << "Number of points in Y0: " << _nm_Y0.size() <<endl;
        _out << "Number of points in Yn: " << _nm_Yn.size() <<endl;
        
#ifdef DEBUG
        std::list<const NOMAD::Eval_Point *>::const_iterator itY0 = _nm_Y0.begin();
        _out << "The list Y0 contains: " << std::endl;
        for ( itY0 =_nm_Y0.begin(); itY0 != _nm_Y0.end(); ++itY0)
            _out << *(*itY0) ;
        _out << endl;
        
        std::list<const NOMAD::Eval_Point *>::const_iterator itYn = _nm_Yn.begin();
        _out << "The list Yn contains: " << std::endl;
        for ( itYn =_nm_Yn.begin(); itYn != _nm_Yn.end(); ++itYn)
            _out << *(*itYn) ;
        _out << endl;
        
#endif
    }
    
}
