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
 \file   Stats.hpp
 \brief  Algorithm stats (headers)
 \author Sebastien Le Digabel
 \date   2010-04-22
 \see    Stats.cpp
 */
#ifndef __STATS__
#define __STATS__

#include "Clock.hpp"
#include "Double.hpp"
#include "Model_Stats.hpp"

namespace NOMAD {
    
    /// Algorithm stats.
    class Stats {
        
    private:
        
        /// Number of evaluations.
        /**
         Blackbox evaluations + cache hits.
         */
        int           _eval;
        
        /// Number of simulated blackbox evaluations.
        /**
         Blackbox evaluations + initial cache hits.
         */
        int           _sim_bb_eval;
        
        int           _sgte_eval;        ///< Number of surrogate evaluations.
        int           _sgte_cost;        ///< Surrogate cost.
        int           _bb_eval;          ///< Number of blackbox evaluations.
        int           _block_eval;       ///< Number of block of evaluations.
        int           _failed_eval;      ///< Number of failed evaluations.
        int           _cache_hits;       ///< Number of cache hits.
        int           _interrupted_eval; ///< Number of interrupted sequence of evaluations.
        int           _iterations;       ///< Number of iterations.
        NOMAD::Double _stat_sum;         ///< Sum stat (output of type NOMAD::STAT_SUM).
        NOMAD::Double _stat_avg;         ///< Average stat (output of type NOMAD::STAT_AVG).
        int           _cnt_avg;          ///< Number of values for the \c avg stat.
        int           _p1_iterations;    ///< Iterations in MADS phase one.
        int           _p1_bbe;           ///< Blackbox evaluations in phase one.
        NOMAD::Clock  _clock;            ///< Wall-clock time.
        int           _mads_runs;        ///< Number of MADS runs (for multi-objective runs).
        
        // Polls:
        int           _nb_poll_searches; ///< Number of poll searches.
        int           _poll_pts;         ///< Number of poll points.
        int           _poll_success;     ///< Number of poll successes.
        
        // Extended polls:
        int           _nb_ext_polls;      ///< Number of extended polls.
        int           _ext_poll_pts;      ///< Number of extended poll points.
        int           _ext_poll_succ;     ///< Number of extended poll successes.
        int           _ext_poll_bb_eval;  ///< Number of extended poll blackbox evaluations.
        int           _ext_poll_descents; ///< Number of extended poll descents.
        
        // Speculative searches:
        int           _nb_spec_searches;  ///< Number of speculative searches.
        int           _spec_pts;          ///< Number of speculative search points.
        int           _spec_success;      ///< Number of speculative search successes.
        
#ifdef USE_MPI
        int           _asynchronous_success; ///< Number of asynchronous successes.
        int           MPI_MPI_data_size;        ///< Size of MPI messages.
#endif
        
        // LH searches:
        int           _nb_LH_searches;    ///< Number of Latin-Hypercube (LH) searches.
        int           LH_LH_pts;            ///< Number of Latin-Hypercube (LH) search points.  //zhenghua
        int           LH_LH_success;        ///< Number of LH search successes. //zhenghua
        
        // NM searches:
        int           _nb_NM_searches;    ///< Number of Nelder-Meads (NM) searches.
        int           NM_NM_pts;            ///< Number of NM search points.  //zhenghua
        int           NM_NM_step ;          ///< Number of NM steps.
        int           NM_NM_initial_step ;          ///< Number of NM initial steps.
        int           NM_NM_success;        ///< Number of NM search successes.
        int           NM_NM_bb_eval;        ///< Number of NM search blackbox evaluations.
        int           NM_NM_reflect_step;    ///< Number of NM reflect step.
        int           NM_NM_expand_step;    ///< Number of NM expand step.
        int           NM_NM_inside_contraction_step; ///< Number of NM inside contraction step.
        int           NM_NM_outside_contraction_step;///< Number of NM outside contraction step.
        int           NM_NM_shrink_step;    ///< Number of NM shrink step.
        
        // TM line searches:
        int           _nb_TM_searches;    ///< Number of Trend matrix (TM) line searches.
        int           TM_TM_pts;            ///< Number of TM search points.  //zhenghua
        int           TM_TM_success;        ///< Number of TM line search successes.
        int           TM_TM_bb_eval;        ///< Number of TM search blackbox evaluations.

        
        // Cache searches:
        int           _nb_cache_searches; ///< Number of cache searches (CS).
        int           CS_CS_pts;            ///< Number of CS search points.  //zhenghua
        int           CS_CS_success;        ///< Number of CS search successes.
        
        /// Model stats.
        NOMAD::Model_Stats _model_stats;
        
        // VNS searches:
        int           _nb_VNS_searches;   ///< Number of VNS searches.
        int           VNS_VNS_pts;           ///< Number of VNS search points.
        int           VNS_VNS_success;       ///< Number of VNS search successes.
        int           VNS_VNS_bb_eval;       ///< Number of VNS blackbox evaluations.
        int           VNS_VNS_sgte_eval;     ///< Number of VNS surrogate evaluations.
        
        // User searches:
        int           _nb_usr_searches;   ///< Number of user searches.
        int           _usr_srch_pts;      ///< Number of user search points.
        int           _usr_srch_success;  ///< Number of user search successes.
        
        // Dynamic management of poll directions
        int           _nb_success_dyn_dir;  ///< Number of successfull polling in the direction added dynamically
        
        
    public:
        
        /// Constructor.
        /**
         \param sgte_cost Surrogate cost -- \b IN -- \b optional (default = \c -1).
         */
        explicit Stats ( int sgte_cost = -1 ) : _sgte_cost(sgte_cost) { reset(); }
        
        /// Copy constructor.
        /**
         \param s The copied object -- \b IN.
         */
        explicit Stats ( const Stats & s )
        : _eval                 ( s._eval                 ) ,
        _sim_bb_eval          ( s._sim_bb_eval          ) ,
        _sgte_eval            ( s._sgte_eval            ) ,
        _sgte_cost            ( s._sgte_cost            ) ,
        _bb_eval              ( s._bb_eval              ) ,
        _block_eval           ( s._block_eval        ) ,
        _failed_eval          ( s._failed_eval          ) ,
        _cache_hits           ( s._cache_hits           ) ,
        _interrupted_eval     ( s._interrupted_eval     ) ,
        _iterations           ( s._iterations           ) ,
        _stat_sum             ( s._stat_sum             ) ,
        _stat_avg             ( s._stat_avg             ) ,
        _cnt_avg              ( s._cnt_avg              ) ,
        _p1_iterations        ( s._p1_iterations        ) ,
        _p1_bbe               ( s._p1_bbe               ) ,
        _clock                ( s._clock                ) ,
        _mads_runs            ( s._mads_runs            ) ,
        _nb_poll_searches     ( s._nb_poll_searches     ) ,
        _poll_pts             ( s._poll_pts             ) ,
        _poll_success         ( s._poll_success         ) ,
        _nb_ext_polls         ( s._nb_ext_polls         ) ,
        _ext_poll_pts         ( s._ext_poll_pts         ) ,
        _ext_poll_succ        ( s._ext_poll_succ        ) ,
        _ext_poll_bb_eval     ( s._ext_poll_bb_eval     ) ,
        _ext_poll_descents    ( s._ext_poll_descents    ) ,
        _nb_spec_searches     ( s._nb_spec_searches     ) ,
        _spec_pts             ( s._spec_pts             ) ,
        _spec_success         ( s._spec_success         ) ,
#ifdef USE_MPI
        _asynchronous_success ( s._asynchronous_success ) ,
        MPI_MPI_data_size        ( s.MPI_MPI_data_size        ) ,
#endif
        _nb_LH_searches       ( s._nb_LH_searches       ) ,
        LH_LH_pts               ( s.LH_LH_pts               ) ,
        LH_LH_success           ( s.LH_LH_success           ) ,
        _nb_NM_searches       ( s._nb_NM_searches       ) ,
        NM_NM_pts               ( s.NM_NM_pts               ) ,
        NM_NM_step              ( s.NM_NM_step              ) ,
        NM_NM_initial_step      ( s.NM_NM_initial_step      ) ,
        NM_NM_success           ( s.NM_NM_success           ) ,
        NM_NM_bb_eval           ( s.NM_NM_bb_eval           ) ,
        NM_NM_reflect_step      ( s.NM_NM_reflect_step      ) ,
        NM_NM_expand_step       ( s.NM_NM_expand_step       ) ,
        NM_NM_inside_contraction_step ( s.NM_NM_inside_contraction_step ),
        NM_NM_outside_contraction_step ( s.NM_NM_outside_contraction_step ),
        NM_NM_shrink_step       ( s.NM_NM_shrink_step       ) ,
        _nb_TM_searches       ( s._nb_TM_searches       ) ,
        TM_TM_pts               ( s.TM_TM_pts               ) ,
        TM_TM_success           ( s.TM_TM_success           ) ,
        _nb_cache_searches    ( s._nb_cache_searches    ) ,
        CS_CS_pts               ( s.CS_CS_pts               ) ,
        CS_CS_success           ( s.CS_CS_success           ) ,
        _nb_VNS_searches      ( s._nb_VNS_searches      ) ,
        VNS_VNS_pts              ( s.VNS_VNS_pts              ) ,
        VNS_VNS_success          ( s.VNS_VNS_success          ) ,
        VNS_VNS_bb_eval          ( s.VNS_VNS_bb_eval          ) ,
        VNS_VNS_sgte_eval        ( s.VNS_VNS_sgte_eval        ) ,
        _nb_usr_searches      ( s._nb_usr_searches      ) ,
        _usr_srch_pts         ( s._usr_srch_pts         ) ,
        _usr_srch_success     ( s._usr_srch_success     ) ,
        _nb_success_dyn_dir   ( s._nb_success_dyn_dir   ) {}
        
        /// Affectation operator.
        /**
         \param s The right-hand side object -- \b IN.
         */
        Stats & operator = ( const Stats & s );
        
        /// Destructor.
        virtual ~Stats ( void ) {}
        
        /// Reset the stats.
        void reset ( void );
        
        /// Update stats from another NOMAD::Stats object.
        /**
         \param s The other NOMAD::Stats object -- \b IN.
         \param for_search A flag equal to \c true if the method
         has been invoked within a search step
         -- \b IN.
         */
        void update ( const Stats & s , bool for_search );
        
        /// Update the \c sum stat.
        /**
         \param d New \c sum element -- \b IN.
         */
        void update_stat_sum ( const NOMAD::Double & d );
        
        /// Update the \c avg stat.
        /**
         \param d New \c avg element -- \b IN.
         */
        void update_stat_avg ( const NOMAD::Double & d );
        
        /// Add \c 1 to stat \c _mads_run.
        void add_mads_run ( void ) { ++_mads_runs; }
        
        /// Set the number of MADS runs.
        void set_mads_runs ( int mads_runs ) { _mads_runs = mads_runs; }
        
        /// Add \c 1 to stat \c _eval.
        void add_eval ( void ) { ++_eval; }
        
        /// Add \c 1 to stat \c _sim_bb_eval.
        void add_sim_bb_eval ( void ) { ++_sim_bb_eval; }
        
        /// Add \c 1 to stat \c _sgte_eval.
        void add_sgte_eval ( void ) { ++_sgte_eval; }
        
        /// Add \c count_eval to stat \c _sgte_eval.
        void add_sgte_eval ( int count_eval ) { _sgte_eval+=count_eval; }
        
        /// Add \c 1 to stat \c _bb_eval.
        void add_bb_eval ( void ) { ++_bb_eval; }
        
        /// Add \c 1 to stat \c _block_eval.
        void add_one_block_eval ( void ) { ++_block_eval; }
        
        /// Add \c count_eval to stat \c _bb_eval.
        void add_bb_eval ( int count_eval ) { _bb_eval+=count_eval; }
        
        /// Add \c 1 to stat \c _failed_eval.
        void add_failed_eval ( void ) { ++_failed_eval; }
        
        /// Add \c 1 to stat \c _cache_hits.
        void add_cache_hit ( void ) { ++_cache_hits; }

        /// Add an integer to stat \c _cache_hits.
        /**
         \param i The integer -- \b IN.
         */
        void add_cache_hits ( int i ) { _cache_hits+=i ;}
        
        /// Add \c 1 to stat \c _interrupted_eval.
        void add_interrupted_eval ( void ) { ++_interrupted_eval; }
        
        /// Add \c 1 to stat \c _iterations.
        void add_iteration ( void ) { ++_iterations; }
        
        /// Add \c 1 to stat \c _nb_poll_searches.
        void add_nb_poll_searches ( void ) { ++_nb_poll_searches; }
        
        /// Add \c 1 to stat \c _poll_success.
        void add_poll_success ( void ) { ++_poll_success; }
        
        /// Add \c 1 to stat \c _nb_ext_polls.
        void add_nb_ext_polls ( void ) { ++_nb_ext_polls; }
        
        /// Add \c 1 to stat \c _ext_poll_succ.
        void add_ext_poll_succ ( void ) { ++_ext_poll_succ; }
        
        /// Add \c 1 to stat \c _ext_poll_descents.
        void add_ext_poll_descent ( void ) { ++_ext_poll_descents; }
        
        /// Add \c 1 to stat \c _nb_spec_searches.
        void add_nb_spec_searches ( void ) { ++_nb_spec_searches; }
        
        /// Add \c 1 to stat \c _spec_success.
        void add_spec_success ( void ) { ++_spec_success; }
        
        /// Add \c 1 to stat \c _LH_success.
        void add_LH_success ( void ) { ++LH_LH_success; }
        
        /// Add \c 1 to stat \c _NM_success.
        void add_NM_success ( void ) { ++NM_NM_success; }

        /// Add \c 1 to stat \c _NM_step.
        void add_NM_step ( void ) { ++NM_NM_step; }

        /// Add \c 1 to stat \c _NM_initial_step.
        void add_NM_initial_step ( void ) { ++NM_NM_initial_step; }
        
        /// Add \c 1 to stat \c _NM_expand_step.
        void add_NM_expand_step ( void ) { ++NM_NM_expand_step; }

        /// Add \c 1 to stat \c _NM_expand_step.
        void add_NM_reflect_step ( void ) { ++NM_NM_reflect_step; }
        
        /// Add \c 1 to stat \c _NM_inside_contraction_step.
        void add_NM_inside_contraction_step ( void ) { ++NM_NM_inside_contraction_step; }

        /// Add \c 1 to stat \c _NM_outside_contraction_step.
        void add_NM_outside_contraction_step ( void ) { ++NM_NM_outside_contraction_step; }

        /// Add \c 1 to stat \c _NM_shrink_step.
        void add_NM_shrink_step ( void ) { ++NM_NM_shrink_step; }
        
        /// Add \c 1 to stat \c _TM_success.
        void add_TM_success ( void ) { ++TM_TM_success; }
        
        /// Add \c 1 to stat \c _nb_cache_searches.
        void add_nb_cache_searches ( void ) { ++_nb_cache_searches; }
        
        /// Add \c 1 to stat \c _nb_LH_searches.
        void add_nb_LH_searches ( void ) { ++_nb_LH_searches; }
        
        /// Add \c 1 to stat \c _nb_NM_searches.
        void add_nb_NM_searches ( void ) { ++_nb_NM_searches; }

        /// Add \c 1 to stat \c _nb_TM_searches.
        void add_nb_TM_searches ( void ) { ++_nb_TM_searches; }

        /// Add \c 1 to stat \c _CS_success.
        void add_CS_success ( void ) { ++CS_CS_success; }
        
        /// Add \c 1 to stat \c _nb_VNS_searches.
        void add_nb_VNS_searches ( void ) { ++_nb_VNS_searches; }
        
        /// Add \c 1 to stat \c _VNS_success.
        void add_VNS_success ( void ) { ++VNS_VNS_success; }
        
        /// Add \c 1 to stat \c _nb_usr_searches.
        void add_nb_usr_searches ( void ) { ++_nb_usr_searches; }
        
        /// Add \c 1 to stat \c _usr_srch_success.
        void add_usr_srch_success ( void ) { ++_usr_srch_success; }
        
        /// Add \c 1 to stat \c _nb_success_dyn_dir.
        void add_nb_success_dyn_dir(void) {++_nb_success_dyn_dir;}
        
        /// Add an integer to stat \c _p1_iterations.
        /**
         \param i The integer -- \b IN.
         */
        void add_p1_iterations ( int i ) { _p1_iterations += i; }
        
        /// Add an integer to stat \c _p1_bbe.
        /**
         \param i The integer -- \b IN.
         */
        void add_p1_bbe ( int i ) { _p1_bbe += i; }
        
        /// Add an integer to stat \c _poll_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_poll_pts ( int i ) { _poll_pts += i; }
        
        /// Add an integer to stat \c _ext_poll_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_ext_poll_pts ( int i ) { _ext_poll_pts += i; }
        
        /// Add an integer to stat \c _ext_poll_bb_eval.
        /**
         \param i The integer -- \b IN.
         */
        void add_ext_poll_bb_eval ( int i ) { _ext_poll_bb_eval += i;}
        
        /// Add an integer to stat \c _spec_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_spec_pts ( int i ) { _spec_pts += i; }
        
        /// Add an integer to stat \c _LH_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_LH_pts ( int i ) { LH_LH_pts += i; }
        
        /// Add an integer to stat \c _NM_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_NM_pts ( int i ) { NM_NM_pts += i; }

        /// Add an integer to stat \c _TM_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_TM_pts ( int i ) { TM_TM_pts += i; }

        
        /// Add an integer to stat \c _CS_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_CS_pts ( int i ) { CS_CS_pts += i; }
        
        /// Update model stats.
        void update_model_stats ( const NOMAD::Model_Stats & ms )
        {
            _model_stats.update ( ms );
        }
        
        /// Add an integer to stat \c _VNS_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_VNS_pts ( int i ) { VNS_VNS_pts += i; }
        
        /// Add an integer to stat \c _VNS_bb_eval.
        /**
         \param i The integer -- \b IN.
         */
        void add_VNS_bb_eval ( int i ) { VNS_VNS_bb_eval += i; }

        /// Add an integer to stat \c _NM_bb_eval.
        /**
         \param i The integer -- \b IN.
         */
        void add_NM_bb_eval ( int i ) { NM_NM_bb_eval += i; }

        /// Add an integer to stat \c _TM_bb_eval.
        /**
         \param i The integer -- \b IN.
         */
        void add_TM_bb_eval ( int i ) { TM_TM_bb_eval += i; }

        
        /// Add an integer to stat \c _VNS_sgte_eval.
        /**
         \param i The integer -- \b IN.
         */
        void add_VNS_sgte_eval ( int i ) { VNS_VNS_sgte_eval += i; }
        
        /// Add an integer to stat \c _usr_srch_pts.
        /**
         \param i The integer -- \b IN.
         */
        void add_usr_srch_pts ( int i ) { _usr_srch_pts += i; }
        
#ifdef USE_MPI
        
        /// Add \c 1 to stat \c _asynchronous_success.
        void add_asynchronous_success ( void ) { ++_asynchronous_success; }
        
        /// Add an integer to stat \c _MPI_data_size.
        /**
         \param i The integer -- \b IN.
         */
        void set_MPI_data_size ( int i ) { MPI_MPI_data_size = i; }
#endif
        
        /// Access to the stat \c _eval.
        /**
         \return The stat \c _eval.
         */
        int get_eval ( void ) const { return _eval; }
        
        /// Access to the stat \c _sim_bb_eval.
        /**
         \return The stat \c _sim_bb_eval.
         */
        int get_sim_bb_eval ( void ) const { return _sim_bb_eval; }
        
        /// Access to the stat \c _sgte_eval.
        /**
         \return The stat \c _sgte_eval.
         */
        int get_sgte_eval ( void ) const { return _sgte_eval; }
        
        /// Access to the real time stat.
        /**
         \return The real time stat.
         */
        int get_real_time ( void ) const { return _clock.get_real_time(); }
        
        /// Access to the stat \c _iterations.
        /**
         \return The stat \c _iterations.
         */
        int get_iterations ( void ) const { return _iterations; }
        
        /// Access to the stat \c _failed_eval.
        /**
         \return The stat \c _failed_eval.
         */
        int get_failed_eval ( void ) const { return _failed_eval; }
        
        /// Access to the stat \c _mads_runs.
        /**
         \return The stat \c _mads_runs.
         */
        int get_mads_runs ( void ) const { return _mads_runs; }
        
        /// Access to the stat \c _LH_pts.
        /**
         \return The stat \c _LH_pts.
         */
        int get_LH_pts ( void ) const { return LH_LH_pts; }
        
        /// Access to the stat \c _NM_pts.
        /**
         \return The stat \c _NM_pts.
         */
        int get_NM_pts ( void ) const { return NM_NM_pts; }

        /// Access to the stat \c _NM_success.
        /**
         \return The stat \c _NM_success.
         */
        int get_NM_success ( void ) const { return NM_NM_success; }

        /// Access to the stat \c _NM_expand_step.
        /**
         \return The stat \c _NM_expand_step.
         */
        int get_NM_expand_step ( void ) const { return NM_NM_expand_step; }

        /// Access to the stat \c _NM_initial_step.
        /**
         \return The stat \c _NM_initial_step.
         */
        int get_NM_initial_step ( void ) const { return NM_NM_initial_step; }

        
        /// Access to the stat \c _NM_reflect_step.
        /**
         \return The stat \c _NM_reflect_step.
         */
        int get_NM_reflect_step ( void ) const { return NM_NM_reflect_step; }

        /**
         \return The stat \c _NM_inside_contraction_step.
         */
        int get_NM_inside_contraction_step ( void ) const { return NM_NM_inside_contraction_step; }
        
        /**
         \return The stat \c _NM_outside_contraction_step.
         */
        int get_NM_outside_contraction_step ( void ) const { return NM_NM_outside_contraction_step; }
        
        /**
         \return The stat \c _NM_shrink_step.
         */
        int get_NM_shrink_step ( void ) const { return NM_NM_shrink_step; }
        
        /**
         \return The stat \c _NM_bb_eval.
         */
        int get_NM_bb_eval ( void ) const { return NM_NM_bb_eval; }

        
        /// Access to the stat \c _nb_TM_searches.
        /**
         \return The stat \c _nb_TM_searches.
         */
        int get_nb_TM_searches ( void ) const { return _nb_TM_searches; }
        
        /// Access to the stat \c _TM_pts.
        /**
         \return The stat \c _TM_pts.
         */
        int get_TM_pts ( void ) const { return TM_TM_pts; }
        
        /// Access to the stat \c _TM_success.
        /**
         \return The stat \c _TM_success.
         */
        int get_TM_success ( void ) const { return TM_TM_success; }

        /**
         \return The stat \c _TM_bb_eval.
         */
        int get_TM_bb_eval ( void ) const { return TM_TM_bb_eval; }
        
        /// Access to the stat \c _CS_pts.
        /**
         \return The stat \c _CS_pts.
         */
        int get_CS_pts ( void ) const { return CS_CS_pts; }
        
        /// Access to the stat \c _VNS_bb_eval.
        /**
         \return The stat \c _VNS_bb_eval.
         */
        int get_VNS_bb_eval ( void ) const { return VNS_VNS_bb_eval; }
        
        /// Access to the stat \c _nb_VNS_SEARCHES.
        /**
         \return The stat \c _nb_VNS_SEARCHES.
         */
        int get_nb_VNS_SEARCHES ( void ) const { return _nb_VNS_searches; }

        /// Access to the stat \c _nb_NM_SEARCHES.
        /**
         \return The stat \c _nb_NM_SEARCHES.
         */
        int get_nb_NM_SEARCHES ( void ) const { return _nb_NM_searches; }

        /// Access to the stat \c _VNS_sgte_eval.
        /**
         \return The stat \c _VNS_sgte_eval.
         */
        int get_VNS_sgte_eval ( void ) const { return VNS_VNS_sgte_eval; }
        
        /// Access to the number of cache hits.
        /**
         \return The number of cache hits.
         */
        int get_cache_hits ( void ) const { return _cache_hits; }
        
        /// Access to the number of blackbox evaluations (includes surrogate cost).
        /**
         \return The number of blackbox evaluations.
         */
        int get_bb_eval ( void ) const
        {
            return ( _sgte_cost > 0 ) ? _bb_eval + _sgte_eval / _sgte_cost : _bb_eval;
        }
        
        
        /// Access to the number of block of evaluations (includes bb and surrogates).
        /**
         \return The number of blackbox evaluations.
         */
        int get_block_eval ( void ) const
        {
            return _block_eval;
        }
        
        /// Access to the \c sum stat.
        /**
         \return The \c sum stat.
         */
        NOMAD::Double get_stat_sum ( void ) const
        {
            return _stat_sum;
        }
        
        /// Access to the \c avg stat.
        /**
         \return The \c avg stat.
         */
        NOMAD::Double get_stat_avg ( void ) const
        {
            return ( _cnt_avg > 0 ) ? _stat_avg/_cnt_avg : NOMAD::Double();
        }
        
        /// Access to the stat \c _nb_ext_polls.
        /**
         \return The stat \c _nb_ext_polls.
         */
        int get_nb_ext_poll ( void ) const
        {
            return _nb_ext_polls;
        }
        
        /// Access to the model stats.
        /**
         \return The model stats.
         */
        const NOMAD::Model_Stats & get_model_stats ( void ) const { return _model_stats; }
        
        /// Display.
        /**
         \param out The NOMAD::Display object -- \b IN.
         */
        void display ( const NOMAD::Display & out ) const;
    };
    
    /// Display a NOMAD::Stats object.
    /**
     \param out The NOMAD::Display object               -- \b IN.
     \param s   The NOMAD::Stats object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
     */
    inline const NOMAD::Display & operator << ( const NOMAD::Display & out ,
                                               const NOMAD::Stats   & s     )
    {
        s.display ( out );
        return out;
    }
}

#endif
