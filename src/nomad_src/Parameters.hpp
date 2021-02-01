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
 \file   Parameters.hpp
 \brief  NOMAD Parameters (headers)
 \author Sebastien Le Digabel
 \date   2010-04-21
 \see    Parameters.cpp
 */
#ifndef __PARAMETERS__
#define __PARAMETERS__

#include <algorithm>
#include <fstream>
#include <ctime>
#include <map>
#include "Parameter_Entries.hpp"
#include "Signature.hpp"
#include "Algo_Parameters.hpp"

namespace NOMAD {
    
    /// Structure to represent the models parameters.
    struct model_params_type {
        
        NOMAD::model_type search1;        ///< First model search (MS) type.
        NOMAD::model_type search2;        ///< Second model search (MS) type.
        
        NOMAD::model_type eval_sort;      ///< Use models to sort trial pts (MES).
        
        bool search_optimistic;           ///< If the MS is optimistic or not.
        bool search_proj_to_mesh;         ///< Model solution projected or not to mesh.
        int  search_max_trial_pts;        ///< Limit on the number of MS trial pts.
        bool eval_sort_cautious;          ///< Cautious strategy for MES.
        
        NOMAD::Double quad_radius_factor; ///< Quadratic model radius factor \c r.
        bool          quad_use_WP;        ///< Use or not a well-poisedness strategy.
        int           quad_min_Y_size;    ///< Limit inf on the size of \c Y.
        int           quad_max_Y_size;    ///< Limit sup on the size of \c Y.
        NOMAD::Double model_np1_quad_epsilon;///< Ortho n+1 quadratic model epsilon (scaling used for (n+1)th dynamic direction with Ortho n+1)
        
    };
    
    /// Class for the NOMAD parameters.
    class DLL_API Parameters : public Algo_Parameters, private NOMAD::Uncopyable {
        
    private:
        
        bool           _to_be_checked; ///< Access control to the parameters.
        static  bool   _warning_has_been_displayed;
        NOMAD::Display _out;           ///< Display.
        
        /*----------------------------------------------------------------------*/
        
        /// Initializations.
        void init ( void );
        
        /// Interpretation of the entries for the starting point.
        /**
         \param entries Parameter entries -- \b IN.
         */
        void interpret_x0 ( const NOMAD::Parameter_Entries & entries );
        
        /// Interpretation of the entries for mesh and poll sizes.
        /**
         \param entries Parameter entries -- \b IN.
         \param param_name Parameter name -- \b IN.
         */
        void interpret_mesh_sizes ( const NOMAD::Parameter_Entries & entries    ,
                                   const std::string              & param_name   );
        
        /// Interpretation of the entries for granular variables.
        /**
         \param entries Parameter entries -- \b IN.
         \param param_name Parameter name -- \b IN.
         */
        void interpret_granularity ( const NOMAD::Parameter_Entries & entries    ,
                                    const std::string              & param_name   );
        
        
        /// Interpretation of the entries for bounds and fixed variables.
        /**
         \c BFVS stands for "Bounds and Fixed VariableS".
         \param entries Parameter entries -- \b IN.
         \param param_name Parameter name -- \b IN.
         */
        void interpret_BFVS ( const NOMAD::Parameter_Entries & entries    ,
                             const std::string              & param_name   );
        
        /// Interpretation of the entries for parameter \c BB_INPUT_TYPE.
        /**
         \param entries Parameter entries -- \b IN.
         */
        void interpret_bb_input_type ( const NOMAD::Parameter_Entries & entries );
        
        /// Interpretation of the entries for parameter \c PERIODIC_VARIABLE.
        /**
         \param entries Parameter entries -- \b IN.
         */
        void interpret_periodic_var ( const NOMAD::Parameter_Entries & entries );
        
        /// Interpretation of the entries for parameter \c VARIABLE_GROUP.
        /**
         \param entries Parameter entries -- \b IN.
         */
        void interpret_var_groups ( const NOMAD::Parameter_Entries & entries );
        
        // Interpretation of the entries for parameter \c F_TARGET.
        /**
         \param entries Parameter entries -- \b IN.
         */
        void interpret_f_target ( const NOMAD::Parameter_Entries & entries );
        
        /// Interpretation of the entries for parameter \c TREND_MATRIX.
        /**
         \param entries Parameter entries -- \b IN.
         */
        void interpret_trend_matrix ( const NOMAD::Parameter_Entries & entries );
        
        /// Delete the list of starting points.
        void delete_x0s ( void );
        
        /// Check a directory name.
        /**
         \param dir Directory name -- \b IN.
         \return A boolean equal to \c true if the directory is valid.
         */
        static bool check_directory ( std::string & dir );
        
        /// Add a seed to a file name.
        /**
         - Transforms \c file_name.ext into \c file_name.seed.ext.
         \param n_seed Length of the seed \c s_seed     -- \b IN.
         \param s_seed A string representating the seed -- \b IN.
         \param file_name File name                     -- \b IN/OUT.
         */
        static void add_seed_to_file_name ( int                 n_seed    ,
                                           const std::string & s_seed    ,
                                           std::string       & file_name   );
        
        /// Check if a specified direction type is set.
        /**
         \param dt The specified direction type -- \b IN.
         \return true if the specified direction type is used.
         */
        bool has_direction_type(NOMAD::direction_type dt) const;
        
        
        /*----------------------------------------------------------------------*/
        
    public:
        
        /// Exception class for an invalid parameter.
        class Invalid_Parameter : public NOMAD::Exception {
        public:
            /// Constructor.
            Invalid_Parameter ( const std::string & file ,
                               int                 line ,
                               const std::string & msg    )
            : NOMAD::Exception ( file , line , msg ) {
                _type_msg = "Invalid Parameter";
            }
            
        };
        
        /// Exception class for a bad access.
        class Bad_Access : public NOMAD::Exception {
        public:
            /// Constructor.
            Bad_Access ( const std::string & file ,
                        int                 line ,
                        const std::string & msg    )
            : NOMAD::Exception ( file , line , msg ) {}
        };
        
        /*----------------------------------------------------------------------*/
        
        /// Constructor #1.
        /**
         \param out Display -- \b IN.
         */
        explicit Parameters ( const NOMAD::Display & out )
        : _out              ( out  ) ,
        _std_signature    ( NULL ) ,
        _extern_signature ( NULL )   { init(); }
        
        /// Constructor #2.
        /**
         From an extern signature.
         \param extern_signature A pointer to the extern signature -- \b IN.
         \param out Display -- \b IN.
         */
        explicit Parameters ( NOMAD::Signature     * extern_signature ,
                             const NOMAD::Display & out                )
        : _out              ( out              ) ,
        _std_signature    ( NULL             ) ,
        _extern_signature ( extern_signature )   { init(); }
        
        /// Destructor.
        virtual ~Parameters ( void );
        
        // Added for Algo_Parameters
        /// Check if algo parameters are compatible with current parameters
        /**
         \param p Parameters for compatibility  -- \b IN.
         */
        bool is_algo_compatible ( const NOMAD::Parameters & p) const ;
        
        
        // Added for Algo_Parameters
        /// Check if algo parameters are compatible with current parameters
        /**
         \param ap Algorithmic parameters for compatibility  -- \b IN.
         */
        virtual bool is_algo_compatible ( const Algo_Parameters & ap) const ;
        
        
        // Added for Algo_Parameters
        /// Return algo name
        virtual std::string get_algo_name() const { return std::string("NOMAD");}
        
        // Added for Algo_Parameters
        /// Return algo version
        virtual std::string get_algo_version() const { return BASE_VERSION;}
        
        /// Display parameter help.
        /**
         - Help option obtained with command \c nomad \c -h.
         - A list of parameters is given as the program arguments.
         \param argc      Number of command line arguments              -- \b IN.
         \param argv      The command line arguments                    -- \b IN.
         \param developer Bool to request developer help (defaut=false) -- \b IN.
         */
        void help ( int argc , char ** argv , bool developer=false) const;
        
        /// Display parameter help.
        /**
         For a specific parameter.
         \param param_name Name of the parameter                         -- \b IN.
         \param developer  Bool to request developer help (defaut=false) -- \b IN.
         */
        void help ( const std::string & param_name, bool developer=false ) const;
        
        /// Display parameter help.
        /**
         For a list of parameters.
         \param param_names_list    List of parameter names              -- \b IN.
         \param developer           Bool to request developer help (defaut=false) -- \b IN.
         */
        void help ( const std::list<std::string> & param_names_list , bool developer=false) const;
        
        /// Reset.
        void reset ( void ) { init(); }
        
        /// Read parameters from a file. 1/2
        /**
         \param param_file Name of the parameters file -- \b IN.
         */
        void read ( const std::string & param_file );
        
        /// Read parameters from a Parameter_Entries 2/2
        /**
         \param entries -- \b IN.
         */
        void read ( const NOMAD::Parameter_Entries & entries );

	/// Read parameters from iostream. This is for snomadr,  zhenghua
	/**
	  \param fin Name of the io stream -- \b IN.
					*/
	void read (std::iostream & fin);
	
        /// Check the parameters.
        /**
         \param remove_history_file A boolean equal to \c true if the history file
         has to be deleted
         -- \b IN -- \b optional (default = \c true).
         \param remove_solution_file A boolean equal to \c true if the solution file
         has to be deleted
         -- \b IN -- \b optional (default = \c true).
         \param remove_stats_file A boolean equal to \c true if the stats file
         has to be deleted
         -- \b IN -- \b optional (default = \c true).
         */
        void check ( bool remove_history_file  = true ,
                    bool remove_solution_file = true ,
                    bool remove_stats_file    = true   );
        
        /// Force the check flag to be satisfied.
        /**
         This advanced function should be used with care.
         */
        void force_check_flag ( void ) { _to_be_checked = false; }
        
        /// Display.
        /**
         \param out The NOMAD::Display object -- \b IN.
         */
        void display ( const NOMAD::Display & out ) const;
        
        /// Display.
        /**
         Uses \c this->_out as NOMAD::Display object.
         */
        void display ( void ) const { display ( _out ); }
        
        /// Reset the display for warning.
        /**
         */
        static void reset_display_warning ( ){ _warning_has_been_displayed = false ;}
        
        
        /*--------------------------------------------------------------*/
        /*           Attributes and methods listed by categories        */
        /*--------------------------------------------------------------*/
        
        // Algorithm and miscellaneous parameters:
        // ---------------------------------------
    private:
        std::string            _problem_dir;   ///< Problem directory.
        std::string            _tmp_dir;       ///< Temporary directory.
        bool                   _reject_unknown_parameters;  ///< Error if an unknown parameter is encountered
        int                    _seed;          ///< Seed.
        int                    _max_eval;      ///< Maximum number of evaluations.
        int                    _max_bb_eval;   ///< Maximum number of blackbox evaluations.
        
        // Report blackbox evaluation number for sub-algos (VNS,ExtendedPoll) on display stat
        // (if false the number of eval is not reported and a symbol "+" is displayed
        bool                   _report_bbe_value; ///< Parameter \c REPORT_BBE_VALUE
        
        // Report surrogate evaluation number for sub-algos (VNS,ExtendedPoll) on display stat
        // (if false the number of eval is not reported and a symbol "+" is displayed
        bool                   _report_sgte_value; ///< Parameter \c REPORT_SGTE_VALUE
        
        
        // Report block evaluation number for sub-algos (VNS, ExtendedPoll) on display stat
        // (if false the number of eval is not reported and a symbol "+" is displayed
        bool                   _report_blk_eva_value; ///< Parameter \c REPORT_BLK_EVA_VALUE
        
        int                    _max_block_eval;  ///< Maximum number of block evaluations.
        
        std::list<std::string> _display_stats;    ///< Stats keywords for \c DISPLAY_STATS.
        bool                   _display_all_eval; ///< Parameter \c DISPLAY_ALL_EVAL
        
        /// Equal to \c true if \c _max_bb_eval has been entered.
        bool                   _max_bbe_decided;
        
        /// Maximum number of simulated evaluations.
        int                    _max_sim_bb_eval;
        
        
        /// Maximum number of evaluation during intensification of poll and/or search
        int                    _max_eval_intensification;
        
        /// Intensification on poll or on search or on both
        NOMAD::intensification_type _intensification_type;
        
        int                    _max_time;           ///< Maximum time.
        int                    _max_iterations;     ///< Maximum number of iterations.
        int                    _max_cons_failed_it; ///< Maximum number of consecutive failed iterations.
        float                  _max_cache_memory;   ///< Maximum cache memory.
        bool                   _stop_if_feasible;   ///< Stop if a feasible solution is found.
        NOMAD::Point           _f_target;           ///< Target for the objective function.
        NOMAD::Double          _stat_sum_target;    ///< Target for the STAT_SUM stat.
        NOMAD::Double          L_L_curve_target;     ///< Target for the L_CURVE criterion.  // zhenghua, for compiling on Solaris.
        bool                   _snap_to_bounds;     ///< Snap or not the points to the bounds.
        bool                   _user_calls_enabled; ///< Enable calls to user functions.
        bool                   _asynchronous;       ///< Asynchronous version for parallelism.
        
    public:
        
        /// Test to know if parameters have to be checked.
        /**
         \return A boolean equal to \c true if parameters have
         to be checked with function NOMAD::Parameters::check().
         */
        bool to_be_checked ( void ) const { return _to_be_checked; }
        
        /// When an unknown parameter is read, should an error occur.
        /**
         \return true if unknown parameters trigger an error (default) and false if unknown parameters are ignored.
         */
        bool get_reject_unknown_parameters() { return _reject_unknown_parameters; }
        
        /// Access to the display degrees.
        /**
         - Use NOMAD::Parameters::out().get_X_dd()
         to access other specific display degrees.
         - \see Display.hpp.
         \param d The 4 display degrees as a string -- \b OUT.
         */
        void get_display_degree ( std::string & d ) const;
        
        /// Access to the display degree.
        /**
         \return General display degree.
         */
        int get_display_degree ( void ) const;
        
        /// Access to display.
        /**
         \return The NOMAD::Display object.
         */
        const NOMAD::Display & out ( void ) const;
        
        /// Access to the \c DISPLAY_STATS parameter.
        /**
         \return The \c DISPLAY_STATS parameter.
         */
        const std::list<std::string> & get_display_stats ( void ) const;
        
        /// Access to the \c DISPLAY_ALL_EVAL parameter.
        /**
         \return The \c DISPLAY_ALL_EVAL parameter.
         */
        bool get_display_all_eval ( void ) const;
        
        /// Access to the \c POINT_DISPLAY_LIMIT parameter.
        /**
         \return The \c POINT_DISPLAY_LIMIT parameter.
         */
        int get_point_display_limit ( void ) const;
        
        /// Access to the seed.
        /**
         \return The seed.
         */
        int get_seed ( void ) const;
        
        /// Access to the maximum number of evaluations.
        /**
         \return The maximum number of evaluations.
         */
        int get_max_eval ( void ) const;
        
        /// Access to the \c MAX_SIM_BB_EVAL parameter.
        /**
         \return The \c MAX_SIM_BB_EVAL parameter.
         */
        int get_max_sim_bb_eval ( void ) const;
        
        /// Access to the \c MAX_BB_EVAL parameter.
        /**
         \return The \c MAX_BB_EVAL parameter.
         */
        int get_max_bb_eval ( void ) const;
        
        /// Access to the \c MAX_BLOCk_EVAL parameter.
        /**
         \return The \c MAX_BLOCK_EVAL parameter.
         */
        int get_max_block_eval ( void ) const;
        
        /// Access to the \c INTENSITIFICATION_TYPE parameter.
        /**
         \return The \c INTENSIFICATION_TYPE parameter.
         */
        NOMAD::intensification_type get_intensification_type ( void ) const { return _intensification_type ; }
        
        /// Access to the \c MAX_EVAL_INTENSIFICATION parameter.
        /**
         \return The \c MAX_EVAL_INTENSIFICATION parameter.
         */
        int get_max_eval_intensification ( void ) const { return _max_eval_intensification ; }
        
        
        
        
        /// Access to the \c MAX_TIME parameter.
        /**
         \return The \c MAX_TIME parameter.
         */
        int get_max_time ( void ) const;
        
        /// Access to the maximum number of iterations.
        /**
         \return The maximum number of iterations.
         */
        int get_max_iterations ( void ) const;
        
        /// Access to the maximum number of consecutive failed iterations.
        /**
         \return The maximum number of consecutive failed iterations.
         */
        int get_max_consecutive_failed_iterations ( void ) const;
        
        /// Access to the maximum cache memory.
        /**
         \return The maximum cache memory.
         */
        float get_max_cache_memory ( void ) const;
        
        /// Access to the \c STOP_IF_FEASIBLE parameter.
        /**
         \return The \c STOP_IF_FEASIBLE parameter.
         */
        bool get_stop_if_feasible ( void ) const;
        
        /// Access to the \c F_TARGET parameter.
        /**
         \return The \c F_TARGET parameter.
         */
        const NOMAD::Point & get_f_target ( void ) const;
        
        /// Access to the \c STAT_SUM_TARGET parameter.
        /**
         \return The \c STAT_SUM_TARGET parameter.
         */
        const NOMAD::Double & get_stat_sum_target ( void ) const;
        
        /// Access to the \c L_CURVE_TARGET parameter.
        /**
         \return The \c L_CURVE_TARGET parameter.
         */
        const NOMAD::Double & get_L_curve_target ( void ) const;
        
        /// Access to the problem directory.
        /**
         \return The problem directory.
         */
        const std::string & get_problem_dir ( void ) const;
        
        /// Access to the temporary directory.
        /**
         \return The temporary directory.
         */
        const std::string & get_tmp_dir ( void ) const;
        
        /// Access to the \c SNAP_TO_BOUNDS parameter.
        /**
         \return The \c SNAP_TO_BOUNDS parameter.
         */
        bool get_snap_to_bounds ( void ) const;
        
        /// Access to the \c USER_CALLS_ENABLED parameter.
        /**
         \return The \c USER_CALLS_ENABLED parameter.
         */
        bool get_user_calls_enabled ( void ) const;
        
        /// Access to the \c ASYNCHRONOUS parameter.
        /**
         \return The \c ASYNCHRONOUS parameter.
         */
        bool get_asynchronous ( void ) const;
        
        
        /// Access to the \c REPORT_BBE_VALUE parameter.
        /**
         \return The \c REPORT_BBE_VALUE parameter.
         */
        bool get_report_bbe_value ( void ) const;
        
        /// Access to the \c REPORT_SGTE_VALUE parameter.
        /**
         \return The \c REPORT_SGTE_VALUE parameter.
         */
        bool get_report_sgte_value ( void ) const;
        
        
        /// Access to the \c REPORT_BLK_EVA_VALUE parameter.
        /**
         \return The \c REPORT_BLK_EVA_VALUE parameter.
         */
        bool get_report_blk_eva_value ( void ) const;
        
        
        /// Set the \c REJECT_UNKNOWN_PARAMETERS parameter.
        /**
         \param reject The \c REJECT_UNKNOWN_PARAMETERS value.
         */
        void set_REJECT_UNKNOWN_PARAMETERS(const bool reject);
        
        /// Set the \c POINT_DISPLAY_LIMIT parameter.
        /**
         \param dl The \c POINT_DISPLAY_LIMIT parameter -- \b IN.
         */
        void set_POINT_DISPLAY_LIMIT ( int dl );
        
        /// Set the \c DISPLAY_STATS parameter.
        /**
         - From a list of strings.
         \param ds The \c DISPLAY_STATS parameter -- \b IN.
         */
        void set_DISPLAY_STATS ( const std::list<std::string> & ds );
        
        /// Set the \c DISPLAY_STATS parameter.
        /**
         - From a string.
         \param ds The \c DISPLAY_STATS parameter -- \b IN.
         */
        void set_DISPLAY_STATS ( const std::string & ds );
        
        /// Set the \c DISPLAY_ALL_EVAL parameter.
        /**
         \param dae The \c DISPLAY_ALL_EVAL parameter -- \b IN.
         */
        void set_DISPLAY_ALL_EVAL ( bool dae );
        
        /// Set the display degree.
        /**
         - Accepts also NOMAD::dd_type arguments.
         \param dd Display degree -- \b IN.
         \return \c true if the operation succeeded.
         */
        bool set_DISPLAY_DEGREE ( int dd );
        
        /// Set the display degrees.
        /**
         - From a string with the 4 degrees.
         \param dd Display degree -- \b IN.
         \return \c true if the operation succeeded.
         */
        bool set_DISPLAY_DEGREE ( const std::string & dd );
        
        /// Set the display degrees.
        /*
         \param gen_dd    General display degree   -- \b IN.
         \param search_dd Search display degree    -- \b IN.
         \param poll_dd   Poll display degree      -- \b IN.
         \param iter_dd   Iterative display degree -- \b IN.
         */
        void set_DISPLAY_DEGREE ( int gen_dd    ,
                                 int search_dd ,
                                 int poll_dd   ,
                                 int iter_dd     );
        
        /// Set the \c OPEN_BRACE parameter.
        /**
         \param s The \c OPEN_BRACE parameter -- \b IN.
         */
        void set_OPEN_BRACE ( const std::string & s );
        
        /// Set the \c CLOSED_BRACE parameter.
        /**
         \param s The \c CLOSED_BRACE parameter -- \b IN.
         */
        void set_CLOSED_BRACE ( const std::string & s );
        
        /// Set the seed.
        /**
         \param seed The seed -- \b IN.
         */
        void set_SEED ( int seed );
        
        /// Set the \c MAX_EVAL parameter.
        /**
         \param me The \c MAX_EVAL parameter -- \b IN.
         */
        void set_MAX_EVAL ( int me );
        
        /// Set the \c MAX_SIM_BB_EVAL parameter.
        /**
         \param msbe The \c MAX_SIM_BB_EVAL parameter -- \b IN.
         */
        void set_MAX_SIM_BB_EVAL ( int msbe );
        
        /// Set the \c MAX_BB_EVAL parameter.
        /**
         \param bbe The \c MAX_BB_EVAL parameter -- \b IN.
         */
        void set_MAX_BB_EVAL ( int bbe );
        
        /// Set the \c MAX_BLOCK_EVAL parameter.
        /**
         \param blk The \c MAX_BLOCK__EVAL parameter -- \b IN.
         */
        void set_MAX_BLOCK_EVAL ( int blk );
        
        /// Set the \c MAX_EVAL_INTENSIFICATION parameter.
        /**
         \param bbe The \c MAX_EVAL_INTENSIFICATION parameter -- \b IN.
         */
        void set_MAX_EVAL_INTENSIFICATION ( int bbe );
        
        /// Set the \c INTENSIFICATION_TYPE parameter.
        /**
         \param it The \c INTENSIFICATION_TYPE parameter -- \b IN.
         */
        void set_INTENSIFICATION_TYPE ( NOMAD::intensification_type it );
        
        /// Set the \c POLL_INTENSIFICATION_DIRECTION parameter.
        /**
         \param pid The \c POLL_INTENSIFICATION_DIRECTION parameter -- \b IN.
         */
        void set_POLL_INTENSIFICATION_DIRECTION ( NOMAD::direction_type pid );
        
        /// Set the \c MAX_TIME parameter.
        /**
         \param mt The \c MAX_TIME parameter -- \b IN.
         */
        void set_MAX_TIME ( int mt );
        
        /// Set the \c MAX_ITERATIONS parameter.
        /**
         \param mi The \c MAX_ITERATIONS parameter -- \b IN.
         */
        void set_MAX_ITERATIONS ( int mi );
        
        /// Set the \c MAX_CONSECUTIVE_FAILED_ITERATIONS parameter.
        /**
         \param mcsi The \c MAX_CONSECUTIVE_FAILED_ITERATIONS parameter -- \b IN.
         */
        void set_MAX_CONSECUTIVE_FAILED_ITERATIONS ( int mcsi );
        
        /// Set the \c MAX_CACHE_MEMORY parameter.
        /**
         \param mcm The \c MAX_CACHE_MEMORY parameter -- \b IN.
         */
        void set_MAX_CACHE_MEMORY ( float mcm );
        
        /// Set the \c STOP_IF_FEASIBLE parameter.
        /**
         \param sif The \c STOP_IF_FEASIBLE parameter -- \b IN.
         */
        void set_STOP_IF_FEASIBLE ( bool sif );
        
        /// Set the \c STAT_SUM_TARGET parameter.
        /**
         \param target The \c STAT_SUM_TARGET parameter -- \b IN.
         */
        void set_STAT_SUM_TARGET ( const NOMAD::Double & target );
        
        /// Set the \c L_CURVE_TARGET parameter.
        /**
         \param target The \c L_CURVE_TARGET parameter -- \b IN.
         */
        void set_L_CURVE_TARGET ( const NOMAD::Double & target );
        
        /// Set the \c PROBLEM_DIR parameter.
        /**
         \param dir The \c PROBLEM_DIR parameter -- \b IN.
         */
        void set_PROBLEM_DIR ( const std::string & dir );
        
        /// Set the \c TMP_DIR parameter.
        /**
         \param dir The \c TMP_DIR parameter -- \b IN.
         */
        void set_TMP_DIR ( const std::string & dir );
        
        /// Set the \c SNAP_TO_BOUNDS parameter.
        /**
         \param stb The \c SNAP_TO_BOUNDS parameter -- \b IN.
         */
        void set_SNAP_TO_BOUNDS ( bool stb );
        
        /// Set the \c USER_CALLS_ENABLED parameter.
        /**
         \param uce The \c USER_CALLS_ENABLED parameter -- \b IN.
         */
        void set_USER_CALLS_ENABLED ( bool uce );
        
        /// Set the \c ASYNCHRONOUS parameter.
        /**
         \param as The \c ASYNCHRONOUS parameter -- \b IN.
         */
        void set_ASYNCHRONOUS ( bool as );
        
        /// Set the \c F_TARGET parameter.
        /**
         - Multi-objective version.
         \param target The \c F_TARGET parameter -- \b IN.
         */
        void set_F_TARGET ( const NOMAD::Point & target );
        
        /// Set the F_TARGET \c parameter.
        /**
         - Single-objective version.
         \param target The \c F_TARGET parameter -- \b IN.
         */
        void set_F_TARGET ( const NOMAD::Double & target );
        
        /// Reset the \c F_TARGET parameter.
        void reset_f_target ( void ) { _f_target.clear(); }
        
        /// Set the \c EPSILON parameter.
        /**
         \param eps The \c EPSILON parameter -- \b IN.
         */
        void set_EPSILON ( const NOMAD::Double & eps )
        {
            NOMAD::Double::set_epsilon ( eps.value() );
        }
        
        /// Set the \c UNDEF_STR parameter.
        /**
         \param undef_str The \c UNDEF_STR parameter -- \b IN.
         */
        void set_UNDEF_STR ( const std::string & undef_str )
        {
            NOMAD::Double::set_undef_str ( undef_str );
        }
        
        /// Set the \c INF_STR parameter.
        /**
         \param inf_str The \c INF_STR parameter -- \b IN.
         */
        void set_INF_STR ( const std::string & inf_str )
        {
            NOMAD::Double::set_inf_str ( inf_str );
        }
        
        /// Access to the \c EPSILON parameter.
        /**
         \return The \c EPSILON parameter.
         */
        const NOMAD::Double get_epsilon ( void ) const
        {
            return NOMAD::Double::get_epsilon();
        }
        
        /// Access to the \c UNDEF_STR parameter.
        /**
         \return The \c UNDEF_STR parameter.
         */
        const std::string get_undef_str ( void ) const
        {
            return NOMAD::Double::get_undef_str();
        }
        
        /// Access to the \c INF_STR parameter.
        /**
         \return The \c INF_STR parameter.
         */
        const std::string get_inf_str ( void ) const
        {
            return NOMAD::Double::get_inf_str();
        }
        
        // Output files:
        // -------------
    private:
        
        /// List of stats for parameter \c STATS_FILE.
        std::list<std::string> _stats_file;
        
        /// Name of the stats file.
        std::string _stats_file_name;
        
        /// Parameter \c ADD_SEED_TO_FILE_NAME.
        bool _add_seed_to_file_names;
        
        std::string _solution_file;     ///< Parameter \c SOLUTION_FILE.
        std::string _history_file;      ///< Parameter \c HISTORY_FILE.
        std::string _cache_file;        ///< Parameter \c CACHE_FILE.
        int         _cache_save_period; ///< Parameter \c CACHE_SAVE_PERIOD.
        
    public:
        
        /// Access to the list of stats for the \c STATS_FILE parameter.
        /**
         \return The list of stats.
         */
        const std::list<std::string> & get_stats_file ( void ) const;
        
        /// Access to the name of the stats file.
        /**
         \return The name of the stats file.
         */
        const std::string & get_stats_file_name ( void ) const;
        
        /// Access to the \c SOLUTION_FILE parameter.
        /**
         \return The \c SOLUTION_FILE parameter.
         */
        const std::string & get_solution_file ( void ) const;
        
        /// Access to the \c HISTORY_FILE parameter.
        /**
         \return The \c HISTORY_FILE parameter.
         */
        const std::string & get_history_file ( void ) const;
        
        /// Access to the \c CACHE_FILE parameter.
        /**
         \return The \c CACHE_FILE parameter.
         */
        const std::string & get_cache_file ( void ) const;
        
        /// Access to the \c CACHE_SAVE_PERIOD parameter.
        /**
         \return The \c CACHE_SAVE_PERIOD parameter.
         */
        int get_cache_save_period ( void ) const;
        
        /// Access to the \c ADD_SEED_TO_FILE_NAME parameter.
        /**
         \return The \c ADD_SEED_TO_FILE_NAME parameter.
         */
        bool get_add_seed_to_file_names ( void ) const;
        
        /// Reset the \c STATS_FILE parameter.
        void reset_stats_file ( void );
        
        /// Set the \c STATS_FILE parameter.
        /**
         \param file_name Name of the stats file -- \b IN.
         \param stats     List of stats          -- \b IN.
         */
        void set_STATS_FILE ( const std::string            & file_name ,
                             const std::list<std::string> & stats       );
        
        /// Set the \c STATS_FILE parameter.
        /**
         \param file_name Name of the stats file -- \b IN.
         \param stats     The stats              -- \b IN.
         */
        void set_STATS_FILE ( const std::string & file_name ,
                             const std::string & stats       );
        
        /// Set the \c ADD_SEED_TO_FILE_NAME parameter.
        /**
         \param astfn The \c ADD_SEED_TO_FILE_NAME parameter -- \b IN.
         */
        void set_ADD_SEED_TO_FILE_NAMES ( bool astfn );
        
        /// Set the \c SOLUTION_FILE parameter.
        /**
         \param sf The \c SOLUTION_FILE parameter -- \b IN.
         */
        void set_SOLUTION_FILE ( const std::string & sf );
        
        /// Set the \c HISTORY_FILE parameter.
        /**
         \param hf The \c HISTORY_FILE parameter -- \b IN.
         */
        void set_HISTORY_FILE ( const std::string & hf );
        
        /// Set the \c CACHE_FILE parameter.
        /**
         \param cf The \c CACHE_FILE parameter -- \b IN.
         */
        void set_CACHE_FILE ( const std::string & cf );
        
        /// Set the \c CACHE_SAVE_PERIOD parameter.
        /**
         \param csp The \c CACHE_SAVE_PERIOD parameter -- \b IN.
         */
        void set_CACHE_SAVE_PERIOD ( int csp );
        
        // Searches:
        // ---------
    private:
        
        /// Speculative search (MADS search).
        bool _speculative_search;
        
        bool _disable_models;  ///< Models disablement
        
        NOMAD::model_params_type _model_params; ///< Models parameters.
        
        // VNS search parameters:
        bool          VNS_VNS_search;  ///< Flag for the VNS search.   //zhenghua
        NOMAD::Double VNS_VNS_trigger; ///< VNS trigger.
        
        
        // NelderMead search simplex parameters:
        bool          NM_NM_search;  ///< Flag for the NM search.   //zhenghua
        
        // NelderMead simplex update parameters
        Double        NM_NM_gamma ;     ///< Shrink parameter
        Double        NM_NM_delta_ic ;  ///< Inside contraction parameter
        Double        NM_NM_delta_oc ;  ///< Outside contraction parameter
        Double        NM_NM_delta_e ;   ///< Expansion parameter
        
        // NelderMead search parameters:
        bool          NM_NM_search_intensive;  ///< Flag for the NM search intensive mode.
        bool          NM_NM_search_opportunistic;  ///< Flag for the NM search opportunistic mode.
        int           NM_NM_search_max_trial_pts;  ///< NM max number of trial points at each iteration.
        Double        NM_NM_search_min_simplex_vol; ///< NM min simplex volume for stopping.
        Double        NM_NM_search_include_factor; ///< NM search initial simplex inclusion factor.
        Double        NM_NM_search_rank_eps; ///< NM search epsilon for rank calculation of NM simplex.
        int           NM_NM_search_max_trial_pts_nfactor; ///< NM max number of trial points at each iteration : nfactor * dim.
        bool          NM_NM_search_use_only_Y; ///< Flag for the NM search using only Y to establish dominance of new points.
        bool          NM_NM_search_scaled_DZ; ///< Flag for the NM search using scaled DZ (Delta) for simplex characteristics.
        bool          NM_NM_search_init_Y_iter; ///< Flag for the NM search initial simplex Y to be obtained iteratively.
        bool          NM_NM_search_use_short_Y0; ///< Flag for the NM search using only two points for Y0.
        bool          NM_NM_search_init_Y_best_von; ///< Flag for the NM search picking point for init Y having best normalized volume.
        
        
        
        // Latin-Hypercube (LH) search:
        int  LH_LH_search_p0;      ///< Number of initial LH search points.   //zhenghua
        int  LH_LH_search_pi;      ///< LH search points at each iteration.   //zhenghua
        bool _opportunistic_LH;  ///< Parameter \c OPPORTUNISTIC_LH.
        bool _opp_LH_is_defined; ///< A boolean equal to \c true if a LH has been defined.
        
        bool _cache_search;               ///< Cache search.
        bool _opportunistic_cache_search; ///< Parameter \c OPPORTUNISTIC_CACHE_SEARCH.
        
        /// A boolean equal to \c true if \c OPPORTUNISTIC_CACHE_SEARCH has been defined.
        bool _opp_CS_is_defined;
        
        bool _random_eval_sort;            ///< Use random order to sort trial pts.
        
        
    public:
        
        /// Access to the \c SPECULATIVE_SEARCH parameter.
        /**
         \return The \c SPECULATIVE_SEARCH parameter.
         */
        bool get_speculative_search ( void ) const;
        
        /// Check if a model search is specified.
        /**
         \return A boolean equal to \c true if a model search is defined.
         */
        bool has_model_search ( void ) const;
        
        /// Access to the \c MODEL_SEARCH parameter.
        /**
         \param  i Index of the model search (1 or 2) -- \b IN.
         \return The \c MODEL_SEARCH parameter.
         */
        NOMAD::model_type get_model_search ( int i ) const;
        
        /// Access to the \c MODEL_SEARCH_OPTIMISTIC parameter.
        /**
         \return The \c MODEL_SEARCH_OPTIMISTIC parameter.
         */
        bool get_model_search_optimistic ( void ) const;
        
        /// Access to the \c MODEL_SEARCH_PROJ_TO_MESH parameter.
        /**
         \return The \c MODEL_SEARCH_PROJ_TO_MESH parameter.
         */
        bool get_model_search_proj_to_mesh ( void ) const;
        
        /// Access to the \c MODEL_QUAD_RADIUS_FACTOR parameter.
        /**
         \return The \c MODEL_QUAD_RADIUS_FACTOR parameter.
         */
        const NOMAD::Double & get_model_quad_radius_factor ( void ) const;
        
        /// Access to the \c MODEL_NP1_QUAD_EPSILON parameter.
        /**
         \return The \c MODEL_NP1_QUAD_EPSILON parameter.
         */
        const NOMAD::Double & get_model_np1_quad_epsilon ( void ) const;
        
        
        /// Access to the \c MODEL_QUAD_USE_WP parameter.
        /**
         \return The \c MODEL_QUAD_USE_WP parameter.
         */
        bool get_model_quad_use_WP ( void ) const;
        
        /// Access to the \c MODEL_QUAD_MAX_Y_SIZE parameter.
        /**
         \return The \c MODEL_QUAD_MAX_Y_SIZE parameter.
         */
        int get_model_quad_max_Y_size ( void ) const;
        
        /// Access to the \c MODEL_QUAD_MIN_Y_SIZE parameter.
        /**
         \return The \c MODEL_QUAD_MIN_Y_SIZE parameter.
         */
        int get_model_quad_min_Y_size ( void ) const;
        
        /// Access to the \c MODEL_SEARCH_MAX_TRIAL_PTS parameter.
        /**
         \return The \c MODEL_SEARCH_MAX_TRIAL_PTS parameter.
         */
        int get_model_search_max_trial_pts ( void ) const;
        
        /// Access to the \c MODEL_EVAL_SORT parameter.
        /**
         \return The \c MODEL_EVAL_SORT parameter.
         */
        NOMAD::model_type get_model_eval_sort ( void ) const;
        
        /// Access to the \c MODEL_EVAL_SORT_CAUTIOUS parameter.
        /**
         \return The \c MODEL_EVAL_SORT_CAUTIOUS parameter.
         */
        bool get_model_eval_sort_cautious ( void ) const;
        
        
        /// Access to all the models parameters.
        /**
         \param mp The models parameters -- \b OUT.
         */
        void get_model_parameters ( NOMAD::model_params_type & mp ) const;
        
        /// Access to the \c VNS_SEARCH parameter.
        /**
         \return The \c VNS_SEARCH parameter.
         */
        bool get_VNS_search ( void ) const;
        
        /// Access to the VNS trigger.
        /**
         \return The VNS trigger.
         */
        const NOMAD::Double & get_VNS_trigger ( void ) const;
        
        
        /// Access to the \c NM_SEARCH parameter.
        /**
         \return The \c NM_SEARCH parameter.
         */
        bool get_NM_search ( void ) const;
        
        /// Access to the \c NM_SEARCH_GAMMA parameter.
        /**
         \return The \c simplex shrink parameter.
         */
        const Double & get_NM_gamma ( void ) const;
        
        /// Access to the \c NM_SEARCH_DELTA_IC parameter.
        /**
         \return The \c simplex inside contraction parameter.
         */
        const Double & get_NM_delta_ic ( void ) const;
        
        /// Access to the \c NM_SEARCH_DELTA_OC parameter.
        /**
         \return The \c simplex outside contraction parameter.
         */
        const Double & get_NM_delta_oc ( void ) const;
        
        /// Access to the \c NM_SEARCH_DELTA_E parameter.
        /**
         \return The \c simplex expansion parameter.
         */
        const Double & get_NM_delta_e ( void ) const;
        
        /// Access to the \c NM_SEARCH_INTENSIVE parameter.
        /**
         \return The \c NM_SEARCH_INTENSIVE parameter.
         */
        bool get_NM_search_intensive ( void ) const;
        
        /// Access to the \c NM_SEARCH_OPPORTUNISTIC parameter.
        /**
         \return The \c NM_SEARCH_OPPORTUNISTIC parameter.
         */
        bool get_NM_search_opportunistic ( void ) const;
        
        /// Access to the max number of NM trial points per iteration.
        /**
         \return The max number of NM search trial points per iteration.
         */
        int get_NM_search_max_trial_pts ( void ) const;
        
        /// Access to the minimum volume of the NM simplex.
        /**
         \return The min volume of the NM simplex.
         */
        const NOMAD::Double & get_NM_search_min_simplex_vol ( void ) const;
        
        /// Access to the epsilon for rank calculation of NM simplex.
        /**
         \return The epsilon for rank calculation of NM simplex.
         */
        const NOMAD::Double & get_NM_search_rank_eps ( void ) const;
        
        
        /// Access to the inclusion factor for constructing initial NM simplex.
        /**
         \return The inclusion factor.
         */
        const NOMAD::Double & get_NM_search_include_factor ( void ) const;
        
        
        /// Access to the ratio of feasible points for producing initial NM simplex.
        /**
         \return The feasible ratio.
         */
        const NOMAD::Double & get_NM_search_init_Y_feas_ratio ( void ) const;
        
        /// Access to the \c NM_SEARCH_USE_ONLY_Y parameter.
        /**
         \return The \c NM_SEARCH_USE_ONLY_Y parameter.
         */
        bool get_NM_search_use_only_Y ( void ) const;
        
        /// Access to the \c NM_SEARCH_SCALED_DZ parameter.
        /**
         \return The \c NM_SEARCH_SCALED_DZ parameter.
         */
        bool get_NM_search_scaled_DZ ( void ) const;
        
        /// Access to the \c NM_SEARCH_INIT_Y_ITER parameter.
        /**
         \return The \c NM_SEARCH_INIT_Y_ITER parameter.
         */
        bool get_NM_search_init_Y_iter ( void ) const;
        
        
        /// Access to the \c NM_SEARCH_INIT_Y_BY_TAG parameter.
        /**
         \return The \c NM_SEARCH_INTI_Y_BY_TAG parameter.
         */
        bool get_NM_search_init_Y_by_tag ( void ) const;
        
        
        /// Access to the \c NM_SEARCH_INIT_Y_BEST_VON parameter.
        /**
         \return The \c NM_SEARCH_INTI_Y_BEST_VON parameter.
         */
        bool get_NM_search_init_Y_best_von ( void ) const;
        
        
        /// Access to the \c NM_SEARCH_USE_SHORT_Y0 parameter.
        /**
         \return The \c NM_SEARCH_USE_SHORT_Y0 parameter.
         */
        bool get_NM_search_use_short_Y0 ( void ) const;
        
        
        /// Access to the multiplicative factor to obtain the max number of NM trial points per iteration.
        /**
         \return The multiplicative factor.
         */
        int get_NM_search_max_trial_pts_nfactor ( void ) const;
        
        
        /// Access to the number of initial LH search points.
        /**
         \return The number of initial LH search points.
         */
        int get_LH_search_p0 ( void ) const;
        
        /// Access to the number of LH search points at each iteration.
        /**
         \return The number of LH search points at each iteration.
         */
        int get_LH_search_pi ( void ) const;
        
        /// Access to the \c OPPORTUNISTIC_LH parameter.
        /**
         \return The \c OPPORTUNISTIC_LH parameter.
         */
        bool get_opportunistic_LH ( void ) const;
        
        /// Access to the \c CACHE_SEARCH parameter.
        /**
         \return The \c CACHE_SEARCH parameter.
         */
        bool get_cache_search ( void ) const;
        
        /// Access to the \c DISABLE_EVAL_SORT parameter.
        /**
         \return The \c DISABLE_EVAL_SORT parameter.
         */
        bool get_disable_eval_sort ( void ) const;
        
        /// Access to the \c RANDOM_EVAL_SORT parameter.
        /**
         \return The \c RANDOM_EVAL_SORT parameter.
         */
        bool get_random_eval_sort ( void ) const;
        
        /// Access to the \c OPPORTUNISTIC_CACHE_SEARCH parameter.
        /**
         \return The \c OPPORTUNISTIC_CACHE_SEARCH parameter.
         */
        bool get_opportunistic_cache_search ( void ) const;
        
        /// Set the \c SPECULATIVE_SEARCH parameter.
        /**
         \param ss The \c SPECULATIVE_SEARCH parameter -- \b IN.
         */
        void set_SPECULATIVE_SEARCH ( bool ss );
        
        /// Disable use of models.
        /**
         */
        void set_DISABLE_MODELS ( void );
        
        /// Disable use of sort (lexicographic order used).
        /**
         */
        void set_DISABLE_EVAL_SORT ( void );
        
        /// Set random random eval sort parameter.
        /**
         \param res The \c RANDOM_EVAL_SORT parameter -- \b IN.
         */
        void set_RANDOM_EVAL_SORT ( bool res );
        
        /// Set all the models parameters.
        /**
         \param mp The models parameters -- \b IN.
         */
        void set_model_parameters ( const NOMAD::model_params_type & mp );
        
        /// Set the \c MODEL_SEARCH parameter (1/3).
        /**
         \param i  Index of the model search (1 or 2) -- \b IN.
         \param ms The \c MODEL_SEARCH parameter      -- \b IN.
         */
        void set_MODEL_SEARCH ( int i , NOMAD::model_type ms );
        
        /// Set the \c MODEL_SEARCH parameter (2/3).
        /**
         \param ms The \c MODEL_SEARCH parameter -- \b IN.
         */
        void set_MODEL_SEARCH ( bool ms );
        
        /// Set the \c MODEL_SEARCH parameter (3/3).
        /**
         \param mod_type The \c MODEL_SEARCH parameter      -- \b IN.
         */
        void set_MODEL_SEARCH ( NOMAD::model_type mod_type );
        
        /// Set the \c MODEL_SEARCH_OPTIMISTIC parameter.
        /**
         \param mso The \c MODEL_SEARCH_OPTIMISTIC parameter -- \b IN.
         */
        void set_MODEL_SEARCH_OPTIMISTIC ( bool mso );
        
        /// Set the \c MODEL_SEARCH_PROJ_TO_MESH parameter.
        /**
         \param ptm The \c MODEL_SEARCH_PROJ_TO_MESH parameter -- \b IN.
         */
        void set_MODEL_SEARCH_PROJ_TO_MESH ( bool ptm );
        
        /// Set the \c MODEL_QUAD_RADIUS_FACTOR parameter.
        /**
         \param r The \c MODEL_QUAD_RADIUS_FACTOR parameter -- \b IN.
         */
        void set_MODEL_QUAD_RADIUS_FACTOR ( const NOMAD::Double & r );
        
        /// Set the \c MODEL_QUAD_USE_WP parameter.
        /**
         \param uwp The \c MODEL_QUAD_USE_WP parameter -- \b IN.
         */
        void set_MODEL_QUAD_USE_WP ( bool uwp );
        
        /// Set the \c MODEL_QUAD_MAX_Y_SIZE parameter.
        /**
         \param s The \c MODEL_QUAD_MAX_Y_SIZE parameter -- \b IN.
         */
        void set_MODEL_QUAD_MAX_Y_SIZE ( int s );
        
        /// Set the \c MODEL_QUAD_MIN_Y_SIZE parameter.
        /**
         \param s The \c MODEL_QUAD_MIN_Y_SIZE parameter -- \b IN.
         */
        void set_MODEL_QUAD_MIN_Y_SIZE ( int s );
        
        /// Set the \c MODEL_NP1_QUAD_EPSILON parameter.
        /**
         \param r The \c MODEL_NP1_QUAD_EPSILON parameter -- \b IN.
         */
        void set_MODEL_NP1_QUAD_EPSILON ( const NOMAD::Double & r );
        
        /// Set the \c MODEL_SEARCH_MAX_TRIAL_PTS parameter.
        /**
         \param s The \c MODEL_SEARCH_MAX_TRIAL_PTS parameter -- \b IN.
         */
        void set_MODEL_SEARCH_MAX_TRIAL_PTS ( int s );
        
        /// Set the \c MODEL_EVAL_SORT parameter (1/2).
        /**
         \param mes The \c MODEL_EVAL_SORT parameter -- \b IN.
         */
        void set_MODEL_EVAL_SORT ( NOMAD::model_type mes );
        
        /// Set the \c MODEL_EVAL_SORT parameter (2/2).
        /**
         \param mes The \c MODEL_EVAL_SORT parameter -- \b IN.
         */
        void set_MODEL_EVAL_SORT ( bool mes );
        
        
        /// Set the \c MODEL_EVAL_SORT_CAUTIOUS parameter.
        /**
         \param mesc The \c MODEL_EVAL_SORT_CAUTIOUS parameter -- \b IN.
         */
        void set_MODEL_EVAL_SORT_CAUTIOUS ( bool mesc );
        
        /// Set the \c VNS_SEARCH parameter.
        /**
         \param vns The \c VNS_SEARCH parameter -- \b IN.
         */
        void set_VNS_SEARCH ( bool vns );
        
        
        /// Set the \c NM_SEARCH parameter.
        /**
         \param nm The \c NM_SEARCH parameter -- \b IN.
         */
        void set_NM_SEARCH ( bool nm );
        
        /// Set the \c NM_GAMMA parameter.
        /**
         \param ga The \c NM_GAMMA parameter -- \b IN.
         */
        void set_NM_GAMMA ( const NOMAD::Double & ga );
        
        /// Set the \c NM_DELTA_IC parameter.
        /**
         \param dic The \c NM_DELTA_IC parameter -- \b IN.
         */
        void set_NM_DELTA_IC ( const NOMAD::Double & dic );
        
        /// Set the \c NM_DELTA_OC parameter.
        /**
         \param doc The \c NM_DELTA_OC parameter -- \b IN.
         */
        void set_NM_DELTA_OC ( const NOMAD::Double & doc );
        
        /// Set the \c NM_DELTA_E parameter.
        /**
         \param de The \c NM_DELTA_E parameter -- \b IN.
         */
        void set_NM_DELTA_E ( const NOMAD::Double & de );
        
        /// Set the \c NM_SEARCH_INTENSIVE parameter.
        /**
         \param nm_intens The \c NM_SEARCH_INTENSIVE parameter -- \b IN.
         */
        void set_NM_SEARCH_INTENSIVE ( bool nm_intens );
        
        /// Set the \c NM_SEARCH_SCALED_DZ parameter.
        /**
         \param nm_scaled The \c NM_SEARCH_SCALED_DZ parameter -- \b IN.
         */
        void set_NM_SEARCH_SCALED_DZ ( bool nm_scaled );
        
        /// Set the \c NM_SEARCH_OPPORTUNISTIC parameter.
        /**
         \param nm_oppor The \c NM_SEARCH_OPPORTUNISTIC parameter -- \b IN.
         */
        void set_NM_SEARCH_OPPORTUNISTIC ( bool nm_oppor );
        
        /// Set the \c NM_SEARCH_MIN_SIMPLEX_VOL parameter.
        /**
         \param vol The \c NM_SEARCH_MIN_SIMPLEX_VOL parameter -- \b IN.
         */
        void set_NM_SEARCH_MIN_SIMPLEX_VOL ( NOMAD::Double vol );
        
        /// Set the \c NM_SEARCH_RANK_EPS parameter.
        /**
         \param eps The \c NM_SEARCH_RANK_EPS parameter -- \b IN.
         */
        void set_NM_SEARCH_RANK_EPS ( NOMAD::Double eps );
        
        /// Set the \c NM_SEARCH_INCLUDE_FACTOR parameter.
        /**
         \param f The \c NM_SEARCH_INCLUDE_FACTOR parameter -- \b IN.
         */
        void set_NM_SEARCH_INCLUDE_FACTOR ( NOMAD::Double f );
        
        /// Set the \c NM_SEARCH_INIT_Y_FEAS_RATIO parameter.
        /**
         \param f The \c NM_SEARCH_INIT_Y_FEAS_RATIO parameter -- \b IN.
         */
        void set_NM_SEARCH_INIT_Y_FEAS_RATIO ( NOMAD::Double f );
        
        /// Set the \c NM_SEARCH_INIT_Y_ITER parameter.
        /**
         \param nm_init_y_iter The \c NM_SEARCH_INIT_Y_ITER parameter -- \b IN.
         */
        void set_NM_SEARCH_INIT_Y_ITER ( bool nm_init_y_iter );
        
        /// Set the \c NM_SEARCH_MAX_TRIAL_PTS parameter.
        /**
         \param max_trial_pts The \c NM_SEARCH_MAX_TRIAL_PTS parameter -- \b IN.
         */
        void set_NM_SEARCH_MAX_TRIAL_PTS ( int max_trial_pts );
        
        
        /// Set the \c NM_SEARCH_USE_ONLY_Y parameter.
        /**
         \param nm_use_only_y The \c NM_SEARCH_USE_ONLY_Y parameter.
         */
        void set_NM_SEARCH_USE_ONLY_Y ( bool nm_use_only_y ) ;
        
        
        /// Set the \c NM_SEARCH_USE_SHORT_Y0 parameter.
        /**
         \param nm_use_short_y0 The \c NM_SEARCH_USE_SHORT_Y0 parameter.
         */
        void set_NM_SEARCH_USE_SHORT_Y0 ( bool nm_use_short_y0 ) ;
        
        
        /// Set the \c NM_SEARCH_INIT_Y_BY_TAG parameter.
        /**
         \param nm_init_y_by_tag The \c NM_SEARCH_INIT_Y_BY_TAG parameter.
         */
        void set_NM_SEARCH_INIT_Y_BY_TAG ( bool nm_init_y_by_tag ) ;
        
        
        /// Set the \c NM_SEARCH_INIT_Y_BEST_VON parameter.
        /**
         \param nm_init_y_by_best_von The \c NM_SEARCH_INIT_Y_BEST_VON parameter.
         */
        void set_NM_SEARCH_INIT_Y_BEST_VON ( bool nm_init_y_by_best_von ) ;
        
        
        /// Set the multiplicative factor to obtain the max number of NM trial points per iteration.
        /**
         \param nm_nfactor The \c NM_SEARCH_MAX_TRIAL_PTS_NFACTOR parameter.
         */
        void set_NM_SEARCH_MAX_TRIAL_PTS_NFACTOR ( int nm_nfactor ) ;
        
        
        /// Set the \c VNS_SEARCH parameter.
        /**
         \param trigger The VNS trigger -- \b IN.
         */
        void set_VNS_SEARCH ( const NOMAD::Double & trigger );
        
        /// Set the \c LH_SEARCH parameter.
        /**
         \param p0 Number of initial LH search points -- \b IN.
         \param pi LH search points at each iteration -- \b IN.
         */
        void set_LH_SEARCH ( int p0 , int pi );
        
        /// Set the \c OPPORTUNISTIC_LH parameter.
        /**
         \param olh The \c OPPORTUNISTIC_LH parameter -- \b IN.
         */
        void set_OPPORTUNISTIC_LH ( bool olh );
        
        /// Set the \c CACHE_SEARCH parameter.
        /**
         \param cs The \c CACHE_SEARCH parameter -- \b IN.
         */
        void set_CACHE_SEARCH ( bool cs );
        
        /// Set the \c OPPORTUNISTIC_CACHE_SEARCH parameter.
        /**
         \param ocs The \c OPPORTUNISTIC_CACHE_SEARCH parameter -- \b IN.
         */
        void set_OPPORTUNISTIC_CACHE_SEARCH ( bool ocs );
        
        // Mesh:
        // -----
    private:
        
        
        NOMAD::mesh_type _mesh_type;            ///< The type of mesh used (xmesh [D], gmesh, smesh [old] )
        bool          _anisotropic_mesh;        ///< Anisotropic mesh (gmesh, xmesh only)
        NOMAD::Double _anisotropy_factor; ///< Anisotropic mesh shrink/expand factor (gmesh only)
        NOMAD::Double _mesh_update_basis;        ///< Mesh update basis (tau).
        NOMAD::Double _poll_update_basis;        ///< Poll update basis (beta).
        int           _mesh_coarsening_exponent; ///< Mesh coarsening exponent.
        int           _mesh_refining_exponent;   ///< Mesh refining exponent.
        int           _initial_mesh_index;       ///< Initial mesh index (ell_0).
        NOMAD::Point  _initial_mesh_size; ///< Initial (absolute) mesh size (delta^0).
        NOMAD::Point  _min_mesh_size;     ///< Minimal (absolute) mesh size (delta_min).
        NOMAD::Point  _initial_poll_size; ///< Initial (absolute) poll size (delta^0).
        NOMAD::Point  _min_poll_size;     ///< Minimal (absolute) poll size (Delta_min).
        
        bool          _min_poll_size_defined; ///< \c true if _min_poll_size is user-defined.
        
    public:
        
        /// Access to the \c ANISOTROPIC_MESH parameter.
        /**
         \return The \c ANISOTROPIC_MESH parameter -- \b IN.
         */
        bool get_anisotropic_mesh ( void ) const;
        
        /// Access to the \c ANISOTROPY_FACTOR parameter.
        /**
         \return The \c ANISOTROPY_FACTOR parameter -- \b IN.
         */
        NOMAD::Double get_anisotropy_factor ( void ) const;
        
        /// Access to the \c MESH_TYPE parameter.
        /**
         \return The \c MESH_TYPE parameter -- \b IN.
         */
        const NOMAD::mesh_type & get_mesh_type ( void ) const;
        
        /// Access to the \c POLL_UPDATE_BASIS parameter.
        /**
         \return The \c POLL_UPDATE_BASIS parameter.
         */
        const NOMAD::Double & get_poll_update_basis ( void ) const;
        
        /// Access to the \c MESH_UPDATE_BASIS parameter.
        /**
         \return The \c MESH_UPDATE_BASIS parameter.
         */
        const NOMAD::Double & get_mesh_update_basis ( void ) const;
        
        /// Access to the \c MESH_COARSENING_EXPONENT parameter.
        /**
         \return The \c MESH_COARSENING_EXPONENT parameter.
         */
        int get_mesh_coarsening_exponent ( void ) const;
        
        /// Access to the \c MESH_REFINING_EXPONENT parameter.
        /**
         \return The \c MESH_REFINING_EXPONENT parameter.
         */
        int get_mesh_refining_exponent ( void ) const;
        
        /// Access to the \c INITIAL_MESH_INDEX parameter.
        /**
         \return The \c INITIAL_MESH_INDEX parameter.
         */
        int get_initial_mesh_index ( void ) const;
        
        /// Access to the \c INITIAL_MESH_SIZE parameter.
        /**
         \return The \c INITIAL_MESH_SIZE parameter.
         */
        const NOMAD::Point & get_initial_mesh_size ( void ) const;
        
        /// Access to the \c INITIAL_POLL_SIZE parameter.
        /**
         \return The \c INITIAL_POLL_SIZE parameter.
         */
        const NOMAD::Point & get_initial_poll_size ( void ) const;
        
        /// Access to the \c MIN_MESH_SIZE parameter.
        /**
         \return The \c MIN_MESH_SIZE parameter.
         */
        const NOMAD::Point & get_min_mesh_size ( void ) const;
        
        /// Access to the \c MIN_POLL_SIZE parameter.
        /**
         \return The \c MIN_POLL_SIZE parameter.
         */
        const NOMAD::Point & get_min_poll_size ( void ) const;
        
        /// Access to \c _min_poll_size_defined.
        /**
         \return A boolean equal to \c true if the \c MIN_POLL_SIZE parameter
         has been user-defined.
         */
        bool get_min_poll_size_defined ( void ) const;
        
        /// Set the \c ANISOTROPIC_MESH parameter.
        /**
         \param anis The \c ANISOTROPIC_MESH parameter -- \b IN.
         */
        void set_ANISOTROPIC_MESH ( bool anis );
        
        /// Set the \c ANISOTROPY_FACTOR parameter.
        /**
         \param f The \c ANISOTROPY_FACTOR parameter -- \b IN.
         */
        void set_ANISOTROPY_FACTOR ( const NOMAD::Double & f );
        
        /// Set the \c MESH_TYPE parameter.
        /**
         \param mt The \c MESH_TYPE parameter -- \b IN.
         */
        void set_MESH_TYPE ( NOMAD::mesh_type mt );
        
        /// Set the \c MESH_UPDATE_BASIS parameter.
        /**
         \param tau The \c MESH_UPDATE_BASIS parameter -- \b IN.
         */
        void set_MESH_UPDATE_BASIS ( const NOMAD::Double & tau );
        
        /// Set the \c POLL_UPDATE_BASIS parameter.
        /**
         \param r The \c POLL_UPDATE_BASIS parameter -- \b IN.
         */
        void set_POLL_UPDATE_BASIS ( const NOMAD::Double & r );
        
        /// Set the \c INITIAL_MESH_INDEX parameter.
        /**
         \param ell_0 The \c INITIAL_MESH_INDEX parameter -- \b IN.
         */
        void set_INITIAL_MESH_INDEX ( int ell_0 );
        
        /// Set the \c MESH_REFINING_EXPONENT parameter.
        /**
         \param mre The \c MESH_REFINING_EXPONENT parameter -- \b IN.
         */
        void set_MESH_REFINING_EXPONENT ( int mre );
        
        /// Set the \c MESH_COARSENING_EXPONENT parameter.
        /**
         \param mce The \c MESH_COARSENING_EXPONENT parameter -- \b IN.
         */
        void set_MESH_COARSENING_EXPONENT ( int mce );
        
        /// Set the \c MIN_MESH_SIZE parameter.
        /**
         \param mms      Minimum mesh size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_MIN_MESH_SIZE ( const NOMAD::Double & mms ,
                                bool relative = false );
        
        /// Set the \c MIN_MESH_SIZE parameter.
        /**
         \param index    Index of a variable -- \b IN.
         \param mms      Minimum mesh size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_MIN_MESH_SIZE ( int                   index ,
                                const NOMAD::Double & mms   ,
                                bool                  relative = false );
        
        /// Set the \c MIN_MESH_SIZE parameter.
        /**
         \param mms      Minimum mesh size   -- \b IN.
         \param relative Indicated as relative values
         -- \b IN -- \b optional (default = \c false).
         */
        void set_MIN_MESH_SIZE ( const NOMAD::Point & mms ,
                                bool relative = false );
        
        /// Set the \c MIN_POLL_SIZE parameter.
        /**
         \param mps      Minimum poll size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_MIN_POLL_SIZE ( const NOMAD::Double & mps ,
                                bool relative = false );
        
        /// Set the \c MIN_POLL_SIZE parameter.
        /**
         \param index    Index of a variable -- \b IN.
         \param mps      Minimum poll size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_MIN_POLL_SIZE ( int                   index ,
                                const NOMAD::Double & mps   ,
                                bool                  relative = false );
        
        /// Set the \c MIN_POLL_SIZE parameter.
        /**
         \param mps      Minimum poll size   -- \b IN.
         \param relative Indicated as relative values
         -- \b IN -- \b optional (default = \c false).
         */
        void set_MIN_POLL_SIZE ( const NOMAD::Point & mps ,
                                bool relative = false );
        
        /// Set the \c INITIAL_MESH_SIZE parameter.
        /**
         \param ims      Initial mesh size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_INITIAL_MESH_SIZE ( const NOMAD::Double & ims ,
                                    bool relative = false );
        
        /// Set the \c INITIAL_MESH_SIZE parameter.
        /**
         \param index    Index of a variable -- \b IN.
         \param ims      Initial mesh size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_INITIAL_MESH_SIZE ( int                   index ,
                                    const NOMAD::Double & ims   ,
                                    bool                  relative = false );
        
        /// Set the \c INITIAL_MESH_SIZE parameter.
        /**
         \param ims      Initial mesh size   -- \b IN.
         \param relative Indicated as relative values
         -- \b IN -- \b optional (default = \c false).
         */
        void set_INITIAL_MESH_SIZE ( const NOMAD::Point & ims ,
                                    bool relative = false );
        
        
        /// Set the \c INITIAL_POLL_SIZE parameter.
        /**
         \param ims      Initial poll size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_INITIAL_POLL_SIZE ( const NOMAD::Double & ims ,
                                    bool relative = false );
        
        /// Set the \c INITIAL_POLL_SIZE parameter.
        /**
         \param index    Index of a variable -- \b IN.
         \param ims      Initial poll size   -- \b IN.
         \param relative Indicated as a relative value
         -- \b IN -- \b optional (default = \c false).
         */
        void set_INITIAL_POLL_SIZE ( int                   index ,
                                    const NOMAD::Double & ims   ,
                                    bool                  relative = false );
        
        /// Set the \c INITIAL_POLL_SIZE parameter.
        /**
         \param ims      Initial poll size   -- \b IN.
         \param relative Indicated as relative values
         -- \b IN -- \b optional (default = \c false).
         */
        void set_INITIAL_POLL_SIZE ( const NOMAD::Point & ims , bool relative = false );
        
        // TODO update to prevent access on uncheck attributes
        // Granular variables:
        // -------------------
    private:
        
        NOMAD::Point _granularity;      ///< Granularity.
        
    public:
        
        /// Access to the granular variables.
        /**
         \return The granular variables.
         */
        const NOMAD::Point & get_granularity ( void ) const ;
        
        /// Access to the number of granular variables (granularity !=0).
        /**
         \return The number of real granular variables.
         */
        int get_number_granular_variables ( void ) const;
        
        
        /// Reset the granular variables.
        /**
         This sets all variables with a granularity of zero.
         */
        void reset_granulary ( void );
        
        /// Set the granularity for one variable.
        /**
         \param i     Index of the variable       -- \b IN.
         \param value Granularity of the variable -- \b IN.
         */
        void set_GRANULARITY ( int i , const NOMAD::Double & value );
        
        /// Set the granularity of a series of variables.
        /**
         \param fv The granular variables; This point is of dimension \c n;
         regular variables have a granularity of zero -- \b IN.
         */
        void set_GRANULARITY ( const NOMAD::Point & fv );
        
        
        /// Set the granularity of all variables.
        /**
         \param v The granularity of all variables -- \b IN.
         */
        void set_GRANULARITY ( const NOMAD::Double & v );
        
        
        
        // Directions:
        // -----------
    private:
        
        /// Types of poll directions.
        std::set<NOMAD::direction_type> _direction_types;
        
        /// Types of directions for the secondary poll.
        std::set<NOMAD::direction_type> _sec_poll_dir_types;
        
        /// Types of directions for the intensification poll.
        std::set<NOMAD::direction_type> _int_poll_dir_types;
        
        /// Change direction types to prevent using models for finding the (n+1)th direction
        ///   Apply if  ORTHO N+1 QUAD   ->    ORTHO N+1 NEG.
        void set_DIRECTION_TYPE_NO_MODEL ( void );
        
        
    public:
        
        /// Access to the list of poll direction types.
        /**
         \return The list of poll direction types.
         */
        const std::set<NOMAD::direction_type> & get_direction_types ( void ) const;
        
        /// Access to the list of secondary poll direction types.
        /**
         \return The list of secondary poll direction types.
         */
        const std::set<NOMAD::direction_type> & get_sec_poll_dir_types ( void ) const;
        
        /// Access to the list of intensification poll direction types.
        /**
         \return The list of intensification poll direction types.
         */
        const std::set<NOMAD::direction_type> & get_int_poll_dir_types ( void ) const;
        
        
        /// Check if there are Ortho-MADS directions.
        /**
         \return A boolean equal to \c true if there is at least one
         Ortho-MADS direction in the list of poll directions
         or in the list of secondary poll directions.
         */
        bool has_orthomads_directions ( void ) const;
        
        
        /// Check if a direction type that required dynamic completion for the (n+1)th direction.
        /**
         \return true if a dynamic completion is required.
         */
        bool has_dynamic_direction(void) const;
        
        /// Reset the directions.
        /**
         */
        void reset_directions ( void );
        
        
        /// Add a new direction type.
        /**
         \param dt The new direction type -- \b IN.
         */
        void set_DIRECTION_TYPE ( NOMAD::direction_type dt );
        
        /// Add a set of new direction types.
        /**
         \param dt The set of new direction types -- \b IN.
         */
        void set_DIRECTION_TYPE ( const std::set<NOMAD::direction_type> & dt );
        
        
        /// Add a new direction type for the secondary poll.
        /**
         \param dt The new direction type -- \b IN.
         */
        void set_SEC_POLL_DIR_TYPE ( NOMAD::direction_type dt );
        
        /// Add a set of new direction types for the secondary poll.
        /**
         \param dt The set of new direction types -- \b IN.
         */
        void set_SEC_POLL_DIR_TYPE ( const std::set<NOMAD::direction_type> & dt );
        
        
        /// Add a new direction type for the intensification poll.
        /**
         \param dt The new direction type -- \b IN.
         */
        void set_INT_POLL_DIR_TYPE ( NOMAD::direction_type dt );
        
        /// Add a set of new direction types for the intensification poll.
        /**
         \param dt The set of new direction types -- \b IN.
         */
        void set_INT_POLL_DIR_TYPE ( const std::set<NOMAD::direction_type> & dt );
        
        
        /// Enables use of quad model to determine prospect direction
        /// for Ortho n+1 direction type
        /**
         \param qmpd boolean -- \b IN.
         */
        void set_QUAD_MODEL_PROSPECT_DIR ( bool qmpd ) ;
        
        // Starting point(s):
        // ------------------
    private:
        
        std::vector<NOMAD::Point *> _x0s; ///< List of starting points.
        std::string       _x0_cache_file; ///< Cache file containing starting points.
        
    public:
        
        /// Add a new point in the list of starting points.
        /**
         \param x0 The new point -- \b IN.
         */
        void set_X0 ( const NOMAD::Point & x0 );
        
        /// Indicate a cache file containing starting points.
        /**
         \param file_name Name of the cache file -- \b IN.
         */
        void set_X0 ( const std::string  & file_name );
        
        /// Reset all starting points.
        void reset_X0 ( void );
        
        /// Access to the list of starting points.
        /**
         \return The list of starting points.
         */
        const std::vector<NOMAD::Point *> & get_x0s ( void ) const;
        
        /// Access to the name of a cache file containing starting points.
        /**
         \return The file name.
         */
        const std::string & get_x0_cache_file ( void ) const;
        
        // Signature: standard or extern (only one is != NULL):
        // ----------------------------------------------------
    private:
        
        /// Standard signature.
        /**
         Created and deleted by the Parameters class.
         */
        NOMAD::Signature * _std_signature;
        
        /// Extern signature.
        /**
         Created and deleted outside the class.
         */
        NOMAD::Signature * _extern_signature;
        
    public:
        
        /// Access to the signature.
        /**
         \return The one non-NULL signature.
         */
        NOMAD::Signature * get_signature ( void ) const;
        
        /// Set a new extern signature.
        /**
         Deletes the standard signature.
         \param s A pointer to the extern signature -- \b IN.
         */
        void set_EXTERN_SIGNATURE ( NOMAD::Signature * s );
        
        
        //        /// Reset standard and extern signatures.
        //        /**
        //         Deletes the standard and extern signatures.
        //        */
        //        void reset_signatures ( void );
        
        
        
        // Dimension:
        // ----------
    private:
        
        /// Dimension.
        /**
         - Number of variables.
         - Parameter \c DIMENSION.
         */
        int _dimension;
        
    public:
        
        /// Access to the dimension.
        /**
         \return The dimension.
         */
        int get_dimension ( void ) const;
        
        /// Set the dimension.
        /**
         \param n The dimension -- \b IN.
         \return \c true if the operation succeeded.
         */
        bool set_DIMENSION ( int  n );
        
        // Fixed variables:
        // ----------------
    private:
        
        /// Fixed variables.
        /**
         - This point is of dimension \c n (parameter \c DIMENSION).
         - Undefined values correspond to free variables.
         */
        NOMAD::Point _fixed_variables;
        
        /// Number of free variables.
        int _nb_free_variables;
        
    public:
        
        /// Access to the number of free variables.
        /**
         \return The number of free variables.
         */
        int get_nb_free_variables ( void ) const;
        
        /// Access to the fixed variables.
        /**
         \return The fixed variables.
         */
        const NOMAD::Point & get_fixed_variable ( void ) const
        {
            return get_fixed_variables();
        }
        
        /// Access to the fixed variables.
        /**
         \return The fixed variables.
         */
        const NOMAD::Point & get_fixed_variables ( void ) const;
        
        /// Test if a variable is fixed.
        /**
         \param i Index of the variable -- \b IN.
         \return A boolean equal to \c true if the variable \c i is fixed.
         */
        bool variable_is_fixed ( int i ) const;
        
        /// Reset the fixed variables.
        /**
         This frees all the variables.
         */
        void reset_fixed_variables ( void );
        
        /// Fix one variable.
        /**
         \param i     Index of the variable       -- \b IN.
         \param value Value of the fixed variable -- \b IN.
         */
        void set_FIXED_VARIABLE ( int i , const NOMAD::Double & value );
        
        /// Fix one variable.
        /**
         The value of the variable is based on the starting point.
         \param i Index of the variable -- \b IN.
         */
        void set_FIXED_VARIABLE ( int i );
        
        /// Fix a series of variables.
        /**
         \param fv The fixed variables; This point is of dimension \c n;
         free variables correspond to undefined values -- \b IN.
         */
        void set_FIXED_VARIABLE ( const NOMAD::Point & fv );
        
        /// Free a variable.
        /**
         \param i Index of the variable -- \b IN.
         */
        void set_FREE_VARIABLE ( int i ) { set_FIXED_VARIABLE ( i , NOMAD::Double() ); }
        
        // Periodic variables:
        // -------------------
        
        /// Periodic variables.
        /**
         - This vector is of size \c n.
         - \c _periodic_variables[i] is equal to \c true if
         the variable \c i is periodic.
         */
        std::vector<bool> _periodic_variables;
        
        /// Access to the periodic variables.
        /**
         \return The periodic variables
         */
        const std::vector<bool> & get_periodic_variable ( void ) const
        {
            return get_periodic_variables();
        }
        
        /// Access to the periodic variables.
        /**
         \return The periodic variables
         */
        const std::vector<bool> & get_periodic_variables ( void ) const;
        
        /// Check if there is periodic variables.
        /**
         \return A boolean equal to \c true if there is periodic variables.
         */
        bool has_periodic_variables ( void ) const
        {
            return !_periodic_variables.empty();
        }
        
        /// Reset periodic variables.
        void reset_periodic_variables ( void );
        
        /// Set one variable to be periodic.
        /**
         \param i Index of the variable -- \b IN.
         */
        void set_PERIODIC_VARIABLE ( int i );
        
        /// Set a series of variables to be periodic.
        /**
         \param pv Vector of size \c n indicating the
         periodic variables -- \b IN.
         */
        void set_PERIODIC_VARIABLE ( const std::vector<bool> & pv );
        
        // Categorical variables:
        // ----------------------
    private:
        
        /// Extended poll trigger.
        /**
         - Must be strictly positive.
         - May be relative.
         */
        NOMAD::Double _extended_poll_trigger;
        
        /// Equal to \c true if the extended poll trigger is relative.
        bool _relative_ept;
        
        /// Equal to \c true if the extended poll is enabled.
        bool _extended_poll_enabled;
        
        /// Neighborhood executable for using batch mode with categorical variables.
        std::string _neighbors_exe;
        
    public:
        
        /// Access to the extended poll trigger.
        /**
         \return The extended poll trigger.
         */
        const NOMAD::Double & get_extended_poll_trigger ( void ) const;
        
        /// Check if the extended poll trigger is relative.
        /**
         \return A boolean equal to \c true if the extended poll trigger is relative.
         */
        bool get_relative_ept ( void ) const;
        
        /// Check if the extended poll is enabled.
        /**
         \return A boolean equal to \c true if the extended poll is enabled.
         */
        bool get_extended_poll_enabled ( void ) const;
        
        /// Access to the neighborhood executable.
        /**
         \return The neighborhood executable.
         */
        const std::string & get_neighbors_exe ( void ) const;
        
        /// Set the extended poll trigger.
        /**
         \param ept The extended poll trigger -- \b IN.
         \param rel A boolean equal to \c true if the extended poll
         trigger is relative -- \b IN.
         */
        void set_EXTENDED_POLL_TRIGGER ( const NOMAD::Double & ept , bool rel );
        
        /// Enable or disable the extended poll.
        /**
         \param epe A boolean equal to \c true if the extended poll is enabled -- \b IN.
         */
        void set_EXTENDED_POLL_ENABLED ( bool epe );
        
        /// Set the neighborhood executable.
        /**
         \param ne The neighborhood executable.
         */
        void set_NEIGHBORS_EXE ( const std::string & ne );
        
        // Groups of variables:
        // --------------------
    private:
        
        /// Groups of variables.
        std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> _var_groups;
        
        /// User groups of variables.
        std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> _user_var_groups;
        
        /// Reset a group of variables.
        /**
         \param g Group to reset; may be \c _var_groups or \c _user_var_groups
         -- \b IN/OUT.
         */
        void reset_variable_groups ( std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> & g ) const;
        
    public:
        
        /// Reset groups of variables.
        void reset_variable_groups ( void );
        
        /// Access to the groups of variables.
        /**
         \return The groups of variables.
         */
        const std::set<NOMAD::Variable_Group*,NOMAD::VG_Comp> &
        get_variable_groups ( void ) const;
        
        /// Set one group of variables.
        /**
         Basic version.
         \param var_indexes Indexes of the variables of the group -- \b IN.
         */
        void set_VARIABLE_GROUP ( const std::set<int> & var_indexes );
        
        /// Set one group of variables.
        /**
         Advanced version.
         \param var_indexes         Indexes of the variables of the group  -- \b IN.
         \param prim_poll_dir_types Types of the poll directions           -- \b IN.
         \param sec_poll_dir_types  Types of the secondary poll directions -- \b IN.
         \param int_poll_dir_types  Types of the intensification poll directions -- \b IN.
         */
        void set_VARIABLE_GROUP ( const std::set<int>                   & var_indexes         ,
                                 const std::set<NOMAD::direction_type> & prim_poll_dir_types ,
                                 const std::set<NOMAD::direction_type> & sec_poll_dir_types ,
                                 const std::set<NOMAD::direction_type> & int_poll_dir_types  );
        
        /// Set several groups of variables.
        /**
         \param vg A list of groups of variables -- \b IN.
         */
        void set_VARIABLE_GROUP ( const std::list<NOMAD::Variable_Group*> & vg );
        
        // Bounds and scaling:
        // -------------------
    private:
        
        NOMAD::Point _lb;      ///< Lower bounds.
        NOMAD::Point _ub;      ///< Upper bounds.
        NOMAD::Point _scaling; ///< Scaling.
        
    public:
        
        /// Access to the lower bounds.
        /**
         \return The lower bounds.
         */
        const NOMAD::Point & get_lb ( void ) const;
        
        /// Access to the upper bounds.
        /**
         \return The upper bounds.
         */
        const NOMAD::Point & get_ub ( void ) const;
        
        /// Access to the lower bounds.
        /**
         \return The lower bounds.
         */
        const NOMAD::Point & get_lower_bound ( void ) const { return get_lb(); }
        
        /// Access to the upper bounds.
        /**
         \return The upper bounds.
         */
        const NOMAD::Point & get_upper_bound ( void ) const { return get_ub(); }
        
        /// Access to the scaling.
        /**
         \return The scaling.
         */
        const NOMAD::Point & get_scaling ( void ) const;
        
        /// Reset the bounds.
        void reset_bounds ( void );
        
        /// Reset the scaling.
        void reset_scaling ( void );
        
        /// Set one lower bound.
        /**
         \param i  Index of the variable -- \b IN.
         \param lb Lower bound           -- \b IN.
         */
        void set_LOWER_BOUND ( int i , const NOMAD::Double & lb );
        
        /// Set all lower bounds.
        /**
         Use undefined values for variables without bounds.
         \param lb Lower bounds -- \b IN.
         */
        void set_LOWER_BOUND ( const NOMAD::Point & lb );
        
        /// Set one upper bound.
        /**
         \param i  Index of the variable -- \b IN.
         \param ub Upper bound           -- \b IN.
         */
        void set_UPPER_BOUND ( int i , const NOMAD::Double & ub );
        
        /// Set all upper bounds.
        /**
         Use undefined values for variables without bounds.
         \param ub Upper bounds -- \b IN.
         */
        void set_UPPER_BOUND ( const NOMAD::Point & ub );
        
        /// Set the scaling of one variable.
        /**
         \param i Index of the variable -- \b IN.
         \param s Scaling               -- \b IN.
         */
        void set_SCALING ( int i , const NOMAD::Double & s );
        
        /// Set the scaling for all variables.
        /**
         Use undefined values for variables scaling.
         \param s Scaling -- \b IN.
         */
        void set_SCALING ( const NOMAD::Point & s );
        
        // Blackboxes (inputs and outputs):
        // --------------------------------
    private:
        
        /// Blackbox input types.
        std::vector<NOMAD::bb_input_type> _bb_input_type;
        
        /// Parameter \c BB_INPUT_INCLUDE_TAG.
        bool _bb_input_include_tag;
        
        /// Parameter \c BB_INPUT_INCLUDE_SEED.
        bool _bb_input_include_seed;
        
        /// Blackbox output types.
        /**
         May be modified during the algorithm (PEB constraints).
         */
        mutable std::vector<NOMAD::bb_output_type> _bb_output_type;
        
        /// List of blackbox executables.
        std::list<std::string> _bb_exe;
        
        /// List of objective indexes.
        std::list<int> _index_obj;
        
        /// Index for the output \c STAT_SUM.
        int _index_stat_sum;
        
        /// Index for the output \c STAT_AVG.
        int _index_stat_avg;
        
        /// Index for the output \c CNT_EVAL.
        /**
         - \c 0 or \c 1.
         - \c CNT_EVAL indicates if a blackbox evaluation has to be counted.
         */
        int _index_cnt_eval;
        
        /// Parameter \c BB_REDIRECTION.
        bool _bb_redirection;
        
    public:
        
        /// Access to the number of blackbox outputs.
        /**
         \return The number of blackbox outputs.
         */
        int get_bb_nb_outputs ( void ) const;
        
        /// Access to the \c BB_INPUT_INCLUDE_TAG parameter.
        /**
         \return The \c BB_INPUT_INCLUDE_TAG parameter.
         */
        bool get_bb_input_include_tag ( void ) const;
        
        /// Access to the \c BB_INPUT_INCLUDE_SEED parameter.
        /**
         \return The \c BB_INPUT_INCLUDE_SEED parameter.
         */
        bool get_bb_input_include_seed ( void ) const;
        
        /// Access to the blackbox input types.
        /**
         \return The blackbox input types.
         */
        const std::vector<NOMAD::bb_input_type> & get_bb_input_type ( void ) const;
        
        /// Access to the blackbox output types.
        /**
         \return The blackbox output types.
         */
        const std::vector<NOMAD::bb_output_type> & get_bb_output_type ( void ) const;
        
        /// Access to the objective indexes.
        /**
         \return The list of objective indexes.
         */
        const std::list<int> & get_index_obj ( void ) const;
        
        /// Access to the number of objective functions.
        /**
         \return The number of objective functions.
         */
        int get_nb_obj ( void ) const;
        
        /// Check the display and file stats.
        /**
         \param stats The stats -- \b IN.
         \return A booleam equal to \c true if the stats are valid.
         */
        bool check_display_stats ( const std::list<std::string> & stats ) const;
        
        /// Check if there is a \c STAT_SUM output.
        /**
         \return A boolean equal to \c true if there a \c STAT_SUM output.
         */
        bool check_stat_sum ( void ) const;
        
        /// Check if there is a \c STAT_AVG output.
        /**
         \return A boolean equal to \c true if there is a \c STAT_AVG output.
         */
        bool check_stat_avg ( void ) const;
        
        /// Access the index of output \c CNT_EVAL output.
        /**
         \return Index of the \c CNT_EVAL output.
         \return \c -1 if there is no \c CNT_EVAL output.
         */
        int get_index_cnt_eval ( void ) const;
        
        /// Access the index of output \c STAT_SUM output.
        /**
         \return Index of the \c STAT_SUM output.
         \return \c -1 if there is no \c STAT_SUM output.
         */
        int get_index_stat_sum ( void ) const;
        
        /// Access the index of output \c STAT_AVG output.
        /**
         \return Index of the \c STAT_AVG output.
         \return \c -1 if there is no \c STAT_AVG output.
         */
        int get_index_stat_avg ( void ) const;
        
        /// Access to the list of blackbox executables.
        /**
         \return The list of blackbox executables.
         */
        const std::list<std::string> & get_bb_exe ( void ) const;
        
        /// Access to the \c BB_REDIRECTION parameter.
        /**
         \return The \c BB_REDIRECTION parameter.
         */
        bool get_bb_redirection ( void ) const;
        
        /// Set the \c BB_INPUT_INCLUDE_TAG parameter.
        /**
         \param bbiit The \c BB_INPUT_INCLUDE_TAG parameter -- \b IN.
         */
        void set_BB_INPUT_INCLUDE_TAG ( bool bbiit );
        
        /// Set the \c BB_INPUT_INCLUDE_SEED parameter.
        /**
         \param bbiis The \c BB_INPUT_INCLUDE_SEED parameter -- \b IN.
         */
        void set_BB_INPUT_INCLUDE_SEED ( bool bbiis );
        
        /// Set the blackbox input type of one variable.
        /**
         \param index Index of the variable -- \b IN.
         \param bbit  Type of the variable  -- \b IN.
         */
        void set_BB_INPUT_TYPE ( int index , NOMAD::bb_input_type bbit  );
        
        /// Set the blackbox input types of all variables.
        /**
         \param bbit  Types of the variables -- \b IN.
         */
        void set_BB_INPUT_TYPE ( const std::vector<NOMAD::bb_input_type> & bbit );
        
        /// Set the blackbox input types of all variables.
        /**
         \param bbit Types of the variables -- \b IN.
         */
        void set_BB_INPUT_TYPE ( const std::list<NOMAD::bb_input_type> & bbit );
        
        /// Set the blackbox output types.
        /**
         \param bbot Blackbox output types -- \b IN.
         */
        void set_BB_OUTPUT_TYPE ( const std::list<NOMAD::bb_output_type> & bbot );
        
        /// Set the blackbox output types.
        /**
         \param bbot Blackbox output types -- \b IN.
         */
        void set_BB_OUTPUT_TYPE ( const std::vector<NOMAD::bb_output_type> & bbot );
        
        /// Set a list of blackbox executable names.
        /**
         Must correspond to the blackbox output types.
         \param bbexe The list of blackbox executable names -- \b IN.
         */
        void set_BB_EXE ( const std::list<std::string> & bbexe );
        
        /// Set a list of blackbox executable names.
        /**
         \param m     Number of blackbox outputs and size of \c bbexe -- \b IN.
         \param bbexe The list of blackbox executable names           -- \b IN.
         */
        void set_BB_EXE ( int m , const std::string * bbexe );
        
        /// Set a unique blackbox executable name.
        /**
         \param bbexe The blackbox executable name -- \b IN.
         */
        void set_BB_EXE ( const std::string & bbexe );
        
        /// Set the \c BB_REDIRECTION parameter.
        /**
         \param bbr The \c BB_REDIRECTION parameter -- \b IN.
         */
        void set_BB_REDIRECTION ( bool bbr );
        
        /// Reset the PEB statuses.
        /**
         Set all outputs at NOMAD::PEB_E (PEB constraint in "extreme" status)
         to NOMAD::PEB_P (PEB constraint in "progressive" status).
         */
        void reset_PEB_changes ( void ) const;
        
        /// Change constraints from PEB to PB.
        /**
         Set all outputs at NOMAD::PEB to NOMAD::PB
         */
        void change_PEB_to_PB ( void );
        
        /// Change the status of one PEB constraint.
        /**
         The status is changed from
         NOMAD::PEB_P (PEB constraint in "progressive" status)
         to NOMAD::PEB_E (PEB constraint in "extreme" status).
         \param index Index of the PEB constraint.
         */
        void change_PEB_constraint_status ( int index ) const;
        
        // Trend matrix:
        // ------------
    private:
        
        /// Trend matrix provided as a vector of points.
        std::vector<NOMAD::Point> _trend_matrix;
        
        /// Flag equal to \c true if constraint trend matrix used to sort evaluation points.
        bool _trend_matrix_eval_sort;
        
        /// Flag equal to \c true if trend matrix used to perform a basic line search.
        bool _trend_matrix_basic_line_search;
        
        
    public:
        
        /// Push back a new trend point in the trend matrix. Order must follow bbot order.
        /**
         \param T The new trend point -- \b IN.
         */
        void push_back_trend ( const NOMAD::Point & T ) ;
        
        /// Reset trend matrix.
        void reset_trend_matrix ( void );
        
        /// Get the trend matrix.
        /**
         \return The trend matrix.
         */
        const std::vector<NOMAD::Point> & get_trend_matrix ( void ) const { return _trend_matrix ; }
        
        /// Get the flag for trend matrix basic line search.
        /**
         \return A boolean equal to \c true if trend matrix
         basic line search is enabled.
         */
        bool get_trend_matrix_basic_line_search ( void ) const { return _trend_matrix_basic_line_search ; }
        
        
        /// Set the \c TREND_MATRIX_BASIC_LINE_SEARCH bool parameter.
        /**
         \param t                    -- \b IN.
         */
        void set_TREND_MATRIX_BASIC_LINE_SEARCH ( bool t ) ;
        
        /// Get the flag for trend matrix eval sort.
        /**
         \return A boolean equal to \c true if trend matrix
         eval sort is enabled.
         */
        bool get_trend_matrix_eval_sort ( void ) const { return _trend_matrix_eval_sort ; }
        
        
        /// Set the \c TREND_MATRIX_EVAL_SORT bool parameter.
        /**
         \param t                    -- \b IN.
         */
        void set_TREND_MATRIX_EVAL_SORT ( bool t );
        
        
        // Surrogates:
        // -----------
    private:
        
        /// Surrogate executables.
        /**
         \c _sgte_exe[bb_exe] corresponds to the surrogate associated
         with the blackbox true executable \c bb_exe.
         */
        std::map<std::string,std::string> _sgte_exe;
        
        bool _disable_eval_sort;  ///< Sort disablement
        
        /// Flag equal to \c true if surrogates are used to sort evaluation points.
        bool _sgte_eval_sort;
        
        /// Flag equal to \c true if the problem has a surrogate.
        bool _has_sgte;
        
        /// Flag equal to \c true if NOMAD considers only surrogates.
        bool _opt_only_sgte;
        
        /// Surrogate cost.
        /**
         Number of surrogate evaluations counting as one blackbox evaluation.
         */
        int _sgte_cost;
        
        /// Maximum number of surrogate evaluations.
        int _sgte_max_eval;
        
        /// Surrogate cache file.
        std::string _sgte_cache_file;
        
        
        
    public:
        
        /// Access to the surrogate associated with a truth executable.
        /**
         \param bb_exe The truth executable -- \b IN.
         \return The surrogate executable name.
         \return An empty string if \c bb_exe has no surrogate.
         */
        std::string get_sgte_exe ( const std::string & bb_exe ) const;
        
        /// Access to the \c SGTE_EVAL_SORT parameter.
        /**
         \return The \c SGTE_EVAL_SORT parameter.
         */
        bool get_sgte_eval_sort ( void ) const;
        
        /// Access to the \c SGTE_COST parameter.
        /**
         \return The \c SGTE_COST parameter.
         */
        int get_sgte_cost ( void ) const;
        
        /// Access to the \c MAX_SGTE_EVAL parameter.
        /**
         \return The \c MAX_SGTE_EVAL parameter.
         */
        int get_max_sgte_eval ( void ) const;
        
        /// Access to the \c SGTE_CACHE_FILE parameter.
        /**
         \return The \c SGTE_CACHE_FILE parameter.
         */
        const std::string & get_sgte_cache_file ( void ) const;
        
        /// Access to the \c HAS_SGTEparameter.
        /**
         \return The \c HAS_SGTEparameter.
         */
        bool has_sgte ( void ) const;
        
        /// Access to the \c OPT_ONLY_SGTE parameter.
        /**
         \return The \c OPT_ONLY_SGTE parameter.
         */
        bool get_opt_only_sgte ( void ) const;
        
        /// Check if a surrogate executable has been defined.
        /**
         \return A boolean equal to \c true if
         a surrogate executable has been defined.
         */
        bool has_sgte_exe ( void ) const;
        
        /// Set the \c SGTE_EXE parameter.
        /**
         \param bb_exe   Truth executable                -- \b IN.
         \param sgte_exe Associated surrogate executable -- \b IN.
         */
        void set_SGTE_EXE ( const std::string & bb_exe , const std::string & sgte_exe );
        
        /// Set the \c SGTE_EVAL_SORT parameter.
        /**
         \param ses The \c SGTE_EVAL_SORT parameter -- \b IN.
         */
        void set_SGTE_EVAL_SORT ( bool ses );
        
        /// Set the \c HAS_SGTE parameter.
        /**
         \param hs The \c HAS_SGTE parameter -- \b IN.
         */
        void set_HAS_SGTE ( bool hs );
        
        /// Set the \c OPT_ONLY_SGTE parameter.
        /**
         \param oos The \c OPT_ONLY_SGTE parameter -- \b IN.
         */
        void set_OPT_ONLY_SGTE ( bool oos );
        
        /// Set the \c SGTE_COST parameter.
        /**
         \param sc The \c SGTE_COST parameter -- \b IN.
         */
        void set_SGTE_COST ( int sc );
        
        /// Set the \c MAX_SGTE_EVAL parameter.
        /**
         \param mse The \c MAX_SGTE_EVAL parameter -- \b IN.
         */
        void set_MAX_SGTE_EVAL ( int mse );
        
        /// Set the \c SGTE_CACHE_FILE parameter.
        /**
         \param scf The \c SGTE_CACHE_FILE parameter -- \b IN.
         */
        void set_SGTE_CACHE_FILE ( const std::string & scf );
        
        // Barrier:
        // --------
    private:
        NOMAD::Double     _h_min;   ///< Value of \c h_min.
        NOMAD::Double     _h_max_0; ///< Initial value of \c h_max.
        NOMAD::Double     _rho;     ///< Rho parameter of the progressive barrier.
        NOMAD::hnorm_type _h_norm;  ///< Norm used to compute the feasibility function \c h.
        
        /// Flag equal to \c true if there are constraints.
        bool _has_constraints;
        
        /// Flag equal to \c true if there are filter or progressive barrier constraints.
        bool _has_filter_constraints;
        
        /// Flag equal to \c true if there are extreme barrier constraints.
        bool _has_EB_constraints;
        
        /// Type of the barrier.
        /**
         May be NOMAD::FILTER, NOMAD::PB, NOMAD::PEB_P, or NOMAD::EB.
         */
        NOMAD::bb_output_type _barrier_type;
        
    public:
        
        /// Access to the \c H_MIN parameter.
        /**
         \return The \c H_MIN parameter.
         */
        const NOMAD::Double & get_h_min ( void ) const;
        
        /// Access to the \c H_MAX_0 parameter.
        /**
         \return The \c H_MAX_0 parameter.
         */
        const NOMAD::Double & get_h_max_0 ( void ) const;
        
        /// Access to the \c RHO parameter.
        /**
         \return The \c RHO parameter.
         */
        const NOMAD::Double & get_rho ( void ) const;
        
        /// Access to the \c H_NORM parameter.
        /**
         \return The \c H_NORM parameter.
         */
        NOMAD::hnorm_type get_h_norm ( void ) const;
        
        /// Access to the type of barrier.
        /**
         \return The type of barrier.
         */
        NOMAD::bb_output_type get_barrier_type ( void ) const;
        
        /// Check if there are constraints.
        /**
         \return A boolean equal to \c true if there are constraints.
         */
        bool has_constraints ( void ) const;
        
        /// Check if there are extreme barrier constraints.
        /**
         \return A boolean equal to \c true if there are extreme barrier constraints.
         */
        bool has_EB_constraints ( void ) const;
        
        /// Check if the progressive barrier (PB) is used.
        /**
         \return A boolean equal to \c true if the PB is used.
         */
        bool use_sec_poll_center ( void ) const;
        
        /// Set the \c H_MIN parameter.
        /**
         \param h_min The \c H_MIN parameter -- \b IN.
         */
        void set_H_MIN ( const NOMAD::Double & h_min );
        
        /// Set the \c H_MAX_0 parameter.
        /**
         \param h_max The \c H_MAX_0 parameter -- \b IN.
         */
        void set_H_MAX_0 ( const NOMAD::Double & h_max );
        
        /// Set the \c RHO parameter.
        /**
         \param rho The \c RHO parameter -- \b IN.
         */
        void set_RHO ( const NOMAD::Double & rho );
        
        /// Set the \c H_NORM parameter.
        /**
         \param h_norm The \c H_NORM parameter -- \b IN.
         */
        void set_H_NORM ( NOMAD::hnorm_type h_norm );
        
        // MULTI-MADS parameters:
        // ----------------------
    private:
        
        /// Number of MADS runs in Multi-MADS.
        int _multi_nb_mads_runs;
        
        /// Maximum number of blackbox evaluations in Multi-MADS.
        int _multi_overall_bb_eval;
        
        /// Flag equal to \c true if the delta criterion is used in Multi-MADS.
        bool _multi_use_delta_crit;
        
        /// Bounds on the objective necessary for the display of the \c surf stat.
        NOMAD::Point _multi_f_bounds;
        
        /// Multi-MADS reformulation.
        /**
         May be NOMAD::NORMALIZED, NOMAD::PRODUCT,
         NOMAD::DIST_L1, NOMAD::DIST_L2, or NOMAD::DIST_LINF.
         */
        NOMAD::multi_formulation_type _multi_formulation;
        
    public:
        
        /// Access to the \c MULTI_NB_MADS_RUNS parameter.
        /**
         \return The \c MULTI_NB_MADS_RUNS parameter.
         */
        int get_multi_nb_mads_runs ( void ) const;
        
        /// Access to the \c MULTI_OVERALL_BB_EVAL parameter.
        /**
         \return The \c MULTI_OVERALL_BB_EVAL parameter.
         */
        int get_multi_overall_bb_eval ( void ) const;
        
        /// Access to the \c MULTI_USE_DELTA_CRIT parameter.
        /**
         \return The \c MULTI_USE_DELTA_CRIT parameter.
         */
        bool get_multi_use_delta_crit ( void ) const;
        
        /// Access to the \c MULTI_F_BOUNDS parameter.
        /**
         \return The \c MULTI_F_BOUNDS parameter.
         */
        const NOMAD::Point & get_multi_f_bounds ( void ) const;
        
        /// Access to the \c MULTI_FORMULATION parameter.
        /**
         \return The \c MULTI_FORMULATION parameter.
         */
        NOMAD::multi_formulation_type get_multi_formulation ( void ) const;
        
        /// Set the \c MULTI_NB_MADS_RUNS parameter.
        /**
         \param mads_runs The \c MULTI_NB_MADS_RUNS parameter -- \b IN.
         */
        void set_MULTI_NB_MADS_RUNS ( int mads_runs );
        
        /// Set the \c MULTI_OVERALL_BB_EVAL parameter.
        /**
         \param bbe The \c MULTI_OVERALL_BB_EVAL parameter -- \b IN.
         */
        void set_MULTI_OVERALL_BB_EVAL ( int bbe );
        
        /// Set the \c MULTI_USE_DELTA_CRIT parameter.
        /**
         \param udc The \c MULTI_USE_DELTA_CRIT parameter -- \b IN.
         */
        void set_MULTI_USE_DELTA_CRIT ( bool udc );
        
        /// Set the \c MULTI_F_BOUNDS parameter.
        /**
         \param mfb The \c MULTI_F_BOUNDS parameter -- \b IN.
         */
        void set_MULTI_F_BOUNDS ( const NOMAD::Point & mfb );
        
        /// Set the \c MULTI_FORMULATION parameter.
        /**
         \param mf The \c MULTI_FORMULATION parameter -- \b IN.
         */
        void set_MULTI_FORMULATION ( NOMAD::multi_formulation_type mf );
        
        // Opportunistic strategy parameters:
        // ----------------------------------
    private:
        
        /// Flag equal to \c true if the opportunistic strategy is enabled.
        bool _opportunistic_eval;
        
        /// Minimum number of successes.
        /**
         Parameter \c OPPORTUNISTIC_MIN_NB_SUCCESS.
         */
        int _opportunistic_min_nb_success;
        
        /// Minimum number of evaluations.
        /**
         Parameter \c OPPORTUNISTIC_MIN_EVAL.
         */
        int _opportunistic_min_eval;
        
        /// Minimum (relative) percentage of feasible objective improvement.
        /**
         Parameter \c OPPORTUNISTIC_MIN_F_IMPRVMT.
         */
        NOMAD::Double _opportunistic_min_f_imprvmt;
        
        /// Flag equal to \c true if the lucky evaluation is enabled.
        /**
         Do one more eval "for luck".
         */
        bool _opportunistic_lucky_eval;
        
        
        /// Max block size for list evaluation
        /**
         Parameter \c BB_MAX_BLOCK_SIZE.
         */
        int _bb_max_block_size;
        
        /// Block of points evaluation
        /**
         Parameter \c EVAL_POINTS_AS_BLOCK.
         */
        bool _eval_points_as_block;
        
    public:
        
        /// Access to the \c OPPORTUNISTIC_EVAL parameter.
        /**
         \return The \c OPPORTUNISTIC_EVAL parameter
         */
        bool get_opportunistic_eval ( void ) const;
        
        /// Access to the \c OPPORTUNISTIC_MIN_NB_SUCCESS parameter.
        /**
         \return The \c OPPORTUNISTIC_MIN_NB_SUCCESS parameter.
         */
        int get_opportunistic_min_nb_success ( void ) const;
        
        /// Access to the \c OPPORTUNISTIC_MIN_EVAL parameter.
        /**
         \return The \c OPPORTUNISTIC_MIN_EVAL parameter.
         */
        int get_opportunistic_min_eval ( void ) const;
        
        /// Access to the \c OPPORTUNISTIC_MIN_F_IMPRVMT parameter.
        /**
         \return The \c OPPORTUNISTIC_MIN_F_IMPRVMT parameter.
         */
        const NOMAD::Double & get_opportunistic_min_f_imprvmt ( void ) const;
        
        /// Access to the \c OPPORTUNISTIC_LUCKY_EVAL parameter.
        /**
         \return The \c OPPORTUNISTIC_LUCKY_EVAL parameter.
         */
        bool get_opportunistic_lucky_eval ( void ) const;
        
        /// Access to the \c BB_MAX_BLOCK_SIZE parameter
        /**
         \return the \c BB_MAX_BLOCK_SIZE parameter
         */
        int get_bb_max_block_size ( void ) const ;
        
        /// Set the \c BB_MAX_BLOCK_SIZE parameter.
        /**
         \param bb_block_size The \c BB_MAX_BLOCK_SIZE parameter -- \b IN.
         */
        void set_BB_MAX_BLOCK_SIZE ( int bb_block_size );
        
        /// Access to the \c EVAL_POINTS_AS_BLOCK parameter
        /**
         \return true if points are evaluated as a block (number of elements >=1)
         */
        bool eval_points_as_block() const {return _eval_points_as_block;}
        
        
        /// Set the \c OPPORTUNISTIC_EVAL parameter.
        /**
         \param opp_eval The \c OPPORTUNISTIC_EVAL parameter -- \b IN.
         */
        void set_OPPORTUNISTIC_EVAL ( bool opp_eval );
        
        /// Set the \c OPPORTUNISTIC_MIN_NB_SUCCESS parameter.
        /**
         \param opp_min_nb_succ The \c OPPORTUNISTIC_MIN_NB_SUCCESS parameter -- \b IN.
         */
        void set_OPPORTUNISTIC_MIN_NB_SUCCESS ( int opp_min_nb_succ );
        
        /// Set the \c OPPORTUNISTIC_MIN_EVAL parameter.
        /**
         \param opp_min_eval The \c OPPORTUNISTIC_MIN_EVAL parameter -- \b IN.
         */
        void set_OPPORTUNISTIC_MIN_EVAL ( int opp_min_eval );
        
        /// Set the \c OPPORTUNISTIC_MIN_F_IMPRVMT parameter.
        /**
         \param opp_min_f_imprvt The \c OPPORTUNISTIC_MIN_F_IMPRVMT parameter -- \b IN.
         */
        void set_OPPORTUNISTIC_MIN_F_IMPRVMT ( const NOMAD::Double & opp_min_f_imprvt );
        
        /// Set the \c OPPORTUNISTIC_LUCKY_EVAL parameter.
        /**
         \param opp_lucky_eval The \c OPPORTUNISTIC_LUCKY_EVAL parameter -- \b IN.
         */
        void set_OPPORTUNISTIC_LUCKY_EVAL ( bool opp_lucky_eval );
        
        // RobustMads parameters
        // ----------------------------------
    private:
        bool _robust_mads;                      ///< Run robust mads or regular mads
        NOMAD::Double _robust_mads_standard_dev_factor;    ///< Parameter for smoothing in robust mads
        
    public:
        /// Access to the \c _robust_mads parameter.
        /**
         \return The \c robust_mads parameter.
         */
        bool get_robust_mads ( void ) const
        {
            return _robust_mads;
        }
        
        /// Set the \c ROBUST_MADS parameter.
        /**
         \param robust_mads The \c ROBUST_MADS parameter -- \b IN.
         */
        void set_ROBUST_MADS ( bool robust_mads )
        {
            _robust_mads = robust_mads;
        }
        
        /// Access to the \c _robust_mads_variance parameter.
        /**
         \return The \c robust_mads_variance parameter.
         */
        NOMAD::Double get_robust_mads_standard_dev_factor ( void ) const
        {
            return _robust_mads_standard_dev_factor;
        }
        
        /// Set the \c ROBUST_MADS_STANDARD_DEV_FACTOR parameter.
        /**
         \param beta The \c ROBUST_MADS_STANDARD_DEV_FACTOR parameter -- \b IN.
         */
        void set_ROBUST_MADS_STANDARD_DEV_FACTOR ( NOMAD::Double beta )
        {
            _robust_mads_standard_dev_factor = beta;
        }
        
        
        
        // SGTELIB parameters
        // ----------------------------------
    private:
        
        /// Number of evaluation of the model during each sgtelib_model-search step.
        int _sgtelib_model_eval_nb;
        int _sgtelib_model_candidates_nb;
        int _sgtelib_model_trials;
        
        std::string _sgtelib_model_definition;
        
        std::string _sgtelib_model_display;
        std::string _sgtelib_model_filter;
        
        /// Coefficient of the diversification term in the surrogate problem.
        NOMAD::Double _sgtelib_model_diversification;
        
        /// Coeff to avoid existing points
        NOMAD::Double _sgtelib_model_exclusion_area;
        
        /// Formulation of the surrogate problem in the sgtelib_model search.
        NOMAD::sgtelib_model_formulation_type _sgtelib_model_formulation;
        
        /// Method to compute the probability of feasibility with sgtelib_model.
        NOMAD::sgtelib_model_feasibility_type _sgtelib_model_feasibility;
        
        
    public:
        /// Set the \c SGTELIB_MODEL_EVAL_NB parameter.
        /**
         \param me The \c SGTELIB_MODEL_EVAL_NB parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_EVAL_NB ( const int me );
        
        /// Access to the \c SGTELIB_MODEL_EVAL_NB parameter.
        /**
         \return The \c SGTELIB_MODEL_EVAL_NB parameter.
         */
        int get_SGTELIB_MODEL_EVAL_NB ( void ) const;
        
        
        /// Set the \c SGTELIB_MODEL_TRIALS parameter.
        /**
         \param cn The \c SGTELIB_MODEL_TRIALS parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_TRIALS ( const int cn );
        
        /// Access to the \c SGTELIB_MODEL_TRIALS parameter.
        /**
         \return The \c SGTELIB_MODEL_TRIALS parameter.
         */
        int get_SGTELIB_MODEL_TRIALS ( void ) const;
        
        
        /// Set the \c SGTELIB_MODEL_CANDIDATES_NB parameter.
        /**
         \param me The \c SGTELIB_MODEL_CANDIDATES_NB parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_CANDIDATES_NB ( const int me );
        
        /// Access to the \c SGTELIB_MODEL_CANDIDATES_NB parameter.
        /**
         \return The \c SGTELIB_MODEL_CANDIDATES_NB parameter.
         */
        int get_SGTELIB_MODEL_CANDIDATES_NB ( void ) const;
        
        
        /// Set the \c SGTELIB_MODEL_DEFINITION parameter.
        /**
         \param s The \c SGTELIB_MODEL_DEFINITION parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_DEFINITION( const std::string & s );
        
        /// Access to the \c SGTELIB_MODEL_DEFINITION parameter.
        /**
         \return The \c SGTELIB_MODEL_DEFINITION parameter.
         */
        std::string get_SGTELIB_MODEL_DEFINITION ( void ) const;
        
        
        /// Set the \c SGTELIB_MODEL_DISPLAY parameter.
        /**
         \param s The \c SGTELIB_MODEL_DISPLAY parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_DISPLAY( const std::string & s );
        
        /// Access to the \c SGTELIB_MODEL_DISPLAY parameter.
        /**
         \return The \c SGTELIB_MODEL_DISPLAY parameter.
         */
        std::string get_SGTELIB_MODEL_DISPLAY ( void ) const;
        
        
        /// Set the \c SGTELIB_MODEL_FILTER parameter.
        /**
         \param s The \c SGTELIB_MODEL_FILTER parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_FILTER( const std::string & s );
        
        /// Access to the \c SGTELIB_MODEL_FILTER parameter.
        /**
         \return The \c SGTELIB_MODEL_FILTER parameter.
         */
        std::string get_SGTELIB_MODEL_FILTER ( void ) const;
        
        
        /// Set the \c SGTELIB_MODEL_DIVERSIFICATION parameter.
        /**
         \param dsw The \c SGTELIB_MODEL_DIVERSIFICATION parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_DIVERSIFICATION ( const NOMAD::Double dsw );
        
        /// Access to the \c SGTELIB_MODEL_DIVERSIFICATION parameter.
        /**
         \return The \c SGTELIB_MODEL_DIVERSIFICATION parameter.
         */
        NOMAD::Double get_SGTELIB_MODEL_DIVERSIFICATION ( void ) const;
        
        /// Set the \c SGTELIB_MODEL_EXCLUSION_AREA parameter.
        /**
         \param tc The \c SGTELIB_MODEL_EXCLUSION_AREA parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_EXCLUSION_AREA ( const NOMAD::Double tc );
        
        /// Access to the \c SGTELIB_MODEL_EXCLUSION_AREA parameter.
        /**
         \return The \c SGTELIB_MODEL_EXCLUSION_AREA parameter.
         */
        NOMAD::Double get_SGTELIB_MODEL_EXCLUSION_AREA ( void ) const;
        
        
        /// Set the \c SGTELIB_MODEL_FORMULATION parameter.
        /**
         \param dft The \c SGTELIB_MODEL_FORMULATION parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_FORMULATION ( const NOMAD::sgtelib_model_formulation_type dft );
        
        /// Access to the \c SGTELIB_MODEL_FORMULATION parameter.
        /**
         \return The \c SGTELIB_MODEL_FORMULATION parameter.
         */
        NOMAD::sgtelib_model_formulation_type get_SGTELIB_MODEL_FORMULATION (void) const;
        
        /// Set the \c SGTELIB_MODEL_FEASIBILITY parameter.
        /**
         \param dft The \c SGTELIB_MODEL_FEASIBILITY parameter -- \b IN.
         */
        void set_SGTELIB_MODEL_FEASIBILITY ( const NOMAD::sgtelib_model_feasibility_type dft  );
        
        /// Access to the \c SGTELIB_MODEL_FEASIBILITY parameter.
        /**
         \return The \c SGTELIB_MODEL_FEASIBILITY parameter.
         */
        NOMAD::sgtelib_model_feasibility_type get_SGTELIB_MODEL_FEASIBILITY (void) const;
        
        
        
    };
    
    /*----------------------------------------------------------------------*/
    
    /// Display a NOMAD::Parameters object.
    /**
     \param out The NOMAD::Display object -- \b IN.
     \param p   The NOMAD::Parameters object to be displayed -- \b IN.
     \return    The NOMAD::Display object.
     */
    inline const NOMAD::Display & operator << ( const NOMAD::Display    & out ,
                                               const NOMAD::Parameters & p     )
    {
        p.display ( out );
        return out;
    }
}

#endif
