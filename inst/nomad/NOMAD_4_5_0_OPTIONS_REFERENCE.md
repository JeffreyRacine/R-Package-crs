# NOMAD 4.5.0 Option Reference for `crs`

Generated on: 2026-02-23

Source of truth: vendored NOMAD definition files in `src/nomad4_src/src/Attribute`

Scope:
- This reference is generated from embedded NOMAD 4.5.0 sources used by `crs`.
- `snomadr(opts=...)` passes option names/values through to NOMAD.
- `display.nomad.progress=FALSE` also sets `DISPLAY_DEGREE 0`, `DISPLAY_ALL_EVAL false`, and `DISPLAY_UNSUCCESSFUL false`.
- `nomad.opt` in the working directory is read after `opts` and can override earlier values.
- Option names use NOMAD 4.5.0 terminology directly (for example, `MIN_FRAME_SIZE`).

## Problem Parameters (pbAttributesDefinition.txt)

Total options in this section: 13

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `BB_INPUT_TYPE` | `NOMAD::BBInputTypeList` | `* R` | The variable blackbox input types |
| `DIMENSION` | `size_t` | `0` | Dimension of the optimization problem (required) |
| `FIXED_VARIABLE` | `NOMAD::Point` | `-` | Fix some variables to some specific values |
| `GRANULARITY` | `NOMAD::ArrayOfDouble` | `-` | The granularity of the variables |
| `INITIAL_FRAME_SIZE` | `NOMAD::ArrayOfDouble` | `-` | The initial frame size of MADS |
| `INITIAL_MESH_SIZE` | `NOMAD::ArrayOfDouble` | `-` | The initial mesh size of MADS |
| `LOWER_BOUND` | `NOMAD::ArrayOfDouble` | `-` | The optimization problem lower bounds for each variable |
| `MIN_FRAME_SIZE` | `NOMAD::ArrayOfDouble` | `-` | Termination criterion on minimal frame size of MADS |
| `MIN_MESH_SIZE` | `NOMAD::ArrayOfDouble` | `-` | Termination criterion on minimal mesh size of MADS |
| `UPPER_BOUND` | `NOMAD::ArrayOfDouble` | `-` | The optimization problem upper bounds for each variable |
| `VARIABLE_GROUP` | `NOMAD::ListOfVariableGroup` | `-` | The groups of variables) |
| `X0` | `NOMAD::ArrayOfPoint` | `-` | The initial point(s) |
| `POINT_FORMAT` | `NOMAD::ArrayOfDouble` | `-` | Format of the doubles for trial points |

## Run Parameters (Core) (runAttributesDefinition.txt)

Total options in this section: 35

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `ADD_SEED_TO_FILE_NAMES` | `bool` | `true` | The flag to add seed to the file names |
| `ANISOTROPIC_MESH` | `bool` | `true` | MADS uses anisotropic mesh for generating directions |
| `ANISOTROPY_FACTOR` | `NOMAD::Double` | `0.1` | MADS anisotropy factor for mesh size change |
| `SEARCH_METHOD_MESH_PROJECTION` | `bool` | `true` | Projection on mesh for trial points from a search method |
| `DIRECTION_TYPE` | `NOMAD::DirectionTypeList` | `ORTHO N+1 QUAD` | Direction types for Poll steps |
| `DIRECTION_TYPE_SECONDARY_POLL` | `NOMAD::DirectionTypeList` | `DOUBLE` | Direction types for Mads secondary poll |
| `TRIAL_POINT_MAX_ADD_UP` | `size_t` | `0` | Max number of trial points |
| `ORTHO_MESH_REFINE_FREQ` | `size_t` | `1` | Control mesh refinement frequency |
| `FRAME_CENTER_USE_CACHE` | `bool` | `false` | Find best points in the cache and use them as frame centers |
| `H_MAX_0` | `NOMAD::Double` | `NOMAD::INF` | Initial value of hMax. |
| `HOT_RESTART_FILE` | `std::string` | `hotrestart.txt` | The name of the hot restart file |
| `HOT_RESTART_ON_USER_INTERRUPT` | `bool` | `false` | Flag to perform a hot restart on user interrupt |
| `HOT_RESTART_READ_FILES` | `bool` | `false` | Flag to read hot restart files |
| `HOT_RESTART_WRITE_FILES` | `bool` | `false` | Flag to write hot restart files |
| `MAX_ITERATIONS` | `size_t` | `INF` | The maximum number of iterations of the MADS algorithm |
| `MAX_ITERATION_PER_MEGAITERATION` | `size_t` | `INF` | Maximum number of Iterations to generate for each MegaIteration. |
| `MAX_TIME` | `size_t` | `INF` | Maximum wall-clock time in seconds |
| `MEGA_SEARCH_POLL` | `bool` | `false` | Evaluate points generated from Search and Poll steps all at once |
| `REJECT_UNKNOWN_PARAMETERS` | `bool` | `true` | Flag to reject unknown parameters when checking validity of parameters |
| `RHO` | `NOMAD::Double` | `0.1` | Rho parameter of the progressive barrier |
| `SEED` | `int` | `0` | The seed for the pseudo-random number generator |
| `RNG_ALT_SEEDING` | `bool` | `false` | With this option the seed is used to set xdef |
| `SIMPLE_LINE_SEARCH` | `bool` | `false` | MADS simple line search method complement speculative search |
| `SPECULATIVE_SEARCH` | `bool` | `true` | MADS speculative search method |
| `SPECULATIVE_SEARCH_BASE_FACTOR` | `NOMAD::Double` | `4.0` | Distance of the MADS speculative search method |
| `SPECULATIVE_SEARCH_MAX` | `size_t` | `1` | MADS speculative search method |
| `USER_SEARCH` | `bool` | `false` | MADS user search method provided as callback function |
| `STOP_IF_FEASIBLE` | `bool` | `false` | Stop algorithm once a feasible point is obtained |
| `STOP_IF_PHASE_ONE_SOLUTION` | `bool` | `false` | Stop algorithm once a phase one solution is obtained |
| `USER_CALLS_ENABLED` | `bool` | `true` | Controls the automatic calls to user function |
| `RANDOM_ALGO_SEARCH` | `bool` | `false` | A random search step for Mads using an algo (several iterations) |
| `RANDOM_ALGO_OPTIMIZATION` | `bool` | `false` | A standalone random optimization algo (several iterations) |
| `RANDOM_ALGO_DUMMY_FACTOR` | `size_t` | `1` | Dummy factor for random algo (used as template) |
| `H_NORM` | `NOMAD::HNormType` | `L2` | Norm type for infeasibility measure (h) computation |
| `H_MIN` | `NOMAD::Double` | `0` | Min h value for detecting infeasibility |

## Run Parameters (COOP) (runAttributesDefinitionCOOP.txt)

Total options in this section: 3

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `COOP_MADS_OPTIMIZATION` | `bool` | `false` | COOP-MADS optimization algorithm |
| `COOP_MADS_NB_PROBLEM` | `size_t` | `4` | Number of COOP-MADS problems |
| `COOP_MADS_OPTIMIZATION_CACHE_SEARCH` | `bool` | `true` | COOP-MADS cache search for incumbent synchronization |

## Run Parameters (Coordinate Search) (runAttributesDefinitionCS.txt)

Total options in this section: 1

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `CS_OPTIMIZATION` | `bool` | `false` | Coordinate Search optimization |

## Run Parameters (DMultiMads) (runAttributesDefinitionDMulti.txt)

Total options in this section: 8

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `DMULTIMADS_EXPANSIONINT_LINESEARCH` | `bool` | `false` | DMultiMads Expansion integer linesearch |
| `DMULTIMADS_QUAD_MODEL_STRATEGY` | `NOMAD::DMultiMadsQuadSearchType` | `MULTI` | Quad Model search strategies for DMultiMads |
| `DMULTIMADS_MIDDLEPOINT_SEARCH` | `bool` | `false` | DMultiMads Middle Point search |
| `DMULTIMADS_MIDDLEPOINT_SEARCH_CACHE_MAX` | `size_t` | `50` | DMultiMads middle point search |
| `DMULTIMADS_NM_STRATEGY` | `NOMAD::DMultiMadsNMSearchType` | `DOM` | Nelder-Mead search strategies for DMultiMads |
| `DMULTIMADS_OPTIMIZATION` | `bool` | `false` | DMultiMads solves multiobjective optimization problems |
| `DMULTIMADS_SELECT_INCUMBENT_THRESHOLD` | `size_t` | `1` | Control the choice of the DMultiMads incumbent |
| `DMULTIMADS_QMS_PRIOR_COMBINE_OBJ` | `bool` | `true` | Select compute method for objective of DMultiMads quad model search |

## Run Parameters (DiscoMads) (runAttributesDefinitionDisco.txt)

Total options in this section: 8

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `DISCO_MADS_OPTIMIZATION` | `bool` | `false` | DiscoMads optimization |
| `DISCO_MADS_DETECTION_RADIUS` | `NOMAD::Double` | `1.0` | Radius used to reveal discontinuities in DiscoMads |
| `DISCO_MADS_LIMIT_RATE` | `NOMAD::Double` | `1` | Limit rate of change used to reveal discontinuities in DiscoMads |
| `DISCO_MADS_EXCLUSION_RADIUS` | `NOMAD::Double` | `1` | Radius of exclusion balls around revealing points in DiscoMads |
| `DISCO_MADS_REVEALING_POLL_RADIUS` | `NOMAD::Double` | `2.02` | Revealing poll radius in DiscoMads |
| `DISCO_MADS_REVEALING_POLL_NB_POINTS` | `size_t` | `1` | Number of random points sampled by the revealing poll in DiscoMads |
| `DISCO_MADS_HID_CONST` | `bool` | `false` | Use DiscoMADS to reveal and escape hidden constraints regions |
| `DISCO_MADS_HID_CONST_OUTPUT_VALUE` | `NOMAD::Double` | `1E20` | Value attributed to objective function and PB constraints for failed evaluations |

## Run Parameters (IBEX) (runAttributesDefinitionIBEX.txt)

Total options in this section: 4

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `USE_IBEX` | `bool` | `false` | Boolean to determine if we want to use the functionnalities of IBEX |
| `SYSTEM_FILE_NAME` | `string` | `-` | File with the constraints |
| `SET_FILE` | `bool` | `false` | Boolean to determine if the file of the set is already created |
| `SET_FILE_NAME` | `string` | `-` | File to load with the set |

## Run Parameters (Latin Hypercube) (runAttributesDefinitionLH.txt)

Total options in this section: 2

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `LH_EVAL` | `size_t` | `0` | Latin Hypercube Sampling of points (no optimization) |
| `LH_SEARCH` | `NOMAD::LHSearchType` | `-` | Latin Hypercube Sampling Search method |

## Run Parameters (Nelder-Mead) (runAttributesDefinitionNM.txt)

Total options in this section: 11

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `NM_OPTIMIZATION` | `bool` | `false` | Nelder Mead stand alone optimization for constrained and unconstrained pbs |
| `NM_SEARCH` | `bool` | `true` | Nelder Mead optimization used as a search step for Mads |
| `NM_SIMPLEX_INCLUDE_LENGTH` | `NOMAD::Double` | `INF` | Construct NM simplex using points in cache. |
| `NM_SIMPLEX_INCLUDE_FACTOR` | `size_t` | `8` | Construct NM simplex using points in cache. |
| `NM_DELTA_E` | `NOMAD::Double` | `2` | NM expansion parameter delta_e. |
| `NM_DELTA_IC` | `NOMAD::Double` | `-0.5` | NM inside contraction parameter delta_ic. |
| `NM_DELTA_OC` | `NOMAD::Double` | `0.5` | NM outside contraction parameter delta_oc. |
| `NM_GAMMA` | `NOMAD::Double` | `0.5` | NM shrink parameter gamma. |
| `NM_SEARCH_MAX_TRIAL_PTS_NFACTOR` | `size_t` | `80` | NM-Mads search stopping criterion. |
| `NM_SEARCH_RANK_EPS` | `NOMAD::Double` | `0.01` | NM-Mads epsilon for the rank of DZ. |
| `NM_SEARCH_STOP_ON_SUCCESS` | `bool` | `false` | NM-Mads search stops on success. |

## Run Parameters (PSD-Mads) (runAttributesDefinitionPSDSSD.txt)

Total options in this section: 6

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `PSD_MADS_OPTIMIZATION` | `bool` | `0` | PSD-MADS optimization algorithm |
| `PSD_MADS_NB_VAR_IN_SUBPROBLEM` | `size_t` | `2` | Number of variables in PSD-MADS subproblems |
| `PSD_MADS_NB_SUBPROBLEM` | `size_t` | `INF` | Number of PSD-MADS subproblems |
| `PSD_MADS_ITER_OPPORTUNISTIC` | `bool` | `true` | Opportunistic strategy between the Mads subproblems in PSD-MADS |
| `PSD_MADS_ORIGINAL` | `bool` | `false` | Use NOMAD 3 strategy for mesh update in PSD-MADS |
| `PSD_MADS_SUBPROBLEM_PERCENT_COVERAGE` | `NOMAD::Double` | `70` | Percentage of variables that must be covered in subproblems before updating mesh |

## Run Parameters (QP Solver) (runAttributesDefinitionQPSolver.txt)

Total options in this section: 20

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `QP_OPTIMIZATION` | `bool` | `false` | Quad model stand alone QP optimization for constrained and unconstrained pbs |
| `QP_SEARCH` | `bool` | `false` | A quad model based search step for Mads using a QP solver |
| `QP_SelectAlgo` | `size_t` | `0` | Select the algorithm for QP solver |
| `QP_maxIter` | `size_t` | `20` | QPSolver outter loop iteration limit |
| `QP_tolDistDX` | `NOMAD::Double` | `-1.0` | A quad model based search step for Mads using a QP solver |
| `QP_absoluteTol` | `NOMAD::Double` | `1e-3` | A quad model based search step for Mads using a QP solver |
| `QP_tolCond` | `NOMAD::Double` | `1e-15` | A quad model based search step for Mads using a QP solver |
| `QP_tolMesh` | `NOMAD::Double` | `1.0` | A quad model based search step for Mads using a QP solver |
| `QP_relativeTol` | `NOMAD::Double` | `1e-3` | A quad model based search step for Mads using a QP solver |
| `QP_verbose` | `bool` | `false` | A quad model based search step for Mads using a QP solver |
| `QP_verboseFull` | `bool` | `false` | A quad model based search step for Mads using a QP solver |
| `QP_AugLag_mu0` | `NOMAD::Double` | `0.5` | A quad model based search step for Mads using a QP solver |
| `QP_AugLag_muDecrease` | `NOMAD::Double` | `2` | A quad model based search step for Mads using a QP solver |
| `QP_AugLag_eta0` | `NOMAD::Double` | `1.0` | A quad model based search step for Mads using a QP solver |
| `QP_AugLag_omega0` | `NOMAD::Double` | `1.0` | A quad model based search step for Mads using a QP solver |
| `QP_AugLag_maxIterInner` | `size_t` | `50` | QPSolver inner iteration limit for the subproblem |
| `QP_AugLag_tolDistDXInner` | `NOMAD::Double` | `1e-15` | A quad model based search step for Mads using a QP solver |
| `QP_AugLag_maxSuccessivFail` | `size_t` | `3` | A quad model based search step for Mads using a QP solver |
| `QP_AugLag_successRatio` | `NOMAD::Double` | `0.99` | A quad model based search step for Mads using a QP solver |
| `QP_SEARCH_MODEL_BOX_SIZE_LIMIT` | `NOMAD::Double` | `0` | QP solver generates trial points if bounds box size is above limit |

## Run Parameters (Quadratic Model) (runAttributesDefinitionQuadModel.txt)

Total options in this section: 8

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `QUAD_MODEL_SEARCH` | `bool` | `true` | Quad model search |
| `QUAD_MODEL_SEARCH_SIMPLE_MADS` | `bool` | `false` | Quad model search using a simpler version of Mads |
| `QUAD_MODEL_SEARCH_BOUND_REDUCTION_FACTOR` | `NOMAD::Double` | `1` | Scale the bounds for the quad model search |
| `QUAD_MODEL_DISPLAY` | `std::string` | `-` | Display of a model |
| `QUAD_MODEL_OPTIMIZATION` | `bool` | `false` | Quad model stand alone optimization for constrained and unconstrained pbs |
| `QUAD_MODEL_SEARCH_BOX_FACTOR` | `NOMAD::Double` | `4.0` | Quadratic model search point selection factor |
| `QUAD_MODEL_BOX_FACTOR` | `NOMAD::Double` | `4.0` | Quadratic model points selection box factor |
| `QUAD_MODEL_SEARCH_FORCE_EB` | `bool` | `false` | Quadratic model search optimization using extreme barrier |

## Run Parameters (Sgtelib Model) (runAttributesDefinitionSgtelibModel.txt)

Total options in this section: 14

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `SGTELIB_MODEL_EVAL` | `bool` | `0` | Sgtelib Model Sampling of points |
| `SGTELIB_MODEL_SEARCH` | `bool` | `false` | Model search using Sgtelib |
| `SGTELIB_MODEL_DISPLAY` | `std::string` | `-` | Display of a model |
| `SGTELIB_MODEL_DEFINITION` | `NOMAD::ArrayOfString` | `-` | Definition of the Sgtelib model |
| `SGTELIB_MODEL_SEARCH_TRIALS` | `size_t` | `1` | Max number of sgtelib model search failures before going to the poll step |
| `SGTELIB_MODEL_FORMULATION` | `NOMAD::SgtelibModelFormulationType` | `FS` | Formulation of the sgtelib model problem |
| `SGTELIB_MODEL_FEASIBILITY` | `NOMAD::SgtelibModelFeasibilityType` | `C` | Method used to model the feasibility of a point |
| `SGTELIB_MODEL_DIVERSIFICATION` | `NOMAD::Double` | `0.01` | Coefficient of the exploration term in the sgtelib model problem |
| `SGTELIB_MODEL_SEARCH_EXCLUSION_AREA` | `NOMAD::Double` | `0.0` | Exclusion area for the sgtelib model search around points of the cache |
| `SGTELIB_MODEL_SEARCH_CANDIDATES_NB` | `int` | `-1` | Number of candidates returned by the sgtelib model search |
| `SGTELIB_MIN_POINTS_FOR_MODEL` | `size_t` | `1` | Minimum number of valid points necessary to build a model |
| `SGTELIB_MAX_POINTS_FOR_MODEL` | `size_t` | `500` | Maximum number of valid points used to build a model |
| `SGTELIB_MODEL_SEARCH_FILTER` | `std::string` | `2345` | Methods used in the sgtelib search filter to return several search candidates |
| `SGTELIB_MODEL_RADIUS_FACTOR` | `NOMAD::Double` | `2.0` | Sgtelib model radius factor |

## Run Parameters (VNS) (runAttributesDefinitionVNS.txt)

Total options in this section: 7

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `VNS_MADS_OPTIMIZATION` | `bool` | `false` | VNS MADS stand alone optimization for constrained and unconstrained pbs |
| `VNSMART_MADS_SEARCH` | `bool` | `false` | VNS Mads search under condition of consecutive fails |
| `VNSMART_MADS_SEARCH_THRESHOLD` | `int` | `3` | Threshold for VNS (SMART) Mads search |
| `VNS_MADS_SEARCH` | `bool` | `false` | VNS Mads optimization used as a search step for Mads |
| `VNS_MADS_SEARCH_TRIGGER` | `NOMAD::Double` | `0.75` | VNS Mads search trigger |
| `VNS_MADS_SEARCH_WITH_SURROGATE` | `bool` | `false` | VNS Mads search with surrogate |
| `VNS_MADS_SEARCH_MAX_TRIAL_PTS_NFACTOR` | `size_t` | `100` | VNS-Mads search stopping criterion. |

## Evaluation Parameters (evalAttributesDefinition.txt)

Total options in this section: 5

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `BB_EVAL_FORMAT` | `NOMAD::ArrayOfDouble` | `-` | Format of the doubles sent to the blackbox evaluator |
| `BB_EXE` | `std::string` | `-` | Blackbox executable |
| `BB_REDIRECTION` | `bool` | `true` | Blackbox executable redirection for outputs |
| `BB_OUTPUT_TYPE` | `NOMAD::BBOutputTypeList` | `OBJ` | Type of outputs provided by the blackboxes |
| `SURROGATE_EXE` | `std::string` | `-` | Static surrogate executable |

## Display Parameters (displayAttributesDefinition.txt)

Total options in this section: 15

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `DISPLAY_ALL_EVAL` | `bool` | `false` | Flag to display all evaluations |
| `DISPLAY_DEGREE` | `int` | `2` | Level of verbose during execution |
| `DISPLAY_HEADER` | `size_t` | `40` | Frequency at which the stats header is displayed |
| `DISPLAY_INFEASIBLE` | `bool` | `true` | Flag to display infeasible |
| `DISPLAY_MAX_STEP_LEVEL` | `size_t` | `20` | Depth of the step after which info is not printed |
| `DISPLAY_STATS` | `NOMAD::ArrayOfString` | `BBE OBJ` | Format for displaying the evaluation points |
| `DISPLAY_FAILED` | `bool` | `false` | Flag to display failed evaluation |
| `DISPLAY_UNSUCCESSFUL` | `bool` | `false` | Flag to display unsuccessful |
| `STATS_FILE` | `NOMAD::ArrayOfString` | `-` | The name of the stats file |
| `EVAL_STATS_FILE` | `string` | `-` | The name of the file for stats about evaluations and successes |
| `SOL_FORMAT` | `NOMAD::ArrayOfDouble` | `-` | Internal parameter for format of the solution |
| `OBJ_WIDTH` | `size_t` | `0` | Internal parameter for character width of the objective |
| `HISTORY_FILE` | `std::string` | `-` | The name of the history file |
| `SOLUTION_FILE` | `std::string` | `-` | The name of the file containing the best feasible solution |
| `SOLUTION_FILE_FINAL` | `bool` | `false` | Flag to decide when to write best feasible solution |

## Cache Parameters (cacheAttributesDefinition.txt)

Total options in this section: 2

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `CACHE_FILE` | `std::string` | `""` | Cache file name |
| `CACHE_SIZE_MAX` | `size_t` | `INF` | Maximum number of evaluation points to be stored in the cache |

## Evaluator Control Parameters (evaluatorControlAttributesDefinition.txt)

Total options in this section: 6

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `EVAL_OPPORTUNISTIC` | `bool` | `true` | Opportunistic strategy: Terminate evaluations as soon as a success is found |
| `EVAL_SURROGATE_OPTIMIZATION` | `bool` | `false` | Use static surrogate as blackbox for optimization |
| `EVAL_USE_CACHE` | `bool` | `true` | Use cache in algorithms |
| `EVAL_QUEUE_SORT` | `NOMAD::EvalSortType` | `QUADRATIC_MODEL` | How to sort points before evaluation |
| `PSD_MADS_SUBPROBLEM_MAX_BB_EVAL` | `size_t` | `INF` | Max number of evaluations for each subproblem |
| `SUBPROBLEM_MAX_BB_EVAL` | `size_t` | `INF` | Internal parameter for PSD_MADS_SUBPROBLEM_MAX_BB_EVAL |

## Evaluator Control Global Parameters (evaluatorControlGlobalAttributesDefinition.txt)

Total options in this section: 17

| Option | Type | Default | Short Description |
| --- | --- | --- | --- |
| `BB_MAX_BLOCK_SIZE` | `size_t` | `1` | Size of blocks of points, to be used for parallel evaluations |
| `SURROGATE_MAX_BLOCK_SIZE` | `size_t` | `1` | Size of blocks of points, to be used for parallel evaluations |
| `EVAL_QUEUE_CLEAR` | `bool` | `true` | Opportunistic strategy: Flag to clear EvaluatorControl queue between each run |
| `EVAL_SURROGATE_COST` | `size_t` | `INF` | Cost of the surrogate function versus the true function |
| `MAX_BB_EVAL` | `size_t` | `INF` | Stopping criterion on the number of blackbox evaluations |
| `MAX_BLOCK_EVAL` | `size_t` | `INF` | Stopping criterion on the number of blocks evaluations |
| `MAX_EVAL` | `size_t` | `INF` | Stopping criterion on the number of evaluations (blackbox and cache) |
| `MAX_SURROGATE_EVAL_OPTIMIZATION` | `size_t` | `INF` | Stopping criterion on the number of static surrogate evaluations |
| `MODEL_MAX_BLOCK_SIZE` | `size_t` | `INF` | Internal parameter for QUAD_MODEL_MAX_BLOCK_SIZE and SGTELIB_MODEL_MAX_BLOCK_SIZE |
| `MODEL_MAX_EVAL` | `size_t` | `1000` | Internal parameter for QUAD_MODEL_MAX_EVAL and SGTELIB_MODEL_MAX_EVAL |
| `QUAD_MODEL_MAX_BLOCK_SIZE` | `size_t` | `INF` | Size of blocks of points, to be used for parallel evaluations |
| `QUAD_MODEL_MAX_EVAL` | `size_t` | `2000` | Max number of model evaluations for optimization of the quad model problem |
| `SGTELIB_MODEL_MAX_BLOCK_SIZE` | `size_t` | `INF` | Size of blocks of points, to be used for parallel evaluations |
| `SGTELIB_MODEL_MAX_EVAL` | `size_t` | `2000` | Max number of model evaluations for each optimization of the sgtelib model problem |
| `TMP_DIR` | `std::string` | `-` | Directory where to put temporary files |
| `USE_CACHE_FILE_FOR_RERUN` | `bool` | `false` | Cache file for rerun |
| `NB_THREADS_PARALLEL_EVAL` | `int` | `1` | Max number of threads used for parallel evaluations of each algorithm |

## Totals

- Total option entries across all definition files: 185
- To inspect option help at runtime: `snomadr(information = list("help" = "-h"))`
- To inspect a specific option at runtime: `snomadr(information = list("help" = "-h OPTION_NAME"))`
