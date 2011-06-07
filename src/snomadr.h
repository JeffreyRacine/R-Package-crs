/*
 * =====================================================================================
 *
 *       Filename:  nomadr.hpp
 *
 *    Description:  Header files for nomadr 
 *
 *        Version:  1.0
 *        Created:  15/05/2011 13:51:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zhenghua Nie (ZHN), zhenghua.nie@gmail.com
 *        Company:  McMaster University
 *
 *    Copyright (C) 2011 Zhenghua Nie. All Rights Reserved.
 *    This code is published under GNU GENERAL PUBLIC LICENSE.
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License,  or
 *    (at your option) any later version.
 *      
 *    This program is distributed WITHOUT ANY WARRANTY. See the
 *    GNU General Public License for more details.
 *           
 *    If you do not have a copy of the GNU General Public License,  
 *    write to the Free Software Foundation, Inc., 
 *    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *              
 *        
 * =====================================================================================
 */

#ifndef __NOMADR_HPP__
#define __NOMADR_HPP__

#include <iomanip> 
#include <iostream> 
#include <fstream> 
#include <stdio.h> 
#include <string.h> 
#include <string> 
#include <assert.h> 
#include <sys/types.h>     
#include <sys/stat.h>    
#include <math.h> 
#include <cmath> 

#include "nomad.hpp"

#include <R.h>
#include <Rdefines.h>
#include <Rinterface.h>
#include <R_ext/Utils.h>

typedef union UNION_VARIABLE
{
int ivalue;
double dvalue;
} union_variable;

typedef struct ONE_VAR
{
int var_type;
union_variable l_x;
union_variable u_x;
union_variable x0;
} one_var;


const char *csbbin[]={"continuous", "integer", "categorical", "binary"};
const char *csbbout[]={"obj", "pb", "eb", "peb_p", "peb_e", "filter", "cnt_eval", "stat_avg", "stat_sum", "undefined_bbo"};
const char *stop_message[]={
      "NO_STOP                   ,  ///No stop", 
      "ERROR                     ,  ///< Error", 
      "UNKNOWN_STOP_REASON        ,  ///< Unknown", 
      "CTRL_C                     ,  ///< Ctrl-C", 
      "USER_STOPPED               ,  ///< User-stopped in Evaluator::update_iteration()", 
      "MESH_PREC_REACHED          ,  ///< Mesh minimum precision reached", 
      "X0_FAIL                    ,  ///< Problem with starting point evaluation", 
      "P1_FAIL                    ,  ///< Problem with phase one", 
      "DELTA_M_MIN_REACHED        ,  ///< Min mesh size", 
      "DELTA_P_MIN_REACHED        ,  ///< Min poll size", 
      "L_MAX_REACHED              ,  ///< Max mesh index", 
      "L_LIMITS_REACHED           ,  ///< Mesh index limits", 
      "MAX_TIME_REACHED           ,  ///< Max time", 
      "MAX_BB_EVAL_REACHED        ,  ///< Max number of blackbox evaluations", 
      "MAX_SGTE_EVAL_REACHED      ,  ///< Max number of surrogate evaluations", 
      "MAX_EVAL_REACHED           ,  ///< Max number of evaluations", 
      "MAX_SIM_BB_EVAL_REACHED    ,  ///< Max number of sim bb evaluations", 
      "MAX_ITER_REACHED           ,  ///< Max number of iterations", 
      "FEAS_REACHED               ,  ///< Feasibility", 
      "F_TARGET_REACHED           ,  ///< F_TARGET", 
      "STAT_SUM_TARGET_REACHED    ,  ///< STAT_SUM_TARGET", 
      "L_CURVE_TARGET_REACHED     ,  ///< L_CURVE_TARGET", 
      "MULTI_MAX_BB_REACHED       ,  ///< Max number of blackbox evaluations (multi obj.)", 
      "MULTI_NB_MADS_RUNS_REACHED ,  ///< Max number of MADS runs (multi obj.)", 
      "MULTI_STAGNATION           ,  ///< Stagnation criterion (multi obj.)", 
      "MULTI_NO_PARETO_PTS        ,  ///< No Pareto points (multi obj.)", 
      "MAX_CACHE_MEMORY_REACHED   ,    ///< Max cache memory"
    };

#endif
