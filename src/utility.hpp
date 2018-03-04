/*
 * =====================================================================================
 *
 *       Filename:  utility.hpp
 *
 *    Description:  isnan and isinf.
 *
 *        Version:  1.0
 *        Created:  03/03/2018 18:49:31
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zhenghua Nie (ZHN), 
 *        Company:  McMaster University
 *        
 * =====================================================================================
 */

#ifndef __CRS_UTILITY_HPP__
#define __CRS_UTILITY_HPP__

inline bool crs_isnan ( double x ) { return x != x; }
inline bool crs_isinf (double x) { return !crs_isnan(x) && crs_isnan(x - x); }

#endif
