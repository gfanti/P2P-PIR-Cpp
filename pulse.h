//  Percy++ Copyright 2007,2012 Ian Goldberg <iang@cs.uwaterloo.ca>,
//  Casey Devet <cjdevet@cs.uwaterloo.ca>,
//  Paul Hendry <pshdenry@uwaterloo.ca>,
//  Ryan Henry <rhenry@cs.uwaterloo.ca>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of version 2 of the GNU General Public License as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  There is a copy of the GNU General Public License in the COPYING file
//  packaged with this plugin; if you cannot find it, write to the Free
//  Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//  02110-1301 USA

#ifndef __PULSE_H__
#define __PULSE_H__

// #include <stdint.h>
#include "gf2e.h"

// How many entries needed for the ratio test
const GF216_Element NUM_RATIOS = 3;

// How many bins each entry gets mapped to
const GF216_Element DEGREE = 3;


/* Use the same representation as AES */

extern const GF216_Element GF216_pulse_mtx_4bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_5bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_6bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_7bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_8bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_9bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_10bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_11bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_12bins[10000][3];
extern const GF216_Element GF216_pulse_mtx_1000bins[10000][3];

/* The interpolation matrices for V^-1 */
extern const GF216_Element GF216_V_inv_2servers[2];
extern const GF216_Element GF216_V_inv_3servers[3];
extern const GF216_Element GF216_V_inv_4servers[4];
extern const GF216_Element GF216_V_inv_5servers[5];

#endif
