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

/* This file defines the algorithm that is chosen for the portfolio algorithm.
 * To get the result for (k, t, h), use portfolio_table[k-3][t+1][h-t-1].
 */

#ifndef __PORTFOLIO_H__
#define __PORTFOLIO_H__

#include <string>
#include <map>

// This file is used by the portfolio algorithm in RSDecoder to choose the best
// single-polynomial decoding algorithm to use.  The results are based off of
// timing tests, the results of which are stored in the portfolio_table array.
// The portfolio algorithm only works for k<=20 (as these are the only timing
// results that are stored here).

// A type for the different RS decoding algorithms.  (Maybe TestType should be
// changed to something more suitable).
enum TestType {
    UNDEFINED = 0,
    UNKNOWN,
    EASY,
    BEST,
    BRUTE,
    KOTTER,
    CH_MS,
    BW, 
    CH_MULTI,
    CH_TK1,
    DP, 
    MAX_TESTTYPE
};
const std::string testTypeStrings[MAX_TESTTYPE] = { 
    "undefined",
    "unknown",
    "easy",
    "best",
    "brute",
    "kotter",
    "ch_ms",
    "bw",
    "ch_multi",
    "ch_tk1",
    "dp"
};

// A type for the flavour of the dynamic programming decoding algorithm
enum DPType {
    UNDEFINED_DPTYPE = 0,
    ASSUME_CORRECT,
    ASSUME_WRONG,
    ASSUME_SHARES,
    MAX_DPTYPE
};
const std::string DPTypeStrings[MAX_DPTYPE] = { 
    "undefined",
    "assume_correct",
    "assume_wrong",
    "assume_shares"
};

// Defaults for findpolys
//static TestType defaultTestType = KOTTER;
//static DPType defaultDPType = UNDEFINED_DPTYPE;
//static int defaultGord = 1;


struct Choice {
    TestType testType;
    double time;
    DPType dpType;
    int gord;
};

const int max_k = 25;

extern const Choice portfolio_table[max_k - 2][max_k + 1][max_k + 1];

// Function used by RSDecoder to get the best algorithm to run.
bool portfolioChoice (int k, int t, int h,
        TestType &testType, DPType &dpType, int &gord);

#endif

