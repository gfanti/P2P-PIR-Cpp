//  Percy++ Copyright 2007,2012,2013 Ian Goldberg <iang@cs.uwaterloo.ca>,
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

#ifndef __RECOVER_H__
#define __RECOVER_H__

#include "rsdecoder.h"
#include "gf2e.h"

// The EasyRecover function based on Goldberg's Improving the Robustness of
// Private Information Retrival (2007).

NTL_CLIENT

bool EasyRecover(dbsize_t bytes_per_word, nservers_t t,
	nservers_t h, vector<DecoderResult<ZZ_p> > &results,
	dbsize_t word_number, const vec_ZZ_p &values, const vec_ZZ_p &indices);

bool EasyRecover(dbsize_t bytes_per_word, nservers_t t,
	nservers_t h, vector<DecoderResult<GF2E> > &results, 
	dbsize_t word_number, const vec_GF2E &values_vec, 
	const vec_GF2E &indices_vec);

#endif
