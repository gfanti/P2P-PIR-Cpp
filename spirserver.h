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

#ifndef __SPIR_H__
#define __SPIR_H__

#include <vector>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_vec_ZZ_p.h>
#include "../PolyCommit/PolyCommitCommon.h"
#include "pspir_crypt.h"

#include "percyparams.h"

// This class is used to do the *server-side* work for the SPIR protocol.
class SPIRServerQuery {
public:
    SPIRServerQuery ();
    ~SPIRServerQuery ();

    SPIRServerQuery (const SPIRServerQuery &other);
    SPIRServerQuery& operator= (SPIRServerQuery other);

    // Read the required input from the client
    bool read_spir_input (PercyServerParams& params, std::istream &is);

    // Verify the query
    bool query_verification (PercyServerParams& params, vec_ZZ_p &inputvector);

    // Initialize the randomization
    bool init_randomization (PercyServerParams& params, unsigned char * prng_seed);

    // Add randomization to the reponse
    bool randomize_response (vec_ZZ_p &response);

private:
    // Progress
    enum Progress {
	PROGRESS_NONE, 
	PROGRESS_READ_INPUT, 
	PROGRESS_VERIFICATION, 
	PROGRESS_INIT_RANDOM
    };
    Progress progress;

    vec_G1 commitvec;
    G1 witness, witness_prime;
    uint64_t a_r;
    bool z;
    G1 C_a_prime;
    GT D, D_prime;
    G1 eta_1, eta_2, eta_3;
    GT eta_4, eta_5, eta_6;
    ZZ_p v_0, v_1, v_2, v_3;
    vec_G1 muvec;
    vec_GT nuvec;
    vec_ZZ_p chivec;
    dbsize_t cvec_size;
    uint64_t *cvec;
    vec_ZZ_p ephemeral_block;
    ZZ_p finalshare;
};

#endif
