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

// Holds SPIR information for each client query.
class SPIRClientQuery {
public:
    SPIRClientQuery ();
    ~SPIRClientQuery ();

    SPIRClientQuery (const SPIRClientQuery &other);
    SPIRClientQuery& operator= (SPIRClientQuery other);

    // initializes SPIR parameters
    bool init_parameters (const PercyClientParams &params, 
	    nservers_t num_servers, nservers_t t, const dbsize_t block_number,
	    const vec_ZZ_p &indices, const vec_vec_ZZ_p &shares, 
	    const vec_ZZ_pX &polyvec);

    // sends the infomation to the server
    bool send_to_server (const PercyClientParams &params, 
	    nservers_t server_index, std::ostream &os);

private:
    // Progress
    enum Progress {
	PROGRESS_NONE, 
	PROGRESS_INIT_PARAMS
    };
    Progress progress;

    vec_G1 commitvec;      // Vector of commitments to the polynomials
    ZZ_pX dot_f_a;         // Dot product of `f` vector and `a` vector
    G1 C_a;                // Commitment on dot_f_a
    vec_G1 witnessvec;     // Vector of witnesses to the evaluations

    G1 C_a_prime, witness_prime;
    GT D, D_prime;
    G1 eta_1, eta_2, eta_3;
    GT eta_4, eta_5, eta_6;
    ZZ_p v_0, v_1, v_2, v_3;
    vec_G1 muvec;
    vec_GT nuvec;
    vec_ZZ_p chivec;
    dbsize_t cvec_size;
    uint64_t *cvec;
};

#endif
