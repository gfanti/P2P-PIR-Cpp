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

#ifndef PSPIR_CRYPT_H
#define PSPIR_CRYPT_H

#include <iostream>
#include <gcrypt.h>
//NTL includes
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ.h>
//PolyCommit includes
#include "../PolyCommit/PolyCommitCommon.h"

//Parameter for soundness of ZKPs
#define KAPPA 40

NTL_CLIENT

//NTL_vector_decl(PolyCommitment,vec_PolyCommitment)
typedef vector<PolyCommitment> vec_PolyCommitment;
NTL_vector_decl(G1,vec_G1)
NTL_io_vector_decl(G1,vec_G1)
NTL_vector_decl(G2,vec_G2)
NTL_io_vector_decl(G2,vec_G2)
NTL_vector_decl(GT,vec_GT)
NTL_io_vector_decl(GT,vec_GT)


unsigned char * hash(std::iostream &stream);
ZZ hash_ZZ(std::iostream &stream);
G1 hash_G1(const PolyCommitParams &p, std::iostream &stream);

gcry_md_hd_t * incr_hash_init();
void incr_hash(std::iostream &stream, gcry_md_hd_t *hash_handle_p);
unsigned char * incr_hash_final(gcry_md_hd_t *hash_handle_p);
//G1 hash_G1(const PolyCommitParams &p, std::iostream &stream);
uint64_t hash_long(std::iostream &stream, const unsigned int bits);

/*ZZ_p LagrangeCoefficient(const vec_ZZ_p &indices, const ZZ_p &index, const ZZ_p &point = to_ZZ_p(0));

inline G1 BLSRecombine(const PolyCommitParams &p, const vec_G1 &sigmavec, const vec_ZZ_p &indices);
inline G2 BLSCreateGroupKey(const PolyCommitParams &p, const vec_G2 &pkvec, const vec_ZZ_p &indices);
inline G1 BLSSign(const PolyCommitParams &p, const G1 &message, const ZZ_p &sk);
inline bool BLSVerify(const PolyCommitParams &p, const G1 &message, const G2 &pk, const G1 &sigma);
*/

ZZ_p dot_product(const vec_ZZ_p &rhovec, const uint64_t *avec);

ZZ_pX dot_product(const vec_ZZ_pX &polyvec, const uint64_t *avec);

#endif /* PSPIR_CRYPT_H */
