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

#include <sstream>
#include <algorithm>
#include "spirserver.h"

SPIRServerQuery::SPIRServerQuery () :
    progress(PROGRESS_NONE),
    z(true),
    cvec_size(0),
    cvec(NULL)
{}

SPIRServerQuery::~SPIRServerQuery ()
{
    if (cvec != NULL) {
	delete[] cvec;
    }
}

SPIRServerQuery::SPIRServerQuery (const SPIRServerQuery &other) :
    progress(other.progress),
    commitvec(other.commitvec), witness(other.witness),
    witness_prime(other.witness_prime), a_r(other.a_r), z(other.z),
    C_a_prime(other.C_a_prime), D(other.D), D_prime(other.D_prime),
    eta_1(other.eta_1), eta_2(other.eta_2), eta_3(other.eta_3),
    eta_4(other.eta_4), eta_5(other.eta_5), eta_6(other.eta_6), 
    v_0(other.v_0), v_1(other.v_1), v_2(other.v_2), v_3(other.v_3),
    muvec(other.muvec), nuvec(other.nuvec), chivec(other.chivec),
    cvec_size(other.cvec_size), cvec(NULL),
    ephemeral_block(other.ephemeral_block), finalshare(other.finalshare)
{
    if (cvec_size > 0) {
	cvec = new uint64_t[cvec_size];
    }
    for (dbsize_t i = 0; i < cvec_size; ++i) {
	cvec[i] = other.cvec[i];
    }
}

SPIRServerQuery& SPIRServerQuery::operator= (SPIRServerQuery other)
{
    progress = other.progress;
    commitvec = other.commitvec;
    witness = other.witness;
    witness_prime = other.witness_prime;
    a_r = other.a_r;
    z = other.z;
    C_a_prime = other.C_a_prime;
    D = other.D;
    D_prime = other.D_prime;
    eta_1 = other.eta_1;
    eta_2 = other.eta_2;
    eta_3 = other.eta_3;
    eta_4 = other.eta_4;
    eta_5 = other.eta_5;
    eta_5 = other.eta_6;
    v_0 = other.v_0;
    v_1 = other.v_1;
    v_2 = other.v_2;
    v_3 = other.v_3;
    muvec = other.muvec;
    nuvec = other.nuvec;
    chivec = other.chivec;
    cvec_size = other.cvec_size;
    std::swap(cvec, other.cvec);
    ephemeral_block = other.ephemeral_block;
    finalshare = other.finalshare;
    return *this;
}

bool SPIRServerQuery::read_spir_input (PercyServerParams &params, std::istream &is)
{
    PolyCommitParams * pcparamsp = params.get_pcparamsp();

    // Set the appropriate modulus
    unsigned int modulus_bytes;
    if (params.hybrid()) {
        params.mod_modulussq();
        modulus_bytes = params.modulussq_bytes();
    } else {
        params.mod_modulus();
        modulus_bytes = params.modulus_bytes();
    }

    // Read some values from the params
    dbsize_t num_blocks = params.num_blocks();

    const Pairing &e = pcparamsp->get_pairing();

    // Read in the vector of commitments
    commitvec.SetLength(num_blocks);
    muvec.SetLength(num_blocks);
    nuvec.SetLength(num_blocks);
    for(unsigned int i = 0; i < num_blocks; ++i) {
	commitvec[i] = G1(e);
	muvec[i] = G1(e);
	nuvec[i] = GT(e);
	is >> commitvec[i] >> muvec[i] >> nuvec[i];
	if (is.eof()) {
	    return false;
	}
    }

    // Read in the witness
    witness = G1(e);
    is >> witness;
    if (is.eof()) {
	return false;
    }

    // Read in eta_1,..., eta_6, C_a_prime, witness_prime, D, D_prime
    eta_1 = G1(e);
    eta_2 = G1(e);
    eta_3 = G1(e);
    eta_4 = GT(e);
    eta_5 = GT(e);
    eta_6 = GT(e);
    C_a_prime = G1(e);
    witness_prime = G1(e);
    D = GT(e);
    D_prime = GT(e);
    is >> eta_1 >> eta_2 >> eta_3 >> eta_4 >> eta_5 >> eta_6
	>> C_a_prime >> witness_prime >> D >> D_prime;
    if (is.eof()) {
	return false;
    }

    // Read in v_0,...,v_3
    unsigned char * vbytes = new unsigned char[4 * modulus_bytes];
    is.read((char *)vbytes, 4*modulus_bytes);
    if (is.eof()) {
	return false;
    }
    ZZ inputz[4];
    for (unsigned int k=0; k<4; k++) {
	ZZFromBytes(inputz[k], vbytes + k * modulus_bytes, modulus_bytes);
    }
    v_0 = to_ZZ_p(inputz[0]);
    v_1 = to_ZZ_p(inputz[1]);
    v_2 = to_ZZ_p(inputz[2]);
    v_3 = to_ZZ_p(inputz[3]);
    delete[] vbytes;

    // Read in chivec
    unsigned char * cbytes = new unsigned char[modulus_bytes];
    ZZ cinputz;
    chivec.SetLength(num_blocks);
    for (dbsize_t i = 0; i < num_blocks; ++i) {
	is.read((char *)cbytes, modulus_bytes);
	if (is.eof()) {
	    delete[] cbytes;
	    return false;
	    }
	ZZFromBytes(cinputz, cbytes, modulus_bytes);
	chivec[i] = to_ZZ_p(cinputz);
    }
    delete[] cbytes;

    // Read in cvec
    cvec = new uint64_t[num_blocks];
    cvec_size = num_blocks;
    uint64_t * cvecarr = new uint64_t[num_blocks];
    is.read((char *)cvecarr, num_blocks * sizeof(uint64_t));
    if (is.eof()) {
	delete[] cvecarr;
	return false;
    }
    for (dbsize_t i = 0; i < num_blocks; ++i) {
	cvec[i] = cvecarr[i];
    }
    delete[] cvecarr;

    progress = PROGRESS_READ_INPUT;
    return true;
}

bool SPIRServerQuery::query_verification (PercyServerParams &params, vec_ZZ_p &inputvector)
{
    if (progress < PROGRESS_READ_INPUT) {
	return false;
    }

    PolyCommitParams * pcparamsp = params.get_pcparamsp();

    // Read some values from the params
    dbsize_t num_blocks = params.num_blocks();

    const Pairing &e = pcparamsp->get_pairing();
    const G1 &g = pcparamsp->get_galphai(0);
    const G2 &ghat = pcparamsp->get_ghatalphai(0);
    const G2 &ghatalpha = pcparamsp->get_ghatalphai(1);
    const GT gt = e(g, ghat);
    const G1 &h = pcparamsp->get_h();
    const G1 &halpha = pcparamsp->get_halpha();
    const GT ht = e(h, ghat);
    const GPP<G1> h_pp(e, h);
    const GPP<G1> halpha_pp(e, halpha);
    const GPP<GT> gt_pp(e, gt);
    const GPP<GT> ht_pp(e, ht);
    const GPP<G2> ghat_pp(e, ghat);

    // Query verification phase

    //Compute avec
    uint64_t * avec = new uint64_t[num_blocks];
    stringbuf sb1;
    iostream hash_stream1(&sb1);
    hash_stream1 << commitvec[0];
    avec[0] = hash_long(hash_stream1, KAPPA);
    for (unsigned int i = 1; i < num_blocks; i++) {
	stringbuf sb2;
	iostream hash_stream2(&sb2);
	hash_stream2 << avec[i-1] << commitvec[i];
	avec[i] = hash_long(hash_stream2, KAPPA);
    }
    a_r = avec[num_blocks-1];

    //Compute C_a
    G1 C_a = G1(e, 1);
    uint64_t mask = ((uint64_t) 1)<<(KAPPA-1);
    for (unsigned int k = 0; k < KAPPA; ++k, mask >>= 1) {
	C_a = C_a.square();
	for (dbsize_t i = 0; i < num_blocks; i++) {
	    if ((avec[i] & mask)) {
		C_a *= commitvec[i];
	    }
	}
    }

    //Compute rho_a and verify the evaluation
    ZZ_p rho_a = dot_product(inputvector, avec);
    const Zr ir = to_Zr(e, params.get_sid());
    const Zr fir = to_Zr(e, rho_a);
    const GT LHS = e(C_a, ghat);
    const GT RHS = e(witness, ghatalpha/(ghat^ir)) * e(g, ghat^fir);
    if (!(LHS == RHS)) {
	std::cerr << "Error: Query verification phase failed." << std::endl;
	return false;
    }

    //Verify query proof phase, part I
    const GT Y = e(C_a_prime, ghat) / e(witness_prime, ghatalpha);

    stringbuf zeta_sb;
    iostream zeta_hash_stream(&zeta_sb);
    zeta_hash_stream << Y << eta_1 << eta_2 << eta_3 << eta_4 << eta_5 << eta_6; // << credential;
    const ZZ_p zeta = to_ZZ_p(hash_ZZ(zeta_hash_stream));

    const Zr zeta_Zr = to_Zr(e, zeta);
    const Zr v_0_Zr = to_Zr(e, v_0);

    z = (eta_2 == ((h_pp^v_0_Zr) * (eta_1^zeta_Zr)));
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part I, line 1) failed.\n";
	delete[] avec;
	return false;
    }
    z = (eta_3 == G1::pow2(e, C_a, v_0_Zr, C_a_prime, zeta_Zr));
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part I, line 2) failed.\n";
	delete[] avec;
	return false;
    }
    z = (D == ((gt_pp^to_Zr(e, v_1)) * (ht_pp^to_Zr(e, v_2)) * (eta_4^zeta_Zr)));
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part I, line 3) failed.\n";
	delete[] avec;
	return false;
    }
    z = (eta_5 == GT::pow2(e, D, v_0_Zr, D_prime, zeta_Zr));
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part I, line 4) failed.\n";
	delete[] avec;
	return false;
    }
    z = (eta_6 == ((ht_pp^to_Zr(e, v_3)) * ((D_prime/Y)^zeta_Zr)));
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part I, line 5) failed.\n";
	delete[] avec;
	return false;
    }

    //Verify query proof phase, part II

    uint64_t * bvec = new uint64_t[num_blocks];
    for(dbsize_t i = 0; i < num_blocks; i++) {
	bvec[i] = RandomBits_long(KAPPA);
    }

    G1 lhs1 = G1(e, 1);
    GT lhs2 = GT(e, 1);
    uint64_t lhs3 = 0;

    ZZ_p rhs_exp1 = to_ZZ_p(0);
    ZZ_p rhs1_exp2 = to_ZZ_p(0);
    ZZ_p rhs2_exp2 = to_ZZ_p(0);

    mask = ((uint64_t) 1)<<(KAPPA-1);
    for (unsigned int k = 0; k < KAPPA; ++k, mask >>= 1) {
	lhs1 = lhs1.square();
	lhs2 = lhs2.square();
	for (dbsize_t j = 0; j < num_blocks; ++j)
	{
	    if ((bvec[j] & mask)) {
		lhs1 *= muvec[j];
		lhs2 *= nuvec[j];
	    }
	}
    }

    stringbuf c_sb;
    iostream c_hash_stream(&c_sb);
    for (dbsize_t j = 0; j < num_blocks; ++j)
    {
	rhs_exp1 += chivec[j] * bvec[j];
	rhs1_exp2 += to_ZZ_p(avec[j]) * bvec[j] * cvec[j];
	rhs2_exp2 += to_ZZ_p(bvec[j]) * cvec[j];

	lhs3 += cvec[j];
	c_hash_stream << muvec[j] << nuvec[j];
    }
    const uint64_t c = hash_long(c_hash_stream, KAPPA);
    Zr rhs_exp1_Zr = to_Zr(e, rhs_exp1);

    delete[] avec;
    delete[] bvec;

    z = (lhs1 == ((h_pp^rhs_exp1_Zr) * (eta_1^to_Zr(e, rhs1_exp2))));
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part II, line 1) failed.\n";
	return false;
    }

    z = (lhs2 == ((gt_pp^rhs_exp1_Zr) * (Y^to_Zr(e, rhs2_exp2))));
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part II, line 2) failed.\n";
	return false;
    }

    lhs3 = (lhs3 << (64-KAPPA)) >> (64-KAPPA);
    z = (lhs3 == c);
    if (!z)
    {
	cerr << "Error: verification of query proof phase (part II, line 3) failed.\n";
	return false;
    }

    progress = PROGRESS_VERIFICATION;
    fprintf(stderr, "Query verified\n");
    return true;
}

bool SPIRServerQuery::init_randomization (PercyServerParams &params, unsigned char * prng_seed)
{
    if (progress < PROGRESS_VERIFICATION) {
	return false;
    }

    // Read some values from the params
    dbsize_t words_per_block = params.words_per_block();

    //Set the NTL PRNG seed
    for(unsigned int i=0; i<8; i++) {
	prng_seed[16+i] = 0xff & (a_r >> (8*i));
    }
    NTL::SetSeed(ZZFromBytes(prng_seed, 24));

    //Compute the ephemeral (r+1)th database record for the query
    ephemeral_block.SetLength(words_per_block);
    for (unsigned int c = 0; c < words_per_block; c++) {
	ephemeral_block[c] = random_ZZ_p();
    }

    //Compute the final share, a pseudorandom polynomial g where g(0)=0
    ZZ_pX randpoly = random_ZZ_pX(params.tau());
    SetCoeff(randpoly,0,0);
    finalshare = eval(randpoly, to_ZZ_p(params.get_sid()));

    progress = PROGRESS_INIT_RANDOM;
    return true;
}

bool SPIRServerQuery::randomize_response (vec_ZZ_p &response)
{
    if (progress < PROGRESS_INIT_RANDOM) {
	return false;
    }

    // Add randomness to the response vector
    response += finalshare * ephemeral_block;

    return true;
}

