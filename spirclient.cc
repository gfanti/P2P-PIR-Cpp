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
#include "spirclient.h"

SPIRClientQuery::SPIRClientQuery () :
    cvec_size(0),
    cvec(NULL)
{}

SPIRClientQuery::~SPIRClientQuery ()
{
    if (cvec != NULL) {
	delete[] cvec;
    }
}

SPIRClientQuery::SPIRClientQuery (const SPIRClientQuery &other) :
    progress(other.progress),
    commitvec(other.commitvec), dot_f_a(other.dot_f_a),
    C_a(other.C_a), witnessvec(other.witnessvec), C_a_prime(other.C_a_prime), 
    witness_prime(other.witness_prime), D(other.D), D_prime(other.D_prime),
    eta_1(other.eta_1), eta_2(other.eta_2), eta_3(other.eta_3),
    eta_4(other.eta_4), eta_5(other.eta_5), eta_6(other.eta_6), 
    v_0(other.v_0), v_1(other.v_1), v_2(other.v_2), v_3(other.v_3),
    muvec(other.muvec), nuvec(other.nuvec), chivec(other.chivec),
    cvec_size(other.cvec_size), cvec(NULL)
{
    if (cvec_size > 0) {
	cvec = new uint64_t[cvec_size];
    }
    for (dbsize_t i = 0; i < cvec_size; ++i) {
	cvec[i] = other.cvec[i];
    }
}

SPIRClientQuery& SPIRClientQuery::operator= (SPIRClientQuery other)
{
    progress = other.progress;
    commitvec = other.commitvec;
    dot_f_a = other.dot_f_a;
    C_a = other.C_a;
    witnessvec = other.witnessvec;
    C_a_prime = other.C_a_prime;
    witness_prime = other.witness_prime;
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
    return *this;
}

bool SPIRClientQuery::init_parameters (const PercyClientParams &params, 
	nservers_t num_servers, nservers_t t, const dbsize_t block_number, 
	const vec_ZZ_p &indices, const vec_vec_ZZ_p &shares, 
	const vec_ZZ_pX	&polyvec)
{
    dbsize_t num_blocks = params.num_blocks();

    //Variables used in query proof
    PolyCommitParams * pcparamsp = params.get_pcparamsp();
    const Pairing &e = pcparamsp->get_pairing();
    const G1 &g = pcparamsp->get_galphai(0);
    const G2 &ghat = pcparamsp->get_ghatalphai(0);
    const G2 &ghatalpha = pcparamsp->get_ghatalphai(1);
    const G1 &h = pcparamsp->get_h();
    const G1 &halpha = pcparamsp->get_halpha();
    const GT gt = e(g, ghat);
    const GT ht = e(h, ghat);
    const GPP<G1> h_pp(e, h);
    const GPP<G1> halpha_pp(e, halpha);
    const GPP<GT> gt_pp(e, gt);
    const GPP<GT> ht_pp(e, ht);
    const GPP<G2> ghat_pp(e, ghat);

    commitvec.SetLength(num_blocks);
    for (unsigned int i = 0; i < num_blocks; ++i) {
	// Compute commitment to the polynomial
	PolyCommitment c = PolyCommitment(pcparamsp, polyvec[i]);
	commitvec[i] = c.get_C_fast();
    }

    //Compute the avec hash
    uint64_t * avec = new uint64_t[num_blocks];
    stringbuf sb1;
    iostream hash_stream1(&sb1);
    hash_stream1 << commitvec[0];
    avec[0] = hash_long(hash_stream1, KAPPA);
    for (dbsize_t i = 1; i < num_blocks; ++i) {
	stringbuf sb2;
	iostream hash_stream2(&sb2);
	hash_stream2 << avec[i-1] << commitvec[i];
	avec[i] = hash_long(hash_stream2, KAPPA);
    }

    // Compute dot product of polyvec and avec
    dot_f_a = dot_product(polyvec, avec);

    // Compute a commitment on `dot_f_a`
    PolyCommitment pc_C_a = PolyCommitment(pcparamsp, dot_f_a);
    C_a = pc_C_a.get_C_fast();

    // Compute witnesses
    witnessvec.SetLength(num_servers);
    for (unsigned int j = 0; j < num_servers; ++j) {
	witnessvec[j] = pc_C_a.createwitness(indices[j]);

	ZZ_p rho_a = to_ZZ_p(0);
	for(unsigned int i=0; i<num_blocks; i++) {
	    rho_a += shares[i][j] * avec[i];
	}
    }


    /* Query proof phase */

    // 1. Choose gamma_0 and gamma_1, and compute the PolyCommit commitment (eta_1, eta_2 and eta_3)

    // Choose random gamma_0 and gamma_1
    const ZZ_p gamma_0 = random_ZZ_p();
    const ZZ_p gamma_1 = random_ZZ_p();
    const Zr gamma_0_Zr = to_Zr(e, gamma_0);
    const Zr gamma_1_Zr = to_Zr(e, gamma_1);

    // Compute the PolyCommit commitment and C_a_prime
    eta_1 = h_pp^gamma_0_Zr;
    eta_2 = h_pp^gamma_1_Zr;
    C_a_prime = C_a;
    eta_3 = C_a_prime^gamma_1_Zr;
    C_a_prime ^= gamma_0_Zr;

    // 2. Compute a Pedersen commitment (eta_4, eta_5, and eta_6) and witness (witness_prime)

    // Compute witness
    const PolyCommitment C_a_gamma_0(pcparamsp, dot_f_a * gamma_0);
    witness_prime = C_a_gamma_0.createwitness(to_ZZ_p(0));

    // Compute the Receipt (D) and the zero knowledge proof (D_prime)
    const ZZ_p gamma_2 = random_ZZ_p();
    const GT Y = e(C_a_prime, ghat) / e(witness_prime, ghatalpha);
    D = (gt_pp^to_Zr(e, avec[block_number])) * (ht_pp^to_Zr(e, gamma_2));
    D_prime = D^gamma_0_Zr;

    // Compute eta_4, eta_5 and eta_6
    const ZZ_p gamma_3 = random_ZZ_p();
    const ZZ_p gamma_4 = random_ZZ_p();
    const ZZ_p gamma_5 = random_ZZ_p();
    eta_4 = (gt_pp^to_Zr(e, gamma_3)) * (ht_pp^to_Zr(e, gamma_4));
    eta_5 = D^gamma_1_Zr;
    eta_6 = ht_pp^to_Zr(e, gamma_5);

    // Compute zeta
    stringbuf zeta_sb;
    iostream zeta_hash_stream(&zeta_sb);
    zeta_hash_stream << Y << eta_1 << eta_2 << eta_3 << eta_4 << eta_5 << eta_6; // << credential;
    const ZZ_p zeta = to_ZZ_p(hash_ZZ(zeta_hash_stream));

    // Compute v's
    v_0 = gamma_1 - gamma_0 * zeta;
    v_1 = avec[block_number] - gamma_3 * zeta;
    v_2 = gamma_2 - gamma_4 * zeta;
    v_3 = gamma_5 - gamma_0 * gamma_2 * zeta;

    chivec.SetLength(num_blocks);
    muvec.SetLength(num_blocks);
    nuvec.SetLength(num_blocks);
    // Choose random challenges
    cvec = new uint64_t[num_blocks];
    for (dbsize_t i = 0; i< num_blocks; i++) {
	cvec[i] = (i == block_number) ? 0 : RandomBits_long(KAPPA);
    }

    // Compute mu and nu for each record
    stringbuf c_sb;
    iostream c_hash_stream(&c_sb);
    for (dbsize_t j = 0; j < num_blocks; ++j)
    {
	chivec[j] = random_ZZ_p();
	muvec[j] = h_pp^to_Zr(e, (chivec[j] + to_ZZ_p(avec[j]) * cvec[j] * gamma_0));
	nuvec[j] = gt_pp^to_Zr(e, (chivec[j] + to_ZZ_p(avec[block_number]) * cvec[j] * gamma_0));

	c_hash_stream << muvec[j] << nuvec[j];
    }

    const uint64_t c = hash_long(c_hash_stream, KAPPA);
    //Compute missing challenge
    cvec[block_number] = c;
    for (dbsize_t i = 0; i < num_blocks; i++) {
	cvec[block_number] -= (i != block_number) * cvec[i];
	if (cvec[block_number] < 0) cvec[block_number] += (((uint64_t) 1) << KAPPA);
    } //
    chivec[block_number] -= to_ZZ_p(avec[block_number]) * cvec[block_number] * gamma_0;

    delete[] avec;

    progress = PROGRESS_INIT_PARAMS;
    return true;
}

bool SPIRClientQuery::send_to_server (const PercyClientParams &params,
	nservers_t server_index, std::ostream &os)
{
    if (progress < PROGRESS_INIT_PARAMS) {
	return -1;
    }

    dbsize_t num_blocks = params.num_blocks();
    unsigned long modulus_bytes = params.hybrid() ? params.modulussq_bytes() :
	    params.modulus_bytes();

    //Send the vector of polynomial commitments, muvec, and nuvec
    for (unsigned int i = 0; i < num_blocks; ++i) {
	os << commitvec[i] << muvec[i] << nuvec[i];
    }

    //Send the witness
    os << witnessvec[server_index];

    // Send eta_1,..., eta_6, C_a_prime, witness_prime, D, and D_prime
    os << eta_1 << eta_2 << eta_3 << eta_4 << eta_5 << eta_6
	<< C_a_prime << witness_prime << D << D_prime;

    // Send v_0,...,v_3
    unsigned char * vbytes = new unsigned char[4 * modulus_bytes];
    BytesFromZZ(vbytes, rep(v_0), modulus_bytes);
    BytesFromZZ(vbytes + modulus_bytes, rep(v_1), modulus_bytes);
    BytesFromZZ(vbytes + 2 * modulus_bytes, rep(v_2), modulus_bytes);
    BytesFromZZ(vbytes + 3 * modulus_bytes, rep(v_3), modulus_bytes);
    os.write((char *)vbytes, 4*modulus_bytes);
    delete[] vbytes;

    // Send chivec
    unsigned char * cbytes = new unsigned char[modulus_bytes];
    for (dbsize_t i = 0; i < num_blocks; ++i) {
	BytesFromZZ(cbytes, rep(chivec[i]), modulus_bytes);
	os.write((char *)cbytes, modulus_bytes);
    }
    delete[] cbytes;

    // Send cvec
    uint64_t * cvecarr = new uint64_t[num_blocks];
    for (dbsize_t i = 0; i < num_blocks; ++i) {
	cvecarr[i] = cvec[i];
    }
    os.write((char *)cvecarr, num_blocks * sizeof(uint64_t));
    delete[] cvecarr;

    return true;
}

