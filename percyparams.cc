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

#include <string.h>
#include <fstream>
#include "percyparams.h"
#include "percyio.h"

#ifdef SPIR_SUPPORT
// Open the specified file, read in the PolyCommit parameters contained 
// within, create a new Params object from the parameters, and return it
PolyCommitParams * pcparams_init(const char * filename) {
    if (!filename) {
        return NULL;
    }
    ifstream ifile(filename);
    if(!ifile.is_open()) {
        std::cerr << "Error: Cannot open params file." << endl;
        exit(1);
    }
    PolyCommitParams * pcparamsp = new PolyCommitParams();
    ifile >> *pcparamsp;
    ifile.close();

    return pcparamsp;
}
#endif

PercyParams::PercyParams()
{
    this->version = PERCY_VERSION;
    this->hybrid_protection = false;
    this->_tau = 0;
    this->do_spir = false;
    this->_words_per_block = 1;
    this->_num_blocks = 1;
    this->_max_unsynchronized = 0;
    this->_expansion_factor = 1;
    this->modulus = 257;
    this->mode = MODE_ZZ_P;
    this->pcparams_filename = NULL;
#ifdef SPIR_SUPPORT
    this->pcparamsp = NULL;
#endif
    create_ZZ_pContexts();
}

PercyParams::PercyParams(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	dbsize_t expansion_factor, nservers_t tau, ZZ modulus, PercyMode mode,
	char *pcparams_file, bool do_spir)
{
    this->version = PERCY_VERSION;
    this->hybrid_protection = false;
    this->_tau = tau;
    this->do_spir = do_spir;
    this->_words_per_block = words_per_block;
    this->_num_blocks = num_blocks;
    this->_max_unsynchronized = max_unsynchronized;
    this->_expansion_factor = expansion_factor;
    this->modulus = modulus;
    this->mode = mode;
    this->pcparams_filename = pcparams_file;
#ifdef SPIR_SUPPORT
    this->pcparamsp = pcparams_init(pcparams_file);
#endif
    create_ZZ_pContexts();
}

ostream& operator<<(ostream& os, const PercyParams &params)
{
    unsigned char minibuf[4];

    // Output the magic header
	os.write("PIRC", 4);
	//printf("SENT: PIRC\n");

    // Output the version number and flags
    minibuf[0] = params.version;
    minibuf[1] = (params.hybrid_protection ? 1 : 0) | (params._tau ? 2 : 0);
    os.write((char *)minibuf, 2);

    // Output the SPIR flag
    minibuf[0] = params.do_spir;
    os.write((char *)minibuf, 1);

    // Output the mode
    minibuf[0] = params.mode;
    os.write((char *)minibuf, 1);

    // Output the words_per_block and num_blocks values
    PERCY_WRITE_LE_DBSIZE(os, params._words_per_block);
    PERCY_WRITE_LE_DBSIZE(os, params._num_blocks);

    // Output the max_unsynchronized value
    PERCY_WRITE_LE_DBSIZE(os, params._max_unsynchronized);
    PERCY_WRITE_LE_DBSIZE(os, params._expansion_factor);
    
    // Output the modulus
    percy_write_ZZ(os, params.modulus);

    // Output g, if appropriate
    if (params.hybrid_protection) {
        percy_write_ZZ(os, rep(params.g));
    }

    return os;
}

istream& operator>>(istream& is, PercyParams &params)
{
    unsigned char minibuf[4];

    // Input the magic header
	is.read((char *)minibuf, 4);
	//printf("RECEIVED: %u %u %u %u\n", minibuf[0], minibuf[1], minibuf[2], minibuf[3]);
	
    if (memcmp(minibuf, "PIRC", 4)) {
	std::cerr << "Did not find expected PercyParams header.\n";
	return is;
    }

    // Input the version number and flags
    is.read((char *)minibuf, 2);
    params.version = minibuf[0];
    params.hybrid_protection = ( (minibuf[1] & 1) == 1 );
    params._tau = ( (minibuf[1] & 2) == 2 );
    if (params.version != PERCY_VERSION) {
	std::cerr << "Did not find expected PercyParams version number " <<
	    PERCY_VERSION << ".\n";
	return is;
    }

    // Input the SPIR flag
    is.read((char *)minibuf, 1);
    params.do_spir = (PercyMode) minibuf[0];

    // Input the mode
    is.read((char *)minibuf, 1);
    params.mode = (PercyMode) minibuf[0];

    // Input the words_per_block and num_blocks values
    PERCY_READ_LE_DBSIZE(is, params._words_per_block);
    PERCY_READ_LE_DBSIZE(is, params._num_blocks);

    // Input the max_unsynchronized value
    PERCY_READ_LE_DBSIZE(is, params._max_unsynchronized);

    // Input the expansion factor value
    PERCY_READ_LE_DBSIZE(is, params._expansion_factor);
    
    // Input the modulus
    percy_read_ZZ(is, params.modulus);
    params.create_ZZ_pContexts();

    // Input g, if appropriate
    if (params.hybrid_protection) {
	ZZ_pContext savectx;
	savectx.save();
	params.modsqctx.restore();
	ZZ gz;
	percy_read_ZZ(is, gz);
	params.g = to_ZZ_p(gz);
	savectx.restore();
    }

    return is;
}

PercyClientParams::PercyClientParams(dbsize_t words_per_block, dbsize_t num_blocks, 
        dbsize_t max_unsynchronized, dbsize_t expansion_factor, nservers_t tau, ZZ p, ZZ q)
{
    init_hybrid(words_per_block, num_blocks, max_unsynchronized, expansion_factor, tau, p, q);
}

void PercyClientParams::init_hybrid(dbsize_t words_per_block, dbsize_t num_blocks, 
        dbsize_t max_unsynchronized, dbsize_t expansion_factor, nservers_t tau, ZZ p, ZZ q)
{
    this->version = PERCY_VERSION;
    this->hybrid_protection = true;
    this->_tau = tau;
    this->_words_per_block = words_per_block;
    this->_num_blocks = num_blocks;
    this->_max_unsynchronized = max_unsynchronized;
    this->_expansion_factor = expansion_factor;
    ZZ modulus;
    modulus = p * q;
    this->modulus = modulus;
    this->mode = MODE_ZZ_P;
    this->p1 = p;
    this->p2 = q;
#ifdef SPIR_SUPPORT
    this->pcparamsp = NULL;
#endif
    create_ZZ_pContexts();

    // Generate the Paillier public and private parts
    ZZ pm1, qm1;
    pm1 = p - 1;
    qm1 = q - 1;
    this->lambda = pm1 * qm1 / GCD(pm1, qm1);

    ZZ_pContext savectx;
    savectx.save();
    mod_modulussq();
    random(this->g);
    ZZ muinv = rep(power(this->g, this->lambda) - 1) / modulus;
    mod_modulus();
    this->mu = inv(to_ZZ_p(muinv));
    savectx.restore();
}

void PercyParams::create_ZZ_pContexts()
{
    // Create the ZZ_pContexts
    ZZ_pContext modctx(modulus);
    ZZ_pContext modsqctx(modulus * modulus);
    this->modctx = modctx;
    this->modsqctx = modsqctx;
}

PercyClientParams::PercyClientParams(dbsize_t words_per_block, dbsize_t num_blocks, 
        dbsize_t max_unsynchronized, dbsize_t expansion_factor, nservers_t tau, unsigned long modulus_bits)
{
    // Pick the sizes for the primes
    unsigned long qsize = modulus_bits / 2;
    unsigned long psize = modulus_bits - qsize;

    // Generate random primes of the appropriate size.  We ensure the
    // top two bits are set so that their product is of the right
    // bitlength.
    ZZ pbase, qbase, p, q;
    RandomBits(pbase, psize);
    if (psize >= 1) SetBit(pbase, psize-1);
    if (psize >= 2) SetBit(pbase, psize-2);
    NextPrime(p, pbase);
    RandomBits(qbase, qsize);
    if (qsize >= 1) SetBit(qbase, qsize-1);
    if (qsize >= 2) SetBit(qbase, qsize-2);
    NextPrime(q, qbase);

    init_hybrid(words_per_block, num_blocks, max_unsynchronized, expansion_factor, tau, p, q);
}
