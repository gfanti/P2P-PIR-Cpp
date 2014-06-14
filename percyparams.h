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

#ifndef __PERCYPARAMS_H__
#define __PERCYPARAMS_H__

#include <iostream>
#include <NTL/ZZ_p.h>
#include <vector>
#ifdef SPIR_SUPPORT
#include <PolyCommitCommon.h>
#endif
#include "percytypes.h"

#define PERCY_VERSION 1

NTL_CLIENT

enum PercyMode { MODE_ZZ_P, MODE_GF28, MODE_GF216, MODE_CHOR, MODE_RS_SYNC}; //, MODE_PULSE_SYNC };

enum PercyThreadMethod {THREAD_METHOD_NONE, THREAD_METHOD_PTHREAD, THREAD_METHOD_FORK};

enum PercyThreadingType {THREADING_ROWS, THREADING_QUERIES};

class PercyParams {
    friend ostream& operator<<(ostream& os, const PercyParams &params);
    friend istream& operator>>(istream& is, PercyParams &params);

    //friend void mpi_write_params (int dest, const PercyParams& params, MPI_Comm comm);
    //friend void mpi_read_params (int src, PercyParams& params, MPI_Comm comm);

public:
    PercyParams();

    // Use the given modulus directly; only encryption will be possible
    PercyParams(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	    nservers_t tau, ZZ modulus, PercyMode mode,
	    char *pcparams_file = NULL, bool do_spir = false);

    unsigned long modulus_bytes() const {
	return NumBytes(modulus);
    }
    unsigned long modulussq_bytes() const {
	return NumBytes(modulus*modulus);
    }
    dbsize_t bytes_per_word() const {
	unsigned long r = NumBytes(modulus);
	if (r > 1) return r - 1;
	return 1;
    }
    dbsize_t words_per_block() const {
	return _words_per_block;
    }
    dbsize_t num_blocks() const {
	return _num_blocks;
    }
    dbsize_t max_unsynchronized() const {
	return _max_unsynchronized;
    }
    dbsize_t bytes_per_block() const {
	return bytes_per_word() * _words_per_block;
    }
    nservers_t tau() const {
	return _tau;
    }
    void mod_modulus() const {
	modctx.restore();
    }
    void mod_modulussq() const {
	modsqctx.restore();
    }
    const ZZ &get_modulus() const {
	return modulus;
    }
    const ZZ_p &get_g() const {
	return g;
    }
    bool hybrid() const {
	return hybrid_protection;
    }
	bool spir() const {
	return do_spir;
	}
    bool modulus_match(ZZ testmod) const {
	return modulus == testmod;
    }
    PercyMode get_mode() const {
        return mode;
    }
    char * get_pcparams_filename() const {
	return pcparams_filename;
    }	
#ifdef SPIR_SUPPORT
    PolyCommitParams * get_pcparamsp() const {
        return pcparamsp;
    }
#endif
    //bool is_gf28() const {
	//return modulus == 256;
    //}
    // Encrypt the given plaintext.  The current ZZ_p context must be
    // modsqctx.
    ZZ_p encrypt(ZZ plaintext) const {
	ZZ_p r;
	random(r);
	return power(g, plaintext) * power(r, modulus);
    }

protected:
    void create_ZZ_pContexts();
    unsigned char version;
    bool hybrid_protection;
    bool do_spir;
    nservers_t _tau;
    dbsize_t _words_per_block, _num_blocks, _max_unsynchronized;
    ZZ_pContext modctx, modsqctx;
    // Paillier public key
    ZZ modulus;
    ZZ_p g;   // mod modulus^2
    PercyMode mode;
    char * pcparams_filename;
#ifdef SPIR_SUPPORT
    PolyCommitParams * pcparamsp;
#endif
};

class PercyClientParams : public PercyParams {
public:
    // Use the given modulus directly; only encryption will be possible
    PercyClientParams(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	    nservers_t tau, ZZ modulus, PercyMode mode, char *pcparams_file=NULL, bool do_spir=false)
        : PercyParams(words_per_block, num_blocks, max_unsynchronized, tau, modulus, mode, 
		pcparams_file, do_spir) {}

    // Generate a new public/private key pair of the given keysize
    PercyClientParams(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	    nservers_t tau, unsigned long modulus_bits);

    // Use the given factors to generate a public/private key pair
    PercyClientParams(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	    nservers_t tau, ZZ p, ZZ q);// : PercyParams(words_per_block, num_blocks, tau, p, q) {}

    // Decrypt the given ciphertext.  This routine will change the
    // current ZZ_p context to modctx.
    ZZ_p decrypt(ZZ_p ciphertext) const {
	modsqctx.restore();
	ZZ Lval = rep(power(ciphertext, lambda) - 1) / modulus;
	modctx.restore();
	ZZ_p ret = to_ZZ_p(Lval) * mu;
	return ret;
    }		
    ZZ get_p1() const { return p1; }
    ZZ get_p2() const { return p2; }
protected:
    void init_hybrid(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	    nservers_t tau, ZZ p, ZZ q);
    // Paillier private key
    ZZ lambda;
    ZZ_p mu;  // mod modulus
    ZZ p1, p2;
};


struct serverinfo {
    nservers_t sid;
    char * addr;
    uint16_t port;
};

class PercyServerParams : public PercyParams {
public:
    // Use the given modulus directly; only encryption will be possible
    PercyServerParams(dbsize_t words_per_block, dbsize_t num_blocks, dbsize_t max_unsynchronized,
	    nservers_t tau, ZZ modulus, PercyMode mode, bool be_byzantine,
	    char *pcparams_filename, bool do_spir, nservers_t sid,
	    dbsize_t num_threads = 0, PercyThreadingType ttype = THREADING_ROWS,
	    PercyThreadMethod tmethod = THREAD_METHOD_PTHREAD) :
	PercyParams(words_per_block, num_blocks, max_unsynchronized, tau, modulus, mode, 
		pcparams_filename, do_spir), 
	sid(sid),
	be_byzantine(be_byzantine),
	num_threads(num_threads),
	ttype(ttype),
	tmethod(tmethod)
    {}

    bool is_compatible(PercyParams& other)
    {
	bool result = (this->get_mode() == other.get_mode())
		&& (this->words_per_block() == other.words_per_block())
		&& (this->num_blocks() == other.num_blocks())
        && (this->max_unsynchronized() == other.max_unsynchronized())
		&& (this->bytes_per_word() == other.bytes_per_word())
		&& (this->tau() == other.tau())
		&& (this->spir() == other.spir())
		&& (this->hybrid() == other.hybrid());

	if (result && this->hybrid())
	{
	    result &= this->modulus_match(other.get_modulus())
		    && this->get_g() == other.get_g();
	}

        return result;
    }
    nservers_t get_sid() const {
        return sid;
    }

    bool is_byzantine() const {
	return be_byzantine;
    }

    dbsize_t get_num_threads () const {
	return num_threads;
    }

    PercyThreadingType get_ttype () const {
	return ttype;
    }

    PercyThreadMethod get_tmethod () const {
	return tmethod;
    }

protected:
    nservers_t sid;
    bool be_byzantine;

    // Threading parameters
    dbsize_t num_threads;
    PercyThreadingType ttype;
    PercyThreadMethod tmethod;
};


#endif
