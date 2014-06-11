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

//NTL includes
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ.h>
//PolyCommit includes
#include "../PolyCommit/PolyCommitCommon.h"

#include <gcrypt.h>


#include <iostream>
#include <sstream>

#include "pspir_crypt.h"

NTL_CLIENT

NTL_vector_impl(G1,vec_G1)
NTL_io_vector_impl(G1,vec_G1)
NTL_vector_impl(G2,vec_G2)
NTL_io_vector_impl(G2,vec_G2)
NTL_vector_impl(GT,vec_GT)
NTL_io_vector_impl(GT,vec_GT)


const int SHA256_LEN = gcry_md_get_algo_dlen(GCRY_MD_SHA256);

//hash(): hash a stream using SHA256
//  Note: Creates a buffer for the entire stream. For large, incrementally-created streams, 
//  you should use incr_hash instead and cumulatively hash to avoid creating a huge buffer.
unsigned char * hash(std::iostream &stream)
{
  stream.flush();
  unsigned char *obuf = new unsigned char[ SHA256_LEN + 1 ];
  obuf[SHA256_LEN] = '\0';
  stream.seekg(0, ios::end);
  const unsigned int length = (unsigned int) stream.tellg();
  stream.seekg(0, ios::beg);
  char *buffer = new char[length];
  stream.read(buffer, length);
  gcry_md_hash_buffer(GCRY_MD_SHA256, obuf, buffer, length);
  delete [] buffer;
  return obuf;
}

//hash_ZZ(): hash a stream using SHA256, and convert to ZZ
ZZ hash_ZZ(std::iostream &stream)
{
  return ZZFromBytes(hash(stream), SHA256_LEN);
}

//hash_G1(): hash a stream using SHA256, and convert to G1
G1 hash_G1(const PolyCommitParams &p, std::iostream &stream)
{
  const Pairing &e = p.get_pairing();
  unsigned char * hashv = hash(stream);
  return G1(e, hashv, SHA256_LEN);
}

//incr_hash_init(): Initialize a gcrypt hash handle for SHA256 hashing.
gcry_md_hd_t * incr_hash_init()
{
  gcry_md_hd_t * hash_handle_p = new gcry_md_hd_t;
  gcry_md_open(hash_handle_p, GCRY_MD_SHA256, 0);
  return hash_handle_p;
}

//incr_hash(): Incrementally hash data into an existing gcrypt hash handle.
//  Note: Maybe further efficiency gains could be made by re-using a buffer across all
//  incremental hashings, rather than allocate one each time.
void incr_hash(std::iostream &stream, gcry_md_hd_t *hash_handle_p)
{
  stream.flush();
  stream.seekg(0, ios::end);
  const unsigned int length = (unsigned int) stream.tellg();
  stream.seekg(0, ios::beg);
  char *buffer = new char[length];
  stream.read(buffer, length);
  gcry_md_write(*hash_handle_p, buffer, length);
  delete [] buffer;
}

//incr_hash_final(): Obtain the final hash from a gcrypt hash handle, and destroy it.
unsigned char * incr_hash_final(gcry_md_hd_t *hash_handle_p)
{
  unsigned char *obuf = new unsigned char[ SHA256_LEN + 1 ];
  obuf[SHA256_LEN] = '\0';
  const unsigned char *hash_val = gcry_md_read(*hash_handle_p, GCRY_MD_SHA256);
  //Copy the hash value, because as far as I can tell, the value returned by gcry_md_read 
  //is not valid after gcry_md_close() is called.
  strncpy((char *)obuf, (char *)hash_val, SHA256_LEN);
  gcry_md_close(*hash_handle_p);

  return obuf;
}

uint64_t hash_long(std::iostream &stream, const unsigned int bits)
{
  return trunc_long(hash_ZZ(stream), bits);
}

ZZ_p dot_product(const vec_ZZ_p &rhovec, const uint64_t *avec)
{
  const unsigned int blocks = rhovec.length();
  ZZ_p result;
  clear(result);
  for (uint64_t i = 0; i < blocks; i++)
  {
    result += rhovec[i] * avec[i];
  }
  return result;
}

ZZ_pX dot_product(const vec_ZZ_pX &polyvec, const uint64_t *avec)
{
  const unsigned int blocks = polyvec.length();
  ZZ_pX result;
  clear(result);
  for (uint64_t i = 0; i < blocks; i++)
  {
    result += polyvec[i] * avec[i];
  }
  return result;
}

