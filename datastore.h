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

#ifndef __DATASTORE_H__
#define __DATASTORE_H__
#include <ZZ.h>
#include <set>
#include <map>
#include "percyparams.h"

NTL_CLIENT

class DataStore {
public:
    virtual ~DataStore() = 0;

    // Retrieve the cth word of the jth block of the database as a ZZ
    // Both c and j are 0-based.
    virtual ZZ get_word(dbsize_t c, dbsize_t j) = 0;
};

class MemoryDataStore : public DataStore {
protected:
    unsigned char *database;
    dbsize_t bytes_per_word, words_per_block, num_blocks;

public:
    MemoryDataStore (unsigned char * database, const PercyParams &params,
	    bool tau_independent = false);
    virtual ~MemoryDataStore ();
    
    // Retrieve the cth word of the jth block of the database as a ZZ
    // Both c and j are 0-based.
    virtual ZZ get_word(dbsize_t c, dbsize_t j);

    // Get the pointer to the data
    const unsigned char *get_data() const {return database;}

    // Get the size params
    dbsize_t get_bytes_per_word() { return bytes_per_word; }
    dbsize_t get_words_per_block() { return words_per_block; }
    dbsize_t get_num_blocks() { return num_blocks; }
};

class FileDataStore : public MemoryDataStore {
protected:
    int dbfd;
    unsigned char *mapptr;
    dbbits_t totbytes;
    dboffset_t offset;

public:
    // Create a DataStore backed by a file.
    FileDataStore(const char *filename, const PercyParams &params,
	    bool tau_independent = false);
    virtual ~FileDataStore();
};

#endif
