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

#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include "datastore.h"
#include "percyio.h"

#define FILENAME_SIZE 120

// You have to define even pure virtual destructors in C++.  You just do.
DataStore::~DataStore() {}


MemoryDataStore::MemoryDataStore (unsigned char * database, const PercyParams &params,
	bool tau_independent) :
    database(database)
{
    bytes_per_word = tau_independent ? params.modulus_bytes()
	: params.bytes_per_word();
    words_per_block = params.words_per_block();
    num_blocks = params.num_blocks();

    /*
    std::cerr << "bytes_per_word = " << bytes_per_word << "\n";
    std::cerr << "words_per_block = " << words_per_block << "\n";
    std::cerr << "num_blocks = " << num_blocks << "\n";
    */
}

MemoryDataStore::~MemoryDataStore () {}

// Retrieve the cth word of the jth block of the database as a ZZ.
// Both c and j are 0-based.
ZZ MemoryDataStore::get_word(dbsize_t c, dbsize_t j)
{
    ZZ ret;
    ZZFromBytes(ret, database +
	    (j*words_per_block+c)*(bytes_per_word), bytes_per_word);
    return ret;
}


FileDataStore::FileDataStore(const char *filename, const PercyParams &params,
	bool tau_independent) :
    MemoryDataStore(0, params)
{
    offset = 0;

    if (tau_independent) {
	// Process the PIRD header
	ifstream dbhead(filename);

	// We'll turn these into exceptions later
	char headbuf[6];
	dbhead.read(headbuf, 6);
	if (memcmp(headbuf, "PIRD\x01\x00", 6)) {
	    std::cerr << "Split database not in PIRD format.\n";
	    exit(1);
	}
	ZZ dbmodulus;
	percy_read_ZZ(dbhead, dbmodulus);

	if (!params.modulus_match(dbmodulus)) {
	    std::cerr << "Incorrect modulus for split database.\n";
	    exit(1);
	}

	offset = dbhead.tellg();
	dbhead.close();
    }

    // Open the file so we can mmap it
    dbfd = open(filename, O_RDONLY);

    if (dbfd < 0) {
	fprintf(stderr, "Could not open database %s\n", filename);
	exit(1);
    }
    totbytes = bytes_per_word * words_per_block * num_blocks;
    struct stat st;
    fstat(dbfd, &st);

    /*
    std::cerr << "database = " << filename << "\n";
    std::cerr << "st.st_size = " << st.st_size << "\n";
    std::cerr << "totbytes = " << totbytes << "\n";
    std::cerr << "offset = " << offset << "\n";
    */
    if (st.st_size < (off_t)totbytes+offset) {
        fprintf(stderr, "Database too small!\n");
        exit(1);
    }
    mapptr = (unsigned char *)MMAP(NULL, totbytes+offset, PROT_READ,
	    MAP_SHARED, dbfd, 0);
    if (mapptr == MAP_FAILED) {
	perror("mmap");
	exit(1);
    }
    database = mapptr + offset;
}

FileDataStore::~FileDataStore()
{
    munmap(mapptr, totbytes + offset);
    close(dbfd);
}

