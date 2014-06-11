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

#ifndef __THREADEDSERVER_H__
#define __THREADEDSERVER_H__

#include <vector>
#include <streambuf>
#include <iostream>
#include <stdio.h>
#include <NTL/vec_ZZ_p.h>
#include "datastore.h"
#include "percyparams.h"
#include "percyserver.h"
#include "percytypes.h"
#include "datastore.h"

class ThreadedDataStore : public FileDataStore {
public:
    ThreadedDataStore (const char * filename, const PercyServerParams &params,
	    bool tau_independent = false);
    virtual ~ThreadedDataStore ();

    MemoryDataStore * get_subds (dbsize_t index);
    dbsize_t get_num_threads ();

protected:
    dbsize_t num_threads;
    PercyThreadingType ttype;
    std::vector<MemoryDataStore*> subdatastores;
};

struct BufferInfo {
    char * addr;
    dbsize_t size;
};

class ThreadStreamBuf : public std::streambuf {
public:
    ThreadStreamBuf ();
    virtual ~ThreadStreamBuf ();

    void add_inbuffer (char * addr, dbsize_t size);
    void add_outbuffer (char * addr, dbsize_t size);

    virtual std::streamsize showmanyc ();
    virtual int underflow ();
    virtual int uflow ();
    virtual int overflow (int c = EOF);

    bool in_eof ();
    bool out_eof ();

private:
    std::vector<BufferInfo> inbuffers;
    std::vector<BufferInfo> outbuffers;
    std::streamsize inbufferindex;
    std::streamsize outbufferindex;
    dbsize_t incharindex;
    dbsize_t outcharindex;
};


class PercyThreadedServer : public PercyServer {
public:
    // Initialize a server with the given datastore
    PercyThreadedServer(ThreadedDataStore * ds);
    ~PercyThreadedServer();

    virtual bool handle_request(PercyServerParams &params, std::istream &is, std::ostream &os);

private:
    friend void * PThreadWork (void * arg);
    bool thread_work (dbsize_t thread_index, PercyServerParams &params, 
	    nqueries_t num_queries, unsigned char * const queries, unsigned char * const responses);

    dbsize_t num_threads;

    std::vector<PercyServer*> subservers;
};

#endif
