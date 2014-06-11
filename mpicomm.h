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

#ifndef __MPICOMM_H__
#define __MPICOMM_H__

#include <mpi.h>
#include <iostream>
#include <vector>
#include "percytypes.h"
#include "percyparams.h"
#include "percycomm.h"

class MPIStreamBuf : public std::streambuf
{
public:
    MPIStreamBuf (int other_rank, MPI_Comm comm = MPI_COMM_WORLD);
    ~MPIStreamBuf ();

protected:
    virtual int underflow ();
    virtual int overflow (int c = EOF);
    virtual std::streamsize xsputn (const char * s, std::streamsize n);

private:
    char * buffer;
    int buffersize;

    int other_rank;
    int my_rank;
    MPI_Comm comm;
};

#endif
