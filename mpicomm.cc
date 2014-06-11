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

#include <mpi.h>
#include "mpicomm.h"

// MPI Streaming implementation

MPIStreamBuf::MPIStreamBuf (int other_rank, MPI_Comm comm) :
    buffer(NULL),
    buffersize(0),
    other_rank(other_rank),
    comm(comm)
{
    MPI_Comm_rank(comm, &my_rank);
}

MPIStreamBuf::~MPIStreamBuf ()
{
    if (buffer != NULL) {
	delete[] buffer;
    }
}

int MPIStreamBuf::underflow ()
{
    if (gptr() < egptr()) {
	// The buffer is not all read
	return traits_type::to_int_type(*gptr());	
    }

    // We need to do another MPI receive
    // Get length of next receive.
    MPI_Status s1;
    MPI_Recv(&buffersize, 1, MPI_INT, other_rank, 0, comm, &s1);

    // Set up buffer
    if (buffer != NULL) {
	delete[] buffer;
    }
    buffer = new char[buffersize];

    // Do receive
    MPI_Status s2;
    MPI_Recv(buffer, buffersize, MPI_CHAR, other_rank, 0, comm, &s2);

#ifdef VERBOSE_MPI
    fprintf(stderr, "[MPI] %d <- %d - Received %d character(s)\n", my_rank, other_rank, buffersize);
#endif

    // Change pointers
    setg(buffer, buffer, buffer + buffersize);

    // Return current character
    return traits_type::to_int_type(*gptr());
}

int MPIStreamBuf::overflow (int c)
{
    if (c != EOF) {
	// Send the character
	char cc = traits_type::to_char_type(c);
	if (xsputn(&cc, 1) < 1) {
	    return EOF;
	}
    }
    return c;
}

std::streamsize MPIStreamBuf::xsputn (const char * s, std::streamsize n)
{
    // Send the size
    int size = (int)n;
    MPI_Send(&size, 1, MPI_INT, other_rank, 0, comm);

    // Send the string
    MPI_Send((char *)s, size, MPI_CHAR, other_rank, 0, comm);

#ifdef VERBOSE_MPI
    fprintf(stderr, "[MPI] %d -> %d - Sent %d character(s)\n", my_rank, other_rank, size);
#endif

    return size;
}

