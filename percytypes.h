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

#ifndef __PERCYTYPES_H__
#define __PERCYTYPES_H__

//  Percy++
//  Copyright 2007 Ian Goldberg <iang@cs.uwaterloo.ca>
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
//  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
//  02111-1307  USA

#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <stddef.h>	// size_t
#include <sys/types.h>  // off_t

// This header contains typedefs for seamless switching between 32-
// and 64-bit builds of Percy++.
//
// The macros that we use to detect the build type are adapted from
// the discussion at: "http://stackoverflow.com/questions/5272825/detecting-64bit-compile-in-c"

// LEGEND
// sid_t      : server ids
// dbsize_t   : blocks; words/block; bytes/word
// dbbits_t   : total bits in database
// dboffset_t : offset into database (in bytes)
// nqueries_t : number of simultanous queries
// nservers_t : number of servers; privacy threshold; tau-independence

//#if UINTPTR_MAX == 0xffffffff /* 32-bit */
#if defined(__i386__) /* IA-32 */

//	#warning "Detected 32-bit build."

	#define BUILD_32
	typedef uint32_t sid_t;
	typedef uint32_t dbsize_t;
	typedef uint64_t dbbits_t;
	typedef off_t dboffset_t;
	typedef uint16_t nqueries_t;
	typedef uint16_t nservers_t;

	#define MMAP(u,v,w,x,y,z) mmap(u,v,w,x,y,z)
	#define PERCY_READ_LE_DBSIZE(a,b) percy_read_le_uint32(a,b)
	#define PERCY_WRITE_LE_DBSIZE(a,b) percy_write_le_uint32(a,b)

//#elif UINTPTR_MAX == 0xffffffffffffffff /* 64-bit */
#elif defined(__x86_64__)|defined(__IA64__) /* AMD64|IA-64 */
//#elif defined(__LP64__) /* (long int,pointer)->64-bit; int->32-bit */

//	#warning "Detected 64-bit build."

	#define BUILD_64
	#define _W160_OPT

	typedef uint32_t sid_t;
	typedef uint64_t dbsize_t;
	typedef uint64_t dbbits_t;
	typedef off64_t dboffset_t;
	typedef uint16_t nqueries_t;
	typedef uint16_t nservers_t;

	#define MMAP(u,v,w,x,y,z) mmap64(u,v,w,x,y,z)
	#define PERCY_READ_LE_DBSIZE(a,b) percy_read_le_uint64(a,b)
	#define PERCY_WRITE_LE_DBSIZE(a,b) percy_write_le_uint64(a,b)

#else

	#error "Unsupported architecture."

#endif

#endif
