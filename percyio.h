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

#ifndef __PERCYIO_H__
#define __PERCYIO_H__

#include <ZZ.h>
#include <stdint.h>

NTL_CLIENT

// Read and write little-endian XX-bit "uintXX_t"s and ZZs.
void percy_write_le_uint64(ostream &is, uint64_t val);
void percy_write_le_uint32(ostream &os, uint32_t val);
void percy_write_le_uint16(ostream &os, uint16_t val);
void percy_write_ZZ(ostream &os, ZZ val);

void percy_read_le_uint64(istream &is, uint64_t &val);
void percy_read_le_uint32(istream &is, uint32_t &val);
void percy_read_le_uint16(istream &is, uint16_t &val);
void percy_read_ZZ(istream &is, ZZ &val);

#endif
