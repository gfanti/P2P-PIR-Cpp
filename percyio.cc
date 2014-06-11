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

#include "percyio.h"

// Read and write little-endian XX-bit "uintXX_t"s and ZZs.
void percy_write_le_uint64(ostream &os, uint64_t val)
{
    unsigned char buf[8];

    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;
    buf[2] = (val >> 16) & 0xff;
    buf[3] = (val >> 24) & 0xff;
	buf[4] = (val >> 32) & 0xff;
	buf[5] = (val >> 40) & 0xff;
	buf[6] = (val >> 48) & 0xff;
	buf[7] = (val >> 56) & 0xff;

    os.write((char *)buf, 8);
}

void percy_write_le_uint32(ostream &os, uint32_t val)
{
    unsigned char buf[4];

    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;
    buf[2] = (val >> 16) & 0xff;
    buf[3] = (val >> 24) & 0xff;

    os.write((char *)buf, 4);
}

void percy_write_le_uint16(ostream &os, uint16_t val)
{
    unsigned char buf[2];

    buf[0] = val & 0xff;
    buf[1] = (val >> 8) & 0xff;

    os.write((char *)buf, 2);
}

void percy_write_ZZ(ostream &os, ZZ val)
{ 
    // Output the length
    long len = NumBytes(val);
    if (len > 65535) {
	std::cerr << "ZZ is too long (" << len << " bytes)\n";
	len = 65535;
    }
    unsigned char *buf = new unsigned char[len];
    BytesFromZZ(buf, val, len);
    percy_write_le_uint16(os, len);
    os.write((char *)buf, len);
    delete[] buf;
}

void percy_read_le_uint64(istream &is, uint64_t &val)
{
    unsigned char buf[8];

    is.read((char *)buf, 8);

    val = buf[0] | (buf[1]<<8) | (buf[2]<<16) | (buf[3]<<24)
	     | ((uint64_t)buf[4]<<32) | ((uint64_t)buf[5]<<40)
		 | ((uint64_t)buf[6]<<48) | ((uint64_t)buf[7]<<56);
}

void percy_read_le_uint32(istream &is, uint32_t &val)
{
    unsigned char buf[4];

    is.read((char *)buf, 4);

    val = buf[0] | (buf[1] << 8) | (buf[2] << 16) | (buf[3] << 24);
}

void percy_read_le_uint16(istream &is, uint16_t &val)
{
    unsigned char buf[2];

    is.read((char *)buf, 2);

    val = buf[0] | (buf[1] << 8);
}

void percy_read_ZZ(istream &is, ZZ &val)
{
    // Input the length
    unsigned short len;
    percy_read_le_uint16(is, len);
    unsigned char *buf = new unsigned char[len];
    is.read((char *)buf, len);
    ZZFromBytes(val, buf, len);
    delete[] buf;
}
