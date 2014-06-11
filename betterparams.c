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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* For a given n,h,l, we must find integers k,m with k<m such that the
   following is negative. */
static long latticecond(long n, long h, long l, long m, long k)
{
    return l*m*(m-1)/2 + n*(k+1)*k/2 - h*k*m;
}

static long oldm(long n, long h, long l)
{
    return (2*h*n+h*h-l*n-1)/(h*h-l*n); /* implicit floor */
}

static long oldk(long n, long h, long l, long m)
{
    return h*m/n - 1;  /* implicit floor */
}

static long goodm(long n, long h, long l)
{
    return ((h-l)*n)/(h*h-l*n)+1;  /* implicit floor */
}

static long goodk(long n, long h, long l, long m)
{
    return (h*m)/n;  /* implicit floor */
}

static long long cost(long n, long h, long l, long m, long k)
{
    return (long long)m*m*m*n*k*(n-h)*k;
}

int main(int argc, char **argv)
{
    long n,h,l,m,k,v,om,ok,ov;
    long long cst, ocst;

    long startn = 3;

    if (argc > 1) {
	startn = atol(argv[1]);
    }

    for (n=startn;;++n) {
	for (l=1;l<=n;++l) {
	    for (h=(long)(floor(sqrt((double)n*l)))+1;h<=n;++h) {
		m = goodm(n,h,l);
		k = goodk(n,h,l,m);
		v = latticecond(n,h,l,m,k);
		cst = cost(n,h,l,m,k);
		om = oldm(n,h,l);
		ok = oldk(n,h,l,om);
		ov = latticecond(n,h,l,om,ok);
		ocst = cost(n,h,l,om,ok);
		if (v >= 0) {
		    exit(1);
		}
		while (v < 0) {
		    printf("%ld %ld %ld: %ld %ld %ld [%lld] / %ld %ld %ld [%lld]\n",
			    n,h,l,m,k,v,cst,om,ok,ov,ocst);
		    --k;
		    v = latticecond(n,h,l,m,k);
		}
	    }
	}
    }
}
