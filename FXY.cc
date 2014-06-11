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

#include "FXY.h"

NTL_OPEN_NNS

ZZ comb(long n, long k)
{
  //
  // compute binomial coefficients
  //
  // comb(n,k) = n! / k!(n-k)!
  //

  if( k==0 ) return to_ZZ(1);
  if( k==1 ) return to_ZZ(n);

  ZZ ONE = to_ZZ(1);
  ZZ N = to_ZZ(n);
  ZZ K = to_ZZ(k);
  ZZ c = to_ZZ(1);
  if( k < n-k ){
    for(ZZ i = N; i >= N-K+ONE; i--){
      c = ((c*i) / (N-i+ONE));
    }
  }else{
    for(ZZ i = N; i >= K+ONE; i--){
      c = ((c*i) / (N-i+ONE));
    }
  } // if-then-else
  
  return c;
}

NTL_CLOSE_NNS
