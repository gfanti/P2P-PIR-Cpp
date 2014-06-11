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
//
//  The FXY.cc and FXY.h files are adapted from files in version
//  5.4 of the NTL package by Victor Shoup <victor@shoup.net>, with
//  modifications by M. Jason Hinek <mjhinek@alumni.uwaterloo.ca> and
//  Ian Goldberg <iang@cs.uwaterloo.ca>.

#ifndef __FXY_H__
#define __FXY_H__

#include <NTL/ZZ_pX.h>
#include <NTL/GF2EX.h>

NTL_OPEN_NNS

#define FTXY FXY<F,FX,vec_FX>



/************************************************************

                         FXY

The class FXY implements bivariate polynomial arithmetic over a field F.
Polynomials are represented as vec_FX's.
If f is a FXY, then f.rep is a vec_FX.
The zero polynomial is represented as a zero length vector.
Otherwise. f.rep[0] is the constant-term, and f.rep[f.rep.length()-1]
is the leading coefficient, which is always non-zero.
The member f.rep is public, so the vector representation is fully
accessible.
Use the member function normalize() to strip leading zeros.

**************************************************************/

template <class F, class FX, class vec_FX>
class FXY {

public:

typedef vec_FX VectorBaseType; 


vec_FX rep;


/***************************************************************

          Constructors, Destructors, and Assignment

****************************************************************/


FXY()
//  initial value 0

   { }


FXY(INIT_SIZE_TYPE, long n) { rep.SetMaxLength(n); }

FXY(const FXY& a) : rep(a.rep) { }
// initial value is a


FXY& operator=(const FXY& a) 
   { rep = a.rep; return *this; }

~FXY() { }

// strip leading zeros
void normalize(){
   long n;
   FX* p;

   // First normalize the coefficients
   n = rep.length();
   if (n == 0) return;
   p = rep.elts() + n;
   while (n > 0) {
       --p;
       p->normalize();
       --n;
   }
   
   n = rep.length();
   p = rep.elts() + n;
   while (n > 0 && IsZero(*--p)) {
      n--;
   }
   rep.SetLength(n);
}

void SetMaxLength(long n) 
// pre-allocate space for n coefficients.
// Value is unchanged

   { rep.SetMaxLength(n); }


void kill() 
// free space held by this polynomial.  Value becomes 0.

   { rep.kill(); }

static const FXY& zero() {
   static FXY z;
   return z;
}


FXY(FXY& x, INIT_TRANS_TYPE) : rep(x.rep, INIT_TRANS) { }

inline FXY(long i, const FX& c);
inline FXY(long i, const F& c);
inline FXY(long i, long c);

FXY& operator=(long a) { conv(*this, a); return *this; }
FXY& operator=(const F& a) { conv(*this, a); return *this; }
FXY& operator=(const FX& a) { conv(*this, a); return *this; }


};



/********************************************************************

                           input and output

I/O format:

   [a_0 a_1 ... a_n],

represents the polynomial a_0 + a_1*Y + ... + a_n*Y^n.

On output, all coefficients will be polynomials of the type FX.
amd a_n not zero (the zero polynomial is [ ]).
On input, the coefficients are arbitrary polynomials which are
then reduced modulo p, and leading zeros stripped.

*********************************************************************/


template <class F, class FX, class vec_FX>
NTL_SNS istream& operator>>(NTL_SNS istream& s, FXY<F,FX,vec_FX>& x)
{
   s >> x.rep;
   x.normalize();
   return s;
}

template <class F, class FX, class vec_FX>
NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const FXY<F,FX,vec_FX>& a)
{
   return s << a.rep;
}




/**********************************************************

                   Some utility routines

***********************************************************/



template <class F, class FX, class vec_FX>
inline long deg(const FXY<F,FX,vec_FX>& a) { return a.rep.length() - 1; }
// degree of a polynomial in Z_p[X][Y] (ie, degree wrt y).
// note that the zero polynomial has degree -1.

template <class F, class FX, class vec_FX>
inline long degY(const FXY<F,FX,vec_FX>& a) { return a.rep.length() - 1; }
// degree of a polynomial in Z_p[X][Y] (ie, degree wrt y) 
// note that the zero polynomial has degree -1.

// degree of a polynomial in Z_p[Y][X] instead of Z_p[X][Y] (ie, degree wrt x)
// note that the zero polynomial has degree -1.
template <class F, class FX, class vec_FX>
long degX(const FXY<F,FX,vec_FX>& a)
{
  //
  // degree of polynomial in Z_p[Y][X] instead of Z_p[X][Y] 
  // 
  //
  
  long degY = a.rep.length() - 1;
  long degX = -1;

  if(degY > -1){
    degX = deg(a.rep[0]);
    for(long i = 1; i <= degY; i++){
      if( deg(a.rep[i]) > degX )
	degX = deg(a.rep[i]);
    }
    
  }
  
  return degX;
       
}



// zero if i not in range
template <class F, class FX, class vec_FX>
const FX& coeff(const FXY<F,FX,vec_FX>& a, long i)
{
   if (i < 0 || i > deg(a))
      return FX::zero();
   else
      return a.rep[i];
}

// x = a[i], or zero if i not in range
template <class F, class FX, class vec_FX>
void GetCoeff(FX& x, const FXY<F,FX,vec_FX>& a, long i)
{
   if (i < 0 || i > deg(a))
      clear(x);
   else
      x = a.rep[i];
}

// zero if a == 0
template <class F, class FX, class vec_FX>
const FX& LeadCoeff(const FXY<F,FX,vec_FX>& a)
{
   if (IsZero(a))
      return FX::zero();
   else
      return a.rep[deg(a)];
}

// zero if a == 0
template <class F, class FX, class vec_FX>
const FX& ConstTerm(const FXY<F,FX,vec_FX>& a)
{
   if (IsZero(a))
      return FX::zero();
   else
      return a.rep[0];
}

// x[i][j] = a, error is raised if i < 0
template <class F, class FX, class vec_FX>
void SetCoeff(FXY<F,FX,vec_FX>& x, long i, long j, const F& a)
{
   long k, m;

   if (i < 0) 
      Error("SetCoeff: negative index");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m) {
      /* careful: a may alias a coefficient of x */

#if 0
      long alloc = x.rep.allocated();

      if (alloc > 0 && i >= alloc) {
         FXTemp aa_tmp;  FX& aa = aa_tmp.val();
         aa = a;
         x.rep.SetLength(i+1);
         x.rep[i] = aa;
      }
      else {
#endif
         x.rep.SetLength(i+1);
         SetCoeff(x.rep[i], j, a);
#if 0
      }
#endif

      for (k = m+1; k < i; k++)
         clear(x.rep[k]);
   }
   else
	SetCoeff(x.rep[i], j, a);

   x.normalize();
}

// x[i] = a, error is raised if i < 0
template <class F, class FX, class vec_FX>
void SetCoeff(FXY<F,FX,vec_FX>& x, long i, const FX& a)
{
   long j, m;

   if (i < 0) 
      Error("SetCoeff: negative index");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m) {
      /* careful: a may alias a coefficient of x */

#if 0
      long alloc = x.rep.allocated();

      if (alloc > 0 && i >= alloc) {
         FXTemp aa_tmp;  FX& aa = aa_tmp.val();
         aa = a;
         x.rep.SetLength(i+1);
         x.rep[i] = aa;
      }
      else {
#endif
         x.rep.SetLength(i+1);
         x.rep[i] = a;
#if 0
      }
#endif

      for (j = m+1; j < i; j++)
         clear(x.rep[j]);
   }
   else
      x.rep[i] = a;

   x.normalize();
}

template <class F, class FX, class vec_FX>
void SetCoeff(FXY<F,FX,vec_FX>& x, long i, long a)
{
   if (a == 1) 
      SetCoeff(x, i);
   else {
      FX T;
      conv(T, a);
      SetCoeff(x, i, T);
   }
}

// x[i] = 1, error is raised if i < 0
template <class F, class FX, class vec_FX>
void SetCoeff(FXY<F,FX,vec_FX>& x, long i)
{
   long j, m;

   if (i < 0) 
      Error("coefficient index out of range");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m) {
      x.rep.SetLength(i+1);
      for (j = m+1; j < i; j++)
         clear(x.rep[j]);
   }
   set(x.rep[i]);
   x.normalize();
}

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>::FXY(long i, const FX& a)
   { SetCoeff(*this, i, a); } 

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>::FXY(long i, long a)
   { SetCoeff(*this, i, a); } 

// y is set to the monomial Y
template <class F, class FX, class vec_FX>
void SetY(FXY<F,FX,vec_FX>& y)
{
   clear(y);
   SetCoeff(y, 1);
}

// test if a = Y
template <class F, class FX, class vec_FX>
long IsY(const FXY<F,FX,vec_FX>& a)
{
   return deg(a) == 1 && IsOne(LeadCoeff(a)) && IsZero(ConstTerm(a));
}

template <class F, class FX, class vec_FX>
inline void clear(FXY<F,FX,vec_FX>& x) 
// x = 0

   { x.rep.SetLength(0); }

template <class F, class FX, class vec_FX>
inline void set(FXY<F,FX,vec_FX>& x)
// x = 1

   { x.rep.SetLength(1); set(x.rep[0]); }

template <class F, class FX, class vec_FX>
inline void swap(FXY<F,FX,vec_FX>& x, FXY<F,FX,vec_FX>& y)
// swap x & y (only pointers are swapped)

   { swap(x.rep, y.rep); }

template <class F, class FX, class vec_FX>
void random(FXY<F,FX,vec_FX>& x, long n);
template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> random_FXY(long n)
   { FXY<F,FX,vec_FX> x; random(x, n); NTL_OPT_RETURN(FTXY, x); }
// generate a random polynomial of degree < n 

template <class F, class FX, class vec_FX>
void trunc(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, long m);
// x = a % Y^m

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> trunc(const FXY<F,FX,vec_FX>& a, long m)
   { FXY<F,FX,vec_FX> x; trunc(x, a, m); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
void RightShift(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, long n);
// x = a/Y^n

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> RightShift(const FXY<F,FX,vec_FX>& a, long n)
   { FXY<F,FX,vec_FX> x; RightShift(x, a, n); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
void LeftShift(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, long n);
// x = a*Y^n

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> LeftShift(const FXY<F,FX,vec_FX>& a, long n)
   { FXY<F,FX,vec_FX> x; LeftShift(x, a, n); NTL_OPT_RETURN(FTXY, x); }

// a =  X^n * P(X) 
template <class F, class FX, class vec_FX>
FX shiftX(const FX& P, long n)
{
  // a =  X^n * P(X) 
  // n >= 0 for now

  long deg = P.rep.length()-1;
  
  FX x;
  x.rep.SetLength( deg+1+n );

  for(long i = 0; i < n; i++){
    x.rep[i] = 0;
  }
  for(long i = 0; i <= deg; i++){
    x.rep[i+n] = P.rep[i];
  }
  return x;
}

// returns  P(X) / X^n
template <class F, class FX, class vec_FX>
FXY<F,FX,vec_FX> backShiftX(const FXY<F,FX,vec_FX>& P, long n)
{
  // returns  P(X) / X^n

  long deg = P.rep.length() - 1; 
  long degX;

  FXY<F,FX,vec_FX> x;
  x.rep.SetLength( deg+1 );

  for(long i = 0; i <= deg; i++){ // work on each coefficient of Y
    if( IsZero(P.rep[i]) ){
      x.rep[i] = P.rep[i];
    }else{
      degX = (P.rep[i].rep.length()-1)-n;
      x.rep[i].rep.SetLength( degX + 1);
      for(long j = 0; j <= degX; j++){
	x.rep[i].rep[j] = P.rep[i].rep[j+n];
      }
    }
  }
  x.normalize();
  return x;
}

#ifndef NTL_TRANSITION

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator>>(const FXY<F,FX,vec_FX>& a, long n)
   { FXY<F,FX,vec_FX> x; RightShift(x, a, n); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator<<(const FXY<F,FX,vec_FX>& a, long n)
   { FXY<F,FX,vec_FX> x; LeftShift(x, a, n); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator<<=(FXY<F,FX,vec_FX>& x, long n)
   { LeftShift(x, x, n); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator>>=(FXY<F,FX,vec_FX>& x, long n)
   { RightShift(x, x, n); return x; }

#endif



template <class F, class FX, class vec_FX>
void diff(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a);
// x = derivative of a wrt Y

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> diff(const FXY<F,FX,vec_FX>& a)
   { FXY<F,FX,vec_FX> x; diff(x, a); NTL_OPT_RETURN(FTXY, x); }


template <class F, class FX, class vec_FX>
void MakeMonic(FXY<F,FX,vec_FX>& x)
{
    if (IsZero(x))
	return;

    F lc = LeadCoeff(LeadCoeff(x));
    if (lc == 1)
	return;

    F lcinv;
    inv(lcinv, lc);
    x *= lcinv;
}

template <class F, class FX, class vec_FX>
void reverse(FXY<F,FX,vec_FX>& c, const FXY<F,FX,vec_FX>& a, long hi);

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> reverse(const FXY<F,FX,vec_FX>& a, long hi)
   { FXY<F,FX,vec_FX> x; reverse(x, a, hi); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline void reverse(FXY<F,FX,vec_FX>& c, const FXY<F,FX,vec_FX>& a)
{  reverse(c, a, deg(a)); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> reverse(const FXY<F,FX,vec_FX>& a)
   { FXY<F,FX,vec_FX> x; reverse(x, a); NTL_OPT_RETURN(FTXY, x); }

// Returns the degree of the rightmost term whose coefficient has the highest
// degree.  Returns -1 if the given polynomial is 0;
template <class F, class FX, class vec_FX>
long pivotIndex (const FXY<F,FX,vec_FX>& a) {
	long max_index = -1;
	long max_degree = -1;
	for (long i = 0; i < a.rep.length(); ++i) {
		long d = deg(a.rep[i]);
		if (d >= max_degree) {
			max_index = i;
			max_degree = d;
		}
	}
	return max_index;
}

template <class F, class FX, class vec_FX>
FXY<F,FX,vec_FX> swapVarOrder (const FXY<F,FX,vec_FX>& a) {
	FXY<F,FX,vec_FX> new_a;
	for (long i = 0; i < a.rep.length(); ++i) {
		for (long j = 0; j <= deg(a.rep[i]); ++j) {
			if (j > deg(new_a)) {
				SetCoeff(new_a, j, FX(i, coeff(a.rep[i], j)));
			} else {
				SetCoeff(new_a.rep[j], i, coeff(a.rep[i], j));
			}
		}
	}
	return new_a;
}

template <class F, class FX, class vec_FX>
FX eval(const FXY<F,FX,vec_FX>& a, F x) {
	FX result;
	for (long i = 0; i < a.rep.length(); ++i) {
		result += power(x, i) * a.rep[i];
	}
	return result;
}


/*******************************************************************

                        conversion routines

********************************************************************/



template <class F, class FX, class vec_FX>
void conv(FXY<F,FX,vec_FX>& x, long a)
{
   if (a == 0)
      clear(x);
   else if (a == 1)
      set(x);
   else {
      FX T;
      conv(T, a);
      conv(x, T);
   }
}

template <class F, class FX, class vec_FX>
void conv(FXY<F,FX,vec_FX>& x, const F& a)
{
   if (IsZero(a))
      clear(x);
   else {
      FX T;
      conv(T, a);
      conv(x, T);
   }
}

template <class F, class FX, class vec_FX>
void conv(FXY<F,FX,vec_FX>& x, const FX& a)
{
   if (IsZero(a))
      x.rep.SetLength(0);
   else {
      x.rep.SetLength(1);
      x.rep[0] = a;

      // note: if a aliases x.rep[i], i > 0, this code
      //       will still work, since is is assumed that
      //       SetLength(1) will not relocate or destroy x.rep[i]
   }
}
template <class F, class FX, class vec_FX>
void conv(FXY<F,FX,vec_FX>& x, const vec_FX& a)
{
   x.rep = a;
   x.normalize();
}

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> to_FXY(long a)
   { FXY<F,FX,vec_FX> x; conv(x, a); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> to_FXY(const F& a)
   { FXY<F,FX,vec_FX> x; conv(x, a); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> to_FXY(const FX& a)
   { FXY<F,FX,vec_FX> x; conv(x, a); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> to_FXY(const vec_FX& a)
   { FXY<F,FX,vec_FX> x; conv(x, a); NTL_OPT_RETURN(FTXY, x); }



/*************************************************************

                        Comparison

**************************************************************/

template <class F, class FX, class vec_FX>
long IsZero(const FXY<F,FX,vec_FX>& a)
{
   return a.rep.length() == 0;
}

template <class F, class FX, class vec_FX>
long IsOne(const FXY<F,FX,vec_FX>& a)
{
    return a.rep.length() == 1 && IsOne(a.rep[0]);
}

template <class F, class FX, class vec_FX>
inline long operator==(const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
{
   return a.rep == b.rep;
}

template <class F, class FX, class vec_FX>
inline long operator!=(const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
{
   return !(a == b);
}

template <class F, class FX, class vec_FX>
long operator==(const FXY<F,FX,vec_FX>& a, long b)
{
   if (b == 0)
      return IsZero(a);

   if (b == 1)
      return IsOne(a);

   long da = deg(a);

   if (da > 0)
      return 0;

   F bb;
   bb = b;

   if (da < 0)
      return IsZero(bb);

   return a.rep[0] == bb;
}

template <class F, class FX, class vec_FX>
long operator==(const FXY<F,FX,vec_FX>& a, const F& b)
{
   if (IsZero(b))
      return IsZero(a);

   long da = deg(a);

   if (da != 0)
      return 0;

   return a.rep[0] == b;
}

template <class F, class FX, class vec_FX>
long operator==(const FXY<F,FX,vec_FX>& a, const FX& b)
{
   if (IsZero(b))
      return IsZero(a);

   long da = deg(a);

   if (da != 0)
      return 0;

   return a.rep[0] == b;
}

template <class F, class FX, class vec_FX>
inline long operator==(long a, const FXY<F,FX,vec_FX>& b) { return b == a; }
template <class F, class FX, class vec_FX>
inline long operator==(const F& a, const FXY<F,FX,vec_FX>& b) { return b == a; }
template <class F, class FX, class vec_FX>
inline long operator==(const FX& a, const FXY<F,FX,vec_FX>& b) { return b == a; }

template <class F, class FX, class vec_FX>
inline long operator!=(const FXY<F,FX,vec_FX>& a, long b) { return !(a == b); }
template <class F, class FX, class vec_FX>
inline long operator!=(const FXY<F,FX,vec_FX>& a, const F& b) { return !(a == b); }
template <class F, class FX, class vec_FX>
inline long operator!=(const FXY<F,FX,vec_FX>& a, const FX& b) { return !(a == b); }

template <class F, class FX, class vec_FX>
inline long operator!=(long a, const FXY<F,FX,vec_FX>& b) { return !(a == b); }
template <class F, class FX, class vec_FX>
inline long operator!=(const F& a, const FXY<F,FX,vec_FX>& b) { return !(a == b); }
template <class F, class FX, class vec_FX>
inline long operator!=(const FX& a, const FXY<F,FX,vec_FX>& b) { return !(a == b); }


/***************************************************************

                         Addition

****************************************************************/

//
// x = a + b
//
template <class F, class FX, class vec_FX>
void add(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
{
   long degYa = deg(a);
   long degYb = deg(b);
   long minab = min(degYa, degYb);
   long maxab = max(degYa, degYb);
   x.rep.SetLength(maxab+1);

   long i;
   const FX *ap, *bp; 
   FX* xp;

   for (i = minab+1, ap = a.rep.elts(), bp = b.rep.elts(), xp = x.rep.elts();
        i; i--, ap++, bp++, xp++)
      add(*xp, (*ap), (*bp));

   if (degYa > minab && &x != &a)
      for (i = degYa-minab; i; i--, xp++, ap++)
         *xp = *ap;
   else if (degYb > minab && &x != &b)
      for (i = degYb-minab; i; i--, xp++, bp++)
         *xp = *bp;
   else
      x.normalize();
}

template <class F, class FX, class vec_FX>
void add(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const FX& b)
{
  long n = a.rep.length(); 
  if (n == 0) {
    conv(x, b);
  }
  else if (&x == &a) {
    add(x.rep[0], a.rep[0], b);
    x.normalize();
  }
  else if (x.rep.MaxLength() == 0) {
    x = a;
    add(x.rep[0], a.rep[0], b);
    x.normalize();
  }
  else {
    // ugly...b could alias a coeff of x
    
    FX *xp = x.rep.elts();
    add(xp[0], a.rep[0], b);
    x.rep.SetLength(n);
    xp = x.rep.elts();
    const FX *ap = a.rep.elts();
    long i;
    for (i = 1; i < n; i++)
      xp[i] = ap[i];
    x.normalize();
  }
}

template <class F, class FX, class vec_FX>
void add(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const F& b)
{
   if (a.rep.length() == 0) {
      conv(x, b);
   }
   else {
      if (&x != &a) x = a;
      add(x.rep[0], x.rep[0], b);
      x.normalize();
   }
}

template <class F, class FX, class vec_FX>
void add(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, long b)
{
   if (a.rep.length() == 0) {
      conv(x, b);
   }
   else {
      if (&x != &a) x = a;
      add(x.rep[0], x.rep[0], b);
      x.normalize();
   }
}

template <class F, class FX, class vec_FX>
inline void add(FXY<F,FX,vec_FX>& x, const FX& a, const FXY<F,FX,vec_FX>& b) { add(x, b, a); }
template <class F, class FX, class vec_FX>
inline void add(FXY<F,FX,vec_FX>& x, const F&  a, const FXY<F,FX,vec_FX>& b) { add(x, b, a); }
template <class F, class FX, class vec_FX>
inline void add(FXY<F,FX,vec_FX>& x, long         a, const FXY<F,FX,vec_FX>& b) { add(x, b, a); }


//
// x = a - b
//
template <class F, class FX, class vec_FX>
void sub(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
{
   long da = deg(a);
   long db = deg(b);
   long minab = min(da, db);
   long maxab = max(da, db);
   x.rep.SetLength(maxab+1);

   long i;
   const FX *ap, *bp; 
   FX* xp;

   for (i = minab+1, ap = a.rep.elts(), bp = b.rep.elts(), xp = x.rep.elts();
        i; i--, ap++, bp++, xp++)
      sub(*xp, (*ap), (*bp));

   if (da > minab && &x != &a)
      for (i = da-minab; i; i--, xp++, ap++)
         *xp = *ap;
   else if (db > minab)
      for (i = db-minab; i; i--, xp++, bp++)
         negate(*xp, *bp);
   else
      x.normalize();

}

template <class F, class FX, class vec_FX>
void sub(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const FX& b)
{
   long n = a.rep.length();
   if (n == 0) {
      conv(x, b);
      negate(x, x);
   }
   else if (&x == &a) {
      sub(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else if (x.rep.MaxLength() == 0) {
      x = a;
      sub(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else {
      // ugly...b could alias a coeff of x

      FX *xp = x.rep.elts();
      sub(xp[0], a.rep[0], b);
      x.rep.SetLength(n);
      xp = x.rep.elts();
      const FX *ap = a.rep.elts();
      long i;
      for (i = 1; i < n; i++)
         xp[i] = ap[i];
      x.normalize();
   }
}

template <class F, class FX, class vec_FX>
void sub(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const F& b)
{
   if (IsZero(b)) {
      x = a;
      return;
   }

   if (a.rep.length() == 0) {
      x.rep.SetLength(1);
      x.rep[0] = b;
      negate(x.rep[0], x.rep[0]);
   }
   else {
      if (&x != &a) x = a;
      sub(x.rep[0], x.rep[0], b);
   }
   x.normalize();
}

template <class F, class FX, class vec_FX>
void sub(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, long b)
{
   if (b == 0) {
      x = a;
      return;
   }

   if (a.rep.length() == 0) {
      x.rep.SetLength(1);
      x.rep[0] = b;
      negate(x.rep[0], x.rep[0]);
   }
   else {
      if (&x != &a) x = a;
      sub(x.rep[0], x.rep[0], b);
   }
   x.normalize();
}

template <class F, class FX, class vec_FX>
void sub(FXY<F,FX,vec_FX>& x, const FX& a, const FXY<F,FX,vec_FX>& b)
{
  FXY<F,FX,vec_FX> temp;
  conv(temp,a);
  
  negate(x, b);
  add(x, x, temp);
}

template <class F, class FX, class vec_FX>
void sub(FXY<F,FX,vec_FX>& x, const F&  a, const FXY<F,FX,vec_FX>& b)
{
  FXY<F,FX,vec_FX> temp;
  conv(temp,a);
  
  negate(x, b);
  add(x, x, temp);
}

template <class F, class FX, class vec_FX>
void sub(FXY<F,FX,vec_FX>& x, long a, const FXY<F,FX,vec_FX>& b)
{
  FXY<F,FX,vec_FX> temp;
  conv(temp,a);
  
  negate(x, b);
  add(x, x, temp);
}


//
// x = -a
//
template <class F, class FX, class vec_FX>
void negate(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a)
{
   long n = a.rep.length();
   x.rep.SetLength(n);

   const FX* ap = a.rep.elts();
   FX* xp = x.rep.elts();
   long i;

   for (i = n; i; i--, ap++, xp++)
      negate((*xp), (*ap));

}


//
// operator versions for add, sub and negate
//
template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator+(const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; add(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator+(const FXY<F,FX,vec_FX>& a, const FX& b)
   { FXY<F,FX,vec_FX> x; add(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator+(const FXY<F,FX,vec_FX>& a, const F& b)
   { FXY<F,FX,vec_FX> x; add(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator+(const FXY<F,FX,vec_FX>& a, long b)
   { FXY<F,FX,vec_FX> x; add(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator+(const FX& a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; add(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator+(const F&  a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; add(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator+(long a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; add(x, a, b); NTL_OPT_RETURN(FTXY, x); }


template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; sub(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(const FXY<F,FX,vec_FX>& a, const FX& b)
   { FXY<F,FX,vec_FX> x; sub(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(const FXY<F,FX,vec_FX>& a, const F& b)
   { FXY<F,FX,vec_FX> x; sub(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(const FXY<F,FX,vec_FX>& a, long b)
   { FXY<F,FX,vec_FX> x; sub(x, a, b); NTL_OPT_RETURN(FTXY, x); }


template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(const FX& a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; sub(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(const F& a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; sub(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(long a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; sub(x, a, b); NTL_OPT_RETURN(FTXY, x); }


template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator+=(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& b)
   { add(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator+=(FXY<F,FX,vec_FX>& x, const FX& b)
   { add(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator+=(FXY<F,FX,vec_FX>& x, const F& b)
   { add(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator+=(FXY<F,FX,vec_FX>& x, long b)
   { add(x, x, b); return x; }


template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator-=(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& b)
   { sub(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator-=(FXY<F,FX,vec_FX>& x, const FX& b)
   { sub(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator-=(FXY<F,FX,vec_FX>& x, const F& b)
   { sub(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator-=(FXY<F,FX,vec_FX>& x, long b)
   { sub(x, x, b); return x; }


template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator-(const FXY<F,FX,vec_FX>& a) 
   { FXY<F,FX,vec_FX> x; negate(x, a); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator++(FXY<F,FX,vec_FX>& x) 
   { add(x, x, 1); return x; }

template <class F, class FX, class vec_FX>
inline void operator++(FXY<F,FX,vec_FX>& x, int) 
   { add(x, x, 1); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator--(FXY<F,FX,vec_FX>& x) 
   { sub(x, x, 1); return x; }

template <class F, class FX, class vec_FX>
inline void operator--(FXY<F,FX,vec_FX>& x, int) 
   { sub(x, x, 1); }



/*************************************************************

                      Computing y-roots

**************************************************************/

ZZ comb(long n, long k);

template <class F, class FX, class vec_FX>
void mul(FX &x, const FX &a, const ZZ &b)
{  
  F T;
  conv(T, b);  
  mul(x, a, T);  
}  

//
// returns Q(x,xy+a)/X^m, for largest m>=0 such that division is exact 
//
template <class F, class FX, class vec_FX>
FXY<F,FX,vec_FX> mapit(const FXY<F,FX,vec_FX>& Q, const F& a)
{

  // 
  // if Q(x,y) has no y terms (ie deg(Q) < 1) 
  // there is nothing to do here
  //

  if( deg(Q) < 1 ) return Q;  


  //
  // Q(x,y) has some y terms, so here comes the work
  //

  FXY<F,FX,vec_FX> temp;          // we build up the new 
                                  // polynomial in temp

  long degY = deg(Q);
  temp.rep.SetLength( degY + 1 );      

  FX coef1, coef2;             // we build up each coefficient with coef

  F aa;

  //
  // first we let y --> xy + a
  //

  //std::cout << "--- mapit ---" << endl;

  for(long i = 0; i <= degY; i++){
    coef1 = 0;
    
    for(long j = i; j <= degY; j++){
      coef2 = 0;

      power(aa,a,j-i);
      
      mul(coef2, Q.rep[j], aa);
      //std::cout << "." ;
      mul<F,FX,vec_FX>(coef2, coef2, comb(j,i));
      // std::cout << "." ;
      add(coef1,coef1,coef2);
      //std::cout << "." << endl;
      //      std::cout << ".." << coef2 << " " << aa << endl;

    } // for j

    //std::cout << "-shift " << coef1 << endl;
    coef1 = shiftX<F,FX,vec_FX>(coef1, i);
    //    std::cout << "." << coef1 << endl;
    temp.rep[i] = coef1;

  } // for i
  
  //std::cout << "start back-shift" << endl;
  //std::cout << temp << " -- " << minX(temp) << endl;
  temp.normalize();
  temp = backShiftX(temp, minX(temp));
  //std::cout << "end back-shift" << temp << endl;

  return temp;
}

// returns all roots of polynomial P(x)
template <class F, class vec_F, class FX, class vec_FX, class vec_pair_FX_long>
vec_F findroots_FX(const FX &P)
{
  

  FX Q;             //  
  Q = P;               // Q = monic version of P
  MakeMonic(Q);        //

  FX      factor;   // one square free factor

  vec_FX  factors;  // all factors (of factor) from Berlekamp
  long      nfactors;  //

  vec_F  roots;     // list of roots of Q
  long     nroots = 0; //

  F root;           // a single root of Q

  
  // first do square free decomposition 
  vec_pair_FX_long sqrfreeQ = SquareFreeDecomp(Q);
  //  std::cout << "sqrfree = " << sqrfreeQ << endl;


  // work on each factor from square free decomposition
  for(long i = 0; i < sqrfreeQ.length(); i++){
    factor = sqrfreeQ[i].a;
    //    std::cout << "factor = " << sqrfreeQ[i] << endl;

    factors = SFBerlekamp(factor);
    nfactors = factors.length();

    //    std::cout << ".." << factors << " " << nfactors << endl;

    for(long j = 0; j < nfactors; j++){
      //    std::cout << factors[j] << deg( factors[j] ) << endl;
      if( deg(factors[j]) == 1 ){
	root = FindRoot( factors[j] );
	//	std::cout << "root " << root << endl; 
	append(roots, root); 
	//	std::cout << roots[nroots] << endl; 
	nroots++; 
      }
      
    }
    
  }
  
  return roots;

}


template <class F, class FX, class vec_FX>
long minDegX(const FX& P, long max)
{
  //
  // return m >= 0 such that X^m divides P(X)
  //

  if( IsZero(P) ) return max;
  
  long min = 0;
  while( IsZero(P.rep[min]) ){
    min++;
  }
  return min;
}

template <class F, class FX, class vec_FX>
long minX(const FXY<F,FX,vec_FX>& Q)
{
  //
  // return largest m >= 0 such that X^m divides Q(x,y)
  //
  
  if( IsZero(Q) ) return 0;

  long max = degX(Q);
  long minX = minDegX<F,FX,vec_FX>(Q.rep[0],max); 
  long degreeY = deg(Q);
  long currentMinX;

  for( long i = 1; i <= degreeY; i++){
    currentMinX = minDegX<F,FX,vec_FX>(Q.rep[i], max);
    if( currentMinX < minX)
      minX = currentMinX;
  }
  
  return minX;
}

template <class F, class FX, class vec_FX>
void multiply(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& Q, long a)
{
  //
  // x = aQ (ie, Q added to itself a times)
  // 
  
  // do the easy stuff 
  if (a < 1) {
    Error("multiply: multiplier smaller than 1");
  }
  
  if ( IsOne(Q) || a == 1) {  // a1 = 1, 1Q = Q
    x = Q;
    return;
  }
  
  if (a == 2){
    add(x,Q,Q);
    return;
  }
  
  long dQ = deg(Q);
  
/*   if (dQ == 0) { */
/*     x = multiply(ConstTerm(Q), a); */
/*     return; */
/*   } */
  
  FXY<F,FX,vec_FX> res;
  res.SetMaxLength(dQ + 1);
  res = 0;
  
  long k = NumBits(a);
  long i;

  // double & add 
  for (i = k - 1; i >= 0; i--) {
    add(res,res,res); // double
    if (bit(a, i))
      add(res,res,Q); // add
  }
  
  x = res;
}

template <class F, class FX, class vec_FX>
void multiply(FX& x, const FX& Q, long a)
{
  //
  // x = aQ (ie, Q added to itself a times)
  // 
  
  // do the easy stuff 
  if (a < 1) {
    Error("multiply: multiplier smaller than 1");
  }
  
  if ( IsOne(Q) || a == 1) {  // a1 = 1, 1Q = Q
    x = Q;
    return;
  }
  
  if (a == 2){
    add(x,Q,Q);
    return;
  }
  
  long dQ = deg(Q);
  
/*   if (dQ == 0) { */
/*     x = multiply(ConstTerm(Q), a); */
/*     return; */
/*   } */
  
  FX res;
  res.SetMaxLength(dQ + 1);
  res = 0;
  
  long k = NumBits(a);
  long i;
  
  for (i = k - 1; i >= 0; i--) {
    add(res,res,res); // double
    if (bit(a, i))
      add(res,res,Q); // add
  }
  
  x = res;
}

template <class F, class FX, class vec_FX>
void multiply(FX& x, const FX& Q, ZZ a)
{
  //
  // x = aQ (ie, Q added to itself a times)
  // 
  
  // do the easy stuff 
  if( a < to_ZZ(1) ) {
    Error("multiply: multiplier smaller than 1");
  }
  
  if( IsOne(Q) || IsOne(a) ) {  // a1 = 1, 1Q = Q
    x = Q;
    return;
  }
  
  if (a == to_ZZ(2)){
    add(x,Q,Q);
    return;
  }
  
  long dQ = deg(Q);
  
/*   if (dQ == 0) { */
/*     x = multiply(ConstTerm(Q), a); */
/*     return; */
/*   } */
  
  FX res;
  res.SetMaxLength(dQ + 1);
  res = 0;
  
  long k = NumBits(a);
  long i;
  
  for (i = k - 1; i >= 0; i--) {
    add(res,res,res); // double
    if (bit(a, i))
      add(res,res,Q); // add
  }
  
  x = res;
}

/*****************************************************************

                        Multiplication

******************************************************************/

//
// x = a * b  (and associated operations)
//

// x = a * b
template <class F, class FX, class vec_FX>
void mul(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
{
  long da = deg(a);
  long db = deg(b);
  
  // do the easy cases first

  if (da < 0 || db < 0){  // a or b is zero
    clear(x); 
    return; 
  }   

  if (da == 0 && db == 0){  // a and b are FX
    x.rep.SetLength(1);     // x is FX too then
    mul( x.rep[0], a.rep[0], b.rep[0] );
    return; 
  }

  if (da == 0){ // a is FX and b is FXY
    mul(x,b,a.rep[0]);
    return;
  }

  if (db == 0){ // b is FX and a is FXY
    mul(x,a,b.rep[0]);
    return;
  }
  
  // da,db >= 1 
  
  PlainMul(x,a,b);
}

template <class F, class FX, class vec_FX>
void mul(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const FX& b)
{
  long da = deg(a);
  long db = deg(b);
  
  // do the easy cases first

  if (da < 0 || db < 0){  // a or b is zero
    clear(x); 
    return; 
  }   

  if (da == 0){  // a and b are both FX
    x.rep.SetLength(da+1);
    mul( x.rep[0], a.rep[0], b );
    return; 
  }  
  
  if (db == 0){  // b is F
    mul(x, a, b.rep[0]); 
    return; 
  }  
  
  // da >= 1 (simply multiply each coefficient of a by b)
  
  x = a;
  for(long i = 0; i <= da; i++){
    mul(x.rep[i],a.rep[i],b);
  }
  x.normalize();

}

template <class F, class FX, class vec_FX>
void mul(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const F& b)
{
   if (IsZero(b)) {
      clear(x);
      return;
   }

   if (IsOne(b)) {
      x = a;
      return;
   }

   long i, da;

   const FX *ap;
   FX* xp;

   F t;
   t = b;

   da = deg(a);
   x.rep.SetLength(da+1);
   ap = a.rep.elts();
   xp = x.rep.elts();

   for (i = 0; i <= da; i++) 
      mul(xp[i], ap[i], t);

   x.normalize();
}

template <class F, class FX, class vec_FX>
void mul(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, long b)
{
   F T;
   conv(T, b);
   mul(x, a, T);
}

// x = a^2
template <class F, class FX, class vec_FX>
void sqr(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a)
{
  mul(x,a,a);  // I'll take the easy way out for now!
}

template <class F, class FX, class vec_FX>
void sqr(FXY<F,FX,vec_FX>& x, const FX& a)
{
  FX temp;
  sqr(temp,a);  
  conv(x,temp);
}

template <class F, class FX, class vec_FX>
void sqr(FXY<F,FX,vec_FX>& x, const F& a)
{
  F temp;
  sqr(temp,a);  
  conv(x,temp);
}

template <class F, class FX, class vec_FX>
void sqr(FXY<F,FX,vec_FX>& x, long a)
{
  F ap;
  ap = a;
  F temp;
  sqr(temp, ap);
  conv(x,temp);
}

// x = a^e (e >= 0)
template <class F, class FX, class vec_FX>
void power(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, long e)
{
  // do the easy stuff 
  if (e < 0) {
    Error("power: negative exponent");
  }
  
  if (e == 0) {  // a^0 = 1
    x = 1;
    return;
  }
  
  if (a == 0 || a == 1 || e == 1) { 
    x = a;
    return;
  }
  
  if (e == 2){
    sqr(x,a);
    return;
  }

   long da = deg(a);

   if (da == 0) {
      x = power(ConstTerm(a), e);
      return;
   }

   if (da > (NTL_MAX_LONG-1)/e)
      Error("overflow in power");

   FXY<F,FX,vec_FX> res;
   res.SetMaxLength(da*e + 1);
   res = 1;
   
   long k = NumBits(e);
   long i;

   for (i = k - 1; i >= 0; i--) {
      sqr(res, res);
      if (bit(e, i))
         mul(res, res, a);
   }

   x = res;
}



// inline versions

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> sqr(const FXY<F,FX,vec_FX>& a) 
   { FXY<F,FX,vec_FX> x; sqr(x, a); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> power(const FXY<F,FX,vec_FX>& a, long e)
   { FXY<F,FX,vec_FX> x; power(x, a, e); NTL_OPT_RETURN(FTXY, x); }




// other variations (of input arguments)

template <class F, class FX, class vec_FX>
inline void mul(FXY<F,FX,vec_FX>& x, const FX& a, const FXY<F,FX,vec_FX>& b) 
   { mul(x, b, a); }

template <class F, class FX, class vec_FX>
inline void mul(FXY<F,FX,vec_FX>& x, const F& a, const FXY<F,FX,vec_FX>& b) 
   { mul(x, b, a); }

template <class F, class FX, class vec_FX>
inline void mul(FXY<F,FX,vec_FX>& x, long a, const FXY<F,FX,vec_FX>& b) 
   { mul(x, b, a); }


// operator versions

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator*(const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b) 
   { FXY<F,FX,vec_FX> x; mul(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator*(const FXY<F,FX,vec_FX>& a, const FX& b) 
   { FXY<F,FX,vec_FX> x; mul(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator*(const FXY<F,FX,vec_FX>& a, const F& b)
   { FXY<F,FX,vec_FX> x; mul(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator*(const FXY<F,FX,vec_FX>& a, long b)
   { FXY<F,FX,vec_FX> x; mul(x, a, b); NTL_OPT_RETURN(FTXY, x); }


template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator*(const FX& a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; mul(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator*(const F& a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; mul(x, a, b); NTL_OPT_RETURN(FTXY, x); }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX> operator*(long a, const FXY<F,FX,vec_FX>& b)
   { FXY<F,FX,vec_FX> x; mul(x, a, b); NTL_OPT_RETURN(FTXY, x); }



template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator*=(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& b)
   { mul(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator*=(FXY<F,FX,vec_FX>& x, const FX& b)
   { mul(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator*=(FXY<F,FX,vec_FX>& x, const F& b)
   { mul(x, x, b); return x; }

template <class F, class FX, class vec_FX>
inline FXY<F,FX,vec_FX>& operator*=(FXY<F,FX,vec_FX>& x, long b)
   { mul(x, x, b); return x; }






//
// Classical (plain) n^2 algorithms
//

template <class F, class FX, class vec_FX>
void PlainMul(FXY<F,FX,vec_FX>& x, const FXY<F,FX,vec_FX>& a, const FXY<F,FX,vec_FX>& b)
{
  long da = deg(a);
  long db = deg(b);
  
  // do the easy cases first
  
  if (da < 0 || db < 0){  // a or b is zero
    clear(x); 
    return; 
  }   
  
  if (da == 0 && db == 0){  // a and b are FX
    x.rep.SetLength(1);     // x is FX too then
    mul( x.rep[0], a.rep[0], b.rep[0] );
    return; 
  }
  
  if (da == 0){ // a is FX and b is FXY
    mul(x,b,a.rep[0]);
    return;
  }
  
  if (db == 0){ // b is FX and a is FXY
    mul(x,a,b.rep[0]);
    return;
  }
  
  // da,db >= 1 
  
  long d = da+db;  // degree of new polynomial
  
  FXY<F,FX,vec_FX> la, lb;    // create new copies of a and b
  la = a;           // just in case &x == &a or &x == &b
  lb = b;           // 

  x.rep.SetLength(d+1);
  
  long i, j, jmin, jmax;
  static FX t, accum;
  
  for (i = 0; i <= d; i++) {
    jmin = max(0, i-db);
    jmax = min(da, i);
    clear(accum);
    for (j = jmin; j <= jmax; j++) {
      mul(t, la.rep[j], lb.rep[i-j]);
      //mul(t, rep(ap[j]), rep(bp[i-j]));
      add(accum, accum, t);
    }
    SetCoeff( x, i, accum);
  }
  x.normalize();
}

typedef FXY<ZZ_p, ZZ_pX, vec_ZZ_pX> ZZ_pXY;
typedef FXY<GF2E, GF2EX, vec_GF2EX> GF2EXY;

NTL_CLOSE_NNS

#endif
