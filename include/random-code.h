//---------------------------------------------------------------------------------------------------------------//
// Copyright (C) 2012 Tyson Whitehead
//
// This code is free software; you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation; either version 2, or (at your option) any later
// version.
//
//---------------------------------------------------------------------------------------------------------------//

#if _POSIX_C_SOURCE < 200112
#error _POSIX_C_SOURCE 200112 or greater required
#endif


#ifndef RANDOM_CODE
#define RANDOM_CODE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "random-type.h"
#include "random-func.h"
#include "random-data.h"

#include "error-func.h"
#include "error-code.h"


//---------------------------------------------------------------------------------------------------------------//
// MersenneTwister_begin: () -> (MersenneTwister)
//
// MersenneTwister initialized from /dev/urandom.
//
MersenneTwister MersenneTwister_begin() {
  MersenneTwister mersennetwister;

  // Allocate space for state
  if ( !(mersennetwister = malloc(sizeof(struct _MersenneTwister))) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for MersenneTwister", 
                   sizeof(struct _MersenneTwister));

  // Initial state with /dev/urandom
  mersennetwister = MersenneTwister_urandom(mersennetwister);

  return mersennetwister;
}


// MersenneTwister_end: (*MersenneTwister) -> ()
//
// Release resources associated with mersennetwister.
//
void MersenneTwister_end(MersenneTwister const mersennetwister) {
  free(mersennetwister);
}


// MersenneTwister_urandom: (*MersenneTwister) -> (MersenneTwister)
//
// Seed the mersennetwister state using /dev/urandom.
//
MersenneTwister MersenneTwister_urandom(MersenneTwister const mersennetwister) {
  FILE *file;
  UInt index;

  // Open /dev/urandom
  if ( !(file = fopen("/dev/urandom","r")) )
    Error_dieErrNo(1, "unable to open /dev/urandom for reading");

  // Fill the state array (avoiding all all zeros)
  do {
    // Fill the state array
    if ( fread(mersennetwister->state, sizeof(UInt32), MERSENNETWISTER_N, file) != MERSENNETWISTER_N )
      Error_dieErrNo(1, "unable to read %tu bytes from /dev/urandom",
                     sizeof(mersennetwister->state)/sizeof(mersennetwister->state[0]));

    // Find first non-zero entry (implies not all zero if it exists)
    for (index=0; index < MERSENNETWISTER_N && mersennetwister->state[index]==0; index+=1);
  } while (index == MERSENNETWISTER_N);

  // Close /dev/urandom
  if (fclose(file))
    Error_dieErrNo(1, "unable to close /dev/urandom");

  // Reset the index
  mersennetwister->index = 0;

  return mersennetwister;
}


// MersenneTwister_seed: (*MersenneTwister,UInt32) -> (MersenneTwister)
//
// Seed the mersennertwister state using the Matsumoto and Nishimura's initialization generator
//
//   x_k = 1812433253 * (x_{k-1} ^ x_{k-1}>>30) + k
//
// where x_0 = seed.  This cannot generate a zero state as x_0 = 0 gives x_1 = 1.
//
// The 1812433253 multiplier is from Brosh and Niederreiter paper on optimal multipliers for linear congruential
// generators (LCG).[2,3]  This is not a linear congruential generator though as Matsumoto has added extra bits.
//
// The shifted bit (added in 2002) ensures the most significant bits of the seed also effect the least
// significant bits of x_k.  The addition of the index is presumably to make a zero state impossible.
//
// As there are only 2^32-1 unique seed (compared to the 32^624-1 unique seeds for the Meresenne Twister),
// this routine should only be used for cases where it is necessary to regenerate a stream, such as debugging.
//
MersenneTwister MersenneTwister_seed(MersenneTwister const mersennetwister, UInt32 const seed) {
  // Generate state from the seed (can't give an entire state of zero)
  mersennetwister->state[0] = seed;

  for (UInt index = 1; index < MERSENNETWISTER_N; index+= 1)
    mersennetwister->state[index] =
      1812433253*(mersennetwister->state[index-1] ^ mersennetwister->state[index-1]>>30) + index;

  // Reset the index
  mersennetwister->index = 0;

  return mersennetwister;
}


// MersenneTwister_next: (*MersenneTwister) -> (MersenneTwister)
//
// Generate the next chunk of state via the recursion equation
//
//   x_{k+n} = x_{k+m} ^ t(x_k & 0x80000000 | x_{k+1} & 0x7fffffff)
//
// where t(x) = x>>1 ^ (x&1 ? A : 0).
//
MersenneTwister MersenneTwister_next(MersenneTwister const mersennetwister) {
  // Generate next bit of state (unwound into three loops to avoid having to modulus the index)
  UInt index;
  UInt32 value;

  for (index=0; index<MERSENNETWISTER_N-MERSENNETWISTER_M; index+=1){
    value = (mersennetwister->state[index]&0x80000000) | (mersennetwister->state[index+1]&0x7fffffff);
    mersennetwister->state[index] =
      mersennetwister->state[index+MERSENNETWISTER_M] ^ value>>1 ^ (value&1 ? MERSENNETWISTER_A : 0);
  }
  for (; index<MERSENNETWISTER_N-1; index+=1) {
    value = (mersennetwister->state[index]&0x80000000) | (mersennetwister->state[index+1]&0x7fffffff);
    mersennetwister->state[index] =
      mersennetwister->state[index-(MERSENNETWISTER_N-MERSENNETWISTER_M)] ^
      value>>1 ^ (value&1 ? MERSENNETWISTER_A : 0);
  }
  value = (mersennetwister->state[MERSENNETWISTER_N-1]&0x80000000) | (mersennetwister->state[0]&0x7fffffff);
  mersennetwister->state[MERSENNETWISTER_N-1] =
    mersennetwister->state[MERSENNETWISTER_M-1] ^ value>>1 ^ (value&1 ? MERSENNETWISTER_A : 0);

  // Reset extraction point
  mersennetwister->index = 0;

  return mersennetwister;
}


// MersenneTwister_extract_UInt32: (*MersenneTwister) -> (MersenneTwister, UInt32)
//
// Return 32b of the current state, after tempering with
//
//   y_k = x_k ^  (x_k >> U)
//   y_k = y_k ^ ((y_k << S) & B)
//   y_k = y_k ^ ((y_k << T) & C)
//   y_k = y_k ^  (y_k >> L)
//
// and advance to next bit of state (generating more if requried).
//
MersenneTwister_UInt32 MersenneTwister_extract_UInt32(MersenneTwister mersennetwister) {
  // Remix the state if we've used it all
  if (mersennetwister->index >= MERSENNETWISTER_N)
    mersennetwister = MersenneTwister_next(mersennetwister);

  // Generate the random value (see paper)
  UInt32 value;

  value = mersennetwister->state[mersennetwister->index];
  mersennetwister->index += 1;

  value ^= value>>MERSENNETWISTER_U;
  value ^= value<<MERSENNETWISTER_S & MERSENNETWISTER_B;
  value ^= value<<MERSENNETWISTER_T & MERSENNETWISTER_C;
  value ^= value>>MERSENNETWISTER_L;
  
  return pack_MersenneTwister_UInt32(mersennetwister, value);
}


//---------------------------------------------------------------------------------------------------------------//
// Random_uniform_UInt32: (*MersenneTwister,UInt32) -> (MersenneTwister, UInt32)
//
// Uniform random integer on [0,n) (via rejecting last tail and clipping to range)
//
MersenneTwister_UInt32 Random_uniform_UInt32(MersenneTwister mersennetwister, UInt32 const n) {
  UInt32 value;

  do
    unpack_MersenneTwister_UInt32( &mersennetwister, &value, 
                                   MersenneTwister_extract_UInt32(mersennetwister) );
  while (value/n >= UINT32_MAX/n);

  value %= n;

  return pack_MersenneTwister_UInt32(mersennetwister, value);
}


// Random_binomial_UInt32: (*MersenneTwister, UInt32, Float32) -> (MersenneTwister, UInt32)
//
// A binomial integer on [0,n] (via Kachitvichyanukul and Schmesier's accept/reject criteria)[4]
//
MersenneTwister_UInt32 Random_binomial_UInt32(MersenneTwister mersennetwister,
                                                    UInt32 const n, Float32 const p) {
  // Calculate parameter constants (caching these can speed things up)
  Float32 r,q, nr,nrq;
  Float32 s,a, f0;
  Float32 x_left,x_med,x_right, l_left,l_right;
  Float32 c, prob1,prob2,prob3,prob4;
  UInt32 med;

  r = fminf(p, 1.0-p);
  q = 1.0-r;
  nr = n*r;

  s = r/q;
  a = s*(n+1.0);
  f0 = powf(q,n);

  if(nr >= 10.0) {
    Float32 temp1,temp2;

    nrq = nr*q;

    temp1 = nr+r;
    med = (UInt32)floorf(temp1);
    prob1 = floorf(2.195*sqrtf(nrq)-4.6*q)+0.5;
    x_med = med+0.5;
    x_left = x_med-prob1;
    x_right = x_med+prob1;

    temp2 = (temp1-x_left)/(temp1-x_left*r);
    l_left = temp2*(1.0+0.5*temp2);
    temp2 = (x_right-temp1)/(x_right*q);
    l_right = temp2*(1.0+0.5*temp2);

    c = 0.134+20.5/(15.3+med);
    prob2 = prob1*(1.0+2.0*c);
    prob3 = prob2+c/l_left;
    prob4 = prob3+c/l_right;
  }

  // Generate the binomial RV
  UInt32 y;

  // Use BINV (Binomial CDF Inverse) for nr < 10.0
  if (nr < 10.0) {
    Float32 u;
    Float32 prob = f0;

    unpack_MersenneTwister_Float32( &mersennetwister, &u,
                                    Random_uniform_Float32(mersennetwister) );

    for (y = 0; u > prob && y < n; ++y){
      u -= prob;
      prob *= (a/(y+1)-s);
    }
  }
  // Use BTPE (Binomial, Triangle, Parallelogram, Exponential) for nr >= 10.0
  else {
    Float32 u;
    Float32 v;

    unpack_MersenneTwister_Float32( &mersennetwister,&u,
                                    Random_uniform_Float32(mersennetwister) );
    u *= prob4;
    unpack_MersenneTwister_Float32( &mersennetwister,&v,
                                    Random_uniform_Float32(mersennetwister) );

    // This region falls under the scaled PDF (see diagram in paper)
    if (u <= prob1) {
      // Region 1 -- triangular
      y = (UInt32)floor(x_med-prob1*v+u);
    }
    // These regions may not fall under the scaled PDF (see diagram in paper)
    else {
      // Region 2 -- parallelograms
      if (u <= prob2) {
        Float32 const x = x_left+(u-prob1)/c;
        v = v*c+1.0-fabs(med-x+0.5)/prob1;
        if (v > 1.0)
          return Random_binomial_UInt32(mersennetwister, n, p);
        y = (UInt32)floorf(x);
      }
      // Region 3 -- left exponential tail
      else if (u <= prob3) {
        Float32 const x = x_left+logf(v)/l_left;
        if (x < 0)
          return Random_binomial_UInt32(mersennetwister, n, p);
        v *= (u-prob2)*l_left;
        y = (UInt32)floorf(x);

      }
      // Region 4 -- right exponential tail
      else {
        Float32 const x = x_right-logf(v)/l_right;
        if (x > n+1.0)
          return Random_binomial_UInt32(mersennetwister, n, p);
        v *= (u-prob3)*l_right;
        y = (UInt32)floorf(x);
      }

      // Accept/reject tests
      Float32 const ydiff = fabs(y-med);

      if (ydiff <= 20 || ydiff >= nrq/2.0-1.0) {
        // Evaluate f(y) via f(y) = f(y-1)(a/x-s) starting at the mode
        Float32 f = 1.0;

        if (med < y)
          for (UInt32 i = med+1; i <= y; ++i)
            f *= a/i-s;
        else if (med > y)
          for (UInt32 i = y+1; i <= med; ++i)
            f /= a/i-s;

        if (v > f)
          return Random_binomial_UInt32(mersennetwister, n, p);
      }
      else {
        // Check the value of ln(v) against upper and lower bounds of ln(f(y))
        Float32 const bound = (ydiff/nrq)*((ydiff*(ydiff/3.0+0.625)+1.0/6.0)/nrq+0.5);
        Float32 const y2norm = -ydiff*ydiff/(2.0*nrq);
        Float32 const vln = logf(v);

        if (vln > y2norm+bound)
          return Random_binomial_UInt32(mersennetwister, n, p);
        if (vln >= y2norm-bound) {
          // Stirling's approximation (within machine accuracy)
          Float32 const y1 = y+1, y2 = y1*y1;
          Float32 const f1 = med+1, f2 = f1*f1;
          Float32 const z1 = n+1-med, z2 = z1*z1;
          Float32 const w1 = n-y+1, w2 = w1*w1;

          if (vln > (x_med*logf(f1/y1)+(n-med+0.5)*logf(z1/w1)+(y-med)*logf(w1*r/(y1*q))
                     +(13860.0-(462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0
                     +(13860.0-(462.0-(132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z1/166320.0
                     +(13860.0-(462.0-(132.0-(99.0-140.0/y2)/y2)/y2)/y2)/y1/166320.0
                     +(13860.0-(462.0-(132.0-(99.0-140.0/w2)/w2)/w2)/w2)/w1/166320.0))
            return Random_binomial_UInt32(mersennetwister, n, p);
        }
      }
    }
  }

  if (p > 0.5)
    y = n-y;

  return pack_MersenneTwister_UInt32(mersennetwister, y);
}


// Random_uniform_Float32: (*MersenneTwister) -> (MersenneTwister, Float32)
//
// Uniform random floating point number on [0,1) (by shifting in uniform integer to precision).
//
MersenneTwister_Float32 Random_uniform_Float32(MersenneTwister mersennetwister) {
  // Shift in bits after the decimal place until mantissa is filled
  Float32 x, p;

  x = 0.0;
  p = 1.0;

  do {
    // Generate random 32b
    UInt32 u;

    unpack_MersenneTwister_UInt32( &mersennetwister, &u,
                                   MersenneTwister_extract_UInt32(mersennetwister) );

    // Shift in leading 23b (this fills all the bits in the mantissa)
    p *= 0x1p-23;
    x += ( u >> 9 & (1<<23)-1 )*p;

    // Shift in next 9b (exponent will absorbs leading zeros leaving room for more)
    p *= 0x1p-9;
    x += ( u >> 0 & (1<<9)-1 )*p;
  } while ( x + p*0x1p-1 != x );

  return pack_MersenneTwister_Float32(mersennetwister, x);
}


// Random_normal2_Float32: (*MersenneTwister) -> (MersenneTwister, Float32,Float32)
//
// Two floating point N(0,1) numbers (via the polar Box-Muller transformation).
//
MersenneTwister_Float32_Float32 Random_normal2_Float32(MersenneTwister mersennetwister) {
  // Generate a uniform draw on an open unit circle (excluding the origin)
  Float32 x,y, r;

  do {
    unpack_MersenneTwister_Float32( &mersennetwister, &x,
                                    Random_uniform_Float32(mersennetwister) );
    unpack_MersenneTwister_Float32( &mersennetwister, &y,
                                    Random_uniform_Float32(mersennetwister) );

    x = 2.0*x-1.0;
    y = 2.0*y-1.0;

    r = x*x+y*y;
  } while (r == 0.0 || r > 1.0);

  // Convert into two N(0,1) numbers
  r = sqrtf(-2.0*logf(r)/r);

  return pack_MersenneTwister_Float32_Float32(mersennetwister, x*r,y*r);
}


//---------------------------------------------------------------------------------------------------------------//
// Tuples
MersenneTwister_UInt32 pack_MersenneTwister_UInt32(MersenneTwister const first, UInt32 const second) {
  MersenneTwister_UInt32 const value = { .first = first, .second = second };
  return value;
}

MersenneTwister_Float32 pack_MersenneTwister_Float32(MersenneTwister const first, Float32 const second) {
  MersenneTwister_Float32 const value = { .first = first, .second = second };
  return value;
}

MersenneTwister_Float32_Float32 pack_MersenneTwister_Float32_Float32
(MersenneTwister const first, Float32 const second, Float32 const third) {
  MersenneTwister_Float32_Float32 const value = { .first = first, .second = second, .third = third };
  return value;
}


//
void unpack_MersenneTwister_UInt32(MersenneTwister* const first, UInt32* const second,
                                   MersenneTwister_UInt32 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
}

void unpack_MersenneTwister_Float32(MersenneTwister* const first, Float32* const second,
                                    MersenneTwister_Float32 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
}

void unpack_MersenneTwister_Float32_Float32(MersenneTwister* const first, Float32* const second,
                                            Float32* const third,
                                            MersenneTwister_Float32_Float32 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
}


//---------------------------------------------------------------------------------------------------------------//

#endif // RANDOM_CODE
