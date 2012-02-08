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


#ifndef RANDOM_DATA
#define RANDOM_DATA

#include "random-type.h"

#include "system-type.h"

// Mersenne Twister 19937[1] as described on Wikipedia (presumably based on Matsumato and Nishimura's code)
//
// [1] Makoto Matsumoto and Takuji Nishimura.  Mersenne twister: a 623-dimensionally equidistributed uniform
//     random number generator.  Transactions on Modeling and Computer Simulation, 8(1):3-30, 1998.
//
// [2] Itshak Borosh and Harald Niederreiter.  Optimal Multipliers for Pseudo-Random Number Generation by the
//     Linear Congruential Method.  BIT Numerical Mathematics, 23(1):65-74, 1983.
//
// [3] Donald Knuth.  The Art of Computer Programming.  Addison Wesley, 3rd edition, 1998.
//
// [4] Voratas Kachitvichyanukul and Bruce W. Schmesier.  Binomial Random Variate Generation.  Communications of
//     the ACM, 31(2):216-222, 1988.
//

// Mersenne Twister MT19937 parameters (see paper p.11)[1]
#define MERSENNETWISTER_N 624
#define MERSENNETWISTER_M 397

#define MERSENNETWISTER_A 0x9908b0df

#define MERSENNETWISTER_B 0x9d2c5680
#define MERSENNETWISTER_C 0xefc60000

#define MERSENNETWISTER_U 11
#define MERSENNETWISTER_S 7
#define MERSENNETWISTER_T 15
#define MERSENNETWISTER_L 18


//---------------------------------------------------------------------------------------------------------------//
// Mersenne Twister
struct _MersenneTwister{
  UInt index;
  UInt32 state[MERSENNETWISTER_N];
};


//---------------------------------------------------------------------------------------------------------------//
// Tuples
struct _Tuple_MersenneTwister_UInt32 {
  MersenneTwister first;
  UInt32 second;
};

struct _Tuple_MersenneTwister_Float32 {
  MersenneTwister first;
  Float32 second;
};

struct _Tuple_MersenneTwister_Float32_Float32 {
  MersenneTwister first;
  Float32 second;
  Float32 third;
};


//---------------------------------------------------------------------------------------------------------------//

#endif // RANDOM_DATA

