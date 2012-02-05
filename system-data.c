//---------------------------------------------------------------------------------------------------------------//
// Copyright (C) 2011-12 Tyson Whitehead
//
// This code is free software; you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation; either version 2, or (at your option) any later
// version.
//
//---------------------------------------------------------------------------------------------------------------//

#if _POSIX_C_SOURCE < 200112
#error _POSIX_C_SOURCE 200112 or greater required
#endif


#ifndef SYSTEM_DATA
#define SYSTEM_DATA

#include "system-type.h"

#include "model-data.c"

// Split 64b Z space into components (require CLUSTER=2^N for some N and CLUSTER^DEPTH = 2^64)

#define DEPTH  4
#define CLUSTER 65536

// Buffer size for fatal errors

#define ERRNO_BUFFER 1024

// Z values
//
// Z values form a fractal, which each set of two bits giving a layer of 2x2 connected boxes connected in a Z
// (upper left, upper right, lower left, lower right).  Each layer connects the 2x2 connected boxes from the lower
// layer into a new set of 2x2 connected boxes.  The lowest layer is formed by the two least significant bits.
//
// The X & Y values contained in the Z values can easily be compared by masking.

#define X_MASK 0x5555555555555555u
#define Y_MASK 0xaaaaaaaaaaaaaaaau

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

// Mersenne Twister MT19937 parameters (see paper p.11)
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
// Individuals
struct _SVarieties {
  UInt number;
  Variety variety[];
};


union _SIndividuals_ {
  SIndividuals0 sindividuals0;
  SIndividuals1 sindividuals1;
};


struct _SIndividuals0 {
  UInt64 left_z[CLUSTER];
  UInt64 right_z[CLUSTER];
  SIndividuals_ sindividuals_[CLUSTER];
};


struct _SIndividuals1 {
  UInt64 z[CLUSTER];
  Individual individual[CLUSTER];
};


struct _SIndividuals {
  UInt64 number;
  SIndividuals_ sindividuals_;
};


struct _SSIndividuals {
  UInt number;
  SIndividuals sindividuals[];
};


//---------------------------------------------------------------------------------------------------------------//
// Reduction
struct _SRVarieties {
  UInt number;
  RVariety rvariety[];
};


struct _SSRIndividualsIn {
  UInt number;
  SRIndividualsIn srindividualsin[];
};


struct _SSRIndividualsOut {
  UInt number;
  SRIndividualsOut srindividualsout[];
};


struct _SRIndividualsIn {
  UInt64 number;
  RIndividualIn rindividualin[];
};


struct _SRIndividualsOut {
  UInt64 number;
  RIndividualOut rindividualout[];

};


//---------------------------------------------------------------------------------------------------------------//
// Iterating
struct _IZ {
  Bool valid;
  UInt64 z;
};


struct _IIndividuals {
  Bool valid;

  UInt64 number;
  UInt64 index;

  SIndividuals_ sindividuals_;
  SIndividuals1 sindividuals1;
};


//---------------------------------------------------------------------------------------------------------------//
// Building (_AVarieties and _ASIndividuals are same as _SVarieties and _SSIndividuals)
struct _AIndividuals {
  UInt64 number;

  UInt64 z_min;
  UInt64 z_max;

  SIndividuals_ sindividuals_;
  SIndividuals1 sindividuals1;
};


//---------------------------------------------------------------------------------------------------------------//
// Returned
struct _Box {
  Float32 ul_x;
  Float32 ul_y;
  Float32 lr_x;
  Float32 lr_y;
};


struct _Thread {
};


struct _Space {
  Bool periodic_x;
  Bool periodic_y;

  Float32 scale;

  Float32 size_x;
  Float32 size_y;

  World world;
  SVarieties svarieties;
  SSIndividuals ssindividuals;
};


struct _RSpace {
  RWorld rworld;
  SRVarieties srvarieties;
  SSRIndividualsIn ssrindividualsin;
  SSRIndividualsOut ssrindividualsout;
};


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

struct _Tuple_World_SVarieties_SSIndividuals {
  World first;
  SVarieties second;
  SSIndividuals third;
};

struct _Tuple_Space_World_SVarieties_SSIndividuals {
  Space first;
  World second;
  SVarieties third;
  SSIndividuals fourth;
};

struct _Tuple_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut {
  RWorld first;
  SRVarieties second;
  SSRIndividualsIn third;
  SSRIndividualsOut fourth;
};


//---------------------------------------------------------------------------------------------------------------//

#endif // SYSTEM_DATA
