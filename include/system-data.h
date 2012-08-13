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

#include <stdio.h>

#include "system-type.h"

#include "model-data.h"

// Split 64b Z space into components (require CLUSTER=2^N for some N and CLUSTER^DEPTH = 2^64)

#define DEPTH  4
#define CLUSTER 65536

// Z values
//
// Z values form a fractal, which each set of two bits giving a layer of 2x2 connected boxes connected in a Z
// (upper left, upper right, lower left, lower right).  Each layer connects the 2x2 connected boxes from the lower
// layer into a new set of 2x2 connected boxes.  The lowest layer is formed by the two least significant bits.
//
// The X & Y values contained in the Z values can easily be compared by masking.

#define X_MASK 0x5555555555555555
#define Y_MASK 0xaaaaaaaaaaaaaaaa


//---------------------------------------------------------------------------------------------------------------//
// Space
struct _Space {
  Bool periodic_x;
  Bool periodic_y;

  Float32 scale;

  Float32 size_x;
  Float32 size_y;
};


//---------------------------------------------------------------------------------------------------------------//
// Individuals
struct _SVarieties {
  UInt number;
  Variety variety[];
};


struct _SIndividuals {
  UInt64 number;
  UInt64* z;
  Individual* individual;
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
  SIndividuals sindividuals;
  UInt64 index;
};


//---------------------------------------------------------------------------------------------------------------//
// Building (_AVarieties, _ASIndividuals, and _AIndividuals are same as their constant variants)


//---------------------------------------------------------------------------------------------------------------//
// Tuples
struct _AVarieties_ASIndividuals {
  AVarieties first;
  ASIndividuals second;
};

struct _World_SVarieties_SSIndividuals {
  World first;
  SVarieties second;
  SSIndividuals third;
};

struct _Space_World_SVarieties_SSIndividuals {
  Space first;
  World second;
  SVarieties third;
  SSIndividuals fourth;
};

struct _Space_World_SVarieties_SSIndividuals_FILE_UInt64 {
  Space first;
  World second;
  SVarieties third;
  SSIndividuals fourth;
  FILE* fifth;
  UInt64 sixth;
};

struct _RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut {
  RWorld first;
  SRVarieties second;
  SSRIndividualsIn third;
  SSRIndividualsOut fourth;
};


//---------------------------------------------------------------------------------------------------------------//

#endif // SYSTEM_DATA
