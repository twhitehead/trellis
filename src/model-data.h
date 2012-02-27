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


#ifndef MODEL_DATA
#define MODEL_DATA

#include <stdio.h>

#include "model-type.h"

#include "system-type.h"
#include "system-data.h"

#include "random-type.h"
#include "random-data.h"


//---------------------------------------------------------------------------------------------------------------//
// Thread
struct _Thread {
  MersenneTwister mersennetwister;
};


//---------------------------------------------------------------------------------------------------------------//
// Reduction
struct _RWorld {
};

struct _RVariety {
  UInt64 number_adult;
  UInt64 number_seedling;
};


struct _RIndividualIn {
};


struct _RIndividualOut {
  UInt64 number_lower;
  UInt64 number_higher;
};


//---------------------------------------------------------------------------------------------------------------//
// Model
struct _World {
  UInt year;
  Float32 cell_diameter;
};


//
enum _Variety_Type {
  Variety_Type_PaperBirch,
  Variety_Type_WhitePine,
  Variety_Type_SugarMaple,
  VARIETY_TYPE_INVALID
};


struct _Variety {
  Variety_Type type;

  Float32 height_mature;
  Float32 height_maximum;

  Float32 growth_rate;
  Float32 growth_decay;

  Float32 mortality_intrinsic;
  Float32 mortality_initial;
  Float32 mortality_decay;

  Float32 fecundity_maximum;

  Float32 dispersal_method;
  Float32 dispersal_short;
  Float32 dispersal_long;
};


//
struct _Individual {
  Float32 x;
  Float32 y;
  Float32 height;
};


//---------------------------------------------------------------------------------------------------------------//
// Tuples
struct _World_FILE_UInt64 {
  World first;
  FILE* second;
  UInt64 third;
};

struct _Variety_FILE_UInt64 {
  Variety first;
  FILE* second;
  UInt64 third;
};

struct _Individual_FILE_UInt64 {
  Individual first;
  FILE* second;
  UInt64 third;
};

struct _Float32_Float32_Float32_Float32 {
  Float32 first;
  Float32 second;
  Float32 third;
  Float32 fourth;
};


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_DATA
