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


//---------------------------------------------------------------------------------------------------------------//
// Thread
struct _Thread {
};


//---------------------------------------------------------------------------------------------------------------//
// Reduction
struct _RWorld {
  UInt32 ids;
};

struct _RVariety {
  UInt32 ids;
};


struct _RIndividualIn {
  UInt32 ids;
};


struct _RIndividualOut {
  UInt32 ids;
};


//---------------------------------------------------------------------------------------------------------------//
// Model
struct _World {
  UInt32 id;
};


struct _Variety {
  UInt32 id;
};


struct _Individual {
  UInt32 id;
  Float32 x;
  Float32 y;
  Float32 in;
  Float32 out;
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
