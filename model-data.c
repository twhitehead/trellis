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

#include "model-type.h"

#include "system-type.h"


//---------------------------------------------------------------------------------------------------------------//
// Reduction
struct _RWorld {
};

struct _RVariety {
};


struct _RIndividualIn {
};


struct _RIndividualOut {
};


//---------------------------------------------------------------------------------------------------------------//
// Model
struct _World {
  UInt32 year;
};


struct _Variety {
  Float32 height_mature;
  Float32 height_maximum;

  Float32 growth_rate;
  Float32 growth_competition_lower;
  Float32 growth_competition_higher;

  Float32 mortality_initial;
  Float32 mortality_decay;
  Float32 mortality_intrinsic;

  Float32 fecundity_maximum;

  Float32 masting_time;
  Float32 masting_phase;

  Float32 dispersal_probability_short;
  Float32 dispersal_mode_short;
  Float32 dispersal_mode_long;
};


struct _Individual {
  Float32 x;
  Float32 y;
  Float32 height;
};


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_DATA
