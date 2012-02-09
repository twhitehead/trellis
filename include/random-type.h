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


#ifndef RANDOM_TYPE
#define RANDOM_TYPE


//---------------------------------------------------------------------------------------------------------------//
// Single upper case to distinguish types

typedef struct _MersenneTwister* MersenneTwister;

typedef struct _MersenneTwister_UInt32 MersenneTwister_UInt32;
typedef struct _MersenneTwister_Float32 MersenneTwister_Float32;
typedef struct _MersenneTwister_Float32_Float32 MersenneTwister_Float32_Float32;


//---------------------------------------------------------------------------------------------------------------//

#endif // RANDOM_TYPE
