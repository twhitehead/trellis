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


#ifndef RANDOM_FUNC
#define RANDOM_FUNC

#include "random-type.h"

#include "system-type.h"


//---------------------------------------------------------------------------------------------------------------//
//
MersenneTwister MersenneTwister_begin();
void MersenneTwister_end(MersenneTwister mersennetwister);

MersenneTwister MersenneTwister_urandom(MersenneTwister mersennetwister);
MersenneTwister MersenneTwister_seed(MersenneTwister mersennetwister, UInt32 seed);

MersenneTwister MersenneTwister_next(MersenneTwister mersennetwister);

MersenneTwister_UInt32 MersenneTwister_extract_UInt32(MersenneTwister mersennetwister);

//
MersenneTwister_UInt32 Random_uniform_UInt32(MersenneTwister mersennetwister, UInt32 n);
MersenneTwister_UInt32 Random_binomial_UInt32(MersenneTwister mersennetwister, UInt32 n, Float32 p);
MersenneTwister_Float32 Random_uniform_Float32(MersenneTwister mersennetwister);
MersenneTwister_Float32_Float32 Random_normal2_Float32(MersenneTwister mersennetwister);

//
MersenneTwister_UInt32 pack_MersenneTwister_UInt32(MersenneTwister first, UInt32 second);
MersenneTwister_Float32 pack_MersenneTwister_Float32(MersenneTwister first, Float32 second);
MersenneTwister_Float32_Float32 pack_MersenneTwister_Float32_Float32
(MersenneTwister first, Float32 second, Float32 third);

void unpack_MersenneTwister_UInt32(MersenneTwister* first, UInt32* second,
                                   MersenneTwister_UInt32 tuple);
void unpack_MersenneTwister_Float32(MersenneTwister* first, Float32* second,
                                    MersenneTwister_Float32 tuple);
void unpack_MersenneTwister_Float32_Float32(MersenneTwister* first, Float32* second, Float32* third,
                                            MersenneTwister_Float32_Float32 tuple);


//---------------------------------------------------------------------------------------------------------------//

#endif // RANDOM_FUNC
