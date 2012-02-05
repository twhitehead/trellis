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

Tuple_MersenneTwister_UInt32 MersenneTwister_extract_UInt32(MersenneTwister mersennetwister);
void MersenneTwister_extract_(MersenneTwister* first, UInt32* second,
			      MersenneTwister mersennetwsiter);

//
Tuple_MersenneTwister_UInt32 Random_uniform_UInt32(MersenneTwister mersennetwister, UInt32 n);
Tuple_MersenneTwister_UInt32 Random_binomial_UInt32(MersenneTwister mersennetwister, UInt32 n, Float32 p);
Tuple_MersenneTwister_Float32 Random_uniform_Float32(MersenneTwister mersennetwister);
Tuple_MersenneTwister_Float32_Float32 Random_normal2_Float32(MersenneTwister mersennetwister);

void Random_uniform_UInt32_(MersenneTwister* first, UInt32* second,
			    MersenneTwister mersennetwister, UInt32 n);
void Random_binomial_UInt32_(MersenneTwister* first, UInt32* second,
			     MersenneTwister mersennetwister, UInt32 n, Float32 p);
void Random_uniform_Float32_(MersenneTwister* first, Float32* second,
			     MersenneTwister mersennetwister);
void Random_normal2_Float32_(MersenneTwister* first, Float32* second, Float32* third,
			     MersenneTwister mersennetwister);

//
Tuple_MersenneTwister_UInt32 tuple_MersenneTwister_UInt32(MersenneTwister first, UInt32 second);
Tuple_MersenneTwister_Float32 tuple_MersenneTwister_Float32(MersenneTwister first, Float32 second);
Tuple_MersenneTwister_Float32_Float32 tuple_MersenneTwister_Float32_Float32
(MersenneTwister first, Float32 second, Float32 third);

Tuple_World_SVarieties_SSIndividuals tuple_World_SVarieties_SSIndividuals
(World first, SVarieties second, SSIndividuals third);
Tuple_Space_World_SVarieties_SSIndividuals tuple_Space_World_SVarieties_SSIndividuals
(Space first, World second, SVarieties third, SSIndividuals fourth);
Tuple_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
tuple_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
(RWorld first, SRVarieties second, SSRIndividualsIn third, SSRIndividualsOut fourth);


//---------------------------------------------------------------------------------------------------------------//

#endif // RANDOM_FUNC
