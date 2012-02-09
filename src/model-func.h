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


#ifndef MODEL_FUNC
#define MODEL_FUNC

#include <stdio.h>

#include "model-type.h"

#include "system-type.h"


//---------------------------------------------------------------------------------------------------------------//
//
FILE* World_saveFP(Space space, World world, char const* name, FILE* file);
World_FILE_UInt64 World_loadFP(Space space, char const* name, FILE* file, UInt64 line);

FILE* Variety_saveFP(Space space, World world, Variety variety, char const* name, FILE* file);
Variety_FILE_UInt64 Variety_loadFP(Space space, World world,
                                   char const* name, FILE* file, UInt64 line);

FILE* Individual_saveFP(Space space, World world, Variety variety, Individual individual,
			char const* name, FILE* file);
Individual_FILE_UInt64 Individual_loadFP(Space space, World const world, Variety variety,
                                         char const* name, FILE* file, UInt64 line);

//
RWorld RWorld_first(World world);
RWorld RWorld_rest(World world, Variety variety, Individual individual);
RWorld RWorld_merge(RWorld rworld0, RWorld rworld1);

//
RVariety RVariety_first(World world, Variety variety);
RVariety RVariety_rest(World world, Variety variety, Individual individual);
RVariety RVariety_merge(RVariety rvariety0, RVariety rvariety1);

//
Float32_Float32_Float32_Float32
RIndividualIn_bound(World world, Variety variety0, Individual individual0, Variety variety1);
Bool RIndividualIn_filter(World world, Variety variety0, Individual individual0,
			  Variety variety1, Individual individual1);

RIndividualIn RIndividualIn_first(World world, Variety variety, Individual individual);
RIndividualIn RIndividualIn_rest(World world, Variety variety0, Individual individual0,
				 Variety variety1, Individual individual1);
RIndividualIn RIndividualIn_merge(RIndividualIn rindividualin0, RIndividualIn rindividualin1);

//
Float32_Float32_Float32_Float32
RIndividualOut_bound(World world, Variety variety0, Individual individual0, Variety variety1);
Bool RIndividualOut_filter(World world, Variety variety0, Individual individual0,
			   Variety variety1, Individual individual1);

RIndividualOut RIndividualOut_first(World world, Variety variety, Individual individual);
RIndividualOut RIndividualOut_rest(World world, Variety variety0, Individual individual0,
				   Variety variety1, Individual individual1);
RIndividualOut RIndividualOut_merge(RIndividualOut rindividualin0, RIndividualOut rindividualin1);

//
World World_next(Space space, World world, RWorld rworld, Thread thread);

Variety Variety_next(Space space, World world, Variety variety, RWorld rworld, RVariety rvariety, Thread thread);

AIndividuals Individual_next(AIndividuals aindividuals, Space space,
			     World world, Variety variety, Individual individual,
			     RWorld rworld, RVariety rvariety,
			     RIndividualIn rindividualin, RIndividualOut rindividualout,
			     Thread thread);

//
World_FILE_UInt64 pack_World_FILE_UInt64(World first, FILE* second, UInt64 third);
Variety_FILE_UInt64 pack_Variety_FILE_UInt64(Variety first, FILE* second, UInt64 third);
Individual_FILE_UInt64 pack_Individual_FILE_UInt64(Individual first, FILE* second, UInt64 third);

Float32_Float32_Float32_Float32 pack_Float32_Float32_Float32_Float32
(Float32 first, Float32 second, Float32 third, Float32 fourth);

void unpack_World_FILE_UInt64(World* first, FILE** second, UInt64* third, World_FILE_UInt64 tuple);
void unpack_Variety_FILE_UInt64(Variety* first, FILE** second, UInt64* third, Variety_FILE_UInt64 tuple);
void unpack_Individual_FILE_UInt64(Individual* first, FILE** second, UInt64* third,
                                   Individual_FILE_UInt64 tuple);

void unpack_Float32_Float32_Float32_Float32(Float32* first, Float32* second, Float32* third, Float32* fourth,
                                            Float32_Float32_Float32_Float32 tuple);


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_FUNC
