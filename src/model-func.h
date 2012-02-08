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
Tuple_World_FILE_UInt64 World_loadFP(Space space, char const* name, FILE* file, UInt64 line);
void World_loadFP_(World* first, FILE** second, UInt64* third,
		   Space space, char const* name, FILE* file, UInt64 line);

FILE* Variety_saveFP(Space space, World world, Variety variety, char const* name, FILE* file);
Tuple_Variety_FILE_UInt64 Variety_loadFP(Space space, World world,
					 char const* name, FILE* file, UInt64 line);
void Variety_loadFP_(Variety* first, FILE** second, UInt64* third,
		     Space space, World world, char const* name, FILE* file, UInt64 line);

FILE* Individual_saveFP(Space space, World world, Variety variety, Individual individual,
			char const* name, FILE* file);
Tuple_Individual_FILE_UInt64 Individual_loadFP(Space space, World const world, Variety variety,
					       char const* name, FILE* file, UInt64 line);
void Individual_loadFP_(Individual* first, FILE** second, UInt64* third,
			Space space, World world, Variety variety, char const* name, FILE* file, UInt64 line);

//
RWorld RWorld_first(World world);
RWorld RWorld_rest(World world, Variety variety, Individual individual);
RWorld RWorld_merge(RWorld rworld0, RWorld rworld1);

//
RVariety RVariety_first(World world, Variety variety);
RVariety RVariety_rest(World world, Variety variety, Individual individual);
RVariety RVariety_merge(RVariety rvariety0, RVariety rvariety1);

//
Tuple_Float32_Float32_Float32_Float32
RIndividualIn_bound(World world, Variety variety0, Individual individual0, Variety variety1);
void RIndividualIn_bound_(Float32* first, Float32* second, Float32* third, Float32* fourth,
			  World world, Variety variety0, Individual individual0, Variety variety1);
Bool RIndividualIn_filter(World world, Variety variety0, Individual individual0,
			  Variety variety1, Individual individual1);

RIndividualIn RIndividualIn_first(World world, Variety variety, Individual individual);
RIndividualIn RIndividualIn_rest(World world, Variety variety0, Individual individual0,
				 Variety variety1, Individual individual1);
RIndividualIn RIndividualIn_merge(RIndividualIn rindividualin0, RIndividualIn rindividualin1);

//
Tuple_Float32_Float32_Float32_Float32
RIndividualOut_bound(World world, Variety variety0, Individual individual0, Variety variety1);
void RIndividualOut_bound_(Float32* first, Float32* second, Float32* third, Float32* fourth,
			   World world, Variety variety0, Individual individual0, Variety variety1);
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
Tuple_World_FILE_UInt64 tuple_World_FILE_UInt64(World first, FILE* second, UInt64 third);
Tuple_Variety_FILE_UInt64 tuple_Variety_FILE_UInt64(Variety first, FILE* second, UInt64 third);
Tuple_Individual_FILE_UInt64 tuple_Individual_FILE_UInt64(Individual first, FILE* second, UInt64 third);

Tuple_Float32_Float32_Float32_Float32 tuple_Float32_Float32_Float32_Float32
(Float32 first, Float32 second, Float32 third, Float32 fourth);


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_FUNC
