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


#ifndef MODEL_C_CODE
#define MODEL_C_CODE

#include "model-type.h"
#include "model-func.h"
#include "model-data.h"

#include "system-type.h"
#include "system-func.h"
#include "system-data.h"
#include "system-code.h"


//---------------------------------------------------------------------------------------------------------------//
//
FILE* World_saveFP(Space const space, World const world,
		   char const* const name, FILE* file) {
  if ( fprintf(file, "%"PRIu32"\n",
	       world.year) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


void World_loadFP_(World* const first, FILE** const second, UInt64* const third,
		   Space const space,
		   char const* const name, FILE* const file, UInt64 const line) {
  Tuple_World_FILE_UInt64 const value =
    World_loadFP(space, name, file, line);
  *first = value.first;
  *second = value.second;
  *third = value.third;
}
Tuple_World_FILE_UInt64 World_loadFP(Space const space,
				     char const* const name, FILE* file, UInt64 const line) {
  World world;
  Int records;

  records = fscanf(file,"%"SCNu32"\n",
                   &world.year);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 1 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
	      "YEAR", name, line);

  return tuple_World_FILE_UInt64(world, file, line+1);
}


//---------------------------------------------------------------------------------------------------------------//
//
FILE* Variety_saveFP(Space const space, World const world, Variety const variety,
		     char const* const name, FILE* file) {
  if ( fprintf(file, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	       variety.height_mature, variety.height_maximum,
	       variety.growth_rate, variety.growth_competition_lower, variety.growth_competition_higher,
	       variety.mortality_initial, variety.mortality_decay, variety.mortality_intrinsic,
	       variety.fecundity_maximum,
	       variety.masting_time, variety.masting_phase,
	       variety.dispersal_probability_short, variety.dispersal_mode_short, variety.dispersal_mode_long)
       < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


void Variety_loadFP_(Variety* const first, FILE** const second, UInt64* const third,
		     Space const space, World const world,
		     char const* const name, FILE* const file, UInt64 const line) {
  Tuple_Variety_FILE_UInt64 const value =
    Variety_loadFP(space, world, name, file, line);
  *first = value.first;
  *second = value.second;
  *third = value.third;
}
Tuple_Variety_FILE_UInt64 Variety_loadFP(Space const space, World const world,
					 char const* const name, FILE* file, UInt64 const line) {
  Variety variety;
  Int records;

  records = fscanf(file, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
		   &variety.height_mature, &variety.height_maximum,
		   &variety.growth_rate, &variety.growth_competition_lower, &variety.growth_competition_higher,
		   &variety.mortality_initial, &variety.mortality_decay, &variety.mortality_intrinsic,
		   &variety.fecundity_maximum,
		   &variety.masting_time, &variety.masting_phase,
		   &variety.dispersal_probability_short,
		   &variety.dispersal_mode_short, &variety.dispersal_mode_long);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 14 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
	      "HEIGHT_MATURE HEIGHT_MAXIMUM "
	      "GROWTH_RATE VARIETY_GROWTH_COMPETITION_LOWER VARIETY_GROWTH_COMPETITION_HIGHER "
	      "MORTALITY_INITIAL MORTALITY_DECAY MORTALITY_INTRINSIC "
	      "FECUNDITY_MAXIMUM "
	      "MASTING_TIME MASTING_PHASE "
	      "DISPERSAL_PROBABILITY_SHORT DISPERSAL_MODE_SHORT DISPERSAL_MODE_LONG", name, line);

  return tuple_Variety_FILE_UInt64(variety, file, line+1);
}


//---------------------------------------------------------------------------------------------------------------//
//
FILE* Individual_saveFP(Space const space, World const world, Variety const variety, Individual const individual,
			char const* const name, FILE* file) {
  if ( fprintf(file, "%g %g %g\n",
	       individual.x, individual.y,
	       individual.height) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


void Individual_loadFP_(Individual* const first, FILE** const second, UInt64* const third,
			Space const space, World const world, Variety const variety,
			char const* const name, FILE* const file, UInt64 const line) {
  Tuple_Individual_FILE_UInt64 const value =
    Individual_loadFP(space, world, variety, name, file, line);
  *first = value.first;
  *second = value.second;
  *third = value.third;
}
Tuple_Individual_FILE_UInt64 Individual_loadFP(Space const space, World const world, Variety const variety,
					       char const* const name, FILE* file, UInt64 const line) {
  Individual individual;
  Int records;

  records = fscanf(file, "%g %g %g\n",
		   &individual.x, &individual.y, &individual.height);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 3 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
	      "X Y "
	      "HEIGHT", name, line);
  if (individual.x < 0 || individual.x >= space.size_x)
    Error_die(1, "problem parsing \"%s\":%"PRIu64": "
	      "the constraint 0 <= X=%g < SIZE_X=%g does not hold", name, line,
	      individual.x, space.size_x);
  if (individual.y < 0 || individual.y >= space.size_y)
    Error_die(1, "problem parsing \"%s\":%"PRIu64": "
	      "the constraint 0 <= Y=%g < SIZE_Y=%g does not hold", name, line,
	      individual.y, space.size_y);
  if (individual.height < 0)
    Error_die(1, "problem parsing \"%s\":%"PRIu64": "
	      "the constraint 0 <= HEIGHT=%g", name, line, individual.height);

  return tuple_Individual_FILE_UInt64(individual, file, line+1);
}




//---------------------------------------------------------------------------------------------------------------//
//
RWorld RWorld_first(World const world) {
  RWorld rworld = { };
  return rworld;
}
RWorld RWorld_rest(World const world, Variety const variety, Individual const individual) {
  RWorld rworld = { };
  return rworld;
}
RWorld RWorld_merge(RWorld const rworld0, RWorld const rworld1) {
  RWorld rworld = { };
  return rworld;
}


RVariety RVariety_first(World const world, Variety const variety) {
  RVariety rvariety = { };
  return rvariety;
}
RVariety RVariety_rest(World const world, Variety const variety, Individual const individual) {
  RVariety rvariety = { };
  return rvariety;
}
RVariety RVariety_merge(RVariety const rvariety0, RVariety const rvariety1) {
  RVariety rvariety = { };
  return rvariety;
}


void RIndividualIn_bound_(Float32* const first, Float32* const second, Float32* const third,
			  Float32* const fourth,
			  World const world, Variety const variety0, Individual const individual0,
			  Variety const variety1) {
  Tuple_Float32_Float32_Float32_Float32 value =
    RIndividualIn_bound(world, variety0, individual0, variety1);
  *first = value.first;
  *second = value.second;
  *third = value.third;
  *fourth = value.fourth;
}
Tuple_Float32_Float32_Float32_Float32 RIndividualIn_bound(World const world,
							  Variety const variety0, Individual const individual0,
							  Variety const variety1) {
  Tuple_Float32_Float32_Float32_Float32 box = { };
  return box;
}
void RIndividualOut_bound_(Float32* const first, Float32* const second, Float32* const third,
			   Float32* const fourth,
			   World const world, Variety const variety0, Individual const individual1, 
			   Variety const variety1) {
  Tuple_Float32_Float32_Float32_Float32 value =
    RIndividualOut_bound(world, variety0, individual1, variety1);
  *first = value.first;
  *second = value.second;
  *third = value.third;
  *fourth = value.fourth;
}
Tuple_Float32_Float32_Float32_Float32 RIndividualOut_bound(World const world,
							   Variety const variety0, Individual const individual1, 
							   Variety const variety1) {
  Tuple_Float32_Float32_Float32_Float32 box = { };
  return box;
}


Bool RIndividualIn_filter(World const world, Variety const variety0, Individual const individual0,
			  Variety const variety1, Individual const individual1) {
  return true;
}

Bool RIndividualOut_filter(World const world, Variety const variety0, Individual const individual0,
			   Variety const variety1, Individual const individual1) {
  return true;
}


RIndividualIn RIndividualIn_first(World const world, Variety const variety, Individual const individual) {
  RIndividualIn rindividualin = { };
  return rindividualin;
}
RIndividualIn RIndividualIn_rest(World const world, Variety const variety0, Individual const individual0,
				 Variety const variety1, Individual const individual1) {
  RIndividualIn rindividualin = { };
  return rindividualin;
}
RIndividualIn RIndividualIn_merge(RIndividualIn const rindividualin0, RIndividualIn const rindividualin1) {
  RIndividualIn rindividualin = { };
  return rindividualin;
}


RIndividualOut RIndividualOut_first(World const world, Variety const variety, Individual const individual) {
  RIndividualOut rindividualin = { };
  return rindividualin;
}
RIndividualOut RIndividualOut_rest(World const world, Variety const variety0, Individual const individual0,
				   Variety const variety1, Individual const individual1) {
  RIndividualOut rindividualin = { };
  return rindividualin;
}
RIndividualOut RIndividualOut_merge(RIndividualOut const rindividualin0, RIndividualOut const rindividualin1) {
  RIndividualOut rindividualin = { };
  return rindividualin;
}


//---------------------------------------------------------------------------------------------------------------//
//
World World_next(Space const space, World const world, RWorld const rworld, Thread const thread) {
  return world;
}

Variety Variety_next(Space const space, World const world, Variety const variety,
		     RWorld const rworld, RVariety const rvariety,
		     Thread const thread) {
  return variety;
}

AIndividuals Individual_next(AIndividuals const aindividuals, Space const space,
                             World const world, Variety const variety, Individual const individual,
                             RWorld const rworld, RVariety const rvariety,
                             RIndividualIn const rindividualin, RIndividualOut const rindividualout,
                             Thread const thread) {
  return AIndividuals_append(aindividuals, individual, Z_xy(individual.x, individual.y, space.scale));
}


//---------------------------------------------------------------------------------------------------------------//
// Tuples
Tuple_Float32_Float32_Float32_Float32 tuple_Float32_Float32_Float32_Float32
(Float32 const first, Float32 const second, Float32 const third, Float32 const fourth) {
  Tuple_Float32_Float32_Float32_Float32 const value =
    { .first = first, .second = second, .third = third, .fourth = fourth };
  return value;
}

Tuple_World_FILE_UInt64 tuple_World_FILE_UInt64
(World const first, FILE* const second, UInt64 const third) {
  Tuple_World_FILE_UInt64 const value =
    { .first = first, .second = second, .third = third };
  return value;
}

Tuple_Variety_FILE_UInt64 tuple_Variety_FILE_UInt64
(Variety const first, FILE* const second, UInt64 const third) {
  Tuple_Variety_FILE_UInt64 const value =
    { .first = first, .second = second, .third = third };
  return value;
}

Tuple_Individual_FILE_UInt64 tuple_Individual_FILE_UInt64
(Individual const first, FILE* const second, UInt64 const third) {
  Tuple_Individual_FILE_UInt64 const value =
    { .first = first, .second = second, .third = third };
  return value;
}


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_CODE
