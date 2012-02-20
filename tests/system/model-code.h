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
World World_raw(UInt32 const id) {
  World const world = { .id = id };
  return world;
}


//
FILE* World_saveFP(Space const space, World const world,
                   char const* const name, FILE* file) {
  if ( fprintf(file, "%"PRIu32"\n",
               world.id) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


World_FILE_UInt64 World_loadFP(Space const space,
                               char const* const name, FILE* file, UInt64 const line) {
  World world;
  Int records;

  records = fscanf(file,"%"SCNu32"\n",
                   &world.id);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 1 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
              "ID", name, line);

  return pack_World_FILE_UInt64(world, file, line+1);
}


//
World World_next(Space const space, World const world, RWorld const rworld, Thread const thread) {
  return world;
}


//---------------------------------------------------------------------------------------------------------------//
//
Variety Variety_raw(UInt32 const id) {
  Variety const variety = { .id = id };
  return variety;
}


//
FILE* Variety_saveFP(Space const space, World const world, Variety const variety,
                     char const* const name, FILE* file) {
  if ( fprintf(file, "%"PRIu32"\n",
               variety.id) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


Variety_FILE_UInt64 Variety_loadFP(Space const space, World const world,
                                   char const* const name, FILE* file, UInt64 const line) {
  Variety variety;
  Int records;

  records = fscanf(file, "%"SCNu32"\n",
                   &variety.id);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 1 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
              "ID", name, line);

  return pack_Variety_FILE_UInt64(variety, file, line+1);
}


//
Variety Variety_next(Space const space, World const world, Variety const variety,
                     RWorld const rworld, RVariety const rvariety,
                     Thread const thread) {
  return variety;
}


//---------------------------------------------------------------------------------------------------------------//
//
Individual Individual_raw(UInt32 const id, Float32 const x, Float32 const y,
                          Float32 const in, Float32 const out) {
  Individual const individual = { .id = id, .x = x, .y = y, .in = in, .out = out };
  return individual;
}


//
FILE* Individual_saveFP(Space const space, World const world, Variety const variety, Individual const individual,
                        char const* const name, FILE* file) {
  if ( fprintf(file, "%"PRIu32" %g %g %g %g\n",
               individual.id, individual.x, individual.y, individual.in, individual.out) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


Individual_FILE_UInt64 Individual_loadFP(Space const space, World const world, Variety const variety,
                                               char const* const name, FILE* file, UInt64 const line) {
  Individual individual;
  Int records;

  records = fscanf(file, "%"SCNu32" %g %g %g %g\n",
                   &individual.id, &individual.x, &individual.y, &individual.in, &individual.out);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 5 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
              "ID X Y IN OUT", name, line);
  //  if (individual.x < 0 || individual.x >= space.size_x)
  //    Error_die(1, "problem parsing \"%s\":%"PRIu64": "
  //          "the constraint 0 <= X=%g < SIZE_X=%g does not hold", name, line,
  //          individual.x, space.size_x);
  //  if (individual.y < 0 || individual.y >= space.size_y)
  //    Error_die(1, "problem parsing \"%s\":%"PRIu64": "
  //          "the constraint 0 <= Y=%g < SIZE_Y=%g does not hold", name, line,
  //          individual.y, space.size_y);

  return pack_Individual_FILE_UInt64(individual, file, line+1);
}


//
AIndividuals Individual_next(AIndividuals const aindividuals, Space const space,
                             World const world, Variety const variety, Individual const individual,
                             RWorld const rworld, RVariety const rvariety,
                             RIndividualIn const rindividualin, RIndividualOut const rindividualout,
                             Thread const thread) {
  return AIndividuals_append(aindividuals, individual, Z_xy(individual.x, individual.y, space.scale));
}


//---------------------------------------------------------------------------------------------------------------//
//
RWorld RWorld_raw(UInt32 const ids) {
  RWorld const rworld = { .ids = ids };
  return rworld;
}


//
RWorld RWorld_first(Space const space, World const world) {
  return RWorld_raw(world.id);
}
RWorld RWorld_rest(Space const space, World const world, Variety const variety, Individual const individual) {
  return RWorld_raw(individual.id);
}
RWorld RWorld_merge(RWorld const rworld0, RWorld const rworld1) {
  return RWorld_raw(rworld0.ids ^ rworld1.ids);
}


//---------------------------------------------------------------------------------------------------------------//
//
RVariety RVariety_raw(UInt32 const ids) {
  RVariety const rvariety = { .ids = ids };
  return rvariety;
}


//
RVariety RVariety_first(Space const space, World const world, Variety const variety) {
  return RVariety_raw(variety.id);
}
RVariety RVariety_rest(Space const space, World const world, Variety const variety, Individual const individual) {
  return RVariety_raw(individual.id);
}
RVariety RVariety_merge(RVariety const rvariety0, RVariety const rvariety1) {
  return RVariety_raw(rvariety0.ids ^ rvariety1.ids);
}


//---------------------------------------------------------------------------------------------------------------//
//
RIndividualIn RIndividualIn_raw(UInt32 const ids) {
  RIndividualIn const rindividualin = { .ids = ids };
  return rindividualin;
}


//
Float32_Float32_Float32_Float32 RIndividualIn_bound(Space const space, World const world,
                                                    Variety const variety0, Individual const individual0,
                                                    Variety const variety1) {
  return pack_Float32_Float32_Float32_Float32(individual0.x-individual0.in, individual0.y-individual0.in,
                                              individual0.x+individual0.in, individual0.y+individual0.in);
}

Bool RIndividualIn_filter(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                          Variety const variety0, Individual const individual0,
                          Variety const variety1, Individual const individual1) {
  Float32 const delta_x = individual1.x+mirror_x*space.size_x-individual0.x;
  Float32 const delta_y = individual1.y+mirror_y*space.size_y-individual0.y;
  return delta_x*delta_x + delta_y*delta_y < individual0.in*individual0.in;
}


//
RIndividualIn RIndividualIn_first(Space space, World const world,
                                  Variety const variety0, Individual const individual0) {
  return RIndividualIn_raw(individual0.id);
}

RIndividualIn RIndividualIn_rest(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                                 Variety const variety0, Individual const individual0,
                                 Variety const variety1, Individual const individual1) {
  return RIndividualIn_raw(individual1.id);
}
RIndividualIn RIndividualIn_merge(RIndividualIn const rindividualin0, RIndividualIn const rindividualin1) {
  return RIndividualIn_raw(rindividualin0.ids ^ rindividualin1.ids);
}


//---------------------------------------------------------------------------------------------------------------//
//
RIndividualOut RIndividualOut_raw(UInt32 const ids) {
  RIndividualOut const rindividualout = { .ids = ids };
  return rindividualout;
}


//
Float32_Float32_Float32_Float32 RIndividualOut_bound(Space const space, World const world,
                                                     Variety const variety0, Individual const individual0,
                                                     Variety const variety1) {
  return pack_Float32_Float32_Float32_Float32(individual0.x-individual0.out, individual0.y-individual0.out,
                                              individual0.x+individual0.out, individual0.y+individual0.out);
}

Bool RIndividualOut_filter(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                           Variety const variety0, Individual const individual0,
                           Variety const variety1, Individual const individual1) {
  Float32 const delta_x = individual1.x+mirror_x*space.size_x-individual0.x;
  Float32 const delta_y = individual1.y+mirror_y*space.size_y-individual0.y;
  return delta_x*delta_x + delta_y*delta_y < individual0.out*individual0.out;
}


//
RIndividualOut RIndividualOut_first(Space const space, World const world,
                                    Variety const variety1, Individual const individual1) {
  return RIndividualOut_raw(individual1.id);
}
RIndividualOut RIndividualOut_rest(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                                   Variety const variety0, Individual const individual0,
                                   Variety const variety1, Individual const individual1) {
  return RIndividualOut_raw(individual0.id);
}
RIndividualOut RIndividualOut_merge(RIndividualOut const rindividualout0, RIndividualOut const rindividualout1) {
  return RIndividualOut_raw(rindividualout0.ids ^ rindividualout1.ids);
}


//---------------------------------------------------------------------------------------------------------------//
// Tuples
Float32_Float32_Float32_Float32 pack_Float32_Float32_Float32_Float32
(Float32 const first, Float32 const second, Float32 const third, Float32 const fourth) {
  Float32_Float32_Float32_Float32 const value =
    { .first = first, .second = second, .third = third, .fourth = fourth };
  return value;
}

World_FILE_UInt64 pack_World_FILE_UInt64
(World const first, FILE* const second, UInt64 const third) {
  World_FILE_UInt64 const value =
    { .first = first, .second = second, .third = third };
  return value;
}

Variety_FILE_UInt64 pack_Variety_FILE_UInt64
(Variety const first, FILE* const second, UInt64 const third) {
  Variety_FILE_UInt64 const value =
    { .first = first, .second = second, .third = third };
  return value;
}

Individual_FILE_UInt64 pack_Individual_FILE_UInt64
(Individual const first, FILE* const second, UInt64 const third) {
  Individual_FILE_UInt64 const value =
    { .first = first, .second = second, .third = third };
  return value;
}


//
void unpack_Float32_Float32_Float32_Float32(Float32* const first, Float32* const second,
                                            Float32* const third, Float32* const fourth,
                                            Float32_Float32_Float32_Float32 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
  *fourth = tuple.fourth;
}

void unpack_World_FILE_UInt64(World* const first, FILE** const second, UInt64* const third,
                              World_FILE_UInt64 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
}

void unpack_Variety_FILE_UInt64(Variety* const first, FILE** const second, UInt64* const third,
                                Variety_FILE_UInt64 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
}

void unpack_Individual_FILE_UInt64(Individual* const first, FILE** const second, UInt64* const third,
                                   Individual_FILE_UInt64 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
}


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_CODE
