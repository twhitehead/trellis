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

#include <math.h>
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#include "model-type.h"
#include "model-func.h"
#include "model-data.h"

#include "system-type.h"
#include "system-func.h"
#include "system-data.h"
#include "system-code.h"

#include "random-type.h"
#include "system-func.h"
#include "random-data.h"
#include "random-code.h"


//---------------------------------------------------------------------------------------------------------------//
//
Thread Thread_begin() {
  Thread thread;
  if ( !(thread = malloc(sizeof(struct _Thread)) ) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for Thread",
                   sizeof(struct _Thread));

  thread->mersennetwister = MersenneTwister_begin();

  return thread;
}


void Thread_end(Thread const thread) {
  MersenneTwister_end(thread->mersennetwister);

  free(thread);
}


//---------------------------------------------------------------------------------------------------------------//
//
World World_raw(Float32 const cell_diameter, UInt const year,
                UInt64 const number_seedling, UInt64 const number_adult) {
  World const world = { .cell_diameter = cell_diameter, .year = year,
                        .number_seedling = number_seedling, .number_adult = number_adult };
  return world;
}


//
FILE* World_saveFP(Space const space, World const world,
                   char const* const name, FILE* file) {
  if ( fprintf(file, "%g %u %"PRIu64" %"PRIu64"\n",
               world.cell_diameter, world.year, world.number_seedling, world.number_adult) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


World_FILE_UInt64 World_loadFP(Space const space,
                               char const* const name, FILE* file, UInt64 const line) {
  World world;
  Int records;

  records = fscanf(file,"%g %u %"SCNu64" %"SCNu64"\n",
                   &world.cell_diameter, &world.year, &world.number_seedling, &world.number_adult);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 4 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
              "CELL_DIAMETER "
              "YEAR "
              "NUMBER_SEEDLING NUMBER_ADULT", name, line);

  return pack_World_FILE_UInt64(world, file, line+1);
}


//
World_SVarieties_SSIndividuals World_next(Space const space,
                                          World world, AVarieties avarieties, ASIndividuals asindividuals, 
                                          RWorld const rworld, Thread const thread) {
  // Generate initial variety
  if (world.year == 0) {
    // Initial variety
    Variety const variety = Variety_raw(Variety_TypeOriginal,
                                        2.0, 25.0,
                                        0.1, 30.0, 30.0,
                                        0.01, 0.7, 80.0,
                                        5.0,
                                        10.0, 0.0,
                                        0.95, 10.0, 500.0);
    avarieties = AVarieties_append(avarieties, variety);

    // Initial variety individuals
    AIndividuals aindividuals = AIndividuals_begin();

    UInt const individual_number = space.size_x*space.size_y/(world.cell_diameter*world.cell_diameter);
    for (UInt individual_iterator=0; individual_iterator < individual_number; individual_iterator+=1) {
      Float32 x, y;
      unpack_MersenneTwister_Float32( &thread->mersennetwister, &x,
                                      Random_uniform_Float32(thread->mersennetwister) );
      unpack_MersenneTwister_Float32( &thread->mersennetwister, &y,
                                      Random_uniform_Float32(thread->mersennetwister) );
      x *= space.size_x;
      y *= space.size_y;

      Individual const individual = Individual_raw(x, y, 0.02);

      aindividuals = AIndividuals_append(aindividuals, individual, Z_xy(x, y, space.scale));
    }

    SIndividuals const sindividuals = AIndividuals_end(aindividuals);
    asindividuals = ASIndividuals_append(asindividuals, sindividuals);
  }
  // Generate invader variety
  else if (world.year == 999) {
    // Invader variety
    Variety const variety = Variety_raw(Variety_TypeInvader,
                                        2.0, 25.0,
                                        0.1, 30.0, 30.0,
                                        0.01, 0.7, 80.0,
                                        5.0,
                                        10.0, 0.0,
                                        0.95, 10.0, 500.0);
    avarieties = AVarieties_append(avarieties, variety);

    // Invader variety individuals
    SIndividuals const sindividuals = AIndividuals_end(AIndividuals_begin());
    asindividuals = ASIndividuals_append(asindividuals, sindividuals);
  }

  // Next world
  world.year += 1;
  world.number_seedling += rworld.number_seedling;
  world.number_adult += rworld.number_adult;

  // Package and return
  SVarieties const svarieties = AVarieties_end(avarieties);
  SSIndividuals const ssindividuals = ASIndividuals_end(asindividuals);
  return pack_World_SVarieties_SSIndividuals(world, svarieties, ssindividuals);
}


//---------------------------------------------------------------------------------------------------------------//
//
Variety Variety_raw(Variety_Type const type,
                    Float32 const height_mature, Float32 const height_maximum,
                    Float32 const growth_rate,
                    Float32 const competition_lower, Float32 const competition_higher,
                    Float32 const mortality_intrinsic,
                    Float32 const mortality_initial, Float32 const mortality_decay,
                    Float32 const fecundity_maximum,
                    Float32 const masting_time, Float32 const masting_phase,
                    Float32 const dispersal_method, Float32 const dispersal_short,
                    Float32 const dispersal_long) {
  Variety const variety = { .type = type,
                            .height_mature = height_mature, .height_maximum = height_maximum,
                            .growth_rate = growth_rate,
                            .competition_lower = competition_lower, .competition_higher = competition_higher, 
                            .mortality_intrinsic = mortality_intrinsic,
                            .mortality_initial = mortality_initial, .mortality_decay = mortality_decay,
                            .fecundity_maximum = fecundity_maximum,
                            .masting_time = masting_time, .masting_phase = masting_phase,
                            .dispersal_method = dispersal_method,
                            .dispersal_short = dispersal_short, .dispersal_long = dispersal_long };
  return variety;
}


//
FILE* Variety_saveFP(Space const space, World const world, Variety const variety,
                     char const* const name, FILE* file) {
  if ( fprintf(file, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
               variety.type,
               variety.height_mature, variety.height_maximum,
               variety.growth_rate, variety.competition_lower, variety.competition_higher,
               variety.mortality_intrinsic, variety.mortality_initial, variety.mortality_decay,
               variety.fecundity_maximum,
               variety.masting_time, variety.masting_phase,
               variety.dispersal_method, variety.dispersal_short, variety.dispersal_long)
       < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


Variety_FILE_UInt64 Variety_loadFP(Space const space, World const world,
                                         char const* const name, FILE* file, UInt64 const line) {
  Int type;
  Variety variety;
  Int records;

  records = fscanf(file, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                   &type,
                   &variety.height_mature, &variety.height_maximum,
                   &variety.growth_rate, &variety.competition_lower, &variety.competition_higher,
                   &variety.mortality_intrinsic, &variety.mortality_initial, &variety.mortality_decay,
                   &variety.fecundity_maximum,
                   &variety.masting_time, &variety.masting_phase,
                   &variety.dispersal_method,
                   &variety.dispersal_short, &variety.dispersal_long);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 15 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
              "TYPE "
              "HEIGHT_MATURE HEIGHT_MAXIMUM "
              "GROWTH_RATE "
              "COMPETITION_LOWER COMPETITION_HIGHER "
              "MORTALITY_INTRINSIC MORTALITY_INITIAL MORTALITY_DECAY "
              "FECUNDITY_MAXIMUM "
              "MASTING_TIME MASTING_PHASE "
              "DISPERSAL_METHOD DISPERSAL_SHORT DISPERSAL_LONG", name, line);
  variety.type = type;

  if (variety.type < 0 || variety.type >= VARIETY_TYPE_INVALID)
    Error_die(1, "problem parsing \"%s\":%"PRIu64" "
              "the constraint 0 <= TYPE=%d < VARIETY_TYPE_INVALID=%d does not hold", name, line,
              variety.type, VARIETY_TYPE_INVALID);

  return pack_Variety_FILE_UInt64(variety, file, line+1);
}


//
AVarieties_ASIndividuals Variety_next(AVarieties avarieties, ASIndividuals asindividuals, Space const space,
                                      World const world, Variety const variety, AIndividuals aindividuals,
                                      RWorld const rworld, RVariety const rvariety, Thread const thread) {
  // Seed invaders (from time of variety introduction until 1050)
  if (variety.type == Variety_TypeInvader && world.year < 1050) {
    UInt const individual_number = space.size_x*space.size_y/(world.cell_diameter*world.cell_diameter);
    for (UInt individual_iterator=0; individual_iterator < individual_number; individual_iterator+=1) {
      Float32 x, y;
      unpack_MersenneTwister_Float32( &thread->mersennetwister, &x,
                                      Random_uniform_Float32(thread->mersennetwister) );
      unpack_MersenneTwister_Float32( &thread->mersennetwister, &y,
                                      Random_uniform_Float32(thread->mersennetwister) );
      x *= space.size_x;
      y *= space.size_y;

      Individual const individual = Individual_raw(x, y, 0.02);

      aindividuals = AIndividuals_append(aindividuals, individual, Z_xy(x, y, space.scale));
    }
  }

  // Package and return
  avarieties = AVarieties_append(avarieties, variety);
  asindividuals = ASIndividuals_append(asindividuals, AIndividuals_end(aindividuals));
  return pack_AVarieties_ASIndividuals(avarieties, asindividuals);
}


//---------------------------------------------------------------------------------------------------------------//
//
Individual Individual_raw(Float32 const x, Float32 const y, Float32 const height) {
  Individual individual = { .x = x, .y = y, .height = height };
  return individual;
}


//
FILE* Individual_saveFP(Space const space, World const world, Variety const variety, Individual const individual,
                        char const* const name, FILE* file) {
  if ( fprintf(file, "%g %g %g\n",
               individual.x, individual.y,
               individual.height) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  return file;
}


Individual_FILE_UInt64 Individual_loadFP(Space const space, World const world, Variety const variety,
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

  return pack_Individual_FILE_UInt64(individual, file, line+1);
}


//
AIndividuals Individual_next(AIndividuals aindividuals, Space const space,
                             World const world, Variety const variety, Individual const individual,
                             RWorld const rworld, RVariety const rvariety,
                             RIndividualIn const rindividualin, RIndividualOut const rindividualout,
                             Thread const thread) {
  // Seed
  if (individual.height >= variety.height_mature) {
    UInt32 const seeds = floor( ( variety.fecundity_maximum +
                                  5.0*sin(2.0*M_PI*(world.year + variety.masting_phase)/variety.masting_time) ) *
                                individual.height/variety.height_maximum );

    for (UInt32 seed = 0; seed < seeds; seed += 1) {
      Float32 mean;
      Float32 radius;
      Float32 angle;
      unpack_MersenneTwister_Float32( &thread->mersennetwister, &mean,
                                      Random_uniform_Float32(thread->mersennetwister) );
      unpack_MersenneTwister_Float32( &thread->mersennetwister, &radius,
                                      Random_uniform_Float32(thread->mersennetwister) );
      unpack_MersenneTwister_Float32( &thread->mersennetwister, &angle,
                                      Random_uniform_Float32(thread->mersennetwister) );
      mean = mean < variety.dispersal_method ? variety.dispersal_short : variety.dispersal_long;
      radius = -mean*log1p(-radius);
      angle *= 2.0*M_PI;

      Float const x = remainderf(cosf(angle)*radius, space.size_x) + space.size_x/2.0;
      Float const y = remainderf(sinf(angle)*radius, space.size_y) + space.size_y/2.0;
      Individual const individual = Individual_raw(x, y, 0.02);

      aindividuals = AIndividuals_append(aindividuals, individual, Z_xy(x, y, space.scale));
    }
  }

  // Grow/die
  Float32 const growth = ( variety.growth_rate * individual.height * 
                           logf( variety.height_maximum/individual.height ) *
                           expf( -( rindividualout.number_higher/variety.competition_higher +
                                    rindividualout.number_lower/variety.competition_lower ) ) );
  Float32 const mortality = ( variety.mortality_intrinsic + 
                              ( (1.0-variety.mortality_intrinsic) *
                                variety.mortality_initial*expf(-variety.mortality_decay*growth) ) );

  Float32 survive;
  unpack_MersenneTwister_Float32( &thread->mersennetwister, &survive,
                                  Random_uniform_Float32(thread->mersennetwister) );
  if ( survive >= mortality ) {
    Individual const individual_next = Individual_raw(individual.x, individual.y, individual.height+growth);
    aindividuals = AIndividuals_append(aindividuals, individual_next,
                                       Z_xy(individual_next.x, individual_next.y, space.scale));
  }

  // Package and return
  return aindividuals;
}


//---------------------------------------------------------------------------------------------------------------//
//
RWorld RWorld_raw(UInt64 number_seedling, UInt64 number_adult) {
  RWorld const rworld = { .number_seedling = number_seedling, .number_adult = number_adult };
  return rworld;
}


//
RWorld RWorld_first(Space const space, World const world) {
  return RWorld_raw(0,0);
}
RWorld RWorld_rest(Space const space, World const world, Variety const variety, Individual const individual) {
  return individual.height >= variety.height_mature ? RWorld_raw(0,1) : RWorld_raw(1,0);
}
RWorld RWorld_merge(RWorld const rworld0, RWorld const rworld1) {
  return RWorld_raw(rworld0.number_seedling + rworld1.number_seedling,
                    rworld0.number_adult + rworld1.number_adult);
}


//---------------------------------------------------------------------------------------------------------------//
//
RVariety RVariety_raw(UInt64 number_seedling, UInt64 number_adult) {
  RVariety const rvariety = { .number_seedling = number_seedling, .number_adult = number_adult };
  return rvariety;
}


//
RVariety RVariety_first(Space const space, World const world, Variety const variety) {
  return RVariety_raw(0,0);
}
RVariety RVariety_rest(Space const space, World const world, Variety const variety, Individual const individual) {
  return individual.height >= variety.height_mature ? RVariety_raw(0,1) : RVariety_raw(1,0);
}
RVariety RVariety_merge(RVariety const rvariety0, RVariety const rvariety1) {
  return RVariety_raw(rvariety0.number_seedling + rvariety1.number_seedling,
                      rvariety0.number_adult + rvariety1.number_adult);
}


//---------------------------------------------------------------------------------------------------------------//
//
RIndividualIn RIndividualIn_raw() {
  RIndividualIn const rindividualin = { };
  return rindividualin;
}


//
Float32_Float32_Float32_Float32 RIndividualIn_bound(Space const space, World const world,
                                                    Variety const variety0, Individual const individual0,
                                                    Variety const variety1) {
  return pack_Float32_Float32_Float32_Float32(individual0.x, individual0.y,
                                              individual0.x, individual0.y);
}


Bool RIndividualIn_filter(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                          Variety const variety0, Individual const individual0,
                          Variety const variety1, Individual const individual1) {
  return false;
}


//
RIndividualIn RIndividualIn_first(Space const space, World const world,
                                  Variety const variety, Individual const individual) {
  return RIndividualIn_raw();
}


RIndividualIn RIndividualIn_rest(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                                 Variety const variety0, Individual const individual0,
                                 Variety const variety1, Individual const individual1) {
  return RIndividualIn_raw();
}


RIndividualIn RIndividualIn_merge(RIndividualIn const rindividualin0, RIndividualIn const rindividualin1) {
  return RIndividualIn_raw();
}


//---------------------------------------------------------------------------------------------------------------//
//
RIndividualOut RIndividualOut_raw(UInt64 const number_lower, UInt64 const number_higher) {
  RIndividualOut const rindividualout = { .number_higher = number_higher, .number_lower = number_lower };
  return rindividualout;
}


//
Float32_Float32_Float32_Float32 RIndividualOut_bound(Space const space, World const world,
                                                     Variety const variety0, Individual const individual0,
                                                     Variety const variety1) {
  return pack_Float32_Float32_Float32_Float32
    (individual0.x-world.cell_diameter/2.0, individual0.y-world.cell_diameter/2.0,
     individual0.x+world.cell_diameter/2.0, individual0.y+world.cell_diameter/2.0);
}


Bool RIndividualOut_filter(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                           Variety const variety0, Individual const individual0,
                           Variety const variety1, Individual const individual1) {
  Float32 const delta_x = individual1.x+mirror_x*space.size_x-individual0.x;
  Float32 const delta_y = individual1.y+mirror_y*space.size_y-individual0.y;
  return delta_x*delta_x + delta_y*delta_y <= world.cell_diameter/2.0*world.cell_diameter/2.0;
}


//
RIndividualOut RIndividualOut_first(Space const space, World const world,
                                    Variety const variety, Individual const individual) {
  return RIndividualOut_raw(0, 1);
}


RIndividualOut RIndividualOut_rest(Int const mirror_x, Int const mirror_y, Space const space, World const world,
                                   Variety const variety0, Individual const individual0,
                                   Variety const variety1, Individual const individual1) {
  return individual0.height <= individual1.height ? RIndividualOut_raw(1, 0) : RIndividualOut_raw(0, 1);
}


RIndividualOut RIndividualOut_merge(RIndividualOut const rindividualout0, RIndividualOut const rindividualout1) {
  return RIndividualOut_raw(rindividualout0.number_lower  + rindividualout1.number_lower,
                            rindividualout0.number_higher + rindividualout1.number_higher);
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
