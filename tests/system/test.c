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

#include <stdint.h>
#include <inttypes.h>

#include "system.h"
#include "random.h"
                 

//---------------------------------------------------------------------------------------------------------------//
//
void reduce(Space const space, World const world, SVarieties const svarieties, SSIndividuals const ssindividuals,
            RWorld const rworld, SRVarieties const srvarieties,
            SSRIndividualsIn const ssrindividualsin, SSRIndividualsOut const ssrindividualsout) {
  // Outer variety loop
  RWorld rworld_check = RWorld_first(space, world);

  if ( svarieties->number != srvarieties->number ||
       svarieties->number != ssindividuals->number ||
       srvarieties->number != ssrindividualsin->number ||
       srvarieties->number != ssrindividualsout->number )
    Error_die(1, "disagreement over number of variety components");
  UInt const varieties_number = srvarieties->number;

  for (UInt varieties0_index = 0; varieties0_index < varieties_number; varieties0_index += 1) {
    // Break out variety components
    Variety const variety0 = svarieties->variety[varieties0_index];
    RVariety const rvariety = srvarieties->rvariety[varieties0_index];
    SIndividuals const sindividuals0 = ssindividuals->sindividuals[varieties0_index];
    SRIndividualsIn const srindividualsin = ssrindividualsin->srindividualsin[varieties0_index];
    SRIndividualsOut const srindividualsout = ssrindividualsout->srindividualsout[varieties0_index];

    // Outer individual loop
    RVariety rvariety_check = RVariety_first(space, world, variety0);

    if ( sindividuals0.number != srindividualsin->number ||
         sindividuals0.number != srindividualsout->number )
      Error_die(1, "disagreement over number of individual componentso");

    for ( IIndividuals iindividuals0 = IIndividuals_first(sindividuals0);
          iindividuals0.valid;
          iindividuals0 = IIndividuals_next(iindividuals0) ) {
      // Break out individual components
      Individual const individual0 = IIndividuals_individual(iindividuals0);
      RIndividualIn const rindividualin = srindividualsin->rindividualin[iindividuals0.index];
      RIndividualOut const rindividualout = srindividualsout->rindividualout[iindividuals0.index];

      // Perform reductions
      rvariety_check = RVariety_merge(rvariety_check, RVariety_rest(space, world, variety0, individual0));
      rworld_check = RWorld_merge(rworld_check, RWorld_rest(space, world, variety0, individual0));

      // Inner variety loop
      RIndividualIn rindividualin_check = RIndividualIn_first(space, world, variety0, individual0);
      RIndividualOut rindividualout_check = RIndividualOut_first(space, world, variety0, individual0);

      for (UInt varieties1_index = 0; varieties1_index < varieties_number; varieties1_index += 1) {
        // Break out variety components
        Variety const variety1 = svarieties->variety[varieties1_index];
        SIndividuals const sindividuals1 = ssindividuals->sindividuals[varieties1_index];

        // Periodic in loops
        Float32 box_in_ul_x, box_in_ul_y, box_in_lr_x, box_in_lr_y;
        unpack_Float32_Float32_Float32_Float32
          ( &box_in_ul_x, &box_in_ul_y, &box_in_lr_x, &box_in_lr_y,
            RIndividualIn_bound(space, world, variety0, individual0, variety1) );

        Int const mirror_in_ul_x = space.periodic_x ? floorf(box_in_ul_x/space.size_x) : 0;
        Int const mirror_in_ul_y = space.periodic_y ? floorf(box_in_ul_y/space.size_y) : 0;
        Int const mirror_in_lr_x = space.periodic_x ? floorf(box_in_lr_x/space.size_x) : 0;
        Int const mirror_in_lr_y = space.periodic_y ? floorf(box_in_lr_y/space.size_y) : 0;

        for (Int mirror_in_x = mirror_in_ul_x; mirror_in_x <= mirror_in_lr_x; mirror_in_x += 1) {
          for (Int mirror_in_y = mirror_in_ul_y; mirror_in_y <= mirror_in_lr_y; mirror_in_y += 1) {
            // Inner individual in loop
            for ( IIndividuals iindividuals1 = IIndividuals_first(sindividuals1);
                  iindividuals1.valid;
                  iindividuals1 = IIndividuals_next(iindividuals1) ) {
              // Break out individual components
              Individual const individual1 = IIndividuals_individual(iindividuals1);

              // Perform reduction if in in range
              if ( (varieties0_index != varieties1_index || iindividuals0.index != iindividuals1.index ||
                    mirror_in_x != 0 || mirror_in_y != 0) &&
                   RIndividualIn_filter(mirror_in_x, mirror_in_y, space, world, 
                                        variety0, individual0, variety1, individual1) )
                rindividualin_check = RIndividualIn_merge( rindividualin_check,
                                                           RIndividualIn_rest(mirror_in_x, mirror_in_y,
                                                                              space, world,
                                                                              variety0, individual0,
                                                                              variety1, individual1) );
            }
          }
        }

        // Inner individual out loop (order is reversed for reduction components)
        for ( IIndividuals iindividuals1 = IIndividuals_first(sindividuals1);
              iindividuals1.valid;
              iindividuals1 = IIndividuals_next(iindividuals1) ) {
          // Break out individual components
          Individual const individual1 = IIndividuals_individual(iindividuals1);

          // Periodic out loops
          Float32 box_out_ul_x, box_out_ul_y, box_out_lr_x, box_out_lr_y;
          unpack_Float32_Float32_Float32_Float32
            ( &box_out_ul_x, &box_out_ul_y, &box_out_lr_x, &box_out_lr_y,
              RIndividualOut_bound(space, world, variety1, individual1, variety0) );

          Int const mirror_out_ul_x = space.periodic_x ? floorf(box_out_ul_x/space.size_x) : 0;
          Int const mirror_out_ul_y = space.periodic_y ? floorf(box_out_ul_y/space.size_y) : 0;
          Int const mirror_out_lr_x = space.periodic_x ? floorf(box_out_lr_x/space.size_x) : 0;
          Int const mirror_out_lr_y = space.periodic_y ? floorf(box_out_lr_y/space.size_y) : 0;

          for (Int mirror_out_x = mirror_out_ul_x; mirror_out_x <= mirror_out_lr_x; mirror_out_x += 1) {
            for (Int mirror_out_y = mirror_out_ul_y; mirror_out_y <= mirror_out_lr_y; mirror_out_y += 1) {
              // Perform reduction if in out range
              if ( (varieties0_index != varieties1_index || iindividuals0.index != iindividuals1.index ||
                    mirror_out_x != 0 || mirror_out_y != 0) &&
                   RIndividualOut_filter(mirror_out_x, mirror_out_y, space, world,
                                         variety1, individual1, variety0, individual0) )
                rindividualout_check = RIndividualOut_merge( rindividualout_check,
                                                             RIndividualOut_rest(mirror_out_x, mirror_out_y,
                                                                                 space, world,
                                                                                 variety1, individual1,
                                                                                 variety0, individual0) );
            }
          }
        }
      }

      // Verify individual reductions
      if ( rindividualin.ids != rindividualin_check.ids )
        Error_die(1, "individual in reduction didn't verify for %u %"PRIu64,
                  varieties0_index, iindividuals0.index );
      if ( rindividualout.ids != rindividualout_check.ids )
        Error_die(1, "individual out reduction didn't verify for %u %"PRIu64,
                  varieties0_index, iindividuals0.index );
    }

    // Verify variety reduction
    if ( rvariety.ids != rvariety_check.ids )
      Error_die(1, "variety reduction didn't verify for %u (%"PRIu32" != %"PRIu32")",
                varieties0_index, rvariety.ids, rvariety_check.ids);
  }

  // Verify world reduction
  if ( rworld.ids != rworld_check.ids )
    Error_die(1, "world reduction didn't verify");
}


//---------------------------------------------------------------------------------------------------------------//
//
int main() {
  MersenneTwister mersennetwister = MersenneTwister_begin();
  Space const space = Space_raw(true,true, 30.0, 30.0);

  // Generate a new world with a random id
  printf("Generating state...\n");

  UInt32 id;
  unpack_MersenneTwister_UInt32( &mersennetwister, &id,
                                 MersenneTwister_extract_UInt32(mersennetwister) );
  World const world = World_raw(id);
  
  // Generate 10 random varieties of 900 random individuals
  AVarieties avarieties = AVarieties_begin();
  ASIndividuals asindividuals = ASIndividuals_begin();

  for (UInt variety_iterator=0; variety_iterator<10; variety_iterator+=1) {
    // Generate a new variety with a random id
    UInt32 id;
    unpack_MersenneTwister_UInt32( &mersennetwister, &id,
                                   MersenneTwister_extract_UInt32(mersennetwister) );
    Variety const variety = Variety_raw(id);
    avarieties = AVarieties_append(avarieties, variety);

    // Generate 900 random individuals
    AIndividuals aindividuals = AIndividuals_begin();

    for (UInt individual_iterator=0; individual_iterator<900; individual_iterator+=1) {
      // Generate a new individual with a random id, location, and range
      UInt32 id;
      Float32 x, y, in, out;
      unpack_MersenneTwister_UInt32( &mersennetwister, &id,
                                     MersenneTwister_extract_UInt32(mersennetwister) );
      unpack_MersenneTwister_Float32( &mersennetwister, &x,
                                      Random_uniform_Float32(mersennetwister) );
      unpack_MersenneTwister_Float32( &mersennetwister, &y,
                                      Random_uniform_Float32(mersennetwister) );
      unpack_MersenneTwister_Float32_Float32( &mersennetwister, &in, &out,
                                              Random_normal2_Float32(mersennetwister) );
      x *= space.size_x;
      y *= space.size_y;
      in = fabs(in);
      out = fabs(out);

      Individual const individual = Individual_raw(variety_iterator*900+individual_iterator, x, y, in, out);

      aindividuals = AIndividuals_append(aindividuals, individual, Z_xy(x, y, space.scale));
    }

    SIndividuals const sindividuals = AIndividuals_end(aindividuals);

    asindividuals = ASIndividuals_append(asindividuals, sindividuals);
  }

  SVarieties const svarieties = AVarieties_end(avarieties);
  SSIndividuals const ssindividuals = ASIndividuals_end(asindividuals);
  
  // Save state 
  //printf("Saving state...\n");
  //State_save(space, world, svarieties, ssindividuals, "checkpoint.txt");

  // Load state
  //Space space;
  //World world;
  //SVarieties svarieties;
  //SSIndividuals ssindividuals;

  //printf("Loading state...\n");
  //unpack_Space_World_SVarieties_SSIndividuals( &space, &world, &svarieties, &ssindividuals,
  //                                             State_load("checkpoint.txt") );

  // Reduce state
  printf("Reducing state...\n");

  RWorld rworld;
  SRVarieties srvarieties;
  SSRIndividualsIn ssrindividualsin;
  SSRIndividualsOut ssrindividualsout;

  unpack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
    ( &rworld, &srvarieties, &ssrindividualsin, &ssrindividualsout,
      State_reduce(space, world, svarieties, ssindividuals) );

  // Verify reduction
  printf("Verifying reduction...\n");

  reduce(space, world, svarieties, ssindividuals, rworld, srvarieties, ssrindividualsin, ssrindividualsout);

  // Release resources
  printf("Releasing state...\n");

  SRVarieties_end(srvarieties);
  SSRIndividualsIn_end(ssrindividualsin);
  SSRIndividualsOut_end(ssrindividualsout);

  SVarieties_end(svarieties);
  SSIndividuals_end(ssindividuals);

  MersenneTwister_end(mersennetwister);
    
  return 0;
}
