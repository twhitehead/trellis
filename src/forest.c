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

#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>

#include "system.h"
#include "random.h"
                 

//---------------------------------------------------------------------------------------------------------------//
//
int main() {
  // Setup
  Thread const thread = Thread_begin();
  Space space;
  World world;
  SVarieties svarieties;
  SSIndividuals ssindividuals;

  if ( access("checkpoint.txt", R_OK) == 0 )
    unpack_Space_World_SVarieties_SSIndividuals( &space, &world, &svarieties, &ssindividuals,
                                                 State_load("checkpoint.txt") );
  else {
    space = Space_raw(true,true, 20.0*20.0, 20.0*20.0);
    world = World_raw(0, 20.0);
    svarieties = AVarieties_end(AVarieties_begin());
    ssindividuals = ASIndividuals_end(ASIndividuals_begin());
  }


  // Simulate
  while (world.year <= 2000) {
    printf("Simulating year %u...\n", world.year);

    // Reduce
    RWorld rworld;
    SRVarieties srvarieties;
    SSRIndividualsIn ssrindividualsin;
    SSRIndividualsOut ssrindividualsout;
    unpack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
      ( &rworld, &srvarieties, &ssrindividualsin, &ssrindividualsout,
        State_reduce(space, world, svarieties, ssindividuals) );

    // Stats
    UInt const varieties_number = ( svarieties->number <= srvarieties->number ?
                                    svarieties->number : srvarieties->number );
    for (UInt variety_index = 0; variety_index < varieties_number; variety_index+=1) {
      Variety const variety = svarieties->variety[variety_index];
      RVariety const rvariety = srvarieties->rvariety[variety_index];

      if (variety.type == Variety_TypeOriginal)
        printf("  Original number (adult, seedling) = (%"PRIu64",%"PRIu64")\n",
               rvariety.number_adult, rvariety.number_seedling);
      else if (variety.type == Variety_TypeInvader)
        printf("  Invader number (adult, seedling) = (%"PRIu64",%"PRIu64")\n",
               rvariety.number_adult, rvariety.number_seedling);
    }

    // Advance
    World world_next;
    SVarieties svarieties_next;
    SSIndividuals ssindividuals_next;
    unpack_World_SVarieties_SSIndividuals
      ( &world_next, &svarieties_next, &ssindividuals_next,
        State_next(space, world, svarieties, ssindividuals,
                   rworld, srvarieties, ssrindividualsin, ssrindividualsout,
                   thread) );

    // Clean
    SVarieties_end(svarieties);
    SSIndividuals_end(ssindividuals);

    SRVarieties_end(srvarieties);
    SSRIndividualsIn_end(ssrindividualsin);
    SSRIndividualsOut_end(ssrindividualsout);

    world = world_next;
    svarieties = svarieties_next;
    ssindividuals = ssindividuals_next;

    // Checkpoint
    if (world.year == 1000)
      State_save(space, world, svarieties, ssindividuals, "checkpoint.txt");
  }

  // Cleanup
  printf("Saving final world...\n");
  State_save(space, world, svarieties, ssindividuals, "final.txt");

  SSIndividuals_end(ssindividuals);
  SVarieties_end(svarieties);

  Thread_end(thread);
    
  return 0;
}
