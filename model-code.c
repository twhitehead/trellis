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

#include "model-func.h"
#include "model-data.c"

#include "system-code.c"


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


Box RIndividualIn_bound(World const world, Variety const variety0, Individual const individual0,
			Variety const variety1) {
  Box box = { };
  return box;
}
Box RIndividualOut_bound(World const world, Variety const variety0, Individual const individual1,
			 Variety variety1) {
  Box box = { };
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

int main() {
  Space space;
  World world;
  SVarieties svarieties;
  SSIndividuals ssindividuals;

  State_load_(&space, &world, &svarieties, &ssindividuals, "checkpoint-0.new");
  State_save("checkout-1.txt", space, world, svarieties, ssindividuals);

  SVarieties_end(svarieties);
  SSIndividuals_end(ssindividuals);

  return 0;
}


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_CODE
