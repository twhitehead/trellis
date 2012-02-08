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

#endif // MODEL_CODE
