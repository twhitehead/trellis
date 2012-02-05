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

#include "model-type.h"

#include "system-type.h"


//---------------------------------------------------------------------------------------------------------------//
//
RWorld RWorld_first(World world);
RWorld RWorld_rest(World world, Variety variety, Individual individual);
RWorld RWorld_merge(RWorld rworld0, RWorld rworld1);

//
RVariety RVariety_first(World world, Variety variety);
RVariety RVariety_rest(World world, Variety variety, Individual individual);
RVariety RVariety_merge(RVariety rvariety0, RVariety rvariety1);

//
Box RIndividualIn_bound(World world, Variety variety0, Individual individual0, Variety variety1);
Bool RIndividualIn_filter(World world, Variety variety0, Individual individual0,
			  Variety variety1, Individual individual1);

RIndividualIn RIndividualIn_first(World world, Variety variety, Individual individual);
RIndividualIn RIndividualIn_rest(World world, Variety variety0, Individual individual0,
				 Variety variety1, Individual individual1);
RIndividualIn RIndividualIn_merge(RIndividualIn rindividualin0, RIndividualIn rindividualin1);

//
Box RIndividualOut_bound(World world, Variety variety0, Individual individual0, Variety variety1);
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


//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_FUNC
