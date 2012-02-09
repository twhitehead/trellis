//---------------------------------------------------------------------------------------------------------------//
// Copyright (C) 2011-12 Tyson Whitehead
//
// This code is free software; you can redistribute it and/or modify it under the terms of the GNU General
// Public License as published by the Free Software Foundation; either version 2, or (at your option) any later
// version.
//
//---------------------------------------------------------------------------------------------------------------//

#if _POSIX_C_SOURCE < 200112
#error _POSIX_C_SOURCE 200112 or greater required
#endif


#ifndef SYSTEM_FUNC
#define SYSTEM_FUNC

#include <stdio.h>

#include "system-type.h"

#include "model-type.h"


//---------------------------------------------------------------------------------------------------------------//
//
UInt64 Z_xy(Float32 x, Float32 y, Float32 scale);

//
UInt64 Indices_reverse(UInt64 indices);

//
void SVarieties_end(SVarieties svarieites);

//
void SSIndividuals_end(SSIndividuals ssindividuals);

//
SIndividuals SIndividuals_raw(UInt64 number, SIndividuals_ sindividuals_);
SIndividuals SIndividuals_none();
void SIndividuals_end(SIndividuals sindividuals);

//
SIndividuals_ SIndividuals__sindividuals_(SIndividuals_ sindividuals_, UInt64 index, UInt depth);
SIndividuals1 SIndividuals__sindividuals1(SIndividuals_ sindividuals_, UInt64 index);

//
IZ IZ_valid(UInt64 z);
IZ IZ_invalid();
IZ IZ_zSet(IZ iz, UInt64 z);
IZ IZ_nextBox(IZ iz, UInt64 box_ul_z,UInt64 box_lr_z);

//
IIndividuals IIndividuals_valid(UInt64 number, UInt64 index,
                                SIndividuals_ sindividuals_, SIndividuals1 sindividuals1);
IIndividuals IIndividuals_invalid();
IIndividuals IIndividuals_indexSet(IIndividuals iindividuals, UInt64 index, SIndividuals1 sindividuals1);
UInt64 IIndividuals_z(IIndividuals iindividuals);
Individual IIndividuals_individual(IIndividuals iindividuals);
IIndividuals IIndividuals_first(SIndividuals sindividuals);
IIndividuals IIndividuals_next(IIndividuals iindividuals);
IIndividuals IIndividuals_firstZ(SIndividuals sindividuals, IZ iz);
IIndividuals IIndividuals_nextZ(IIndividuals iindividuals, IZ iz);
IIndividuals IIndividuals_firstBox(SIndividuals sindividuals, UInt64 box_ul_z,UInt64 box_lr_z);
IIndividuals IIndividuals_nextBox(IIndividuals iindividuals, UInt64 box_ul_z,UInt64 box_lr_z);

//
AVarieties AVarieties_begin();
AVarieties AVarieties_append(AVarieties avarieties, Variety variety);
SVarieties AVarieties_end(AVarieties avarieties);

ASIndividuals ASIndividuals_begin();
ASIndividuals ASIndividuals_append(ASIndividuals asindividuals, SIndividuals sindividuals);
SSIndividuals ASIndividuals_end(ASIndividuals asindividuals);

AIndividuals AIndividuals_begin();
AIndividuals AIndividuals_append(AIndividuals aindividuals, Individual individual, UInt64 z);
SIndividuals_ AIndividuals_attach(SIndividuals_ sindividuals_, SIndividuals1 sindividuals1, UInt64 index);
SIndividuals AIndividuals_end(AIndividuals aindividuals);
SIndividuals_ AIndividuals_cache(SIndividuals_ sindividuals_, UInt64 number);
SIndividuals_ AIndividuals_sort(SIndividuals_ sindividuals_, UInt64 number, UInt64 pivot_z);
SIndividuals_ AIndividuals_sortBoth(SIndividuals_ sindividuals_,
				    SIndividuals1 left_sindividuals1_start,
				    SIndividuals1 right_sindividuals1_start,
				    UInt64 left_index_start, UInt64 right_index_start,
				    UInt64 pivot_z);

//
SRVarieties SRVarieties_begin(World world, SVarieties svarieties);
void SRVarieties_end(SRVarieties srvarieties);

SRIndividualsIn SRIndividualsIn_begin(World world, Variety variety, SIndividuals sindividuals);
void SRIndividualsIn_end(SRIndividualsIn srindividualsin);

SRIndividualsOut SRIndividualsOut_begin(World world, Variety variety, SIndividuals sindividuals);
void SRIndividualsOut_end(SRIndividualsOut srindividualsout);

SSRIndividualsIn SSRIndividualsIn_begin(World world, SVarieties svarieties, SSIndividuals ssindividuals);
void SSRIndividualsIn_end(SSRIndividualsIn ssrindividualsin);

SSRIndividualsOut SSRIndividualsOut_begin(World world, SVarieties svarieties, SSIndividuals ssindividuals);
void SSRIndividualsOut_end(SSRIndividualsOut ssrindividualsout);

//
void State_save(Space space, World world, SVarieties svarieties, SSIndividuals ssindividuals, char const* name);
FILE* State_saveFP(Space space, World world, SVarieties svarieties, SSIndividuals ssindividuals,
		   char const* name, FILE* file);

Space_World_SVarieties_SSIndividuals State_load(char const* name);
Space_World_SVarieties_SSIndividuals_FILE_UInt64 State_loadFP(char const* name, FILE* file, UInt64 line);

RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut State_reduce
(Space space, World world, SVarieties svarieties, SSIndividuals ssindividuals);
World_SVarieties_SSIndividuals State_next
(Space space, World world, SVarieties svarieties, SSIndividuals ssindividuals,
 RWorld rworld, SRVarieties srvarieties,
 SSRIndividualsIn ssrindividualsin, SSRIndividualsOut ssrindividualsout, Thread thread);

//
World_SVarieties_SSIndividuals pack_World_SVarieties_SSIndividuals
(World first, SVarieties second, SSIndividuals third);
Space_World_SVarieties_SSIndividuals pack_Space_World_SVarieties_SSIndividuals
(Space first, World second, SVarieties third, SSIndividuals fourth);
Space_World_SVarieties_SSIndividuals_FILE_UInt64 pack_Space_World_SVarieties_SSIndividuals_FILE_UInt64
(Space first, World second, SVarieties third, SSIndividuals fourth, FILE* fifth, UInt64 sixth);
RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
pack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
(RWorld first, SRVarieties second, SSRIndividualsIn third, SSRIndividualsOut fourth);

void unpack_World_SVarieties_SSIndividuals(World* first, SVarieties* second, SSIndividuals* third,
                                           World_SVarieties_SSIndividuals tuple);
void unpack_Space_World_SVarieties_SSIndividuals
(Space* first, World* second, SVarieties* third ,SSIndividuals* fourth,
 Space_World_SVarieties_SSIndividuals tuple);
void unpack_Space_World_SVarieties_SSIndividuals_FILE_UInt64
(Space* first, World* second, SVarieties* third, SSIndividuals* fourth, FILE** fifth, UInt64* sixth,
 Space_World_SVarieties_SSIndividuals_FILE_UInt64 tuple);
void unpack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
(RWorld* first, SRVarieties* second, SSRIndividualsIn* third, SSRIndividualsOut* fourth,
 RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut tuple);


//---------------------------------------------------------------------------------------------------------------//

#endif // SYSTEM_FUNC