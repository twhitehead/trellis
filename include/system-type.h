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


#ifndef SYSTEM_TYPE
#define SYSTEM_TYPE

#include <stdbool.h>
#include <stdint.h>


//---------------------------------------------------------------------------------------------------------------//
// Single upper case to distinguish types
typedef bool Bool;

typedef unsigned int UInt;
typedef uint8_t UInt8;
typedef uint16_t UInt16;
typedef uint32_t UInt32;
typedef uint64_t UInt64;

typedef int Int;
typedef int8_t Int8;
typedef int16_t Int16;
typedef int32_t Int32;
typedef int64_t Int64;

typedef float Float;
typedef float Float32;
typedef double Float64;

typedef struct _Space Space;

typedef struct _SVarieties* SVarieties;
typedef struct _SSIndividuals* SSIndividuals;
typedef struct _SIndividuals SIndividuals;

typedef struct _SRVarieties* SRVarieties;
typedef struct _SSRIndividualsIn* SSRIndividualsIn;
typedef struct _SSRIndividualsOut* SSRIndividualsOut;
typedef struct _SRIndividualsIn* SRIndividualsIn;
typedef struct _SRIndividualsOut* SRIndividualsOut;

typedef struct _IZ IZ;
typedef struct _IIndividuals IIndividuals;

typedef struct _SVarieties* AVarieties;
typedef struct _SSIndividuals* ASIndividuals;
typedef struct _SIndividuals AIndividuals;

typedef struct _AVarieties_ASIndividuals AVarieties_ASIndividuals;
typedef struct _World_SVarieties_SSIndividuals World_SVarieties_SSIndividuals;
typedef struct _Space_World_SVarieties_SSIndividuals Space_World_SVarieties_SSIndividuals;
typedef struct _Space_World_SVarieties_SSIndividuals_FILE_UInt64 
Space_World_SVarieties_SSIndividuals_FILE_UInt64;
typedef struct _RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut;


//---------------------------------------------------------------------------------------------------------------//

#endif // SYSTEM_TYPE
