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


#ifndef MODEL_TYPE
#define MODEL_TYPE


//---------------------------------------------------------------------------------------------------------------//
// Single upper case to distinguish types
typedef struct _Thread* Thread;

typedef enum _Variety_Type Variety_Type;

typedef struct _World World;
typedef struct _Variety Variety;
typedef struct _Individual Individual;

typedef struct _RWorld RWorld;
typedef struct _RVariety RVariety;
typedef struct _RIndividualIn RIndividualIn;
typedef struct _RIndividualOut RIndividualOut;

typedef struct _World_FILE_UInt64 World_FILE_UInt64;
typedef struct _Variety_FILE_UInt64 Variety_FILE_UInt64;
typedef struct _Individual_FILE_UInt64 Individual_FILE_UInt64;
typedef struct _Float32_Float32_Float32_Float32 Float32_Float32_Float32_Float32;

//---------------------------------------------------------------------------------------------------------------//

#endif // MODEL_TYPE
