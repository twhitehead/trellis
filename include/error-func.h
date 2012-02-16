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


#ifndef ERROR_FUNC
#define ERROR_FUNC

#include <stdarg.h>

#include "error-type.h"


//---------------------------------------------------------------------------------------------------------------//
//
void (__attribute__((format (printf,2,3))) Error_dieErrNo)(int value, char const* format, ...);
void Error_vdieErrNo(int value, char const* format, va_list args);

void (__attribute__((format (printf,3,4))) Error_dieErrNoExplict)(int errno_original, int value, 
                                                                  char const* format, ...);
void Error_vdieErrNoExplict(int errno_original, int value, char const* format, va_list args);

void (__attribute__((format (printf,2,3))) Error_die)(int value, char const* format, ...);
void Error_vdie(int value, char const* format, va_list args);

//---------------------------------------------------------------------------------------------------------------//

#endif // ERROR_FUNC
