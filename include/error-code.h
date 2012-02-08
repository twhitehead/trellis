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


#ifndef ERROR_CODE
#define ERROR_CODE

#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "error-type.h"
#include "error-func.h"
#include "error-data.h"


//---------------------------------------------------------------------------------------------------------------//
// Error_dieErrNo:  (int, char const*, ...)     -> ()
// Error_vdieErrNo: (int, char const*, va_args) -> ()
//
// Error_dieErrNoExplict:  (int, int, char const*, ...)     -> ()
// Error_vdieErrNoExplict: (int, int, char const*, va_args) -> ()
//
// If given printf style formatted message, print it to stderr followed by a colon, then print a description of
// the current value of errno to stderr followed by a newline, and finally exit with given error value.
//
// The second versions explicitly take an errno value, the first versions use the global errno value.
//
void Error_dieErrNo(int const value, char const* const format, ...) {
  // Setup va_list for args and invoke that version
  va_list args;

  va_start(args, format);
  Error_vdieErrNo(value, format, args);
  va_end(args);
}

void Error_vdieErrNo(int const value, char const* const format, va_list args) {
  // Invoke explicit version with global errno value
  Error_vdieErrNoExplict(errno, value, format, args);
}


void Error_dieErrNoExplict(int const errno_original, int const value, char const* const format, ...) {
  // Setup va_list for args and invoke that version
  va_list args;

  va_start(args, format);
  Error_vdieErrNoExplict(errno_original, value, format, args);
  va_end(args);
}

void Error_vdieErrNoExplict(int const errno_original, int const value, char const* const format, va_list args) {
  // Output any given printf style formatted message to stderr
  if (format) {
    vfprintf(stderr, format, args);
    fprintf(stderr, ": ");
  }

  // Output description of given errno to stderr
  char buffer[ERRNO_BUFFER];

  if (strerror_r(errno_original, buffer, sizeof(buffer)/sizeof(buffer[0])) == 0)
    fprintf(stderr,"%s\n",buffer);
  else
    fprintf(stderr,"error %d (could not be described due to error %d)\n", errno_original, errno);

  // Quit with given value
  exit(value);
}


// Error_die:  (char const*, ...)     -> ()
// Error_vdie: (char const*, va_args) -> ()
//
// If given printf style formatted message, print it to stderr followed by a newline, and then exit with given
// error value.
//
void Error_die(int const value, char const* const format, ...) {
  // Setup va_list for args and invoke that version
  va_list args;

  va_start(args, format);
  Error_vdie(value, format, args);
  va_end(args);
}


void Error_vdie(int const value, char const* const format, va_list args) {
  // Output any given printf style formatted message to stderr
  if (format) {
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
  }

  // Quit with given value
  exit(value);
}


//---------------------------------------------------------------------------------------------------------------//

#endif // ERROR_CODE
