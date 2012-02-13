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


#ifndef SYSTEM_CODE
#define SYSTEM_CODE

#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#include "system-type.h"
#include "system-func.h"
#include "system-data.h"

#include "error-func.h"
#include "error-code.h"

#include "model-func.h"
#include "model-code.h"


//---------------------------------------------------------------------------------------------------------------//
//
UInt64 Indices_reverse(UInt64 indices);

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

//
SIndividuals_ AIndividuals_attach(SIndividuals_ sindividuals_, SIndividuals1 sindividuals1, UInt64 index);
SIndividuals_ AIndividuals_cache(SIndividuals_ sindividuals_, UInt64 number);
SIndividuals_ AIndividuals_sort(SIndividuals_ sindividuals_, UInt64 number, UInt64 pivot_z);
SIndividuals_ AIndividuals_sortBoth(SIndividuals_ sindividuals_,
				    SIndividuals1 left_sindividuals1_start,
				    SIndividuals1 right_sindividuals1_start,
				    UInt64 left_index_start, UInt64 right_index_start,
				    UInt64 pivot_z);


//---------------------------------------------------------------------------------------------------------------//
// Z_xy: (Float32,Float32, Float32) -> (UInt64)
//
// Put given x and y values in Z format (scale, change to integer, and interleave bits -- y_n,x_n,...,y_0,x_0).
//
// The functional pseudo-code follows.
//
// Z_xy(x,y, scale)
//   let x = round(x*scale)
//       y = round(y*scale)
//   Z_xy_uint(x,y)
//
// Z_xy_uint(x,y)
//   // Non-zero x or y, directly handle bottom x and y bits and obtain rest by recursion
//   x > 0 || y > 0:
//     Z_xy(x/2,y/2)*4 + (y%2)*2 + (x%2)
//   // Zero x and y, done
//   otherwise:
//     0
//
// Using the Z_xy_int code to generate direct conversion tables for 4b values, and then directly converting 32b
// values by treating them as a sequence of eight 4b values gives the following code.
//
UInt64 Z_xy(Float32 const x, Float32 const y, Float32 const scale) {
  Float32 const x_scaled = roundf(x*scale);
  Float32 const y_scaled = roundf(y*scale);
  UInt32 const x_uint = x_scaled < 0 ? 0 : x_scaled > UINT32_MAX ? UINT32_MAX : x_scaled;
  UInt32 const y_uint = y_scaled < 0 ? 0 : y_scaled > UINT32_MAX ? UINT32_MAX : y_scaled;

  uint8_t static const z_table[] = {
    [0x00]=0x00, [0x01]=0x01, [0x02]=0x04, [0x03]=0x05, [0x04]=0x10, [0x05]=0x11, [0x06]=0x14, [0x07]=0x15,
    [0x08]=0x40, [0x09]=0x41, [0x0a]=0x44, [0x0b]=0x45, [0x0c]=0x50, [0x0d]=0x51, [0x0e]=0x54, [0x0f]=0x55,
    [0x10]=0x02, [0x11]=0x03, [0x12]=0x06, [0x13]=0x07, [0x14]=0x12, [0x15]=0x13, [0x16]=0x16, [0x17]=0x17,
    [0x18]=0x42, [0x19]=0x43, [0x1a]=0x46, [0x1b]=0x47, [0x1c]=0x52, [0x1d]=0x53, [0x1e]=0x56, [0x1f]=0x57,
    [0x20]=0x08, [0x21]=0x09, [0x22]=0x0c, [0x23]=0x0d, [0x24]=0x18, [0x25]=0x19, [0x26]=0x1c, [0x27]=0x1d,
    [0x28]=0x48, [0x29]=0x49, [0x2a]=0x4c, [0x2b]=0x4d, [0x2c]=0x58, [0x2d]=0x59, [0x2e]=0x5c, [0x2f]=0x5d,
    [0x30]=0x0a, [0x31]=0x0b, [0x32]=0x0e, [0x33]=0x0f, [0x34]=0x1a, [0x35]=0x1b, [0x36]=0x1e, [0x37]=0x1f,
    [0x38]=0x4a, [0x39]=0x4b, [0x3a]=0x4e, [0x3b]=0x4f, [0x3c]=0x5a, [0x3d]=0x5b, [0x3e]=0x5e, [0x3f]=0x5f,
    [0x40]=0x20, [0x41]=0x21, [0x42]=0x24, [0x43]=0x25, [0x44]=0x30, [0x45]=0x31, [0x46]=0x34, [0x47]=0x35,
    [0x48]=0x60, [0x49]=0x61, [0x4a]=0x64, [0x4b]=0x65, [0x4c]=0x70, [0x4d]=0x71, [0x4e]=0x74, [0x4f]=0x75,
    [0x50]=0x22, [0x51]=0x23, [0x52]=0x26, [0x53]=0x27, [0x54]=0x32, [0x55]=0x33, [0x56]=0x36, [0x57]=0x37,
    [0x58]=0x62, [0x59]=0x63, [0x5a]=0x66, [0x5b]=0x67, [0x5c]=0x72, [0x5d]=0x73, [0x5e]=0x76, [0x5f]=0x77,
    [0x60]=0x28, [0x61]=0x29, [0x62]=0x2c, [0x63]=0x2d, [0x64]=0x38, [0x65]=0x39, [0x66]=0x3c, [0x67]=0x3d,
    [0x68]=0x68, [0x69]=0x69, [0x6a]=0x6c, [0x6b]=0x6d, [0x6c]=0x78, [0x6d]=0x79, [0x6e]=0x7c, [0x6f]=0x7d,
    [0x70]=0x2a, [0x71]=0x2b, [0x72]=0x2e, [0x73]=0x2f, [0x74]=0x3a, [0x75]=0x3b, [0x76]=0x3e, [0x77]=0x3f,
    [0x78]=0x6a, [0x79]=0x6b, [0x7a]=0x6e, [0x7b]=0x6f, [0x7c]=0x7a, [0x7d]=0x7b, [0x7e]=0x7e, [0x7f]=0x7f,
    [0x80]=0x80, [0x81]=0x81, [0x82]=0x84, [0x83]=0x85, [0x84]=0x90, [0x85]=0x91, [0x86]=0x94, [0x87]=0x95,
    [0x88]=0xc0, [0x89]=0xc1, [0x8a]=0xc4, [0x8b]=0xc5, [0x8c]=0xd0, [0x8d]=0xd1, [0x8e]=0xd4, [0x8f]=0xd5,
    [0x90]=0x82, [0x91]=0x83, [0x92]=0x86, [0x93]=0x87, [0x94]=0x92, [0x95]=0x93, [0x96]=0x96, [0x97]=0x97,
    [0x98]=0xc2, [0x99]=0xc3, [0x9a]=0xc6, [0x9b]=0xc7, [0x9c]=0xd2, [0x9d]=0xd3, [0x9e]=0xd6, [0x9f]=0xd7,
    [0xa0]=0x88, [0xa1]=0x89, [0xa2]=0x8c, [0xa3]=0x8d, [0xa4]=0x98, [0xa5]=0x99, [0xa6]=0x9c, [0xa7]=0x9d,
    [0xa8]=0xc8, [0xa9]=0xc9, [0xaa]=0xcc, [0xab]=0xcd, [0xac]=0xd8, [0xad]=0xd9, [0xae]=0xdc, [0xaf]=0xdd,
    [0xb0]=0x8a, [0xb1]=0x8b, [0xb2]=0x8e, [0xb3]=0x8f, [0xb4]=0x9a, [0xb5]=0x9b, [0xb6]=0x9e, [0xb7]=0x9f,
    [0xb8]=0xca, [0xb9]=0xcb, [0xba]=0xce, [0xbb]=0xcf, [0xbc]=0xda, [0xbd]=0xdb, [0xbe]=0xde, [0xbf]=0xdf,
    [0xc0]=0xa0, [0xc1]=0xa1, [0xc2]=0xa4, [0xc3]=0xa5, [0xc4]=0xb0, [0xc5]=0xb1, [0xc6]=0xb4, [0xc7]=0xb5,
    [0xc8]=0xe0, [0xc9]=0xe1, [0xca]=0xe4, [0xcb]=0xe5, [0xcc]=0xf0, [0xcd]=0xf1, [0xce]=0xf4, [0xcf]=0xf5,
    [0xd0]=0xa2, [0xd1]=0xa3, [0xd2]=0xa6, [0xd3]=0xa7, [0xd4]=0xb2, [0xd5]=0xb3, [0xd6]=0xb6, [0xd7]=0xb7,
    [0xd8]=0xe2, [0xd9]=0xe3, [0xda]=0xe6, [0xdb]=0xe7, [0xdc]=0xf2, [0xdd]=0xf3, [0xde]=0xf6, [0xdf]=0xf7,
    [0xe0]=0xa8, [0xe1]=0xa9, [0xe2]=0xac, [0xe3]=0xad, [0xe4]=0xb8, [0xe5]=0xb9, [0xe6]=0xbc, [0xe7]=0xbd,
    [0xe8]=0xe8, [0xe9]=0xe9, [0xea]=0xec, [0xeb]=0xed, [0xec]=0xf8, [0xed]=0xf9, [0xee]=0xfc, [0xef]=0xfd,
    [0xf0]=0xaa, [0xf1]=0xab, [0xf2]=0xae, [0xf3]=0xaf, [0xf4]=0xba, [0xf5]=0xbb, [0xf6]=0xbe, [0xf7]=0xbf,
    [0xf8]=0xea, [0xf9]=0xeb, [0xfa]=0xee, [0xfb]=0xef, [0xfc]=0xfa, [0xfd]=0xfb, [0xfe]=0xfe, [0xff]=0xff
  };

  UInt8 const part[8] = {
    [0] = z_table[(UInt8)(y_uint >>  0 & 0x0f) << 4 | (UInt8)(x_uint >>  0 & 0x0f)],
    [1] = z_table[(UInt8)(y_uint >>  4 & 0x0f) << 4 | (UInt8)(x_uint >>  4 & 0x0f)],
    [2] = z_table[(UInt8)(y_uint >>  8 & 0x0f) << 4 | (UInt8)(x_uint >>  8 & 0x0f)],
    [3] = z_table[(UInt8)(y_uint >> 12 & 0x0f) << 4 | (UInt8)(x_uint >> 12 & 0x0f)],
    [4] = z_table[(UInt8)(y_uint >> 16 & 0x0f) << 4 | (UInt8)(x_uint >> 16 & 0x0f)],
    [5] = z_table[(UInt8)(y_uint >> 20 & 0x0f) << 4 | (UInt8)(x_uint >> 20 & 0x0f)],
    [6] = z_table[(UInt8)(y_uint >> 24 & 0x0f) << 4 | (UInt8)(x_uint >> 24 & 0x0f)],
    [7] = z_table[(UInt8)(y_uint >> 28 & 0x0f) << 4 | (UInt8)(x_uint >> 28 & 0x0f)]
  };

  return ( (UInt64)(part[0]) <<  0 | (UInt64)(part[1]) <<  8 |
           (UInt64)(part[2]) << 16 | (UInt64)(part[3]) << 24 |
           (UInt64)(part[4]) << 32 | (UInt64)(part[5]) << 40 |
           (UInt64)(part[6]) << 48 | (UInt64)(part[7]) << 56 );
}


//---------------------------------------------------------------------------------------------------------------//
// Indices_reverse: (UInt64) -> (UInt64)
//
// For given index, reverse the entries.
//
// This routine works by treating the index as a collection of CLUSTER sized indices.  Shifting the indices off
// of one integer and onto another reverses their order (same as reversing a list).
//
// The functional pseudo-code follows.
//
// Indices_reverse: (UInt64) -> (UInt64)
// Indices_reverse_shift: (UInt64, UInt64, UInt) -> (UInt64)
//
// Indices_reverse(indices)
//   // Start with everything in the forward indices and nothing in the reverse
//   Indices_reverse_shift(indices, 0, DEPTH)
//
// Indices_reverse_shift(forward, reverse, number)
//   // More chunks, shift next one from forward indices to reverse and recurse to handle rest
//   number > 0:
//     Indices_reverse_shift(forward/CLUSTER, reverse*CLUSTER+forward%CLUSTER, number-1)
//   // No more chunks, done
//   otherwise:
//     Indices_reverse
//
// Flattening the recursive calls with loops this becomes the following code.
//
UInt64 Indices_reverse(UInt64 const indices) {
  // For each level, shift the indices off the end of the one integer and onto the other
  UInt64 forward = indices;
  UInt64 reverse = 0;

  for (UInt iterator=0; iterator<DEPTH; iterator+=1) {
    reverse = reverse*CLUSTER + forward%CLUSTER;
    forward /= CLUSTER;
  }

  return reverse;
}


//---------------------------------------------------------------------------------------------------------------//
// Space_raw: (Bool, Bool, Float32, Float32) -> (Space)
//
// Space with given parameters.
//
Space Space_raw(Bool const periodic_x, Bool const periodic_y, Float32 const size_x, Float32 const size_y) {
  Space space = { .periodic_x = periodic_x, .periodic_y = periodic_y,
                  .scale = 0x1p32/fmax(size_x, size_y),
                  .size_x = size_x, .size_y = size_y };
  return space;
}


//---------------------------------------------------------------------------------------------------------------//
// SVarieties_end: (SVarieties*) -> ()
//
// Release resources associated with svarieties.
//
void SVarieties_end(SVarieties const svarieties) {
  free(svarieties);
}


//---------------------------------------------------------------------------------------------------------------//
// SSIndividuals_end: (SSIndividuals*) -> ()
//
// Release resources associated with ssindividuals.
//
void SSIndividuals_end(SSIndividuals const ssindividuals) {
  for (UInt iterator = 0; iterator < ssindividuals->number; iterator+=1)
    SIndividuals_end(ssindividuals->sindividuals[iterator]);

  free(ssindividuals);
}


//---------------------------------------------------------------------------------------------------------------//
// sindividuals: (UInt64, SIndividuals_) -> (SIndividuals)
//
// SIndividuals with for number Individual rooted at sindividuals_.
//
SIndividuals SIndividuals_raw(UInt64 const number, SIndividuals_ const sindividuals_) {
  SIndividuals const sindividuals = { .number = number, .sindividuals_ = sindividuals_ };

  return sindividuals;
}


// SIndividuals_none: () -> (SIndividuals)
//
// SIndividuals with no Individuals.
//
SIndividuals SIndividuals_none() {
  SIndividuals const sindividuals = { .number = 0 };

  return sindividuals;
}


// SIndividuals_end: (*SIndividuals) -> ()
//
// Release resources associated with sindividuals.
//
// The functional pseudo-code follows (just iteration code as pseudo-code system is garbage collected).
//
// SIndividuals_end: (SIndividuals_, UInt64) -> ()
//
// AIndividuals_end_level:  (SIndividuals_, UInt64, UInt64, UInt64) -> ()
// AIndividuals_end_level0: (SIndividuals0, UInt64, UInt64, UInt64) -> ()
// AIndividuals_end_level1: (SIndividuals0, UInt64, UInt64)         -> ()
//
// SIndividuals_end(sindividuals_, number)
//   SIndividuals_end_level(sindividuals_, number, 0, CLUSTER^(DEPTH-1))
//
// SIndividuals_end_level(sindividuals_, number, index, step)
//   // Not at bottom level, handle with level0 routine
//   step > 1:
//     SIndividuals_end_level0(sindividuals_.sindividuals0, number, index, step)
//   // At bottom level, handle with level1 routine
//   otherwise:
//     SIndividuals_end_level1(sindividuals_.sindividuals1, number, index)
//
// SIndividuals_end_level0(sindividuals0, number, index, step):
//   let () =
//         SIndividuals_end_level(sindividuals0.sindividuals_[index/step%CLUSTER], number, index, step/CLUSTER)
//   // Further entires at this level, free them first
//   index+step < number && (index+step)/step%CLUSTER > 0:
//     SIndividuals_end_level0(sindividuals0, number, index+step, step)
//   // No further entries at this level, free it
//     ()
//
// SIndividuals_end_level1(sindividuals1, number, index):
//   // Free it
//   ()
//
// Flattening the recursive calls with loops and a sindividuals_ stack array this becomes the following code.
//
void SIndividuals_end(SIndividuals const sindividuals) {
  SIndividuals_ sindividuals__stack[DEPTH] = { [DEPTH-1] = sindividuals.sindividuals_ };
  UInt64 index = 0;
  UInt bit = 64-64/DEPTH;

  // Iterate through SIndividuals1 levels, removing them and upper SIndividuals0 levels
  while (index < sindividuals.number) {
    // Go down to SIndividuals1 level and free it
    while (bit > 0) {
      sindividuals__stack[bit/(64/DEPTH)-1] =
        sindividuals__stack[bit/(64/DEPTH)].sindividuals0->sindividuals_[(index >> bit) % CLUSTER];
      bit -= 64/DEPTH;
    }

    free(sindividuals__stack[0].sindividuals1);

    // Free all SIndividuals0 levels above that we fall on boundary
    while (bit < 64-64/DEPTH) {
      // Move to level above
      bit += 64/DEPTH;

      // Further entries on this level, free them first to clean up all lower levels
      if (index+CLUSTER < sindividuals.number && (index+CLUSTER >> bit) % CLUSTER > 0)
        break;

      // Free level
      free(sindividuals__stack[bit/(64/DEPTH)].sindividuals0);
    }

    // Advance to next SIndividuals1 level
    index += CLUSTER;
  }
}


//---------------------------------------------------------------------------------------------------------------//
// SIndividuals__sindividuals_: (SIndividuals_, UInt64, UInt) -> (SIndividuals_)
//
// For given index up to given depth, return the SIndividuals_.
//
// The functional pseudo-code follows.
//
// SIndividuals__sindividuals_: (SIndividuals_, UInt64, UInt) -> (SIndividuals_)
// SIndividuals__sindividuals_walk: (SIndividuals_, UInt64, UInt) -> (SIndividuals_)
//
// SIndividuals__sindividuals_(individuals_, index, depth)
//   // Start at top level with reversed index as a list of individual level indices in required order
//   SIndividuals__sindividuals_walk(individuals_, Indices_reverse(index), depth)
//
// SIndividuals__sindividuals_walk(individuals_, indices, depth)
//   // Not at desired depth, use next portion of reversed indices to index next level and recurse
//   depth > 0:
//     SIndividuals__sindividuals_walk(individuals_.individuals0.individuals_[indices%CLUSTER],
//                                     indices/CLUSTER, depth-1)
//   // At desire depth, done
//   otherwise:
//     level
//
// Flattening the recursive calls with loops this becomes the following code.
//
SIndividuals_ SIndividuals__sindividuals_(SIndividuals_ const sindividuals_,
					  UInt64 const index, UInt const depth) {
  // Traverse each level up to the requested depth
  UInt64 indices = Indices_reverse(index);
  SIndividuals_ sindividuals__next = sindividuals_;

  for (UInt iterator=0; iterator<depth; iterator+=1) {
    sindividuals__next = sindividuals__next.sindividuals0->sindividuals_[indices%CLUSTER];
    indices /= CLUSTER;
  }

  return sindividuals__next;
}


// SIndividuals__sindividuals1: (SIndividuals_, UInt64) -> (SIndividuals1)
//
// For given index, return the SIndividuals1.
//
// The functional pseudo-code follows.
//
// SIndividuals__sindividuals1: (SIndividuals_, UInt64) -> (SIndividuals1)
//
// SIndividuals__sindividuals1(sindividuals_, index)
//   SIndividuals__sindividuals_(sindividuals_, index, DEPTH-1).sindividuals1
//
// This becomes the following code.
//
SIndividuals1 SIndividuals__sindividuals1(SIndividuals_ const sindividuals_, UInt64 const index) {
  return SIndividuals__sindividuals_(sindividuals_, index, DEPTH-1).sindividuals1;
}


//---------------------------------------------------------------------------------------------------------------//
// IZ_valid: (UInt64) -> (IZ)
//
// Valid IZ with given values
// 
IZ IZ_valid(UInt64 const z) {
  IZ const iz = { .valid = true, .z = z };

  return iz;
}


// IZ_invalid: (IZ, Bool) -> (IZ)
//
// Invalid IZ.
//
IZ IZ_invalid() {
  IZ const iz = { .valid = false };

  return iz;
}


// IZ_zSet: (IZ, UInt64) -> (IZ)
//
// Duplicate of original IZ with Z set to z.
//
IZ IZ_zSet(IZ const iz, UInt64 const z) {
  // Invalid iz, done
  if (!iz.valid)
    return iz;

  // Make a duplicate with field set
  IZ iz_new = iz;

  iz_new.z = z;
  return iz_new;
}


// IZ_nextBox: (IZ, UInt64,UInt64) -> (IZ)
//
// Smallest Z value greater than given z value inside given box.
//
// This routine works by first checking there are more Z values in the box.  If so, the given z value is a lower
// bound that is increased until the desired Z value (the smallest one greater than the give z value that is
// inside the box) is found.
//
// The given z value is increased by checking the largest 2x2-aligned box immediately following it for overlap
// with given box.  If no overlap is found, the box is skipped by increasing the given z value to the end of the
// box and repeating the procedure.
//
// Once overlap is found, the overlapping box is broken down into 2x2 boxes.  Each 2x2 box without overlap is
// skipped by increasing the given z value to the end of it.  Once overlap is found, the procedure is repeated.
// The Z value is found when the box becomes a singleton.
//
// The functional pseudo-code follows.
//
// IZ_nextBox: (UInt64, UInt64,UInt64) -> (UInt64)
//
// IZ_nextBox(z, box_ul_z,box_lr_z)
//   // Interval [z+1,box_lr_z+1) contains Z value, find it
//   z < box_lr_z:
//     IZ_nextBox_outSize(z, 1, box_ul_z,box_lr_z)
//   // Interval [z+1,box_lr_z+1) does not contain Z value, done
//   otherwise:
//     fail("no more valid z points")
//
// IZ_nextBox_outSize(z, step, box_ul_z,box_lr_z)
//   // On step size boundary, increase step size
//   (z+1) % (step*4) == 0:
//     IZ_nextBox_outSize(z, step*4, box_ul_z,box_lr_z)
//   // Not on step size boundary, check box to next next boundary
//   otherwise:
//     IZ_nextBox_outStep(z, step, box_ul_z,box_lr_z)
//
// IZ_nextBox_outStep(z, step, box_ul_z,box_lr_z)
//   let ul_z = z+1
//   let lr_z = z+step
//   // Next box overlap with target box, locate Z point in it
//   Z_overlap(ul_z,lr_z, box_ul_z,box_lr_z):
//     IZ_nextBox_inSize(z, step, box_ul_z,box_lr_z)
//   // Next box doesn't overlap with target box, skip over it
//   otherwise:
//     IZ_nextBox_outSize(lr_z, step, box_ul_z,box_lr_z)
//
// IZ_nextBox_inSize(z, step, box_ul_z,box_lr_z)
//   // Single point, done
//   step == 1:
//     z+1
//   // Not single point, divide into sub boxes
//   otherwise:
//     IZ_nextBox_inStep(z, step/4, box_ul_z,box_lr_z)
//
// IZ_nextBox_inStep(z, step, box_ul_z,box_lr_z)
//   let ul_z = z+1
//       lr_z = z+step
//   // Next sub box overlaps with target box, locate Z point in it
//   Z_overlap(ul_z,lr_z, box_ul_z,box_lr_z):
//     IZ_nextBox_inSize(z, step, box_ul_z,box_lr_z)
//   // Next sub box doesn't overlap with target box, skip over it
//   otherwise:
//     IZ_nextBox_inStep(lr_z, step, box_ul_z,box_lr_z)
//
// Z_overlap(ul_z0,lr_z0, ul_z1,lr_z1)
//   ( (ul_z0&X_MASK) <= (lr_z1&X_MASK) && (lr_z0&X_MASK) >= (ul_z1&X_MASK) &&
//     (ul_z0&Y_MASK) <= (lr_z1&Y_MASK) && (lr_z0YX_MASK) >= (ul_z1YX_MASK) )
//
// Flattening the recursive calls with loops this becomes the following code.
//
IZ IZ_nextBox(IZ const iz, UInt64 const box_ul_z, UInt64 const box_lr_z) {
  // Interval [z+1,box_lr_z+1) does not contain Z value, done
  if (!iz.valid || iz.z >= box_lr_z)
    return IZ_invalid();

  // Scan forward with progressively larger boxes for overlap with target box (bits below bit are 1)
  UInt64 z = iz.z;
  UInt bit = 0;

  while (1) {
    // While on step size boundary, increase step size
    while ((z+1 & ((UInt64)3<<2*bit)) == 0)
      bit += 1;

    // If previous box overlaps with target box, locate Z point in it
    UInt64 const ul_z = z+1;
    UInt64 const lr_z = z+((UInt64)1<<2*bit);

    if ( (ul_z&X_MASK) <= (box_lr_z&X_MASK) && (lr_z&X_MASK) >= (box_ul_z&X_MASK) &&
         (ul_z&Y_MASK) <= (box_lr_z&Y_MASK) && (lr_z&Y_MASK) >= (box_ul_z&Y_MASK) )
      break;

    // Previous box doesn't overlap with target box, skip over it
    z = lr_z;
  }

  // Scan progressively smaller boxes in overlapping box for first point in target box (bits below bit are 1)
  while (1) {
    // At single point, done
    if (bit == 0)
      break;

    // Decrease step size to divide into sub boxes
    bit -= 1;

    // While previous sub box doesn't overlap with target box, skip over it (could be unrolled to 3 checks)
    while (1) {
      UInt64 const ul_z = z+1;
      UInt64 const lr_z = z+((UInt64)1<<2*bit);

      if ( (ul_z&X_MASK) <= (box_lr_z&X_MASK) && (lr_z&X_MASK) >= (box_ul_z&X_MASK) &&
           (ul_z&Y_MASK) <= (box_lr_z&Y_MASK) && (lr_z&Y_MASK) >= (box_ul_z&Y_MASK) )
        break;

      z = lr_z;
    }
  }

  // Return IZ for located Z
  return IZ_zSet(iz, z+1);
}


//---------------------------------------------------------------------------------------------------------------//
// IIndividuals_valid: (UInt64, UInt64, SIndividuals_, SIndividuals1) -> (IIndividuals)
//
// Valid IIndividuals with given values.
//
IIndividuals IIndividuals_valid(UInt64 const number, UInt64 const index,
                                SIndividuals_ const sindividuals_, SIndividuals1 const sindividuals1) {
  IIndividuals const iindividuals = { .valid = true, .number = number, .index = index,
				      .sindividuals_ = sindividuals_, .sindividuals1 = sindividuals1 };

  return iindividuals;
}


// IIndividuals_invalid: ()) -> (IIndividuals)
//
// Invalid IIndividuals.
//
IIndividuals IIndividuals_invalid() {
  IIndividuals const iindividuals = { .valid = false };

  return iindividuals;
}


// IIndividuals_indexSet: (IIndividuals, UInt64, SIndividuals1) -> (IIndividuals)
//
// Duplicate of original IIndividuals with index and SIndividuals1 set index and sindividuals1.
//
IIndividuals IIndividuals_indexSet(IIndividuals const iindividuals,
				   UInt64 const index, SIndividuals1 const sindividuals1) {
  // Invalid iindividuals, done
  if (!iindividuals.valid)
    return iindividuals;

  // Make a duplicate with field set
  IIndividuals iindividuals_new = iindividuals;

  iindividuals_new.index = index;
  iindividuals_new.sindividuals1 = sindividuals1;
  return iindividuals_new;
}


// IIndividuals_z: (IIndividuals) -> (UInt64)
//
// Z value (of Individual) indexed by iindividuals (must be valid).
//
UInt64 IIndividuals_z(IIndividuals const iindividuals) {
  return iindividuals.sindividuals1->z[iindividuals.index%CLUSTER];
}


// IIndividuals_individual: (IIndividuals) -> (Individual)
//
// Individual indexed by iindividuals (must be valid).
//
Individual IIndividuals_individual(IIndividuals const iindividuals) {
  return iindividuals.sindividuals1->individual[iindividuals.index%CLUSTER];
}


// IIndividuals_first: (SIndividuals) -> (IIndividuals)
//
// First IIndividual for given sindividuals.
//
IIndividuals IIndividuals_first(SIndividuals const sindividuals) {
  // No individuals, done
  if (sindividuals.number == 0)
    return IIndividuals_invalid();

  // Return IIndividuals for first individual
  return IIndividuals_valid(sindividuals.number, 0, sindividuals.sindividuals_,
                            SIndividuals__sindividuals1(sindividuals.sindividuals_,0));

}

// IIndividuals_next: (IIndividuals) -> (IIndividuals)
//
// Next IIndividuals after given iindividuals.
//
IIndividuals IIndividuals_next(IIndividuals const iindividuals) {
  // No more individuals, done
  if (!iindividuals.valid || iindividuals.index+1 == iindividuals.number)
    return IIndividuals_invalid();

  // On different SIndividuals1, advance index and sindividuals1
  if ((iindividuals.index+1)%CLUSTER == 0)
    return IIndividuals_indexSet(iindividuals, iindividuals.index+1,
                                  SIndividuals__sindividuals1(iindividuals.sindividuals_, iindividuals.index+1));
  // On same SIndividuals1, just advance index
  else
    return IIndividuals_indexSet(iindividuals, iindividuals.index+1, iindividuals.sindividuals1);
}


// IIndividuals_firstZ: (SIndividuals, IZ) -> (IIndividuals)
//
// Smallest index with Z value greater or equal to given z value.
//
// This routine works by first checking the interval [0,2^64) contains a valid individual (one with a Z value
// greater or equal to z).  If so, the interval is repeatedly reduced by excluding the right half if the left
// half has a valid individual (as that index is less than all the right half) and the left half otherwise.
//
// The individual is found when the interval becomes a singelton.
//
// The functional pseudo-code follows.  The complexity of SIndividuals indexing is hidden in the [] operator. The
// initial range is fixed at 2^64 to generate fixed length loops (for compile time optimization) and to ensure
// ensure midpoints are 2^N (N=63,...,0) so SIndividuals indexing can progress through the depths in order.
//
// IIndividuals_firstZ: (SIndividuals, UInt64) -> (UInt64)
// IIndividuals_firstZ_size: (SIndividuals, UInt64, UInt64, UInt64) -> (UInt64)
// IIndividuals_firstZ_step: (SIndividuals, UInt64, UInt64, UInt64) -> (UInt64)
//
// IIndividuals_firstZ(sindividuals, z)
//   // In range [0,number), look in [0,number)
//   sindividuals.number > 0 && sindividuals.z[sindividuals.number-1] >= z:
//     IIndividuals_firstZ_step(sindividuals, 0, 2^63, z)
//   // Not in range [0,number), done
//   otherwise:
//     fail("no valid individuals")
//
// IIndividuals_firstZ_size(sindividuals, index, step, z):
//   // At a single point, done
//   step == 1:
//     index
//   // Not at a single point, divide into sub ranges
//   otherwise:
//     IIndividuals_firstZ_step(sindividuals, index, step/2, z)
//
// IIndividuals_firstZ_step(sindividuals, index, step, z)
//   // Not in range [index,index+step), look in [index+step,number)
//   index+step < sindividuals.number && sindividuals.z[index+step-1] < z:
//     IIndividuals_firstZ_size(sindividuals, index+step, step, z)
//   // In range [index,index+step), look in [index,index+step)
//   otherwise:
//     IIndividuals_firstZ_size(sindividuals, index, step, z)
//
// Flattening the recursive calls with loops this becomes the following code.
//
IIndividuals IIndividuals_firstZ(SIndividuals const sindividuals, IZ const iz) {
  // There are no valid individuals with Z >= z, done
  if (sindividuals.number == 0 ||
      sindividuals.sindividuals_.sindividuals0->right_z[sindividuals.number-1 >> 64-64/DEPTH] < iz.z)
    return IIndividuals_invalid();

  // Scan forward in half sized steps 
  SIndividuals_ sindividuals_ = sindividuals.sindividuals_;
  UInt64 const* right_z = sindividuals_.sindividuals0->right_z;
  UInt64 index = 0;
  UInt bit = 63;

  while (1) {
    // Not in range [index,index+step), look in [index+step,number)
    if (((UInt64)1<<bit) < sindividuals.number-index &&
        right_z[(index+((UInt64)1<<bit)-1 >> bit/(64/DEPTH)*(64/DEPTH)) % CLUSTER] < iz.z)
      index += ((UInt64)1<<bit);

    // If at single individual, done
    if (bit == 0)
      break;

    // Decrease step size to divide into sub intervals
    if (bit%(64/DEPTH) == 0) {
      sindividuals_ = sindividuals_.sindividuals0->
        sindividuals_[(index >> bit/(64/DEPTH)*(64/DEPTH)) % CLUSTER];

      if ((bit-1)/(64/DEPTH) > 0)
        right_z = sindividuals_.sindividuals0->right_z;
      else
        right_z = sindividuals_.sindividuals1->z;
    }

    bit -= 1;
  }

  // Return IIndividuals for located individual
  return IIndividuals_valid(sindividuals.number, index, sindividuals.sindividuals_,sindividuals_.sindividuals1);
}


// IIndividuals_nextZ: (IIndividuals, UInt64) -> (IIndividuals)
//
// Smallest index greater than given index with a Z value greater or equal to given z value.
//
// This routine works by first checking the interval [index+1,number) contains a valid individual (one with a Z
// value greater or equal to z).  If so, the given index is a lower bound that is increasing until the desired
// index (the smallest one identifying an individual with a greater or equal Z value to the given one) is found.
//
// The given index is increased by checking the Z range of the largest bit-aligned index interval immediately
// following it for overlap with the given Z value.  If no overlap is found, the index interval is skipped by
// increasing the given index to the start of it and repeating the procedure.
//
// Once overlap is found, the overlapping interval is broken in half into left and right intervals.  If the left
// interval contains a valid individual, the procedure is repeated on it.  Otherwise, it is it repeats on the
// right interval.  The individual is found when the interval becomes a singleton.
//
// The functional pseudo-code follows.  The complexity of SIndividuals indexing is hidden in the [] operator. The
// bit aligning ensures midpoints are 2^N (N=0,...,63) so SIndividuals can progress through depths in order.
//
// IIndividuals_nextZ: (SIndividuals, UInt64, UInt64) -> (UInt64)
//
// IIndividuals_nextZ_outSize: (SIndividuals, UInt64, UInt64, UInt64) -> (UInt64)
// IIndividuals_nextZ_outStep: (SIndividuals, UInt64, UInt64, UInt64) -> (UInt64)
//
// IIndividuals_nextZ_inSize: (SIndividuals, UInt64, UInt64, UInt64) -> (UInt64)
// IIndividuals_nextZ_inStep: (SIndividuals, UInt64, UInt64, UInt64) -> (UInt64)
//
// IIndividuals_nextZ(sindividuals, index, z)
//   // Interval [index+1,number) contains individual, find it
//   sindividuals.number > 0 && sindividuals.z[sindividuals.number-1] >= z:
//     IIndividuals_nextZ_outSize(sindividuals, index, 1, z)
//   // Interval [index+1,number) does not contain individual, done
//   otherwise:
//     fail("no more valid individuals")
//
// IIndividuals_nextZ_outSize(sindividuals, index, step, z)
//   // On step size boundary, increase step size
//   (index+1) % (step*2) == 0:
//     IIndividuals_nextZ_outSize(sindividuals, index, step*2, z)
//   // Not on step size boundary, check interval to next boundary
//   otherwise:
//     IIndividuals_nextZ_outStep(sindividuals, index, step, z)
//
// IIndividuals_nextZ_outStep(sindividuals, index, step, z)
//   // Interval [index+1,index+step+1) does not contain individual, skip over it
//   index+step+1 < sindividuals.number && sindividuals.z[index+step] < z
//     IIndividuals_nextZ_outSize(sindividuals, index+step, step, z)
//   // Interval [index+1,index+step+1) contains individual, locate individual in it
//   otherwise:
//     IIndividuals_nextZ_inSize(sindividuals, index, step, z)
//
// IIndividuals_nextZ_inSize(sindividuals, index, step, z)
//   // Single point, done
//   step == 1:
//     index+1
//   // Not single point, divide into sub intervals
//     IIndividuals_nextZ_inStep(sindividuals, index, step/2, z)
//
// IIndividuals_nextZ_inStep(sindividuals, index, step, z)
//   // Interval [index+1,index+step+1) does not contain individual, skip over it
//   index+step+1 < sindividuals.number && sindividuals.z[index+step] < z
//     IIndividuals_nextZ_inSize(sindividuals, index+step, step, z)
//   // Interval [index+1,index+step+1) contains individual, locate individual in it
//   otherwise:
//     IIndividuals_nextZ_inSize(sindividuals, index, step, z)
//
// Flattening the recursive calls with loops this becomes the following code.
//
IIndividuals IIndividuals_nextZ(IIndividuals const iindividuals, IZ const iz) {
  // There are no valid individuals with Z >= z, done
  if (!iindividuals.valid || !iz.valid ||
      iindividuals.index+1 == iindividuals.number ||
      iindividuals.sindividuals_.sindividuals0->right_z[iindividuals.number-1 >> 64-64/DEPTH] < iz.z)
    return IIndividuals_invalid();

  // Find interval containing suitable individuals (bits below bit are 1)
  SIndividuals_ sindividuals_ = { .sindividuals1 = iindividuals.sindividuals1 };
  UInt64 const* right_z = iindividuals.sindividuals1->z;
  UInt64 index = iindividuals.index;
  UInt bit = 0;

  while (1) {
    // While on step size boundary, increase step size
    while ((index & ((UInt64)1<<bit)) == 1) {
      bit += 1;

      if (bit%(64/DEPTH) == 0) {
        sindividuals_ = SIndividuals__sindividuals_(iindividuals.sindividuals_, index+1, DEPTH-1-bit/(64/DEPTH));
        right_z = sindividuals_.sindividuals0->right_z;
      }
    }

    // Interval [index+1,index+(1<<bit)+1) contains individual, locate individual in it
    if (((UInt64)1<<bit) >= iindividuals.number-index-1 ||
        right_z[(index+((UInt64)1<<bit) >> bit/(64/DEPTH)*(64/DEPTH)) % CLUSTER] >= iz.z)
      break;

    // Interval [index+1,index+(1<<bit)+1) doesn't contain individual, skip over it
    index += ((UInt64)1<<bit);
  }

  // Find first suitable individual in containing interval (bits below bit are 1)
  while(1) {
    // At single individual, done
    if (bit == 0)
      break;

    // Decrease step size to divide into sub intervals
    if (bit%(64/DEPTH) == 0) {
      sindividuals_ = sindividuals_.sindividuals0->
        sindividuals_[(index+1 >> bit/(64/DEPTH)*(64/DEPTH)) % CLUSTER];

      if ((bit-1)/(64/DEPTH) > 0)
        right_z = sindividuals_.sindividuals0->right_z;
      else
        right_z = sindividuals_.sindividuals1->z;
    }

    bit -= 1;

    // If interval [index+1,index+(1<<bit)+1) doesn't contain individual, skip over it
    if (right_z[(index+((UInt64)1<<bit) >> bit/(64/DEPTH)*(64/DEPTH)) % CLUSTER] > iz.z)
      index += ((UInt64)1<<bit);
  }

  // Next index wasn't skipped over
  return IIndividuals_indexSet(iindividuals, index+1, sindividuals_.sindividuals1);
}


// IIndividuals_firstBox: (SIndividuals, UInt64,UInt64) -> (IIndividuals)
//
// Smallest index with Z value inside given box.
//
// This routine works by alternating between finding the next greater Z value that falls in the box and finding
// the next index with a greater or equal Z value.
//
// The first candidate index is one with a Z value greater or equal to the upper left-hand corner.  If the
// candidate index falls in the box it is returned.  If it doesn't, the next greater Z value is found followed by
// the next index value with a greater or equal Z value.  The routine repeats with this as the candidate index.
//
// Consistent scaling and rounding when forming Z values is required to for these functions to correctly iterate
// over all the desired individuals.
//
// The functional pseudo-code follows.  The complexity of SIndividuals indexing is hidden in the [] operator.
// 
// IIndividuals_firstBox(sindividuals, box_ul_z,box_lr_z)
//   // Try one at or after box_ul_z first
//   let index = IIndividuals_firstZ(sindividuals, box_ul_z)
//   IIndividuals_firstBox_possible(sindividuals, index, box_ul_z,box_lr_z)
//
// IIndividuals_firstBox_possible(sindividuals, index, box_ul_z,box_lr_z)
//   // Individual is in box, done 
//   Z_overlap(individuals.z[index],individuals.z[index], box_ul_z,box_lr_z):
//     index
//   // Individual is not in box, try one at or after next Z in box next
//   otherwise:
//     let z_next = IZ_nextBox(sindividuals.z[index], box_ul_z,box_lr_z)
//         index_next = IIndividuals_nextZ(sindividuals, index, z_possible)
//     IIndividuals_firstBox_possible(sindividuals, index_next, box_ul_z,box_lr_z)
//
// Flattening the recursive calls with loops this becomes the following code.
//
IIndividuals IIndividuals_firstBox(SIndividuals const sindividuals,
				   UInt64 const box_ul_z, UInt64 const box_lr_z) {
  // Start with first individual after upper-left hand side of box
  IIndividuals iindividuals_new = IIndividuals_firstZ(sindividuals, IZ_valid(box_ul_z));

  // Iterate if outside of box by alternating between next valid Z followed by next Individual
  while (1) {
    // No more possible individuals, done
    if (!iindividuals_new.valid)
      break;

    // Individual is in box, done
    UInt64 const z = IIndividuals_z(iindividuals_new);
    if ( (z&X_MASK) <= (box_lr_z&X_MASK) && (z&X_MASK) >= (box_ul_z&X_MASK) &&
         (z&Y_MASK) <= (box_lr_z&Y_MASK) && (z&Y_MASK) >= (box_ul_z&Y_MASK) )
      break;

    // Try next individual
    IZ const iz_new = IZ_nextBox(IZ_valid(z), box_ul_z,box_lr_z);
    iindividuals_new = IIndividuals_nextZ(iindividuals_new, iz_new);
  }

  return iindividuals_new;
}


// IIndividuals_nextBox: (IIndividuals, UInt64,UInt64) -> (IIndividuals)
//
// Smallest index greater than given index with a Z value inside given box.
//
// This routine works by alternating between finding the next greater Z value that falls in the box and finding
// the next index with a greater or equal Z value.
//
// The first candidate index is the one one greater than the given one.  If the candidate index falls in the box
// it is returned.  If it doesn't, the next greater Z value is found followed by the next index value with a
// greater or equal Z value.  The routine repeats with this as the candidate index.
//
// Consistent scaling and rounding when forming Z values is required to for these functions to correctly iterate
// over all the desired individuals.
//
// The functional pseudo-code follows.  The complexity of SIndividuals indexing is hidden in the [] operator.
//
// IIndividuals_nextBox(sindividuals, index, box_ul_z,box_lr_z)
//   // Try next individual first first
//   index+1 < sindividuals.number:
//     IIndividuals_nextBox_possible(sindividuals, index+1, box_ul_z,box_lr)
//   otherwise:
//     fail("no more valid individuals")
//
// IIndividuals_nextBox_possible(sindividuals, index, box_ul_z,box_lr_z)
//   // Individual is in box, done 
//   Z_overlap(individuals.z[index],individuals.z[index], box_ul_z,box_lr_z):
//     index
//   // Individual is not in box, try one at or after next Z in box next
//   otherwise:
//     let z_next = IZ_nextBox(sindividuals.z[index], box_ul_z,box_lr_z)
//         index_next = IIndividuals_nextZ(sindividuals, index, z_possible)
//     IIndividuals_nextBox_possible(sindividuals, index_next, box_ul_z,box_lr_z)
//
// Flattening the recursive calls with loops this becomes the following code.
//
IIndividuals IIndividuals_nextBox(IIndividuals const iindividuals,
				  UInt64 const box_ul_z, UInt64 const box_lr_z) {
  // Start with next individual
  IIndividuals iindividuals_new = IIndividuals_next(iindividuals);

  // Iterate if outside of box by alternating between next valid Z followed by next Individual
  while (1) {
    // No more possible individuals, done
    if (!iindividuals_new.valid)
      break;

    // Individual is in box, done
    UInt64 const z = IIndividuals_z(iindividuals_new);
    if ( (z&X_MASK) <= (box_lr_z&X_MASK) && (z&X_MASK) >= (box_ul_z&X_MASK) &&
         (z&Y_MASK) <= (box_lr_z&Y_MASK) && (z&Y_MASK) >= (box_ul_z&Y_MASK) )
      break;

    // Try next individual
    IZ const iz_new = IZ_nextBox(IZ_valid(z), box_ul_z,box_lr_z);
    iindividuals_new = IIndividuals_nextZ(iindividuals_new, iz_new);
  }

  return iindividuals_new;
}


//---------------------------------------------------------------------------------------------------------------//
// AVarieties_begin: () -> (AVarieties)
//
// AVarieties with no Variety.
//
AVarieties AVarieties_begin() {
  AVarieties avarieties;

  // Start with enough space for one
  if ( !(avarieties = malloc(sizeof(struct _SVarieties) + sizeof(Variety))) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for AVarieties", 
		   sizeof(struct _SVarieties) + sizeof(Variety));

  // No entries
  avarieties->number = 0;

  return avarieties;
}


// AVarieties_append: (*AVarieties, Varieties) -> (AVarieties)
//
// Adds svarieties onto end of avarieties.
//
AVarieties AVarieties_append(AVarieties const avarieties, Variety const variety) {
  AVarieties avarieties_new;

  // All space is used up (number == 2^N), expand space by *2
  if ( avarieties->number > 0 && (avarieties->number & (avarieties->number-1)) == 0 ) {
    if ( !(avarieties_new = realloc(avarieties,
                                    sizeof(struct _SVarieties) + sizeof(Variety)*(avarieties->number*2))) )
      Error_dieErrNo(1, "unable to grow Avarieties allocation to %tu bytes",
		     sizeof(struct _SVarieties) + sizeof(Variety)*(avarieties->number*2));
  }
  else
    avarieties_new = avarieties;

  // Add onto end
  avarieties_new->variety[avarieties_new->number] = variety;
  avarieties_new->number += 1;

  return avarieties_new;
}


// AVarieties_end: (*AVarieties) -> (SVarieties)
//
// Convert into SVarieties (assumes AVarieties_ struct is same as SVarieties_).
// 
SVarieties AVarieties_end(AVarieties const avarieties) {
  SVarieties avarieties_new = avarieties;

  // Release extra unused space
  if ( !(avarieties_new = realloc(avarieties,
                                  sizeof(struct _SVarieties) + sizeof(Variety)*avarieties->number)) )
    Error_dieErrNo(1, "unable to shrink AVarieties allocation to %tu bytes",
		   sizeof(struct _SVarieties) + sizeof(Variety)*avarieties->number);

  return avarieties_new;
}


//---------------------------------------------------------------------------------------------------------------//
// ASIndividuals_begin: () -> (ASIndividuals)
//
// AIndividuals with no SIndividuals.
//
ASIndividuals ASIndividuals_begin() {
  ASIndividuals asindividuals;

  // Start with enough space for one
  if ( !(asindividuals = malloc(sizeof(struct _SSIndividuals) + sizeof(SIndividuals))) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for ASIndividuals", 
		   sizeof(struct _SSIndividuals) + sizeof(SIndividuals));

  // No entries
  asindividuals->number = 0;

  return asindividuals;
}


// ASIndividuals_append: (*ASIndividuals, SIndividuals) -> (ASIndividuals)
//
// Adds sindividuals onto end of asindividuals.
//
ASIndividuals ASIndividuals_append(ASIndividuals const asindividuals, SIndividuals const sindividuals) {
  ASIndividuals asindividuals_new;

  // All space is used up (number == 2^N), expand space by *2
  if ( asindividuals->number > 0 && (asindividuals->number & (asindividuals->number-1)) == 0 ) {
    if ( !(asindividuals_new = realloc(asindividuals,
                                       sizeof(struct _SSIndividuals) + 
                                       sizeof(SIndividuals)*(asindividuals->number*2))) )
      Error_dieErrNo(1, "unable to grow ASindividuals allocation to %tu bytes",
		     sizeof(struct _SSIndividuals) + sizeof(SIndividuals)*(asindividuals->number*2));
  }
  else
    asindividuals_new = asindividuals;

  // Add onto end
  asindividuals_new->sindividuals[asindividuals_new->number] = sindividuals;
  asindividuals_new->number += 1;

  return asindividuals_new;
}


// ASIndividuals_end: (*ASIndividuals) -> (SSIndividuals)
//
// Convert into SSIndividuals (assumes ASIndividuals_ struct is same as SSIndividuals_).
// 
SSIndividuals ASIndividuals_end(ASIndividuals const asindividuals) {
  SSIndividuals asindividuals_new;

  // Release extra unused space
  if ( !(asindividuals_new = realloc(asindividuals,
                                     sizeof(struct _SSIndividuals) +
                                     sizeof(SIndividuals)*asindividuals->number)) )
    Error_dieErrNo(1, "unable to shrink ASIndividuals allocation to %tu bytes",
		   sizeof(struct _SSIndividuals) + sizeof(SIndividuals)*asindividuals->number);

  return asindividuals_new;
}


//---------------------------------------------------------------------------------------------------------------//
// AIndividuals_begin: () -> (AIndividuals)
//
// AIndividuals with no Individuals.
//
AIndividuals AIndividuals_begin() {
  AIndividuals const aindividuals = { .number = 0 };

  return aindividuals;
}


// AIndividuals_append: (*AIndividuals, Individual, UInt64) -> (AIndividuals)
//
// Adds individual with Z value z onto end of aindividuals.
//
// The Z value is used for selecting individuals in the IIndividuals_*_box functions.  Consistent scaling and
// rounding when forming Z values is required to for these functions to correctly iterate over all the desired
// individuals.
//
// The functional pseudo-code follows.
//
// AIndividuals_append: (AIndividuals, Individual, UInt64) -> (AIndividuals)
//
// AIndividuals_append(aindividuals, individual, z)
//   // New bottom level required, create one to add individual to
//   aindividuals.number%CLUSTER == 0:
//     // Existing bottom level, attach it and then create new one to add individual to
//     aindividuals.number > 0:
//       aindividuals{ .number = aindividuals.number+1,
//                     .z_min = min(aindividuals.z_min, z), .z_max = max(aindividuals.z_max, z),
//                     .sindividuals_ = AIndividuals_attach(aindividuals.sindividuals_,
//                                                          aindividuals.sindividuals1,
//                                                          aindividuals.number-CLUSTER),
//                     .sindividuals1 = SIndividuals1{ .z[0] = z,
//                                                     .individual[0] = individual } }
//     // No existing bottom level, create new one to add individual to
//     otherwise:
//       aindividuals{ .number = aindividuals.number+1,
//                     .z_min = z, .z_max = z,
//                     .sindividuals1 = SIndividuals1{ .z[0] = z,
//                                                     .individual[0] = individual } }
//   // New bottom level not required, add individual to last one
//   otherwise:
//     aindividuals{ .number = aindividuals.number+1,
//                   .z_min = min(aindividuals.z_min, z), .z_max = max(aindividuals.z_max, z),
//                   .sindividuals1{ .z[number%CLUSTER] = z,
//                                   .individual[number%CLUSTER] = individual } }
//
// This becomes the following code.
//
AIndividuals AIndividuals_append(AIndividuals const aindividuals,
				 Individual const individual, UInt64 const z) {
  AIndividuals aindividuals_new = aindividuals;

  // New bottom level required, create one to add individual to
  if (aindividuals.number%CLUSTER == 0) {
    // Previous bottom level, attach it
    if (aindividuals.number > 0)
      aindividuals_new.sindividuals_ = AIndividuals_attach(aindividuals.sindividuals_, 
                                                           aindividuals.sindividuals1,
                                                           aindividuals.number-CLUSTER);

    // Create new bottom level
    if ( !(aindividuals_new.sindividuals1 = malloc(sizeof(struct _SIndividuals1))) )
      Error_dieErrNo(1, "unable to allocate %tu bytes for new SIndividuals1", sizeof(struct _SIndividuals1));
  }

  // Not first individual, update min and max Z
  if (aindividuals.number > 0) {
    aindividuals_new.z_min = aindividuals_new.z_min <= z ? aindividuals_new.z_min : z;
    aindividuals_new.z_max = aindividuals_new.z_max >= z ? aindividuals_new.z_max : z;
  }
  // First individual, initialize min and max Z
  else {
    aindividuals_new.z_min = z;
    aindividuals_new.z_max = z;
  }

  // Add individual to bottom level
  aindividuals_new.sindividuals1->z[aindividuals.number%CLUSTER] = z;
  aindividuals_new.sindividuals1->individual[aindividuals.number%CLUSTER] = individual;
  aindividuals_new.number += 1;

  return aindividuals_new;
}


// AIndividuals_attach: (*SIndividuals_, SIndividuals1, UInt64) -> (SIndividuals_)
//
// Adds sindividuals1 onto end of sindividuals marked by index (requires index%CLUSTER == 0).
//
// The functional pseudo-code follows.
//
// AIndividuals_attach: (SIndividuals_, SIndividuals1, UInt64) -> (SIndividuals_)
// AIndividuals_attach_old: (SIndividuals_, SIndividuals1, UInt64, UInt) -> (SIndividuals_)
// AIndividuals_attach_new: (SIndividuals1, UInt) -> (SIndividuals_)
//
// AIndividuals_attach(sindividuals_, sindividuals1, index)
//   let indices = Indices_reverse(index)
//   AIndividuals_attach_old(sindividuals_, sindividuals1, indices, DEPTH-1)
//
// AIndividuals_attach_old(sindividuals_, sindividuals1, indices, depth)
//   // Entries below this one, descend into them
//   indices > 0:
//     let sindividuals__next = sindividuals_.sindividuals0.sindividuals_[indices%CLUSTER]
//         sindividuals__next = AIndividuals_attach_old(sindividuals__next, sindividuals1, indices/CLUSTER, depth-1)
//     sindividuals_.sindividuals0.sindividuals_[indices%CLUSTER] = sindividuals__next
//   // No entries below this one, add them
//   otherwise:
//     AIndividuals_attach_new(sindividuals1, depth)
//
// AIndividuals_attach_new(sindividuals1, depth)
//   // SIndividuals0 depth, add SIndividuals0 and continue down
//   depth > 0:
//     let sindividuals_ = AIndividuals_attach_new(sindividuals1, depth-1)
//         sindividuals0 = SIndividuals0.sindividuals_[0] = sindividuals_
//     SIndividuals_.sindividuals0 = sindividuals0
//   // SIndividuals1 depth, add SIndividuals1 and done
//   otherwise:
//     SIndividuals_.sindividuals1 = sindividuals1
//
// Flattening the recursive calls with loops this becomes the following code.
//
SIndividuals_ AIndividuals_attach(SIndividuals_ const sindividuals_,
				  SIndividuals1 const sindividuals1, UInt64 const index) {
  UInt64 indices = Indices_reverse(index);
  SIndividuals_ sindividuals__new = sindividuals_;
  SIndividuals_* sindividuals__link = &sindividuals__new;
  UInt depth = DEPTH-1;

  // While existing intermediate AIndividuals0 entries, descend into them
  while (indices > 0) {
    sindividuals__link = &sindividuals__link->sindividuals0->sindividuals_[indices%CLUSTER];
    indices /= CLUSTER;
    depth -= 1;
  }

  // While needing intermediate AIndividuals0 entries, add them
  while (depth > 0) {
    if ( !(sindividuals__link->sindividuals0 = malloc(sizeof(struct _SIndividuals0))) )
      Error_dieErrNo(1, "unable to allocate %tu bytes for new SIndividuals0", sizeof(struct _SIndividuals0));
    sindividuals__link = &sindividuals__link->sindividuals0->sindividuals_[indices%CLUSTER];
    indices /= CLUSTER;
    depth -= 1;
  }

  // Add sindividuals1
  sindividuals__link->sindividuals1 = sindividuals1;

  return sindividuals__new;
}


// AIndividuals_end: (*AIndividuals) -> (SIndividuals)
//
// Convert into an SIndividuals.
//
// The functional pseudo-code follows.
//
// AIndividuals_end: (AIndividuals) -> (SIndividuals)
//
// AIndividuals_end(aindividuals):
//   // Individuals, attach SIndividuals1, sort them, build the left/right Z cache values, and return SIndividuals
//   aindividuals.number > 0:
//     let sindividuals_ = AIndividuals_attach(aindividuals.sindividuals_, aindividuals.sindividuals1
//                                             (aindividuals.number-1)/CLUSTER*CLUSTER)
//         sindividuals_ = AIndividuals_sort(sindividuals_, aindividuals.number,
//                                           (aindividuals.z_min+aindividuals.z_max)/2)
//         sindividuals_ = AIndividuals_cache(sindividuals_, aindividuals.number)
//     SIndividuals{ .sindividuals_ = sindividuals, .number = aindividuals.number }
//   // No Individuals, return empty SIndividuals
//   otherwise:
//     SIndividuals{ .number = 0 }
//
// This becomes the following code.
//
SIndividuals AIndividuals_end(AIndividuals const aindividuals) {
  // Individuals, attach SIndividuals1, sort them, build the left/right Z cache values, and return SIndividuals
  if (aindividuals.number > 0) {
    SIndividuals_ sindividuals_;

    sindividuals_ = AIndividuals_attach(aindividuals.sindividuals_, aindividuals.sindividuals1,
                                        (aindividuals.number-1)/CLUSTER*CLUSTER);
    sindividuals_ = AIndividuals_sort(sindividuals_, aindividuals.number,
                                      (aindividuals.z_min/2+aindividuals.z_max/2 +
                                       (aindividuals.z_min%2+aindividuals.z_max %2)/2));
    sindividuals_ = AIndividuals_cache(sindividuals_, aindividuals.number);

    return SIndividuals_raw(aindividuals.number, sindividuals_);
  }
  // No Individuals, return empty SIndividuals
  else
    return SIndividuals_none();
}


// AIndividuals_cache: (*SIndividuals_, UInt64) -> (SIndividuals_)
//
// Cache left and right (min and max if sorted) Z values in upper levels.
//
// The functional pseudo-code follows.
//
// AIndividuals_cache: (SIndividuals_, UInt64) -> (SIndividuals_)
//
// AIndividuals_cache_level:  (SIndividuals_, UInt64, UInt64, UInt64) -> (SIndividuals_, UInt64, UInt64)
// AIndividuals_cache_level0: (SIndividuals0, UInt64, UInt64, UInt64) -> (SIndividuals0, UInt64, UInt64)
// AIndividuals_cache_level1: (SIndividuals1, UInt64, UInt64)         -> (SIndividuals1, UInt64, UInt64)
//
// AIndividuals_cache(sindividuals_, number):
//   AIndividuals_cache_level(sindividuals_, number, 0, CLUSTER^(DEPTH-1))
//
// AIndividuals_cache_level(sindividuals_, number, index, step):
//   // Not at bottom level, handle with level0 routine
//   step > 1:
//     let (sindividuals0, left_z, right_z) =
//           AIndividuals_cache_level0(sindividuals_.sindividuals0, number, index, step)
//     (SIndividuals_.sindividuals0 = sindividuals0, left_z, right_z)
//   // At bottom level,  handle with level1 routine
//   otherwise:
//     let (sindividuals1, left_z, right_z) = 
//           AIndividuals_cache_level1(sindividuals_.sindividuals1, number, index)
//     (SIndividuals_.sindividuals1 = sindividuals1, left_z, right_z)
//
// AIndividuals_cache_level0(sindividuals0, number, index, step):
//   let (sindividuals_, left_z, right_z) =
//         AIndividuals_cache_level(sindividuals0.sindividuals_[index/step%CLUSTER], number, index, step/CLUSTER)
//       sindividuals0{ sindividuals_[index/step%CLUSTER] = sindividuals_,
//                      left_z[index/step%CLUSTER] = left_z, right_z[index/step%CLUSTER] = right_z }
//   // Further entries at this level, do them
//   index+step < number && (index+step)/step%CLUSTER > 0:
//     AIndividuals_cache_level0(sindividuals0, number, index+step, step)
//   // No further entries at this level, return left and right values for this level
//   otherwise:
//     (sindividuals0, sindividuals0.left_z[0], sindividuals0.right_z[index/step%CLUSTER])
//
// AIndividuals_cache_level1(sindividuals1, number, index)
//   // Return left and right values for this level
//   (sindividuals1, sindividuals1.z[0], sindividuals1.z[min(index+CLUSTER-1,number-1)%CLUSTER])
//
// Flattening the recursive calls with loops and a sindividuals_ stack array this becomes the following code.
//
SIndividuals_ AIndividuals_cache(SIndividuals_ const sindividuals_, UInt64 const number) {
  SIndividuals_ sindividuals__stack[DEPTH] = { [DEPTH-1] = sindividuals_ };
  UInt64 index = 0;
  UInt bit = 64-64/DEPTH;

  // Iterate through SIndividuals1 levels filling in left and right Z boundaries on upper SIndividuals1 entries
  while (index < number) {
    // Go down to SIndividuals1 level
    while (bit > 0) {
      sindividuals__stack[bit/(64/DEPTH)-1] =
        sindividuals__stack[bit/(64/DEPTH)].sindividuals0->sindividuals_[(index >> bit) % CLUSTER];
      bit -= 64/DEPTH;
    }

    // Left and right most values for SIndividuals1 level are actual Z values
    UInt64 left_z = sindividuals__stack[0].sindividuals1->z[0];
    UInt64 right_z = sindividuals__stack[0].sindividuals1->
      z[(index+CLUSTER-1 < number-1 ? index+CLUSTER-1 : number-1) % CLUSTER];

    // Fill in cached left and right Z values on all SIndividuals0 levels above that we fall on boundary
    while (bit < 64-64/DEPTH) {
      // Move to level above and fill in cached left and right Z values
      bit += 64/DEPTH;
      sindividuals__stack[bit/(64/DEPTH)].sindividuals0->left_z[(index>>bit) % CLUSTER] = left_z;
      sindividuals__stack[bit/(64/DEPTH)].sindividuals0->right_z[(index>>bit) % CLUSTER] = right_z;
      
      // Further entries on this level, process them first to fill in rest of cache
      if (index+CLUSTER < number && (index+CLUSTER >> bit) % CLUSTER > 0)
        break;

      // Load left and right most values for this SIndividauls0 level to update cache above
      left_z = sindividuals__stack[bit/(64/DEPTH)].sindividuals0->left_z[0];
      right_z = sindividuals__stack[bit/(64/DEPTH)].sindividuals0->right_z[(index >> bit) % CLUSTER];
    }

    // Advance to next SIndividuals1 level
    index += CLUSTER;
  }

  return sindividuals__stack[DEPTH-1];
}


// AIndividuals_sort: (*SIndividuals_, UInt64, UInt64) -> (SIndividuals_)
//
// AIndividuals_sortBoth:  (*SIndividuals_, Int64,Int64, Int64) -> (SIndividuals_)
//
// In-place pivot sort.
//
// This routine works by recursively choosing a pivot and breaking elements into a left subset <= pivot and a
// right subset > pivot.  The base case is a subset in which no elements are not equal (maximum and minimum are
// equal).  This implies the subset it fully sorted and will occur for sure once singelton subsets are reached.
// The first pivot is the middle element and later pivots are the average of the minimum and maximum.
//
// The functional pseudo-code follows.  The complexity of SIndividuals indexing is hidden in the [] operator.
//
// AIndividuals_sort: (SIndividuals) -> (SIndividuals)
//
// AIndividuals_sort_both: (SIndividuals, UInt64,UInt64, UInt64) -> (SIndividuals)
// AIndividuals_sort_left:  (SIndividuals, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64, UInt64)
//                            -> (SIndividuals)
// AIndividuals_sort_right: (SIndividuals, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64, UInt64)
//                            -> (SIndividuals)
//
// swap: (SIndividuals, UInt64, UInt64) -> (SIndividuals)
//
// AIndividuals_sort(sindividuals, pivot_z)
//   // If there are individuals, use the centre element as initial pivot
//   sindividuals.number > 0:
//     AIndividuals_sort_both(individuals, 0,sindividuals.number-1, pivot_z)
//   // If there are no individuals, done
//   otherwise:
//     sindividuals
//
// AIndividuals_sort_both(sindividuals, left_start,right_start, pivot_z)
//   // Split into left-hand side <= pivot and right-hand side > pivot
//   let (sindividuals,
//        left_end,  left_min, left_max,
//        right_end, right_min,right_max) = AIndividuals_sort_left(sindividuals, pivot_z,
//                                                                 left_start, left_start,  UINT64_MAX,0,
//                                                                 right_start,right_start, UINT64_MAX,0)
//   // Repeat on each of left-hand and right-hand if they have at least two different elements in them
//       sindividuals = left_min  < left_max:
//           AIndividuals_sort_both(sindividuals, left_start,left_end,    (left_z_min +left_z_max) /2)
//         otherwise:
//           sindividuals
//       sindividuals = right_min < right_max:
//           AIndividuals_sort_both(sindividuals, right_end, right_start, (right_z_min+right_z_max)/2)
//         otherwise:
//           sindividuals
//   // Done
//   sindividuals
//
// AIndividuals_sort_left(sindividuals,
//                        left_start, left_current,  left_z_min, left_z_max,
//                        right_start,right_current, right_z_min,right_z_max,
//                        pivot_z)
//   // Left side hasn't reached right side, process element
//   left_current <= right_current:
//     // Element <= pivot, move it into left-hand side and continue left scan with next element
//     sindividuals.z[left_current] <= pivot_z:
//       AIndividuals_sort_left( sindividuals,
//                               left_start, left_current+1,
//                               min(sindividuals.z[left],left_z_min), max(sindividuals.z[left],left_z_max),
//                               right_start,right_current,   right_z_min,right_z_max,
//                               pivot_z)
//     // Element > pivot, switch to right scan to find element to swap it with
//     sindividuals.z[left_current] >  pivot_z:
//       AIndividuals_sort_right(sindividuals, pivot_z,
//                               left_start, left_current,    left_z_min, left_z_max,
//                               right_start,right_current,   right_z_min,right_z_max,
//                               pivot_z)
//   // Left side reached right side, return sides
//   otherwise:
//     (sindividuals,
//      right_current, left_min, left_max,
//      left_current,  right_min,right_max)
//
// AIndividuals_sort_right(sindividuals, pivot_z,
//            left_start, left_current,  left_z_min, left_z_max,
//            right_start,right_current, right_z_min,right_z_max)
//   // Right side hasn't reached left side, process element
//   left_current <= right_current:
//     // Element > pivot, move it into right-hand side and continue right scan with next element
//     sindividuals[right_current] > pivot_z:
//       AIndividuals_sort_right(sindividuals, pivot_z,
//                               left_start, left_current,    left_z_min,left_z_max,
//                               right_start,right_current-1,
//                               min(sindividuals[right],right_z_min), max(sindividuals[right],right_z_max),
//                               pivot_z)
//     // Element <= pivot, swap it with > element found in left scan and switch back to left scan
//     sindividuals.z[left_current] <= pivot_z:
//       let sindividuals = AIndividuals_swap(sindividuals, left_current, right_current)
//       AIndividuals_sort_left( sindividuals, pivot_z,
//                               left_start, left_current,    left_z_min, left_z_max,
//                               right_start,right_current,   right_z_min,right_z_max,
//                               pivot_z)
//   // Reach side reached left side, return sides
//   otherwise:
//     (sindividuals,
//      right_current, left_min, left_max,
//      left_current,  right_min,right_max)
//
// AIndividuals_swap(sindividuals, index0, index1)
//   let individual0 = sindividuals.individual[index0]
//       individual1 = sindividuals.individual[index1]
//       z0 = sindividuals.z[index0]
//       z1 = sindividuals.z[index1]
//   sindividuals{ .individuals{ [index0] = individuals1, [index1] = individuals0 }
//                 .z{ [index0] = z1, [index1] = z0 } }
//
// Flattening the recursive calls with loops this becomes the following code.
//
SIndividuals_ AIndividuals_sort(SIndividuals_ const sindividuals_, UInt64 const number, UInt64 const pivot_z) {
  // Individuals, pivot sort
  if (number > 0)
    return AIndividuals_sortBoth(sindividuals_,
				 SIndividuals__sindividuals1(sindividuals_, 0),
				 SIndividuals__sindividuals1(sindividuals_, number-1),
				 0, number-1,
				 pivot_z);
  // No individuals, done
  else
    return sindividuals_;
}


SIndividuals_ AIndividuals_sortBoth(SIndividuals_ const sindividuals_,
				    SIndividuals1 const left_sindividuals1_start,
				    SIndividuals1 const right_sindividuals1_start,
				    UInt64 const left_index_start, UInt64 const right_index_start,
				    UInt64 const pivot_z) {
  SIndividuals1 left_sindividuals1_current = left_sindividuals1_start;
  UInt64 left_index_current = left_index_start;
  UInt64 left_z_min = UINT64_MAX;
  UInt64 left_z_max = 0;

  SIndividuals1 right_sindividuals1_current = right_sindividuals1_start;
  UInt64 right_index_current = right_index_start;
  UInt64 right_z_min = UINT64_MAX;
  UInt64 right_z_max = 0;

  // Sweep to centre splitting into left <= pivot and right > pivot by swaping left > pivot and right <= pivot
  while (1) {
    UInt64 left_z;
    UInt64 right_z;

    // Find element from left > pivot
    while (left_index_current <= right_index_current) {
      left_z = left_sindividuals1_current->z[left_index_current%CLUSTER];

      // No element found, move to next and repeat
      if (left_z <= pivot_z) {
        // Update left min and max for new left <= pivot element
        left_z_min = left_z_min <= left_z ? left_z_min : left_z;
        left_z_max = left_z_max >= left_z ? left_z_max : left_z;

        // Advance to next possible left element
        left_index_current += 1;
        if (left_index_current%CLUSTER == 0)
          left_sindividuals1_current = SIndividuals__sindividuals1(sindividuals_, left_index_current);
      }
      // Element found, break
      else
        break;
    }

    // Find element from right <= pivot
    while (left_index_current <= right_index_current) {
      right_z = right_sindividuals1_current->z[right_index_current%CLUSTER];

      // No element found, move to next and repeat
      if (right_z > pivot_z) {
        // Update right min and max for new right > pivot element
        right_z_min = right_z_min <= right_z ? right_z_min : right_z;
        right_z_max = right_z_max >= right_z ? right_z_max : right_z;

        // Advance to next possible right element
        right_index_current -= 1;
        if (right_index_current%CLUSTER == CLUSTER-1)
          right_sindividuals1_current = SIndividuals__sindividuals1(sindividuals_, right_index_current);
      }
      // Element found, break
      else
        break;
    }

    // Not at centre, left > pivot and right <= pivot elements found, swap them
    if (left_index_current <= right_index_current) {
      left_sindividuals1_current->z[left_index_current%CLUSTER] = right_z;
      right_sindividuals1_current->z[right_index_current%CLUSTER] = left_z;

      Individual const left_individual = left_sindividuals1_current->individual[left_index_current%CLUSTER];
      Individual const right_individual = right_sindividuals1_current->individual[right_index_current%CLUSTER];

      left_sindividuals1_current->individual[left_index_current%CLUSTER] = right_individual;
      right_sindividuals1_current->individual[right_index_current%CLUSTER] = left_individual;
    }
    // At centre, seperated into left <= pivot and right > pivot subsets, break
    else
      break;
  }

  // Sort new left <= pivot and right > pivot subsets
  if (left_z_min  < left_z_max)
    AIndividuals_sortBoth(sindividuals_,
			  left_sindividuals1_start,  right_sindividuals1_current,
			  left_index_start,          right_index_current,
			  left_z_min /2+left_z_max /2 + (left_z_min %2+left_z_max %2)/2);
  if (right_z_min < right_z_max)
    AIndividuals_sortBoth(sindividuals_,
			  left_sindividuals1_current,right_sindividuals1_start,
			  left_index_current,        right_index_start,
			  right_z_min/2+right_z_max/2 + (right_z_min%2+right_z_max%2)/2);

  return sindividuals_;
}


//---------------------------------------------------------------------------------------------------------------//
// srvarieties_begin: (World, SVarieties) -> (SRVarieties)
//
// SRVarieties prepared for reduction of svarities (initialized with Reduce_variety_inital).
//
SRVarieties SRVarieties_begin(World const world, SVarieties const svarieties) {
  SRVarieties srvarieties;

  // Create a RVariety for each Variety
  if ( !(srvarieties = malloc(sizeof(struct _SRVarieties) +
                              sizeof(struct _RVariety)*svarieties->number) ) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for SRVarieties",
		   sizeof(struct _SRVarieties)+sizeof(struct _RVariety)*svarieties->number);

  // Initialize each RVariety with RVariety_first
  for (UInt index=0; index<svarieties->number; index+=1) {
    // Break out Variety
    Variety const variety = svarieties->variety[index];

    // Initialize RVariety
    RVariety const rvariety = RVariety_first(world, variety);

    srvarieties->rvariety[index] = rvariety;
  }

  // Number of RVariety
  srvarieties->number = svarieties->number;

  return srvarieties;
}


// SRVarieties_end: (SRVarieties) -> ()
//
// Release resources associated with svarieties.
//
void SRVarieties_end(SRVarieties srvarieties) {
  free(srvarieties);
}


//---------------------------------------------------------------------------------------------------------------//
// SRIndividualsIn_begin: (World, Variety, SIndividuals) -> (SRIndividualsIn)
//
// SRIndividualsIn prepared for reduction of sindividuals (initialized with RIndividualIn_first).
//
SRIndividualsIn SRIndividualsIn_begin(World const world, Variety const variety,
				      SIndividuals const sindividuals) {
  SRIndividualsIn srindividualsin;

  // Create a RIndividualIn for each Individual in sindividuals
  if ( !(srindividualsin = malloc(sizeof(struct _SRIndividualsIn)+
				  sizeof(RIndividualIn)*sindividuals.number) ) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for SRIndividualsIn",
		   sizeof(struct _SRIndividualsIn)+sizeof(RIndividualIn)*sindividuals.number);

  // Initialize each RIndividualIn with RIndividualIn_first
  for ( IIndividuals iindividuals = IIndividuals_first(sindividuals);
        iindividuals.valid;
        iindividuals = IIndividuals_next(iindividuals) ) {
    // Break out Individual
    Individual const individual = IIndividuals_individual(iindividuals);

    // Initialize RIndividualIn
    RIndividualIn const rindividualin = RIndividualIn_first(world, variety, individual);

    srindividualsin->rindividualin[iindividuals.index] = rindividualin;
  }

  // Number of RIndividualsIn
  srindividualsin->number = sindividuals.number;

  return srindividualsin;
}


// SRIndividualsIn_end: (SRIndividualsIn) -> ()
//
// Release resources associated with srindividualsin.
//
void SRIndividualsIn_end(SRIndividualsIn const srindividualsin) {
  free(srindividualsin);
}


// SRIndividualsOut_begin: (World, Variety, SIndividuals) -> (SRIndividualsOut)
//
// SEIndividualsOut prepared for reduction of sindividuals (initialized with RIndividualOut_first).
//
SRIndividualsOut SRIndividualsOut_begin(World const world, Variety const variety,
					SIndividuals const sindividuals) {
  SRIndividualsOut srindividualsout;

  // Create a RIndividualsOut for each Individual
  if ( !(srindividualsout = malloc(sizeof(struct _SRIndividualsOut)+
				   sizeof(RIndividualOut)*sindividuals.number) ) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for SRIndividualsOut",
		   sizeof(struct _SRIndividualsOut)+sizeof(RIndividualOut)*sindividuals.number);

  // Initialize each RIndividualsOut with RIndividualOut_first
  for ( IIndividuals iindividuals = IIndividuals_first(sindividuals);
        iindividuals.valid;
        iindividuals = IIndividuals_next(iindividuals) ) {
    // Break out Individual
    Individual const individual = IIndividuals_individual(iindividuals);

    // Initialize RIndividualOut
    RIndividualOut const rindividualout = RIndividualOut_first(world, variety, individual);

    srindividualsout->rindividualout[iindividuals.index] = rindividualout;
  }

  // Number of RIndividualsOut
  srindividualsout->number = sindividuals.number;

  return srindividualsout;
}


// SRIndividualsOut_end: (SRIndividualsOut) -> ()
//
// Release resources associated with srindividualsout.
//
void SRIndividualsOut_end(SRIndividualsOut const srindividualsout) {
  free(srindividualsout);
}


//---------------------------------------------------------------------------------------------------------------//
// SSRIndividualsIn_begin: (World, SVarieties) -> (SSRIndividualsIn)
//
// SSRIndividualsIn prepared for reduction of ssindividuals (initialized with RIndividualIn_first).
//
SSRIndividualsIn SSRIndividualsIn_begin(World const world, SVarieties const svarieties,
					SSIndividuals const ssindividuals) {
  UInt number = svarieties->number <= ssindividuals->number ? svarieties->number : ssindividuals->number;
  SSRIndividualsIn ssrindividualsin;

  // Create a SRIndivdiualsIn for each (Variety,SIndividuals) pair
  if ( !(ssrindividualsin = malloc(sizeof(struct _SSRIndividualsIn) +
				   sizeof(SRIndividualsIn)*number) ) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for SSRIndividualsIn",
		   sizeof(struct _SSRIndividualsIn)+sizeof(SRIndividualsIn)*number);

  // Initialize each SRIndividualsIn with SRIndividualsIn_begin (uses RIndividualIn_first)
  for (UInt index=0; index<number; index+=1) {
    // Break out (Variety,SIndividuals) pair
    Variety const variety = svarieties->variety[index];
    SIndividuals const sindividuals = ssindividuals->sindividuals[index];

    // Initialize SRIndividualsIn
    SRIndividualsIn const srindividualsin = SRIndividualsIn_begin(world, variety, sindividuals);

    ssrindividualsin->srindividualsin[index] = srindividualsin;
  }

  // Number of SRIndividualsIn
  ssrindividualsin->number = number;

  return ssrindividualsin;
}


// SSRIndividualsIn_end: (SSRIndividualsIn) -> ()
//
// Release resources associated with ssrindividualsin.
//
void SSRIndividualsIn_end(SSRIndividualsIn const ssrindividualsin) {
  // Release SRIndividualsIn resources
  for (UInt index=0; index<ssrindividualsin->number; index+=1) {
    // Break in SRIndividualsIn
    SRIndividualsIn const srindividualsin = ssrindividualsin->srindividualsin[index];

    // Release SRIndividualsIn resources
    SRIndividualsIn_end(srindividualsin);
  }

  // Release SSRIndividualsIn
  free(ssrindividualsin);
}


// SSRIndividualsOut_begin: (World, SVarieties) -> (SSRIndividualsOut)
//
// SSRIndividualsOut prepared for reduction of ssindividuals (initialized with RIndividualOut_first).
//
SSRIndividualsOut SSRIndividualsOut_begin(World const world, SVarieties const svarieties,
					  SSIndividuals const ssindividuals) {
  UInt const number = svarieties->number <= ssindividuals->number ? svarieties->number : ssindividuals->number;
  SSRIndividualsOut ssrindividualsout;

  // Create a SRIndivdiualsOut for each (Variety,SIndividuals) pair
  if ( !(ssrindividualsout = malloc(sizeof(struct _SSRIndividualsOut) +
				    sizeof(SRIndividualsOut)*number) ) )
    Error_dieErrNo(1, "unable to allocate %tu bytes for SSRIndividualsOut",
		   sizeof(struct _SSRIndividualsOut)+sizeof(SRIndividualsOut)*number);

  // Initialize each SRIndividualsOut with SRIndividualsOut_begin (uses RIndividualOut_first)
  for (UInt index=0; index<number; index+=1) {
    // Break out (Variety,SIndividuals) pair
    Variety const variety = svarieties->variety[index];
    SIndividuals const sindividuals = ssindividuals->sindividuals[index];

    // Initialize SRIndividualsOut
    SRIndividualsOut const srindividualsout = SRIndividualsOut_begin(world, variety, sindividuals);

    ssrindividualsout->srindividualsout[index] = srindividualsout;
  }

  // Number of SRIndividualsOut
  ssrindividualsout->number = number;

  return ssrindividualsout;
}


// SSRIndividualsOut_end: (SSRIndividualsOut) -> ()
//
// Release resources associated with ssrindividualsout.
//
void SSRIndividualsOut_end(SSRIndividualsOut const ssrindividualsout) {
  // Release SRIndividualsOut resources
  for (UInt index=0; index<ssrindividualsout->number; index+=1) {
    // Break in SRIndividualsOut
    SRIndividualsOut const srindividualsout = ssrindividualsout->srindividualsout[index];

    // Release SRIndividualsOut resources
    SRIndividualsOut_end(srindividualsout);
  }

  // Release SSRIndividualsOut
  free(ssrindividualsout);
}


//---------------------------------------------------------------------------------------------------------------//
// Tuples
World_SVarieties_SSIndividuals pack_World_SVarieties_SSIndividuals
(World const first, SVarieties const second, SSIndividuals const third) {
  World_SVarieties_SSIndividuals const value = { .first = first, .second = second, .third = third };
  return value;
}

Space_World_SVarieties_SSIndividuals pack_Space_World_SVarieties_SSIndividuals
(Space const first, World const second, SVarieties const third, SSIndividuals const fourth) {
  Space_World_SVarieties_SSIndividuals const value =
    { .first = first, .second = second, .third = third, .fourth = fourth };
  return value;
};

Space_World_SVarieties_SSIndividuals_FILE_UInt64 pack_Space_World_SVarieties_SSIndividuals_FILE_UInt64
(Space const first, World const second, SVarieties const third, SSIndividuals const fourth,
 FILE* const fifth, UInt64 const sixth) {
  Space_World_SVarieties_SSIndividuals_FILE_UInt64 const value =
    { .first = first, .second = second, .third = third, .fourth = fourth, .fifth = fifth, .sixth = sixth };
  return value;
};

RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
pack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
(RWorld const first, SRVarieties const second, SSRIndividualsIn const third, SSRIndividualsOut const fourth) {
  RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut const value =
    { .first = first, .second = second, .third = third, .fourth = fourth };
  return value;
}


//
void unpack_World_SVarieties_SSIndividuals
(World* const first, SVarieties* const second, SSIndividuals* const third,
 World_SVarieties_SSIndividuals const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
}

void unpack_Space_World_SVarieties_SSIndividuals
(Space* const first, World* const second, SVarieties* const third, SSIndividuals* const fourth,
 Space_World_SVarieties_SSIndividuals const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
  *fourth = tuple.fourth;
};

void unpack_Space_World_SVarieties_SSIndividuals_FILE_UInt64
(Space* const first, World* const second, SVarieties* const third, SSIndividuals* const fourth,
 FILE** const fifth, UInt64* const sixth,
 Space_World_SVarieties_SSIndividuals_FILE_UInt64 const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
  *fourth = tuple.fourth;
  *fifth = tuple.fifth;
  *sixth = tuple.sixth;
};

void unpack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut
(RWorld* const first, SRVarieties* const second, SSRIndividualsIn* const third, SSRIndividualsOut* const fourth,
 RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut const tuple) {
  *first = tuple.first;
  *second = tuple.second;
  *third = tuple.third;
  *fourth = tuple.fourth;
}


//---------------------------------------------------------------------------------------------------------------//
// State_save: (char*, Space, World, SVarieties, SSIndividuals) -> ()
// State_saveFP: (FILE*, char*, Space, World, SVarieties, SSIndividuals) -> ()
//
// Create/truncate the given file name and serialize the given space to it.
//
void State_save(Space const space,
		World const world, SVarieties const svarieties, SSIndividuals const ssindividuals,
		char const* const name) {
  // Wrap file handle based routine with opening and close details
  FILE* file;

  if ( !(file = fopen(name,"w")) )
    Error_dieErrNo(1, "unable to create/truncate \"%s\" for writing", name);
  file = State_saveFP(space, world, svarieties, ssindividuals, name, file);
  if (fclose(file))
    Error_dieErrNo(1, "unable to close \"%s\"", name);
}

FILE* State_saveFP(Space const space,
		   World const world, SVarieties const svarieties, SSIndividuals const ssindividuals,
		   char const* const name, FILE* file) {
  // Save space
  if ( fprintf(file, "SPACE %u %u %g %g\n",
	       space.periodic_x, space.periodic_y,
	       space.size_x, space.size_y) < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  // Save world
  if ( fprintf(file, "WORLD ") < 0 )
    Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

  file = World_saveFP(space, world, name, file);

  // Save varieties
  UInt const varieties_number = ( svarieties->number <= ssindividuals->number ?
				  svarieties->number : ssindividuals->number );

  for (UInt varieties_index = 0; varieties_index<varieties_number; varieties_index+=1) {
    // Break out variety components
    Variety const variety = svarieties->variety[varieties_index];
    SIndividuals const sindividuals = ssindividuals->sindividuals[varieties_index];

    // Save variety
    if ( fprintf(file, "VARIETY ") < 0 )
      Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

    file = Variety_saveFP(space, world, variety, name, file);

    // Save individuals
    for ( IIndividuals iindividuals = IIndividuals_first(sindividuals);
          iindividuals.valid;
          iindividuals = IIndividuals_next(iindividuals) ) {
      // Break out individual components
      Individual const individual = IIndividuals_individual(iindividuals);

      // Save individual
      if ( fprintf(file, "INDIVIDUAL ")  < 0 )
        Error_dieErrNo(1, "an error occured while writing to \"%s\"", name);

      file = Individual_saveFP(space, world, variety, individual, name, file);
    }
  }

  return file;
}


// State_load: (char*) -> (Space, World, SVarieties, SSIndividuals)
// State_loadFP: (FILE*, char*) -> (Space, World, SVarieties, SSIndividuals)
//
// Read the given file and restore a serialized forest from it.
//
Space_World_SVarieties_SSIndividuals State_load(char const* const name) {
  FILE* file;
  UInt64 line;
  Space space;
  World world;
  SVarieties svarieties;
  SSIndividuals ssindividuals;

  if ( !(file = fopen(name,"r")) )
    Error_dieErrNo(1, "unable to open \"%s\" for reading", name);
  line = 1;

  unpack_Space_World_SVarieties_SSIndividuals_FILE_UInt64
    ( &space, &world, &svarieties, &ssindividuals, &file, &line,
      State_loadFP(name, file, line) );

  if (fclose(file))
    Error_dieErrNo(1, "unable to close \"%s\"", name);

  return pack_Space_World_SVarieties_SSIndividuals(space, world, svarieties, ssindividuals);
}


Space_World_SVarieties_SSIndividuals_FILE_UInt64
State_loadFP(char const* const name, FILE* file, UInt64 line) {
  char type[11];
  int records;

  // Load type (should be SPACE)
  records = fscanf(file, "%10s", type);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 1 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting \"SPACE\"", name, line);
  
  // Load space data
  Space space;
  UInt periodic_x, periodic_y;

  if ( strcmp(type, "SPACE") != 0 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting \"SPACE\"", name, line);

  records = fscanf(file, "%u %u %g %g\n",
                   &periodic_x, &periodic_y,
                   &space.size_x, &space.size_y);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records < 4 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
	      "PERIODIC_X PERIODIC_Y "
	      "SIZE_X SIZE_Y", name, line);
  line += 1;

  space.scale = 0x1p32/fmax(space.size_x,space.size_y);
  space.periodic_x = periodic_x;
  space.periodic_y = periodic_y;

  // Load next type (should be WORLD)
  records = fscanf(file, "%10s", type);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records == EOF || records < 1 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting \"WORLD\"", name, line);
  
  // Load world data
  World world;

  if ( strcmp(type, "WORLD") != 0 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting \"WORLD\"", name, line);

  unpack_World_FILE_UInt64( &world, &file, &line,
                            World_loadFP(space, name, file, line) );

  // Load type (should be VARIETY or nothing)
  records = fscanf(file, "%10s", type);
  if ( ferror(file) )
    Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
  if ( records != EOF && records < 1 )
    Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting \"VARIETY\" or nothing", name, line);

  // Load varieties
  AVarieties avarieties = AVarieties_begin();
  ASIndividuals asindividuals = ASIndividuals_begin();

  while (records != EOF) {
    // Load variety
    Variety variety;

    if ( strcmp(type, "VARIETY") != 0 )
      Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting \"VARIETY\" or nothing", name, line);

    unpack_Variety_FILE_UInt64( &variety, &file, &line,
                                Variety_loadFP(space, world, name, file, line) );

    // Add variety to varieties
    avarieties = AVarieties_append(avarieties, variety);

    // Load type (should be VARIETY, INDIVIDUAL, or nothing)
    records = fscanf(file, "%10s", type);
    if ( ferror(file) )
      Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
    if ( records != EOF && records < 1 )
      Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
		"\"VARIETY\", \"INDIVIDUAL\", or nothing", name, line);

    // Load variety individuals
    AIndividuals aindividuals = AIndividuals_begin();

    while (records != EOF) {
      // Load individual
      Individual individual;

      if ( strcmp(type, "INDIVIDUAL") != 0 )
	break;
    
      unpack_Individual_FILE_UInt64( &individual, &file, &line,
                                     Individual_loadFP(space, world, variety, name, file, line) );

      // Add individual to individuals
      aindividuals = AIndividuals_append(aindividuals, individual, Z_xy(individual.x, individual.y, space.scale));

      // Load type (should be VARIETY, INDIVIDUAL, or nothing)
      records = fscanf(file, "%10s", type);
      if ( ferror(file) )
	Error_dieErrNo(1, "an error occured while reading from \"%s\"", name);
      if ( records != EOF && records < 1 )
	Error_die(1, "problem parsing \"%s\":%"PRIu64": expecting "
		  "\"VARIETY\", \"INDIVIDUAL\", or nothing", name, line);
    }
    
    // Add individuals to vairety individuals
    SIndividuals sindividuals = AIndividuals_end(aindividuals);

    asindividuals = ASIndividuals_append(asindividuals, sindividuals);
  }

  // Finished additions
  SVarieties const svarieties = AVarieties_end(avarieties);
  SSIndividuals const ssindividuals = ASIndividuals_end(asindividuals);

  return
    pack_Space_World_SVarieties_SSIndividuals_FILE_UInt64(space, world, svarieties, ssindividuals, file, line);
}


// State_reduce: (Space, World, SVarieties, SSIndividuals)
//                 -> (RWorld, SRVarieties, SSRIndividualsIn, SSRIndividualsOut)
//
// Compute the reduction spaces
//
//   RWorld    - fold (variety,individual) by ()
//   SRVariety - fold (variety,individual) by (variety)
//   SSRIndividualsIn  - fold (variety0,individual0,variety1,individual1) by (variety0,individual0) for each
//                       (variety1,individual1) not equal to and in in range of (variety0,individual0)
//   SSRIndividualsOut - fold (variety0,individual0,variety1,individual1) by (variety1,individual1)
//                       (variety1,individual1) not equal to and in out range of (variety0,individual0)
//
// The initial value of the reduction spaces are given by the R*_first functions mapped across the tuples
// indicates above.  These are then folded with the values given by R*_next functions mapped across the the fold
// tuples above via the R*_merge functions.
//
// The individiuals in the in and out ranges are determined in two stages for efficiency.  First the
// RIndividual*_bound functions give a bound which is used as a rough filter.  This set is then refined via the
// RIndividual*_filter functions.
//
// Note that rounding down the Z components of the upper-left corner of the bound and up the lower-right corner
// of the bound ensures that it captures all individuals in the bound whose Z components will also be rounded.
//
// The functional pseudo-code follows.
//
// State_reduce: (World, SVarieties, SSIndividuals) -> (RWorld, SRVarieties, SSRIndividualsIn, SSRIndividualsOut)
//
// State_reduce(world, svarieties, ssindividuals)
//   let // Initial World reduction
//       rworld = RWorld_first(world)
//
//       // Initial Individual (value individual) reduction code
//       State_reduce_individual((variety), individual)
//         // RIndividualIn and RIndividualOut initialization
//         (RIndividualIn_first(world, variety, individual),
//          RIndividualOut_first(world, variety, individual))
//
//       // Initial Variety (value variety) reduction code
//       State_reduce_variety((), (variety, sindividuals))
//         let // RVariety initialization
//             rvariety = RVariety_first(world, variety)
//             // Initial Individual (value individual) reduction map
//             (srindividualsin, srindividualsout) = map(State_reduce_individual, (variety), sindividuals)
//         (rvariety, srindividualsin, srindividualsout)
//
//       // Initial Variety (value variety) reduction map
//       (srvarieties, ssrindividualsin, ssrindividualsout) =
//         map(State_reduce_variety, (), zip(svarieties, ssindividuals))
//
//
//       // Inner Individual (index i1) in reduction code
//       State_reduce_individual1In((ssrindividualsin, v0, i0, v1), i1)
//         // First individual is not second individual and second is in in range, reduce second to first
//         (v0 /= v1 || i0 /= i1) && RIndividualIn_filter(world,
//                                                        svarieties[v0], ssindividuals[v0][i0],
//                                                        svarieties[v0], ssindividuals[v1][i1]):
//           // Inner level SSRIndividualsIn reduction
//           (ssrindividualsin[v0][i0] =
//              RIndividualIn_merge(ssrindividualsin[v0][i0],
//                                  RIndividualIn_next(world,
//                                                     svarieties[v0], ssindividuals[v0][i0],
//                                                     svarieties[v0], ssindividuals[v1][i1])),
//            v0, i0, v1)
//         // First individual is second or second is not in in range, don't reduce second to first
//         otherwise:
//           (ssrindividualin, v0, i0, v1)
//
//       // Inner Individual (index i1) out reduction code
//       State_reduce_individual1Out((ssrindividualsout, v0, i0, v1), i1)
//         // First individual in not second individual and second is in out range, reduce first to second
//         (v0 /= v1 || i0 /= i1) && RIndividualOut_filter(world,
//                                                         svarieties[v0], ssindividuals[v0][i0],
//                                                         svarieties[v0], ssindividuals[v1][i1]):
//           // Inner level SSRIndividualsOut reduction
//           (ssrindividualsout[v1][i1] =
//              RIndividualOut_merge(ssrindividualsout[v1][i1],
//                                   RIndividualOut_next(world,
//                                                       svarieties[v0], ssindividuals[v0][i0],
//                                                       svarieties[v1], ssindividuals[v1][i1])),
//            v0, i0, v1)
//         // First individual is second or second is not in our range, don't reduce first to second
//         otherwise:
//           ssrindividualsout
//
//       // Inner Variety (index v1) reduction code
//       State_reduce_variety1((ssrindividualsin, ssrindividualsout, v0, i0), v1)
//         let // Inner Individual (index i1) reduction fold
//             ssrindividualsin =
//               fold(State_reduce_individual1In,
//                    (ssrindividualsin, v0, i0, v1),
//                    Indicies_filterBox(ssindividuals[v1],
//                                       RIndividualIn_bound(world,
//                                                           svarieties[v0], ssindividuals[v0][i0],
//                                                           svarieties[v1])))
//             // Inner Individual (index i1) reduction fold
//             ssrindividualsout =
//               fold(State_reduce_individual1Out,
//                    (ssrindividualsout, v0, i0, v1),
//                    Indicies_filterBox(ssindividuals[v1],
//                                       RIndividualOut_bound(world,
//                                                            svarieties[v0], ssindividuals[v0][i0],
//                                                            svarieties[v1])))
//         (ssrindividualsin, ssrindividualsout, v0, i0)
//
//       // Outer Individual (index i0) reduction code
//       State_reduce_individual0((rworld, srvarieties, ssrindividualsin, ssrindividualsout, v0), i0)
//         let // Outer level RWorld and SRVariety reductions
//             rworld = RWorld_merge(rworld,
//                                   RWorld_rest(world, svarieties[v0], ssindividuals[v0][i0]))
//             srvarieties[v0] =
//               RVariety_merge(srvarieties[v0],
//                              RVariety_rest(world, svarieties[v0], ssindividuals[v0][i0]))
//             // Outer Individual (index i0) reduction fold
//             (ssrindividualsin, ssrindividualsout) =
//               fold(State_reduce_individual0,
//                    (rworld, srvarieties, ssrindividualsin, ssrindividualsout, v0, i0),
//                    Indicies_all(ssindividuals[v0]))
//         (rworld, srvarieties, ssrindividualsin, ssrindividualsout, v0)
//
//       // Outer Variety (index v0) reduction code
//       State_reduce_variety0((rworld, srvarieties, ssrindividualsin, ssrindividualsout), v0)
//         let // Outer Individual (index i0) reduction fold
//             (rworld, srvarieties, ssrindividualsin, ssrindividualsout) =
//               fold(State_reduce_individual0,
//                    (rworld, srvarieties, ssrindividualsin, ssrindividualsout, v0),
//                    Indicies_all(ssindividuals[v0]))
//         (rworld, srvarieties, ssrindividualsin, ssrindividualsout)
//
//       // Outer Variety (index v0) reduction fold
//       (rworld, srvarieties, ssrindividualsin, ssrindividualsout) =
//         fold(State_reduce_variety0, (rworld, srvarieties, ssrindividualsin, ssrindividualsout),
//              Indicies_all(srvarieties))
//
//   (rworld, srvarieties, ssrindividualsin, ssrindividualsout)
//
// Expressing the folds as loops this becomes the following code.
//
RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut State_reduce
(Space const space, World const world, SVarieties const svarieties, SSIndividuals const ssindividuals) {
  // Reduction setup
  RWorld rworld = RWorld_first(world);
  SRVarieties const srvarieties = SRVarieties_begin(world, svarieties);
  SSRIndividualsIn const ssrindividualsin = SSRIndividualsIn_begin(world, svarieties, ssindividuals);
  SSRIndividualsOut const ssrindividualsout = SSRIndividualsOut_begin(world, svarieties, ssindividuals);

  // Reduction outer loop over varieties
  UInt const varieties_number = ( svarieties->number <= ssindividuals->number ?
				  svarieties->number : ssindividuals->number );

  for (UInt varieties0_index=0; varieties0_index<varieties_number; varieties0_index+=1) {
    // Break out outer variety components
    Variety const variety0 = svarieties->variety[varieties0_index];
    SIndividuals const sindividuals0 = ssindividuals->sindividuals[varieties0_index];
    RVariety rvariety = srvarieties->rvariety[varieties0_index];
    SRIndividualsIn const srindividualsin = ssrindividualsin->srindividualsin[varieties0_index];

    // Outer loop over individuals
    for ( IIndividuals iindividuals0 = IIndividuals_first(sindividuals0);
          iindividuals0.valid;
          iindividuals0 = IIndividuals_next(iindividuals0) ) {
      // Break out outer individual components
      Individual const individual0 = IIndividuals_individual(iindividuals0);
      RIndividualIn rindividualin = srindividualsin->rindividualin[iindividuals0.index];

      // Outer level reductions
      rworld = RWorld_merge(rworld, RWorld_rest(world,variety0,individual0));
      rvariety = RVariety_merge(rvariety, RVariety_rest(world,variety0,individual0));

      // Reduction inner loop over varieties
      for (UInt varieties1_index=0; varieties1_index<varieties_number; varieties1_index+=1) {
        // Break out inner variety components
        Variety const variety1 = svarieties->variety[varieties1_index];
        SIndividuals const sindividuals1 = ssindividuals->sindividuals[varieties1_index];
        SRIndividualsOut const srindividualsout = ssrindividualsout->srindividualsout[varieties1_index];

        // Reduction inner loop over individuals in in loop
	Float32 box_in_ul_x,box_in_ul_y, box_in_lr_x,box_in_lr_y;
        unpack_Float32_Float32_Float32_Float32( &box_in_ul_x,&box_in_ul_y, &box_in_lr_x,&box_in_lr_y,
                                                RIndividualIn_bound(world, variety0,individual0, variety1) );

        UInt64 const box_in_ul_z = Z_xy(box_in_ul_x, box_in_ul_y, space.scale);
        UInt64 const box_in_lr_z = Z_xy(box_in_lr_x, box_in_lr_y, space.scale);

        for ( IIndividuals iindividuals1 = IIndividuals_firstBox(sindividuals1, box_in_ul_z,box_in_lr_z);
              iindividuals1.valid;
              iindividuals1 = IIndividuals_nextBox(iindividuals1, box_in_ul_z,box_in_lr_z) ) {
          // Break out inner individual components
          Individual const individual1 = IIndividuals_individual(iindividuals1);

          // Inner level reductions (individual to self is special case done in reduction setup)
          if ( (varieties0_index != varieties1_index || iindividuals0.index != iindividuals1.index) &&
               RIndividualIn_filter(world, variety0,individual0, variety1,individual1) )
            rindividualin = RIndividualIn_merge(rindividualin,
						RIndividualIn_rest(world, variety0,individual0,
								   variety1,individual1));
        }

        // Reduction inner loop over individuals in out loop
	Float32 box_out_ul_x,box_out_ul_y, box_out_lr_x,box_out_lr_y;
        unpack_Float32_Float32_Float32_Float32( &box_out_ul_x,&box_out_ul_y, &box_out_lr_x,&box_out_lr_y,
                                                RIndividualOut_bound(world, variety0,individual0, variety1) );

        UInt64 const box_out_ul_z = Z_xy(box_out_ul_x, box_out_ul_y, space.scale);
        UInt64 const box_out_lr_z = Z_xy(box_out_lr_x, box_out_lr_y, space.scale);

        for ( IIndividuals iindividuals1 = IIndividuals_firstBox(sindividuals1, box_out_ul_z,box_out_lr_z);
              iindividuals1.valid;
              iindividuals1 = IIndividuals_nextBox(iindividuals1, box_out_ul_z,box_out_lr_z) ) {
          // Break out inner individual level components
          Individual const individual1 = IIndividuals_individual(iindividuals1);
          RIndividualOut rindividualout = srindividualsout->rindividualout[iindividuals1.index];

          // Inner level reductions (individual to self is special case done in reduction setup)
          if ( (varieties0_index != varieties1_index || iindividuals0.index != iindividuals1.index) &&
               RIndividualOut_filter(world, variety0,individual0, variety1,individual1) )
            rindividualout = RIndividualOut_merge(rindividualout,
						  RIndividualOut_rest(world, variety0,individual0,
								      variety1,individual1));

          // Write back reduction
          srindividualsout->rindividualout[iindividuals1.index] = rindividualout;
        }
      }
      // Write back reduction
      srindividualsin->rindividualin[iindividuals0.index] = rindividualin;
    }
  }

  return pack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut(rworld, srvarieties,
                                                                    ssrindividualsin, ssrindividualsout);
}


// State_next: (Space, World, SVarieties, SSIndividuals,
//              RWorld, SRVarieties, SSRIndividualsIn, SSRIndividualsOut)
//               -> (World, SVarieties, SSIndividuals)
//
// Generate the next World, SVarieties, SSIndividuals from the current and reductions
//
//   World - (World, RWorld)
//   SVariety - (World, Variety, RWorld, RVariety)
//   SSIndividuals - (World, Variety, Individual, RWorld, RVariety, RIndividualIn, RIndividualOut)
//
// The new World and SVarieties are given by the *_next functions mapped across the tuples indicated above.  The
// new SSIndividuals is given by the Individual_next function folded across the tuple given above.
//
// The functional pseudo-code follows.
//
// State_next(world, svarieties, ssindividuals, rworld, srvarieties, ssrindividualsin, ssrindividualsout)
//   let // Next world
//       world_new = World_next(space, world, rworld);
//
//       // Individual (value individual, rindividual) construction code
//       State_next_individuals((aindividuals, variety, rvariety), (individual,rindividualin, rindividualout))
//         let // Next individual
//             aindividuals =
//               Individual_next(aindividuals, space, world, variety, individual,
//                               rspace, rvariety, rindividualin, rindividualout)
//         (aindividuals, variety, rvariety)
//
//       // Variety (value variety, rvariety) construction code
//       State_next_variety((), (variety, sindividuals, rvariety, srindividualsin, srindividualsout))
//         let // Next variety
//             variety_new = Variety_next(space, world, variety, rworld, rvariety)
//             // Individual (value individual, rindividual) construction fold
//             aindividuals = AIndividuals_begin()
//             aindividuals = fold(State_next_individuals,
//                                 (aindividuals, variety, rvariety),
//                                 zip(sindividuals, srindividualsin, srindividualsout))
//             sindividuals_new = AIndividuals_end(aindividuals)
//         (variety_new, sindividuals_new)
//
//       // Variety (value variety, rvariety) construction map
//       (svarieties, ssindividuals) =
//         unzip(map(State_next_variety, (), zip5(svarieties, ssindividuals,
//                                                srvarieties, ssrindividualsin, ssrindividualsout)))
//
//   (world_new, svarieties_new, SSIndividuals_new)
//
// Expressing the folds as loops this becomes the following code.
//
World_SVarieties_SSIndividuals State_next
(Space const space, World const world, SVarieties const svarieties, SSIndividuals const ssindividuals,
 RWorld const rworld, SRVarieties const srvarieties,
 SSRIndividualsIn const ssrindividualsin, SSRIndividualsOut const ssrindividualsout, Thread const thread) {
  // Next world
  World const world_new = World_next(space, world, rworld, thread);

  // Construction loop over varieties
  UInt varieties_number;
  AVarieties avarieties = AVarieties_begin();
  ASIndividuals asindividuals = ASIndividuals_begin();

  varieties_number = svarieties->number <= ssindividuals->number ? svarieties->number : ssindividuals->number;
  varieties_number = varieties_number <= srvarieties->number ? varieties_number : srvarieties->number;
  varieties_number = varieties_number <= ssrindividualsin->number ? varieties_number : ssrindividualsin->number;
  varieties_number = varieties_number <= ssrindividualsout->number ? varieties_number : ssrindividualsout->number;

  for (UInt varieties_index=0; varieties_index<varieties_number; varieties_index+=1) {
    // Break out variety components
    Variety const variety = svarieties->variety[varieties_index];
    SIndividuals const sindividuals = ssindividuals->sindividuals[varieties_index];
    RVariety const rvariety = srvarieties->rvariety[varieties_index];
    SRIndividualsIn const srindividualsin = ssrindividualsin->srindividualsin[varieties_index];
    SRIndividualsOut const srindividualsout = ssrindividualsout->srindividualsout[varieties_index];
    
    // Next variety
    Variety const variety_new = Variety_next(space, world, variety, rworld, rvariety, thread);
    avarieties = AVarieties_append(avarieties, variety_new);

    // Construction loop over individuals
    UInt individuals_number;
    AIndividuals aindividuals = AIndividuals_begin();

    individuals_number =
      sindividuals.number <= srindividualsin->number ? sindividuals.number : srindividualsin->number;
    individuals_number =
      individuals_number <= srindividualsout->number ? individuals_number : srindividualsout->number;

    for ( IIndividuals iindividuals = IIndividuals_first(sindividuals);
          iindividuals.valid && iindividuals.index < individuals_number;
          iindividuals = IIndividuals_next(iindividuals) ) {
      // Break out individuals components
      Individual const individual = IIndividuals_individual(iindividuals);
      RIndividualIn const rindividualin = srindividualsin->rindividualin[iindividuals.index];
      RIndividualOut const rindividualout = srindividualsout->rindividualout[iindividuals.index];

      // Next individual
      aindividuals = Individual_next(aindividuals, space, world, variety, individual,
                                     rworld, rvariety, rindividualin, rindividualout,
                                     thread);
    }

    SIndividuals const sindividuals_new = AIndividuals_end(aindividuals);
    asindividuals = ASIndividuals_append(asindividuals, sindividuals_new);
  }

  SVarieties const svarieties_new = AVarieties_end(avarieties);
  SSIndividuals const ssindividuals_new = ASIndividuals_end(asindividuals);

  return pack_World_SVarieties_SSIndividuals(world_new, svarieties_new, ssindividuals_new);
}


//---------------------------------------------------------------------------------------------------------------//

#endif // SYSTEM_CODE
