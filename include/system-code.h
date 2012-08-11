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
UInt64 Bits_oneIffLEMSB(UInt64 const bits);
UInt64 Bits_oneIffLEMSBZ(UInt64 const bits);
UInt64 Bits_zeroIffGTLSB(UInt64 const bits);

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
// Bits_oneIffLEMSB: (UInt64) -> (UInt64)
//
// Returned bits are one if and only if their position is less than or equal to that of the most-significant bit
// set (clear all bits if not bits set, otherwise set all bits up to and including last set bit).
//
UInt64 Bits_oneIffLEMSB(UInt64 const bits) {
  UInt64 reduction = bits;

  for (UInt shift = 1; shift < 64; shift *= 2)
    reduction = reduction | reduction >> shift;

  return reduction;
}


// Bits_oneIffLEMSBZ: (UInt64) -> (UInt64)
//
// Returned bits along each dimension are zero if and only if their position is less than or equal to that of the
// most-significant bit sets along that dimension (clear all bits if no bits set, otherwise set all bits up to
// and including last set bit along each dimension).
//
UInt64 Bits_oneIffLEMSBZ(UInt64 const bits) {
  UInt64 reduction = bits;

  for (UInt shift = 2; shift < 64; shift *= 2)
    reduction = reduction | reduction >> shift;

  return reduction;
}


// Bits_zeroIffGTLSB: (UInt64) -> (UInt64)
//
// Returned bits are zero if and only if their position is greater than that of the least-significant bit set
// (set all bits if no bits set, otherwise clear all bits above last set bit and set all bits below it).
//
UInt64 Bits_zeroIffGTLSB(UInt64 const bits) {
  return bits-1 ^ bits;
}


//---------------------------------------------------------------------------------------------------------------//
// Z_xy: (Float32,Float32, Float32) -> (UInt64)
//
// Put given x and y values in Z format (scale, change to integer, and interleave bits -- y_n,x_n,...,y_0,x_0).
//
// The interleaving in the routine works by separately interspersing the x and y components with zeros.  These
// intersepersed components are then offset and combined.
//
// The intersepersing is done via divide and conqueror.  At each step there is a number of bits to be
// interspersed with zeros.  These bits are divided into a top and bottom group which are themselves interspersed
// via recurision.  These interspersed bits are then combined in order to give the final result.
//
// As the recursive nature of the computation is deterministic, occurs on a bit level, and the combined storage
// requirements of all computations at each level is fixed, it can be expressed as manipulations combined
// independent groups of bits and be performed in parallel.  The details follow from this example
//
//   Z_DIMS = 3
//   Z_BITS = 8
//
//   bits =  8
//   mask =  111111111111 111111111111
//   z_x  =  000000000000 000010011101
//
//   bits =  4
//   mask = (111111111111 111111111111 ^ 111111111111 000000000000)
//        =  000000000000 111111111111
//   z_x  = (000000000000 000010011101 | 000000001001 110100000000) &
//          (000000000000 111111111111 ^ 000000001111 111111110000)
//        =  000000001001 000000001101)
//
//   bits =  2
//   mask = (000000000000 111111111111 ^ 000000111111 111111000000)
//        =  000000111111 000000111111
//   z_x  = (000000001001 000000001101 | 000010010000 000011010000) &
//          (000000111111 000000111111 ^ 000011111100 000011111100)
//        =  000010000001 000011000001
//
//   bits =  1
//   mask = (000000111111 000000111111 ^ 000111111000 000111111000)
//        =  000111000111 000111000111
//   z_x  = (000010000001 000011000001 | 001000000100 001100000100) &
//          (000111000111 000111000111 ^ 001110001110 001110001110)
//        =  001000000001 001001000001
//
// The functional pseudo-code follows.
//
// Z_DIMS = 2
// Z_BITS = 32
//
// Z_xy(x,y, scale)
//   let x = round(x*scale)
//       y = round(y*scale)
//   Z_xy_uint(x,y, 1<<Z_DIMS*Z_BITS-1, Z_BITS)
//
// Z_xy_uint(x,y, mask, bits)
//   // More bits to seperate, split them in half
//   bits > 1:
//     let bits = bits/2
//         mask = mask ^ mask<<Z_DIMS*bits
//         x = (x | x<<(Z_DIMS-1)*bits) & (mask ^ mask<<bits)
//         y = (y | y<<(Z_DIMS-1)*bits) & (mask ^ mask<<bits)
//         Z_xy_uint(x,y, bits, mask)
//   // No more bits to seperate, return results
//   otherwise:
//     y<<1 | x
//
UInt64 Z_xy(Float32 const x, Float32 const y, Float32 const scale) {
  Float32 const x_scaled = roundf(x*scale);
  Float32 const y_scaled = roundf(y*scale);
  UInt32 const x_uint = ( x_scaled < 0                      ? 0 :
                          x_scaled > nextafterf(0x1p32,0.0) ? nextafterf(0x1p32,0.0) :
                          floorf(x_scaled) );
  UInt32 const y_uint = ( y_scaled < 0                      ? 0 :
                          y_scaled > nextafterf(0x1p32,0.0) ? nextafterf(0x1p32,0.0) :
                          floorf(y_scaled) );

  UInt bits = 32;

  UInt64 mask = 0xffffffffffffffff;
  UInt64 z_x = x_uint;
  UInt64 z_y = y_uint;

  while (bits > 1) {
    bits /= 2;
    mask = mask ^ mask<<2*bits;
    z_x = (z_x | z_x<<(2-1)*bits) & (mask ^ mask<<bits);
    z_y = (z_y | z_y<<(2-1)*bits) & (mask ^ mask<<bits);
  }

  return z_y<<1 | z_x;
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
                  .scale = 0x1p32/fmaxf(size_x, size_y),
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
// IZ_iz: (IZ, Bool) -> (IZ)
//
// IZ with given values
//
IZ IZ_iz(bool valid, UInt64 z) {
  IZ const iz = { .valid = valid, .z = z };

  return iz;
}


// IZ_valid: (UInt64) -> (IZ)
//
// Valid IZ with given values
//
IZ IZ_valid(UInt64 const z) {
  return IZ_iz(true, z);
}


// IZ_invalid: () -> (IZ)
//
// Invalid IZ.
//
IZ IZ_invalid() {
  return IZ_iz(false, 0);
}


// IZ_zSet: (IZ, UInt64) -> (IZ)
//
// Duplicate of original IZ with Z set to z.
//
IZ IZ_zSet(IZ const iz, UInt64 const z) {
  return IZ_iz(iz.valid, z);
}


// IZ_nextBox: (IZ, UInt64,UInt64) -> (IZ)
//
// Smallest Z value greater than given Z value inside given box
//
// This routine works by computing a non-overlapping coverage of all subsequent Z values with completely-
// traversed hyperbox.  The next value is then the negative corner of the box clipped to the first of the
// hyperboxes that intersects the box.
//
// The clipped corner must be on the curve as it is in a completely-traversed hyperbox.  It must be the first
// Z value in the box to be traversed as all other values in the box either do not intersect or are greater.
//
// There is one completely-traversed hyperbox for each zero in the given Z value.  The Z values coverered by
// each of these hyperboxes is given by changing the zero into a one and iterating through all possible lower
// bit combinations.
//
// This works because it is precisely what adding ones to the given Z value does.  Eventually a carry sets the
// given zero bit.  Subsequent additions then iterate through all possible lower bit combinations.
//
// The actual computations are as follows.  Consider the given point and box corners (all in Z values) as a
// series of bits in the standard order.  Let these series be known as z and c0 and c1.  For each dimension i
// in the Z values, consider the additional series s0_i and s1_i (same) and o0_i and o1_i (outside).
//
// Let each bit in s0_i and s1_i indicate that c0 and c1, respectively, are the same as z from that bit upwards
// along the appropriate dimension.  Let each bit in o0_i and o1_i indicate that c0 and c1, respectively, are,
// greater and lesser, respectively, than z from that bit upwards (z falls outside the box).
//
// Using + to index the bit on the left and e_i as a mask for dimension i (kronecker delta), these can be written
// recursively as
//
//   ~s0_i = ~+s0_i | (~c0 & z | c0 & ~z) & e_i
//   ~s1_i = ~+s1_i | (~c1 & z | c1 & ~z) & e_i
//
//   o0_i = +o0_i | ~+s0_i &  c0 & ~z & e_i
//   o1_i = +o1_i | ~+s1_i & ~c1 &  z & e_i
//
// Consider an additional series p (possible).  Let each bit in p indicate that setting this bit in z and
// clearing all lower bits would give the negative corner of one of the completely traversed hyperboxes that
// intersects the box.  This can be expressed in terms of s0_i, s1_i, o0_i, and o1_i as
//
//   p = ~z & (\and_i ~+o0_i) & (\and_i ~+o1_i & (~+s1_i | c1 | ~e_i))
//     = ~z & ~(\or_i +o0_i | \or_i +o1_i | \or_i +s1_i & ~c1 & e_i)
//     = ~z & ~(\or_i +o_i | (\or_i +s1_i & e_i) & ~c1)
//     = ~z & ~(+o | ++s1 & ~c1)
//     = ~z & ~+o & (~++s1 | c1)
//
// where ++ indicates indexing the next bit in the same dimension on the left and
//
//   s0 = \or_i s0_i & e_i  (~s0 = ~++s0 | ~c0 & z | c0 & ~z)
//   s1 = \or_i s1_i & e_i  (~s1 = ~++s1 | ~c1 & z | c1 & ~z)
//
//   o_i = o0_i | o1_i
//
//   o0 = \or_i o0_i
//   o1 = \or_i o1_i
//
//   o = \or_i o_i = o0 | o1
//
// The combined recursive equations (those with the ++) come from expanding the non-combined ones once for each
// dimension.  The first hyperbox is indicated by the lowest set bit in p.  No set bits indicate there is no
// further Z values inside the box.
//
// Assume further Z values.  The next one has the same upper bits as z until the lowest set bit in p.  This bit
// is one and the remaning bits are either zero or c0 along each dimension.  The former happens when the negative
// corner of the box is clipped by the hyperbox in that dimension and the latter otherwise.
//
// Consider the series i0_i.  For each hyperbox identified by a set bit in p, let i0_i indicate if c0 is less
// than the negative corner of the hyperbox (c0 has to be clipped).  This can be expressed as
//
//   i0_i = ~s0_i & ~o0_i |_{z=1}
//        = (~+s0_i | (~c0 & z | c0 & ~z) & e_i) & ~(+o0_i | ~+s0_i & c0 & ~z & e_i) |_{z=1}
//        = (~+s0_i | ~c0 & e_i) & ~+o0_i
//        = ~+s0_i | ~c0 & e_i
//
// where the current bit of z is taken as set and p as true (p implying ~+o implying ~+o0_i).
//
// The actual hyperbox of interest is the first of those identified by p.  Let the series r indicate if a bit in
// this hyperbox is fixed at z (occurs at or before the first one in p)
//
//   r = -r & ~-p
//
// Let m0 be a series indicating if the corresponding bit should be zero or c0.  This is created by latching i0_i
// at the hyperbox (p & r)
//
//   m0_i = +m0_i | p & r & i0_i
//
//   m0 = \or_i m0_i & e_i
//      = ++m0 | (... | +p & +r | p & r) & ~++s0 | p & r & ~c0
//
// As before, the combined equation (the one with the ++) comes from expanding the non-combined one once for each
// dimension.
//
// The final solution is then z above r, one where r starts, and a clipped c0 below
//
//   z & ~r | p & r | ~m0 & c0 & +r
//
IZ IZ_nextBox(IZ const iz, UInt64 const corner0, UInt64 const corner1) {
  // Calculate overlap possibility at each level along each dimension
  UInt64 const not_same0 = Bits_oneIffLEMSBZ(corner0 ^ iz.z);
  UInt64 const not_same1 = Bits_oneIffLEMSBZ(corner1 ^ iz.z);

  // Calculate level at which overlap is no longer possible
  UInt64 const outside01 = Bits_oneIffLEMSB(~(not_same1>>2) & ~corner1 & iz.z | ~(not_same0>>2) & corner0 & ~iz.z);

  // Calculate bottom mask for overlap
  UInt64 const possible = ~iz.z & (corner1 | not_same1>>2) & ~(outside01>>1);
  UInt64 const region = Bits_zeroIffGTLSB(possible);

  // Calculate bottom inside masks along each dimension
  UInt64 const mask0 = Bits_oneIffLEMSBZ((region ^ region>>2) & not_same0>>2 | (region ^ region>>1) & ~corner0);

  // Combine top and bottom bits and return IZ
  return IZ_iz(possible && iz.valid,
               ~region & iz.z | region ^ (region>>1) | (region>>1) & ~mask0 & corner0);
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
    while ((index+1 & ((UInt64)1<<bit)) == 0) {
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
    if (((UInt64)1<<bit) < iindividuals.number-index-1 &&
        right_z[(index+((UInt64)1<<bit) >> bit/(64/DEPTH)*(64/DEPTH)) % CLUSTER] < iz.z)
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
SRVarieties SRVarieties_begin(Space const space, World const world, SVarieties const svarieties) {
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
    RVariety const rvariety = RVariety_first(space, world, variety);

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
// SRIndividualsIn_begin: (Space, World, Variety, SIndividuals) -> (SRIndividualsIn)
//
// SRIndividualsIn prepared for reduction of sindividuals (initialized with RIndividualIn_first).
//
SRIndividualsIn SRIndividualsIn_begin(Space const space, World const world, Variety const variety,
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
    RIndividualIn const rindividualin = RIndividualIn_first(space, world, variety, individual);

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
SRIndividualsOut SRIndividualsOut_begin(Space const space, World const world, Variety const variety,
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
    RIndividualOut const rindividualout = RIndividualOut_first(space, world, variety, individual);

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
// SSRIndividualsIn_begin: (Space, World, SVarieties) -> (SSRIndividualsIn)
//
// SSRIndividualsIn prepared for reduction of ssindividuals (initialized with RIndividualIn_first).
//
SSRIndividualsIn SSRIndividualsIn_begin(Space const space, World const world, SVarieties const svarieties,
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
    SRIndividualsIn const srindividualsin = SRIndividualsIn_begin(space, world, variety, sindividuals);

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
SSRIndividualsOut SSRIndividualsOut_begin(Space const space, World const world, SVarieties const svarieties,
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
    SRIndividualsOut const srindividualsout = SRIndividualsOut_begin(space, world, variety, sindividuals);

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
AVarieties_ASIndividuals pack_AVarieties_ASIndividuals(AVarieties const first, ASIndividuals const second) {
  AVarieties_ASIndividuals const value = { .first = first, .second = second };
  return value;
}

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
void unpack_AVarieties_ASIndividuals(AVarieties* const first, ASIndividuals* const second,
                                     AVarieties_ASIndividuals const tuple) {
  *first = tuple.first;
  *second = tuple.second;
}
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

  space.scale = 0x1p32/fmaxf(space.size_x,space.size_y);
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
//   RWorld            - fold (World, Variety, Individual)
//   SRVariety         - fold (World, Variety, Individual) grouped by (World, Variety)
//   SSRIndividualsIn  - fold (World, Variety[0|1], Individual[0|1]) grouped by (World, Variety0, Individual0)
//                       for each (Variety1, Individual1) not equal to and in in range of (Variety0, Individual0)
//   SSRIndividualsOut - fold (World, Variety[0|1], Individual[0|1]) grouped by (World, Variety1, Individual1)
//                       for each (Variety1, Individual1) not equal to and in out range of (Variety0, Individual0)
//
// The initial value of the reduction spaces are given by the R*_first functions mapped across the tuples
// indicates above.  These are then folded with the values given by R*_rest functions mapped across the the fold
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
// Box_clip:   ((Float32, Float32, Float32, Float32), (Float32, Float32, Float32, Float32))
//               -> (Float32, Float32, Float32, Float32)
// Box_offset: ((Float32, Float32, Float32, Float32), (Float32, Float32))
//               -> (Float32, Float32, Float32, Float32)
// Box_grid:   ((Bool, Bool), (Float32, Float32, Float32, Float32), (Float32, Float32))
//               -> (Int32, Int32, Int32, Int32)
//
// State_reduce(space, world, svarieties, ssindividuals)
//   let // Initial World reduction
//       rworld = RWorld_first(space, world)
//
//       // Initial Individual (value individual) reduction code
//       State_reduce_individual((variety), individual)
//         // RIndividualIn and RIndividualOut initialization
//         (RIndividualIn_first(space, world, variety, individual),
//          RIndividualOut_first(space, world, variety, individual))
//
//       // Initial Variety (value variety) reduction code
//       State_reduce_variety((), (variety, sindividuals))
//         let // RVariety initialization
//             rvariety = RVariety_first(space, world, variety)
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
//       State_reduce_individual1In((ssrindividualsin, x, y, v0, i0, v1), i1)
//         // First individual is not second individual and second is in in range, reduce second to first
//         (v0 /= v1 || i0 /= i1) && RIndividualIn_filter(x, y, world,
//                                                        svarieties[v0], ssindividuals[v0][i0],
//                                                        svarieties[v0], ssindividuals[v1][i1]):
//           // Inner level SSRIndividualsIn reduction
//           (ssrindividualsin[v0][i0] =
//              RIndividualIn_merge(ssrindividualsin[v0][i0],
//                                  RIndividualIn_rest(x, y, space, world,
//                                                     svarieties[v0], ssindividuals[v0][i0],
//                                                     svarieties[v0], ssindividuals[v1][i1])),
//            v0, i0, v1)
//         // First individual is second or second is not in in range, don't reduce second to first
//         otherwise:
//           (ssrindividualin, v0, i0, v1)
//
//       // Inner Individual (index i1) out reduction code
//       State_reduce_individual1Out((ssrindividualsout, x, y, v0, i0, v1), i1)
//         // First individual in not second individual and second is in out range, reduce first to second
//         (v0 /= v1 || i0 /= i1) && RIndividualOut_filter(world,
//                                                         svarieties[v0], ssindividuals[v0][i0],
//                                                         svarieties[v0], ssindividuals[v1][i1]):
//           // Inner level SSRIndividualsOut reduction
//           (ssrindividualsout[v1][i1] =
//              RIndividualOut_merge(ssrindividualsout[v1][i1],
//                                   RIndividualOut_rest(x, y, space, world,
//                                                       svarieties[v0], ssindividuals[v0][i0],
//                                                       svarieties[v1], ssindividuals[v1][i1])),
//            v0, i0, v1)
//         // First individual is second or second is not in our range, don't reduce first to second
//         otherwise:
//           ssrindividualsout
//
//       // Periodic in reduction code
//       State_reduce_periodicIn((ssrindividualsin, v0, i0, v1), (x,y))
//         let // Inner Individual (index i1) reduction fold
//             bound = Box_clip(Box_offset(RIndividualIn_bound(x, y, world, svarieties[v0],
//                                                             ssindividuals[v0][i0], svarieties[v1]),
//                                         (x*space.size_x, y*space.size_y))
//                              (0, 0, space.size_x, space.size_y))
//             ssrindividualsin = fold(State_reduce_individual1In, (ssrindividualsin, x, y, v0, i0, v1),
//                                     Indicies_filterBox(ssindividuals[v1], bound))
//         (ssrindividualsin, v0, i0, v1)
//
//       // Periodic out reduction code
//       State_reduce_periodicOut((ssrindividualsout, v0, i0, v1), (x,y))
//         let // Inner Individual (index i1) reduction fold
//             bound = Box_clip(Box_offset(RIndividualOut_bound(x, y, world, svarieties[v0],
//                                                              ssindividuals[v0][i0], svarieties[v1]),
//                                         (x*space.size_x, y*space.size_y)),
//                              (0, 0, space.size_x, space.size_y))
//             ssrindividualsout =  fold(State_reduce_individual1Out, (ssrindividualsout, x, y, v0, i0, v1),
//                                       Indicies_filterBox(ssindividuals[v1], bound))
//         (ssrindividualsout, v0, i0, v1)
//
//       // Inner Variety (index v1) reduction code
//       State_reduce_variety1((ssrindividualsin, ssrindividualsout, v0, i0), v1)
//         let // Periodic in reduction code
//             bound = RIndividualIn_bound(world, svarieties[v0], ssindividuals[v0][i0], svarieties[v1])
//             indicies = Indicies_box((space.periodic_x, space.periodic_y), bound, (space.size_x, space.size_y))
//             ssrindividualsin = fold(State_reduce_periodicIn, (ssrindividualsin, v0, i0, v1), indicies)
//             // Periodic out out reduction code
//             bound = RIndividualOut_bound(world, svarieties[v0], ssindividuals[v0][i0], svarieties[v1])
//             indicies = Indicies_box((space.periodic_x, space.periodic_y), bound, (space.size_x, space.size_y))
//             ssrindividualsout = fold(State_reduce_periodicOut, (ssrindividualsout, v0, i0, v1), indicies)
//         (ssrindividualsin, ssrindividualsout, v0, i0)
//
//       // Outer Individual (index i0) reduction code
//       State_reduce_individual0((rworld, srvarieties, ssrindividualsin, ssrindividualsout, v0), i0)
//         let // Outer level RWorld and SRVariety reductions
//             rworld = RWorld_merge(rworld,
//                                   RWorld_rest(space, world, svarieties[v0], ssindividuals[v0][i0]))
//             srvarieties[v0] =
//               RVariety_merge(srvarieties[v0],
//                              RVariety_rest(space, world, svarieties[v0], ssindividuals[v0][i0]))
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
// Box_clip((box0_ul_x, box0_ul_y, box0_lr_x, box0_lr_y), (box1_ul_x, box1_ul_y, box1_lr_x, box1_lr_y))
//   // Intersection between box0 and box1 (assumes overlap)
//   (max(box0_ul_x, box1_ul_x), max(box0_ul_y, box1_ul_y),
//    min(box0_lr_x, box1_lr_x), min(box0_lr_y, box1_lr_y))
//
// Box_offset((box_ul_x, box_ul_y, box_lr_x, box_lr_y), (offset_x, offset_y))
//   // Shift box by offset
//   (box_ul_x+offset_x, box_ul_y+offset_y, box_lr_x+offset_x, box_lr_y+offset_y)
//
// Box_grid((periodic_x,periodic_y), (box_ul_x, box_ul_y, box_lr_x, box_lr_y), (size_x, size_y))
//   // Grid coverage of box, where grid components are of size and optionally periodic
//   (if(periodic_x, floor(box_ul_x/size_x), 0), if(periodic_y, floor(box_ul_y/size_y), 0),
//    if(periodic_x, floor(box_lr_x/size_y), 0), if(periodic_y, floor(box_lr_y/size_y), 0))
//
// Expressing the folds as loops this becomes the following code.
//
RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut State_reduce
(Space const space, World const world, SVarieties const svarieties, SSIndividuals const ssindividuals) {
  // Reduction setup
  RWorld rworld = RWorld_first(space, world);
  SRVarieties const srvarieties = SRVarieties_begin(space, world, svarieties);
  SSRIndividualsIn const ssrindividualsin = SSRIndividualsIn_begin(space, world, svarieties, ssindividuals);
  SSRIndividualsOut const ssrindividualsout = SSRIndividualsOut_begin(space, world, svarieties, ssindividuals);

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
      rworld = RWorld_merge(rworld, RWorld_rest(space, world,variety0,individual0));
      rvariety = RVariety_merge(rvariety, RVariety_rest(space, world,variety0,individual0));

      // Reduction inner loop over varieties
      for (UInt varieties1_index=0; varieties1_index<varieties_number; varieties1_index+=1) {
        // Break out inner variety components
        Variety const variety1 = svarieties->variety[varieties1_index];
        SIndividuals const sindividuals1 = ssindividuals->sindividuals[varieties1_index];
        SRIndividualsOut const srindividualsout = ssrindividualsout->srindividualsout[varieties1_index];

        // Periodic loops for individuals in in loop
        Float32 bound_in_ul_x,bound_in_ul_y, bound_in_lr_x,bound_in_lr_y;
        unpack_Float32_Float32_Float32_Float32
          ( &bound_in_ul_x,&bound_in_ul_y, &bound_in_lr_x,&bound_in_lr_y,
            RIndividualIn_bound(space, world, variety0,individual0, variety1) );

        Int const mirror_in_ul_x = space.periodic_x ? floorf(bound_in_ul_x/space.size_x) : 0;
        Int const mirror_in_ul_y = space.periodic_y ? floorf(bound_in_ul_y/space.size_y) : 0;
        Int const mirror_in_lr_x = space.periodic_x ? floorf(bound_in_lr_x/space.size_x) : 0;
        Int const mirror_in_lr_y = space.periodic_y ? floorf(bound_in_lr_y/space.size_y) : 0;

        for (Int mirror_in_x = mirror_in_ul_x; mirror_in_x <= mirror_in_lr_x; mirror_in_x +=1 ) {
          for (Int mirror_in_y = mirror_in_ul_y; mirror_in_y <= mirror_in_lr_y; mirror_in_y += 1 ) {
            // Reduction inner loop over individuals in in loop in current mirror
            UInt64 const bound_in_ul_z = Z_xy(fmaxf(0.0, bound_in_ul_x - mirror_in_x*space.size_x),
                                              fmaxf(0.0, bound_in_ul_y - mirror_in_y*space.size_y),
                                              space.scale);
            UInt64 const bound_in_lr_z = Z_xy(fminf(bound_in_lr_x - mirror_in_x*space.size_x, space.size_x),
                                              fminf(bound_in_lr_y - mirror_in_y*space.size_y, space.size_y),
                                              space.scale);

            for ( IIndividuals iindividuals1 = IIndividuals_firstBox(sindividuals1, bound_in_ul_z,bound_in_lr_z);
                  iindividuals1.valid;
                  iindividuals1 = IIndividuals_nextBox(iindividuals1, bound_in_ul_z,bound_in_lr_z) ) {
              // Break out inner individual components
              Individual const individual1 = IIndividuals_individual(iindividuals1);

              // Inner level reductions (individual to self is special case done in reduction setup)
              if ( (varieties0_index != varieties1_index || iindividuals0.index != iindividuals1.index ||
                    mirror_in_x != 0 || mirror_in_y != 0) &&
                   RIndividualIn_filter(mirror_in_x, mirror_in_y, space, world,
                                        variety0,individual0, variety1,individual1) )
                rindividualin = RIndividualIn_merge(rindividualin,
                                                    RIndividualIn_rest(mirror_in_x, mirror_in_y, space, world,
                                                                       variety0,individual0,
                                                                       variety1,individual1));
            }
          }
        }

        // Periodic loops for for individuals in out loop
        Float32 bound_out_ul_x,bound_out_ul_y, bound_out_lr_x,bound_out_lr_y;
        unpack_Float32_Float32_Float32_Float32
          ( &bound_out_ul_x,&bound_out_ul_y, &bound_out_lr_x,&bound_out_lr_y,
            RIndividualOut_bound(space, world, variety0,individual0, variety1) );

        Int const mirror_out_ul_x = space.periodic_x ? floorf(bound_out_ul_x/space.size_x) : 0;
        Int const mirror_out_ul_y = space.periodic_y ? floorf(bound_out_ul_y/space.size_y) : 0;
        Int const mirror_out_lr_x = space.periodic_x ? floorf(bound_out_lr_x/space.size_x) : 0;
        Int const mirror_out_lr_y = space.periodic_y ? floorf(bound_out_lr_y/space.size_y) : 0;

        for (Int mirror_out_x = mirror_out_ul_x; mirror_out_x <= mirror_out_lr_x; mirror_out_x +=1 ) {
          for (Int mirror_out_y = mirror_out_ul_y; mirror_out_y <= mirror_out_lr_y; mirror_out_y += 1 ) {
            // Reduction inner loop over individuals in out loop
            UInt64 const bound_out_ul_z = Z_xy(fmaxf(0.0, bound_out_ul_x - mirror_out_x*space.size_x),
                                               fmaxf(0.0, bound_out_ul_y - mirror_out_y*space.size_y),
                                               space.scale);
            UInt64 const bound_out_lr_z = Z_xy(fminf(bound_out_lr_x - mirror_out_x*space.size_x, space.size_x),
                                               fminf(bound_out_lr_y - mirror_out_y*space.size_y, space.size_y),
                                               space.scale);

            for ( IIndividuals iindividuals1 = IIndividuals_firstBox(sindividuals1, bound_out_ul_z,bound_out_lr_z);
                  iindividuals1.valid;
                  iindividuals1 = IIndividuals_nextBox(iindividuals1, bound_out_ul_z,bound_out_lr_z) ) {
              // Break out inner individual level components
              Individual const individual1 = IIndividuals_individual(iindividuals1);
              RIndividualOut rindividualout = srindividualsout->rindividualout[iindividuals1.index];

              // Inner level reductions (individual to self is special case done in reduction setup)
              if ( (varieties0_index != varieties1_index || iindividuals0.index != iindividuals1.index ||
                    mirror_out_x != 0 || mirror_out_y != 0) &&
                   RIndividualOut_filter(mirror_out_x, mirror_out_y, space, world,
                                         variety0,individual0, variety1,individual1) )
                rindividualout = RIndividualOut_merge(rindividualout,
                                                      RIndividualOut_rest(mirror_out_x, mirror_out_y, space, world,
                                                                          variety0,individual0,
                                                                          variety1,individual1));

              // Write back reduction
              srindividualsout->rindividualout[iindividuals1.index] = rindividualout;
            }
          }
        }
      }
      // Write back reduction
      srindividualsin->rindividualin[iindividuals0.index] = rindividualin;
    }
    // Write back reduction
    srvarieties->rvariety[varieties0_index] = rvariety;
  }

  return pack_RWorld_SRVarieties_SSRIndividualsIn_SSRIndividualsOut(rworld, srvarieties,
                                                                    ssrindividualsin, ssrindividualsout);
}


// State_next: (Space, World, SVarieties, SSIndividuals,
//              RWorld, SRVarieties, SSRIndividualsIn, SSRIndividualsOut)
//               -> (World, SVarieties, SSIndividuals)
//
// Compute the next state by successively replacing/updating the data structure from the outside in.
//
//   (Individual)                     - map (one-to-many) ([R]World, [R]Variety, [R]Individual[In|Out])
//   (Variety,SIndividuals)           - map (one-to-many) ([R]World, [R]Variety) at the (Variety) level
//   (World,SVarieties,SSIndividuals) - map (one-to-one) ([R]World)
//
// The *_next functions are mapped across the tuples given above.  For the one-to-many cases, the results are
// joined/concatenated together.  This allows each of these components to be replaced by an arbitrary number (or
// none) of new ones.  The mapping and joining/concatenating is done from the outside in (Individual, Variety,
// and finally World) and the prior results are used in each new state.  The outside in ordering is required to
// maintain the applicability of the reduction data.  If the inside was mapped first, the outside reduction data
// may no longer map to the new structure.  The updating is required for the outside updates to have any effect.
//
// The problem this update approach addresses is how to abstract the replacement of as much or as little of the
// entire underlying data structure in as easy a way as possible.  The fundamental issue is that only the edges
// of a data structure can be replaced in an arbitrary fashion because any update to internal components are
// contrained to having to maintain the structure decending data is rooted in.  The solution taken here is to
// only update edges by extending the replacement at all stages to reach the edges.
//
// For example, a Variety cannot be replaced with an arbitrary number of other Variety without first addressing
// what will happens to the SIndividuals rooted in the replaced Variety and what SIndividuals will be rooted in
// the new Variety.  How can this be done in a generic fashion?  Include the SIndividuals in the replacement.
// (Variety,SIndividuals) tuples can be replaced by an arbitrary number of other (Variety, SIndividuals) tuples
// because the SIndividuals rooted in the Variety have now been packaged together as a whole.  This effectively
// punts the SIndividuals questions into the replacement routine.  Exactly where it has to be to be generic.
//
// The functional pseudo-code follows.
//
// State_next(world, svarieties, ssindividuals, rworld, srvarieties, ssrindividualsin, ssrindividualsout)
//   let // Individual next construction code (one-to-many)
//       State_next_individuals((world, variety, rvariety, rworld), (individual,rindividualin, rindividualout))
//         let // Next sindividuals (one-to-many)
//             sindividuals = Individual_next(space, world, variety, individual,
//                                            rworld, rvariety, rindividualin, rindividualout)
//         sindividuals
//
//       // Variety next construction code (one-to-many)
//       State_next_variety((world, rworld), (variety, sindividuals, rvariety, srindividualsin, srindividualsout))
//         let // Next sindividuals (many-to-many)
//             sindividuals = join(map(State_next_individuals, (world, variety, rworld, rvariety),
//                                     zip3(sindividuals, srindividualsin, srindividualsout)))
//             // Next svarieties (one-to-many)
//             (svarieties, ssindividuals) = Variety_next(space, world, variety, sindividuals, rworld, rvariety)
//         (svarieties, ssindividuals)
//
//       // World next construction code (one-to-one)
//       State_next_world((), (world, svariety, ssindividuals,
//                             rworld, srvariety, ssrindividualsin, ssrindividualsout))
//         let // Next svarieties (many-to-many)
//             (svarieties,ssindividuals) unzip(join(map(State_next_variety, (world, rworld),
//                                                   zip5(svarieties, ssindividuals,
//                                                        srvarieties, ssrindividualsin, ssrindividualsout))))
//             // Next world (one-to-one)
//             (world,svarieties,ssindividuals) = World_next(space, world, svarieties, ssindividuals, rworld)
//         (world,svarieties,ssindividuals)
//
//   State_next_world((), (world, svariety, ssindividuals,
//                         rworld, srvariety, ssrindividualsin, ssrindividualsout))
//
// Expressing the map as loops this becomes the following code.
//
World_SVarieties_SSIndividuals State_next
(Space const space, World const world, SVarieties const svarieties, SSIndividuals const ssindividuals,
 RWorld const rworld, SRVarieties const srvarieties,
 SSRIndividualsIn const ssrindividualsin, SSRIndividualsOut const ssrindividualsout, Thread const thread) {
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
    
    // Construction loop over individuals
    AIndividuals aindividuals = AIndividuals_begin();

    UInt individuals_number;
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

      // Next individuals
      aindividuals = Individual_next(aindividuals, space, world, variety, individual,
                                     rworld, rvariety, rindividualin, rindividualout,
                                     thread);
    }

    // Next varieties
    unpack_AVarieties_ASIndividuals(&avarieties, &asindividuals,
                                    Variety_next(avarieties, asindividuals, space, world, variety, aindividuals,
                                                 rworld, rvariety, thread));
  }

  // Next world
  return World_next(space, world, avarieties, asindividuals, rworld, thread);
}


//---------------------------------------------------------------------------------------------------------------//

#endif // SYSTEM_CODE
