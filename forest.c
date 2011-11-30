#include <stdint.h>

// Split 64b Z space into components (require CLUSTER=2^N for some N and CLUSTER^DEPTH = 2^64)

#define DEPTH  4
#define CLUSTER 65536

// Buffer size for fatal errors

#define ERRNO_BUFFER 1024

// Z values
//
// Z values form a fractal, which each set of two bits giving a layer of 2x2 connected boxes connected in a Z
// (upper left, upper right, lower left, lower right).  Each layer connects the 2x2 connected boxes from the lower
// layer into a new set of 2x2 connected boxes.  The lowest layer is formed by the two least significant bits.
//
// The X & Y values contained in the Z values can easily be compared by masking.

#define X_MASK 0x5555555555555555u
#define Y_MASK 0xaaaaaaaaaaaaaaaau


//---------------------------------------------------------------------------------------------------------------//
// Single upper case to distinguish types
typedef unsigned int UInt;
typedef uint32_t UInt32;
typedef uint64_t UInt64;

typedef int Int;
typedef int32_t Int32;
typedef int64_t Int64;

typedef struct _Forest Forest;

typedef struct _Varieties Varieties;
typedef struct _Variety Variety;

typedef struct _Individuals Individuals;
typedef union _Level Level;
typedef struct _Level0 Level0;
typedef struct _Level1 Level1;
typedef struct _Individual Individual;


//---------------------------------------------------------------------------------------------------------------//
// Individuals
struct _Individual {
  float x;
  float y;
  float height;
};

union _Level {
  Level0* level0;
  Level1* level1;
};

struct _Level1 {
  UInt64 z[CLUSTER];
  Individual individual[CLUSTER];
};

struct _Level0 {
  UInt64 left_z[CLUSTER];
  UInt64 right_z[CLUSTER];
  Level level[CLUSTER];
};

struct _Individuals {
  UInt64 number;
  Level level;
};

// Varieties (individuals grouped by model parameters)
struct _Variety {
  float height_mature;
  float height_maximum;

  float growth_rate;
  float growth_competition_lower;
  float growth_competition_higher;

  float mortality_initial;
  float mortality_decay;
  float mortality_intrinsic;

  float fecundity_maximum;

  float masting_time;
  float masting_phase;

  float dispersal_probability_short;
  float dispersal_mode_short;
  float dispersal_mode_long;

  Individuals individuals;
};

struct _Varieties {
  UInt number;
  Variety* varieties;
};

// Forest (all individuals and global model parameters)
struct _Forest {
  float scale;

  bool periodic_x;
  bool periodic_y;

  float size_x;
  float size_y;

  Varieties varieties;
};


//---------------------------------------------------------------------------------------------------------------//
// z_xy: (UInt32, UInt32) -> (UInt64)
//
// Put given X and Y values in Z format (interleave the X and Y bits -- y_n,x_n,...,y_0,x_0).
//
// The functional pseudo-code follows.
//
// z_xy(x,y)
//   // Non-zero x or y, directly handle bottom x and y bits and obtain rest by recursion
//   x > 0 || y > 0:
//     z_xy(x/2,y/2)*4 + (y%2)*2 + (x%2)
//   // Zero x and y, done
//   otherwise:
//     0
//
// Using this code the generate direct conversion tables for 4bit values, and then directly converting 32bit
// values by treating them as a sequence of eight 4bit values gives the following code.
//
UInt64 z_xy(UInt32 x, UInt32 y) {
  static const uint8_t z_table[] = {
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
  uint8_t part[8];

  part[0] = z_table[(uint8_t)(y >>  0 & 0x0f) << 4 | (uint8_t)(x >>  0 & 0x0f)];
  part[1] = z_table[(uint8_t)(y >>  4 & 0x0f) << 4 | (uint8_t)(x >>  4 & 0x0f)];
  part[2] = z_table[(uint8_t)(y >>  8 & 0x0f) << 4 | (uint8_t)(x >>  8 & 0x0f)];
  part[3] = z_table[(uint8_t)(y >> 12 & 0x0f) << 4 | (uint8_t)(x >> 12 & 0x0f)];
  part[4] = z_table[(uint8_t)(y >> 16 & 0x0f) << 4 | (uint8_t)(x >> 16 & 0x0f)];
  part[5] = z_table[(uint8_t)(y >> 20 & 0x0f) << 4 | (uint8_t)(x >> 20 & 0x0f)];
  part[6] = z_table[(uint8_t)(y >> 24 & 0x0f) << 4 | (uint8_t)(x >> 24 & 0x0f)];
  part[7] = z_table[(uint8_t)(y >> 28 & 0x0f) << 4 | (uint8_t)(x >> 28 & 0x0f)];

  return ( (uint64_t)(part[0]) <<  0 | (uint64_t)(part[1]) <<  8 |
           (uint64_t)(part[2]) << 16 | (uint64_t)(part[3]) << 24 |
           (uint64_t)(part[4]) << 32 | (uint64_t)(part[5]) << 40 |
           (uint64_t)(part[6]) << 48 | (uint64_t)(part[7]) << 56 );
}


//---------------------------------------------------------------------------------------------------------------//
// indices_reverse: (UInt64, UInt) -> (UInt64)
//
// For given index up to given depth, reverse the entries.
//
// This routine works by treating the index as a collection of CLUSTER sized indicies.  Shifting the indices off
// of one integer and onto another reverses their order (same as reversing a list).
//
// The functional pseudo-code follows.
//
// indices_reverse: (UInt64, UInt) -> (UInt64)
// indices_reverse_shift: (UInt64, UInt64, UInt) -> (UInt64)
//
// indices_reverse(indices, depth)
//   // Start with everything in the forward indices and nothing in the reverse
//   indices_reverse_shift(indicies, 0, depth)
//
// indices_reverse_shift(forward, reverse, number)
//   // More chunks, shift next one from forward indices to reverse and recurse to handle rest
//   number > 0:
//     indices_reverse_shift(forward/CLUSTER, reverse*CLUSTER+forward%CLUSTER, number-1)
//   // No more chunks, done
//   otherwise:
//     indices_reverse
//
// Flattening the recursive calls with loops this becomes the following code.
//
UInt64 indices_reverse(UInt64 indices, UInt number) {
  // For each level, shift the indices off the end of the one integer and onto the other
  UInt64 forward = indices;
  UInt64 reverse = 0;

  for (UInt iterator=number; iterator>0; iterator-=1) {
    reverse = reverse*CLUSTER + forward%CLUSTER;
    forward /= CLUSTER;
  }

  return reverse;
}


// individuals_level: (Individuals, UInt64, UInt) -> (Level)
//
// For given index up to given depth, return the level.
//
// The functional pseudo-code follows.
//
// individuals_level: (Individuals, UInt64, UInt) -> (Level)
// individuals_level_walk: (Level, UInt64, UInt) -> (Level)
//
// individuals_level(individuals, index, depth)
//   // Start at top level with reversed index as a list of individual level indices in required order
//   individuals_level_walk(level(individuals), indices_reverse(index), depth)
//
// individuals_level_walk(level, indices, depth)
//   // Not at desired depth, use next porition of reversed indices to index next level and recurse
//   depth > 0:
//     individuals_level_walk(level(level0(level))[indices%CLUSTER], indices/CLUSTER, depth-1)
//   // At desire depth, done
//   otherwise:
//     level
//
// Flattening the recursive calls with loops this becomes the following code.
//
Level individuals_level(Individuals* individuals, UInt64 index, UInt depth) {
  // Traverse each level up to the requested depth
  UInt64 indices = indices_reverse(index, depth+1);
  Level level = individuals->level;

  for (UInt depth_iterator=depth; depth_iterator>0; depth_iterator-=1) {
    level = level.level0->level[indices%CLUSTER];
    indices /= CLUSTER;
  }

  return level;
}


// individuals_level1: (Individuals, UInt64) -> (Level1)
//
// For given index, return the level1.
//
// The functional pseudo-code follows.
//
// individuals_level1: (Individuals, UInt64) -> (Level1)
//
// individuals_level1(individuals, index)
//   // Use general routine to get level1
//   level1(individuals_level(individuals, index, DEPTH-1))
//
// This becomes the following code.
//
Level1* individuals_level1(Individuals* individuals, UInt64 index) {
  // Use general routine to get level1
  return individuals_level(individuals, index, DEPTH-1).level1;
}


// individuals_append: (Individuals) -> (Level1)
//
// Add a new Level1 onto end of individuals (assumes number%CLUSTER == 0).
//
// The functional pseudo-code follows.
//
// individuals_append: (Individuals) -> (Individuals)
// individuals_existing: (Level, UIn64, UInt) -> (Level)
// individuals_add: (UInt) -> (Level)
//
// individuals_append(individuals)
//   let indices = indices_reverse(number(individuals), DEPTH)
//   // Entries below this one, descend into them
//   indices > 0:
//     individuals { .level = individuals_existing(level(individuals)[indices%CLUSTER], indices/CLUSTER, DEPTH-1) }
//   // No entries below this one, add them
//   otherwise:
//     individuals { .level = individuals_add(DEPTH) }
//
// individuals_walk(level, indices, depth)
//   // Entries below this one, descend into them
//   indices > 0:
//     individuals_existing(level(level0(level))[indicies%CLUSTER], indicies/CLUSTER, depth-1)
//   // No entries below this one, add them
//   otherwise:
//     individuals_add(depth)
//
// individuals_add(depth)
//   // Level0 depth, add Level0 and continue down
//   depth > 0:
//     level(level0() { .level[0] = individuals_add(depth-1) })
//   // Level1 depth, add Level1 and done
//   otherwise:
//     level(level1())
//
// Flattening the recursive calls with loops this becomes the following code.
//
Level1* individuals_append(Individuals* individuals) {
  UInt64 indices = indices_reverse(individuals->number, DEPTH);
  Level* level = &individuals->level;
  UInt depth = DEPTH;

  // While existing intermediate Level0 entries, descend into them
  while (indices > 0) {
    level = &level->level0->level[indices%CLUSTER];
    indices /= CLUSTER;
    depth -= 1;
  }

  // While needing intermediate Level0 entries, add them
  while (depth > 1) {
    if ( !(level->level0 = malloc(sizeof(Level0))) )
      die_errno(1, "unable to allocate %tu bytes for new Level0", sizeof(Level0));
    level = &level->level0->level[indices%CLUSTER];
    indices /= CLUSTER;
    depth -= 1;
  }

  // Add final Level1 entry
  if ( !(level->level1 = malloc(sizeof(Level1))) )
    die_errno(1, "unable to allocate %tu bytes for new Level1", sizeof(Level1));

  return level->level1;
}


// individuals_z: (Individuals, UInt64) -> (UInt64)
//
// For given index, return the Z value.
//
// The functional pseudo-code follows.
//
// individuals_z: (Individuals, UInt64) -> (UInt64)
//
// individuals_z(individuals, index)
//   z(individuals_level1(individuals, index))[index%CLUSTER]
//
// This becomes the following code.
//
UInt64 individuals_z(Individuals* individuals, UInt64 index) {
  return individuals_level1(individuals, index)->z[index%CLUSTER];
}


// individuals_left_z: (Individuals, UInt64, UInt) -> (UInt64)
//
// For given index up to given depth, return the left Z bound array.
//
// The functional pseudo-code follows.
//
// individuals_left_z: (Individuals, UInt64, UInt) -> (UInt64)
//
// individuals_left_z(individuals, index, depth)
//   // Given depth is Level0
//   depth > DEPTH-1:
//     z_left(level0(individuals_level(individuals, index, depth)))
//   // Given depth is Level1
//   otherwise:
//     z(level1(individuals_level(individuals, index, depth)))
//
// This becomes the following code
//
UInt64* individuals_left_z(Individuals* individuals, UInt64 index, UInt depth) {
  // If depth is a level0 level, the left Z bound array is left_z
  if (depth < DEPTH-1)
    return individuals_level(individuals, index, depth).level0->left_z;
  // If depth is a level1 level, the left Z bound array is z (the actual Z values)
  else
    return individuals_level(individuals, index, depth).level1->z;
}


//---------------------------------------------------------------------------------------------------------------//
// individual_left: (Individuals, UInt64, UInt64) -> ()
//
// Individual with greatest index and Z value less than or equal to given index and Z value.
//
// This routine works by treating the given index and an upper bound and reducing it to the point the desired
// index (the greatest one identifying a individual with a lesser or equal Z value to the given one) is
// identified.
//
// The given index is reduced by checking the Z range in the largest bit-aligned index-range immediately
// proceeding it for overlap with the given Z value.  If no overlap is found, the index range is skipped by
// reducing the given index to the start of it and repeating the procedure.  If the index reaches zero, then
// there is no index with a less than or equal Z value.
//
// Once overlap is found, the overlapping region is broken down into sub regions.  Each sub region without
// overlap is skipped by reducing the given index to the start of it.  Once an overlapping region is found, the
// procedure is repeated.
//
// The functional pseudo-code follows.
//
// individual_left_out:      (UInt64, UInt64, UInt64) -> (UInt64)
// individual_left_out_size: (UInt64, UInt64, UInt64) -> (UInt64)
// individual_left_out_step: (UInt64, UInt64, UInt64) -> (UInt64)
//
// individual_left_in_size: (UInt64, UInt64, UInt64) -> (UInt64)
// individual_left_in_step: (UInt64, UInt64, UInt64) -> (UInt64)
//
// individual_left_out(index, step, z)
//   // At first individual, done
//   index == 0:
//     fail("no more valid individuals")
//   // Not at first individual, determine step size
//   otherwise:
//     individual_left_out_size(index, step, z)
//
// individual_left_out_size(index, step, z)
//   // On step size boundary, increase step size
//   index % (step*2) == 0:
//     individual_left_out_size(index, step*2, z)
//   // Not on step size boundary, check interval to next boundary
//   otherwise:
//     individual_left_out_step(index, step, z)
//
// individual_left_out_step(index, step, z)
//   // Previous interval contains individual, locate individual in it
//   individual[index-step] <= z:
//     individual_left_in_size(index, step, z)
//   // Previous interval doesn't contain individual, skip over it
//   otherwise:
//     individual_left_out(index-step, step, z)
//
// individual_left_in_size(index, step, z)
//   // Single point, done
//   step == 1:
//     index-1
//   // Not single point, divide into sub regions
//   otherwise:
//     individual_left_in_step(index, step/2, z)
//
// individual_left_in_step(index, step, z)
//   // Previous sub region contains individual, locate z point in it
//   individual[index-step] <= z:
//     individual_left_in_size(index, step, z)
//   // Previous sub region doesn't contain individual, skip over it
//   otherwise:
//     individual_left_in_step(index-step, step, z)
//
// Flattening the recursive calls with loops this becomes the following code.
//
void individual_left(Individuals* individuals, UInt64 index, UInt64 z) {
  UInt64* left_z = individuals_left_z(individuals, index-1, DEPTH-1);
  UInt64 step = 1;
  UInt depth = DEPTH-1;

  // Find interval containing suitable individuals
  while (1) {
    // If at first individual, done
    if (index == 0) {
      abort(); // ...
    }

    // While on step size boundary, increase step size
    while ((index & (step*2-1)) == 0)
      step *= 2;

    if (step >= CLUSTER) {
      do {
        index /= CLUSTER;
        step /= CLUSTER;
        depth -= 1;
      } while (step >= CLUSTER);

      left_z = individuals_left_z(individuals, index-1, depth);
    }

    // If previous interval contains individual, locate individual in it
    if (left_z[(index-step)%CLUSTER] <= z)
      break;

    // Previous interval doesn't contain individual, skip over it
    index -= step;
  }

  // Find first suitable individual in containing interval
  while(1) {
    // If at single individual, done
    if (step == 1) {
      if (depth == DEPTH-1) {
        abort(); // ... index-1 ...
      }

      index *= CLUSTER;
      step *= CLUSTER;
      depth += 1;

      left_z = individuals_left_z(individuals, index-1, depth);
    }

    // Decrease step size to divide into sub intervals
    step /= 2;

    // If previous sub interval doesn't contain individual, skip over it
    if (left_z[(index-step)%CLUSTER] > z)
      index -= step;
  }
}


//---------------------------------------------------------------------------------------------------------------//
// z_left: (UInt64, UInt64, UInt64,UInt64) -> ()
//
// Greatest Z value lesser than a given Z value falling into a target box.
//
// This routine works by treating the given Z value as an upper bound and reducing it to the point the desired
// Z value (the greatest one falling into the target box and lesser than the given value) is identified.
//
// The given Z value is reduced by checking the largest 2x2-aligned box immediately proceeding it for overlap
// with the target region.  If no overlap is found, the box is skipped by reducing the given Z value to the
// start of the box and repeating the procedure.  If the given Z value reaches zero, then there is no lesser Z
// value.
//
// If overlap is found, the overlapping box is broken down into 2x2 boxes.  Each 2x2 box without overlap is
// skipped by reducing the given Z value start of it.  Once an overlapping box is found, the procedure is
// repeated.
//
// The functional pseudo-code follows.
//
// z_left_out:      (UInt64, UInt64, UInt64,UInt64) -> (UInt64)
// z_left_out_size: (UInt64, UInt64, UInt64,UInt64) -> (UInt64)
// z_left_out_step: (UInt64, UInt64, UInt64,UInt64) -> (UInt64)
//
// z_left_in_size: (UInt64, UInt64, UInt64,UInt64) -> (UInt64)
// z_left_in_step: (UInt64, UInt64, UInt64,UInt64) -> (UInt64)
//
// z_overlap: (UInt64,UInt64, UInt64,UInt64) -> (Bool)
//
// z_left_out(z, step, box_ul_z,box_lr_z)
//   // At first z point, done
//   z == 0:
//     fail("no more valid z points")
//   // Not at first z point, determine step size
//   otherwise:
//     z_left_out_size(z, step, box_ul_z,box_lr_z)
//
// z_left_out_size(z, step, box_ul_z,box_lr_z)
//   // On step size boundary, increase step size
//   z % (step*4) == 0:
//     z_left_out_size(z, step*4, box_ul_z,box_lr_z)
//   // Not on step size boundary, check box to next boundary
//   otherwise:
//     z_left_out_step(z, step, box_ul_z,box_lr_z)
//
// z_left_out_step(z, step, box_ul_z,box_lr_z)
//   let ul_z = z-step
//       lr_z = z-1
//   // Previous box overlap with target box, locate z point in it
//   z_overlap(ul_z,lr_z, box_ul_z,box_lr_z):
//     z_left_in_size(z, step, box_ul_z,box_lr_z)
//   // Previous box doesn't overlap with target box, skip over it
//   otherwise:
//     z_left_out(ul_z, step, box_ul_z,box_lr_z)
//
// z_left_in_size(z, step, box_ul_z,box_lr_z)
//   // Single point, done
//   step == 1:
//     z-1
//   // Not single point, divide into sub boxes
//   otherwise:
//     z_left_in_step(z, step/4, box_ul_z,box_lr_z)
//
// z_left_in_step(z, step, box_ul_z,box_lr_z)
//   let ul_z = z-step
//       lr_z = z-1
//   // Previous sub box overlaps with target box, locate z point in it
//   z_overlap(ul_z,lr_z, box_ul_z,box_lr_z):
//     z_left_in_size(z, step, box_ul_z,box_lr_z)
//   otherwise:
//   // Previous sub box doesn't overlap with target box, skip over it
//     z_left_in_step(ul_z, step, box_ul_z,box_lr_z)
//
//  z_overlap(ul_z0,lr_z0, ul_z1,lr_z1)
//    ( (ul_z0&X_MASK) <= (lr_z1&X_MASK) && (lr_z0&X_MASK) >= (ul_z1&X_MASK) &&
//      (ul_z0&Y_MASK) <= (lr_z1&Y_MASK) && (lr_z0YX_MASK) >= (ul_z1YX_MASK) )
//
// Flattening the recursive calls with loops this becomes the following code.
//
void z_left(UInt64 z, UInt64 step, UInt64 box_ul_z,UInt64 box_lr_z) {
  // Scan backwards with progressively larger boxes for overlap with target box
  while (1) {
    // At first z point, done
    if (z == 0) {
      abort(); // ...
    }

    // While on step size boundary, increase step size
    while ((z & (step*4-1)) == 0)
      step *= 4;

    // If previous box overlaps with target box, locate z point in it
    UInt64 ul_z = z-step;
    UInt64 lr_z = z-1;

    if ( (ul_z&X_MASK) <= (box_lr_z&X_MASK) && (lr_z&X_MASK) >= (box_ul_z&X_MASK) &&
         (ul_z&Y_MASK) <= (box_lr_z&Y_MASK) && (lr_z&Y_MASK) >= (box_ul_z&Y_MASK) )
      break;

    // Previous box doesn't overlap with target box, skip over it
    z = ul_z;
  }

  // Scan progressively smaller boxes in overlapping box for first point in target box
  while (1) {
    // At single point, done
    if (step == 1) {
      abort(); // ... z-1 ...
    }

    // Decrease step size to divide into sub boxes
    step /= 4;

    // While previous sub box doesn't overlap with target box, skip over it
    while (1) {
      UInt64 ul_z = z-step;
      UInt64 lr_z = z-1;

      if ( (ul_z&X_MASK) <= (box_lr_z&X_MASK) && (lr_z&X_MASK) >= (box_ul_z&X_MASK) &&
           (ul_z&Y_MASK) <= (box_lr_z&Y_MASK) && (lr_z&Y_MASK) >= (box_ul_z&Y_MASK) )
        break;

      z = ul_z;
    }
  }
}


// sort: (Individuals) -> ()
//
// sort_both:  (Individuals, Int64,Int64, Int64) -> ()
//
// In-place pivot sort.
//
// This routine works by recursively choosing a pivot and breaking elements into a left subset <= pivot and a
// right subset > pivot.  The base case is a subset in which no elements are not equal (maximum and minimum are
// equal).  This implies the subset it fully sorted and will occur for sure once singelton subsets are reached.
// The first pivot is the middle element and later pivots are the average of the minimum and maximum.
//
// The functional pseudo-code follows.
//
// sort: (Individuals) -> ()
//
// sort_both:  (Individuals, UInt64,UInt64, UInt64) -> ()
// sort_left:  (Individuals, UInt64, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64) -> ()
// sort_right: (Individuals, UInt64, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64, UInt64,UInt64) -> ()
//
// sort(individuals)
//   // If there are individuals, use the centre element as initial pivot
//   size(individuals) > 0:
//     sort_both(individuals, 0,size(individuals)-1, individuals[size(individuals)/2])
//   // If there are no individuals, done
//   otherwise:
//     0
//
// sort_both(individuals, left_start,right_start, pivot_z)
//   // Split into left-hand side <= pivot and right-hand side > pivot
//   let (left_end,  left_min, left_max,
//        right_end, right_min,right_max) = sort_left(individuals, pivot_z,
//                                                    left_start, left_start,  UINT64_MAX,0,
//                                                    right_start,right_start, UINT64_MAX,0)
//   // Repeat on each of left-hand and right-hand if they have at least two different elements in them
//       () = left_min  < left_max:
//              sort_both(individuals, left_start,left_end,    (left_min_z +left_max_z) /2)
//       () = right_min < right_max:
//              sort_both(individuals, right_end, right_start, (right_min_z+right_max_z)/2)
//   // Done
//   ()
//
// sort_left(individuals, pivot_z,
//           left_start, left_current,  left_min_z, left_max_z,
//           right_start,right_current, right_min_z,right_max_z)
//   // Left side hasn't reached right side, process element
//   left_current <= right_current:
//     // Element <= pivot, move it into left-hand side and continue left scan with next element
//     individuals[left_current] <= pivot_z:
//       sort_left( individuals, pivot_z,
//                  left_start, left_current+1,
//                  min(individuals[left],left_min_z), max(individuals[left],left_max_z),
//                  right_start,right_current, right_min_z,right_max_z)
//     // Element > pivot, switch to right scan to find element to swap it with
//     individuals[left_current] >  pivot_z:
//       sort_right(individuals, pivot_z,
//                  left_start, left_current,  left_min_z, left_max_z,
//                  right_start,right_current, right_min_z,right_max_z)
//   // Left side reached right side, return sides
//   otherwise:
//     (right_current, left_min, left_max,
//      left_current,  right_min,right_max)
//
// sort_right(individuals, pivot_z,
//            left_start, left_current,  left_min_z, left_max_z,
//            right_start,right_current, right_min_z,right_max_z)
//   // Right side hasn't reached left side, process element
//   left_current <= right_current:
//     // Element > pivot, move it into right-hand side and continue right scan with next element
//     individuals[right_current] > pivot_z:
//       sort_right(individuals, pivot_z,
//                  left_start, left_current,    left_min_z,left_max_z,
//                  right_start,right_current-1,
//                  min(individuals[right],right_min_z),max(individuals[right],right_max_z))
//     // Element <= pivot, swap it with > element found in left scan and switch back to left scan
//     individuals[left_current] <= pivot_z:
//       swap(individuals, left_current, right_current)
//       sort_left( individuals, pivot_z,
//                  left_start, left_current,    left_min_z, left_max_z,
//                  right_start,right_current,   right_min_z,right_max_z)
//   // Reach side reached left side, return sides
//   otherwise:
//     (right_current, left_min, left_max,
//      left_current,  right_min,right_max)
//
// Flattening the recursive calls with do loops this becomes the following code.
//
void sort(Individuals* individuals) {
  // If there are individuals, use the centre element for initial pivot
  if (individuals->number)
    return sort_both(individuals, 0,individuals->number-1, individuals_z(individuals, individuals->number/2));
  // If there are no individuals, done
  else
    return;
}

static void sort_both(Individuals* individuals, UInt64 left_start,UInt64 right_start, UInt64 pivot_z) {
  UInt64 left_current = left_start;
  UInt64 left_min_z = UINT64_MAX;
  UInt64 left_max_z = 0;

  UInt64 right_current = right_start;
  UInt64 right_min_z = UINT64_MAX;
  UInt64 right_max_z = 0;

  Level1* left_level1 = individuals_level1(individuals, left_current);
  Level1* right_level1 = individuals_level1(individuals, right_current);

  // Sweep to centre splitting into left <= pivot and right > pivot by swaping left > pivot and right <= pivot
  while (1) {
    UInt64 left_z;
    UInt64 right_z;

    // Find element from left > pivot
    while (left_current <= right_current) {
      left_z = left_level1->z[left_current%CLUSTER];

      // No element found, move to next and repeat
      if (left_z <= pivot_z) {
        // Update left min and max for new left <= pivot element
        left_min_z = left_min_z <= left_z ? left_min_z : left_z;
        left_max_z = left_max_z >= left_z ? left_max_z : left_z;

        // Advance to next possible left element
        left_current += 1;
        if (left_current%CLUSTER == 0)
          left_level1 = individuals_level1(individuals, left_current);
      }
      // Element found, break
      else
        break;
    }

    // Find element from right <= pivot
    while (left_current <= right_current) {
      right_z = right_level1->z[right_current%CLUSTER];

      // No element found, move to next and repeat
      if (right_z > pivot_z) {
        // Update right min and max for new right > pivot element
        right_min_z = right_min_z <= right_z ? right_min_z : right_z;
        right_max_z = right_max_z >= right_z ? right_max_z : right_z;

        // Advance to next possible right element
        right_current -= 1;
        if (right_current%CLUSTER == CLUSTER-1)
          right_level1 = individuals_level1(individuals, right_current);
      }
      // Element found, break
      else
        break;
    }

    // Not at centre, left > pivot and right <= pivot elements found, swap them
    if (left_current <= right_current) {
      left_level1->z[left_current%CLUSTER] = right_z;
      right_level1->z[right_current%CLUSTER] = left_z;

      Individual left_individual = left_level1->individual[left_current%CLUSTER];
      Individual right_individual = right_level1->individual[right_current%CLUSTER];

      left_level1->individual[left_current%CLUSTER] = right_individual;
      right_level1->individual[right_current%CLUSTER] = left_individual;
    }
    // At centre, seperated into left <= pivot and right > pivot subsets, break
    else
      break;
  }

  // Sort new left <= pivot and right > pivot subsets
  if (left_min_z  < left_max_z)
    sort_both(individuals, left_start,  right_current,
              left_min_z /2+left_max_z /2 + (left_min_z %2+left_max_z %2)/2);
  if (right_min_z < right_max_z)
    sort_both(individuals, left_current,right_start,
              right_min_z/2+right_max_z/2 + (right_min_z%2+right_max_z%2)/2);
}


//---------------------------------------------------------------------------------------------------------------//
// die_errno:  (int, char*, ...)     -> ()
// vdie_errno: (int, char*, va_args) -> ()
//
// die_errno_explicit:  (int, int, char*, ...)     -> ()
// vdie_errno_explicit: (int, int, char*, va_args) -> ()
//
// If given printf style formatted message, print it to stderr followed by a colon, then print a description of
// the current value of errno to stderr followed by a newline, and finally exit with given error value.
//
// The second versions explicitly take an errno value, the first versions use the global errno value.
//
void die_errno(int value, char* format, ...) {
  // Setup va_list for args and invoke that version
  va_list args;

  va_start(args, format);
  vdie_errno(value, format, args);
  va_end(args);
}

void vdie_errno(int value, char* format, va_list args) {
  // Invoke explicit version with global errno value
  vdie_errno_explicit(errno, value, format, args);
}


void die_errno_explicit(int errno_original, int value, char* format, ...) {
  // Setup va_list for args and invoke that version
  va_list args;

  va_start(args, format);
  vdie_errno_explicit(errno_original, value, format, args);
  va_end(args);
}

void vdie_errno_explicit(int errno_original, int value, char* format, va_list args) {
  // Output any given printf style formatted message to stderr
  if (format) {
    vfprintf(stderr, format, args);
    fprintf(stderr, ": ");
  }

  // Output description of given errno to stderr
  char buffer[ERRNO_BUFFER];

  if (strerror_r(errno_original, buffer, sizeof(buffer)/sizeof(buffer[0])))
    fprintf(stderr,"%s\n",buffer);
  else
    fprintf(stderr,"error %d (could not be described due to error %d)\n", errno_original, errno);

  // Quit with given value
  exit(value);
}


// die:  (char*, ...)     -> ()
// vdie: (char*, va_args) -> ()
//
// If given printf style formatted message, print it to stderr followed by a newline, and then exit with given
// error value.
//
void die(int value, char* format, ...) {
  // Setup va_list for args and invoke that version
  va_list args;

  va_start(args, format);
  vdie(value, format, args);
  va_end(args);
}


void vdie(int value, char* format, va_list args) {
  // Output any given printf style formatted message to stderr
  if (format) {
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
  }

  // Quit with given value
  exit(value);
}


//---------------------------------------------------------------------------------------------------------------//
// save: (Forest*, char*) -> ()
// save_forest: (FILE*, Forest*, char*) -> ()
//
// Create/truncate the given file name and serialize the given forest to it.
//
void save(Forest* forest, char* name) {
  // Wrap file handle based routine with opening and close details
  FILE* file;

  if ( !(file = fopen(name,"w")) )
    die_errno(1, "unable to create/truncate \"%s\" for writing", name);
  save_forest(file, forest, name);
  if (fclose(file))
    die_errno(1, "unable to close \"%s\"", name);
}

static void save_forest(FILE* file, Forest* forest, char* name) {
  // Save forest level data
  if ( fprintf(file, "FOREST %d %d %g %g\n",
               forest->periodic_x, forest->periodic_y,
               forest->size_x, forest->size_y) < 0 )
    die_errno(1, "an error occured while writing to \"%s\"", name);

  // Save varieties
  Varieties* varieties = &forest->varieties;

  for (UInt varieties_iterator=0; varieties_iterator<varieties->number; varieties_iterator+=1) {
    // Save variety
    Variety *variety = &varieties->varieties[varieties_iterator];

    if ( fprintf(file, "VARIETY %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                 variety->height_mature, variety->height_maximum,
                 variety->growth_rate, variety->growth_competition_lower, variety->growth_competition_higher,
                 variety->mortality_initial, variety->mortality_decay, variety->mortality_intrinsic,
                 variety->fecundity_maximum,
                 variety->masting_time, variety->masting_phase,
                 variety->dispersal_probability_short, variety->dispersal_mode_short, variety->dispersal_mode_long)
         < 0 )
      die_errno(1, "an error occured while writing to \"%s\"", name);

    // Save individuals
    Individuals* individuals = &variety->individuals;
    Level1* level1;

    for (UInt64 individuals_iterator=0; individuals_iterator<individuals->number; individuals_iterator+=1) {
      // Switch levels if required
      if (individuals_iterator%CLUSTER == 0)
        level1 = individuals_level1(individuals, individuals_iterator);

      // Save individual
      if ( fprintf(file, "INDIVIDUAL %g %g %g\n",
                   level1->individual[individuals_iterator%CLUSTER].x,
                   level1->individual[individuals_iterator%CLUSTER].y,
                   level1->individual[individuals_iterator%CLUSTER].height) < 0 )
        die_errno(1, "an error occured while writing to \"%s\"", name);
    }
  }
}


//---------------------------------------------------------------------------------------------------------------//
// load: (char*) -> (Forest*)
// load_forest: (FILE*, char*) -> (Forest*)
//
// Read the given file and restore a serialized forest from it.
//
Forest* load(char* name) {
  Forest* forest;
  FILE* file;

  if ( !(file = fopen(name,"r")) )
    die_errno(1, "unable to open \"%s\" for reading", name);
  forest = load_forest(file, name);
  if (fclose(file))
    die_errno(1, "unable to close \"%s\"", name);

  return forest;
}


static Forest* load_forest(FILE* file, char* name) {
  UInt64 line = 1;
  int records;

  // Load forest level data
  Forest forest;
  UInt periodic_x, periodic_y;

  records = fscanf(file, "FOREST %u %u %g %g\n",
                   &periodic_x, &periodic_y,
                   &forest.size_x, &forest.size_y);
  forest.scale = 0x1p32/fmax(forest.size_x,forest.size_y);
  forest.periodic_x = periodic_x;
  forest.periodic_y = periodic_y;
  if ( ferror(file) )
    die_errno(1, "an error occured while reading from \"%s\"", name);
  if ( records < 4 )
    die(1, "problem parsing \"%s\":%"PRIu64": expecting \"FOREST\" "
        "PERIODIC_X PERIODIC_Y "
        "SIZE_X SIZE_Y", name, line);
  else
    line += 1;

  // Load varieties
  Varieties varieties = { .number = 0, .varieties = 0 };

  while (1) {
    // Load variety
    Variety variety;

    records = fscanf(file, "VARIETY %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                     &variety.height_mature, &variety.height_maximum,
                     &variety.growth_rate, &variety.growth_competition_lower, &variety.growth_competition_higher,
                     &variety.mortality_initial, &variety.mortality_decay, &variety.mortality_intrinsic,
                     &variety.fecundity_maximum,
                     &variety.masting_time, &variety.masting_phase,
                     &variety.dispersal_probability_short,
                     &variety.dispersal_mode_short, &variety.dispersal_mode_long);
    if ( ferror(file) )
      die_errno(1, "an error occured while reading from \"%s\"", name);
    else if (records == 0 || records == EOF)
      break;
    else if (records < 14)
      die(1, "problem parsing \"%s\":%"PRIu64": expecting \"VARIETY\" "
          "HEIGHT_MATURE HEIGHT_MAXIMUM "
          "GROWTH_RATE VARIETY_GROWTH_COMPETITION_LOWER VARIETY_GROWTH_COMPETITION_HIGHER "
          "MORTALITY_INITIAL MORTALITY_DECAY MORTALITY_INTRINSIC "
          "FECUNDITY_MAXIMUM "
          "MASTING_TIME MASTING_PHASE "
          "DISPERSAL_PROBABILITY_SHORT DISPERSAL_MODE_SHORT DISPERSAL_MODE_LONG", name, line);
    else
      line += 1;

    // Load individuals
    Individuals individuals = { .number = 0 };

    while (1) {
      // Load individual
      Individual individual;
      Level1* level1;

      records = fscanf(file, "INDIVIDUAL %g %g %g\n",
                       &individual.x, &individual.y, &individual.height);
      if ( ferror(file) )
        die_errno(1, "an error occured while reading from \"%s\"", name);
      else if (records == 0 || records == EOF)
        break;
      else if (records < 3)
        die(1, "problem parsing \"%s\":%"PRIu64": expecting \"INDIVIDUAL\" "
            "X Y "
            "HEIGHT", name, line);
      else if (individual.x < 0 || individual.x >= forest.size_x)
        die(1, "problem parsing \"%s\":%"PRIu64": expecting \"INDIVIDUAL\" "
            "the constraint 0 <= X=%g < SIZE_X=%g does not hold", name, line, individual.x, forest.size_x);
      else if (individual.y < 0 || individual.y >= forest.size_y)
        die(1, "problem parsing \"%s\":%"PRIu64": expecting \"INDIVIDUAL\" "
            "the constraint 0 <= Y=%g < SIZE_Y=%g does not hold", name, line, individual.y, forest.size_y);
      else if (individual.height < 0)
        die(1, "problem parsing \"%s\":%"PRIu64": expecting \"INDIVIDUAL\" "
            "the constraint 0 <= HEIGHT=%g", name, line, individual.height);
      else
        line += 1;

      // Add individual
      if ( !(individuals.number%CLUSTER) )
        level1 = individuals_append(&individuals);

      level1->z[individuals.number%CLUSTER] = z_xy((UInt32)(forest.scale*individual.x),
                                                   (UInt32)(forest.scale*individual.y));
      level1->individual[individuals.number%CLUSTER] = individual;
      individuals.number += 1;
    }

    // Add individuals
    variety.individuals = individuals;

    // Add variety (memory is allocated in chunks of 2^n-1: 0,1,3,7,15,...)
    if ( !(varieties.number & (varieties.number+1)) )  // Need to a reallocate if (varieties.number+1) == 2^n
      if ( !(varieties.varieties = realloc(varieties.varieties, sizeof(Variety)*((varieties.number+1)*2-1))) )
        die_errno(1, "unable to grow varieties allocation to %tu bytes",
                  sizeof(Variety)*((varieties.number+1)*2-1));
    varieties.varieties[varieties.number] = variety;
    varieties.number += 1;
  }

  // Add varieties (memory is shrunk to an exact fit)
  if ( !(varieties.varieties = realloc(varieties.varieties, sizeof(Variety)*varieties.number)) )
    die_errno(1, "unable to shrink varieties allocation to %tu bytes", sizeof(Variety)*varieties.number);
  forest.varieties = varieties;

  // Allocate forest for return
  Forest* forest_allocated;

  if ( !(forest_allocated = malloc(sizeof(Forest))) )
    die_errno(1, "unable to allocate %tu bytes for new Forest", sizeof(Forest));
  *forest_allocated = forest;

  return forest_allocated;
}


//---------------------------------------------------------------------------------------------------------------//
// Program entrance point

int main() {

  return 0;
}
