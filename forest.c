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
// Program entrance point

int main() {

  return 0;
}
