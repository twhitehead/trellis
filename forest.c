#include <stdint.h>

// Split of 64b address space into 16b components (can't be changed without changing grove walking code)

#define HEIGHT  4
#define CLUSTER 65536


// Single upper case to distinguish types

typedef unsigned int UInt;
typedef uint32_t UInt32;
typedef uint64_t UInt64;

typedef int Int;
typedef int32_t Int32;
typedef int64_t Int64;

typedef struct _Forest Forest;
typedef struct _Groves Groves;
typedef struct _Grove Grove;


// Basic forest structure

struct _Forest {
  UInt64 size;
  Groves* groves;
};

struct _Groves {
  UInt64 z_min[CLUSTER];
  UInt64 z_max[CLUSTER];
  union {
    Groves* groves[CLUSTER];
    Grove* grove[CLUSTER];
  };
};

struct _Grove {
  UInt64 z[CLUSTER];
};


// In-place pivot sort.
//
// This routine works by recursively choosing a pivot and breaking elements into a left subset <= pivot and a
// right subset > pivot.  The base case is a subset in which no elements are not equal (maximum and minimum are
// equal).  This implies the subset it fully sorted and will occur for sure once singelton subsets are reached.
// The first pivot is the middle element and later pivots are the average of the minimum and maximum.
//
// The functional pseudo-code follows.
//
// sort(trees, left_start,right_start, pivot_z)
//   // Split into left-hand side <= pivot and right-hand side > pivot
//   let (left_end,  left_min, left_max,
//        right_end, right_min,right_max) = sort_left(trees, pivot_z,
//                                                    left_start, left_start,  UINT64_MAX,0,
//                                                    right_start,right_start, UINT64_MAX,0)
//   // Repeat on each of left-hand and right-hand if they have at least two different elements in them
//       () = left_min  < left_max:
//              sort(trees, left_start,left_end,    (left_min_z +left_max_z) /2)
//       () = right_min < right_max:
//              sort(trees, right_end, right_start, (right_min_z+right_max_z)/2)
//   // Done
//   ()
//
// sort_left(trees, pivot_z,
//           left_start, left_current,  left_min_z, left_max_z,
//           right_start,right_current, right_min_z,right_max_z)
//   // Left side hasn't reached right side, process element
//   left_current <= right_current:
//     // Element <= pivot, move it into left-hand side and continue left scan with next element
//     trees[left_current] <= pivot_z:
//       sort_left( trees, pivot_z,
//                  left_start, left_current+1, min(trees[left], left_min_z), max(trees[left], left_max_z),
//                  right_start,right_current,                   right_min_z,                  right_max_z)
//     // Element > pivot, switch to right scan to find element to swap it with
//     trees[left_current] >  pivot_z:
//       sort_right(trees, pivot_z,
//                  left_start, left_current,                    left_min_z,                   left_max_z,
//                  right_start,right_current,                   right_min_z,                  right_max_z)
//   // Left side reached right side, return sides
//   otherwise:
//     (right_current, left_min, left_max,
//      left_current,  right_min,right_max)
//
// sort_right(trees, pivot_z,
//            left_start, left_current,  left_min_z, left_max_z,
//            right_start,right_current, right_min_z,right_max_z)
//   // Right side hasn't reached left side, process element
//   left_current <= right_current:
//     // Element > pivot, move it into right-hand side and continue right scan with next element
//     trees[right_current] > pivot_z:
//       sort_right(trees, pivot_z,
//                  left_start, left_current,                    left_min_z,                   left_max_z),
//                  right_start,right_current-1,min(trees[right],right_min_z),max(trees[right],right_max_z))
//     // Element <= pivot, swap it with > element found in left scan and switch back to left scan
//     trees[left_current] <= pivot_z:
//       swap(trees, left_current, right_current)
//       sort_left( trees, pivot_z,
//                  left_start, left_current,                    left_min_z,                   left_max_z,
//                  right_start,right_current,                   right_min_z,                  right_max_z)
//   // Reach side reached left side, return sides
//   otherwise:
//     (right_current, left_min, left_max,
//      left_current,  right_min,right_max)
//
// Flattening the recursive calls with do loops this becomes the following code.
//
void sort(Forest* forest) {
  // Use centre element for initial pivot
  sort_both(forest, 0,forest->size-1,
            forest->groves
            ->groves[forest->size/2/CLUSTER/CLUSTER/CLUSTER%CLUSTER]
            ->groves[forest->size/2/CLUSTER/CLUSTER%CLUSTER]
            ->grove[forest->size/2/CLUSTER%CLUSTER]
            ->z[forest->size/2%CLUSTER]);
}

void sort_both(Forest* forest, UInt64 left_start,UInt64 right_start, UInt64 pivot_z) {
  UInt64 left_current = left_start;
  UInt64 left_min_z = UINT64_MAX;
  UInt64 left_max_z = 0;

  UInt64 right_current = right_start;
  UInt64 right_min_z = UINT64_MAX;
  UInt64 right_max_z = 0;

  left_grove = forest->groves
    ->groves[left_current/CLUSTER/CLUSTER/CLUSTER%CLUSTER]
    ->groves[left_current/CLUSTER/CLUSTER%CLUSTER]
    ->groves[left_current/CLUSTER%CLUSTER];
  right_grove = forest->groves
    ->groves[right_current/CLUSTER/CLUSTER/CLUSTER%CLUSTER]
    ->groves[right_current/CLUSTER/CLUSTER%CLUSTER]
    ->groves[right_current/CLUSTER%CLUSTER];

  // Sweep to centre splitting into left <= pivot and right > pivot by swaping left > pivot and right <= pivot
  while (1) {
    UInt64 left_z;
    UInt64 right_z;

    // Find element from left > pivot
    while (left_current <= right_current) {
      left_z = left_grove->z[left_current%CLUSTER];

      // No element found, move to next and repeat
      if (left_z <= pivot_z) {
        // Update left min and max for new left <= pivot element
        left_min_z = min(left_min_z,left_z);
        left_max_z = max(left_z,left_max_z);

        // Advance to next possible left element
        left_current += 1;
        if (left_current%CLUSTER == 0)
          left_grove = forest->groves
            ->groves[left_current/CLUSTER/CLUSTER/CLUSTER%CLUSTER]
            ->groves[left_current/CLUSTER/CLUSTER%CLUSTER]
            ->grove[left_current/CLUSTER%CLUSTER];
      }
      // Element found, break
      else
        break;
    }

    // Find element from right <= pivot
    while (left_current <= right_current) {
      right_z = right_grove->z[right_current%CLUSTER];

      // No element found, move to next and repeat
      if (right_z > pivot_z) {
        // Update right min and max for new right > pivot element
        right_min_z = min(right_min_z,right_z);
        right_max_z = max(right_z,right_max_z);

        // Advance to next possible right element
        right_current -= 1;
        if (right_current%CLUSTER == CLUSTER-1)
          right_grove = forest->groves
            ->groves[right_current/CLUSTER/CLUSTER/CLUSTER%CLUSTER]
            ->groves[right_current/CLUSTER/CLUSTER%CLUSTER]
            ->grove[right_current/CLUSTER%CLUSTER];
      }
      // Element found, break
      else
        break;
    }

    // Not at centre, left > pivot and right <= pivot elements found, swap them
    if (left_current <= right_current) {
      left_grove->z  = right_z;
      right_grove->z = left_z;
    }
    // At centre, seperated into left <= pivot and right > pivot subsets, break
    else
      break;
  }

  // Sort new left <= pivot and right > pivot subsets
  if (left_min_z  < left_max_z)
    sort_both(forest, left_start,  right_current, left_min_z /2+left_max_z /2 + (left_min_z %2+left_max_z %2)/2);
  if (right_min_z < right_max_z)
    sort_both(forest, left_current,right_start,   right_min_z/2+right_max_z/2 + (right_min_z%2+right_max_z%2)/2);
}


// Program entrance point

int main() {

  return 0;
}
