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

#include <stdint.h>
#include <inttypes.h>

#include "random.h"


//---------------------------------------------------------------------------------------------------------------//
// Stats (sample mean and variance) on draws for the various distributions
MersenneTwister uniform_UInt32(MersenneTwister mersennetwister, UInt32 const n) {
  Float32 const mean = n*(n-1.0)/2.0/n;
  Float32 const variance = (n-1.0)*n*(2.0*n-1.0)/6.0/n-mean*mean;

  UInt32 max = 0;
  UInt32 min = UINT32_MAX;
  Float64 mean_sum = 0;
  Float64 variance_sum = 0;

  for ( UInt iter=0; iter<1000000; iter+=1 ) {
    UInt32 value;
    unpack_MersenneTwister_UInt32( &mersennetwister, &value,
                                   Random_uniform_UInt32(mersennetwister, n) );
    max = max >= value ? max : value;
    min = min <= value ? min : value;
    mean_sum = mean_sum + value;
    variance_sum = variance_sum + value*value;
  }

  Float32 mean_ = mean_sum/1000000.0;
  Float32 variance_ = (variance_sum-mean_*mean_sum)/(1000000.0-1.0);

  printf("1000000 uniform_UInt32(%"PRIu32"):\n"
         "  min = %"PRIu32"\n"
         "  max = %"PRIu32"\n"
         "  mean = %g (%g)\n"
         "  variance = %g (%g)\n",
         n, min, max, mean_, mean, variance_, variance);

  return mersennetwister;
}


MersenneTwister binomial_UInt32(MersenneTwister mersennetwister, UInt32 const n, Float32 const p) {
  Float32 const mean = n*p;
  Float32 const variance = n*p*(1.0-p);

  UInt32 max = 0;
  UInt32 min = UINT32_MAX;
  Float64 mean_sum = 0;
  Float64 variance_sum = 0;

  for ( UInt iter=0; iter<1000000; iter+=1 ) {
    UInt32 value;
    unpack_MersenneTwister_UInt32( &mersennetwister, &value,
                                   Random_binomial_UInt32(mersennetwister, n, p) );
    max = max >= value ? max : value;
    min = min <= value ? min : value;
    mean_sum = mean_sum + value;
    variance_sum = variance_sum + value*value;
  }

  Float32 mean_ = mean_sum/1000000.0;
  Float32 variance_ = (variance_sum-mean_*mean_sum)/(1000000.0-1.0);

  printf("1000000 binomial_UInt32(%"PRIu32",%g):\n"
         "  min = %"PRIu32"\n"
         "  max = %"PRIu32"\n"
         "  mean = %g (%g)\n"
         "  variance = %g (%g)\n",
         n, p, min, max, mean_, mean, variance_, variance);

  return mersennetwister;
}


MersenneTwister uniform_Float32(MersenneTwister mersennetwister) {
  Float32 const mean = 1.0/2.0;
  Float32 const variance = 1.0/12.0;

  Float32 max = -INFINITY;
  Float32 min = INFINITY;
  Float64 mean_sum = 0;
  Float64 variance_sum = 0;

  for ( UInt iter=0; iter<1000000; iter+=1 ) {
    Float32 value;
    unpack_MersenneTwister_Float32( &mersennetwister, &value,
                                    Random_uniform_Float32(mersennetwister) );
    max = max >= value ? max : value;
    min = min <= value ? min : value;
    mean_sum = mean_sum + value;
    variance_sum = variance_sum + value*value;
  }

  Float32 mean_ = mean_sum/1000000.0;
  Float32 variance_ = (variance_sum-mean_*mean_sum)/(1000000.0-1.0);

  printf("1000000 uniform_Float32():\n"
         "  min = %g\n"
         "  max = %g\n"
         "  mean = %g (%g)\n"
         "  variance = %g (%g)\n",
         min, max, mean_, mean, variance_, variance);

  return mersennetwister;
}


MersenneTwister normal2_Float32(MersenneTwister mersennetwister) {
  Float32 const mean = 0.0;
  Float32 const variance = 1.0;

  Float32 max = -INFINITY;
  Float32 min = INFINITY;
  Float64 mean_sum = 0;
  Float64 variance_sum = 0;

  for ( UInt iter=0; iter<1000000/2; iter+=1 ) {
    Float32 value0, value1;
    unpack_MersenneTwister_Float32_Float32( &mersennetwister, &value0, &value1,
                                            Random_normal2_Float32(mersennetwister) );
    max = max >= value0 ? max : value0;
    min = min <= value0 ? min : value0;
    mean_sum = mean_sum + value0;
    variance_sum = variance_sum + value0*value0;

    max = max >= value1 ? max : value1;
    min = min <= value1 ? min : value1;
    mean_sum = mean_sum + value1;
    variance_sum = variance_sum + value1*value1;
  }

  Float32 mean_ = mean_sum/1000000.0;
  Float32 variance_ = (variance_sum-mean_*mean_sum)/(1000000.0-1.0);

  printf("1000000 normal2_Float32():\n"
         "  min = %g\n"
         "  max = %g\n"
         "  mean = %g (%g)\n"
         "  variance = %g (%g)\n",
         min, max, mean_, mean, variance_, variance);

  return mersennetwister;
}


//---------------------------------------------------------------------------------------------------------------//
//
int main() {
  MersenneTwister mersennetwister = MersenneTwister_begin();

  mersennetwister = uniform_UInt32(mersennetwister, 10);
  mersennetwister = binomial_UInt32(mersennetwister, 10, 0.25);  // nr < 10.0 -- BINV region
  mersennetwister = binomial_UInt32(mersennetwister, 40, 0.75);  // nr >= 10.0 -- BTPE region
  mersennetwister = uniform_Float32(mersennetwister);
  mersennetwister = normal2_Float32(mersennetwister);

  MersenneTwister_end(mersennetwister);

  return 0;
}
