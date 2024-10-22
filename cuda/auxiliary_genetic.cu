/*

OmegaDocker, a GPU-accelerated docking program developed for larger ligands, such as oligopeptides, oligonucleotides, and other oligomers.
Copyright (C) 2024 Gao Quan-Ze. All rights reserved.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.

*/


// -------------------------------------------------------
//
// -------------------------------------------------------
inline __device__ uint32_t gpu_rand(uint32_t* prng_states, int i)
// The GPU device function generates a random int
// with a linear congruential generator.
// Each thread (supposing num_of_runs*pop_size blocks and NUM_OF_THREADS_PER_BLOCK threads per block)
// has its own state which is stored in the global memory area pointed by
// prng_states (thread with ID tx in block with ID bx stores its state in prng_states[bx*NUM_OF_THREADS_PER_BLOCK+$
// The random number generator uses the gcc linear congruential generator constants.
{
  uint state;
  // Current state of the threads own PRNG
  state = prng_states[i];
  // Calculating next state
  state = (RAND_A*state+RAND_C);
  // Saving next state to memory
  prng_states[i] = state;
  return state;
}

// -------------------------------------------------------
//
// -------------------------------------------------------
inline __device__ float gpu_randf(uint32_t* prng_states, int i)
// The GPU device function generates a
// random float greater than (or equal to) 0 and less than 1.
// It uses gpu_rand() function.
{
  float state;
  // State will be between 0 and 1
  state =  ((float)gpu_rand(prng_states, i) / (float)MAX_UINT)*0.999999f;
  return state;
}

// -------------------------------------------------------
//
// -------------------------------------------------------
inline __device__ void map_angle(float& angle)
// The GPU device function maps
// the input parameter to the interval 0...360
// (supposing that it is an angle).
{
  while (angle >= 360.0f) {
    angle -= 360.0f;
  }

  while (angle < 0.0f) {
    angle += 360.0f;
  }
}

