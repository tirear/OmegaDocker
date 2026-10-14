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


#define ADADELTA_AUTOSTOP

// Gradient-based adadelta minimizer
// https://arxiv.org/pdf/1212.5701.pdf
// Alternative to Solis-Wets / Steepest-Descent / FIRE

// "rho": controls degree of memory of previous gradients
//        ranges between [0, 1[
//        "rho" = 0.9 most popular value
// "epsilon":  to better condition the square root

// Adadelta parameters (TODO: to be moved to header file?)
//#define RHO             0.9f
//#define EPSILON         1e-6
#define RHO             0.8f
#define EPSILON         1e-2f

// Enabling "DEBUG_ENERGY_ADADELTA" requires
// manually enabling "DEBUG_ENERGY_KERNEL" in calcenergy.cl
//#define DEBUG_ENERGY_ADADELTA
//#define PRINT_ADADELTA_ENERGIES
//#define PRINT_ADADELTA_GENES_AND_GRADS
//#define PRINT_ADADELTA_ATOMIC_COORDS
//#define DEBUG_SQDELTA_ADADELTA

// Enable DEBUG_ADADELTA_MINIMIZER for a seeing a detailed ADADELTA evolution
// If only PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION is enabled,
// then a only a simplified ADADELTA evolution will be shown
//#define DEBUG_ADADELTA_MINIMIZER
//#define PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION

// Enable this for debugging ADADELTA from a defined initial genotype
//#define DEBUG_ADADELTA_INITIAL_2BRT


__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_gradient_minAD_kernel(
                          float* pMem_conformations_next,
                          float* pMem_energies_next,
                          // sans ==>
                          int   *pMem_energyToConforIndex,
                          int   *pMem_evals_of_new_entities,
                 // volatile float *gFloatBuff,
                          float *gFloatBuff,
                          size_t elemNum_gFloatBuffBase,
                          int    run_win_id,
                          int    run_id

//                           float    *atom_charges_const,
//                           uint32_t *atom_types_const,
//                           uint32_t *atom_types_map_const,
//                           char     *ignore_inter_const,
//
//                           uint32_t *intraE_contributors_const,

//                           int      *rotlist_const,
//
//                           float    *ref_coords_const,
//                           float    *rotbonds_moving_vectors_const,
//                           float    *rotbonds_unit_vectors_const,
//
//                           int      *rotbonds,
//                           int      *rotbonds_atoms,
//                           int      *num_rotating_atoms_per_rotbond
                          // <== sans
                         )
// The GPU global function performs gradient-based minimization on (some) entities of conformations_next.
// The number of OpenCL compute units (CU) which should be started equals to num_of_minEntities*num_of_runs.
// This way the first num_of_lsentities entity of each population will be subjected to local search
// (and each CU carries out the algorithm for one entity).
// Since the first entity is always the best one in the current population,
// it is always tested according to the ls probability, and if it not to be
// subjected to local search, the entity with ID num_of_lsentities is selected instead of the first one (with ID 0).
{
//   long long int time1 = clock64(); // sans clock
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------

  // Determining entity, and its run, energy, and genotype
//   int   run_id = blockIdx.x / cData.dockpars.num_of_lsentities; // sans
  float energy;
  float energy2; // sans check12
  int   tooShortDistNum; // sans
  int   tooFarDistNum;   // sans
  // Energy may go up, so we keep track of the best energy ever calculated.
  // Then, we return the genotype corresponding
  // to the best observed energy, i.e. "best_genotype"
  __shared__ int entity_id;
  __shared__ float best_energy;
  // __shared__ float previousEnergy; // sans check12
  __shared__ float bestDrugEnergy; // sans check12
  __shared__ float sFloatAccumulator;
  __shared__ int sIntAccumulator; // sans

//   volatile __shared__ float sFloatBuff[6 * MAX_NUM_OF_ATOMS + 5 * ACTUAL_GENOTYPE_LENGTH]; // sans
// calc_coords cartesian_gradient genotype best_genotype gradient square_delta square_gradient
// GENOTYPE_LENGTH_IN_GLOBMEM 64
// MAX_NUM_OF_ROTBONDS 57
// ACTUAL_GENOTYPE_LENGTH (MAX_NUM_OF_ROTBONDS+6)

  // Ligand-atom position and partial energies
  // sans ==>
  int actual_genotype_length = cData.dockpars.num_of_genes;
  int genotype_length_in_globmem = actual_genotype_length + 1;

  // sans check ==>
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next1a[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
//   __syncthreads();
  // <== sans check
  // sans check ==>
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next1a[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
//   __syncthreads();
  // <== sans check
  float3 *calc_coords;
  float  *genotype;
#ifdef FLOAT_GRADIENTS
    float3* cartesian_gradient;
#else
    int3*   cartesian_gradient;
#endif
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (6 * NUM_OF_THREADS_PER_BLOCK + 5 * actual_genotype_length));
#ifdef FLOAT_GRADIENTS
    cartesian_gradient = (float3*)(calc_coords + NUM_OF_THREADS_PER_BLOCK);
#else
    cartesian_gradient = (int3*)  (calc_coords + NUM_OF_THREADS_PER_BLOCK);
#endif
    genotype = (float*)(cartesian_gradient + NUM_OF_THREADS_PER_BLOCK); // so far used 3*2*NUM_OF_THREADS_PER_BLOCK
  }
  else
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (6 * cData.dockpars.num_of_atoms + 5 * actual_genotype_length));
#ifdef FLOAT_GRADIENTS
    cartesian_gradient = (float3*)(calc_coords + cData.dockpars.num_of_atoms);
#else
    cartesian_gradient = (int3*)  (calc_coords + cData.dockpars.num_of_atoms);
#endif
    genotype = (float*)(cartesian_gradient + cData.dockpars.num_of_atoms); // so far used 3*2*cData.dockpars.num_of_atoms
  }
  // <== sans


//   float3* calc_coords = (float3*)sFloatBuff; // sans
  // Gradient of the intermolecular energy per each ligand atom
  // Also used to store the accummulated gradient per each ligand atom
// #ifdef FLOAT_GRADIENTS
//   float3* cartesian_gradient = (float3*)(calc_coords + MAX_NUM_OF_ATOMS); // sans
// #else
//   int3* cartesian_gradient = (int3*)(calc_coords + MAX_NUM_OF_ATOMS); // sans
// #endif
  // Genotype pointers
//   float* genotype = (float*)(cartesian_gradient + MAX_NUM_OF_ATOMS); // so far used 3*2*MAX_NUM_OF_ATOMS // sans
//   float* best_genotype = (float*) (genotype + ACTUAL_GENOTYPE_LENGTH); // sans
  float* best_genotype = (float*) (genotype + actual_genotype_length); // sans

  // Partial results of the gradient step
//   float* gradient = best_genotype + ACTUAL_GENOTYPE_LENGTH; // sans
  float* gradient = best_genotype + actual_genotype_length; // sans

  // Squared updates E[dx^2]
//   float* square_delta = gradient + ACTUAL_GENOTYPE_LENGTH; // sans
  float* square_delta = gradient + actual_genotype_length; // sans

  // Vector for storing squared gradients E[g^2]
//   float* square_gradient = square_delta + ACTUAL_GENOTYPE_LENGTH; // so far used 5*ACTUAL_GENOTYPE_LENGTH // sans
  float* square_gradient = square_delta + actual_genotype_length; // so far used 5*actual_genotype_length // sans

  // Iteration counter for the minimizer
  uint32_t iteration_cnt = 0;

  if (threadIdx.x == 0)
  {
    // Since entity 0 is the best one due to elitism,
    // it should be subjected to random selection
    entity_id = blockIdx.x % cData.dockpars.num_of_lsentities; // sans this should always be equal to blockIdx.x
    if (entity_id == 0) {
      // If entity 0 is not selected according to LS-rate,
      // choosing another entity
//       if (100.0f*gpu_randf(cData.pMem_prng_states) > cData.dockpars.lsearch_rate) { // sans check
      if (100.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x) > cData.dockpars.lsearch_rate) { // sans check
        entity_id = cData.dockpars.num_of_lsentities; // AT - Should this be (uint)(dockpars_pop_size * gpu_randf(dockpars_prng_states))?
      }
    }

    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION) // sans
//     printf("\n");
    printf("-------> Start of ADADELTA minimization cycle\n");
    printf("%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n", "threadIdx.x: ", threadIdx.x, "blockDim.x: ", blockDim.x, "blockIdx.x: ", blockIdx.x, "gridDim.x: ", gridDim.x, "run_id: ", run_id, "entity_id: ", entity_id, "num_of_lsentities: ", cData.dockpars.num_of_lsentities, "pop_size: ", cData.dockpars.pop_size, "num_of_genes", cData.dockpars.num_of_genes);
// sans ./host/src/performdocking.cpp.Cuda:483: cData.dockpars.num_of_lsentities        = (unsigned int) (mypars->lsearch_rate/100.0*mypars->pop_size + 0.5);
// sans ./host/src/performdocking.cpp.Cuda:466: cData.dockpars.num_of_genes                 = myligand_reference.num_of_rotbonds + 6;
// sans blockNum = blocksPerGridForEachGradMinimizerEntity = cData.dockpars.num_of_lsentities * mypars->num_of_runs;
// sans num_of_lsentities = ga_pop_size
// sans num_of_runs = ga_run
//     printf("%20s %6u\n", "run_id: ", run_id);
//     printf("%20s %6u\n", "entity_id: ", entity_id);
//     printf("%20s %6u\n", "num_of_lsentities: ", cData.dockpars.num_of_lsentities);
//     printf("\n");
//     printf("%20s \n", "LGA genotype: ");
//     printf("%20s %.6f\n", "initial energy: ", energy);
    #endif // sans
  }
  __syncthreads();
//   energy = pMem_energies_next[run_id * cData.dockpars.pop_size + entity_id]; // sans
  // energy = pMem_energies_next[entity_id]; // sans
// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy1: %d %f\n", blockIdx.x, pMem_energies_next[blockIdx.x]);
// }

//   int offset = (run_id * cData.dockpars.pop_size + entity_id) * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
  int offset = entity_id * genotype_length_in_globmem; // sans
//   long long int time2 = clock64(); // sans clock
  for (int i = threadIdx.x ; i < cData.dockpars.num_of_genes; i += blockDim.x)
// sans cData.dockpars.num_of_genes = myligand_reference.num_of_rotbonds + 6
  {
    genotype[i] = pMem_conformations_next[offset + i];
  }
//   long long int time3 = clock64(); // sans clock

  // sans check ==>
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next1[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
  // <== sans check

  // -------------------------------------------------------------------
  // Calculate gradients (forces) for intermolecular energy
  // Derived from autodockdev/maps.py
  // -------------------------------------------------------------------

  #if defined (DEBUG_ENERGY_KERNEL)
  float interE;
  float intraE;
  #endif

  // Update vector, i.e., "delta".
  // It is added to the genotype to create the next genotype.
  // E.g. in steepest descent "delta" is -1.0 * stepsize * gradient

  // Asynchronous copy should be finished by here
  __syncthreads();

  // Initializing vectors
  for(uint32_t i = threadIdx.x;
               i < cData.dockpars.num_of_genes;
               i+= blockDim.x)
  {
    gradient[i]        = 0.0f;
    square_gradient[i] = 0.0f;
    square_delta[i]    = 0.0f;
    best_genotype[i] = genotype[i];
  }
//   long long int time4 = clock64(); // sans clock

  // Initializing best energy
  if (threadIdx.x == 0) {
    best_energy = INFINITY;
    // previousEnergy = 0.0; // sans check12
    bestDrugEnergy = INFINITY; // sans check12
  }

#ifdef ADADELTA_AUTOSTOP
  // __shared__ float rho; // sans
  __shared__ int rhoL; // sans
  __shared__ int rho;  // sans
  __shared__ int cons_succ;
  __shared__ int cons_fail;
  if (threadIdx.x == 0) {
    // rho = 1.0f; // sans
    rho = 0;   // sans
    rhoL = cData.dockpars.ligaRhoLimit; // sans
    cons_succ = 0;
    cons_fail = 0;
  }
#endif
  __syncthreads(); // sans check12

  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // Perform adadelta iterations

  // The termination criteria is based on
  // a maximum number of iterations, and
  // the minimum step size allowed for single-floating point numbers
  // (IEEE-754 single float has a precision of about 6 decimal digits)
//   long long int time5-1 = clock64(); // sans clock
// sans ligand ==>
  uint halfMaxIter;
  uint maxIter;
  maxIter = cData.dockpars.max_num_of_iters;
  halfMaxIter = maxIter / 2;
  if(halfMaxIter == 0)
  {
    halfMaxIter++;
    maxIter++;
  }

  do {
    // Printing number of ADADELTA iterations
    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
    if (threadIdx.x == 0) {
      #if defined (DEBUG_ADADELTA_MINIMIZER)
      printf("%s\n", "----------------------------------------------------------");
      #endif
      printf("%-15s %-3u ", "# ADADELTA iteration: ", iteration_cnt);
    }
    #endif

    // =============================================================
    // =============================================================
    // =============================================================
    // Calculating energy & gradient
    __syncthreads();

// for(int i = threadIdx.x;
//             i < cData.dockpars.num_of_genes;
//             i+= blockDim.x)
// {
//   if(isnan(genotype[i]))
//   {
//     printf("sans check13 NAN kernel_ad1 %d %f run_id:%d\n", i, genotype[i], run_id);
//   }
// }
    if(rhoL > 0 && blockIdx.x >= cData.pMem_reservedIdvNum[run_id]) // sans
    { // sans
      gpu_calc_energrad_ligand
      (
                      // Some OpenCL compilers don't allow declaring
                      // local variables within non-kernel functions.
                      // These local variables must be declared in a kernel,
                      // and then passed to non-kernel functions.
                      genotype,
                      energy,
                      energy2, // sans check12
                      tooShortDistNum, // sans
                      tooFarDistNum,   // sans
//                       run_id, // sans useless
                      calc_coords,
                      #if defined (DEBUG_ENERGY_KERNEL)
                      interE,
                      intraE,
                      #endif
                      // Gradient-related arguments
                      // Calculate gradients (forces) for intermolecular energy
                      // Derived from autodockdev/maps.py
                      cartesian_gradient,
                      gradient,
                      &sFloatAccumulator,
                      &sIntAccumulator
                      // sans ==>
//                       atom_charges_const,
//                       atom_types_const,
//                       atom_types_map_const,
//                       ignore_inter_const,
//
//                       intraE_contributors_const,

//                       rotlist_const,
//
//                       ref_coords_const,
//                       rotbonds_moving_vectors_const,
//                       rotbonds_unit_vectors_const,
//
//                       rotbonds,
//                       rotbonds_atoms,
//                       num_rotating_atoms_per_rotbond
                      // <== sans
      );
    // sans ==>
    }
    else
    {
      rhoL = -1;
      tooShortDistNum = 0;
      tooFarDistNum = 0;
    }
    // <== sans
    // sans check12 ==>
// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy2: %d %f run_id:%d\n", blockIdx.x, energy, run_id);
// }
    // if(previousEnergy != 0.0)
    // {
    //   printf("sans check12 %d %d %f %f %f\n", blockIdx.x, threadIdx.x, energy, previousEnergy, energy - previousEnergy);
    // }
    // previousEnergy = energy;
    // <== sans check12

    // =============================================================
    // =============================================================
    // =============================================================
    #if defined (DEBUG_ENERGY_ADADELTA)
    if (threadIdx.x == 0) {
      #if defined (PRINT_ADADELTA_ENERGIES)
      printf("\n");
      printf("%-10s %-10.6f \n", "intra: ",  intraE);
      printf("%-10s %-10.6f \n", "grids: ",  interE);
      printf("%-10s %-10.6f \n", "Energy: ", intraE + interE);
      #endif

      #if defined (PRINT_ADADELTA_GENES_AND_GRADS)
      for(uint i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %13s %5s %15s %15s\n", "gene_id", "gene.value", "|", "gene.grad", "(autodockdevpy units)");
        }
        printf("%13u %13.6f %5s %15.6f %15.6f\n", i, genotype[i], "|", gradient[i], (i<3)? (gradient[i]/0.375f):(gradient[i]*180.0f/PI_FLOAT));
      }
      #endif

      #if defined (PRINT_ADADELTA_ATOMIC_COORDS)
      for(uint i = 0; i < cData.dockpars.num_of_atoms; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%s\n", "Coordinates calculated by calcenergy.cl");
          printf("%12s %12s %12s %12s\n", "atom_id", "coords.x", "coords.y", "coords.z");
        }
        printf("%12u %12.6f %12.6f %12.6f\n", i, calc_coords_x[i], calc_coords_y[i], calc_coords_z[i]);
      }
      printf("\n");
      #endif
    }
    __syncthreads();
    #endif // DEBUG_ENERGY_ADADELTA

    for(int i = threadIdx.x;
            i < cData.dockpars.num_of_genes;
            i+= blockDim.x)
    {
      // if (energy < best_energy) // we need to be careful not to change best_energy until we had a chance to update the whole array
      //   best_genotype[i] = genotype[i];

      // Accumulating gradient^2 (eq.8 in the paper)
      // square_gradient corresponds to E[g^2]
      square_gradient[i] = RHO * square_gradient[i] + (1.0f - RHO) * gradient[i] * gradient[i];

      // Computing update (eq.9 in the paper)
      // float delta = -1.0f * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON)); // sans check12
      // sans check12 ==>
      // add step size
      float delta = cData.dockpars.negaLigaStepSize * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON));
      // <== sans check12

      // Accumulating update^2
      // square_delta corresponds to E[dx^2]
      square_delta[i] = RHO * square_delta[i] + (1.0f - RHO) * delta * delta;

      // Applying update
// if(isnan(genotype[i]))
// {
//   printf("sans check13 NAN kernel_ad2 %d %f\n", i, genotype[i]);
// }
// float ggg = genotype[i];
      genotype[i] += delta;
// if(isnan(genotype[i]))
// {
//   printf("sans check13 NAN kernel_ad %d ori%f %f d%f s%f g%f\n", i, ggg, genotype[i], delta, square_gradient[i], gradient[i]);
// }
    }
    __syncthreads();

    #if defined (DEBUG_SQDELTA_ADADELTA)
    if (/*(get_group_id(0) == 0) &&*/ (threadIdx.x == 0)) {
      for(int i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %20s %15s %15s %15s\n", "gene", "sq_grad", "delta", "sq_delta", "new.genotype");
        }
        printf("%13u %20.6f %15.6f %15.6f %15.6f\n", i, square_gradient[i], delta[i], square_delta[i], genotype[i]);
      }
    }
    __syncthreads();
    #endif

    // Updating number of ADADELTA iterations (energy evaluations)
    iteration_cnt = iteration_cnt + 1;
    if (threadIdx.x == 0) {
      // // sans check12 ==>
      // if (blockIdx.x == 0)
      // {
      //   printf("Previous drug energy: %f, drug energy: %f, best drug energy: %f, rho: %d\n", previousEnergy, energy2, bestDrugEnergy, rho);
      //   previousEnergy = energy2;
      // }
      // // <== sans check12
      if(energy2 < bestDrugEnergy) // sans check12
      {
        bestDrugEnergy = energy2; // sans check12
#ifdef ADADELTA_AUTOSTOP
        cons_succ++;
        cons_fail = 0;
#endif
    }
#ifdef ADADELTA_AUTOSTOP
      else
      {
        cons_succ = 0;
        cons_fail++;
      }
#endif

      #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
      printf("%20s %10.6f\n", "new.energy: ", energy);
      #endif

      #if defined (DEBUG_ENERGY_ADADELTA)
      printf("%-18s [%-5s]---{%-5s}   [%-10.7f]---{%-10.7f}\n", "-ENERGY-KERNEL7-", "GRIDS", "INTRA", partial_interE[0], partial_intraE[0]);
      #endif
#ifdef ADADELTA_AUTOSTOP
      // // sans check12 ==>
      // if (blockIdx.x == 0)
      // {
      //   printf("cons_succ: %d, cons_fail: %d\n", cons_succ, cons_fail);
      // }
      // // <== sans check12
      if (cons_succ >= 4)
      {
        // rho *= LS_EXP_FACTOR; // sans
        rho -= 1; // sans
        cons_succ = 0;
      }
      else
      {
        if (cons_fail >= 4)
        {
          // rho *= LS_CONT_FACTOR; // sans
          rho += 1; // sans
          cons_fail = 0;
        }
      }
#endif
    }
      __syncthreads(); // making sure that iteration_cnt is up-to-date
#ifdef ADADELTA_AUTOSTOP
    if((tooShortDistNum != 0) && (rho > cData.dockpars.collRhoLimit))
    {
      break;
    }
//   } while ((iteration_cnt < cData.dockpars.max_num_of_iters)  && (rho > 0.01f)); // sans
  } while (( tooShortDistNum != 0 || tooFarDistNum != 0) && iteration_cnt < halfMaxIter && rho < rhoL); // sans
#else
  } while (iteration_cnt < halfMaxIter);
#endif

  __syncthreads();
  for(uint32_t i = threadIdx.x;
               i < cData.dockpars.num_of_genes;
               i+= blockDim.x)
  {
    gradient[i]        = 0.0f;
    square_gradient[i] = 0.0f;
    square_delta[i]    = 0.0f;
  }

#ifdef ADADELTA_AUTOSTOP
  if (threadIdx.x == 0) {
    // rho = 1.0f; // sans
    // sans ==>
    rho  = 0;
    if(tooShortDistNum == 0 && tooFarDistNum == 0)
    {
      rhoL = cData.dockpars.compRhoLimit;
    }
    else
    {
      rhoL = -1;
    }
    // <== sans
    cons_succ = 0;
    cons_fail = 0;
  }
#endif
  __syncthreads();
// <== sans ligand
//   long long int time5-2 = clock64(); // sans clock
  do {
    // Printing number of ADADELTA iterations
    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
    if (threadIdx.x == 0) {
      #if defined (DEBUG_ADADELTA_MINIMIZER)
      printf("%s\n", "----------------------------------------------------------");
      #endif
      printf("%-15s %-3u ", "# ADADELTA iteration: ", iteration_cnt);
    }
    #endif

    // =============================================================
    // =============================================================
    // =============================================================
    // Calculating energy & gradient
    __syncthreads();

// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy3: %d %f run_id:%d\n", blockIdx.x, energy, run_id);
// }
    gpu_calc_energrad(
                      // Some OpenCL compilers don't allow declaring
                      // local variables within non-kernel functions.
                      // These local variables must be declared in a kernel,
                      // and then passed to non-kernel functions.
                      genotype,
                      energy,
//                       run_id, // sans useless
                      calc_coords,
                      #if defined (DEBUG_ENERGY_KERNEL)
                      interE,
                      intraE,
                      #endif
                      // Gradient-related arguments
                      // Calculate gradients (forces) for intermolecular energy
                      // Derived from autodockdev/maps.py
                      cartesian_gradient,
                      gradient,
                      &sFloatAccumulator
                      // sans ==>
//                       atom_charges_const,
//                       atom_types_const,
//                       atom_types_map_const,
//                       ignore_inter_const,
//
//                       intraE_contributors_const,

//                       rotlist_const,
//
//                       ref_coords_const,
//                       rotbonds_moving_vectors_const,
//                       rotbonds_unit_vectors_const,
//
//                       rotbonds,
//                       rotbonds_atoms,
//                       num_rotating_atoms_per_rotbond
                      // <== sans
                     );
    // sans check12 ==>
// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy4: %d %f run_id:%d\n", blockIdx.x, energy, run_id);
// }
    // if(previousEnergy != 0.0)
    // {
    //   printf("sans check12 %d %d %f %f %f\n", blockIdx.x, threadIdx.x, energy, previousEnergy, energy - previousEnergy);
    // }
    // previousEnergy = energy;
    // <== sans check12

    // =============================================================
    // =============================================================
    // =============================================================
    #if defined (DEBUG_ENERGY_ADADELTA)
    if (threadIdx.x == 0) {
      #if defined (PRINT_ADADELTA_ENERGIES)
      printf("\n");
      printf("%-10s %-10.6f \n", "intra: ",  intraE);
      printf("%-10s %-10.6f \n", "grids: ",  interE);
      printf("%-10s %-10.6f \n", "Energy: ", intraE + interE);
      #endif

      #if defined (PRINT_ADADELTA_GENES_AND_GRADS)
      for(uint i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %13s %5s %15s %15s\n", "gene_id", "gene.value", "|", "gene.grad", "(autodockdevpy units)");
        }
        printf("%13u %13.6f %5s %15.6f %15.6f\n", i, genotype[i], "|", gradient[i], (i<3)? (gradient[i]/0.375f):(gradient[i]*180.0f/PI_FLOAT));
      }
      #endif

      #if defined (PRINT_ADADELTA_ATOMIC_COORDS)
      for(uint i = 0; i < cData.dockpars.num_of_atoms; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%s\n", "Coordinates calculated by calcenergy.cl");
          printf("%12s %12s %12s %12s\n", "atom_id", "coords.x", "coords.y", "coords.z");
        }
        printf("%12u %12.6f %12.6f %12.6f\n", i, calc_coords_x[i], calc_coords_y[i], calc_coords_z[i]);
      }
      printf("\n");
      #endif
    }
    __syncthreads();
    #endif // DEBUG_ENERGY_ADADELTA

    for(int i = threadIdx.x;
            i < cData.dockpars.num_of_genes;
            i+= blockDim.x)
    {
      if (energy < best_energy) // we need to be careful not to change best_energy until we had a chance to update the whole array
        best_genotype[i] = genotype[i];

      // Accumulating gradient^2 (eq.8 in the paper)
      // square_gradient corresponds to E[g^2]
      square_gradient[i] = RHO * square_gradient[i] + (1.0f - RHO) * gradient[i] * gradient[i];

      // Computing update (eq.9 in the paper)
      // float delta = -1.0f * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON)); // sans check12
      // sans check12 ==>
      // add step size
      float delta = cData.dockpars.negaCompStepSize * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON));
      // <== sans check12

      // Accumulating update^2
      // square_delta corresponds to E[dx^2]
      square_delta[i] = RHO * square_delta[i] + (1.0f - RHO) * delta * delta;

      // Applying update
      genotype[i] += delta;
    }
    __syncthreads();

    #if defined (DEBUG_SQDELTA_ADADELTA)
    if (/*(get_group_id(0) == 0) &&*/ (threadIdx.x == 0)) {
      for(int i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %20s %15s %15s %15s\n", "gene", "sq_grad", "delta", "sq_delta", "new.genotype");
        }
        printf("%13u %20.6f %15.6f %15.6f %15.6f\n", i, square_gradient[i], delta[i], square_delta[i], genotype[i]);
      }
    }
    __syncthreads();
    #endif

    // Updating number of ADADELTA iterations (energy evaluations)
    iteration_cnt = iteration_cnt + 1;
    if (threadIdx.x == 0) {
      // // sans check12 ==>
      // if (blockIdx.x == 0)
      // {
      //   printf("complex energy: %f, best complex energy: %f, rho: %d\n", energy, best_energy, rho);
      // }
      // // <== sans check12
      if (energy < best_energy)
      {
        best_energy = energy;
#ifdef ADADELTA_AUTOSTOP
        cons_succ++;
        cons_fail = 0;
#endif
      }
#ifdef ADADELTA_AUTOSTOP
      else
      {
        // if(best_energy == INFINITY)
        // {
        //   printf("sans check13 INFINITY 1 entity_id: %d, en %f, b_en %f\n", entity_id, energy, best_energy);
        // }
        // if(energy == INFINITY)
        // {
        //   printf("sans check13 INFINITY 2 entity_id: %d, en %f, b_en %f\n", entity_id, energy, best_energy);
        // }
        cons_succ = 0;
        cons_fail++;
      }
#endif

      #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
      printf("%20s %10.6f\n", "new.energy: ", energy);
      #endif

      #if defined (DEBUG_ENERGY_ADADELTA)
      printf("%-18s [%-5s]---{%-5s}   [%-10.7f]---{%-10.7f}\n", "-ENERGY-KERNEL7-", "GRIDS", "INTRA", partial_interE[0], partial_intraE[0]);
      #endif
#ifdef ADADELTA_AUTOSTOP
      if (cons_succ >= 4)
      {
        // rho *= LS_EXP_FACTOR; // sans
        rho -= 1; // sans
        cons_succ = 0;
      }
      else
      {
        if (cons_fail >= 4)
        {
          // rho *= LS_CONT_FACTOR; // sans
          rho += 1; // sans
          cons_fail = 0;
        }
      }
#endif
    }
      __syncthreads(); // making sure that iteration_cnt is up-to-date
#ifdef ADADELTA_AUTOSTOP
  // } while ((iteration_cnt < cData.dockpars.max_num_of_iters)  && (rho > 0.01f)); // sans
  } while ((iteration_cnt < maxIter)  && (rho < rhoL)); // sans
#else
  } while (iteration_cnt < maxIter); // sans
#endif
//   long long int time6 = clock64(); // sans clock
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
// __syncthreads(); // sans check
  // Mapping torsion angles
  for (uint32_t gene_counter = threadIdx.x+3;
                gene_counter < cData.dockpars.num_of_genes;
                gene_counter += blockDim.x)
  {
    map_angle(best_genotype[gene_counter]);
  }

  // Updating old offspring in population
  __syncthreads();

//   offset = (run_id * cData.dockpars.pop_size + entity_id) * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
  offset = entity_id * genotype_length_in_globmem; // sans
  for (uint gene_counter = threadIdx.x;
            gene_counter < cData.dockpars.num_of_genes;
            gene_counter+= blockDim.x)
  {
    pMem_conformations_next[gene_counter + offset] = best_genotype[gene_counter];
  }
//   __syncthreads(); // sans check

  // Updating eval counter and energy
//   long long int time7 = clock64(); // sans clock
  if (threadIdx.x == 0) {
    /* sans
    cData.pMem_evals_of_new_entities[run_id * cData.dockpars.pop_size + entity_id] += iteration_cnt;
    pMem_energies_next[run_id * cData.dockpars.pop_size + entity_id] = best_energy;
    sans */
    // sans ==>
    pMem_evals_of_new_entities[entity_id] += iteration_cnt;
    // if(best_energy == INFINITY)
    // {
    //   printf("sans check13 INFINITY 3 entity_id: %d\n", entity_id);
    // }
    pMem_energies_next[entity_id] = best_energy;
    // <== sans

    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
    printf("\n");
    printf("Termination criteria: ( #adadelta-iters >= %-3u )\n", dockpars_max_num_of_iters);
    printf("-------> End of ADADELTA minimization cycle, num of energy evals: %u, final energy: %.6f\n", iteration_cnt, best_energy);
    #endif
  }

  // sans ==>
  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
  // <== sans

//   if (threadIdx.x == 0 && blockIdx.x == 0)
//   {
//     printf("\nSANS TIME\n");                      // sans clock
//     printf("TIME1:   %lld\n", time2 - time1);     // sans clock
//     printf("TIME2:   %lld\n", time3 - time2);     // sans clock
//     printf("TIME3:   %lld\n", time4 - time3);     // sans clock
//     printf("TIME4-1: %lld\n", time5-1 -   time4); // sans clock
//     printf("TIME4-2: %lld\n", time5-2 - time5-1); // sans clock
//     printf("TIME5:   %lld\n", time6 -   time5-2); // sans clock
//     printf("TIME6:   %lld\n", time7 - time6);     // sans clock
//   }
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next2[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
  // <== sans check
}


// sans ==>
__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_gradient_minAD_ligand_kernel(
                          float* pMem_conformations_next,
                          float* pMem_energies_next,
                          // sans ==>
                          int   *pMem_energyToConforIndex,
                          int   *pMem_evals_of_new_entities,
                 // volatile float *gFloatBuff,
                          float *gFloatBuff,
                          size_t elemNum_gFloatBuffBase,
                          int    run_win_id,
                          int    run_id

//                           float    *atom_charges_const,
//                           uint32_t *atom_types_const,
//                           uint32_t *atom_types_map_const,
//                           char     *ignore_inter_const,
//
//                           uint32_t *intraE_contributors_const,

//                           int      *rotlist_const,
//
//                           float    *ref_coords_const,
//                           float    *rotbonds_moving_vectors_const,
//                           float    *rotbonds_unit_vectors_const,
//
//                           int      *rotbonds,
//                           int      *rotbonds_atoms,
//                           int      *num_rotating_atoms_per_rotbond
                          // <== sans
                         )
// The GPU global function performs gradient-based minimization on (some) entities of conformations_next.
// The number of OpenCL compute units (CU) which should be started equals to num_of_minEntities*num_of_runs.
// This way the first num_of_lsentities entity of each population will be subjected to local search
// (and each CU carries out the algorithm for one entity).
// Since the first entity is always the best one in the current population,
// it is always tested according to the ls probability, and if it not to be
// subjected to local search, the entity with ID num_of_lsentities is selected instead of the first one (with ID 0).
{
//   long long int time1 = clock64(); // sans clock
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------

  // Determining entity, and its run, energy, and genotype
//   int   run_id = blockIdx.x / cData.dockpars.num_of_lsentities; // sans
  float energy;
  float energy2; // sans check12
  int   tooShortDistNum; // sans
  int   tooFarDistNum;   // sans
  // Energy may go up, so we keep track of the best energy ever calculated.
  // Then, we return the genotype corresponding
  // to the best observed energy, i.e. "best_genotype"
  __shared__ int entity_id;
  __shared__ float best_energy;
  // __shared__ float previousEnergy; // sans check12
  __shared__ float bestDrugEnergy; // sans check12
  __shared__ float sFloatAccumulator;
  __shared__ int sIntAccumulator; // sans

//   volatile __shared__ float sFloatBuff[6 * MAX_NUM_OF_ATOMS + 5 * ACTUAL_GENOTYPE_LENGTH]; // sans
// calc_coords cartesian_gradient genotype best_genotype gradient square_delta square_gradient
// GENOTYPE_LENGTH_IN_GLOBMEM 64
// MAX_NUM_OF_ROTBONDS 57
// ACTUAL_GENOTYPE_LENGTH (MAX_NUM_OF_ROTBONDS+6)

  // Ligand-atom position and partial energies
  // sans ==>
  int actual_genotype_length = cData.dockpars.num_of_genes;
  int genotype_length_in_globmem = actual_genotype_length + 1;

  // sans check ==>
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next1a[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
//   __syncthreads();
  // <== sans check
  // sans check ==>
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next1a[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
//   __syncthreads();
  // <== sans check
  float3 *calc_coords;
  float  *genotype;
#ifdef FLOAT_GRADIENTS
    float3* cartesian_gradient;
#else
    int3*   cartesian_gradient;
#endif
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (6 * NUM_OF_THREADS_PER_BLOCK + 5 * actual_genotype_length));
#ifdef FLOAT_GRADIENTS
    cartesian_gradient = (float3*)(calc_coords + NUM_OF_THREADS_PER_BLOCK);
#else
    cartesian_gradient = (int3*)  (calc_coords + NUM_OF_THREADS_PER_BLOCK);
#endif
    genotype = (float*)(cartesian_gradient + NUM_OF_THREADS_PER_BLOCK); // so far used 3*2*NUM_OF_THREADS_PER_BLOCK
  }
  else
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (6 * cData.dockpars.num_of_atoms + 5 * actual_genotype_length));
#ifdef FLOAT_GRADIENTS
    cartesian_gradient = (float3*)(calc_coords + cData.dockpars.num_of_atoms);
#else
    cartesian_gradient = (int3*)  (calc_coords + cData.dockpars.num_of_atoms);
#endif
    genotype = (float*)(cartesian_gradient + cData.dockpars.num_of_atoms); // so far used 3*2*cData.dockpars.num_of_atoms
  }
  // <== sans


//   float3* calc_coords = (float3*)sFloatBuff; // sans
  // Gradient of the intermolecular energy per each ligand atom
  // Also used to store the accummulated gradient per each ligand atom
// #ifdef FLOAT_GRADIENTS
//   float3* cartesian_gradient = (float3*)(calc_coords + MAX_NUM_OF_ATOMS); // sans
// #else
//   int3* cartesian_gradient = (int3*)(calc_coords + MAX_NUM_OF_ATOMS); // sans
// #endif
  // Genotype pointers
//   float* genotype = (float*)(cartesian_gradient + MAX_NUM_OF_ATOMS); // so far used 3*2*MAX_NUM_OF_ATOMS // sans
//   float* best_genotype = (float*) (genotype + ACTUAL_GENOTYPE_LENGTH); // sans
  float* best_genotype = (float*) (genotype + actual_genotype_length); // sans

  // Partial results of the gradient step
//   float* gradient = best_genotype + ACTUAL_GENOTYPE_LENGTH; // sans
  float* gradient = best_genotype + actual_genotype_length; // sans

  // Squared updates E[dx^2]
//   float* square_delta = gradient + ACTUAL_GENOTYPE_LENGTH; // sans
  float* square_delta = gradient + actual_genotype_length; // sans

  // Vector for storing squared gradients E[g^2]
//   float* square_gradient = square_delta + ACTUAL_GENOTYPE_LENGTH; // so far used 5*ACTUAL_GENOTYPE_LENGTH // sans
  float* square_gradient = square_delta + actual_genotype_length; // so far used 5*actual_genotype_length // sans

  // Iteration counter for the minimizer
  uint32_t iteration_cnt = 0;

  if (threadIdx.x == 0)
  {
    // Since entity 0 is the best one due to elitism,
    // it should be subjected to random selection
    entity_id = blockIdx.x % cData.dockpars.num_of_lsentities; // sans this should always be equal to blockIdx.x
    if (entity_id == 0) {
      // If entity 0 is not selected according to LS-rate,
      // choosing another entity
//       if (100.0f*gpu_randf(cData.pMem_prng_states) > cData.dockpars.lsearch_rate) { // sans check
      if (100.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x) > cData.dockpars.lsearch_rate) { // sans check
        entity_id = cData.dockpars.num_of_lsentities; // AT - Should this be (uint)(dockpars_pop_size * gpu_randf(dockpars_prng_states))?
      }
    }

    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION) // sans
//     printf("\n");
    printf("-------> Start of ADADELTA minimization cycle\n");
    printf("%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n%20s %6u\n", "threadIdx.x: ", threadIdx.x, "blockDim.x: ", blockDim.x, "blockIdx.x: ", blockIdx.x, "gridDim.x: ", gridDim.x, "run_id: ", run_id, "entity_id: ", entity_id, "num_of_lsentities: ", cData.dockpars.num_of_lsentities, "pop_size: ", cData.dockpars.pop_size, "num_of_genes", cData.dockpars.num_of_genes);
// sans ./host/src/performdocking.cpp.Cuda:483: cData.dockpars.num_of_lsentities        = (unsigned int) (mypars->lsearch_rate/100.0*mypars->pop_size + 0.5);
// sans ./host/src/performdocking.cpp.Cuda:466: cData.dockpars.num_of_genes                 = myligand_reference.num_of_rotbonds + 6;
// sans blockNum = blocksPerGridForEachGradMinimizerEntity = cData.dockpars.num_of_lsentities * mypars->num_of_runs;
// sans num_of_lsentities = ga_pop_size
// sans num_of_runs = ga_run
//     printf("%20s %6u\n", "run_id: ", run_id);
//     printf("%20s %6u\n", "entity_id: ", entity_id);
//     printf("%20s %6u\n", "num_of_lsentities: ", cData.dockpars.num_of_lsentities);
//     printf("\n");
//     printf("%20s \n", "LGA genotype: ");
//     printf("%20s %.6f\n", "initial energy: ", energy);
    #endif // sans
  }
  __syncthreads();
//   energy = pMem_energies_next[run_id * cData.dockpars.pop_size + entity_id]; // sans
  // energy = pMem_energies_next[entity_id]; // sans
// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy1: %d %f\n", blockIdx.x, pMem_energies_next[blockIdx.x]);
// }

//   int offset = (run_id * cData.dockpars.pop_size + entity_id) * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
  int offset = entity_id * genotype_length_in_globmem; // sans
//   long long int time2 = clock64(); // sans clock
  for (int i = threadIdx.x ; i < cData.dockpars.num_of_genes; i += blockDim.x)
// sans cData.dockpars.num_of_genes = myligand_reference.num_of_rotbonds + 6
  {
    genotype[i] = pMem_conformations_next[offset + i];
  }
//   long long int time3 = clock64(); // sans clock

  // sans check ==>
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next1[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
  // <== sans check

  // -------------------------------------------------------------------
  // Calculate gradients (forces) for intermolecular energy
  // Derived from autodockdev/maps.py
  // -------------------------------------------------------------------

  #if defined (DEBUG_ENERGY_KERNEL)
  float interE;
  float intraE;
  #endif

  // Update vector, i.e., "delta".
  // It is added to the genotype to create the next genotype.
  // E.g. in steepest descent "delta" is -1.0 * stepsize * gradient

  // Asynchronous copy should be finished by here
  __syncthreads();

  // Initializing vectors
  for(uint32_t i = threadIdx.x;
               i < cData.dockpars.num_of_genes;
               i+= blockDim.x)
  {
    gradient[i]        = 0.0f;
    square_gradient[i] = 0.0f;
    square_delta[i]    = 0.0f;
    best_genotype[i] = genotype[i];
  }
//   long long int time4 = clock64(); // sans clock

  // Initializing best energy
  if (threadIdx.x == 0) {
    best_energy = INFINITY;
    // previousEnergy = 0.0; // sans check12
    bestDrugEnergy = INFINITY; // sans check12
  }

#ifdef ADADELTA_AUTOSTOP
  // __shared__ float rho; // sans
  __shared__ int rhoL; // sans
  __shared__ int rho;  // sans
  __shared__ int cons_succ;
  __shared__ int cons_fail;
  if (threadIdx.x == 0) {
    // rho = 1.0f; // sans
    rho = 0;   // sans
    rhoL = cData.dockpars.ligaRhoLimit; // sans
    cons_succ = 0;
    cons_fail = 0;
  }
#endif
  __syncthreads(); // sans check12

  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // Perform adadelta iterations

  // The termination criteria is based on
  // a maximum number of iterations, and
  // the minimum step size allowed for single-floating point numbers
  // (IEEE-754 single float has a precision of about 6 decimal digits)
//   long long int time5-1 = clock64(); // sans clock
// sans ligand ==>
  uint halfMaxIter;
  uint maxIter;
  maxIter = cData.dockpars.max_num_of_iters;
  halfMaxIter = maxIter / 2;
  if(halfMaxIter == 0)
  {
    halfMaxIter++;
    maxIter++;
  }

  do {
    // Printing number of ADADELTA iterations
    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
    if (threadIdx.x == 0) {
      #if defined (DEBUG_ADADELTA_MINIMIZER)
      printf("%s\n", "----------------------------------------------------------");
      #endif
      printf("%-15s %-3u ", "# ADADELTA iteration: ", iteration_cnt);
    }
    #endif

    // =============================================================
    // =============================================================
    // =============================================================
    // Calculating energy & gradient
    __syncthreads();

// for(int i = threadIdx.x;
//             i < cData.dockpars.num_of_genes;
//             i+= blockDim.x)
// {
//   if(isnan(genotype[i]))
//   {
//     printf("sans check13 NAN kernel_ad1 %d %f run_id:%d\n", i, genotype[i], run_id);
//   }
// }
    if(blockIdx.x >= cData.pMem_reservedIdvNum[run_id]) // sans
    { // sans
      gpu_calc_energrad_ligand
      (
                      // Some OpenCL compilers don't allow declaring
                      // local variables within non-kernel functions.
                      // These local variables must be declared in a kernel,
                      // and then passed to non-kernel functions.
                      genotype,
                      energy,
                      energy2, // sans check12
                      tooShortDistNum, // sans
                      tooFarDistNum,   // sans
//                       run_id, // sans useless
                      calc_coords,
                      #if defined (DEBUG_ENERGY_KERNEL)
                      interE,
                      intraE,
                      #endif
                      // Gradient-related arguments
                      // Calculate gradients (forces) for intermolecular energy
                      // Derived from autodockdev/maps.py
                      cartesian_gradient,
                      gradient,
                      &sFloatAccumulator,
                      &sIntAccumulator
                      // sans ==>
//                       atom_charges_const,
//                       atom_types_const,
//                       atom_types_map_const,
//                       ignore_inter_const,
//
//                       intraE_contributors_const,

//                       rotlist_const,
//
//                       ref_coords_const,
//                       rotbonds_moving_vectors_const,
//                       rotbonds_unit_vectors_const,
//
//                       rotbonds,
//                       rotbonds_atoms,
//                       num_rotating_atoms_per_rotbond
                      // <== sans
      );
    // sans ==>
    }
    else
    {
      rhoL = -1;
      tooShortDistNum = 0;
      tooFarDistNum = 0;
    }
    // <== sans
    // sans check12 ==>
// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy2: %d %f run_id:%d\n", blockIdx.x, energy, run_id);
// }
    // if(previousEnergy != 0.0)
    // {
    //   printf("sans check12 %d %d %f %f %f\n", blockIdx.x, threadIdx.x, energy, previousEnergy, energy - previousEnergy);
    // }
    // previousEnergy = energy;
    // <== sans check12

    // =============================================================
    // =============================================================
    // =============================================================
    #if defined (DEBUG_ENERGY_ADADELTA)
    if (threadIdx.x == 0) {
      #if defined (PRINT_ADADELTA_ENERGIES)
      printf("\n");
      printf("%-10s %-10.6f \n", "intra: ",  intraE);
      printf("%-10s %-10.6f \n", "grids: ",  interE);
      printf("%-10s %-10.6f \n", "Energy: ", intraE + interE);
      #endif

      #if defined (PRINT_ADADELTA_GENES_AND_GRADS)
      for(uint i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %13s %5s %15s %15s\n", "gene_id", "gene.value", "|", "gene.grad", "(autodockdevpy units)");
        }
        printf("%13u %13.6f %5s %15.6f %15.6f\n", i, genotype[i], "|", gradient[i], (i<3)? (gradient[i]/0.375f):(gradient[i]*180.0f/PI_FLOAT));
      }
      #endif

      #if defined (PRINT_ADADELTA_ATOMIC_COORDS)
      for(uint i = 0; i < cData.dockpars.num_of_atoms; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%s\n", "Coordinates calculated by calcenergy.cl");
          printf("%12s %12s %12s %12s\n", "atom_id", "coords.x", "coords.y", "coords.z");
        }
        printf("%12u %12.6f %12.6f %12.6f\n", i, calc_coords_x[i], calc_coords_y[i], calc_coords_z[i]);
      }
      printf("\n");
      #endif
    }
    __syncthreads();
    #endif // DEBUG_ENERGY_ADADELTA

    for(int i = threadIdx.x;
            i < cData.dockpars.num_of_genes;
            i+= blockDim.x)
    {
      // if (energy < best_energy) // we need to be careful not to change best_energy until we had a chance to update the whole array
      //   best_genotype[i] = genotype[i];

      // Accumulating gradient^2 (eq.8 in the paper)
      // square_gradient corresponds to E[g^2]
      square_gradient[i] = RHO * square_gradient[i] + (1.0f - RHO) * gradient[i] * gradient[i];

      // Computing update (eq.9 in the paper)
      // float delta = -1.0f * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON)); // sans check12
      // sans check12 ==>
      // add step size
      float delta = cData.dockpars.negaLigaStepSize * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON));
      // <== sans check12

      // Accumulating update^2
      // square_delta corresponds to E[dx^2]
      square_delta[i] = RHO * square_delta[i] + (1.0f - RHO) * delta * delta;

      // Applying update
// if(isnan(genotype[i]))
// {
//   printf("sans check13 NAN kernel_ad2 %d %f\n", i, genotype[i]);
// }
// float ggg = genotype[i];
      genotype[i] += delta;
// if(isnan(genotype[i]))
// {
//   printf("sans check13 NAN kernel_ad %d ori%f %f d%f s%f g%f\n", i, ggg, genotype[i], delta, square_gradient[i], gradient[i]);
// }
    }
    __syncthreads();

    #if defined (DEBUG_SQDELTA_ADADELTA)
    if (/*(get_group_id(0) == 0) &&*/ (threadIdx.x == 0)) {
      for(int i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %20s %15s %15s %15s\n", "gene", "sq_grad", "delta", "sq_delta", "new.genotype");
        }
        printf("%13u %20.6f %15.6f %15.6f %15.6f\n", i, square_gradient[i], delta[i], square_delta[i], genotype[i]);
      }
    }
    __syncthreads();
    #endif

    // Updating number of ADADELTA iterations (energy evaluations)
    iteration_cnt = iteration_cnt + 1;
    if (threadIdx.x == 0) {
      // // sans check12 ==>
      // if (blockIdx.x == 0)
      // {
      //   printf("Previous drug energy: %f, drug energy: %f, best drug energy: %f, rho: %d\n", previousEnergy, energy2, bestDrugEnergy, rho);
      //   previousEnergy = energy2;
      // }
      // // <== sans check12
      if(energy2 < bestDrugEnergy) // sans check12
      {
        bestDrugEnergy = energy2; // sans check12
#ifdef ADADELTA_AUTOSTOP
        cons_succ++;
        cons_fail = 0;
#endif
      }
#ifdef ADADELTA_AUTOSTOP
      else
      {
        cons_succ = 0;
        cons_fail++;
      }
#endif

      #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
      printf("%20s %10.6f\n", "new.energy: ", energy);
      #endif

      #if defined (DEBUG_ENERGY_ADADELTA)
      printf("%-18s [%-5s]---{%-5s}   [%-10.7f]---{%-10.7f}\n", "-ENERGY-KERNEL7-", "GRIDS", "INTRA", partial_interE[0], partial_intraE[0]);
      #endif
#ifdef ADADELTA_AUTOSTOP
      // // sans check12 ==>
      // if (blockIdx.x == 0)
      // {
      //   printf("cons_succ: %d, cons_fail: %d\n", cons_succ, cons_fail);
      // }
      // // <== sans check12
      if (cons_succ >= 4)
      {
        // rho *= LS_EXP_FACTOR; // sans
        rho -= 1; // sans
        cons_succ = 0;
      }
      else
      {
        if (cons_fail >= 4)
        {
          // rho *= LS_CONT_FACTOR; // sans
          rho += 1; // sans
          cons_fail = 0;
        }
      }
#endif
    }
      __syncthreads(); // making sure that iteration_cnt is up-to-date
#ifdef ADADELTA_AUTOSTOP
    if((tooShortDistNum != 0) && (rho > cData.dockpars.collRhoLimit))
    {
      break;
    }
//   } while ((iteration_cnt < cData.dockpars.max_num_of_iters)  && (rho > 0.01f)); // sans
  } while (( tooShortDistNum != 0 || tooFarDistNum != 0) && iteration_cnt < halfMaxIter && rho < rhoL); // sans
#else
  } while (iteration_cnt < halfMaxIter);
#endif

  __syncthreads();
  for(uint32_t i = threadIdx.x;
               i < cData.dockpars.num_of_genes;
               i+= blockDim.x)
  {
    gradient[i]        = 0.0f;
    square_gradient[i] = 0.0f;
    square_delta[i]    = 0.0f;
  }

#ifdef ADADELTA_AUTOSTOP
  if (threadIdx.x == 0) {
    // rho = 1.0f; // sans
    // sans ==>
    rho  = 0;
    if(tooShortDistNum == 0 && tooFarDistNum == 0)
    {
      rhoL = cData.dockpars.compRhoLimit;
    }
    else
    {
      rhoL = -1;
    }
    // <== sans
    cons_succ = 0;
    cons_fail = 0;
  }
#endif
  __syncthreads();
// <== sans ligand
//   long long int time5-2 = clock64(); // sans clock
  do {
    // Printing number of ADADELTA iterations
    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
    if (threadIdx.x == 0) {
      #if defined (DEBUG_ADADELTA_MINIMIZER)
      printf("%s\n", "----------------------------------------------------------");
      #endif
      printf("%-15s %-3u ", "# ADADELTA iteration: ", iteration_cnt);
    }
    #endif

    // =============================================================
    // =============================================================
    // =============================================================
    // Calculating energy & gradient
    __syncthreads();

// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy3: %d %f run_id:%d\n", blockIdx.x, energy, run_id);
// }
    gpu_calc_energrad_ligand2(
                      // Some OpenCL compilers don't allow declaring
                      // local variables within non-kernel functions.
                      // These local variables must be declared in a kernel,
                      // and then passed to non-kernel functions.
                      genotype,
                      energy,
//                       run_id, // sans useless
                      calc_coords,
                      #if defined (DEBUG_ENERGY_KERNEL)
                      interE,
                      intraE,
                      #endif
                      // Gradient-related arguments
                      // Calculate gradients (forces) for intermolecular energy
                      // Derived from autodockdev/maps.py
                      cartesian_gradient,
                      gradient,
                      &sFloatAccumulator
                      // sans ==>
//                       atom_charges_const,
//                       atom_types_const,
//                       atom_types_map_const,
//                       ignore_inter_const,
//
//                       intraE_contributors_const,

//                       rotlist_const,
//
//                       ref_coords_const,
//                       rotbonds_moving_vectors_const,
//                       rotbonds_unit_vectors_const,
//
//                       rotbonds,
//                       rotbonds_atoms,
//                       num_rotating_atoms_per_rotbond
                      // <== sans
                     );
    // sans check12 ==>
// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
//   printf("sans kernel_ad energy4: %d %f run_id:%d\n", blockIdx.x, energy, run_id);
// }
    // if(previousEnergy != 0.0)
    // {
    //   printf("sans check12 %d %d %f %f %f\n", blockIdx.x, threadIdx.x, energy, previousEnergy, energy - previousEnergy);
    // }
    // previousEnergy = energy;
    // <== sans check12

    // =============================================================
    // =============================================================
    // =============================================================
    #if defined (DEBUG_ENERGY_ADADELTA)
    if (threadIdx.x == 0) {
      #if defined (PRINT_ADADELTA_ENERGIES)
      printf("\n");
      printf("%-10s %-10.6f \n", "intra: ",  intraE);
      printf("%-10s %-10.6f \n", "grids: ",  interE);
      printf("%-10s %-10.6f \n", "Energy: ", intraE + interE);
      #endif

      #if defined (PRINT_ADADELTA_GENES_AND_GRADS)
      for(uint i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %13s %5s %15s %15s\n", "gene_id", "gene.value", "|", "gene.grad", "(autodockdevpy units)");
        }
        printf("%13u %13.6f %5s %15.6f %15.6f\n", i, genotype[i], "|", gradient[i], (i<3)? (gradient[i]/0.375f):(gradient[i]*180.0f/PI_FLOAT));
      }
      #endif

      #if defined (PRINT_ADADELTA_ATOMIC_COORDS)
      for(uint i = 0; i < cData.dockpars.num_of_atoms; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%s\n", "Coordinates calculated by calcenergy.cl");
          printf("%12s %12s %12s %12s\n", "atom_id", "coords.x", "coords.y", "coords.z");
        }
        printf("%12u %12.6f %12.6f %12.6f\n", i, calc_coords_x[i], calc_coords_y[i], calc_coords_z[i]);
      }
      printf("\n");
      #endif
    }
    __syncthreads();
    #endif // DEBUG_ENERGY_ADADELTA

    for(int i = threadIdx.x;
            i < cData.dockpars.num_of_genes;
            i+= blockDim.x)
    {
      if (energy < best_energy) // we need to be careful not to change best_energy until we had a chance to update the whole array
        best_genotype[i] = genotype[i];

      // Accumulating gradient^2 (eq.8 in the paper)
      // square_gradient corresponds to E[g^2]
      square_gradient[i] = RHO * square_gradient[i] + (1.0f - RHO) * gradient[i] * gradient[i];

      // Computing update (eq.9 in the paper)
      // float delta = -1.0f * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON)); // sans check12
      // sans check12 ==>
      // add step size
      float delta = cData.dockpars.negaCompStepSize * gradient[i] * sqrt((float)(square_delta[i] + EPSILON) / (float)(square_gradient[i] + EPSILON));
      // <== sans check12

      // Accumulating update^2
      // square_delta corresponds to E[dx^2]
      square_delta[i] = RHO * square_delta[i] + (1.0f - RHO) * delta * delta;

      // Applying update
      genotype[i] += delta;
    }
    __syncthreads();

    #if defined (DEBUG_SQDELTA_ADADELTA)
    if (/*(get_group_id(0) == 0) &&*/ (threadIdx.x == 0)) {
      for(int i = 0; i < cData.dockpars.num_of_genes; i++) {
        if (i == 0) {
          printf("\n%s\n", "----------------------------------------------------------");
          printf("%13s %20s %15s %15s %15s\n", "gene", "sq_grad", "delta", "sq_delta", "new.genotype");
        }
        printf("%13u %20.6f %15.6f %15.6f %15.6f\n", i, square_gradient[i], delta[i], square_delta[i], genotype[i]);
      }
    }
    __syncthreads();
    #endif

    // Updating number of ADADELTA iterations (energy evaluations)
    iteration_cnt = iteration_cnt + 1;
    if (threadIdx.x == 0) {
      // // sans check12 ==>
      // if (blockIdx.x == 0)
      // {
      //   printf("complex energy: %f, best complex energy: %f, rho: %d\n", energy, best_energy, rho);
      // }
      // // <== sans check12
      if (energy < best_energy)
      {
        best_energy = energy;
#ifdef ADADELTA_AUTOSTOP
        cons_succ++;
        cons_fail = 0;
#endif
      }
#ifdef ADADELTA_AUTOSTOP
      else
      {
        // if(best_energy == INFINITY)
        // {
        //   printf("sans check13 INFINITY 1 entity_id: %d, en %f, b_en %f\n", entity_id, energy, best_energy);
        // }
        // if(energy == INFINITY)
        // {
        //   printf("sans check13 INFINITY 2 entity_id: %d, en %f, b_en %f\n", entity_id, energy, best_energy);
        // }
        cons_succ = 0;
        cons_fail++;
      }
#endif

      #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
      printf("%20s %10.6f\n", "new.energy: ", energy);
      #endif

      #if defined (DEBUG_ENERGY_ADADELTA)
      printf("%-18s [%-5s]---{%-5s}   [%-10.7f]---{%-10.7f}\n", "-ENERGY-KERNEL7-", "GRIDS", "INTRA", partial_interE[0], partial_intraE[0]);
      #endif
#ifdef ADADELTA_AUTOSTOP
      if (cons_succ >= 4)
      {
        // rho *= LS_EXP_FACTOR; // sans
        rho -= 1; // sans
        cons_succ = 0;
      }
      else
      {
        if (cons_fail >= 4)
        {
          // rho *= LS_CONT_FACTOR; // sans
          rho += 1; // sans
          cons_fail = 0;
        }
      }
#endif
    }
      __syncthreads(); // making sure that iteration_cnt is up-to-date
#ifdef ADADELTA_AUTOSTOP
  // } while ((iteration_cnt < cData.dockpars.max_num_of_iters)  && (rho > 0.01f)); // sans
  } while ((iteration_cnt < maxIter)  && (rho < rhoL)); // sans
#else
  } while (iteration_cnt < maxIter); // sans
#endif
//   long long int time6 = clock64(); // sans clock
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
// __syncthreads(); // sans check
  // Mapping torsion angles
  for (uint32_t gene_counter = threadIdx.x+3;
                gene_counter < cData.dockpars.num_of_genes;
                gene_counter += blockDim.x)
  {
    map_angle(best_genotype[gene_counter]);
  }

  // Updating old offspring in population
  __syncthreads();

//   offset = (run_id * cData.dockpars.pop_size + entity_id) * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
  offset = entity_id * genotype_length_in_globmem; // sans
  for (uint gene_counter = threadIdx.x;
            gene_counter < cData.dockpars.num_of_genes;
            gene_counter+= blockDim.x)
  {
    pMem_conformations_next[gene_counter + offset] = best_genotype[gene_counter];
  }
//   __syncthreads(); // sans check

  // Updating eval counter and energy
//   long long int time7 = clock64(); // sans clock
  if (threadIdx.x == 0) {
    /* sans
    cData.pMem_evals_of_new_entities[run_id * cData.dockpars.pop_size + entity_id] += iteration_cnt;
    pMem_energies_next[run_id * cData.dockpars.pop_size + entity_id] = best_energy;
    sans */
    // sans ==>
    pMem_evals_of_new_entities[entity_id] += iteration_cnt;
    // if(best_energy == INFINITY)
    // {
    //   printf("sans check13 INFINITY 3 entity_id: %d\n", entity_id);
    // }
    pMem_energies_next[entity_id] = best_energy;
    // <== sans

    #if defined (DEBUG_ADADELTA_MINIMIZER) || defined (PRINT_ADADELTA_MINIMIZER_ENERGY_EVOLUTION)
    printf("\n");
    printf("Termination criteria: ( #adadelta-iters >= %-3u )\n", dockpars_max_num_of_iters);
    printf("-------> End of ADADELTA minimization cycle, num of energy evals: %u, final energy: %.6f\n", iteration_cnt, best_energy);
    #endif
  }

  // sans ==>
  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
  // <== sans

//   if (threadIdx.x == 0 && blockIdx.x == 0)
//   {
//     printf("\nSANS TIME\n");                      // sans clock
//     printf("TIME1:   %lld\n", time2 - time1);     // sans clock
//     printf("TIME2:   %lld\n", time3 - time2);     // sans clock
//     printf("TIME3:   %lld\n", time4 - time3);     // sans clock
//     printf("TIME4-1: %lld\n", time5-1 -   time4); // sans clock
//     printf("TIME4-2: %lld\n", time5-2 - time5-1); // sans clock
//     printf("TIME5:   %lld\n", time6 -   time5-2); // sans clock
//     printf("TIME6:   %lld\n", time7 - time6);     // sans clock
//   }
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel_ad run_id/e/g/pMem_conformations_next2[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
  // <== sans check
}
// sans <==


__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
checkcheckKernelAD(
                          float* pMem_conformations_next,
                          float* pMem_energies_next,
                          // sans ==>
                          int   *pMem_energyToConforIndex,
                          int   *pMem_evals_of_new_entities,
                 // volatile float *gFloatBuff,
                          float *gFloatBuff,
                          size_t elemNum_gFloatBuffBase,
                          int    run_win_id,
                          int    run_id

//                           float    *atom_charges_const,
//                           uint32_t *atom_types_const,
//                           uint32_t *atom_types_map_const,
//                           char     *ignore_inter_const,
//
//                           uint32_t *intraE_contributors_const,

//                           int      *rotlist_const,
//
//                           float    *ref_coords_const,
//                           float    *rotbonds_moving_vectors_const,
//                           float    *rotbonds_unit_vectors_const,
//
//                           int      *rotbonds,
//                           int      *rotbonds_atoms,
//                           int      *num_rotating_atoms_per_rotbond
                          // <== sans
                         )
// The GPU global function performs gradient-based minimization on (some) entities of conformations_next.
// The number of OpenCL compute units (CU) which should be started equals to num_of_minEntities*num_of_runs.
// This way the first num_of_lsentities entity of each population will be subjected to local search
// (and each CU carries out the algorithm for one entity).
// Since the first entity is always the best one in the current population,
// it is always tested according to the ls probability, and if it not to be
// subjected to local search, the entity with ID num_of_lsentities is selected instead of the first one (with ID 0).
{

// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
//   for(int i = 0; i < cData.dockpars.pop_size; i++)
//   {
//     printf("sans kernel_ad run_win_id blockIdx.x energy2 %d %d %f\n", run_win_id, i, pMem_energies_next[i]);
//   }
// }

}

void gpu_gradient_minAD(
                        uint32_t blocks,
                        uint32_t threads,
                        /* sans
                        float*   pMem_conformations_next,
                        float*   pMem_energies_next
                        sans */
                        // sans ==>
                        bool     ligandOnly,
                        float  **pMem_conformations_next,
                        float  **pMem_energies_next,
                        int    **pMem_energyToConforIndex,
                        int    **pMem_evals_of_new_entities,
                        float   *gFloatBuff,
                        size_t   elemNum_gFloatBuffBase,
                        int      run_id_init,
                        int      run_id_end
                        // <== sans
                       )
{
//   gpu_gradient_minAD_kernel<<<blocks, threads>>>(pMem_conformations_next, pMem_energies_next); // sans
  // sans ==>
  int run_win_id = 0;
  if(ligandOnly == false)
  {
    for(int i = run_id_init; i < run_id_end; i++)
    {
      gpu_gradient_minAD_kernel<<<blocks, threads>>>(
        pMem_conformations_next[i],
        pMem_energies_next[i],
        pMem_energyToConforIndex[i],
        pMem_evals_of_new_entities[i],
        gFloatBuff,
        elemNum_gFloatBuffBase,
        run_win_id,
        i
      );
    }
  }
  else
  {
    for(int i = run_id_init; i < run_id_end; i++)
    {
      gpu_gradient_minAD_ligand_kernel<<<blocks, threads>>>(
        pMem_conformations_next[i],
        pMem_energies_next[i],
        pMem_energyToConforIndex[i],
        pMem_evals_of_new_entities[i],
        gFloatBuff,
        elemNum_gFloatBuffBase,
        run_win_id,
        i
      );
    }
  }
  // <== sans
  LAUNCHERROR("gpu_gradient_minAD_kernel");
#if 0
  cudaError_t status;
  status = cudaDeviceSynchronize();
  RTERROR(status, "gpu_gradient_minAD_kernel");
  status = cudaDeviceReset();
  RTERROR(status, "failed to shut down");
  exit(0);
#endif
}
