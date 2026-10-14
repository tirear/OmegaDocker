/*

OmegaDocker, a GPU-accelerated docking program developed for larger ligands, such as oligopeptides, oligonucleotides, and other oligomers.
Copyright (C) 2024 Gao Quan-Ze.

OmegaDocker is modified from AutoDock-GPU v1.5.2:
AutoDock-GPU, an OpenCL implementation of AutoDock 4.2 running a Lamarckian Genetic Algorithm
Copyright (C) 2017 TU Darmstadt, Embedded Systems and Applications Group, Germany. All rights reserved.
For some of the code, Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

AutoDock is a Trade Mark of the Scripps Research Institute.

SPDX-License-Identifier: LGPL-2.1-or-later

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



//#define DEBUG_ENERGY_KERNEL4

__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_gen_and_eval_newpops_kernel(
                                float* pMem_conformations_current,
                                float* pMem_energies_current,
                                float* pMem_conformations_next,
                                float* pMem_energies_next,
                                // sans ==>
                                int   *pMem_energyToConforIndex,
                                int   *pMem_evals_of_new_entities,
                                float *gFloatBuff,
                                size_t elemNum_gFloatBuffBase,
                                int    run_win_id,
                                int    run_id // sans check

//                                 float    *atom_charges_const,
//                                 uint32_t *atom_types_const,
//                                 uint32_t *atom_types_map_const,
//                                 char     *ignore_inter_const,
//
//                                 uint32_t *intraE_contributors_const,

//                                 int      *rotlist_const,
//
//                                 float    *ref_coords_const,
//                                 float    *rotbonds_moving_vectors_const,
//                                 float    *rotbonds_unit_vectors_const
                                // <== sans
                               )
// The GPU global function
{
//   __shared__ float  offspring_genotype[ACTUAL_GENOTYPE_LENGTH]; // sans
  __shared__ int    parent_candidates [4];
  __shared__ float  candidate_energies[4];
  __shared__ int    parents           [2];
  __shared__ int    covr_point        [2];
  __shared__ float  randnums          [10];
  // __shared__ float  sBestEnergy       [32]; // sans
  // __shared__ int    sBestID           [32]; // sans
//   __shared__ float3 calc_coords       [MAX_NUM_OF_ATOMS]; // sans
  // sans ==>
  int actual_genotype_length = cData.dockpars.num_of_genes;
  int genotype_length_in_globmem = actual_genotype_length + 1;
  unsigned int offsetN;
  unsigned int offsetC;
  float3 *calc_coords;
  float  *offspring_genotype;
  __shared__ int reserved;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * NUM_OF_THREADS_PER_BLOCK + actual_genotype_length));
    offspring_genotype = (float *) (calc_coords + NUM_OF_THREADS_PER_BLOCK);;
  }
  else
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * cData.dockpars.num_of_atoms + actual_genotype_length));
    offspring_genotype = (float *) (calc_coords + cData.dockpars.num_of_atoms);
  }
  if(threadIdx.x == 0)
  {
    if(pMem_energies_current[blockIdx.x] < pMem_energies_current[0] + cData.dockpars.addonEng)
    {
      if(fabsf(pMem_energies_current[blockIdx.x] - pMem_energies_next[pMem_energyToConforIndex[blockIdx.x]]) > cData.dockpars.convergedThre)
      {
        if(blockIdx.x < cData.dockpars.maxReservedIdvNum)
        {
          reserved = 1;
        }
      }
      else if(blockIdx.x < cData.dockpars.topN4rescore)
      {
        reserved = 1;
      }
    }
    else
    {
      reserved = 0;
    }
  }
  __syncthreads();
  // <== sans
  __shared__ float  sFloatAccumulator;
//   int run_id; // sans
  int temp_covr_point;
  float energy;
  // int bestID;

  // In this case this compute-unit is responsible for elitist selection
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel4 run_id/e/g/pMem_conformations_current[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_current[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
  // <== sans check
// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
//   for(int i = 0; i < cData.dockpars.pop_size; i++)
//   {
//     printf("sans kernel4 run_win_id blockIdx.x energy1 %d %d %f\n", run_win_id, i, pMem_energies_current[i]);
//   }
// }
//   if ((blockIdx.x % cData.dockpars.pop_size) == 0) { // sans
  if(blockIdx.x == 0) { // sans
    /* sans
    // Find and copy best member of population to position 0
    if (threadIdx.x < cData.dockpars.pop_size)
    {
      // bestID = blockIdx.x + threadIdx.x; // sans
      // energy = pMem_energies_current[blockIdx.x + threadIdx.x]; // sans
      // sans ==>
      bestID = threadIdx.x;
      energy = pMem_energies_current[threadIdx.x];
      // <== sans
    }
    else
    {
      bestID = -1;
      energy = FLT_MAX;
    }

    // Scan through population (we already picked up a blockDim's worth above so skip)
//     for (int i = blockIdx.x + blockDim.x + threadIdx.x; i < blockIdx.x + cData.dockpars.pop_size; i += blockDim.x) // sans
    for (int i = blockDim.x + threadIdx.x; i < cData.dockpars.pop_size; i += blockDim.x) // sans
    {
      float e = pMem_energies_current[i];
      if (e < energy)
      {
        bestID = i;
        energy = e;
      }
    }

    // Reduce to shared memory by warp
    int tgx = threadIdx.x & cData.warpmask;
    WARPMINIMUM2(tgx, energy, bestID);
    int warpID = threadIdx.x >> cData.warpbits;
    if (tgx == 0)
    {
      sBestID[warpID] = bestID;
      sBestEnergy[warpID] = fminf(MAXENERGY, energy);
    }
    __threadfence();
    __syncthreads();

    // Perform final reduction in warp 0
    if (warpID == 0)
    {
      int blocks = blockDim.x / 32;
      if (tgx < blocks)
      {
        bestID = sBestID[tgx];
        energy = sBestEnergy[tgx];
      }
      else
      {
        bestID = -1;
        energy = FLT_MAX;
      }
      WARPMINIMUM2(tgx, energy, bestID);

      if (tgx == 0)
      {
        pMem_energies_next[blockIdx.x] = energy; // sans
        cData.pMem_evals_of_new_entities[blockIdx.x] = 0; // sans
        // sans ==>
        pMem_energies_next[0] = energy;
        pMem_evals_of_new_entities[0] = 0;
        // <== sans
        sBestID[0] = bestID;
      }
    }
    __threadfence();
    __syncthreads();
    sans */
    // sans ==>
    if(threadIdx.x == 0)
    {
      // pMem_energies_next[0] = pMem_energies_current[0];
      pMem_evals_of_new_entities[0] = 0;
    }
    // <== sans

    // Copy best genome to next generation
//     int dOffset = blockIdx.x * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
//     int sOffset = sBestID[0] * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
    // int sOffset = sBestID[0] * genotype_length_in_globmem; // sans
    offsetC = pMem_energyToConforIndex[0] * genotype_length_in_globmem;
    for (int i = threadIdx.x ; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
//       pMem_conformations_next[dOffset + i] = pMem_conformations_current[sOffset + i]; // sans
      // pMem_conformations_next[i] = pMem_conformations_current[sOffset + i]; // sans
      pMem_conformations_next[i] = pMem_conformations_current[offsetC + i]; // sans
//       printf("sans kernel4-3 run_id/i/pMem_conformations_next[i]: %d/%d/%f\n", run_id, i, pMem_conformations_next[i]); // sans check
    }
  }
  // sans ==>
  else if(reserved == 1)
  {
    if(threadIdx.x == 0)
    {
      cData.pMem_reservedIdvNum[run_id]++;
      // pMem_energies_next[blockIdx.x] = pMem_energies_current[blockIdx.x];
      pMem_evals_of_new_entities[blockIdx.x] = 0;
    }
    offsetN = blockIdx.x * genotype_length_in_globmem;
    offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem;
    for(int i = threadIdx.x; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
      pMem_conformations_next[offsetN + i] = pMem_conformations_current[offsetC + i];
    }
  }
  // <== sans
  else
  {
    // sans ==>
    uint32_t evalNum = 0;
    float ligandEnergy = 65536.0f;
    float ligandEnergyTarget = 65500.0f;
    while(ligandEnergy >= ligandEnergyTarget)
    {
    // <== sans
      // Generating the following random numbers:
      // [0..3] for parent candidates,
      // [4..5] for binary tournaments, [6] for deciding crossover,
      // [7..8] for crossover points, [9] for local search
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < 10;
                    gene_counter += blockDim.x)
      {
        // sans check ==>
//         uint old;
//         if(threadIdx.x == 0)
//         {
//           old = cData.pMem_prng_states[(gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x];
//         }
        // <= sans check
//         randnums[gene_counter] = gpu_randf(cData.pMem_prng_states); // sans check
        randnums[gene_counter] = gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x); // sans check
        // sans check ==>
//         if(threadIdx.x == 0)
//         {
// //           printf("sans check kernel4 threadIdx.x/gene_counter/i/cData.pMem_prng_states/old/blockIdx.x/blockDim.x/gridDim.x: %d/%u/%u/%u/%u/%u/%u/%u\n", threadIdx.x, gene_counter, blockIdx.x * blockDim.x + threadIdx.x, cData.pMem_prng_states[blockIdx.x * blockDim.x + threadIdx.x], old, blockIdx.x, blockDim.x, gridDim.x);
//           printf("sans check kernel4 threadIdx.x/gene_counter/i/cData.pMem_prng_states/old/blockIdx.x/blockDim.x/gridDim.x: %d/%u/%u/%u/%u/%u/%u/%u\n", threadIdx.x, gene_counter, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x, cData.pMem_prng_states[(gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x], old, blockIdx.x + gridDim.x * run_id, blockDim.x, gridDim.x * run_id); // sans check
//         }
        // <= sans check
      }
#if 0
      if ((threadIdx.x == 0) && (blockIdx.x == 1))
      {
        printf("%06d ", blockIdx.x);
        for (int i = 0; i < 10; i++)
          printf("%12.6f ", randnums[i]);
        printf("\n");
      }
#endif
      // Determining run ID
//       run_id = blockIdx.x / cData.dockpars.pop_size; // sans useless
      __threadfence();
      __syncthreads();

      if (threadIdx.x < 4) //it is not ensured that the four candidates will be different...
      {
        /* sans
        parent_candidates[threadIdx.x]  = (int) (cData.dockpars.pop_size*randnums[threadIdx.x]); //using randnums[0..3]
        candidate_energies[threadIdx.x] = pMem_energies_current[run_id*cData.dockpars.pop_size+parent_candidates[threadIdx.x]];
        sans */
        // sans ==>
        parent_candidates[threadIdx.x]  = pMem_energyToConforIndex[(int) (cData.dockpars.alivePopNum * randnums[threadIdx.x])]; //using randnums[0..3]
        candidate_energies[threadIdx.x] = pMem_energies_current[parent_candidates[threadIdx.x]];
        // <== sans
      }
      __threadfence();
      __syncthreads();

      if (threadIdx.x < 2)
      {
        // Notice: dockpars_tournament_rate was scaled down to [0,1] in host
        // to reduce number of operations in device
        if (candidate_energies[2*threadIdx.x] < candidate_energies[2*threadIdx.x+1])
        {
          if (/*100.0f**/randnums[4+threadIdx.x] < cData.dockpars.tournament_rate) { // using randnum[4..5]
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x];
          }
          else {
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x+1];
          }
        }
        else
        {
          if (/*100.0f**/randnums[4+threadIdx.x] < cData.dockpars.tournament_rate) {
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x+1];
          }
          else {
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x];
          }
        }
      }
      __threadfence();
      __syncthreads();

      if(blockIdx.x < cData.dockpars.alivePopNum) // sans
      { // sans
        // Performing crossover
        // Notice: dockpars_crossover_rate was scaled down to [0,1] in host
        // to reduce number of operations in device
        if (/*100.0f**/randnums[6] < cData.dockpars.crossover_rate) // Using randnums[6]
        {
          if (threadIdx.x < 2) {
            // Using randnum[7..8]
            covr_point[threadIdx.x] = (int) ((cData.dockpars.num_of_genes-1)*randnums[7+threadIdx.x]);
          }
          __threadfence();
          __syncthreads();

          // covr_point[0] should store the lower crossover-point
          if (threadIdx.x == 0) {
            if (covr_point[1] < covr_point[0]) {
              temp_covr_point = covr_point[1];
              covr_point[1]   = covr_point[0];
              covr_point[0]   = temp_covr_point;
            }
          }

          __threadfence();
          __syncthreads();

          unsigned int offset0 = parents[0] * genotype_length_in_globmem; // sans
          unsigned int offset1 = parents[1] * genotype_length_in_globmem; // sans
          for (uint32_t gene_counter = threadIdx.x;
                        gene_counter < cData.dockpars.num_of_genes;
                        gene_counter+= blockDim.x)
          {
            // Two-point crossover
            if (covr_point[0] != covr_point[1])
            {
              if ((gene_counter <= covr_point[0]) || (gene_counter > covr_point[1]))
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter]; // sans
              else
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[1])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter]; // sans
            }
            // Single-point crossover
            else
            {
              if (gene_counter <= covr_point[0])
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter]; // sans
              else
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[1])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter]; // sans
            }
          }
        }
        else //no crossover
        {
          offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem; // sans
          for (uint32_t gene_counter = threadIdx.x;
                        gene_counter < cData.dockpars.num_of_genes;
                        gene_counter+= blockDim.x)
          {
//             offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM + gene_counter]; // sans
            // // offspring_genotype[gene_counter] = pMem_conformations_current[parents[0] * genotype_length_in_globmem + gene_counter]; // sans
            offspring_genotype[gene_counter] = pMem_conformations_current[offsetC + gene_counter]; // sans
          }
        } // End of crossover

        __threadfence();
        __syncthreads();

        // Performing mutation
        unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x; // sans
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Notice: dockpars_mutation_rate was scaled down to [0,1] in host
          // to reduce number of operations in device
//           if (/*100.0f**/gpu_randf(cData.pMem_prng_states) < cData.dockpars.mutation_rate) // sans check
          if (/*100.0f**/gpu_randf(cData.pMem_prng_states, overallID) < cData.dockpars.mutation_rate) // sans check
          {
            // Translation genes
            if (gene_counter < 3) {
//               offspring_genotype[gene_counter] += cData.dockpars.abs_max_dmov*(2.0f*gpu_randf(cData.pMem_prng_states)-1.0f); // sans check
              offspring_genotype[gene_counter] += cData.dockpars.abs_max_dmov*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f); // sans check
            }
            // Orientation and torsion genes
            else {
//               offspring_genotype[gene_counter] += cData.dockpars.abs_max_dang*(2.0f*gpu_randf(cData.pMem_prng_states)-1.0f); // sans check
              offspring_genotype[gene_counter] += cData.dockpars.abs_max_dang*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f); // sans check
              map_angle(offspring_genotype[gene_counter]);
            }
          }
        } // End of mutation
      // sans ==>
      }
      else
      {
        // regenerate entities at the same position
        unsigned int offset0 = parents[0] * genotype_length_in_globmem; // sans
        unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x; // sans
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          if(gene_counter < 6)
          {
            offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter];
          }
          else
          {
            offspring_genotype[gene_counter] = gpu_randf(cData.pMem_prng_states, overallID) * 360;
          }
        }
      }
      // <== sans

      // Calculating energy of new offspring
      __threadfence();
      __syncthreads();
      // gpu_calc_energy( // sans
      gpu_calc_all_and_ligand_energy(
                      offspring_genotype,
                      energy,
                      ligandEnergy, // sans
//                       run_id, // sans useless
                      calc_coords,
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
//                       rotbonds_unit_vectors_const
                      // <== sans
                     );
      // sans ==>
      evalNum++;
      ligandEnergyTarget *= 1 + 0.01 * evalNum;
      // if (threadIdx.x == 0)
      // {
      //   printf("sans blockIdx.x energy ligandEnergy///ligandEnergyTarget evalNum: %d %f %f///%f %d\n", blockIdx.x, energy, ligandEnergy, ligandEnergyTarget, evalNum);
      // }
      // <== sans
    }
    if (threadIdx.x == 0) {
      pMem_energies_next[blockIdx.x] = energy;
//       cData.pMem_evals_of_new_entities[run_id][blockIdx.x] = 1; // sans
      // pMem_evals_of_new_entities[blockIdx.x] = 1; // sans
      pMem_evals_of_new_entities[blockIdx.x] = evalNum; // sans

      #if defined (DEBUG_ENERGY_KERNEL4)
      printf("%-18s [%-5s]---{%-5s}   [%-10.8f]---{%-10.8f}\n", "-ENERGY-KERNEL4-", "GRIDS", "INTRA", interE, intraE);
      #endif
    }


    // Copying new offspring to next generation
    offsetN = blockIdx.x * genotype_length_in_globmem; // sans
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
//       pMem_conformations_next[blockIdx.x * GENOTYPE_LENGTH_IN_GLOBMEM + gene_counter] = offspring_genotype[gene_counter]; // sans
      pMem_conformations_next[offsetN + gene_counter] = offspring_genotype[gene_counter]; // sans
// if(isnan(pMem_conformations_next[blockIdx.x * genotype_length_in_globmem + gene_counter]))
// {
//   printf("sans check13 NAN kernel4 %d %d %d %d %f %f\n", blockIdx.x * genotype_length_in_globmem + gene_counter, blockIdx.x, genotype_length_in_globmem, gene_counter, pMem_conformations_next[blockIdx.x * genotype_length_in_globmem + gene_counter], offspring_genotype[gene_counter]);
// }
    }
  }

// __threadfence();
// __syncthreads();
// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
//   for(int i = 0; i < cData.dockpars.pop_size; i++)
//   {
//     printf("sans kernel4 run_win_id blockIdx.x energy2 %d %d %f cData.dockpars.alivePopNum%d\n", run_win_id, i, pMem_energies_next[i], cData.dockpars.alivePopNum);
//   }
// }
// __threadfence();
// __syncthreads();

  // sans check ==>
//   __threadfence();
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel4 run_id/e/g/pMem_conformations_next[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
//   __threadfence();
//   __syncthreads();
  // <== sans check
#if 0
  if ((threadIdx.x == 0) && (blockIdx.x == 0))
  {
    printf("%06d %16.8f ", blockIdx.x, pMem_energies_next[blockIdx.x]);
    for (int i = 0; i < cData.dockpars.num_of_genes; i++)
      printf("%12.6f ", pMem_conformations_next[GENOTYPE_LENGTH_IN_GLOBMEM*blockIdx.x + i]);
  }
#endif
}


// sans ==>
__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_gen_and_eval_newpops_ligand_kernel(
                                float* pMem_conformations_current,
                                float* pMem_energies_current,
                                float* pMem_conformations_next,
                                float* pMem_energies_next,
                                // sans ==>
                                int   *pMem_energyToConforIndex,
                                int   *pMem_evals_of_new_entities,
                                float *gFloatBuff,
                                size_t elemNum_gFloatBuffBase,
                                int    run_win_id,
                                int    run_id // sans check

//                                 float    *atom_charges_const,
//                                 uint32_t *atom_types_const,
//                                 uint32_t *atom_types_map_const,
//                                 char     *ignore_inter_const,
//
//                                 uint32_t *intraE_contributors_const,

//                                 int      *rotlist_const,
//
//                                 float    *ref_coords_const,
//                                 float    *rotbonds_moving_vectors_const,
//                                 float    *rotbonds_unit_vectors_const
                                // <== sans
                               )
// The GPU global function
{
//   __shared__ float  offspring_genotype[ACTUAL_GENOTYPE_LENGTH]; // sans
  __shared__ int    parent_candidates [4];
  __shared__ float  candidate_energies[4];
  __shared__ int    parents           [2];
  __shared__ int    covr_point        [2];
  __shared__ float  randnums          [10];
  // __shared__ float  sBestEnergy       [32]; // sans
  // __shared__ int    sBestID           [32]; // sans
//   __shared__ float3 calc_coords       [MAX_NUM_OF_ATOMS]; // sans
  // sans ==>
  int actual_genotype_length = cData.dockpars.num_of_genes;
  int genotype_length_in_globmem = actual_genotype_length + 1;
  unsigned int offsetN;
  unsigned int offsetC;
  float3 *calc_coords;
  float  *offspring_genotype;
  __shared__ int reserved;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * NUM_OF_THREADS_PER_BLOCK + actual_genotype_length));
    offspring_genotype = (float *) (calc_coords + NUM_OF_THREADS_PER_BLOCK);;
  }
  else
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * cData.dockpars.num_of_atoms + actual_genotype_length));
    offspring_genotype = (float *) (calc_coords + cData.dockpars.num_of_atoms);
  }
  if(threadIdx.x == 0)
  {
    if(pMem_energies_current[blockIdx.x] < pMem_energies_current[0] + cData.dockpars.addonEng)
    {
      if(fabsf(pMem_energies_current[blockIdx.x] - pMem_energies_next[pMem_energyToConforIndex[blockIdx.x]]) > cData.dockpars.convergedThre)
      {
        if(blockIdx.x < cData.dockpars.maxReservedIdvNum)
        {
          reserved = 1;
        }
      }
      else if(blockIdx.x < cData.dockpars.topN4rescore)
      {
        reserved = 1;
      }
    }
    else
    {
      reserved = 0;
    }
  }
  __syncthreads();
  // <== sans
  __shared__ float  sFloatAccumulator;
//   int run_id; // sans
  int temp_covr_point;
  float energy;
  // int bestID;

  // In this case this compute-unit is responsible for elitist selection
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel4 run_id/e/g/pMem_conformations_current[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_current[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
  // <== sans check
// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
//   for(int i = 0; i < cData.dockpars.pop_size; i++)
//   {
//     printf("sans kernel4 run_win_id blockIdx.x energy1 %d %d %f\n", run_win_id, i, pMem_energies_current[i]);
//   }
// }
//   if ((blockIdx.x % cData.dockpars.pop_size) == 0) { // sans
  if(blockIdx.x == 0) { // sans
    /* sans
    // Find and copy best member of population to position 0
    if (threadIdx.x < cData.dockpars.pop_size)
    {
      // bestID = blockIdx.x + threadIdx.x; // sans
      // energy = pMem_energies_current[blockIdx.x + threadIdx.x]; // sans
      // sans ==>
      bestID = threadIdx.x;
      energy = pMem_energies_current[threadIdx.x];
      // <== sans
    }
    else
    {
      bestID = -1;
      energy = FLT_MAX;
    }

    // Scan through population (we already picked up a blockDim's worth above so skip)
//     for (int i = blockIdx.x + blockDim.x + threadIdx.x; i < blockIdx.x + cData.dockpars.pop_size; i += blockDim.x) // sans
    for (int i = blockDim.x + threadIdx.x; i < cData.dockpars.pop_size; i += blockDim.x) // sans
    {
      float e = pMem_energies_current[i];
      if (e < energy)
      {
        bestID = i;
        energy = e;
      }
    }

    // Reduce to shared memory by warp
    int tgx = threadIdx.x & cData.warpmask;
    WARPMINIMUM2(tgx, energy, bestID);
    int warpID = threadIdx.x >> cData.warpbits;
    if (tgx == 0)
    {
      sBestID[warpID] = bestID;
      sBestEnergy[warpID] = fminf(MAXENERGY, energy);
    }
    __threadfence();
    __syncthreads();

    // Perform final reduction in warp 0
    if (warpID == 0)
    {
      int blocks = blockDim.x / 32;
      if (tgx < blocks)
      {
        bestID = sBestID[tgx];
        energy = sBestEnergy[tgx];
      }
      else
      {
        bestID = -1;
        energy = FLT_MAX;
      }
      WARPMINIMUM2(tgx, energy, bestID);

      if (tgx == 0)
      {
        pMem_energies_next[blockIdx.x] = energy; // sans
        cData.pMem_evals_of_new_entities[blockIdx.x] = 0; // sans
        // sans ==>
        pMem_energies_next[0] = energy;
        pMem_evals_of_new_entities[0] = 0;
        // <== sans
        sBestID[0] = bestID;
      }
    }
    __threadfence();
    __syncthreads();
    sans */
    // sans ==>
    if(threadIdx.x == 0)
    {
      // pMem_energies_next[0] = pMem_energies_current[0];
      pMem_evals_of_new_entities[0] = 0;
    }
    // <== sans

    // Copy best genome to next generation
//     int dOffset = blockIdx.x * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
//     int sOffset = sBestID[0] * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
    // int sOffset = sBestID[0] * genotype_length_in_globmem; // sans
    offsetC = pMem_energyToConforIndex[0] * genotype_length_in_globmem;
    for (int i = threadIdx.x ; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
//       pMem_conformations_next[dOffset + i] = pMem_conformations_current[sOffset + i]; // sans
      // pMem_conformations_next[i] = pMem_conformations_current[sOffset + i]; // sans
      pMem_conformations_next[i] = pMem_conformations_current[offsetC + i]; // sans
//       printf("sans kernel4-3 run_id/i/pMem_conformations_next[i]: %d/%d/%f\n", run_id, i, pMem_conformations_next[i]); // sans check
    }
  }
  // sans ==>
  else if(reserved == 1)
  {
    if(threadIdx.x == 0)
    {
      cData.pMem_reservedIdvNum[run_id]++;
      // pMem_energies_next[blockIdx.x] = pMem_energies_current[blockIdx.x];
      pMem_evals_of_new_entities[blockIdx.x] = 0;
    }
    offsetN = blockIdx.x * genotype_length_in_globmem;
    offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem;
    for(int i = threadIdx.x; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
      pMem_conformations_next[offsetN + i] = pMem_conformations_current[offsetC + i];
    }
  }
  // <== sans
  else
  {
    // sans ==>
    uint32_t evalNum = 0;
    float ligandEnergy = 65536.0f;
    float ligandEnergyTarget = 65500.0f;
    while(ligandEnergy >= ligandEnergyTarget)
    {
    // <== sans
      // Generating the following random numbers:
      // [0..3] for parent candidates,
      // [4..5] for binary tournaments, [6] for deciding crossover,
      // [7..8] for crossover points, [9] for local search
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < 10;
                    gene_counter += blockDim.x)
      {
        // sans check ==>
//         uint old;
//         if(threadIdx.x == 0)
//         {
//           old = cData.pMem_prng_states[(gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x];
//         }
        // <= sans check
//         randnums[gene_counter] = gpu_randf(cData.pMem_prng_states); // sans check
        randnums[gene_counter] = gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x); // sans check
        // sans check ==>
//         if(threadIdx.x == 0)
//         {
// //           printf("sans check kernel4 threadIdx.x/gene_counter/i/cData.pMem_prng_states/old/blockIdx.x/blockDim.x/gridDim.x: %d/%u/%u/%u/%u/%u/%u/%u\n", threadIdx.x, gene_counter, blockIdx.x * blockDim.x + threadIdx.x, cData.pMem_prng_states[blockIdx.x * blockDim.x + threadIdx.x], old, blockIdx.x, blockDim.x, gridDim.x);
//           printf("sans check kernel4 threadIdx.x/gene_counter/i/cData.pMem_prng_states/old/blockIdx.x/blockDim.x/gridDim.x: %d/%u/%u/%u/%u/%u/%u/%u\n", threadIdx.x, gene_counter, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x, cData.pMem_prng_states[(gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x], old, blockIdx.x + gridDim.x * run_id, blockDim.x, gridDim.x * run_id); // sans check
//         }
        // <= sans check
      }
#if 0
      if ((threadIdx.x == 0) && (blockIdx.x == 1))
      {
        printf("%06d ", blockIdx.x);
        for (int i = 0; i < 10; i++)
          printf("%12.6f ", randnums[i]);
        printf("\n");
      }
#endif
      // Determining run ID
//       run_id = blockIdx.x / cData.dockpars.pop_size; // sans useless
      __threadfence();
      __syncthreads();

      if (threadIdx.x < 4) //it is not ensured that the four candidates will be different...
      {
        /* sans
        parent_candidates[threadIdx.x]  = (int) (cData.dockpars.pop_size*randnums[threadIdx.x]); //using randnums[0..3]
        candidate_energies[threadIdx.x] = pMem_energies_current[run_id*cData.dockpars.pop_size+parent_candidates[threadIdx.x]];
        sans */
        // sans ==>
        parent_candidates[threadIdx.x]  = pMem_energyToConforIndex[(int) (cData.dockpars.alivePopNum * randnums[threadIdx.x])]; //using randnums[0..3]
        candidate_energies[threadIdx.x] = pMem_energies_current[parent_candidates[threadIdx.x]];
        // <== sans
      }
      __threadfence();
      __syncthreads();

      if (threadIdx.x < 2)
      {
        // Notice: dockpars_tournament_rate was scaled down to [0,1] in host
        // to reduce number of operations in device
        if (candidate_energies[2*threadIdx.x] < candidate_energies[2*threadIdx.x+1])
        {
          if (/*100.0f**/randnums[4+threadIdx.x] < cData.dockpars.tournament_rate) { // using randnum[4..5]
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x];
          }
          else {
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x+1];
          }
        }
        else
        {
          if (/*100.0f**/randnums[4+threadIdx.x] < cData.dockpars.tournament_rate) {
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x+1];
          }
          else {
            parents[threadIdx.x] = parent_candidates[2*threadIdx.x];
          }
        }
      }
      __threadfence();
      __syncthreads();

      if(blockIdx.x < cData.dockpars.alivePopNum) // sans
      { // sans
        // Performing crossover
        // Notice: dockpars_crossover_rate was scaled down to [0,1] in host
        // to reduce number of operations in device
        if (/*100.0f**/randnums[6] < cData.dockpars.crossover_rate) // Using randnums[6]
        {
          if (threadIdx.x < 2) {
            // Using randnum[7..8]
            covr_point[threadIdx.x] = (int) ((cData.dockpars.num_of_genes-1)*randnums[7+threadIdx.x]);
          }
          __threadfence();
          __syncthreads();

          // covr_point[0] should store the lower crossover-point
          if (threadIdx.x == 0) {
            if (covr_point[1] < covr_point[0]) {
              temp_covr_point = covr_point[1];
              covr_point[1]   = covr_point[0];
              covr_point[0]   = temp_covr_point;
            }
          }

          __threadfence();
          __syncthreads();

          unsigned int offset0 = parents[0] * genotype_length_in_globmem; // sans
          unsigned int offset1 = parents[1] * genotype_length_in_globmem; // sans
          for (uint32_t gene_counter = threadIdx.x;
                        gene_counter < cData.dockpars.num_of_genes;
                        gene_counter+= blockDim.x)
          {
            // Two-point crossover
            if (covr_point[0] != covr_point[1])
            {
              if ((gene_counter <= covr_point[0]) || (gene_counter > covr_point[1]))
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter]; // sans
              else
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[1])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter]; // sans
            }
            // Single-point crossover
            else
            {
              if (gene_counter <= covr_point[0])
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter]; // sans
              else
//                 offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[1])*GENOTYPE_LENGTH_IN_GLOBMEM+gene_counter]; // sans
                offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter]; // sans
            }
          }
        }
        else //no crossover
        {
          offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem; // sans
          for (uint32_t gene_counter = threadIdx.x;
                        gene_counter < cData.dockpars.num_of_genes;
                        gene_counter+= blockDim.x)
          {
//             offspring_genotype[gene_counter] = pMem_conformations_current[(run_id*cData.dockpars.pop_size+parents[0])*GENOTYPE_LENGTH_IN_GLOBMEM + gene_counter]; // sans
            // // offspring_genotype[gene_counter] = pMem_conformations_current[parents[0] * genotype_length_in_globmem + gene_counter]; // sans
            offspring_genotype[gene_counter] = pMem_conformations_current[offsetC + gene_counter]; // sans
          }
        } // End of crossover

        __threadfence();
        __syncthreads();

        // Performing mutation
        unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x; // sans
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Notice: dockpars_mutation_rate was scaled down to [0,1] in host
          // to reduce number of operations in device
//           if (/*100.0f**/gpu_randf(cData.pMem_prng_states) < cData.dockpars.mutation_rate) // sans check
          if (/*100.0f**/gpu_randf(cData.pMem_prng_states, overallID) < cData.dockpars.mutation_rate) // sans check
          {
            // Translation genes
            if (gene_counter < 3) {
//               offspring_genotype[gene_counter] += cData.dockpars.abs_max_dmov*(2.0f*gpu_randf(cData.pMem_prng_states)-1.0f); // sans check
              offspring_genotype[gene_counter] += cData.dockpars.abs_max_dmov*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f); // sans check
            }
            // Orientation and torsion genes
            else {
//               offspring_genotype[gene_counter] += cData.dockpars.abs_max_dang*(2.0f*gpu_randf(cData.pMem_prng_states)-1.0f); // sans check
              offspring_genotype[gene_counter] += cData.dockpars.abs_max_dang*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f); // sans check
              map_angle(offspring_genotype[gene_counter]);
            }
          }
        } // End of mutation
      // sans ==>
      }
      else
      {
        // regenerate entities at the same position
        unsigned int offset0 = parents[0] * genotype_length_in_globmem; // sans
        unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x; // sans
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          if(gene_counter < 6)
          {
            offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter];
          }
          else
          {
            offspring_genotype[gene_counter] = gpu_randf(cData.pMem_prng_states, overallID) * 360;
          }
        }
      }
      // <== sans

      // Calculating energy of new offspring
      __threadfence();
      __syncthreads();
      gpu_calc_ligand_energy(
                      offspring_genotype,
                      energy,
//                       run_id, // sans useless
                      calc_coords,
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
//                       rotbonds_unit_vectors_const
                      // <== sans
                     );
      // sans ==>
      evalNum++;
      ligandEnergyTarget *= 1 + 0.01 * evalNum;
      // <== sans
    }
    if (threadIdx.x == 0) {
      pMem_energies_next[blockIdx.x] = energy;
//       cData.pMem_evals_of_new_entities[run_id][blockIdx.x] = 1; // sans
      // pMem_evals_of_new_entities[blockIdx.x] = 1; // sans
      pMem_evals_of_new_entities[blockIdx.x] = evalNum; // sans

      #if defined (DEBUG_ENERGY_KERNEL4)
      printf("%-18s [%-5s]---{%-5s}   [%-10.8f]---{%-10.8f}\n", "-ENERGY-KERNEL4-", "GRIDS", "INTRA", interE, intraE);
      #endif
    }


    // Copying new offspring to next generation
    offsetN = blockIdx.x * genotype_length_in_globmem; // sans
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
//       pMem_conformations_next[blockIdx.x * GENOTYPE_LENGTH_IN_GLOBMEM + gene_counter] = offspring_genotype[gene_counter]; // sans
      pMem_conformations_next[offsetN + gene_counter] = offspring_genotype[gene_counter]; // sans
// if(isnan(pMem_conformations_next[blockIdx.x * genotype_length_in_globmem + gene_counter]))
// {
//   printf("sans check13 NAN kernel4 %d %d %d %d %f %f\n", blockIdx.x * genotype_length_in_globmem + gene_counter, blockIdx.x, genotype_length_in_globmem, gene_counter, pMem_conformations_next[blockIdx.x * genotype_length_in_globmem + gene_counter], offspring_genotype[gene_counter]);
// }
    }
  }

// __threadfence();
// __syncthreads();
// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
//   for(int i = 0; i < cData.dockpars.pop_size; i++)
//   {
//     printf("sans kernel4 run_win_id blockIdx.x energy2 %d %d %f cData.dockpars.alivePopNum%d\n", run_win_id, i, pMem_energies_next[i], cData.dockpars.alivePopNum);
//   }
// }
// __threadfence();
// __syncthreads();

  // sans check ==>
//   __threadfence();
//   __syncthreads();
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     int pop_size=50; // sans check
//     for(int e = 0; e < pop_size; e++) // sans check
//     {
//       for(int g = 0; g < genotype_length_in_globmem; g++)
//       {
//         printf("sans kernel4 run_id/e/g/pMem_conformations_next[], %d/%d/%d/%f\n", run_id, e, g, pMem_conformations_next[e*genotype_length_in_globmem+g]);
//       }
//     }
//   }
//   __threadfence();
//   __syncthreads();
  // <== sans check
#if 0
  if ((threadIdx.x == 0) && (blockIdx.x == 0))
  {
    printf("%06d %16.8f ", blockIdx.x, pMem_energies_next[blockIdx.x]);
    for (int i = 0; i < cData.dockpars.num_of_genes; i++)
      printf("%12.6f ", pMem_conformations_next[GENOTYPE_LENGTH_IN_GLOBMEM*blockIdx.x + i]);
  }
#endif
}
// <== sans
__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
checkcheckKernel4(
                                float* pMem_conformations_current,
                                float* pMem_energies_current,
                                float* pMem_conformations_next,
                                float* pMem_energies_next,
                                // sans ==>
                                int   *pMem_energyToConforIndex,
                                int   *pMem_evals_of_new_entities,
                                float *gFloatBuff,
                                size_t elemNum_gFloatBuffBase,
                                int    run_win_id,
                                int    run_id
                               )
// The GPU global function
{
if(threadIdx.x == 0 && blockIdx.x == 0)
{
  for(int i = 0; i < cData.dockpars.pop_size; i++)
  {
    printf("sans kernel4 run_win_id blockIdx.x energy2 %d %d %f cData.dockpars.alivePopNum%d\n", run_win_id, i, pMem_energies_next[i], cData.dockpars.alivePopNum);
  }
}
}


void gpu_gen_and_eval_newpops(
                              uint32_t blocks,
                              uint32_t threadsPerBlock,
                              /* sans
                              float*   pMem_conformations_current,
                              float*   pMem_energies_current,
                              float*   pMem_conformations_next,
                              float*   pMem_energies_next,
                              sans */
                              // sans ==>
                              bool     ligandOnly,
                              float  **pMem_conformations_current,
                              float  **pMem_energies_current,
                              float  **pMem_conformations_next,
                              float  **pMem_energies_next,
                              int    **pMem_energyToConforIndex,
                              int    **pMem_evals_of_new_entities,
                              float   *gFloatBuff,
                              size_t   elemNum_gFloatBuffBase,
                              int      run_id_init,
                              int      run_id_end
//                               GpuData& cData
                              // <== sans
                             )
{
//   gpu_gen_and_eval_newpops_kernel<<<blocks, threadsPerBlock>>>(pMem_conformations_current, pMem_energies_current, pMem_conformations_next, pMem_energies_next); // sans
  // sans ==>
  int run_win_id = 0;
  if(ligandOnly == false)
  {
    for(int i = run_id_init; i < run_id_end; i++)
    {
//       printf("sans gpu_gen_and_eval_newpops blocks/threadsPerBlock: %zu/%zu/%d\n", blocks, threadsPerBlock, i);fflush(stdout); // sans check
      gpu_gen_and_eval_newpops_kernel<<<blocks, threadsPerBlock>>>
      (
        pMem_conformations_current[i],
        pMem_energies_current[i],
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
//       printf("sans gpu_gen_and_eval_newpops blocks/threadsPerBlock: %zu/%zu/%d\n", blocks, threadsPerBlock, i);fflush(stdout); // sans check
      gpu_gen_and_eval_newpops_ligand_kernel<<<blocks, threadsPerBlock>>>
      (
        pMem_conformations_current[i],
        pMem_energies_current[i],
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
  LAUNCHERROR("gpu_gen_and_eval_newpops_kernel");
#if 0
  cudaError_t status;
  status = cudaDeviceSynchronize();
  RTERROR(status, "gpu_gen_and_eval_newpops_kernel");
  status = cudaDeviceReset();
  RTERROR(status, "failed to shut down");
  exit(0);
#endif
}
