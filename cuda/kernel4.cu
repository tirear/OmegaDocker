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


//#define DEBUG_ENERGY_KERNEL4

__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_gen_and_eval_newpops_kernel(
                                float* pMem_conformations_current,
                                float* pMem_energies_current,
                                float* pMem_conformations_next,
                                float* pMem_energies_next,
                                int   *pMem_energyToConforIndex,
                                int   *pMem_evals_of_new_entities,
                                float *gFloatBuff,
                                size_t elemNum_gFloatBuffBase,
                                int    run_win_id,
                                int    run_id
                               )
// The GPU global function
{
//   __shared__ float  offspring_genotype[ACTUAL_GENOTYPE_LENGTH];
  __shared__ int    parent_candidates [4];
  __shared__ float  candidate_energies[4];
  __shared__ int    parents           [2];
  __shared__ int    covr_point        [2];
  __shared__ float  randnums          [10];
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
  __shared__ float  sFloatAccumulator;
  int temp_covr_point;
  float energy;

  if(blockIdx.x == 0) {
    if(threadIdx.x == 0)
    {
      pMem_evals_of_new_entities[0] = 0;
    }

    // Copy best genome to next generation
    offsetC = pMem_energyToConforIndex[0] * genotype_length_in_globmem;
    for (int i = threadIdx.x ; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
      pMem_conformations_next[i] = pMem_conformations_current[offsetC + i];
    }
  }
  else if(reserved == 1)
  {
    if(threadIdx.x == 0)
    {
      cData.pMem_reservedIdvNum[run_id]++;
      pMem_evals_of_new_entities[blockIdx.x] = 0;
    }
    offsetN = blockIdx.x * genotype_length_in_globmem;
    offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem;
    for(int i = threadIdx.x; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
      pMem_conformations_next[offsetN + i] = pMem_conformations_current[offsetC + i];
    }
  }
  else
  {
    // Generating the following random numbers:
    // [0..3] for parent candidates,
    // [4..5] for binary tournaments, [6] for deciding crossover,
    // [7..8] for crossover points, [9] for local search
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < 10;
                  gene_counter += blockDim.x)
    {
      randnums[gene_counter] = gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x);
    }
    // Determining run ID
    __threadfence();
    __syncthreads();

    if (threadIdx.x < 4) //it is not ensured that the four candidates will be different...
    {
      parent_candidates[threadIdx.x]  = pMem_energyToConforIndex[(int) (cData.dockpars.alivePopNum * randnums[threadIdx.x])]; //using randnums[0..3]
      candidate_energies[threadIdx.x] = pMem_energies_current[parent_candidates[threadIdx.x]];
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

    if(blockIdx.x < cData.dockpars.alivePopNum)
    {
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

        unsigned int offset0 = parents[0] * genotype_length_in_globmem;
        unsigned int offset1 = parents[1] * genotype_length_in_globmem;
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Two-point crossover
          if (covr_point[0] != covr_point[1])
          {
            if ((gene_counter <= covr_point[0]) || (gene_counter > covr_point[1]))
              offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter];
            else
              offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter];
          }
          // Single-point crossover
          else
          {
            if (gene_counter <= covr_point[0])
              offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter];
            else
              offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter];
          }
        }
      }
      else //no crossover
      {
        offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem;
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          offspring_genotype[gene_counter] = pMem_conformations_current[offsetC + gene_counter];
        }
      } // End of crossover

      __threadfence();
      __syncthreads();

      // Performing mutation
      unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x;
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < cData.dockpars.num_of_genes;
                    gene_counter+= blockDim.x)
      {
        // Notice: dockpars_mutation_rate was scaled down to [0,1] in host
        // to reduce number of operations in device
        if (/*100.0f**/gpu_randf(cData.pMem_prng_states, overallID) < cData.dockpars.mutation_rate)
        {
          // Translation genes
          if (gene_counter < 3) {
            offspring_genotype[gene_counter] += cData.dockpars.abs_max_dmov*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f); 
          }
          // Orientation and torsion genes
          else {
            offspring_genotype[gene_counter] += cData.dockpars.abs_max_dang*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f); 
            map_angle(offspring_genotype[gene_counter]);
          }
        }
      } // End of mutation
    }
    else
    {
      // regenerate entities at the same position
      unsigned int offset0 = parents[0] * genotype_length_in_globmem;
      unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x;
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

    // Calculating energy of new offspring
    __threadfence();
    __syncthreads();
    gpu_calc_energy(
                    offspring_genotype,
                    energy,
                    calc_coords,
                    &sFloatAccumulator
                   );
    if (threadIdx.x == 0) {
      pMem_energies_next[blockIdx.x] = energy;
      pMem_evals_of_new_entities[blockIdx.x] = 1;

      #if defined (DEBUG_ENERGY_KERNEL4)
      printf("%-18s [%-5s]---{%-5s}   [%-10.8f]---{%-10.8f}\n", "-ENERGY-KERNEL4-", "GRIDS", "INTRA", interE, intraE);
      #endif
    }


    // Copying new offspring to next generation
    offsetN = blockIdx.x * genotype_length_in_globmem;
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
      pMem_conformations_next[offsetN + gene_counter] = offspring_genotype[gene_counter];
    }
  }
}


__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_gen_and_eval_newpops_ligand_kernel(
                                float* pMem_conformations_current,
                                float* pMem_energies_current,
                                float* pMem_conformations_next,
                                float* pMem_energies_next,
                                int   *pMem_energyToConforIndex,
                                int   *pMem_evals_of_new_entities,
                                float *gFloatBuff,
                                size_t elemNum_gFloatBuffBase,
                                int    run_win_id,
                                int    run_id
                               )
// The GPU global function
{
  __shared__ int    parent_candidates [4];
  __shared__ float  candidate_energies[4];
  __shared__ int    parents           [2];
  __shared__ int    covr_point        [2];
  __shared__ float  randnums          [10];
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
  __shared__ float  sFloatAccumulator;
  int temp_covr_point;
  float energy;

  if(blockIdx.x == 0) {
    if(threadIdx.x == 0)
    {
      pMem_evals_of_new_entities[0] = 0;
    }

    // Copy best genome to next generation
    offsetC = pMem_energyToConforIndex[0] * genotype_length_in_globmem;
    for (int i = threadIdx.x ; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
      pMem_conformations_next[i] = pMem_conformations_current[offsetC + i];
    }
  }
  else if(reserved == 1)
  {
    if(threadIdx.x == 0)
    {
      cData.pMem_reservedIdvNum[run_id]++;
      pMem_evals_of_new_entities[blockIdx.x] = 0;
    }
    offsetN = blockIdx.x * genotype_length_in_globmem;
    offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem;
    for(int i = threadIdx.x; i < cData.dockpars.num_of_genes; i += blockDim.x)
    {
      pMem_conformations_next[offsetN + i] = pMem_conformations_current[offsetC + i];
    }
  }
  else
  {
    // Generating the following random numbers:
    // [0..3] for parent candidates,
    // [4..5] for binary tournaments, [6] for deciding crossover,
    // [7..8] for crossover points, [9] for local search
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < 10;
                  gene_counter += blockDim.x)
    {
      randnums[gene_counter] = gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x);
    }
    // Determining run ID
    __threadfence();
    __syncthreads();

    if (threadIdx.x < 4) //it is not ensured that the four candidates will be different...
    {
      parent_candidates[threadIdx.x]  = pMem_energyToConforIndex[(int) (cData.dockpars.alivePopNum * randnums[threadIdx.x])]; //using randnums[0..3]
      candidate_energies[threadIdx.x] = pMem_energies_current[parent_candidates[threadIdx.x]];
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

    if(blockIdx.x < cData.dockpars.alivePopNum)
    {
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

        unsigned int offset0 = parents[0] * genotype_length_in_globmem;
        unsigned int offset1 = parents[1] * genotype_length_in_globmem;
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Two-point crossover
          if (covr_point[0] != covr_point[1])
          {
            if ((gene_counter <= covr_point[0]) || (gene_counter > covr_point[1]))
              offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter];
            else
              offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter];
          }
          // Single-point crossover
          else
          {
            if (gene_counter <= covr_point[0])
              offspring_genotype[gene_counter] = pMem_conformations_current[offset0 + gene_counter];
            else
              offspring_genotype[gene_counter] = pMem_conformations_current[offset1 + gene_counter];
          }
        }
      }
      else //no crossover
      {
        offsetC = pMem_energyToConforIndex[blockIdx.x] * genotype_length_in_globmem;
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          offspring_genotype[gene_counter] = pMem_conformations_current[offsetC + gene_counter];
        }
      } // End of crossover

      __threadfence();
      __syncthreads();

      // Performing mutation
      unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x;
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < cData.dockpars.num_of_genes;
                    gene_counter+= blockDim.x)
      {
        // Notice: dockpars_mutation_rate was scaled down to [0,1] in host
        // to reduce number of operations in device
        if (/*100.0f**/gpu_randf(cData.pMem_prng_states, overallID) < cData.dockpars.mutation_rate)
        {
          // Translation genes
          if (gene_counter < 3) {
            offspring_genotype[gene_counter] += cData.dockpars.abs_max_dmov*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f);
          }
          // Orientation and torsion genes
          else {
            offspring_genotype[gene_counter] += cData.dockpars.abs_max_dang*(2.0f*gpu_randf(cData.pMem_prng_states, overallID)-1.0f); 
            map_angle(offspring_genotype[gene_counter]);
          }
        }
      } // End of mutation
    }
    else
    {
      // regenerate entities at the same position
      unsigned int offset0 = parents[0] * genotype_length_in_globmem;
      unsigned int overallID = (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x;
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

    // Calculating energy of new offspring
    __threadfence();
    __syncthreads();
    gpu_calc_ligand_energy(
                    offspring_genotype,
                    energy,
                    calc_coords,
                    &sFloatAccumulator
                   );
    if (threadIdx.x == 0) {
      pMem_energies_next[blockIdx.x] = energy;
      pMem_evals_of_new_entities[blockIdx.x] = 1;

      #if defined (DEBUG_ENERGY_KERNEL4)
      printf("%-18s [%-5s]---{%-5s}   [%-10.8f]---{%-10.8f}\n", "-ENERGY-KERNEL4-", "GRIDS", "INTRA", interE, intraE);
      #endif
    }


    // Copying new offspring to next generation
    offsetN = blockIdx.x * genotype_length_in_globmem;
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
      pMem_conformations_next[offsetN + gene_counter] = offspring_genotype[gene_counter];
    }
  }
}


__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
checkcheckKernel4(
                                float* pMem_conformations_current,
                                float* pMem_energies_current,
                                float* pMem_conformations_next,
                                float* pMem_energies_next,
                                int   *pMem_energyToConforIndex,
                                int   *pMem_evals_of_new_entities,
                                float *gFloatBuff,
                                size_t elemNum_gFloatBuffBase,
                                int    run_win_id,
                                int    run_id
                               )
// The GPU global function
{}


void gpu_gen_and_eval_newpops(
                              uint32_t blocks,
                              uint32_t threadsPerBlock,
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
                             )
{
  int run_win_id = 0;
  if(ligandOnly == false)
  {
    for(int i = run_id_init; i < run_id_end; i++)
    {
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
