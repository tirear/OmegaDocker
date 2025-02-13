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


// if defined, new (experimental) SW genotype moves that are dependent
// on nr of atoms and nr of torsions of ligand are used
#define SWAT3 // Third set of Solis-Wets hyperparameters by Andreas Tillack

__global__ void
#if (__CUDA_ARCH__ == 750)
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
#else
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1408 / NUM_OF_THREADS_PER_BLOCK)
#endif
gpu_perform_LS_kernel(
                      float* pMem_conformations_next,
                      float* pMem_energies_next,
                      int   *pMem_energyToConforIndex,
                      int   *pMem_evals_of_new_entities,
                      float *gFloatBuff,
                      size_t elemNum_gFloatBuffBase,
                      int    run_win_id,
                      int    run_id
                     )
// The GPU global function performs local search on the pre-defined entities of conformations_next.
// The number of blocks which should be started equals to num_of_lsentities*num_of_runs.
// This way the first num_of_lsentities entity of each population will be subjected to local search
// (and each block carries out the algorithm for one entity).
// Since the first entity is always the best one in the current population,
// it is always tested according to the ls probability, and if it not to be
// subjected to local search, the entity with ID num_of_lsentities is selected instead of the first one (with ID 0).
{
  __shared__ float rho;
  __shared__ int   cons_succ;
  __shared__ int   cons_fail;
  __shared__ int   iteration_cnt;
  __shared__ int   evaluation_cnt;

  __shared__ float offspring_energy;
  __shared__ float sFloatAccumulator;
  __shared__ int entity_id;


  int actual_genotype_length = cData.dockpars.num_of_genes;
  int genotype_length_in_globmem = actual_genotype_length + 1;
  float candidate_energy;
  float3 *calc_coords;
  float  *genotype_candidate;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * NUM_OF_THREADS_PER_BLOCK + 4 * actual_genotype_length));
    genotype_candidate = (float*)(calc_coords + NUM_OF_THREADS_PER_BLOCK);
  }
  else
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * cData.dockpars.num_of_atoms + 4 * actual_genotype_length));
    genotype_candidate = (float*)(calc_coords + cData.dockpars.num_of_atoms);
  }

  // Genotype pointers
  float* genotype_deviate = (float*)(genotype_candidate + actual_genotype_length);
  float* genotype_bias = (float*)(genotype_deviate + actual_genotype_length);
  float* offspring_genotype = (float*)(genotype_bias + actual_genotype_length);

  // Determining run ID and entity ID
  // Initializing offspring genotype
  if (threadIdx.x == 0)
  {
    entity_id = blockIdx.x % cData.dockpars.num_of_lsentities;
    // Since entity 0 is the best one due to elitism,
    // it should be subjected to random selection
    if (entity_id == 0) {
      // If entity 0 is not selected according to LS-rate,
      // choosing an other entity
      if (100.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x) > cData.dockpars.lsearch_rate) {
        entity_id = cData.dockpars.num_of_lsentities;
      }
    }

    offspring_energy = pMem_energies_next[entity_id];
    rho = 1.0f;
    cons_succ = 0;
    cons_fail = 0;
    iteration_cnt = 0;
    evaluation_cnt = 0;
  }
  __syncthreads();

  size_t offset = entity_id * genotype_length_in_globmem;
  for (uint32_t gene_counter = threadIdx.x;
                gene_counter < cData.dockpars.num_of_genes;
                gene_counter+= blockDim.x)
  {
    offspring_genotype[gene_counter] = pMem_conformations_next[offset + gene_counter];
    genotype_bias[gene_counter] = 0.0f;
  }
  __syncthreads();


#ifdef SWAT3
  float lig_scale = 1.0f/sqrt((float)cData.dockpars.num_of_atoms);
  float gene_scale = 1.0f/sqrt((float)cData.dockpars.num_of_genes);
#endif
  while ((iteration_cnt < cData.dockpars.max_num_of_iters) && (rho > cData.dockpars.rho_lower_bound))
  {
    // New random deviate
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
#ifdef SWAT3
      genotype_deviate[gene_counter] = rho*(2.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x)-1.0f)*(gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x) < gene_scale);

      // Translation genes
      if (gene_counter < 3) {
        genotype_deviate[gene_counter] *= cData.dockpars.base_dmov_mul_sqrt3;
      }
      // Orientation and torsion genes
      else {
        if (gene_counter < 6) {
          genotype_deviate[gene_counter] *= cData.dockpars.base_dang_mul_sqrt3 * lig_scale;
        } else {
          genotype_deviate[gene_counter] *= cData.dockpars.base_dang_mul_sqrt3 * gene_scale;
        }
      }
#else
      genotype_deviate[gene_counter] = rho*(2.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x)-1.0f)*(gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x)<0.3f);

      // Translation genes
      if (gene_counter < 3) {
        genotype_deviate[gene_counter] *= cData.dockpars.base_dmov_mul_sqrt3;
      }
      // Orientation and torsion genes
      else {
        genotype_deviate[gene_counter] *= cData.dockpars.base_dang_mul_sqrt3;
      }
#endif
    }

    // Generating new genotype candidate
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
      genotype_candidate[gene_counter] = offspring_genotype[gene_counter] +
                                         genotype_deviate[gene_counter]   +
                                         genotype_bias[gene_counter];
    }
    // Evaluating candidate
    __syncthreads();

    // =================================================================
    gpu_calc_energy(
                    genotype_candidate,
                    candidate_energy,
                    calc_coords,
                    &sFloatAccumulator
                   );
    // =================================================================
    if (threadIdx.x == 0) {
      evaluation_cnt++;
    }
    __syncthreads();

    if (candidate_energy < offspring_energy)  // If candidate is better, success
    {
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < cData.dockpars.num_of_genes;
                    gene_counter+= blockDim.x)
      {
        // Updating offspring_genotype
        offspring_genotype[gene_counter] = genotype_candidate[gene_counter];
        // Updating genotype_bias
        genotype_bias[gene_counter] = 0.6f*genotype_bias[gene_counter] + 0.4f*genotype_deviate[gene_counter];
      }

      // Work-item 0 will overwrite the shared variables
      // used in the previous if condition
      __syncthreads();

      if (threadIdx.x == 0)
      {
        offspring_energy = candidate_energy;
        cons_succ++;
        cons_fail = 0;
      }
    }
    else // If candidate is worse, check the opposite direction
    {
      // Generating the other genotype candidate
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < cData.dockpars.num_of_genes;
                    gene_counter+= blockDim.x)
      {
        genotype_candidate[gene_counter] = offspring_genotype[gene_counter] -
                                           genotype_deviate[gene_counter] -
                                           genotype_bias[gene_counter];
      }

      // Evaluating candidate
      __syncthreads();

      // =================================================================
      gpu_calc_energy(
                      genotype_candidate,
                      candidate_energy,
                      calc_coords,
                      &sFloatAccumulator
                     );
      // =================================================================

      if (threadIdx.x == 0) {
        evaluation_cnt++;

        #if defined (DEBUG_ENERGY_KERNEL)
        printf("%-18s [%-5s]---{%-5s}   [%-10.8f]---{%-10.8f}\n", "-ENERGY-KERNEL3-", "GRIDS", "INTRA", partial_interE[0], partial_intraE[0]);
        #endif
      }
      __syncthreads();

      if (candidate_energy < offspring_energy) // If candidate is better, success
      {
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Updating offspring_genotype
          offspring_genotype[gene_counter] = genotype_candidate[gene_counter];
          // Updating genotype_bias
          genotype_bias[gene_counter] = 0.6f*genotype_bias[gene_counter] - 0.4f*genotype_deviate[gene_counter];
        }

        // Work-item 0 will overwrite the shared variables
        // used in the previous if condition
              __syncthreads();

        if (threadIdx.x == 0)
        {
          offspring_energy = candidate_energy;
          cons_succ++;
          cons_fail = 0;
        }
      }
      else  // Failure in both directions
      {
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Updating genotype_bias
          genotype_bias[gene_counter] = 0.5f*genotype_bias[gene_counter];
        }
        if (threadIdx.x == 0)
        {
          cons_succ = 0;
          cons_fail++;
        }
      }
    }

    // Changing rho if needed
    if (threadIdx.x == 0)
    {
      iteration_cnt++;
      if (cons_succ >= cData.dockpars.cons_limit)
      {
        rho *= LS_EXP_FACTOR;
        cons_succ = 0;
      }
      else
        if (cons_fail >= cData.dockpars.cons_limit)
        {
          rho *= LS_CONT_FACTOR;
          cons_fail = 0;
        }
    }
    __syncthreads();
  }

  // Updating eval counter and energy
  if (threadIdx.x == 0) {
    pMem_evals_of_new_entities[entity_id] += evaluation_cnt;
    pMem_energies_next[entity_id] = offspring_energy;
  }

  // Mapping torsion angles and writing out results
  offset = entity_id * genotype_length_in_globmem;
  for (uint32_t gene_counter = threadIdx.x;
                gene_counter < cData.dockpars.num_of_genes;
                gene_counter+= blockDim.x)
  {
    if (gene_counter >= 3) {
      map_angle(offspring_genotype[gene_counter]);
    }
    pMem_conformations_next[offset + gene_counter] = offspring_genotype[gene_counter];
  }

  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
}


__global__ void
#if (__CUDA_ARCH__ == 750)
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
#else
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1408 / NUM_OF_THREADS_PER_BLOCK)
#endif
gpu_perform_ligand_LS_kernel(
                      float* pMem_conformations_next,
                      float* pMem_energies_next,
                      int   *pMem_energyToConforIndex,
                      int   *pMem_evals_of_new_entities,
                      float *gFloatBuff,
                      size_t elemNum_gFloatBuffBase,
                      int    run_win_id,
                      int    run_id
                     )
// The GPU global function performs local search on the pre-defined entities of conformations_next.
// The number of blocks which should be started equals to num_of_lsentities*num_of_runs.
// This way the first num_of_lsentities entity of each population will be subjected to local search
// (and each block carries out the algorithm for one entity).
// Since the first entity is always the best one in the current population,
// it is always tested according to the ls probability, and if it not to be
// subjected to local search, the entity with ID num_of_lsentities is selected instead of the first one (with ID 0).
{
  __shared__ float rho;
  __shared__ int   cons_succ;
  __shared__ int   cons_fail;
  __shared__ int   iteration_cnt;
  __shared__ int   evaluation_cnt;

  __shared__ float offspring_energy;
  __shared__ float sFloatAccumulator;
  __shared__ int entity_id;


  int actual_genotype_length = cData.dockpars.num_of_genes;
  int genotype_length_in_globmem = actual_genotype_length + 1;
  float candidate_energy;
  float3 *calc_coords;
  float  *genotype_candidate;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * NUM_OF_THREADS_PER_BLOCK + 4 * actual_genotype_length));
    genotype_candidate = (float*)(calc_coords + NUM_OF_THREADS_PER_BLOCK);
  }
  else
  {
    calc_coords = (float3*) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * (3 * cData.dockpars.num_of_atoms + 4 * actual_genotype_length));
    genotype_candidate = (float*)(calc_coords + cData.dockpars.num_of_atoms);
  }

  // Genotype pointers
  float* genotype_deviate = (float*)(genotype_candidate + actual_genotype_length);
  float* genotype_bias = (float*)(genotype_deviate + actual_genotype_length);
  float* offspring_genotype = (float*)(genotype_bias + actual_genotype_length);

  // Determining run ID and entity ID
  // Initializing offspring genotype
  if (threadIdx.x == 0)
  {
    entity_id = blockIdx.x % cData.dockpars.num_of_lsentities;
    // Since entity 0 is the best one due to elitism,
    // it should be subjected to random selection
    if (entity_id == 0) {
      // If entity 0 is not selected according to LS-rate,
      // choosing an other entity
      if (100.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x) > cData.dockpars.lsearch_rate) {
        entity_id = cData.dockpars.num_of_lsentities;
      }
    }

    offspring_energy = pMem_energies_next[entity_id];
    rho = 1.0f;
    cons_succ = 0;
    cons_fail = 0;
    iteration_cnt = 0;
    evaluation_cnt = 0;
  }
  __syncthreads();

  size_t offset = entity_id * genotype_length_in_globmem;
  for (uint32_t gene_counter = threadIdx.x;
                gene_counter < cData.dockpars.num_of_genes;
                gene_counter+= blockDim.x)
  {
    offspring_genotype[gene_counter] = pMem_conformations_next[offset + gene_counter];
    genotype_bias[gene_counter] = 0.0f;
  }
  __syncthreads();


#ifdef SWAT3
  float lig_scale = 1.0f/sqrt((float)cData.dockpars.num_of_atoms);
  float gene_scale = 1.0f/sqrt((float)cData.dockpars.num_of_genes);
#endif
  while ((iteration_cnt < cData.dockpars.max_num_of_iters) && (rho > cData.dockpars.rho_lower_bound))
  {
    // New random deviate
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
#ifdef SWAT3
      genotype_deviate[gene_counter] = rho*(2.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x)-1.0f)*(gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x) < gene_scale);

      // Translation genes
      if (gene_counter < 3) {
        genotype_deviate[gene_counter] *= cData.dockpars.base_dmov_mul_sqrt3;
      }
      // Orientation and torsion genes
      else {
        if (gene_counter < 6) {
          genotype_deviate[gene_counter] *= cData.dockpars.base_dang_mul_sqrt3 * lig_scale;
        } else {
          genotype_deviate[gene_counter] *= cData.dockpars.base_dang_mul_sqrt3 * gene_scale;
        }
      }
#else
      genotype_deviate[gene_counter] = rho*(2.0f*gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x)-1.0f)*(gpu_randf(cData.pMem_prng_states, (gridDim.x * run_id + blockIdx.x) * blockDim.x + threadIdx.x)<0.3f);

      // Translation genes
      if (gene_counter < 3) {
        genotype_deviate[gene_counter] *= cData.dockpars.base_dmov_mul_sqrt3;
      }
      // Orientation and torsion genes
      else {
        genotype_deviate[gene_counter] *= cData.dockpars.base_dang_mul_sqrt3;
      }
#endif
    }

    // Generating new genotype candidate
    for (uint32_t gene_counter = threadIdx.x;
                  gene_counter < cData.dockpars.num_of_genes;
                  gene_counter+= blockDim.x)
    {
      genotype_candidate[gene_counter] = offspring_genotype[gene_counter] +
                                         genotype_deviate[gene_counter]   +
                                         genotype_bias[gene_counter];
    }
    // Evaluating candidate
    __syncthreads();

    // =================================================================
    gpu_calc_ligand_energy(
                    genotype_candidate,
                    candidate_energy,
                    calc_coords,
                    &sFloatAccumulator
                   );
    // =================================================================
    if (threadIdx.x == 0) {
      evaluation_cnt++;
    }
    __syncthreads();

    if (candidate_energy < offspring_energy)  // If candidate is better, success
    {
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < cData.dockpars.num_of_genes;
                    gene_counter+= blockDim.x)
      {
        // Updating offspring_genotype
        offspring_genotype[gene_counter] = genotype_candidate[gene_counter];
        // Updating genotype_bias
        genotype_bias[gene_counter] = 0.6f*genotype_bias[gene_counter] + 0.4f*genotype_deviate[gene_counter];
      }

      // Work-item 0 will overwrite the shared variables
      // used in the previous if condition
      __syncthreads();

      if (threadIdx.x == 0)
      {
        offspring_energy = candidate_energy;
        cons_succ++;
        cons_fail = 0;
      }
    }
    else // If candidate is worse, check the opposite direction
    {
      // Generating the other genotype candidate
      for (uint32_t gene_counter = threadIdx.x;
                    gene_counter < cData.dockpars.num_of_genes;
                    gene_counter+= blockDim.x)
      {
        genotype_candidate[gene_counter] = offspring_genotype[gene_counter] -
                                           genotype_deviate[gene_counter] -
                                           genotype_bias[gene_counter];
      }

      // Evaluating candidate
      __syncthreads();

      // =================================================================
      gpu_calc_ligand_energy(
                      genotype_candidate,
                      candidate_energy,
                      calc_coords,
                      &sFloatAccumulator
                     );
      // =================================================================

      if (threadIdx.x == 0) {
        evaluation_cnt++;

        #if defined (DEBUG_ENERGY_KERNEL)
        printf("%-18s [%-5s]---{%-5s}   [%-10.8f]---{%-10.8f}\n", "-ENERGY-KERNEL3-", "GRIDS", "INTRA", partial_interE[0], partial_intraE[0]);
        #endif
      }
      __syncthreads();

      if (candidate_energy < offspring_energy) // If candidate is better, success
      {
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Updating offspring_genotype
          offspring_genotype[gene_counter] = genotype_candidate[gene_counter];
          // Updating genotype_bias
          genotype_bias[gene_counter] = 0.6f*genotype_bias[gene_counter] - 0.4f*genotype_deviate[gene_counter];
        }

        // Work-item 0 will overwrite the shared variables
        // used in the previous if condition
              __syncthreads();

        if (threadIdx.x == 0)
        {
          offspring_energy = candidate_energy;
          cons_succ++;
          cons_fail = 0;
        }
      }
      else  // Failure in both directions
      {
        for (uint32_t gene_counter = threadIdx.x;
                      gene_counter < cData.dockpars.num_of_genes;
                      gene_counter+= blockDim.x)
        {
          // Updating genotype_bias
          genotype_bias[gene_counter] = 0.5f*genotype_bias[gene_counter];
        }
        if (threadIdx.x == 0)
        {
          cons_succ = 0;
          cons_fail++;
        }
      }
    }

    // Changing rho if needed
    if (threadIdx.x == 0)
    {
      iteration_cnt++;
      if (cons_succ >= cData.dockpars.cons_limit)
      {
        rho *= LS_EXP_FACTOR;
        cons_succ = 0;
      }
      else
        if (cons_fail >= cData.dockpars.cons_limit)
        {
          rho *= LS_CONT_FACTOR;
          cons_fail = 0;
        }
    }
    __syncthreads();
  }

  // Updating eval counter and energy
  if (threadIdx.x == 0) {
    pMem_evals_of_new_entities[entity_id] += evaluation_cnt;
    pMem_energies_next[entity_id] = offspring_energy;
  }

  // Mapping torsion angles and writing out results
  offset = entity_id * genotype_length_in_globmem;
  for (uint32_t gene_counter = threadIdx.x;
                gene_counter < cData.dockpars.num_of_genes;
                gene_counter+= blockDim.x)
  {
    if (gene_counter >= 3) {
      map_angle(offspring_genotype[gene_counter]);
    }
    pMem_conformations_next[offset + gene_counter] = offspring_genotype[gene_counter];
  }

  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
}


void gpu_perform_LS(
                    uint32_t blocks,
                    uint32_t threads,
                    bool     ligandOnly,
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
      gpu_perform_LS_kernel<<<blocks, threads>>>
      (
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
      gpu_perform_ligand_LS_kernel<<<blocks, threads>>>
      (
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
  LAUNCHERROR("gpu_perform_LS_kernel");
#if 0
  cudaError_t status;
  status = cudaDeviceSynchronize();
  RTERROR(status, "gpu_perform_LS_kernel");
  status = cudaDeviceReset();
  RTERROR(status, "failed to shut down");
  exit(0);
#endif
}
