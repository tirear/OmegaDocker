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

__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_calc_initpop_kernel(
                        float* pMem_conformations_current,
                        float* pMem_energies_current,
                        int   *pMem_energyToConforIndex,
                        int   *pMem_evals_of_new_entities,
                        float *gFloatBuff,
                        size_t elemNum_gFloatBuffBase,
                        int    run_win_id
                       )
{
  __shared__ float sFloatAccumulator;
  float3 *calc_coords;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * NUM_OF_THREADS_PER_BLOCK);
  }
  else
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * cData.dockpars.num_of_atoms);
  }
  float  energy = 0.0f;
  int genotype_length_in_globmem = cData.dockpars.num_of_genes + 1;
  float* pGenotype = pMem_conformations_current + blockIdx.x * genotype_length_in_globmem;

  // =============================================================
  gpu_calc_energy(
                  pGenotype,
                  energy,
                  calc_coords,
                  &sFloatAccumulator
                 );
  // =============================================================

  // Write out final energy
  if (threadIdx.x == 0)
  {
    pMem_energies_current[blockIdx.x] = energy;
    pMem_evals_of_new_entities[blockIdx.x] = 1; // [run_num][ndv_num] individual (ndv)
  }

  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
}


__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_calc_ligand_initpop_kernel(
                        float* pMem_conformations_current,
                        float* pMem_energies_current,
                        int   *pMem_energyToConforIndex,
                        int   *pMem_evals_of_new_entities,
                        float *gFloatBuff,
                        size_t elemNum_gFloatBuffBase,
                        int    run_win_id
                       )
{
  __shared__ float sFloatAccumulator;
  float3 *calc_coords;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * NUM_OF_THREADS_PER_BLOCK);
  }
  else
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * cData.dockpars.num_of_atoms);
  }
  float  energy = 0.0f;
  int genotype_length_in_globmem = cData.dockpars.num_of_genes + 1;
  float* pGenotype = pMem_conformations_current + blockIdx.x * genotype_length_in_globmem;

  // =============================================================
  gpu_calc_ligand_energy(
                  pGenotype,
                  energy,
                  calc_coords,
                  &sFloatAccumulator
                 );
  // =============================================================

  // Write out final energy
  if (threadIdx.x == 0)
  {
    pMem_energies_current[blockIdx.x] = energy;
    pMem_evals_of_new_entities[blockIdx.x] = 1; // [run_num][ndv_num] individual (ndv)
  }

  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
}


void gpu_calc_initpop(
                      uint32_t blocks,
                      uint32_t threadsPerBlock,
                      bool     ligandOnly,
                      float  **pConformations_current,
                      float  **pEnergies_current,
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
      gpu_calc_initpop_kernel<<<blocks, threadsPerBlock>>>
      (
        pConformations_current[i],
        pEnergies_current[i],
        pMem_energyToConforIndex[i],
        pMem_evals_of_new_entities[i],
        gFloatBuff,
        elemNum_gFloatBuffBase,
        run_win_id
      );
    }
  }
  else
  {
    for(int i = run_id_init; i < run_id_end; i++)
    {
      gpu_calc_ligand_initpop_kernel<<<blocks, threadsPerBlock>>>
      (
        pConformations_current[i],
        pEnergies_current[i],
        pMem_energyToConforIndex[i],
        pMem_evals_of_new_entities[i],
        gFloatBuff,
        elemNum_gFloatBuffBase,
        run_win_id
      );
    }
  }
  LAUNCHERROR("gpu_calc_initpop_kernel");
#if 0
  cudaError_t status;
  status = cudaDeviceSynchronize();
  RTERROR(status, "gpu_calc_initpop_kernel");
  status = cudaDeviceReset();
  RTERROR(status, "failed to shut down");
  exit(0);
#endif
}
