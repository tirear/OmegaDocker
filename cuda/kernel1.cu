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


// sans check ==>
// __global__ void sansTest0(int * aa)
// {
//   if(threadIdx.x == 0 && blockIdx.x == 0)
//   {
// //     for(uint i = 0; i < cData.dockpars.rotbondlist_length; i++)
//     for(uint i = 0; i < 2432; i++)
//       printf("sans sansTest0 i/aa[i]: %u/%d\n", i, aa[i]);
//   }
// }
//
// void sansGpuTest0(uint32_t a, uint32_t b, int * aa)
// {
//   sansTest0<<<a, b>>>(aa);
// }
//
// __global__ void sansTest1()
// {
//   if(threadIdx.x == 0 && blockIdx.x == 0)
//   {
//     for(uint i = 0; i < cData.dockpars.rotbondlist_length; i++)
//       printf("sans sansTest1 i/cData.pKerconst_rotlist_rotlist_const[i]: %u/%u/%d\n", i, cData.dockpars.rotbondlist_length, cData.pKerconst_rotlist_rotlist_const[i]);
//   }
// }
//
// void sansGpuTest1(uint32_t a, uint32_t b)
// {
//   sansTest1<<<a, b>>>();
// }
// <= sans check

__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_calc_initpop_kernel(
                        float* pMem_conformations_current,
                        float* pMem_energies_current,
                        // sans ==>
                        int   *pMem_energyToConforIndex,
                        int   *pMem_evals_of_new_entities,
                        float *gFloatBuff,
                        size_t elemNum_gFloatBuffBase,
                        int    run_win_id

//                         float    *atom_charges_const,
//                         uint32_t *atom_types_const,
//                         uint32_t *atom_types_map_const,
//                         char     *ignore_inter_const,
//
//                         uint32_t *intraE_contributors_const,

//                         int      *rotlist_const,
//
//                         float    *ref_coords_const,
//                         float    *rotbonds_moving_vectors_const,
//                         float    *rotbonds_unit_vectors_const
                        // <== sans
                       )
{
//   __shared__ float3 calc_coords[MAX_NUM_OF_ATOMS]; // sans
  __shared__ float sFloatAccumulator;
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-1\n");
//   }
  // <== sans check
  // sans ==>
  float3 *calc_coords;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * NUM_OF_THREADS_PER_BLOCK);
  }
  else
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * cData.dockpars.num_of_atoms);
  }
  // <== sans
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-2\n");
//   }
  // <== sans check
  float  energy = 0.0f;
//   int    run_id = blockIdx.x / cData.dockpars.pop_size; // sans useless
//   float* pGenotype = pMem_conformations_current + blockIdx.x * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
  int genotype_length_in_globmem = cData.dockpars.num_of_genes + 1;      // sans
  float* pGenotype = pMem_conformations_current + blockIdx.x * genotype_length_in_globmem; // sans

  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-3\n");
//   }
  // <== sans check
  // =============================================================
  gpu_calc_energy(
                  pGenotype,
                  energy,
//                   run_id, // sans useless
                  calc_coords,
                  &sFloatAccumulator
                  // sans ==>
//                   atom_charges_const,
//                   atom_types_const,
//                   atom_types_map_const,
//                   ignore_inter_const,
//
//                   intraE_contributors_const,

//                   rotlist_const,
//
//                   ref_coords_const,
//                   rotbonds_moving_vectors_const,
//                   rotbonds_unit_vectors_const
                  // <== sans
                 );
  // =============================================================
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-4\n");
//   }
  // <== sans check
//     calc_coords = NULL; // sans check
//     pGenotype   = NULL; // sans check

  // Write out final energy
  if (threadIdx.x == 0)
  {
    pMem_energies_current[blockIdx.x] = energy;
//     cData.pMem_evals_of_new_entities[run_id][blockIdx.x] = 1; // sans
    pMem_evals_of_new_entities[blockIdx.x] = 1; // sans [run_num][ndv_num] individual (ndv)
  }

  // sans ==>
  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
  // <== sans

  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-5\n");
//   }
  // <== sans check
}


// sans ==>
__global__ void
__launch_bounds__(NUM_OF_THREADS_PER_BLOCK, 1024 / NUM_OF_THREADS_PER_BLOCK)
gpu_calc_ligand_initpop_kernel(
                        float* pMem_conformations_current,
                        float* pMem_energies_current,
                        // sans ==>
                        int   *pMem_energyToConforIndex,
                        int   *pMem_evals_of_new_entities,
                        float *gFloatBuff,
                        size_t elemNum_gFloatBuffBase,
                        int    run_win_id

//                         float    *atom_charges_const,
//                         uint32_t *atom_types_const,
//                         uint32_t *atom_types_map_const,
//                         char     *ignore_inter_const,
//
//                         uint32_t *intraE_contributors_const,

//                         int      *rotlist_const,
//
//                         float    *ref_coords_const,
//                         float    *rotbonds_moving_vectors_const,
//                         float    *rotbonds_unit_vectors_const
                        // <== sans
                       )
{
//   __shared__ float3 calc_coords[MAX_NUM_OF_ATOMS]; // sans
  __shared__ float sFloatAccumulator;
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-1\n");
//   }
  // <== sans check
  // sans ==>
  float3 *calc_coords;
  if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * NUM_OF_THREADS_PER_BLOCK);
  }
  else
  {
    calc_coords = (float3 *) (gFloatBuff + run_win_id * elemNum_gFloatBuffBase + blockIdx.x * 3 * cData.dockpars.num_of_atoms);
  }
  // <== sans
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-2\n");
//   }
  // <== sans check
  float  energy = 0.0f;
//   int    run_id = blockIdx.x / cData.dockpars.pop_size; // sans useless
//   float* pGenotype = pMem_conformations_current + blockIdx.x * GENOTYPE_LENGTH_IN_GLOBMEM; // sans
  int genotype_length_in_globmem = cData.dockpars.num_of_genes + 1;      // sans
  float* pGenotype = pMem_conformations_current + blockIdx.x * genotype_length_in_globmem; // sans

  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-3\n");
//   }
  // <== sans check
  // =============================================================
  gpu_calc_ligand_energy(
                  pGenotype,
                  energy,
//                   run_id, // sans useless
                  calc_coords,
                  &sFloatAccumulator
                  // sans ==>
//                   atom_charges_const,
//                   atom_types_const,
//                   atom_types_map_const,
//                   ignore_inter_const,
//
//                   intraE_contributors_const,

//                   rotlist_const,
//
//                   ref_coords_const,
//                   rotbonds_moving_vectors_const,
//                   rotbonds_unit_vectors_const
                  // <== sans
                 );
  // =============================================================
  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-4\n");
//   }
  // <== sans check
//     calc_coords = NULL; // sans check
//     pGenotype   = NULL; // sans check

  // Write out final energy
  if (threadIdx.x == 0)
  {
    pMem_energies_current[blockIdx.x] = energy;
//     cData.pMem_evals_of_new_entities[run_id][blockIdx.x] = 1; // sans
    pMem_evals_of_new_entities[blockIdx.x] = 1; // sans [run_num][ndv_num] individual (ndv)
  }

  // sans ==>
  if(threadIdx.x == 0)
  {
    pMem_energyToConforIndex[blockIdx.x] = blockIdx.x;
  }
  // <== sans

  // sans check ==>
//   if(blockIdx.x == 0 && threadIdx.x == 0)
//   {
//     printf("sans kernel1-5\n");
//   }
  // <== sans check
}
// <== sans


void gpu_calc_initpop(
                      uint32_t blocks,
                      uint32_t threadsPerBlock,
                      /* sans
                      float*   pConformations_current,
                      float*   pEnergies_current,
                      sans */
                      // sans ==>
                      bool     ligandOnly,
                      float  **pConformations_current,
                      float  **pEnergies_current,
                      int    **pMem_energyToConforIndex,
                      int    **pMem_evals_of_new_entities,
                      float   *gFloatBuff,
                      size_t   elemNum_gFloatBuffBase,
                      int      run_id_init,
                      int      run_id_end
//                       GpuData& cData
                      // <== sans
                     )
{
//   gpu_calc_initpop_kernel<<<blocks, threadsPerBlock>>>(pConformations_current, pEnergies_current); // sans
  // sans ==>
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
      // run_win_id++; // sans force run_window = 1 from v.0.17
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
      // run_win_id++; // sans force run_window = 1 from v.0.17
    }
  }
  // <== sans
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
