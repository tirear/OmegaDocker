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
gpu_sum_evals_kernel(
                     int *pMem_evals_of_new_entities,
                     int  run_id
                    )
// The GPU global function sums the evaluation counter states
// which are stored in evals_of_new_entities array foreach entity,
// calculates the sums for each run and stores it in evals_of_runs array.
// The number of blocks which should be started equals to num_of_runs,
// since each block performs the summation for one run.
{
  __shared__ int sSum_evals;
  int partsum_evals = 0;
  int* pEvals_of_new_entities = pMem_evals_of_new_entities;
  for (int entity_counter = threadIdx.x;
           entity_counter < cData.dockpars.pop_size;
           entity_counter += blockDim.x) 
  {
    partsum_evals += pEvals_of_new_entities[entity_counter];
  }
//   partsum_evals = 0;
  
  // Perform warp-wise reduction
  REDUCEINTEGERSUM(partsum_evals, &sSum_evals);
  if (threadIdx.x == 0)
  {
      cData.pMem_gpu_evals_of_runs[run_id] += sSum_evals;
  }
}

void gpu_sum_evals(uint32_t blocks, uint32_t threadsPerBlock,
                   int    **pMem_evals_of_new_entities,
                   int      run_id)
{
  gpu_sum_evals_kernel<<<blocks, threadsPerBlock>>>(pMem_evals_of_new_entities[run_id], run_id);
  LAUNCHERROR("gpu_sum_evals_kernel");
#if 0
  cudaError_t status;
  status = cudaDeviceSynchronize();
  RTERROR(status, "gpu_sum_evals_kernel");
  status = cudaDeviceReset();
  RTERROR(status, "failed to shut down");
  exit(0);
#endif
}
