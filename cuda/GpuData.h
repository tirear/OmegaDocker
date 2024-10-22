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


#ifndef GPUDATADOTH
#define GPUDATADOTH
#include <float.h>


static const int   TERMBITS         = 10;
static const float MAXTERM          = (float)(1 << (31 - TERMBITS - 8)); // 2^(31 - 10 - 8) = 2^13 = 8192
static const float TERMSCALE        = (float)(1 << TERMBITS);            // 2^10 = 1024
static const float ONEOVERTERMSCALE = 1.0f / TERMSCALE;                  // 1 / 1024 = 0.000977
static const float MAXREDUCE        = (float)(1 << (31 - TERMBITS - 4)); // 2^(31 - 10 - 4) = 2^17 = 131072

static const float MAXENERGY        = FLT_MAX / 100.0f; // Used to cap absurd energies so placeholder energy is always skipped in sorts
static const float MAXFORCE         = FLT_MAX / 100.0f; // Used to cap absurd gradients

#define RTERROR(status, s) \
  if (status != cudaSuccess) { \
    printf("%s %s\n", s, cudaGetErrorString(status)); \
    assert(0); \
    cudaDeviceReset(); \
    exit(-1); \
  }


#define SYNCHRONOUS
#ifdef SYNCHRONOUS
#define LAUNCHERROR(s) \
  { \
    cudaError_t status = cudaGetLastError(); \
    if (status != cudaSuccess) { \
      printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
      cudaDeviceReset(); \
      exit(-1); \
    } \
    status = cudaDeviceSynchronize(); \
    RTERROR(status, s); \
  }
#else
#define LAUNCHERROR(s) \
  { \
    cudaError_t status = cudaGetLastError(); \
    if (status != cudaSuccess) { \
      printf("Error: %s launching kernel %s\n", cudaGetErrorString(status), s); \
      cudaDeviceReset(); \
      exit(-1); \
    } \
  }
#endif

typedef struct
{
  int             num_of_atoms;
  int             true_ligand_atoms;
  int             num_of_atypes;
  int             num_of_map_atypes;
  int             num_of_intraE_contributors;
  int             gridsize_x;
  int             gridsize_y;
  int             gridsize_z;
  int             gridsize_x_times_y;
  int             gridsize_x_times_y_times_z;
  float           grid_spacing;
  int             rotbondlist_length;
  float           coeff_elec;
  float           elec_min_distance;
  float           coeff_desolv;
  int             pop_size;
  int             num_of_genes;
  float           tournament_rate;
  float           crossover_rate;
  float           mutation_rate;
  float           abs_max_dmov;
  float           abs_max_dang;
  float           lsearch_rate;
  float           smooth;
  unsigned int    num_of_lsentities;
  float           rho_lower_bound;
  float           base_dmov_mul_sqrt3;
  float           base_dang_mul_sqrt3;
  unsigned int    cons_limit;
  unsigned int    max_num_of_iters;
  float           qasp;
  float           adam_beta1;
  float           adam_beta2;
  float           adam_epsilon;
  unsigned long num_of_runs;
  bool  miniCoorOut;
  float negaLigaStepSize;
  float negaCompStepSize;
  int   ligaRhoLimit;
  int   compRhoLimit;
  int   collRhoLimit;
  float collDistance;
  int   alivePopNum;
  bool  miniOnly;
  float addonEng;
  int   maxReservedIdvNum;
  float convergedThre;
  unsigned int topN4rescore;
} GpuDockparameters;

struct GpuData {
  int                             devnum;
  int                             devid;
  int                             preallocated_gridsize;
  GpuDockparameters               dockpars;

  // Consolidated constants and memory pointers to reduce kernel launch overhead
  float    *pKerconst_interintra_atom_charges_const;
  uint32_t *pKerconst_interintra_atom_types_const;
  uint32_t *pKerconst_interintra_atom_types_map_const;
  char     *pKerconst_interintra_ignore_inter_const;

  uint32_t *pKerconst_intracontrib_intraE_contributors_const;

  unsigned int       *pKerconst_intra_atom_types_reqm_const;
  unsigned short int *pKerconst_intra_VWpars_exp_const;
  float              *pKerconst_intra_reqm_AB_const;
  float              *pKerconst_intra_VWpars_AC_const;
  float              *pKerconst_intra_VWpars_BD_const;
  float              *pKerconst_intra_dspars_S_const;
  float              *pKerconst_intra_dspars_V_const;
  int  *pKerconst_rotlist_rotlist_const;
  int  *pKerconst_rotlist_rotlist_atomId_const;
  int  *pKerconst_rotlist_rotlist_rotatableBondId_const;

  float *pKerconst_conform_ref_coords_const;
  float *pKerconst_conform_rotbonds_moving_vectors_const;
  float *pKerconst_conform_rotbonds_unit_vectors_const;

  int *pKerconst_grads_rotbonds;
  int *pKerconst_grads_rotbonds_atoms;
  int *pKerconst_grads_num_rotating_atoms_per_rotbond;
  float*                          pMem_fgrids;
  int                           **pMem_evals_of_new_entities; // [run_num][ndv_num] individual (ndv)
  int                            *pMem_reservedIdvNum;        // [run_num]
  int*                            pMem_gpu_evals_of_runs;
  uint32_t*                       pMem_prng_states;
  float*                          pMem_angle_const;
  float*                          pMem_dependence_on_theta_const;
  float*                          pMem_dependence_on_rotangle_const;

  // CUDA-specific constants
  unsigned int                    warpmask;
  unsigned int                    warpbits;
};

struct GpuTempData {
  float*      pMem_fgrids;
  float     **pMem_conformations1;
  float     **pMem_conformations2;
  float     **pMem_energies1;
  float     **pMem_energies2;
  int       **pMem_energyToConforIndex;   // [run_num][ndv_num] individual (ndv)
  int       **pMem_evals_of_new_entities; // [run_num][ndv_num] individual (ndv)
  int        *pMem_reservedIdvNum;        // [run_num]
  int*        pMem_gpu_evals_of_runs;
  uint32_t*   pMem_prng_states;
  char*       device_name;
  bool        device_busy;
  float*   gFloatBuff;
};
#endif
