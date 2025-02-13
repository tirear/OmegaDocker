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


#ifndef CALCENERGY_H_
#define CALCENERGY_H_

#include <math.h>
#include <stdio.h>
#include <cstdint>

#include "calcenergy_basic.h"
#include "miscellaneous.h"
#include "processligand.h"
#include "getparameters.h"

// This struct is passed to the GPU global functions (OpenCL kernels) as input.
// Its members are parameters related to the ligand, the grid
// and the genetic algorithm, or they are pointers of GPU (ADM FPGA) memory areas
// used for storing different data such as the current
// and the next population genotypes and energies, the grids,
// the evaluation counters and the random number generator states.
typedef struct
{
  int           num_of_atoms;
  int           true_ligand_atoms;
  int           num_of_atypes;
  int           num_of_map_atypes;
  int           num_of_intraE_contributors;
  int           gridsize_x;
  int           gridsize_y;
  int           gridsize_z;
  float         grid_spacing;
  float*        fgrids;
  int           rotbondlist_length;
  float         coeff_elec;
  float         elec_min_distance;
  float         coeff_desolv;
  float*        conformations_current;
  float*        energies_current;
  float*        conformations_next;
  float*        energies_next;
  int*          evals_of_new_entities;
  unsigned int* prng_states;
  int           pop_size;
  int           num_of_genes;
  float         tournament_rate;
  float         crossover_rate;
  float         mutation_rate;
  float         abs_max_dmov;
  float         abs_max_dang;
  float         lsearch_rate;
  float         smooth;
  unsigned int  num_of_lsentities;
  float         rho_lower_bound;
  float         base_dmov_mul_sqrt3;
  float         base_dang_mul_sqrt3;
  unsigned int  cons_limit;
  unsigned int  max_num_of_iters;
  float         qasp;
  float         adam_beta1;
  float         adam_beta2;
  float         adam_epsilon;
} Dockparameters;

typedef struct
{
  float    *atom_charges_const;
  uint32_t *atom_types_const;
  uint32_t *atom_types_map_const;
  char     *ignore_inter_const;
} kernelconstant_interintra;

typedef struct
{
  uint32_t *intraE_contributors_const;
} kernelconstant_intracontrib;

typedef struct
{
  unsigned int       *atom_types_reqm_const;
  unsigned short int *VWpars_exp_const;
  float              *reqm_AB_const;
  float              *VWpars_AC_const;
  float              *VWpars_BD_const;
  float              *dspars_S_const;
  float              *dspars_V_const;
} kernelconstant_intra;

typedef struct
{
  int *rotlist_const;
  int  *rotlist_atomId;
  int  *rotlist_rotatableBondId;
} kernelconstant_rotlist;

typedef struct
{
  float *ref_coords_const;
  float *rotbonds_moving_vectors_const;
  float *rotbonds_unit_vectors_const;
} kernelconstant_conform;

typedef struct
{
  // Added for calculating torsion-related gradients.
  // Passing list of rotbond-atoms ids to the GPU.
  // Contains the same information as processligand.h/Liganddata->rotbonds 

  // Contains the same information as processligand.h/Liganddata->atom_rotbonds
  // "atom_rotbonds": array that contains the rotatable bonds - atoms assignment.
  // If the element atom_rotbonds[atom index][rotatable bond index] is equal to 1,
  // it means,that the atom must be rotated if the bond rotates. A 0 means the opposite.
  int *rotbonds;
  int *rotbonds_atoms;
  int *num_rotating_atoms_per_rotbond;
} kernelconstant_grads;

int prepare_const_fields_for_gpu(
                                 Liganddata*                  myligand_reference,
                                 Dockpars*                    mypars,
                                 kernelconstant_interintra*   KerConst_interintra,
                                 kernelconstant_intracontrib* KerConst_intracontrib,
                                 kernelconstant_intra*        KerConst_intra,
                                 kernelconstant_rotlist*      KerConst_rotlist,
                                 kernelconstant_conform*      KerConst_conform,
                                 kernelconstant_grads*        KerConst_grads
                                );

void make_reqrot_ordering(
                          int *number_of_req_rotations,
                          int *atom_id_of_numrots,
                          int num_of_atoms
                         );

int gen_rotlist(
                Liganddata* myligand,
                int * rotlist,
                int   *rotlist_atomId,
                int   *rotlist_rotatableBondId
               );

#endif /* CALCENERGY_H_ */
