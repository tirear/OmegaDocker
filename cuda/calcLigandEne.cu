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


//#define DEBUG_ENERGY_KERNEL


// All related pragmas are in defines.h (accesible by host and device code)

__device__ void gpu_calc_ligand_energy(
                                float*  pGenotype,
                                float&  energy,
                                float3* calc_coords,
                                float*  pFloatAccumulator
                               )
// The GPU device function calculates the energy of the entity described by genotype, dockpars and the liganddata
// arrays in constant memory and returns it in the energy parameter. The parameter run_id has to be equal to the ID
// of the run whose population includes the current entity (which can be determined with blockIdx.x), since this
// determines which reference orientation should be used.
{
  energy = 0.0f;
  float inRangeDist = cData.dockpars.collDistance * 0.5;
#if defined (DEBUG_ENERGY_KERNEL)
  float interE = 0.0f;
  float intraE = 0.0f;
#endif

  // Initializing gradients (forces)
  // Derived from autodockdev/maps.py
  for (uint atom_id = threadIdx.x;
      atom_id < cData.dockpars.num_of_atoms;
      atom_id+= blockDim.x) {
    // Initialize coordinates
    calc_coords[atom_id].x = cData.pKerconst_conform_ref_coords_const[3*atom_id];
    calc_coords[atom_id].y = cData.pKerconst_conform_ref_coords_const[3*atom_id+1];
    calc_coords[atom_id].z = cData.pKerconst_conform_ref_coords_const[3*atom_id+2];
  }

  // General rotation moving vector
  float4 genrot_movingvec;
  genrot_movingvec.x = pGenotype[0];
  genrot_movingvec.y = pGenotype[1];
  genrot_movingvec.z = pGenotype[2];
  genrot_movingvec.w = 0.0f;
  // Convert orientation genes from sex. to radians
  float phi         = pGenotype[3] * DEG_TO_RAD;
  float theta       = pGenotype[4] * DEG_TO_RAD;
  float genrotangle = pGenotype[5] * DEG_TO_RAD;

  float4 genrot_unitvec;
  float sin_angle = sin(theta);
  float s2 = sin(genrotangle * 0.5f);
  genrot_unitvec.x = s2*sin_angle*cos(phi);
  genrot_unitvec.y = s2*sin_angle*sin(phi);
  genrot_unitvec.z = s2*cos(theta);
  genrot_unitvec.w = cos(genrotangle*0.5f);

  // uint g1 = cData.dockpars.gridsize_x;
  // uint g2 = cData.dockpars.gridsize_x_times_y;
  // uint g3 = cData.dockpars.gridsize_x_times_y_times_z;

  __syncthreads();

  // ================================================
  // CALCULATING ATOMIC POSITIONS AFTER ROTATIONS
  // ================================================
  for (uint rotation_counter  = threadIdx.x;
            rotation_counter  < cData.dockpars.rotbondlist_length;
            rotation_counter += blockDim.x)
  {
    int rotation_list_element = cData.pKerconst_rotlist_rotlist_const[rotation_counter];

    if ((rotation_list_element & RLIST_DUMMY_MASK) == 0) // If not dummy rotation
    {
      uint atom_id = cData.pKerconst_rotlist_rotlist_atomId_const[rotation_counter];

      // Capturing atom coordinates
      float4 atom_to_rotate;
      atom_to_rotate.x = calc_coords[atom_id].x;
      atom_to_rotate.y = calc_coords[atom_id].y;
      atom_to_rotate.z = calc_coords[atom_id].z;
      atom_to_rotate.w = 0.0f;

      // initialize with general rotation values
      float4 rotation_unitvec;
      float4 rotation_movingvec;
      if (atom_id < cData.dockpars.true_ligand_atoms){
        rotation_unitvec = genrot_unitvec;
        rotation_movingvec = genrot_movingvec;
      } else{
        rotation_unitvec.x = 0.0f; rotation_unitvec.y = 0.0f; rotation_unitvec.z = 0.0f;
        rotation_unitvec.w = 1.0f;
        rotation_movingvec.x = 0.0f; rotation_movingvec.y = 0.0f; rotation_movingvec.z = 0.0f;
        rotation_movingvec.w = 0.0f;
      }

      if ((rotation_list_element & RLIST_GENROT_MASK) == 0) // If rotating around rotatable bond
      {
        uint rotbond_id = cData.pKerconst_rotlist_rotlist_rotatableBondId_const[rotation_counter];

        float rotation_angle = pGenotype[6+rotbond_id]*DEG_TO_RAD*0.5f;
        float s = sin(rotation_angle);
        rotation_unitvec.x = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id];
        rotation_unitvec.y = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id+1];
        rotation_unitvec.z = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id+2];
        rotation_unitvec.w = cos(rotation_angle);
        rotation_movingvec.x = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id];
        rotation_movingvec.y = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id+1];
        rotation_movingvec.z = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id+2];
        // Performing additionally the first movement which
        // is needed only if rotating around rotatable bond
        atom_to_rotate.x -= rotation_movingvec.x;
        atom_to_rotate.y -= rotation_movingvec.y;
        atom_to_rotate.z -= rotation_movingvec.z;
      }

      // Performing rotation and final movement
      float4 qt = quaternion_rotate(atom_to_rotate, rotation_unitvec);
      calc_coords[atom_id].x = qt.x + rotation_movingvec.x;
      calc_coords[atom_id].y = qt.y + rotation_movingvec.y;
      calc_coords[atom_id].z = qt.z + rotation_movingvec.z;
    } // End if-statement not dummy rotation

      __syncthreads();

  } // End rotation_counter for-loop


  // In paper: intermolecular and internal energy calculation
  // are independent from each other, -> NO BARRIER NEEDED
  // but require different operations,
  // thus, they can be executed only sequentially on the GPU.
  float delta_distance = 0.5f * cData.dockpars.smooth;
  float smoothed_distance;

  // ================================================
  // CALCULATING INTRAMOLECULAR ENERGY
  // ================================================
  for (uint contributor_counter = threadIdx.x;
            contributor_counter < cData.dockpars.num_of_intraE_contributors;
            contributor_counter += blockDim.x)
  {
    // Getting atom IDs
    uint32_t atom1_id = cData.pKerconst_intracontrib_intraE_contributors_const[2*contributor_counter];
    uint32_t atom2_id = cData.pKerconst_intracontrib_intraE_contributors_const[2*contributor_counter+1];

    // Calculating vector components of vector going
    // from first atom's to second atom's coordinates
    float subx = calc_coords[atom1_id].x - calc_coords[atom2_id].x;
    float suby = calc_coords[atom1_id].y - calc_coords[atom2_id].y;
    float subz = calc_coords[atom1_id].z - calc_coords[atom2_id].z;

    // Calculating atomic_distance
    float dist = sqrt(subx*subx + suby*suby + subz*subz);
    float atomic_distance = dist * cData.dockpars.grid_spacing;

    // Getting type IDs
    uint32_t atom1_typeid = cData.pKerconst_interintra_atom_types_const[atom1_id];
    uint32_t atom2_typeid = cData.pKerconst_interintra_atom_types_const[atom2_id];

    uint32_t atom1_type_vdw_hb = cData.pKerconst_intra_atom_types_reqm_const[atom1_typeid];
    uint32_t atom2_type_vdw_hb = cData.pKerconst_intra_atom_types_reqm_const[atom2_typeid];

    // ------------------------------------------------
    // Required only for flexrings
    // Checking if this is a CG-G0 atomic pair.
    // If so, then adding energy term (E = G * distance).
    // Initial specification required NON-SMOOTHED distance.
    // This interaction is evaluated at any distance,
    // so no cuttoffs considered here!
    // vbond is G when calculating flexrings, 0.0 otherwise
    if
    (
      (
        (
          (atom1_type_vdw_hb == ATYPE_CG_IDX) ||
          (atom1_type_vdw_hb == ATYPE_HG_IDX) ||
          (atom1_type_vdw_hb == ATYPE_NG_IDX) ||
          (atom1_type_vdw_hb == ATYPE_OG_IDX)
        )
        &&
        (
          atom2_type_vdw_hb == ATYPE_G0_IDX
        )
      )
      ||
      (
        (
          atom1_type_vdw_hb == ATYPE_G0_IDX
        )
        &&
        (
          (atom2_type_vdw_hb == ATYPE_CG_IDX) ||
          (atom2_type_vdw_hb == ATYPE_HG_IDX) ||
          (atom2_type_vdw_hb == ATYPE_NG_IDX) ||
          (atom2_type_vdw_hb == ATYPE_OG_IDX)
        )
      )
    )
    {
      // float vbond = G;
      // energy += vbond * atomic_distance * atomic_distance;
      // The potential energy of restrained atom-atom distance should not be used to judge the pose is good or not. Therefore the above 2 lines should be marked.

      // Any distance between restrainted atoms should be not larger than the "inRangeDist". If this is the case, the energy value is forced to increased, allowing for a higher chance of rebirth.
      if(atomic_distance > inRangeDist)
      {
        energy += 65536.0f;
      }
    }
    else
    {
      
      // ------------------------------------------------
      // Calculating energy contributions
      // Cuttoff1: internuclear-distance at 8A only for vdw and hbond.
      if (atomic_distance < 8.0f)
      {
        uint32_t idx = atom1_typeid * cData.dockpars.num_of_atypes + atom2_typeid;
        ushort exps = cData.pKerconst_intra_VWpars_exp_const[idx];
        char m=(exps & 0xFF00)>>8;
        char n=(exps & 0xFF);
        // Getting optimum pair distance (opt_distance)
        float opt_distance = cData.pKerconst_intra_reqm_AB_const[idx];
      
        // Getting smoothed distance
        // smoothed_distance = function(atomic_distance, opt_distance)
        float opt_dist_delta = opt_distance - atomic_distance;
        if(fabs(opt_dist_delta)>=delta_distance){
          smoothed_distance = atomic_distance + copysign(delta_distance,opt_dist_delta);
        } else smoothed_distance = opt_distance;
        // Calculating van der Waals / hydrogen bond term
        energy += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra_VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        #if defined (DEBUG_ENERGY_KERNEL)
        intraE += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra_VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        #endif
      } // if cuttoff1 - internuclear-distance at 8A
      
      // Calculating energy contributions
      // Cuttoff2: internuclear-distance at 20.48A only for el and sol.
      if (atomic_distance < 20.48f)
      {
        if(atomic_distance<cData.dockpars.elec_min_distance)
          atomic_distance=cData.dockpars.elec_min_distance;
        float q1 = cData.pKerconst_interintra_atom_charges_const[atom1_id];
        float q2 = cData.pKerconst_interintra_atom_charges_const[atom2_id];
//        float exp_el = native_exp(DIEL_B_TIMES_H*atomic_distance);
        float dist2 = atomic_distance*atomic_distance;
        // Calculating desolvation term
        // 1/25.92 = 0.038580246913580245
        float desolv_energy =  ((cData.pKerconst_intra_dspars_S_const[atom1_typeid] +
               cData.dockpars.qasp*fabs(q1)) * cData.pKerconst_intra_dspars_V_const[atom2_typeid] +
              (cData.pKerconst_intra_dspars_S_const[atom2_typeid] +
               cData.dockpars.qasp*fabs(q2)) * cData.pKerconst_intra_dspars_V_const[atom1_typeid]) *
               (
                cData.dockpars.coeff_desolv*(12.96f-0.1063f*dist2*(1.0f-0.001947f*dist2)) /
                (12.96f+dist2*(0.4137f+dist2*(0.00357f+0.000112f*dist2)))
               );
        // Calculating electrostatic term
        // float dist_shift=atomic_distance+1.26366f;
        // dist2=dist_shift*dist_shift;
        // float diel = (1.10859f / dist2)+0.010358f;
        #ifndef DIEL_FIT_ABC
        float dist_shift=atomic_distance+1.26366f;
        dist2=dist_shift*dist_shift;
        float diel = 1.10859f / dist2 + 0.010358f;
        #else
        float dist_shift=atomic_distance+1.588f;
        dist2=dist_shift*dist_shift;
        float disth_shift=atomic_distance+0.794f;
        float disth4=disth_shift*disth_shift;
        disth4*=disth4;
        float diel = 1.404f / dist2 + 0.072f / disth4 + 0.00831f;
        #endif
        float es_energy = cData.dockpars.coeff_elec * q1 * q2 / atomic_distance;
        energy += diel * es_energy + desolv_energy;
      
        #if defined (DEBUG_ENERGY_KERNEL)
        intraE += diel * es_energy + desolv_energy;
        #endif
      } // if cuttoff2 - internuclear-distance at 20.48A
    }
  } // End contributor_counter for-loop (INTRAMOLECULAR ENERGY)

  // reduction to calculate energy
  REDUCEFLOATSUM(energy, pFloatAccumulator)
#if defined (DEBUG_ENERGY_KERNEL)
  REDUCEFLOATSUM(intraE, pFloatAccumulator)
#endif
}




__device__ void gpu_calc_ligand_energy_final(
                                float*  pGenotype,
                                float&  energy,
                                float*  energyx,
                                float3* calc_coords,
                                float*  pFloatAccumulator
                               )
// The GPU device function calculates the energy of the entity described by genotype, dockpars and the liganddata
// arrays in constant memory and returns it in the energy parameter. The parameter run_id has to be equal to the ID
// of the run whose population includes the current entity (which can be determined with blockIdx.x), since this
// determines which reference orientation should be used.
{
  for(int i = 0; i< 5; i++)
  {
    energyx[i] = 0.0f;
  }
  energy = 0.0f;
  float inRangeDist = cData.dockpars.collDistance * 0.5;
#if defined (DEBUG_ENERGY_KERNEL)
  float interE = 0.0f;
  float intraE = 0.0f;
#endif

  // Initializing gradients (forces)
  // Derived from autodockdev/maps.py
  for (uint atom_id = threadIdx.x;
      atom_id < cData.dockpars.num_of_atoms;
      atom_id+= blockDim.x) {
    // Initialize coordinates
    calc_coords[atom_id].x = cData.pKerconst_conform_ref_coords_const[3*atom_id];
    calc_coords[atom_id].y = cData.pKerconst_conform_ref_coords_const[3*atom_id+1];
    calc_coords[atom_id].z = cData.pKerconst_conform_ref_coords_const[3*atom_id+2];
  }

  // General rotation moving vector
  float4 genrot_movingvec;
  genrot_movingvec.x = pGenotype[0];
  genrot_movingvec.y = pGenotype[1];
  genrot_movingvec.z = pGenotype[2];
  genrot_movingvec.w = 0.0f;
  // Convert orientation genes from sex. to radians
  float phi         = pGenotype[3] * DEG_TO_RAD;
  float theta       = pGenotype[4] * DEG_TO_RAD;
  float genrotangle = pGenotype[5] * DEG_TO_RAD;

  float4 genrot_unitvec;
  float sin_angle = sin(theta);
  float s2 = sin(genrotangle * 0.5f);
  genrot_unitvec.x = s2*sin_angle*cos(phi);
  genrot_unitvec.y = s2*sin_angle*sin(phi);
  genrot_unitvec.z = s2*cos(theta);
  genrot_unitvec.w = cos(genrotangle*0.5f);

  // uint g1 = cData.dockpars.gridsize_x;
  // uint g2 = cData.dockpars.gridsize_x_times_y;
  // uint g3 = cData.dockpars.gridsize_x_times_y_times_z;

  __syncthreads();

  // ================================================
  // CALCULATING ATOMIC POSITIONS AFTER ROTATIONS
  // ================================================
  for (uint rotation_counter  = threadIdx.x;
            rotation_counter  < cData.dockpars.rotbondlist_length;
            rotation_counter += blockDim.x)
  {
    int rotation_list_element = cData.pKerconst_rotlist_rotlist_const[rotation_counter];


    if ((rotation_list_element & RLIST_DUMMY_MASK) == 0) // If not dummy rotation
    {
      uint atom_id = cData.pKerconst_rotlist_rotlist_atomId_const[rotation_counter];

      // Capturing atom coordinates
      float4 atom_to_rotate;
      atom_to_rotate.x = calc_coords[atom_id].x;
      atom_to_rotate.y = calc_coords[atom_id].y;
      atom_to_rotate.z = calc_coords[atom_id].z;
      atom_to_rotate.w = 0.0f;

      // initialize with general rotation values
      float4 rotation_unitvec;
      float4 rotation_movingvec;
      if (atom_id < cData.dockpars.true_ligand_atoms){
        rotation_unitvec = genrot_unitvec;
        rotation_movingvec = genrot_movingvec;
      } else{
        rotation_unitvec.x = 0.0f; rotation_unitvec.y = 0.0f; rotation_unitvec.z = 0.0f;
        rotation_unitvec.w = 1.0f;
        rotation_movingvec.x = 0.0f; rotation_movingvec.y = 0.0f; rotation_movingvec.z = 0.0f;
        rotation_movingvec.w = 0.0f;
      }

      if ((rotation_list_element & RLIST_GENROT_MASK) == 0) // If rotating around rotatable bond
      {
        uint rotbond_id = cData.pKerconst_rotlist_rotlist_rotatableBondId_const[rotation_counter];

        float rotation_angle = pGenotype[6+rotbond_id]*DEG_TO_RAD*0.5f;
        float s = sin(rotation_angle);
        rotation_unitvec.x = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id];
        rotation_unitvec.y = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id+1];
        rotation_unitvec.z = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id+2];
        rotation_unitvec.w = cos(rotation_angle);
        rotation_movingvec.x = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id];
        rotation_movingvec.y = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id+1];
        rotation_movingvec.z = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id+2];
        // Performing additionally the first movement which
        // is needed only if rotating around rotatable bond
        atom_to_rotate.x -= rotation_movingvec.x;
        atom_to_rotate.y -= rotation_movingvec.y;
        atom_to_rotate.z -= rotation_movingvec.z;
      }

      // Performing rotation and final movement
      float4 qt = quaternion_rotate(atom_to_rotate, rotation_unitvec);
      calc_coords[atom_id].x = qt.x + rotation_movingvec.x;
      calc_coords[atom_id].y = qt.y + rotation_movingvec.y;
      calc_coords[atom_id].z = qt.z + rotation_movingvec.z;
    } // End if-statement not dummy rotation

      __syncthreads();

  } // End rotation_counter for-loop

  // In paper: intermolecular and internal energy calculation
  // are independent from each other, -> NO BARRIER NEEDED
  // but require different operations,
  // thus, they can be executed only sequentially on the GPU.
  float delta_distance = 0.5f * cData.dockpars.smooth;
  float smoothed_distance;

  // ================================================
  // CALCULATING INTRAMOLECULAR ENERGY
  // ================================================
  for (uint contributor_counter = threadIdx.x;
            contributor_counter < cData.dockpars.num_of_intraE_contributors;
            contributor_counter += blockDim.x)
  {
    // Getting atom IDs
    uint32_t atom1_id = cData.pKerconst_intracontrib_intraE_contributors_const[2*contributor_counter];
    uint32_t atom2_id = cData.pKerconst_intracontrib_intraE_contributors_const[2*contributor_counter+1];

    // Calculating vector components of vector going
    // from first atom's to second atom's coordinates
    float subx = calc_coords[atom1_id].x - calc_coords[atom2_id].x;
    float suby = calc_coords[atom1_id].y - calc_coords[atom2_id].y;
    float subz = calc_coords[atom1_id].z - calc_coords[atom2_id].z;

    // Calculating atomic_distance
    float dist = sqrt(subx*subx + suby*suby + subz*subz);
    float atomic_distance = dist * cData.dockpars.grid_spacing;

    // Getting type IDs
    uint32_t atom1_typeid = cData.pKerconst_interintra_atom_types_const[atom1_id];
    uint32_t atom2_typeid = cData.pKerconst_interintra_atom_types_const[atom2_id];

    uint32_t atom1_type_vdw_hb = cData.pKerconst_intra_atom_types_reqm_const[atom1_typeid];
    uint32_t atom2_type_vdw_hb = cData.pKerconst_intra_atom_types_reqm_const[atom2_typeid];

    // ------------------------------------------------
    // Required only for flexrings
    // Checking if this is a CG-G0 atomic pair.
    // If so, then adding energy term (E = G * distance).
    // Initial specification required NON-SMOOTHED distance.
    // This interaction is evaluated at any distance,
    // so no cuttoffs considered here!
    // vbond is G when calculating flexrings, 0.0 otherwise
    if
    (
      (
        (
          (atom1_type_vdw_hb == ATYPE_CG_IDX) ||
          (atom1_type_vdw_hb == ATYPE_HG_IDX) ||
          (atom1_type_vdw_hb == ATYPE_NG_IDX) ||
          (atom1_type_vdw_hb == ATYPE_OG_IDX)
        )
        &&
        (
          atom2_type_vdw_hb == ATYPE_G0_IDX
        )
      )
      ||
      (
        (
          atom1_type_vdw_hb == ATYPE_G0_IDX
        )
        &&
        (
          (atom2_type_vdw_hb == ATYPE_CG_IDX) ||
          (atom2_type_vdw_hb == ATYPE_HG_IDX) ||
          (atom2_type_vdw_hb == ATYPE_NG_IDX) ||
          (atom2_type_vdw_hb == ATYPE_OG_IDX)
        )
      )
    )
    {
      // float vbond = G;
      // energy += vbond * atomic_distance * atomic_distance;
      // The potential energy of restrained atom-atom distance should not be used to judge the pose is good or not. Therefore the above 2 lines should be marked.

      // Any distance between restrainted atoms should be not larger than the "inRangeDist". If this is the case, the energy value is forced to increased, allowing for a higher chance of rebirth.
      if(atomic_distance > inRangeDist)
      {
        energy += 65536.0f;
      }
    }
    else
    {
      
      // ------------------------------------------------
      // Calculating energy contributions
      // Cuttoff1: internuclear-distance at 8A only for vdw and hbond.
      if (atomic_distance < 8.0f)
      {
        uint32_t idx = atom1_typeid * cData.dockpars.num_of_atypes + atom2_typeid;
        ushort exps = cData.pKerconst_intra_VWpars_exp_const[idx];
        char m=(exps & 0xFF00)>>8;
        char n=(exps & 0xFF);
        // Getting optimum pair distance (opt_distance)
        float opt_distance = cData.pKerconst_intra_reqm_AB_const[idx];
      
        // Getting smoothed distance
        // smoothed_distance = function(atomic_distance, opt_distance)
        float opt_dist_delta = opt_distance - atomic_distance;
        if(fabs(opt_dist_delta)>=delta_distance){
          smoothed_distance = atomic_distance + copysign(delta_distance,opt_dist_delta);
        } else smoothed_distance = opt_distance;
        // Calculating van der Waals / hydrogen bond term
        energyx[3] += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra_VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        energy += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra_VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        #if defined (DEBUG_ENERGY_KERNEL)
        intraE += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra_VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        #endif
      } // if cuttoff1 - internuclear-distance at 8A
      
      // Calculating energy contributions
      // Cuttoff2: internuclear-distance at 20.48A only for el and sol.
      if (atomic_distance < 20.48f)
      {
        if(atomic_distance<cData.dockpars.elec_min_distance)
          atomic_distance=cData.dockpars.elec_min_distance;
        float q1 = cData.pKerconst_interintra_atom_charges_const[atom1_id];
        float q2 = cData.pKerconst_interintra_atom_charges_const[atom2_id];
//        float exp_el = native_exp(DIEL_B_TIMES_H*atomic_distance);
        float dist2 = atomic_distance*atomic_distance;
        // Calculating desolvation term
        // 1/25.92 = 0.038580246913580245
        float desolv_energy =  ((cData.pKerconst_intra_dspars_S_const[atom1_typeid] +
               cData.dockpars.qasp*fabs(q1)) * cData.pKerconst_intra_dspars_V_const[atom2_typeid] +
              (cData.pKerconst_intra_dspars_S_const[atom2_typeid] +
               cData.dockpars.qasp*fabs(q2)) * cData.pKerconst_intra_dspars_V_const[atom1_typeid]) *
               (
                cData.dockpars.coeff_desolv*(12.96f-0.1063f*dist2*(1.0f-0.001947f*dist2)) /
                (12.96f+dist2*(0.4137f+dist2*(0.00357f+0.000112f*dist2)))
               );
        // Calculating electrostatic term
        // float dist_shift=atomic_distance+1.26366f;
        // dist2=dist_shift*dist_shift;
        // float diel = (1.10859f / dist2)+0.010358f;
        #ifndef DIEL_FIT_ABC
        float dist_shift=atomic_distance+1.26366f;
        dist2=dist_shift*dist_shift;
        float diel = 1.10859f / dist2 + 0.010358f;
        #else
        float dist_shift=atomic_distance+1.588f;
        dist2=dist_shift*dist_shift;
        float disth_shift=atomic_distance+0.794f;
        float disth4=disth_shift*disth_shift;
        disth4*=disth4;
        float diel = 1.404f / dist2 + 0.072f / disth4 + 0.00831f;
        #endif
        float es_energy = cData.dockpars.coeff_elec * q1 * q2 / atomic_distance;
        energyx[4] += diel * es_energy;
        energyx[5] += desolv_energy;
        energy += diel * es_energy + desolv_energy;
      
        #if defined (DEBUG_ENERGY_KERNEL)
        intraE += diel * es_energy + desolv_energy;
        #endif
      } // if cuttoff2 - internuclear-distance at 20.48A
    }
  } // End contributor_counter for-loop (INTRAMOLECULAR ENERGY)

  // reduction to calculate energy
// BBB = energy;
  REDUCEFLOATSUM(energy, pFloatAccumulator)
  REDUCEFLOATSUM(energyx[0], pFloatAccumulator)
  REDUCEFLOATSUM(energyx[1], pFloatAccumulator)
  REDUCEFLOATSUM(energyx[2], pFloatAccumulator)
  REDUCEFLOATSUM(energyx[3], pFloatAccumulator)
  REDUCEFLOATSUM(energyx[4], pFloatAccumulator)
  REDUCEFLOATSUM(energyx[5], pFloatAccumulator)
#if defined (DEBUG_ENERGY_KERNEL)
  REDUCEFLOATSUM(intraE, pFloatAccumulator)
#endif
}
