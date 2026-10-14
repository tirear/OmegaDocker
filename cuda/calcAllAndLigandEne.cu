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

__device__ void gpu_calc_all_and_ligand_energy(
                                float*  pGenotype,
                                float&  energy,
                                float&  ligandEnergy, // sans
//                                 int&    run_id, // sans useless
                                float3* calc_coords,
                                float*  pFloatAccumulator
                                // sans ==>
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
// The GPU device function calculates the energy of the entity described by genotype, dockpars and the liganddata
// arrays in constant memory and returns it in the energy parameter. The parameter run_id has to be equal to the ID
// of the run whose population includes the current entity (which can be determined with blockIdx.x), since this
// determines which reference orientation should be used.
{
  energy = 0.0f;
  ligandEnergy = 0.0f; // sans
  float energyTmp;
  float inRangeDist = cData.dockpars.collDistance * 0.5; // sans The good enough distance between atom and dummy atom with restraints.
// float EEE;
// float BBB;
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
    /* sans
    calc_coords[atom_id].x = cData.pKerconst_conform->ref_coords_const[3*atom_id];
    calc_coords[atom_id].y = cData.pKerconst_conform->ref_coords_const[3*atom_id+1];
    calc_coords[atom_id].z = cData.pKerconst_conform->ref_coords_const[3*atom_id+2];
    sans */
    // sans ==>
    calc_coords[atom_id].x = cData.pKerconst_conform_ref_coords_const[3*atom_id];
    calc_coords[atom_id].y = cData.pKerconst_conform_ref_coords_const[3*atom_id+1];
    calc_coords[atom_id].z = cData.pKerconst_conform_ref_coords_const[3*atom_id+2];
// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
//   printf("sans check13 COORCOOR !! x%f y%f z%f %d\n", calc_coords[atom_id].x, calc_coords[atom_id].y, calc_coords[atom_id].z, atom_id);
// }
    // <== sans
    // sans check ==>
//     if(threadIdx.x == 0)
//     {
//       printf("sans calcenergy 1 atom_id/calc_coords[atom_id].x/y/z: %d/%f/%f/%f\n", atom_id, calc_coords[atom_id].x, calc_coords[atom_id].y, calc_coords[atom_id].z);
//     }
    // <== sans check
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

  uint g1 = cData.dockpars.gridsize_x;
  uint g2 = cData.dockpars.gridsize_x_times_y;
  uint g3 = cData.dockpars.gridsize_x_times_y_times_z;

  __syncthreads();

  // ================================================
  // CALCULATING ATOMIC POSITIONS AFTER ROTATIONS
  // ================================================
  for (uint rotation_counter  = threadIdx.x;
            rotation_counter  < cData.dockpars.rotbondlist_length;
            rotation_counter += blockDim.x)
  {
//     int rotation_list_element = cData.pKerconst_rotlist->rotlist_const[rotation_counter]; // sans
    int rotation_list_element = cData.pKerconst_rotlist_rotlist_const[rotation_counter]; // sans
    // sans check ==>
//     if(threadIdx.x == 0)
//     {
//       printf("sans calcenergy 1-1 cData.dockpars.true_ligand_atoms/rotation_counter/cData.dockpars.rotbondlist_length/cData.pKerconst_rotlist_rotlist_const[rotation_counter]: %u/%u/%u/%d\n", cData.dockpars.true_ligand_atoms, rotation_counter, cData.dockpars.rotbondlist_length, cData.pKerconst_rotlist_rotlist_const[rotation_counter]);
//     }
    // <== sans check

    if ((rotation_list_element & RLIST_DUMMY_MASK) == 0) // If not dummy rotation
    {
//       uint atom_id = rotation_list_element & RLIST_ATOMID_MASK; // sans
      uint atom_id = cData.pKerconst_rotlist_rotlist_atomId_const[rotation_counter]; // sans

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
//         uint rotbond_id = (rotation_list_element & RLIST_RBONDID_MASK) >> RLIST_RBONDID_SHIFT; // sans
        uint rotbond_id = cData.pKerconst_rotlist_rotlist_rotatableBondId_const[rotation_counter]; // sans

        float rotation_angle = pGenotype[6+rotbond_id]*DEG_TO_RAD*0.5f;
        float s = sin(rotation_angle);
        /* sans
        rotation_unitvec.x = s*cData.pKerconst_conform->rotbonds_unit_vectors_const[3*rotbond_id];
        rotation_unitvec.y = s*cData.pKerconst_conform->rotbonds_unit_vectors_const[3*rotbond_id+1];
        rotation_unitvec.z = s*cData.pKerconst_conform->rotbonds_unit_vectors_const[3*rotbond_id+2];
        sans */
        // sans ==>
        rotation_unitvec.x = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id];
        rotation_unitvec.y = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id+1];
        rotation_unitvec.z = s*cData.pKerconst_conform_rotbonds_unit_vectors_const[3*rotbond_id+2];
        // <== sans
        rotation_unitvec.w = cos(rotation_angle);
        /* sans
        rotation_movingvec.x = cData.pKerconst_conform->rotbonds_moving_vectors_const[3*rotbond_id];
        rotation_movingvec.y = cData.pKerconst_conform->rotbonds_moving_vectors_const[3*rotbond_id+1];
        rotation_movingvec.z = cData.pKerconst_conform->rotbonds_moving_vectors_const[3*rotbond_id+2];
        sans */
        // sans ==>
        rotation_movingvec.x = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id];
        rotation_movingvec.y = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id+1];
        rotation_movingvec.z = cData.pKerconst_conform_rotbonds_moving_vectors_const[3*rotbond_id+2];
        // <== sans
        // Performing additionally the first movement which
        // is needed only if rotating around rotatable bond
        atom_to_rotate.x -= rotation_movingvec.x;
        atom_to_rotate.y -= rotation_movingvec.y;
        atom_to_rotate.z -= rotation_movingvec.z;
        // sans check ==>
//         if(threadIdx.x == 0)
//         {
//           printf("sans calcenergy 1-1-1 rotbond_id/rotation_angle/rotation_unitvec.x/rotation_unitvec.y/rotation_unitvec.z: %u/%f/%f/%f/%f\n", rotbond_id, rotation_angle, rotation_unitvec.x, rotation_unitvec.y, rotation_unitvec.z);
//         }
        // <== sans check
      }

      // Performing rotation and final movement
      float4 qt = quaternion_rotate(atom_to_rotate, rotation_unitvec);
      calc_coords[atom_id].x = qt.x + rotation_movingvec.x;
      calc_coords[atom_id].y = qt.y + rotation_movingvec.y;
      calc_coords[atom_id].z = qt.z + rotation_movingvec.z;
      // sans check ==>
//       if(threadIdx.x == 0)
//       {
//         printf("sans calcenergy 2 atom_id/calc_coords[atom_id].x/y/z: %d/%f/%f/%f\n", atom_id, calc_coords[atom_id].x, calc_coords[atom_id].y, calc_coords[atom_id].z);
//       }
      // <== sans check
    } // End if-statement not dummy rotation

      __syncthreads();

  } // End rotation_counter for-loop

  // ================================================
  // CALCULATING INTERMOLECULAR ENERGY
  // ================================================
  float weights[8];
  float cube[8];
  for (uint atom_id = threadIdx.x;
            atom_id < cData.dockpars.num_of_atoms;
            atom_id+= blockDim.x)
  {
//     if (cData.pKerconst_interintra->ignore_inter_const[atom_id]>0) // first two atoms of a flex res are to be ignored here // sans
    if (cData.pKerconst_interintra_ignore_inter_const[atom_id]>0) // first two atoms of a flex res are to be ignored here // sans
      continue;
    float x = calc_coords[atom_id].x;
    float y = calc_coords[atom_id].y;
    float z = calc_coords[atom_id].z;
    /* sans
    float q = cData.pKerconst_interintra->atom_charges_const[atom_id];
    uint atom_typeid = cData.pKerconst_interintra->atom_types_map_const[atom_id];
    sans */
    // sans ==>
    float q = cData.pKerconst_interintra_atom_charges_const[atom_id];
    uint atom_typeid = cData.pKerconst_interintra_atom_types_map_const[atom_id];
    // <== sans
// if(threadIdx.x == 0 && blockIdx.x < 1)
// {
  // int ccc = 0;
  // if(atom_id == 0)
  // {
  //   printf("sans COORCOOR x y z: %f %f %f\n", x, y, z);
  //   for (int aa = 0; aa < cData.dockpars.num_of_atoms; aa++)
  //   {
  //     if ((x < 0) || (y < 0) || (z < 0) || (x >= cData.dockpars.gridsize_x-1)
  //                                     || (y >= cData.dockpars.gridsize_y-1)
  //                                     || (z >= cData.dockpars.gridsize_z-1))
  //     {
  //       printf("sans COORCOOR x y z: %d %f %f %f WAOH\n", aa, calc_coords[aa].x, calc_coords[aa].y, calc_coords[aa].z);
  //       ccc++;
  //     }
  //     else
  //     {
  //       printf("sans COORCOOR x y z: %d %f %f %f\n", aa, calc_coords[aa].x, calc_coords[aa].y, calc_coords[aa].z);
  //     }
  //   }
  // }
  // printf("sans COORCOOR cData.dockpars.gridsize: %d %d %d WAOH%d\n", cData.dockpars.gridsize_x, cData.dockpars.gridsize_y, cData.dockpars.gridsize_z, ccc);
  // printf("sans COORCOOR cData.dockpars.gridsize: %d %d %d\n", cData.dockpars.gridsize_x, cData.dockpars.gridsize_y, cData.dockpars.gridsize_z);
// }
    if ((x < 0) || (y < 0) || (z < 0) || (x >= cData.dockpars.gridsize_x-1)
                                      || (y >= cData.dockpars.gridsize_y-1)
                                      || (z >= cData.dockpars.gridsize_z-1)){
// printf("sans COORCOOR outgrid x y z: %f %f %f id%d %d\n", x, y, z, atom_id, cData.dockpars.gridsize_x);
// sans ==>
#ifdef RESTORING_MAP_GRADIENT
      x -= 0.5f * cData.dockpars.gridsize_x;
      y -= 0.5f * cData.dockpars.gridsize_y;
      z -= 0.5f * cData.dockpars.gridsize_z;
      energy += 21.0f * (x*x+y*y+z*z); //100000.0f;
      #if defined (DEBUG_ENERGY_KERNEL)
      interE += 21.0f * (x*x+y*y+z*z);
      #endif
#else
// <== sans
// EEE = energy;
      energy += 65536.0f; //100000.0f;
// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
// printf("sans COORCOOR EEE: %f %f\n", EEE, energy);
// }
// sans ==>
      #if defined (DEBUG_ENERGY_KERNEL)
      interE += 16777216.0f; //100000.0f;
      #endif
#endif
// <== sans
      continue; // get on with loop as our work here is done (we crashed into the walls)
    }
// printf("sans COORCOOR ingrid x y z: %f %f %f id%d %d\n", x, y, z, atom_id, cData.dockpars.gridsize_x);
    // Getting coordinates
    float x_low  = floor(x);
    float y_low  = floor(y);
    float z_low  = floor(z);

    // Grid value at 000
    float* grid_value_000 = cData.pMem_fgrids + ((ulong)(x_low  + y_low*g1  + z_low*g2)<<2);

    float dx = x - x_low;
    float omdx = 1.0f - dx;
    float dy = y - y_low;
    float omdy = 1.0f - dy;
    float dz = z - z_low;
    float omdz = 1.0f - dz;

    // Calculating interpolation weights
    weights [idx_000] = omdx*omdy*omdz;
    weights [idx_010] = omdx*dy*omdz;
    weights [idx_001] = omdx*omdy*dz;
    weights [idx_011] = omdx*dy*dz;
    weights [idx_100] = dx*omdy*omdz;
    weights [idx_110] = dx*dy*omdz;
    weights [idx_101] = dx*omdy*dz;
    weights [idx_111] = dx*dy*dz;

    ulong mul_tmp = atom_typeid*g3<<2;
    cube[0] = *(grid_value_000+mul_tmp+0);
    cube[1] = *(grid_value_000+mul_tmp+1);
    cube[2] = *(grid_value_000+mul_tmp+2);
    cube[3] = *(grid_value_000+mul_tmp+3);
    cube[4] = *(grid_value_000+mul_tmp+4);
    cube[5] = *(grid_value_000+mul_tmp+5);
    cube[6] = *(grid_value_000+mul_tmp+6);
    cube[7] = *(grid_value_000+mul_tmp+7);

    // Calculating affinity energy
    energy += cube[0]*weights[0] + cube[1]*weights[1] + cube[2]*weights[2] + cube[3]*weights[3] + cube[4]*weights[4] + cube[5]*weights[5] + cube[6]*weights[6] + cube[7]*weights[7];
    #if defined (DEBUG_ENERGY_KERNEL)
    interE += cube[0]*weights[0] + cube[1]*weights[1] + cube[2]*weights[2] + cube[3]*weights[3] + cube[4]*weights[4] + cube[5]*weights[5] + cube[6]*weights[6] + cube[7]*weights[7];
    #endif

    // Capturing electrostatic values
    atom_typeid = cData.dockpars.num_of_map_atypes;

    mul_tmp = atom_typeid*g3<<2; // different atom type id to get charge IA
    cube[0] = *(grid_value_000+mul_tmp+0);
    cube[1] = *(grid_value_000+mul_tmp+1);
    cube[2] = *(grid_value_000+mul_tmp+2);
    cube[3] = *(grid_value_000+mul_tmp+3);
    cube[4] = *(grid_value_000+mul_tmp+4);
    cube[5] = *(grid_value_000+mul_tmp+5);
    cube[6] = *(grid_value_000+mul_tmp+6);
    cube[7] = *(grid_value_000+mul_tmp+7);

    // Calculating affinity energy
    energy += q * (cube[0]*weights[0] + cube[1]*weights[1] + cube[2]*weights[2] + cube[3]*weights[3] + cube[4]*weights[4] + cube[5]*weights[5] + cube[6]*weights[6] + cube[7]*weights[7]);
    #if defined (DEBUG_ENERGY_KERNEL)
    interE += q *(cube[0]*weights[0] + cube[1]*weights[1] + cube[2]*weights[2] + cube[3]*weights[3] + cube[4]*weights[4] + cube[5]*weights[5] + cube[6]*weights[6] + cube[7]*weights[7]);
    #endif

    // Need only magnitude of charge from here on down
    q = fabs(q);
    // Capturing desolvation values (atom_typeid+1 compared to above => mul_tmp + g3*4)
    mul_tmp += g3<<2;
    cube[0] = *(grid_value_000+mul_tmp+0);
    cube[1] = *(grid_value_000+mul_tmp+1);
    cube[2] = *(grid_value_000+mul_tmp+2);
    cube[3] = *(grid_value_000+mul_tmp+3);
    cube[4] = *(grid_value_000+mul_tmp+4);
    cube[5] = *(grid_value_000+mul_tmp+5);
    cube[6] = *(grid_value_000+mul_tmp+6);
    cube[7] = *(grid_value_000+mul_tmp+7);

    // Calculating affinity energy
    energy += q * (cube[0]*weights[0] + cube[1]*weights[1] + cube[2]*weights[2] + cube[3]*weights[3] + cube[4]*weights[4] + cube[5]*weights[5] + cube[6]*weights[6] + cube[7]*weights[7]);
    #if defined (DEBUG_ENERGY_KERNEL)
    interE += q *(cube[0]*weights[0] + cube[1]*weights[1] + cube[2]*weights[2] + cube[3]*weights[3] + cube[4]*weights[4] + cube[5]*weights[5] + cube[6]*weights[6] + cube[7]*weights[7]);
    #endif
  } // End atom_id for-loop (INTERMOLECULAR ENERGY)
#if defined (DEBUG_ENERGY_KERNEL)
  REDUCEFLOATSUM(interE, pFloatAccumulator)
#endif

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
    /* sans
    uint32_t atom1_id = cData.pKerconst_intracontrib->intraE_contributors_const[2*contributor_counter];
    uint32_t atom2_id = cData.pKerconst_intracontrib->intraE_contributors_const[2*contributor_counter+1];
    sans */
    // sans =>
    uint32_t atom1_id = cData.pKerconst_intracontrib_intraE_contributors_const[2*contributor_counter];
    uint32_t atom2_id = cData.pKerconst_intracontrib_intraE_contributors_const[2*contributor_counter+1];
    // <== sans

    // Calculating vector components of vector going
    // from first atom's to second atom's coordinates
    float subx = calc_coords[atom1_id].x - calc_coords[atom2_id].x;
    float suby = calc_coords[atom1_id].y - calc_coords[atom2_id].y;
    float subz = calc_coords[atom1_id].z - calc_coords[atom2_id].z;

    // Calculating atomic_distance
    float dist = sqrt(subx*subx + suby*suby + subz*subz);
    float atomic_distance = dist * cData.dockpars.grid_spacing;

    // Getting type IDs
    /* sans
    uint32_t atom1_typeid = cData.pKerconst_interintra->atom_types_const[atom1_id];
    uint32_t atom2_typeid = cData.pKerconst_interintra->atom_types_const[atom2_id];
    sans */
    // sans ==>
    uint32_t atom1_typeid = cData.pKerconst_interintra_atom_types_const[atom1_id];
    uint32_t atom2_typeid = cData.pKerconst_interintra_atom_types_const[atom2_id];
    // <== sans

    /* sans
    uint32_t atom1_type_vdw_hb = cData.pKerconst_intra->atom_types_reqm_const [atom1_typeid];
    uint32_t atom2_type_vdw_hb = cData.pKerconst_intra->atom_types_reqm_const [atom2_typeid];
    sans */
    // sans ==>
    uint32_t atom1_type_vdw_hb = cData.pKerconst_intra_atom_types_reqm_const[atom1_typeid];
    uint32_t atom2_type_vdw_hb = cData.pKerconst_intra_atom_types_reqm_const[atom2_typeid];
    // <== sans

    // ------------------------------------------------
    // Required only for flexrings
    // Checking if this is a CG-G0 atomic pair.
    // If so, then adding energy term (E = G * distance).
    // Initial specification required NON-SMOOTHED distance.
    // This interaction is evaluated at any distance,
    // so no cuttoffs considered here!
    // vbond is G when calculating flexrings, 0.0 otherwise
    /* sans dna/rna
    float vbond = G * (float)(((atom1_type_vdw_hb == ATYPE_CG_IDX) && (atom2_type_vdw_hb == ATYPE_G0_IDX)) ||
            ((atom1_type_vdw_hb == ATYPE_G0_IDX) && (atom2_type_vdw_hb == ATYPE_CG_IDX)));
    energy += vbond * atomic_distance;
    sans dna/rna */
    // sans dna/rna ==>
    if
    (
      (atom1_type_vdw_hb > ATYPE_G0_IDX && atom2_type_vdw_hb == ATYPE_G0_IDX)
      ||
      (atom1_type_vdw_hb == ATYPE_G0_IDX && atom2_type_vdw_hb > ATYPE_G0_IDX)
    )
    {
      // float vbond = G;
      // energy += vbond * atomic_distance * atomic_distance;
      // sans The potential energy of restrained atom-atom distance should not be used to judge the pose is good or not. Therefore the above 2 lines should be marked.

      // sans Any distance between restrainted atoms should be not larger than the "inRangeDist". If this is the case, the energy value is forced to increased, allowing for a higher chance of rebirth.
      if(atomic_distance > inRangeDist)
      {
        // printf("sans atomDist: %d %d %d %f %f %f\n",blockIdx.x, atom1_id+1, atom2_id+1, atomic_distance, calc_coords[atom1_id].x, calc_coords[atom2_id].x);
        energy += 65536.0f;
        ligandEnergy += 65536.0f;
      }
    }
    else
    {
    // <== sans dna/rna

      // ------------------------------------------------
      // Calculating energy contributions
      // Cuttoff1: internuclear-distance at 8A only for vdw and hbond.
      if (atomic_distance < 8.0f)
      {
        uint32_t idx = atom1_typeid * cData.dockpars.num_of_atypes + atom2_typeid;
//         ushort exps = cData.pKerconst_intra->VWpars_exp_const[idx]; // sans
        ushort exps = cData.pKerconst_intra_VWpars_exp_const[idx]; // sans
        char m=(exps & 0xFF00)>>8;
        char n=(exps & 0xFF);
        // Getting optimum pair distance (opt_distance)
//         float opt_distance = cData.pKerconst_intra->reqm_AB_const[idx]; // sans
        float opt_distance = cData.pKerconst_intra_reqm_AB_const[idx]; // sans

        // Getting smoothed distance
        // smoothed_distance = function(atomic_distance, opt_distance)
        float opt_dist_delta = opt_distance - atomic_distance;
        if(fabs(opt_dist_delta)>=delta_distance){
          smoothed_distance = atomic_distance + copysign(delta_distance,opt_dist_delta);
        } else smoothed_distance = opt_distance;
        // Calculating van der Waals / hydrogen bond term
        /* sans
        energy += (cData.pKerconst_intra->VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra->VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        sans */
        // sans ==>
        energyTmp = (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra_VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        energy += energyTmp;
        ligandEnergy += energyTmp;
        // <== sans
        #if defined (DEBUG_ENERGY_KERNEL)
        /* sans
        intraE += (cData.pKerconst_intra->VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra->VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        sans */
        // sans ==>
        intraE += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -__powf(smoothed_distance,m-n)*cData.pKerconst_intra_VWpars_BD_const[idx])
                   *__powf(smoothed_distance,-m);
        // <== sans
        #endif
      } // if cuttoff1 - internuclear-distance at 8A

      // Calculating energy contributions
      // Cuttoff2: internuclear-distance at 20.48A only for el and sol.
      if (atomic_distance < 20.48f)
      {
        if(atomic_distance<cData.dockpars.elec_min_distance)
          atomic_distance=cData.dockpars.elec_min_distance;
        /* sans
        float q1 = cData.pKerconst_interintra->atom_charges_const[atom1_id];
        float q2 = cData.pKerconst_interintra->atom_charges_const[atom2_id];
        sans */
        // sans ==>
        float q1 = cData.pKerconst_interintra_atom_charges_const[atom1_id];
        float q2 = cData.pKerconst_interintra_atom_charges_const[atom2_id];
        // <== sans
//        float exp_el = native_exp(DIEL_B_TIMES_H*atomic_distance);
        float dist2 = atomic_distance*atomic_distance;
        // Calculating desolvation term
        // 1/25.92 = 0.038580246913580245
        /* sans
        float desolv_energy =  ((cData.pKerconst_intra->dspars_S_const[atom1_typeid] +
               cData.dockpars.qasp*fabs(q1)) * cData.pKerconst_intra->dspars_V_const[atom2_typeid] +
              (cData.pKerconst_intra->dspars_S_const[atom2_typeid] +
               cData.dockpars.qasp*fabs(q2)) * cData.pKerconst_intra->dspars_V_const[atom1_typeid]) *
               (
                cData.dockpars.coeff_desolv*(12.96f-0.1063f*dist2*(1.0f-0.001947f*dist2)) /
                (12.96f+dist2*(0.4137f+dist2*(0.00357f+0.000112f*dist2)))
               );
        sans */
        // sans ==>
        float desolv_energy =  ((cData.pKerconst_intra_dspars_S_const[atom1_typeid] +
               cData.dockpars.qasp*fabs(q1)) * cData.pKerconst_intra_dspars_V_const[atom2_typeid] +
              (cData.pKerconst_intra_dspars_S_const[atom2_typeid] +
               cData.dockpars.qasp*fabs(q2)) * cData.pKerconst_intra_dspars_V_const[atom1_typeid]) *
               (
                cData.dockpars.coeff_desolv*(12.96f-0.1063f*dist2*(1.0f-0.001947f*dist2)) /
                (12.96f+dist2*(0.4137f+dist2*(0.00357f+0.000112f*dist2)))
               );
        // <== sans
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
        // energy += diel * es_energy + desolv_energy; // sans
        // sans ==>
        energyTmp = diel * es_energy + desolv_energy;
        energy += energyTmp;
        ligandEnergy += energyTmp;
        // <== sans

        #if defined (DEBUG_ENERGY_KERNEL)
        intraE += diel * es_energy + desolv_energy;
        #endif
      } // if cuttoff2 - internuclear-distance at 20.48A
    } // sans
  } // End contributor_counter for-loop (INTRAMOLECULAR ENERGY)

  // reduction to calculate energy
// BBB = energy;
  REDUCEFLOATSUM(energy, pFloatAccumulator)
  REDUCEFLOATSUM(ligandEnergy, pFloatAccumulator)
// if(threadIdx.x == 0 && blockIdx.x == 0)
// {
// printf("sans COORCOOR FINAL: %f %f\n", BBB, energy);
// }
#if defined (DEBUG_ENERGY_KERNEL)
  REDUCEFLOATSUM(intraE, pFloatAccumulator)
#endif
}
