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


// IMPORTANT: The following block contains definitions
// already made either in energy or gradient calculation files.
// For that reason, these are commented here.

// IMPORTANT: the code of gradient calculation was the initial template.
// Then, statements corresponding to enery calculations were added gradually.
// The latter can be distinguised this way: they are place within lines without indentation.

#define CONVERT_INTO_ANGSTROM_RADIAN  // DO NOT UNDEFINE, NO REALLY! DO!!! NOT!!! UNDEFINE!!! SML 200608
#define SCFACTOR_ANGSTROM_RADIAN (1.0f/(DEG_TO_RAD * DEG_TO_RAD))

// Enables full floating point gradient calculation.
// Use is not advised as:
// - the determinism gradients (aka integer gradients) are much faster *and*
// - speed up the local search convergence
// Please only use for debugging
// #define FLOAT_GRADIENTS

// Enable restoring map gradient
// Currently, this is not a good idea
// #define RESTORING_MAP_GRADIENT

__device__ void gpu_calc_energrad_ligand2(
                                  float*  genotype,
                                  float&  global_energy,
                                  float3* calc_coords,
#if defined (DEBUG_ENERGY_KERNEL)
                                  float&  interE,
                                  float&  pintraE,
#endif
#ifdef FLOAT_GRADIENTS
                                  float3* gradient,
#else
                                  int3*   gradient,
#endif
                                  float*  fgradient_genotype,
                                  float*  pFloatAccumulator
                                 )
{
  float energy = 0.0f;
  float energyTmp;
  float inRangeDist = cData.dockpars.collDistance * 0.5;
#if defined (DEBUG_ENERGY_KERNEL)
  interE = 0.0f;
  intraE = 0.0f;
#endif

  // Initializing gradients (forces)
  // Derived from autodockdev/maps.py
  for (uint32_t atom_id = threadIdx.x;
                atom_id < cData.dockpars.num_of_atoms; // makes sure that gradient sum reductions give correct results if dockpars_num_atoms < NUM_OF_THREADS_PER_BLOCK
                atom_id+= blockDim.x)
  {
    // Initialize coordinates
    calc_coords[atom_id].x = cData.pKerconst_conform_ref_coords_const[3*atom_id];
    calc_coords[atom_id].y = cData.pKerconst_conform_ref_coords_const[3*atom_id+1];
    calc_coords[atom_id].z = cData.pKerconst_conform_ref_coords_const[3*atom_id+2];

    // Intermolecular gradients
    gradient[atom_id].x = 0;
    gradient[atom_id].y = 0;
    gradient[atom_id].z = 0;
  }

  // Initializing gradient genotypes
  for (uint32_t gene_cnt = threadIdx.x;
                gene_cnt < cData.dockpars.num_of_genes;
                gene_cnt+= blockDim.x)
  {
    fgradient_genotype[gene_cnt] = 0;
  }

  // General rotation moving vector
  float4 genrot_movingvec;
  genrot_movingvec.x = genotype[0];
  genrot_movingvec.y = genotype[1];
  genrot_movingvec.z = genotype[2];
  genrot_movingvec.w = 0.0f;

  // Convert orientation genes from sex. to radians
  float phi         = genotype[3] * DEG_TO_RAD;
  float theta       = genotype[4] * DEG_TO_RAD;
  float genrotangle = genotype[5] * DEG_TO_RAD;

  float4 genrot_unitvec;
  float sin_angle = sin(theta);
  float s2 = sin(genrotangle*0.5f);
  genrot_unitvec.x = s2*sin_angle*cos(phi);
  genrot_unitvec.y = s2*sin_angle*sin(phi);
  genrot_unitvec.z = s2*cos(theta);
  genrot_unitvec.w = cos(genrotangle*0.5f);
  float is_theta_gt_pi = 1.0f-2.0f*(float)(sin_angle < 0.0f);

  // uint32_t  g1 = cData.dockpars.gridsize_x;
  // uint32_t  g2 = cData.dockpars.gridsize_x_times_y;
  // uint32_t  g3 = cData.dockpars.gridsize_x_times_y_times_z;

  __syncthreads();

  // ================================================
  // CALCULATING ATOMIC POSITIONS AFTER ROTATIONS
  // ================================================
  uint32_t rotbond_id;
  float rotation_angle;
  float s;
  for (uint32_t rotation_counter = threadIdx.x;
                rotation_counter < cData.dockpars.rotbondlist_length;
                rotation_counter+=blockDim.x)
  {
    int rotation_list_element = cData.pKerconst_rotlist_rotlist_const[rotation_counter];

    if ((rotation_list_element & RLIST_DUMMY_MASK) == 0) // If not dummy rotation
    {
      uint32_t atom_id = cData.pKerconst_rotlist_rotlist_atomId_const[rotation_counter];

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
        rotbond_id = cData.pKerconst_rotlist_rotlist_rotatableBondId_const[rotation_counter];

        rotation_angle = genotype[6+rotbond_id]*DEG_TO_RAD*0.5f;
        s = sin(rotation_angle);
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
      float4 qt = quaternion_rotate(atom_to_rotate,rotation_unitvec);
      calc_coords[atom_id].x = qt.x + rotation_movingvec.x;
      calc_coords[atom_id].y = qt.y + rotation_movingvec.y;
      calc_coords[atom_id].z = qt.z + rotation_movingvec.z;
    } // End if-statement not dummy rotation
      __syncthreads();
  } // End rotation_counter for-loop

  // Inter- and intra-molecular energy calculation
  // are independent from each other, so NO barrier is needed here.
  // As these two require different operations,
  // they can be executed only sequentially on the GPU.
  float delta_distance = 0.5f * cData.dockpars.smooth;
  float smoothed_distance;

  // ================================================
  // CALCULATING INTRAMOLECULAR GRADIENTS
  // ================================================
#ifdef REPRO
  // Simplest way to ensure random order of atomic addition doesn't make answers irreproducible: use only 1 thread
  if (threadIdx.x==0) for (uint32_t contributor_counter = 0; contributor_counter < cData.dockpars.num_of_intraE_contributors; contributor_counter+= 1) {
#else
  for (uint32_t contributor_counter = threadIdx.x;
                contributor_counter < cData.dockpars.num_of_intraE_contributors;
                contributor_counter+= blockDim.x) {
#endif
    // Storing in a private variable
    // the gradient contribution of each contributing atomic pair
    float priv_gradient_per_intracontributor= 0.0f;

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
      float vbond = G;
      energyTmp = vbond * atomic_distance;
      priv_gradient_per_intracontributor += energyTmp;
      // energyTmp *= atomic_distance;
      // energy += energyTmp;
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
        float rmn = __powf(smoothed_distance,m-n);
        float rm = __powf(smoothed_distance,-m);
        energy += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -rmn*cData.pKerconst_intra_VWpars_BD_const[idx])*rm;
        priv_gradient_per_intracontributor += (n*cData.pKerconst_intra_VWpars_BD_const[idx]*rmn
                                              -m*cData.pKerconst_intra_VWpars_AC_const[idx])*rm
                                              /smoothed_distance;
        #if defined (DEBUG_ENERGY_KERNEL)
        intraE += (cData.pKerconst_intra_VWpars_AC_const[idx]
                   -rmn*cData.pKerconst_intra_VWpars_BD_const[idx])*rm
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
      
        // http://www.wolframalpha.com/input/?i=1%2F(x*(A%2B(B%2F(1%2BK*exp(-h*B*x)))))
/*        float exp_el_DIEL_K = exp_el + DIEL_K;
        float upper = DIEL_A * exp_el_DIEL_K*exp_el_DIEL_K +
                      DIEL_B * exp_el * (DIEL_B_TIMES_H_TIMES_K*atomic_distance + exp_el_DIEL_K);
        float lower = atomic_distance * (DIEL_A * exp_el_DIEL_K + DIEL_B * exp_el);
        lower *= lower;*/
      
//        priv_gradient_per_intracontributor +=  -dockpars_coeff_elec * q1 * q2 * native_divide (upper, lower) -
//                                               0.0771605f * atomic_distance * desolv_energy;
        priv_gradient_per_intracontributor += -(es_energy / atomic_distance) * diel
                                              #ifndef DIEL_FIT_ABC
                                              -es_energy * 2.21718f / (dist2*dist_shift)
                                              #else
                                              -es_energy * ((2.808f / (dist2*dist_shift)) + (0.288f / (disth4*disth_shift)))
                                              #endif
                                              -0.0771605f * atomic_distance * desolv_energy; // 1/3.6^2 = 1/12.96 = 0.0771605
      } // if cuttoff2 - internuclear-distance at 20.48A
    }

    // Decomposing "priv_gradient_per_intracontributor"
    // into the contribution of each atom of the pair.
    // Distances in Angstroms of vector that goes from
    // "atom1_id"-to-"atom2_id", therefore - subx, - suby, and - subz are used
    float grad_div_dist = -priv_gradient_per_intracontributor / dist;
#ifdef FLOAT_GRADIENTS
    float priv_intra_gradient_x = subx * grad_div_dist;
    float priv_intra_gradient_y = suby * grad_div_dist;
    float priv_intra_gradient_z = subz * grad_div_dist;
#else
    int priv_intra_gradient_x = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * subx * grad_div_dist)));
    int priv_intra_gradient_y = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * suby * grad_div_dist)));
    int priv_intra_gradient_z = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * subz * grad_div_dist)));
#endif

    // Calculating gradients in xyz components.
    // Gradients for both atoms in a single contributor pair
    // have the same magnitude, but opposite directions
#ifdef FLOAT_GRADIENTS
    ATOMICSUBF32(&gradient[atom1_id].x, priv_intra_gradient_x);
    ATOMICSUBF32(&gradient[atom1_id].y, priv_intra_gradient_y);
    ATOMICSUBF32(&gradient[atom1_id].z, priv_intra_gradient_z);

    ATOMICADDF32(&gradient[atom2_id].x, priv_intra_gradient_x);
    ATOMICADDF32(&gradient[atom2_id].y, priv_intra_gradient_y);
    ATOMICADDF32(&gradient[atom2_id].z, priv_intra_gradient_z);
#else
    ATOMICSUBI32(&gradient[atom1_id].x, priv_intra_gradient_x);
    ATOMICSUBI32(&gradient[atom1_id].y, priv_intra_gradient_y);
    ATOMICSUBI32(&gradient[atom1_id].z, priv_intra_gradient_z);

    ATOMICADDI32(&gradient[atom2_id].x, priv_intra_gradient_x);
    ATOMICADDI32(&gradient[atom2_id].y, priv_intra_gradient_y);
    ATOMICADDI32(&gradient[atom2_id].z, priv_intra_gradient_z);

#endif
  } // End contributor_counter for-loop (INTRAMOLECULAR ENERGY)
  __syncthreads();

  // Transform gradients_inter_{x|y|z}
  // into local_gradients[i] (with four quaternion genes)
  // Derived from autodockdev/motions.py/forces_to_delta_genes()

  // Transform local_gradients[i] (with four quaternion genes)
  // into local_gradients[i] (with three Shoemake genes)
  // Derived from autodockdev/motions.py/_get_cube3_gradient()
  // ------------------------------------------

  // start by populating "gradient_intra_*" with torque values
  float4 torque_rot;
  torque_rot.x = 0.0f;
  torque_rot.y = 0.0f;
  torque_rot.z = 0.0f;
  float gx = 0.0f;
  float gy = 0.0f;
  float gz = 0.0f;
  // overall rotation is only for the moving ligand
  for (uint32_t atom_cnt = threadIdx.x;
                atom_cnt < cData.dockpars.true_ligand_atoms;
                atom_cnt+= blockDim.x) {
    float3 r;
    r.x = (calc_coords[atom_cnt].x - genrot_movingvec.x) * cData.dockpars.grid_spacing;
    r.y = (calc_coords[atom_cnt].y - genrot_movingvec.y) * cData.dockpars.grid_spacing;
    r.z = (calc_coords[atom_cnt].z - genrot_movingvec.z) * cData.dockpars.grid_spacing;

    // Re-using "gradient_inter_*" for total gradient (inter+intra)
    float3 force;
#ifdef FLOAT_GRADIENTS
    force.x = gradient[atom_cnt].x;
    force.y = gradient[atom_cnt].y;
    force.z = gradient[atom_cnt].z;
#else
    force.x = ONEOVERTERMSCALE * (float)gradient[atom_cnt].x;
    force.y = ONEOVERTERMSCALE * (float)gradient[atom_cnt].y;
    force.z = ONEOVERTERMSCALE * (float)gradient[atom_cnt].z;
#endif
    gx += force.x;
    gy += force.y;
    gz += force.z;
    float4 tr = cross(r, force);
    torque_rot.x += tr.x;
    torque_rot.y += tr.y;
    torque_rot.z += tr.z;
  }

  // Do a reduction over the total gradient containing prepared "gradient_intra_*" values
  REDUCEFLOATSUM(torque_rot.x, pFloatAccumulator);
  REDUCEFLOATSUM(torque_rot.y, pFloatAccumulator);
  REDUCEFLOATSUM(torque_rot.z, pFloatAccumulator);
  // TODO
  // -------------------------------------------------------
  // Obtaining energy and translation-related gradients
  // -------------------------------------------------------
  // reduction over partial energies and prepared "gradient_intra_*" values
  REDUCEFLOATSUM(energy, pFloatAccumulator);
#if defined (DEBUG_ENERGY_KERNEL)
  REDUCEFLOATSUM(intraE, pFloatAccumulator);
#endif
  REDUCEFLOATSUM(gx, pFloatAccumulator);
  REDUCEFLOATSUM(gy, pFloatAccumulator);
  REDUCEFLOATSUM(gz, pFloatAccumulator);
  global_energy = energy;
#ifndef FLOAT_GRADIENTS
  int* gradient_genotype = (int*)fgradient_genotype;
#endif
  if (threadIdx.x == 0) {
    // Scaling gradient for translational genes as
    // their corresponding gradients were calculated in the space
    // where these genes are in Angstrom,
    // but OmegaDocker translational genes are within in grids
#ifdef FLOAT_GRADIENTS
    fgradient_genotype[0] = gx * cData.dockpars.grid_spacing;
    fgradient_genotype[1] = gy * cData.dockpars.grid_spacing;
    fgradient_genotype[2] = gz * cData.dockpars.grid_spacing;

    #if defined (PRINT_GRAD_TRANSLATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("gradient_x:%f\n", fgradient_genotype [0]);
    printf("gradient_y:%f\n", fgradient_genotype [1]);
    printf("gradient_z:%f\n", fgradient_genotype [2]);
    #endif
#else
    gradient_genotype[0] = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * gx * cData.dockpars.grid_spacing)));
    gradient_genotype[1] = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * gy * cData.dockpars.grid_spacing)));
    gradient_genotype[2] = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * gz * cData.dockpars.grid_spacing)));

    #if defined (PRINT_GRAD_TRANSLATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("gradient_x:%f\n", gradient_genotype [0]);
    printf("gradient_y:%f\n", gradient_genotype [1]);
    printf("gradient_z:%f\n", gradient_genotype [2]);
    #endif
#endif
  }
  __syncthreads();

  // ------------------------------------------
  // Obtaining rotation-related gradients
  // ------------------------------------------
  if (threadIdx.x == 0) {
    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-20s %-10.6f %-10.6f %-10.6f\n", "final torque: ", torque_rot.x, torque_rot.y, torque_rot.z);
    #endif

    // Derived from rotation.py/axisangle_to_q()
    // genes[3:7] = rotation.axisangle_to_q(torque, rad)
    float torque_length = norm3df(torque_rot.x, torque_rot.y, torque_rot.z);
    torque_length += (torque_length<1e-20f)*1e-20f;

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-20s %-10.6f\n", "torque length: ", torque_length);
    #endif

    /*
    // Infinitesimal rotation in radians
    const float infinitesimal_radian = 1E-5;
    */

    // Finding the quaternion that performs
    // the infinitesimal rotation around torque axis
    float4 quat_torque;
    quat_torque.x = torque_rot.x * SIN_HALF_INFINITESIMAL_RADIAN / torque_length;
    quat_torque.y = torque_rot.y * SIN_HALF_INFINITESIMAL_RADIAN / torque_length;
    quat_torque.z = torque_rot.z * SIN_HALF_INFINITESIMAL_RADIAN / torque_length;
    quat_torque.w = COS_HALF_INFINITESIMAL_RADIAN;

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-20s %-10.6f\n", "INFINITESIMAL_RADIAN: ", INFINITESIMAL_RADIAN);
    printf("%-20s %-10.6f %-10.6f %-10.6f %-10.6f\n", "quat_torque (w,x,y,z): ", quat_torque.w, quat_torque.x, quat_torque.y, quat_torque.z);
    #endif

    // Converting quaternion gradients into orientation gradients
    // Derived from autodockdev/motion.py/_get_cube3_gradient
    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s %-10.6f %-10.6f %-10.6f %-10.6f\n", "current_q (w,x,y,z): ", genrot_unitvec.w, genrot_unitvec.x, genrot_unitvec.y, genrot_unitvec.z);
    #endif

    // This is where we want to be in quaternion space
    // target_q = rotation.q_mult(q, current_q)
    float4 target_q = quaternion_multiply(quat_torque, genrot_unitvec);

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s %-10.6f %-10.6f %-10.6f %-10.6f\n", "target_q (w,x,y,z): ", target_q.w, target_q.x, target_q.y, target_q.z);
    #endif

    // This is where we are in the orientation axis-angle space
    // Equivalent to "current_oclacube" in autodockdev/motions.py
    float current_phi      = fmod_pi2(PI_TIMES_2 + phi);
    float current_theta    = fmod_pi2(PI_TIMES_2 + theta);
    float current_rotangle = fmod_pi2(PI_TIMES_2 + genrotangle);

    // This is where we want to be in the orientation axis-angle space
    float target_phi, target_theta, target_rotangle;

    // target_oclacube = quaternion_to_oclacube(target_q, theta_larger_than_pi)
    // Derived from autodockdev/motions.py/quaternion_to_oclacube()
    // In our terms means quaternion_to_oclacube(target_q{w|x|y|z}, theta_larger_than_pi)
    target_rotangle = 2.0f * fast_acos(target_q.w); // = 2.0f * ang;
    float sin_ang = sqrt(1.0f-target_q.w*target_q.w); // = native_sin(ang);

    target_theta = PI_TIMES_2 + is_theta_gt_pi * fast_acos(target_q.z / sin_ang );
    target_phi   = fmod_pi2((atan2( is_theta_gt_pi*target_q.y, is_theta_gt_pi*target_q.x) + PI_TIMES_2));

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s %-10.6f %-10.6f %-10.6f\n", "target_axisangle (1,2,3): ", target_phi, target_theta, target_rotangle);
    #endif

    // The infinitesimal rotation will produce an infinitesimal displacement
    // in shoemake space. This is to guarantee that the direction of
    // the displacement in shoemake space is not distorted.
    // The correct amount of displacement in shoemake space is obtained
    // by multiplying the infinitesimal displacement by shoemake_scaling:
    //float shoemake_scaling = native_divide(torque_length, INFINITESIMAL_RADIAN/*infinitesimal_radian*/);
    float orientation_scaling = torque_length * INV_INFINITESIMAL_RADIAN;

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s %-10.6f\n", "orientation_scaling: ", orientation_scaling);
    #endif

    // Derivates in cube3
    float grad_phi, grad_theta, grad_rotangle;
    grad_phi      = orientation_scaling * (fmod_pi2(target_phi      - current_phi      + PI_FLOAT) - PI_FLOAT);
    grad_theta    = orientation_scaling * (fmod_pi2(target_theta    - current_theta    + PI_FLOAT) - PI_FLOAT);
    grad_rotangle = orientation_scaling * (fmod_pi2(target_rotangle - current_rotangle + PI_FLOAT) - PI_FLOAT);

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s \n", "grad_axisangle (1,2,3) - before empirical scaling: ");
    printf("%-13s %-13s %-13s \n", "grad_phi", "grad_theta", "grad_rotangle");
    printf("%-13.6f %-13.6f %-13.6f\n", grad_phi, grad_theta, grad_rotangle);
    #endif

    // Correcting theta gradients interpolating
    // values from correction look-up-tables
    // (X0,Y0) and (X1,Y1) are known points
    // How to find the Y value in the straight line between Y0 and Y1,
    // corresponding to a certain X?
    /*
      | dependence_on_theta_const
      | dependence_on_rotangle_const
      |
      |
      |                        Y1
      |
      |             Y=?
      |    Y0
      |_________________________________ angle_const
           X0         X        X1
    */

    // Finding the index-position of "grad_delta" in the "angle_const" array
    //uint index_theta    = floor(native_divide(current_theta    - angle_const[0], angle_delta));
    //uint index_rotangle = floor(native_divide(current_rotangle - angle_const[0], angle_delta));
    uint index_theta    = floor((current_theta    - cData.pMem_angle_const[0]) * inv_angle_delta);
    uint index_rotangle = floor((current_rotangle - cData.pMem_angle_const[0]) * inv_angle_delta);

    // Interpolating theta values
    // X0 -> index - 1
    // X1 -> index + 1
    // Expresed as weighted average:
    // Y = [Y0 * ((X1 - X) / (X1-X0))] +  [Y1 * ((X - X0) / (X1-X0))]
    // Simplified for GPU (less terms):
    // Y = [Y0 * (X1 - X) + Y1 * (X - X0)] / (X1 - X0)
    // Taking advantage of constant:
    // Y = [Y0 * (X1 - X) + Y1 * (X - X0)] * inv_angle_delta

    float X0, Y0;
    float X1, Y1;
    float dependence_on_theta;    //Y = dependence_on_theta

    // Using interpolation on out-of-bounds elements results in hang
    if ((index_theta <= 0) || (index_theta >= 999))
    {
      dependence_on_theta = cData.pMem_dependence_on_theta_const[stick_to_bounds(index_theta,0,999)];
    } else
    {
      X0 = cData.pMem_angle_const[index_theta];
      X1 = cData.pMem_angle_const[index_theta+1];
      Y0 = cData.pMem_dependence_on_theta_const[index_theta];
      Y1 = cData.pMem_dependence_on_theta_const[index_theta+1];
      dependence_on_theta = (Y0 * (X1-current_theta) + Y1 * (current_theta-X0)) * inv_angle_delta;
    }

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s %-10.6f\n", "dependence_on_theta: ", dependence_on_theta);
    #endif

    // Interpolating rotangle values
    float dependence_on_rotangle;   //Y = dependence_on_rotangle
    // Using interpolation on previous and/or next elements results in hang
    // Using interpolation on out-of-bounds elements results in hang
    if ((index_rotangle <= 0) || (index_rotangle >= 999))
    {
      dependence_on_rotangle = cData.pMem_dependence_on_rotangle_const[stick_to_bounds(index_rotangle,0,999)];
    } else
    {
      X0 = cData.pMem_angle_const[index_rotangle];
      X1 = cData.pMem_angle_const[index_rotangle+1];
      Y0 = cData.pMem_dependence_on_rotangle_const[index_rotangle];
      Y1 = cData.pMem_dependence_on_rotangle_const[index_rotangle+1];
      dependence_on_rotangle = (Y0 * (X1-current_rotangle) + Y1 * (current_rotangle-X0)) * inv_angle_delta;
    }

    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s %-10.6f\n", "dependence_on_rotangle: ", dependence_on_rotangle);
    #endif

    // Setting gradient rotation-related genotypes in cube
    // Multiplicating by DEG_TO_RAD is to make it uniform to DEG (see torsion gradients)
#ifdef FLOAT_GRADIENTS
    fgradient_genotype[3] = (grad_phi / (dependence_on_theta * dependence_on_rotangle)) * DEG_TO_RAD;
    fgradient_genotype[4] = (grad_theta / dependence_on_rotangle) * DEG_TO_RAD;
    fgradient_genotype[5] = grad_rotangle * DEG_TO_RAD;
    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s \n", "grad_axisangle (1,2,3) - after empirical scaling: ");
    printf("%-13s %-13s %-13s \n", "grad_phi", "grad_theta", "grad_rotangle");
    printf("%-13.6f %-13.6f %-13.6f\n", fgradient_genotype[3], fgradient_genotype[4], fgradient_genotype[5]);
    #endif
#else
    gradient_genotype[3] = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * (grad_phi / (dependence_on_theta * dependence_on_rotangle)) * DEG_TO_RAD)));
    gradient_genotype[4] = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * (grad_theta / dependence_on_rotangle) * DEG_TO_RAD)));
    gradient_genotype[5] = lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * grad_rotangle * DEG_TO_RAD)));
    #if defined (PRINT_GRAD_ROTATION_GENES)
    printf("\n%s\n", "----------------------------------------------------------");
    printf("%-30s \n", "grad_axisangle (1,2,3) - after empirical scaling: ");
    printf("%-13s %-13s %-13s \n", "grad_phi", "grad_theta", "grad_rotangle");
    printf("%-13.6f %-13.6f %-13.6f\n", gradient_genotype[3], gradient_genotype[4], gradient_genotype[5]);
    #endif
#endif
  }
  __syncthreads();

  // ------------------------------------------
  // Obtaining torsion-related gradients
  // ------------------------------------------
  uint32_t num_torsion_genes = cData.dockpars.num_of_genes-6;
  for (uint32_t idx = threadIdx.x; idx < num_torsion_genes * cData.dockpars.num_of_atoms; idx += blockDim.x) {
    uint32_t rotable_atom_cnt = idx / num_torsion_genes;
    rotbond_id = idx - rotable_atom_cnt * num_torsion_genes; // this is a bit cheaper than % (modulo)

    if (rotable_atom_cnt >= cData.pKerconst_grads_num_rotating_atoms_per_rotbond[rotbond_id])
      continue; // Nothing to do

    // Querying ids of atoms belonging to the rotatable bond in question
    int atom1_id = cData.pKerconst_grads_rotbonds[2*rotbond_id];
    int atom2_id = cData.pKerconst_grads_rotbonds[2*rotbond_id+1];

    float3 atomRef_coords;
    atomRef_coords.x = calc_coords[atom1_id].x;
    atomRef_coords.y = calc_coords[atom1_id].y;
    atomRef_coords.z = calc_coords[atom1_id].z;
    float3 rotation_unitvec;

    rotation_unitvec.x = calc_coords[atom2_id].x - atomRef_coords.x;
    rotation_unitvec.y = calc_coords[atom2_id].y - atomRef_coords.y;
    rotation_unitvec.z = calc_coords[atom2_id].z - atomRef_coords.z;
    float l = rnorm3df(rotation_unitvec.x, rotation_unitvec.y, rotation_unitvec.z);
    rotation_unitvec.x *= l;
    rotation_unitvec.y *= l;
    rotation_unitvec.z *= l;

    // Torque of torsions
    uint lig_atom_id;
    if(cData.dockpars.num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
    {
      lig_atom_id = cData.pKerconst_grads_rotbonds_atoms[NUM_OF_THREADS_PER_BLOCK * rotbond_id + rotable_atom_cnt];
    }
    else
    {
      lig_atom_id = cData.pKerconst_grads_rotbonds_atoms[cData.dockpars.num_of_atoms * rotbond_id + rotable_atom_cnt];
    }
    float4 torque_tor;
    float3 r, atom_force;

    // Calculating torque on point "A"
    // They are converted back to Angstroms here
    r.x = (calc_coords[lig_atom_id].x - atomRef_coords.x);
    r.y = (calc_coords[lig_atom_id].y - atomRef_coords.y);
    r.z = (calc_coords[lig_atom_id].z - atomRef_coords.z);

    // Re-using "gradient_inter_*" for total gradient (inter+intra)
#ifdef FLOAT_GRADIENTS
    atom_force.x = gradient[lig_atom_id].x;
    atom_force.y = gradient[lig_atom_id].y;
    atom_force.z = gradient[lig_atom_id].z;
#else
    atom_force.x = ONEOVERTERMSCALE * gradient[lig_atom_id].x;
    atom_force.y = ONEOVERTERMSCALE * gradient[lig_atom_id].y;
    atom_force.z = ONEOVERTERMSCALE * gradient[lig_atom_id].z;
#endif
    torque_tor = cross(r, atom_force);
    float torque_on_axis = (rotation_unitvec.x * torque_tor.x  +
          rotation_unitvec.y * torque_tor.y  +
          rotation_unitvec.z * torque_tor.z) * cData.dockpars.grid_spacing;

    // Assignment of gene-based gradient
    // - this works because a * (a_1 + a_2 + ... + a_n) = a*a_1 + a*a_2 + ... + a*a_n
#ifdef FLOAT_GRADIENTS
    ATOMICADDF32(&fgradient_genotype[rotbond_id+6], torque_on_axis * DEG_TO_RAD); /*(M_PI / 180.0f)*/;
#else
    ATOMICADDI32(&gradient_genotype[rotbond_id+6], lrintf(fminf(MAXTERM, fmaxf(-MAXTERM, TERMSCALE * torque_on_axis * DEG_TO_RAD)))); /*(M_PI / 180.0f)*/;
#endif
  }
  __syncthreads();

#ifndef FLOAT_GRADIENTS
  for (uint32_t gene_cnt = threadIdx.x;
                gene_cnt < cData.dockpars.num_of_genes;
                gene_cnt+= blockDim.x) {
    fgradient_genotype[gene_cnt] = ONEOVERTERMSCALE * (float)gradient_genotype[gene_cnt];
  }
  __syncthreads();
#endif
  #if defined (CONVERT_INTO_ANGSTROM_RADIAN)
  for (uint32_t gene_cnt = threadIdx.x + 3; // Only for gene_cnt > 2 means start gene_cnt at 3
                gene_cnt < cData.dockpars.num_of_genes;
                gene_cnt+= blockDim.x)
  {
    fgradient_genotype[gene_cnt] *= cData.dockpars.grid_spacing * cData.dockpars.grid_spacing * SCFACTOR_ANGSTROM_RADIAN;
  }
  __syncthreads();
  #endif
}