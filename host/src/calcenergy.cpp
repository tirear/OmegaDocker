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


// Output interaction pairs
// #define INTERACTION_PAIR_INFO

#include "calcenergy.h"

int prepare_const_fields_for_gpu(
                                 Liganddata*                  myligand_reference,
                                 Dockpars*                    mypars,
                                 kernelconstant_interintra*   KerConst_interintra,
                                 kernelconstant_intracontrib* KerConst_intracontrib,
                                 kernelconstant_intra*        KerConst_intra,
                                 kernelconstant_rotlist*      KerConst_rotlist,
                                 kernelconstant_conform*      KerConst_conform,
                                 kernelconstant_grads*        KerConst_grads
                                )
// The function fills the constant memory field of the GPU
// based on the parameters describing ligand, flexres, and
// docking parameters as well as reference orientation angles.
{
  // Some variables
  int i, j;
  float* floatpoi;
  int *intpoi;

// --------------------------------------------------
// Allocating memory on the heap (not stack) with new
// --------------------------------------------------
// atom_charges:            Stores the ligand atom charges.
//                          Element i corresponds to atom with atom ID i in myligand_reference.
  float* atom_charges            = new float[myligand_reference->num_of_atoms];
// atom_types:              Stores the ligand atom type IDs according to myligand_reference.
//                          Element i corresponds to atom with ID i in myligand_reference.
  int*   atom_types              = new int[myligand_reference->num_of_atoms];
// intraE_contributors:     Each three contiguous items describe an intramolecular contributor.
//                          The first two elements store the atom ID of the contributors according to myligand_reference.
//                          The third element is 0, if no H-bond can occur between the two atoms, and 1, if it can.
  int max_num_of_intraE_contributors = myligand_reference->num_of_atoms * myligand_reference->num_of_atoms;
  int*   intraE_contributors     = new int[2 * max_num_of_intraE_contributors];
  size_t num_used_atom_types      = myligand_reference->num_of_atypes;
  size_t num_used_atom_type_pairs = myligand_reference->num_of_atypes * myligand_reference->num_of_atypes;
  unsigned short* VWpars_exp     = new unsigned short[num_used_atom_type_pairs];
  float* reqm_AB                 = new float[num_used_atom_type_pairs];
// VWpars_AC_const:         Stores the A or C van der Waals parameters.
//                          The element i*MAX_NUM_OF_ATYPES+j and j*MAX_NUM_OF_ATYPES+i corresponds to A or C in case of
//                          H-bond for atoms with type ID i and j (according to myligand_reference).
  float* VWpars_AC               = new float[num_used_atom_type_pairs];
// VWpars_BD:               Stores the B or D van der Waals parameters similar to VWpars_AC.
  float* VWpars_BD               = new float[num_used_atom_type_pairs];
// dspars_S:                Stores the S desolvation parameters.
//                          The element i corresponds to the S parameter of atom with type ID i
//                          according to myligand_reference.
  float* dspars_S                = new float[myligand_reference->num_of_atypes];
  float* dspars_V                = new float[myligand_reference->num_of_atypes];
  int max_num_rotations;
  int windowSize;
  if(myligand_reference->num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    windowSize = NUM_OF_THREADS_PER_BLOCK;
  }
  else
  {
    windowSize = myligand_reference->num_of_atoms;
  }
  max_num_rotations = windowSize * (myligand_reference->num_of_rotbonds + 1);
  int   *rotlist                 = new int[max_num_rotations];
  int   *rotlist_atomId          = (int *)  calloc(max_num_rotations, sizeof(int));
  int   *rotlist_rotatableBondId = (int *)  calloc(max_num_rotations, sizeof(int));
// ref_coords_x:            Stores the x coordinates of the reference ligand atoms.
//                          Element i corresponds to the x coordinate of the atom with atom ID i (according to myligand_reference).
  float* ref_coords_x            = new float[myligand_reference->num_of_atoms];
// ref_coords_y:            Stores the y coordinates of the reference ligand atoms similarly to ref_coords_x.
  float* ref_coords_y            = new float[myligand_reference->num_of_atoms];
// ref_coords_z:            Stores the z coordinates of the reference ligand atoms similarly to ref_coords_x.
  float* ref_coords_z            = new float[myligand_reference->num_of_atoms];
// rotbonds_moving_vectors: Stores the coordinates of rotatable bond moving vectors. Element i, i+1 and i+2 (where i%3=0)
//                          correspond to the moving vector coordinates x, y and z of rotbond ID i, respectively
//                          (according to myligand_reference).
  float* rotbonds_moving_vectors = new float[3*myligand_reference->num_of_rotbonds];
// rotbonds_unit_vectors:   Stores the coordinates of rotatable bond unit vectors similarly to rotbonds_moving_vectors.
  float* rotbonds_unit_vectors = new float[3*myligand_reference->num_of_rotbonds];

// rot_bonds:               Added for calculating torsion-related gradients.
//                          Passing list of rotbond-atom ids to the GPU.
//                          Contains the same information as processligand.h/Liganddata->rotbonds
//                          Each row corresponds to one rotatable bond of the ligand.
//                          The rotatable bond is described with the indices of the two atoms which are connected
//                          to each other by the bond.
//                          The row index is equal to the index of the rotatable bond.
  int*   rotbonds                = new int[2*myligand_reference->num_of_rotbonds];

// rotbonds_atoms:          Contains the same information as processligand.h/Liganddata->atom_rotbonds:
//                                  Matrix that contains the rotatable bonds - atoms assignment.
//                                  If the element [atom index][rotatable bond index] is equal to 1,
//                                  then the atom must be rotated if the bond rotates. A 0 means the opposite.
  int*   rotbonds_atoms          = new int[max_num_rotations];

// num_rotating_atoms_per_rotbond:
//                          Each entry corresponds to a rotbond_id
//                          The value of an entry indicates the number of atoms that rotate along with that rotbond_id
  int*   num_rotating_atoms_per_rotbond = new int[myligand_reference->num_of_rotbonds];
// --------------------------------------------------

  // charges and type id-s
  floatpoi = atom_charges;
  intpoi = atom_types;

  for (i=0; i < myligand_reference->num_of_atoms; i++)
  {
    *floatpoi = (float) myligand_reference->atom_idxyzq[i][4];
    *intpoi  = (int) myligand_reference->atom_idxyzq[i][0];
    floatpoi++;
    intpoi++;
  }

  // intramolecular energy contributors
  myligand_reference->num_of_intraE_contributors = 0;
  for (i=0; i<myligand_reference->num_of_atoms-1; i++)
    for (j=i+1; j<myligand_reference->num_of_atoms; j++)
    {
      if (myligand_reference->intraE_contributors[i][j]){
#ifdef INTERACTION_PAIR_INFO
        printf("Pair interaction between: %i <-> %i\n",i+1,j+1);
#endif
        myligand_reference->num_of_intraE_contributors++;
      }
    }

  if (myligand_reference->num_of_intraE_contributors > max_num_of_intraE_contributors)
  {
    printf("Error: Number of intramolecular energy contributor is larger than maximum (%d).\n",max_num_of_intraE_contributors);
    fflush(stdout);
    return 1;
  }

  intpoi = intraE_contributors;
  for (i=0; i<myligand_reference->num_of_atoms-1; i++)
    for (j=i+1; j<myligand_reference->num_of_atoms; j++)
    {
      if (myligand_reference->intraE_contributors[i][j] == 1)
      {
        *intpoi = (int) i;
        intpoi++;
        *intpoi = (int) j;
        intpoi++;
      }
    }

  // van der Waals parameters
  for (i=0; i<myligand_reference->num_of_atypes; i++)
    for (j=0; j<myligand_reference->num_of_atypes; j++)
    {
      *(VWpars_exp + i*myligand_reference->num_of_atypes + j) = myligand_reference->VWpars_exp[i][j];
      floatpoi = reqm_AB + i*myligand_reference->num_of_atypes + j;
      *floatpoi = (float) myligand_reference->reqm_AB[i][j];

      if (is_H_bond(myligand_reference->base_atom_types[i], myligand_reference->base_atom_types[j]) &&
          (!is_mod_pair(myligand_reference->atom_types[i], myligand_reference->atom_types[j], mypars->nr_mod_atype_pairs, mypars->mod_atype_pairs)))
      {
        floatpoi = VWpars_AC + i*myligand_reference->num_of_atypes + j;
        *floatpoi = (float) myligand_reference->VWpars_C[i][j];
        floatpoi = VWpars_AC + j*myligand_reference->num_of_atypes + i;
        *floatpoi = (float) myligand_reference->VWpars_C[j][i];

        floatpoi = VWpars_BD + i*myligand_reference->num_of_atypes + j;
        *floatpoi = (float) myligand_reference->VWpars_D[i][j];
        floatpoi = VWpars_BD + j*myligand_reference->num_of_atypes + i;
        *floatpoi = (float) myligand_reference->VWpars_D[j][i];
      }
      else
      {
        floatpoi = VWpars_AC + i*myligand_reference->num_of_atypes + j;
        *floatpoi = (float) myligand_reference->VWpars_A[i][j];
        floatpoi = VWpars_AC + j*myligand_reference->num_of_atypes + i;
        *floatpoi = (float) myligand_reference->VWpars_A[j][i];

        floatpoi = VWpars_BD + i*myligand_reference->num_of_atypes + j;
        *floatpoi = (float) myligand_reference->VWpars_B[i][j];
        floatpoi = VWpars_BD + j*myligand_reference->num_of_atypes + i;
        *floatpoi = (float) myligand_reference->VWpars_B[j][i];
      }
    }

  // desolvation parameters
  for (i=0; i<myligand_reference->num_of_atypes; i++)
  {
    dspars_S[i] = myligand_reference->solpar[i];
    dspars_V[i] = myligand_reference->volume[i];
  }

  // generate rotation list
  if(gen_rotlist(myligand_reference,
                 rotlist,
                 rotlist_atomId,
                 rotlist_rotatableBondId) != 0)
  {
    printf("Error: Number of required rotations is larger than maximum (%d).\n",max_num_rotations);
    return 1;
  }

  // coordinates of reference ligand
  for (i=0; i < myligand_reference->num_of_atoms; i++)
  {
    ref_coords_x[i] = myligand_reference->atom_idxyzq[i][1];
    ref_coords_y[i] = myligand_reference->atom_idxyzq[i][2];
    ref_coords_z[i] = myligand_reference->atom_idxyzq[i][3];
  }

  // rotatable bond vectors
  for (i=0; i < myligand_reference->num_of_rotbonds; i++){
    for (j=0; j<3; j++)
    {
      rotbonds_moving_vectors[3*i+j] = myligand_reference->rotbonds_moving_vectors[i][j];
      rotbonds_unit_vectors[3*i+j] = myligand_reference->rotbonds_unit_vectors[i][j];
    }
  }

  // Added for calculating torsion-related gradients.
  // Passing list of rotbond-atoms ids to the GPU.
  // Contains the same information as processligand.h/Liganddata->rotbonds
  for (i=0; i < myligand_reference->num_of_rotbonds; i++)
  {
    rotbonds [2*i]   = myligand_reference->rotbonds[i][0]; // id of first-atom
    rotbonds [2*i+1] = myligand_reference->rotbonds[i][1]; // id of second atom
  }

  // Contains the same information as processligand.h/Liganddata->atom_rotbonds
  // "atom_rotbonds": array that contains the rotatable bonds - atoms assignment.
  // If the element atom_rotbonds[atom index][rotatable bond index] is equal to 1,
  // it means,that the atom must be rotated if the bond rotates. A 0 means the opposite.
  for (i=0; i<myligand_reference->num_of_rotbonds; i++)
  {
    num_rotating_atoms_per_rotbond [i] = 0;
  }

  for (i=0; i < myligand_reference->num_of_rotbonds; i++)
  {
    // Pointing to the mem area corresponding to a given rotbond
    intpoi = rotbonds_atoms + windowSize*i;

    for (j=0; j < myligand_reference->num_of_atoms; j++)
    {
      /*
      rotbonds_atoms [MAX_NUM_OF_ATOMS*i+j] = myligand_reference->atom_rotbonds [j][i];
      */

      // If an atom rotates with a rotbond, then
      // add its atom-id to the entry corresponding to the rotbond-id.
      // Also, count the number of atoms that rotate with a certain rotbond
      if (myligand_reference->atom_rotbonds [j][i] == 1){
        *intpoi = j;
        intpoi++;
        num_rotating_atoms_per_rotbond [i] ++;
      }

    }
  }

  int m;
  size_t m2;
  
  for (m=0;m<myligand_reference->num_of_atoms;m++){
    KerConst_interintra->ignore_inter_const[m] = (char)myligand_reference->ignore_inter[m];
  }
  for (m=0;m<myligand_reference->num_of_atoms;m++){ KerConst_interintra->atom_charges_const[m]    = atom_charges[m];    }
  for (m=0;m<myligand_reference->num_of_atoms;m++){ KerConst_interintra->atom_types_const[m]      = atom_types[m];      }
  for (m=0;m<myligand_reference->num_of_atoms;m++){ KerConst_interintra->atom_types_map_const[m]  = myligand_reference->atom_map_to_fgrids[m]; }

  for (m=0;m<2*max_num_of_intraE_contributors;m++){ KerConst_intracontrib->intraE_contributors_const[m] = intraE_contributors[m]; }

  for (m2=0;m2<num_used_atom_types;m2++)      { KerConst_intra->atom_types_reqm_const[m2] = myligand_reference->atom_types_reqm[m2]; }
  for (m2=0;m2<num_used_atom_type_pairs;m2++) { KerConst_intra->VWpars_exp_const[m2]      = VWpars_exp[m2]; }
  for (m2=0;m2<num_used_atom_type_pairs;m2++) { KerConst_intra->reqm_AB_const[m2]         = reqm_AB[m2]; }
  for (m2=0;m2<num_used_atom_type_pairs;m2++) { KerConst_intra->VWpars_AC_const[m2]       = VWpars_AC[m2]; }
  for (m2=0;m2<num_used_atom_type_pairs;m2++) { KerConst_intra->VWpars_BD_const[m2]       = VWpars_BD[m2]; }
  for (m2=0;m2<num_used_atom_types;m2++)      { KerConst_intra->dspars_S_const[m2]        = dspars_S[m2]; }
  for (m2=0;m2<num_used_atom_types;m2++)      { KerConst_intra->dspars_V_const[m2]        = dspars_V[m2]; }

  for (m=0;m<max_num_rotations;m++) {
    KerConst_rotlist->rotlist_const[m]  = rotlist[m];
    KerConst_rotlist->rotlist_atomId[m]          = rotlist_atomId[m];
    KerConst_rotlist->rotlist_rotatableBondId[m] = rotlist_rotatableBondId[m];
  }

  for (m=0;m<myligand_reference->num_of_atoms;m++) {
    KerConst_conform->ref_coords_const[3*m]     = ref_coords_x[m];
    KerConst_conform->ref_coords_const[3*m+1]   = ref_coords_y[m];
    KerConst_conform->ref_coords_const[3*m+2]   = ref_coords_z[m];
  }

  for (m=0;m<3*myligand_reference->num_of_rotbonds;m++)
  {
    KerConst_conform->rotbonds_moving_vectors_const[m] = rotbonds_moving_vectors[m];
    KerConst_conform->rotbonds_unit_vectors_const[m]   = rotbonds_unit_vectors[m];
  }

  // Added for calculating torsion-related gradients.
  // Passing list of rotbond-atoms ids to the GPU.
  // Contains the same information as processligand.h/Liganddata->rotbonds
  for (m=0;m<2*myligand_reference->num_of_rotbonds;m++)       { KerConst_grads->rotbonds[m]           = rotbonds[m]; }
  for (m=0;m<max_num_rotations;m++)   { KerConst_grads->rotbonds_atoms[m]                 = rotbonds_atoms[m]; }
  for (m=0;m<myligand_reference->num_of_rotbonds;m++)       { KerConst_grads->num_rotating_atoms_per_rotbond[m] = num_rotating_atoms_per_rotbond[m]; }

  delete[] atom_charges;
  delete[] atom_types;
  delete[] intraE_contributors;
  delete[] VWpars_exp;
  delete[] reqm_AB;
  delete[] VWpars_AC;
  delete[] VWpars_BD;
  delete[] dspars_S;
  delete[] dspars_V;
  delete[] rotlist;
  delete[] ref_coords_x;
  delete[] ref_coords_y;
  delete[] ref_coords_z;
  delete[] rotbonds_moving_vectors;
  delete[] rotbonds_unit_vectors;
  delete[] rotbonds;
  delete[] rotbonds_atoms;
  delete[] num_rotating_atoms_per_rotbond;
  free(rotlist_atomId);
  free(rotlist_rotatableBondId);

  return 0;
}



void make_reqrot_ordering(
                          int *number_of_req_rotations,
                          int *atom_id_of_numrots,
                          int num_of_atoms
                         )
// The function puts the first array into a descending order and
// performs the same operations on the second array (since element i of
// number_or_req_rotations and element i of atom_id_of_numrots correspond to each other).
// Element i of the former array stores how many rotations have to be perfomed on the atom
// whose atom ID is stored by element i of the latter one. The third parameter has to be equal
// to the number of ligand atoms
{
  int i, j;
  int temp;

  for (j=0; j<num_of_atoms-1; j++)
    for (i=num_of_atoms-2; i>=j; i--)
      if (number_of_req_rotations[i+1] > number_of_req_rotations[i])
      {
        temp = number_of_req_rotations[i];
        number_of_req_rotations[i] = number_of_req_rotations[i+1];
        number_of_req_rotations[i+1] = temp;

        temp = atom_id_of_numrots[i];
        atom_id_of_numrots[i] = atom_id_of_numrots[i+1];
        atom_id_of_numrots[i+1] = temp;
      }
}



int gen_rotlist(
                Liganddata* myligand,
                int        *rotlist,
                int   *rotlist_atomId,
                int   *rotlist_rotatableBondId
               )
// The function generates the rotation list which will be stored in the constant memory field rotlist_const by
// prepare_const_fields_for_gpu(). The structure of this array is described at that function.
{
  int atom_id, rotb_id, parallel_rot_id, rotlist_id;
  int  *number_of_req_rotations;
  int  *atom_id_of_numrots;
  bool *atom_wasnt_rotated_yet;
  int max_num_rotations;
  if(myligand->num_of_atoms < NUM_OF_THREADS_PER_BLOCK)
  {
    number_of_req_rotations = (int *)  calloc(NUM_OF_THREADS_PER_BLOCK, sizeof(int));
    atom_id_of_numrots      = (int *)  calloc(NUM_OF_THREADS_PER_BLOCK, sizeof(int));
    atom_wasnt_rotated_yet  = (bool *) calloc(NUM_OF_THREADS_PER_BLOCK, sizeof(bool));
    max_num_rotations = NUM_OF_THREADS_PER_BLOCK * (myligand->num_of_rotbonds + 1);
  }
  else
  {
    number_of_req_rotations = (int *)  calloc(myligand->num_of_atoms, sizeof(int));
    atom_id_of_numrots      = (int *)  calloc(myligand->num_of_atoms, sizeof(int));
    atom_wasnt_rotated_yet  = (bool *) calloc(myligand->num_of_atoms, sizeof(bool));
    max_num_rotations = myligand->num_of_atoms * (myligand->num_of_rotbonds + 1);
  }
  if(number_of_req_rotations == NULL)
  {
    printf("ERROR: Cannot allocate memory for number_of_req_rotations of gen_rotlist() on main memory.\n");
    exit(1);
  }
  if(atom_id_of_numrots == NULL)
  {
    printf("ERROR: Cannot allocate memory for atom_id_of_numrots of gen_rotlist() on main memory.\n");
    exit(1);
  }
  if(atom_wasnt_rotated_yet == NULL)
  {
    printf("ERROR: Cannot allocate memory for atom_wasnt_rotated_yet of gen_rotlist() on main memory.\n");
    exit(1);
  }
  
  int new_rotlist_element;
  bool rotbond_found;
  int rotbond_candidate;
  int remaining_rots_around_rotbonds;
  
  myligand->num_of_rotcyc = 0;
  myligand->num_of_rotations_required = 0;


  for (atom_id=0; atom_id<NUM_OF_THREADS_PER_BLOCK; atom_id++)
    number_of_req_rotations[atom_id] = 0;

  for (atom_id=0; atom_id<myligand->num_of_atoms; atom_id++)
  {
    atom_id_of_numrots[atom_id] = atom_id;
    atom_wasnt_rotated_yet[atom_id] = true;

    number_of_req_rotations[atom_id] = 1;

    for (rotb_id=0; rotb_id<myligand->num_of_rotbonds; rotb_id++)
    {
      if (myligand->atom_rotbonds[atom_id][rotb_id] != 0)
        (number_of_req_rotations[atom_id])++;
    }

    myligand->num_of_rotations_required += number_of_req_rotations[atom_id];
  }

  rotlist_id=0;
  make_reqrot_ordering(number_of_req_rotations, atom_id_of_numrots, myligand->num_of_atoms);
  while (number_of_req_rotations[0] != 0) // if the atom with the most remaining rotations has to be rotated 0 times, done
  {
    if (rotlist_id == max_num_rotations)
    {
      printf("sans calcenergy gen_rotlist 2-3-3 WhAt?!!!\n");fflush(stdout); // sans
      return 1;
    }

    // putting the NUM_OF_THREADS_PER_BLOCK pieces of most important rotations to the list
    for (parallel_rot_id=0; parallel_rot_id<NUM_OF_THREADS_PER_BLOCK; parallel_rot_id++)
    {
      if (number_of_req_rotations[parallel_rot_id] == 0) // if the atom has not to be rotated anymore, dummy rotation
      {
        new_rotlist_element = RLIST_DUMMY_MASK;
      }
      else
      {
        atom_id = atom_id_of_numrots[parallel_rot_id];
        rotlist_atomId[rotlist_id] = atom_id;
        new_rotlist_element = ((int) atom_id) & RLIST_ATOMID_MASK;

        if (number_of_req_rotations[parallel_rot_id] == 1)
        {
          new_rotlist_element |= RLIST_GENROT_MASK;
        }
        else
        {
          rotbond_found = false;
          rotbond_candidate = myligand->num_of_rotbonds - 1;
          remaining_rots_around_rotbonds = number_of_req_rotations[parallel_rot_id] - 1; // -1 because of genrot
          while (!rotbond_found)
          {
            if (myligand->atom_rotbonds[atom_id][rotbond_candidate] != 0) // if the atom has to be rotated around current candidate
            {
              if (remaining_rots_around_rotbonds == 1) // if current value of remaining rots is 1, the proper rotbond is found
                rotbond_found = true;
              else
                remaining_rots_around_rotbonds--; // if not, decresing remaining rots (that is, skipping rotations which have to be performed later
            }

            if (!rotbond_found)
              rotbond_candidate--;

            if (rotbond_candidate < 0)
            {
              printf("sans calcenergy gen_rotlist 2-3-3 OMG!!!\n");fflush(stdout); // sans
              return 1;
            }
          }

          rotlist_rotatableBondId[rotlist_id] = rotbond_candidate;
          new_rotlist_element |= (((int) rotbond_candidate) << RLIST_RBONDID_SHIFT) & RLIST_RBONDID_MASK;
        }

        if (atom_wasnt_rotated_yet[atom_id])
        {
          new_rotlist_element |= RLIST_FIRSTROT_MASK;
        }

        // put atom_id's next rotation to rotlist
        atom_wasnt_rotated_yet[atom_id] = false;
        (number_of_req_rotations[parallel_rot_id])--;
      }


      rotlist[rotlist_id] = new_rotlist_element;
      rotlist_id++;
    }

    make_reqrot_ordering(number_of_req_rotations, atom_id_of_numrots, myligand->num_of_atoms);
    (myligand->num_of_rotcyc)++;
  }

  free(number_of_req_rotations);
  free(atom_id_of_numrots);
  free(atom_wasnt_rotated_yet);
  return 0;
}
