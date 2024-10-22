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


#include <stdio.h>
#include <errno.h>
#include "processresult.h"


void write_basic_info(
                            FILE*       fp,
                            Liganddata* ligand_ref,
                      const Dockpars*   mypars,
                      const Gridinfo*   mygrid,
                      const int*        argc,
                            char**      argv
                     )
// The function writes basic information (such as docking parameters) to the file whose file pointer is the first parameter of the function.
{
  int i;

  fprintf(fp, "***********************************\n");
  fprintf(fp, "**    OMEGADOCKER REPORT FILE    **\n");
  fprintf(fp, "***********************************\n\n\n");

  // Writing out docking parameters

  fprintf(fp, "         DOCKING PARAMETERS        \n");
  fprintf(fp, "===================================\n\n");

  fprintf(fp, "Ligand file:                               %s\n", mypars->ligandfile);
  bool flexres = false;
  if (mypars->flexresfile!=NULL){
      if ( strlen(mypars->flexresfile)>0 ) {
        fprintf(fp, "Flexible residue file:                     %s\n", mypars->flexresfile);
        flexres = true;
      }
  }
  fprintf(fp, "Grid fld file:                             %s\n", mypars->fldfile);

  fprintf(fp, "Number of energy evaluations:              %ld\n", mypars->num_of_energy_evals);
  fprintf(fp, "Number of generations:                     %ld\n", mypars->num_of_generations);
  fprintf(fp, "Size of population:                        %ld\n", mypars->pop_size);
  fprintf(fp, "Rate of crossover:                         %lf%%\n", (double) mypars->crossover_rate);
  fprintf(fp, "Tournament selection probability limit:    %lf%%\n", (double) mypars->tournament_rate);
  fprintf(fp, "Rate of mutation:                          %lf%%\n", (double) mypars->mutation_rate);
  fprintf(fp, "Maximal allowed delta movement:            +/- %lfA\n", (double) mypars->abs_max_dmov*mygrid->spacing);
  fprintf(fp, "Maximal allowed delta angle:               +/- %lf\n\n", (double) mypars->abs_max_dang);

  fprintf(fp, "Rate of local search:                      %lf%%\n", mypars->lsearch_rate);

  fprintf(fp, "Maximal number of local search iterations: %ld\n", mypars->max_num_of_iters);
  fprintf(fp, "Rho lower bound:                           %lf\n", (double) mypars->rho_lower_bound);
  fprintf(fp, "Spread of local search delta movement:     %lfA\n", (double) mypars->base_dmov_mul_sqrt3*mygrid->spacing/sqrt(3.0));
  fprintf(fp, "Spread of local search delta angle:        %lf\n", (double) mypars->base_dang_mul_sqrt3/sqrt(3.0));
  fprintf(fp, "Limit of consecutive successes/failures:   %ld\n\n", mypars->cons_limit);

  fprintf(fp, "Unbound model:                             ");
  if (mypars->unbound_model == 0)
    fprintf(fp, "BOUND\n");
  else
    if (mypars->unbound_model == 1)
      fprintf(fp, "EXTENDED\n");
    else
      fprintf(fp, "COMPACT\n");

  fprintf(fp, "Number of pdb files to be generated:       %d\n", mypars->gen_pdbs);

  fprintf(fp, "Initial population:                        ");
  if (!mypars->load_xml)
    fprintf(fp, "GENERATE\n");
  else
    fprintf(fp, "LOAD FROM FILE (%s)\n",mypars->load_xml);

#ifndef TOOLMODE
  fprintf(fp, "\n\nProgram call in command line was:          ");
  for (i=0; i<*argc; i++){
    fprintf(fp, "%s ", argv [i]);
    if (argcmp("filelist", argv[i], 'B')){
      if(mypars->filelist_files>1){
        fprintf(fp, "%s ", mypars->ligandfile);
        i+=mypars->filelist_files; // skip ahead in case there are multiple entries here
      }
    }
    if (argcmp("xml2dlg", argv[i], 'X')){
      if(mypars->xml_files>1){
        fprintf(fp, "%s ", mypars->load_xml);
        i+=mypars->xml_files; // skip ahead in case there are multiple entries here
      }
    }
  }
  fprintf(fp, "\n\n");
#endif
  fprintf(fp, "\n");

  // Writing out receptor parameters

  fprintf(fp, "        RECEPTOR PARAMETERS        \n");
  fprintf(fp, "===================================\n\n");

  fprintf(fp, "Receptor name:                             %s\n", mygrid->receptor_name.c_str());
  fprintf(fp, "Number of grid points (x, y, z):           %d, %d, %d\n", mygrid->size_xyz [0], mygrid->size_xyz [1], mygrid->size_xyz [2]);
  fprintf(fp, "Grid size (x, y, z):                       %lf, %lf, %lfA\n", mygrid->size_xyz_angstr [0], mygrid->size_xyz_angstr [1], mygrid->size_xyz_angstr [2]);
  fprintf(fp, "Grid spacing:                              %lfA\n", mygrid->spacing);
  fprintf(fp, "\n\n");

  // Writing out ligand parameters
  if(flexres)
    fprintf(fp, "     LIGAND+FLEXRES PARAMETERS     \n");
  else
    fprintf(fp, "         LIGAND PARAMETERS         \n");
  fprintf(fp, "===================================\n\n");

  fprintf(fp, "Ligand name:                               ");
  int len = strlen(mypars->ligandfile) - 6;
  for(i=0; i<len; i++) fputc(mypars->ligandfile[i], fp);
  fputc('\n', fp);
  if(flexres){
    fprintf(fp, "Flexres name:                              ");
    int len = strlen(mypars->flexresfile) - 6;
    for(i=0; i<len; i++) fputc(mypars->flexresfile[i], fp);
    fputc('\n', fp);
    fprintf(fp, "Number of ligand atoms:                    %d\n", ligand_ref->true_ligand_atoms);
    fprintf(fp, "Number of flexres atoms:                   %d\n", ligand_ref->num_of_atoms-ligand_ref->true_ligand_atoms);
    fprintf(fp, "Number of ligand rotatable bonds:          %d\n", ligand_ref->true_ligand_rotbonds);
    fprintf(fp, "Number of flexres rotatable bonds:         %d\n", ligand_ref->num_of_rotbonds-ligand_ref->true_ligand_rotbonds);
  }
  fprintf(fp, "Number of atoms:                           %d\n", ligand_ref->num_of_atoms);
  fprintf(fp, "Number of rotatable bonds:                 %d\n", ligand_ref->num_of_rotbonds);
  fprintf(fp, "Number of atom types:                      %d\n", ligand_ref->num_of_atypes);

  fprintf(fp, "Number of intraE contributors:             %d\n", ligand_ref->num_of_intraE_contributors);
  fprintf(fp, "Number of required rotations:              %d\n", ligand_ref->num_of_rotations_required);
  fprintf(fp, "Number of rotation cycles:                 %d\n", ligand_ref->num_of_rotcyc);

  fprintf(fp, "\n\n");
}

void write_basic_info_dlg(
                                FILE*       fp,
                                Liganddata* ligand_ref,
                          const Dockpars*   mypars,
                          const Gridinfo*   mygrid,
                          const int*        argc,
                                char**      argv
                         )
// The function writes basic information (such as docking parameters) to the file whose file pointer is the first parameter of the function.
{
  int i;

  if(mypars->xml2dlg && mypars->dlg2stdout) fprintf(fp, "\nXML2DLG: %s\n", mypars->load_xml);
  fprintf(fp, "OmegaDocker version: %s\n\n", VERSION);

  fprintf(fp, "**********************************************************\n");
  fprintf(fp, "**    OMEGADOCKER  AUTODOCKTOOLS-COMPATIBLE DLG FILE    **\n");
  fprintf(fp, "**********************************************************\n\n\n");

  // Writing out docking parameters

  fprintf(fp, "    DOCKING PARAMETERS\n");
  fprintf(fp, "    ________________________\n\n\n");

  fprintf(fp, "Ligand file:                               %s\n", mypars->ligandfile);
  bool flexres = false;
  if (mypars->flexresfile!=NULL){
      if ( strlen(mypars->flexresfile)>0 ) {
        fprintf(fp, "Flexible residue file:                     %s\n", mypars->flexresfile);
        flexres = true;
      }
  }
  fprintf(fp, "Grid fld file:                             %s\n\n", mypars->fldfile);

  fprintf(fp, "Random seed:                               %u", mypars->seed[0]);
  if(mypars->seed[1]>0) fprintf(fp,", %u",mypars->seed[1]);
  if(mypars->seed[2]>0) fprintf(fp,", %u",mypars->seed[2]);
  fprintf(fp, "\n");
  fprintf(fp, "Number of runs:                            %lu\n", mypars->num_of_runs);

  if(!mypars->xml2dlg){
    fprintf(fp, "Number of energy evaluations:              %ld\n", mypars->num_of_energy_evals);
    fprintf(fp, "Number of generations:                     %ld\n", mypars->num_of_generations);
    fprintf(fp, "Size of population:                        %ld\n", mypars->pop_size);
    fprintf(fp, "Rate of crossover:                         %lf%%\n", (double) mypars->crossover_rate);
    fprintf(fp, "Tournament selection probability limit:    %lf%%\n", (double) mypars->tournament_rate);
    fprintf(fp, "Rate of mutation:                          %lf%%\n", (double) mypars->mutation_rate);
    fprintf(fp, "Maximal allowed delta movement:            +/- %lfA\n", (double) mypars->abs_max_dmov*mygrid->spacing);
    fprintf(fp, "Maximal allowed delta angle:               +/- %lf\n\n", (double) mypars->abs_max_dang);

    fprintf(fp, "Rate of local search:                      %lf%%\n", mypars->lsearch_rate);

    fprintf(fp, "Maximal number of local search iterations: %ld\n", mypars->max_num_of_iters);
    fprintf(fp, "Rho lower bound:                           %lf\n", (double) mypars->rho_lower_bound);
    fprintf(fp, "Spread of local search delta movement:     %lfA\n", (double) mypars->base_dmov_mul_sqrt3*mygrid->spacing/sqrt(3.0));
    fprintf(fp, "Spread of local search delta angle:        %lf\n", (double) mypars->base_dang_mul_sqrt3/sqrt(3.0));
    fprintf(fp, "Limit of consecutive successes/failures:   %ld\n\n", mypars->cons_limit);
  }

    fprintf(fp, "Handle symmetry during clustering:         ");
  if (mypars->handle_symmetry)
    fprintf(fp, "YES\n");
  else
    fprintf(fp, "NO\n");

  fprintf(fp, "RMSD tolerance:                            %lfA\n\n", mypars->rmsd_tolerance);

#ifndef TOOLMODE
  fprintf(fp, "Program call in command line was:          ");
  for (i=0; i<*argc; i++){
    fprintf(fp, "%s ", argv [i]);
    if (argcmp("filelist", argv[i], 'B')){
      if(mypars->filelist_files>1){
        fprintf(fp, "%s ", mypars->ligandfile);
        i+=mypars->filelist_files; // skip ahead in case there are multiple entries here
      }
    }
    if (argcmp("xml2dlg", argv[i], 'X')){
      if(mypars->xml_files>1){
        fprintf(fp, "%s ", mypars->load_xml);
        i+=mypars->xml_files; // skip ahead in case there are multiple entries here
      }
    }
  }
  fprintf(fp, "\n\n");
#endif
  fprintf(fp, "\n");

  // Writing out receptor parameters

  fprintf(fp, "    GRID PARAMETERS\n");
  fprintf(fp, "    ________________________\n\n\n");

  fprintf(fp, "Receptor name:                             %s\n", mygrid->receptor_name.c_str());
  fprintf(fp, "Number of grid points (x, y, z):           %d, %d, %d\n", mygrid->size_xyz [0],
      mygrid->size_xyz [1], mygrid->size_xyz [2]);
  fprintf(fp, "Grid size (x, y, z):                       %lf, %lf, %lfA\n", mygrid->size_xyz_angstr [0],
      mygrid->size_xyz_angstr [1], mygrid->size_xyz_angstr [2]);
  fprintf(fp, "Grid spacing:                              %lfA\n", mygrid->spacing);
  fprintf(fp, "\n\n");

  // Writing out ligand parameters
  if(flexres)
  fprintf(fp, "    LIGAND+FLEXRES PARAMETERS\n");
  else
    fprintf(fp, "    LIGAND PARAMETERS\n");
  fprintf(fp, "    ________________________\n\n\n");

  fprintf(fp, "Ligand name:                               ");
  int len = strlen(mypars->ligandfile) - 6;
  for(i=0; i<len; i++) fputc(mypars->ligandfile[i], fp);
  fputc('\n', fp);
  if(flexres){
    fprintf(fp, "Flexres name:                              ");
    int len = strlen(mypars->flexresfile) - 6;
    for(i=0; i<len; i++) fputc(mypars->flexresfile[i], fp);
    fputc('\n', fp);
    fprintf(fp, "Number of ligand atoms:                    %d\n", ligand_ref->true_ligand_atoms);
    fprintf(fp, "Number of flexres atoms:                   %d\n", ligand_ref->num_of_atoms-ligand_ref->true_ligand_atoms);
    fprintf(fp, "Number of ligand rotatable bonds:          %d\n", ligand_ref->true_ligand_rotbonds);
    fprintf(fp, "Number of flexres rotatable bonds:         %d\n", ligand_ref->num_of_rotbonds-ligand_ref->true_ligand_rotbonds);
  }
  fprintf(fp, "Number of atoms:                           %d\n", ligand_ref->num_of_atoms);
  fprintf(fp, "Number of rotatable bonds:                 %d\n", ligand_ref->num_of_rotbonds);
  fprintf(fp, "Number of atom types:                      %d\n", ligand_ref->num_of_atypes);
  fprintf(fp, "\n\n");

  if(!mypars->xml2dlg){
    fprintf(fp, "    DUMMY DATA (only for ADT-compatibility)\n");
    fprintf(fp, "    ________________________\n\n\n");
    fprintf(fp, "DPF> outlev 1\n");
    fprintf(fp, "DPF> ga_run %lu\n", mypars->num_of_runs);
    fprintf(fp, "DPF> fld %s\n", mygrid->fld_name.c_str());
    fprintf(fp, "DPF> move %s\n", mypars->ligandfile);
    if(flexres) fprintf(fp, "DPF> flexres %s\n", mypars->flexresfile);
    fprintf(fp, "\n\n");
  }
}

void make_resfiles(
                         float*        final_population,
                         float*        energies,
                         IntraTables*  tables,
                         Liganddata*   ligand_ref,
                         Liganddata*   ligand_from_pdb,
                   const Liganddata*   ligand_xray,
                   const Dockpars*     mypars,
                         int           evals_performed,
                         mpz_t         generations_used,
                   const Gridinfo*     mygrid,
                   const int*          argc,
                         char**        argv,
                         int           debug,
                         int           run_cnt,
                         float&        best_energy_of_all,
                         Ligandresult* best_result
                  )
// The function writes out final_population generated by get_result
// as well as different parameters about the docking, the receptor and the ligand to a file called fdock_report.txt in a
// readable and understandable format. The ligand_from_pdb parametere must be the Liganddata which includes the original
// ligand conformation as the result conformations will be compared to this one. The structs containing the grid informations
// and docking parameters are required as well as the number and values of command line arguments. The ligand_ref parameter
// describes the ligand with the reference orientation (gene values of final_population refer to this one, that is, this can
// be moved and rotated according to the genotype values). The function returns some information about the best result wich
// was found with the best_result parameter.
{
  FILE* fp = stdout; // takes care of compile warning down below (and serves as a visual bug tracker in case fp is written to accidentally)
  int i,j;
  double entity_rmsds;
  char    *mpz_str;
  double **init_atom_idxyzq;
  init_atom_idxyzq = (double **) malloc(ligand_ref->num_of_atoms * sizeof(double *));
  if(init_atom_idxyzq == NULL)
  {
    printf("ERROR: Cannot allocate memory for init_atom_idxyzq of make_resfiles() on main memory.\n");
    exit(1);
  }
  size_t idxyzq_size = 5 * sizeof(double);
  for(i = 0; i < ligand_ref->num_of_atoms; i++)
  {
    init_atom_idxyzq[i] = (double *) malloc(idxyzq_size);
    if(init_atom_idxyzq[i] == NULL)
    {
      printf("ERROR: Cannot allocate memory for init_atom_idxyzq[%d] of make_resfiles() on main memory.\n", i);
      exit(1);
    }
    memcpy(init_atom_idxyzq[i], ligand_ref->atom_idxyzq[i], idxyzq_size);
  }
  int len = strlen(mypars->ligandfile) - 6 + 24 + 10 + 10; // length with added bits for things below (numbers below 11 digits should be a safe enough threshold)
  char* temp_filename = (char*)malloc((len+1)*sizeof(char)); // +\0 at the end
  char* name_ext_start;
  float accurate_interE;
  float accurate_intraflexE;
  float accurate_intraE;
  float accurate_interflexE;

  int rescorePopSize;
  if(mypars->num_of_runs > 1)
  {
    rescorePopSize = mypars->pop_size;
  }
  else if (mypars->num_of_runs == 1)
  {
    rescorePopSize = mypars->topN4rescore;
  }
  else
  {
    printf("ERROR: Bug? mypars->num_of_runs (%d) is not larger than zero!?\n", mypars->num_of_runs);
    exit(1);
  }
  printf("Recalculate energies by CPU for %d poses and rerank.\n", rescorePopSize);fflush(stdout);

  sprintf(temp_filename, "final_population_run%d.txt", run_cnt+1);
  
  if (mypars->gen_finalpop) // if final population files are not required, no file will be opened.
  {
    fp = fopen(temp_filename, "w");
    if(fp==NULL){
      printf("Error: Cannot create file %s for output of final population: %s\n",temp_filename,strerror(errno));
      exit(5);
    }

    write_basic_info(fp, ligand_ref, mypars, mygrid, argc, argv); // Write basic information about docking and molecule parameters to file

    fprintf(fp, "           COUNTER STATES           \n");
    fprintf(fp, "===================================\n\n");
    fprintf(fp, "Number of energy evaluations performed:    %d\n", evals_performed);
    mpz_str = mpz_get_str(NULL, 10, generations_used);
    fprintf(fp, "Number of generations used:                %s\n", mpz_str);
    free(mpz_str);
    fprintf(fp, "\n\n");
    fprintf(fp, "     STATE OF FINAL POPULATION     \n");
    fprintf(fp, "===================================\n\n");

    fprintf(fp, " Entity |      dx [A]      |      dy [A]      |      dz [A]      |      phi []      |     theta []     |  alpha_genrot [] |");
    for (i=0; i<ligand_from_pdb->num_of_rotbonds; i++)
      fprintf(fp, "  alpha_rotb%2d [] |", i);
    fprintf(fp, " intramolecular energy | intermolecular energy |     total energy calculated by CPU / calculated by GPU / difference    | RMSD [A] | \n");

    fprintf(fp, "--------+------------------+------------------+------------------+------------------+------------------+------------------+");
    for (i=0; i<ligand_from_pdb->num_of_rotbonds; i++)
      fprintf(fp, "------------------+");
    fprintf(fp, "-----------------------+-----------------------+------------------------------------------------------------------------+----------+ \n");
  }

  // Writing out state of final population

  strcpy(temp_filename, mypars->ligandfile);
  name_ext_start = temp_filename + strlen(mypars->ligandfile) - 6; // without .pdbqt
  
  bool rmsd_valid = true;
  if (mypars->given_xrayligandfile == true) {
    if(!((ligand_xray->num_of_atoms == ligand_ref->num_of_atoms) || (ligand_xray->num_of_atoms == ligand_ref->true_ligand_atoms))){
      printf("Warning: RMSD can't be calculated, atom number mismatch %d (ref) vs. %d!\n",ligand_xray->true_ligand_atoms,ligand_ref->true_ligand_atoms);
      rmsd_valid = false;
    }
  }
  else {
    if(ligand_from_pdb->true_ligand_atoms != ligand_ref->true_ligand_atoms){
      printf("Warning: RMSD can't be calculated, atom number mismatch %d (ref) vs. %d!\n",ligand_xray->true_ligand_atoms,ligand_ref->true_ligand_atoms);
      rmsd_valid = false;
    }
  }

  int actual_genotype_length = ligand_ref->num_of_rotbonds + 6;
  int genotype_length_in_globmem = actual_genotype_length + 1;
  printf("  Pose serial:");
  for(i = 0; i < rescorePopSize; i++)
  {
    printf(" [%d]", i+1);fflush(stdout);
    // start from original coordinates
    
    for(int atomI = 0; atomI < ligand_ref->num_of_atoms; atomI++)
    {
      memcpy(ligand_ref->atom_idxyzq[atomI], init_atom_idxyzq[atomI], idxyzq_size);
    }

    if(mypars->xml2dlg)
    {
      double axisangle[4];
      double *genotype;
      genotype = (double *) malloc(genotype_length_in_globmem * sizeof(double));
      if(genotype == NULL)
      {
        printf("ERROR: Cannot allocate memory for genotype of make_resfiles() on main memory.\n");
        exit(1);
      }
      
      for (int j=0; j<actual_genotype_length; j++)
      {
        genotype [j] = (final_population+i*genotype_length_in_globmem)[j];
      }
      genotype[genotype_length_in_globmem-1] = (final_population+i*genotype_length_in_globmem)[genotype_length_in_globmem-1];
      axisangle[0] = genotype[3];
      axisangle[1] = genotype[4];
      axisangle[2] = genotype[5];
      axisangle[3] = genotype[genotype_length_in_globmem-1];
      change_conform(ligand_ref, mygrid, genotype, axisangle, debug);
      free(genotype);
    }
    else
    {
      change_conform_f(ligand_ref, mygrid, final_population+i*genotype_length_in_globmem, debug);
    }
    
    if(mypars->num_of_runs > 1)
    {
      // the map interaction of flex res atoms is stored in accurate_intraflexE
      if (i == 0)
      {
        accurate_interE = calc_interE_f(mygrid, ligand_ref, mypars->ligandOnly, mypars->cautionOn, mypars->shakingRadius, best_result->shakingGridE, 0.0005, debug, accurate_intraflexE, &(best_result->interE_elec), best_result->peratom_vdw, best_result->peratom_elec); // calculate intermolecular and per atom energies
      }
      else
      {
        accurate_interE = calc_interE_f(mygrid, ligand_ref, mypars->ligandOnly, mypars->cautionOn, mypars->shakingRadius, best_result->shakingGridE, 0.0005, debug, accurate_intraflexE); // calculating the intermolecular energy
      }
      
      if (mypars->contact_analysis && (i==0)){
        best_result->analysis = analyze_ligand_receptor(mygrid, ligand_ref, mypars->receptor_atoms.data(), mypars->receptor_map, mypars->receptor_map_list, 0.0005, debug, mypars->H_cutoff, mypars->V_cutoff);
      }
    }
    else if(mypars->num_of_runs == 1)
    {
      accurate_interE = calc_interE_f(mygrid, ligand_ref, mypars->ligandOnly, mypars->cautionOn, mypars->shakingRadius, best_result[i].shakingGridE, 0.0005, debug, accurate_intraflexE, &(best_result[i].interE_elec), best_result[i].peratom_vdw, best_result[i].peratom_elec); // calculate intermolecular and per atom energies
      
      if (mypars->contact_analysis){
        best_result[i].analysis = analyze_ligand_receptor(mygrid, ligand_ref, mypars->receptor_atoms.data(), mypars->receptor_map, mypars->receptor_map_list, 0.0005, debug, mypars->H_cutoff, mypars->V_cutoff);
      }
    }
    else
    {
      printf("ERROR: Bug? mypars->num_of_runs (%d) is not larger than zero!?\n", mypars->num_of_runs);
      exit(1);
    }

    scale_ligand(ligand_ref, mygrid->spacing);

    if(mypars->num_of_runs > 1)
    {
      // the interaction between flex res and ligand is stored in accurate_interflexE
      if(mypars->contact_analysis && (i==0))
      {
        accurate_intraE = calc_intraE_f(ligand_ref, mygrid, mypars->shakingRadius, best_result->shakingGridE, best_result->shakingGridEL, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, debug, accurate_interflexE, &(best_result->analysis), mypars->receptor_atoms.data() + mypars->nr_receptor_atoms, mypars->R_cutoff, mypars->H_cutoff, mypars->V_cutoff);
      }
      else
      {
        accurate_intraE = calc_intraE_f(ligand_ref, mygrid, mypars->shakingRadius, best_result->shakingGridE, best_result->shakingGridEL, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, debug, accurate_interflexE);
      }
      
      move_ligand(ligand_ref, mygrid->origo_real_xyz, mygrid->origo_real_xyz); //moving it according to grid location
      
      if ((mypars->gen_finalpop) || (i==0)){ // rmsd value is only needed in either of those cases
        if (rmsd_valid){
          if (mypars->given_xrayligandfile)
            entity_rmsds = calc_rmsd(ligand_xray->atom_idxyzq, ligand_ref->atom_idxyzq, ligand_xray->num_of_atoms, mypars->handle_symmetry); //calculating rmds compared to original xray file
          else
            entity_rmsds = calc_rmsd(ligand_from_pdb->atom_idxyzq, ligand_ref->atom_idxyzq, ligand_from_pdb->true_ligand_atoms, mypars->handle_symmetry); //calculating rmds compared to original pdb file
        } else entity_rmsds = 100000;
      }
      
      // copying best result to output parameter
      if (i == 0) // assuming this is the best one (final_population is arranged), however, the
      {           // arrangement was made according to the inaccurate values calculated by FPGA
        best_result->genotype = final_population+i*genotype_length_in_globmem;
        best_result->interE = accurate_interE;
        best_result->interflexE = accurate_interflexE;
        best_result->intraE = accurate_intraE;
        best_result->intraflexE = accurate_intraflexE;
        calcEntropy(best_result, mygrid, mypars, ligand_ref);
        for(int i = 0; i < ligand_ref->num_of_atoms; i++)
        {
          memcpy(best_result->atom_idxyzq[i], ligand_ref->atom_idxyzq[i], 5 * sizeof(double));
        }
        best_result->rmsd_from_ref = entity_rmsds;
        best_result->run_number = run_cnt+1;
        if(mypars->contact_analysis){
          // sort by analysis type
          for(unsigned int j=0; j<best_result->analysis.size(); j++)
            for(unsigned int k=0; k<best_result->analysis.size()-j-1; k++)
              if(best_result->analysis[k].type>best_result->analysis[k+1].type) // percolate larger types numbers up
                std::swap(best_result->analysis[k], best_result->analysis[k+1]);
        }
      }
    }
    else if(mypars->num_of_runs == 1)
    {
      if(mypars->contact_analysis)
      {
        accurate_intraE = calc_intraE_f(ligand_ref, mygrid, mypars->shakingRadius, best_result[i].shakingGridE, best_result[i].shakingGridEL, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, debug, accurate_interflexE, &(best_result[i].analysis), mypars->receptor_atoms.data() + mypars->nr_receptor_atoms, mypars->R_cutoff, mypars->H_cutoff, mypars->V_cutoff);
      }
      else
      {
        accurate_intraE = calc_intraE_f(ligand_ref, mygrid, mypars->shakingRadius, best_result[i].shakingGridE, best_result[i].shakingGridEL, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, debug, accurate_interflexE);
      }
      
      move_ligand(ligand_ref, mygrid->origo_real_xyz, mygrid->origo_real_xyz); //moving it according to grid location
      
      if (rmsd_valid){
        if (mypars->given_xrayligandfile)
          entity_rmsds = calc_rmsd(ligand_xray->atom_idxyzq, ligand_ref->atom_idxyzq, ligand_xray->num_of_atoms, mypars->handle_symmetry); //calculating rmds compared to original xray file
        else
          entity_rmsds = calc_rmsd(ligand_from_pdb->atom_idxyzq, ligand_ref->atom_idxyzq, ligand_from_pdb->true_ligand_atoms, mypars->handle_symmetry); //calculating rmds compared to original pdb file
      } else entity_rmsds = 100000;
      
      best_result[i].genotype = final_population+i*genotype_length_in_globmem;
      best_result[i].interE = accurate_interE;
      best_result[i].interflexE = accurate_interflexE;
      best_result[i].intraE = accurate_intraE;
      best_result[i].intraflexE = accurate_intraflexE;
      calcEntropy(best_result + i, mygrid, mypars, ligand_ref);
      for(int j = 0; j < ligand_ref->num_of_atoms; j++)
      {
        memcpy(best_result[i].atom_idxyzq[j], ligand_ref->atom_idxyzq[j], 5 * sizeof(double));
      }
      best_result[i].rmsd_from_ref = entity_rmsds;
      best_result[i].run_number = run_cnt+1;
      if(mypars->contact_analysis){
        // sort by analysis type
        for(unsigned int j=0; j<best_result[i].analysis.size(); j++)
          for(unsigned int k=0; k<best_result[i].analysis.size()-j-1; k++)
            if(best_result[i].analysis[k].type>best_result[i].analysis[k+1].type) // percolate larger types numbers up
              std::swap(best_result[i].analysis[k], best_result[i].analysis[k+1]);
      }
    }
    else
    {
      printf("ERROR: Bug? mypars->num_of_runs (%d) is not larger than zero!?\n", mypars->num_of_runs);
      exit(1);
    }
    
    if (i < mypars->gen_pdbs) //if it is necessary, making new pdbqts for best entities
    {
      sprintf(name_ext_start, "_docked_run%d_entity%d.pdbqt", run_cnt+1, i+1); //name will be <original pdb filename>_docked_<number starting from 1>.pdb
      gen_new_pdbfile(mypars->ligandfile, temp_filename, ligand_ref);
    }
    if (mypars->gen_finalpop)
    {
      fprintf(fp, "  %3d   |", i+1);

      for (j=0; j<3; j++)
      {
        fprintf(fp, "    %10.3f    |", final_population [i*genotype_length_in_globmem+j]*(mygrid->spacing));
      }
      for (j=3; j<6+ligand_from_pdb->num_of_rotbonds; j++)
      {
        fprintf(fp, "    %10.3f    |", final_population [i*genotype_length_in_globmem+j]);
      }

      fprintf(fp, " %21.3f |", accurate_intraE);
      fprintf(fp, " %21.3f |", accurate_interE);
      fprintf(fp, "  %21.3f / %21.3f / %21.3f |", accurate_intraE + accurate_interE, energies[i], energies[i] - (accurate_intraE + accurate_interE));

      fprintf(fp, " %8.3lf | \n", entity_rmsds);
    }
    if(mypars->num_of_runs > 1 && mypars->savetime) break;
  }
  printf("\n\n");fflush(stdout);
  
  for(i = 0; i < ligand_ref->num_of_atoms; i++)
  {
    memcpy(ligand_ref->atom_idxyzq[i], init_atom_idxyzq[i], idxyzq_size);
  }
  if (mypars->gen_finalpop) fclose(fp);
  free(temp_filename);
  for(i = 0; i < ligand_ref->num_of_atoms; i++)
  {
    free(init_atom_idxyzq[i]);
  }
  free(init_atom_idxyzq);
}

void ligand_calc_output(
                              FILE*         fp,
                        const char*         prefix,
                              IntraTables*  tables,
                        const Liganddata*   ligand,
                        const Dockpars*     mypars,
                        const Gridinfo*     mygrid,
                              bool          output_analysis,
                              bool          output_energy
                       )
{
  if(output_analysis || output_energy)
  {
    Liganddata *calc_lig;
    calc_lig = (Liganddata *) malloc(sizeof(Liganddata));
    if(calc_lig == NULL)
    {
      printf("ERROR: Cannot allocate memory for calc_lig of ligand_calc_output() on main memory.\n");
      exit(1);
    }
    LiganddataCopy(calc_lig, ligand);

    Ligandresult calc;
    calc.atom_idxyzq = (double **) malloc(ligand->num_of_atoms * sizeof(double *));
    if(calc.atom_idxyzq == NULL)
    {
      printf("ERROR: Cannot allocate memory for calc.atom_idxyzq\n");
      exit(1);
    }
    calc.shakingGridE = (float **) malloc(ligand->num_of_atoms * sizeof(float *));
    if(calc.shakingGridE == NULL)
    {
      printf("ERROR: Cannot allocate memory for calc.shakingGridE\n");
      exit(1);
    }
    for(int i = 0; i < ligand->num_of_atoms; i++)
    {
      calc.atom_idxyzq[i] = (double *) malloc(5 * sizeof(double));
      if(calc.atom_idxyzq[i] == NULL)
      {
        printf("ERROR: Cannot allocate memory for calc.atom_idxyzq[i]\n");
        exit(1);
      }
      
      calc.shakingGridE[i] = NULL;
    }
    
    calc.peratom_vdw = (float *) malloc(ligand->num_of_atoms * sizeof(float));
    if(calc.peratom_vdw == NULL)
    {
      printf("ERROR: Cannot allocate memory for calc.peratom_vdw\n");
      exit(1);
    }
    
    calc.peratom_elec = (float *) malloc(ligand->num_of_atoms * sizeof(float));
    if(calc.peratom_elec == NULL)
    {
      printf("ERROR: Cannot allocate memory for calc.peratom_elec\n");
      exit(1);
    }
    double orig_vec[3];
    for (unsigned int i=0; i<3; i++)
      orig_vec [i] = -mygrid->origo_real_xyz [i];
    move_ligand(/*&*/calc_lig, orig_vec, orig_vec);
    scale_ligand(/*&*/calc_lig, 1.0/mygrid->spacing);
    calc.interE = calc_interE_f(mygrid, /*&*/calc_lig, mypars->ligandOnly, mypars->cautionOn, mypars->shakingRadius, calc.shakingGridE, 0.0005, 0, calc.intraflexE, &(calc.interE_elec), calc.peratom_vdw, calc.peratom_elec);
    if (output_analysis){
      calc.analysis = analyze_ligand_receptor(mygrid, /*&*/calc_lig, mypars->receptor_atoms.data(), mypars->receptor_map, mypars->receptor_map_list, 0.0005, 0, mypars->H_cutoff, mypars->V_cutoff);
    }
    scale_ligand(/*&*/calc_lig, mygrid->spacing);
    // the interaction between flex res and ligand is stored in accurate_interflexE
    if(output_analysis)
      calc.intraE = calc_intraE_f(/*&*/calc_lig, mygrid, mypars->shakingRadius, calc.shakingGridE, calc.shakingGridEL, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, 0, calc.interflexE, &(calc.analysis), mypars->receptor_atoms.data() + mypars->nr_receptor_atoms, mypars->R_cutoff, mypars->H_cutoff, mypars->V_cutoff);
    else
      calc.intraE = calc_intraE_f(/*&*/calc_lig, mygrid, mypars->shakingRadius, calc.shakingGridE, calc.shakingGridEL, 8, mypars->smooth, 0, mypars->elec_min_distance, tables, 0, calc.interflexE);
    move_ligand(/*&*/calc_lig, mygrid->origo_real_xyz, mygrid->origo_real_xyz); //moving it according to grid location
    calcEntropy(&calc, mygrid, mypars, calc_lig);
    
    if (output_analysis){
      // sort by analysis type
      for(unsigned int j=0; j<calc.analysis.size(); j++)
        for(unsigned int k=0; k<calc.analysis.size()-j-1; k++)
          if(calc.analysis[k].type>calc.analysis[k+1].type) // percolate larger types numbers up
            std::swap(calc.analysis[k], calc.analysis[k+1]);
      if(calc.analysis.size()>0){
        fprintf(fp, "ANALYSIS: COUNT %lu\n", calc.analysis.size());
        std::string types    = "TYPE    {";
        std::string lig_id   = "LIGID   {";
        std::string ligname  = "LIGNAME {";
        std::string rec_id   = "RECID   {";
        std::string rec_name = "RECNAME {";
        std::string residue  = "RESIDUE {";
        std::string res_id   = "RESID   {";
        std::string chain    = "CHAIN   {";
        char item[8], pad[8];
        for(unsigned int j=0; j<calc.analysis.size(); j++){
          if(j>0){
            types    += ",";
            lig_id   += ",";
            ligname  += ",";
            rec_id   += ",";
            rec_name += ",";
            residue  += ",";
            res_id   += ",";
            chain    += ",";
          }
          switch(calc.analysis[j].type){
            case 0: types += "   \"R\"";
                    break;
            case 1: types += "   \"H\"";
                    break;
            default:
            case 2: types += "   \"V\"";
                    break;
          }
          sprintf(item, "%5d ", calc.analysis[j].lig_id);   lig_id+=item;
          sprintf(item, "\"%s\"", calc.analysis[j].lig_name); sprintf(pad, "%6s", item); ligname+=pad;
          sprintf(item, "%5d ", calc.analysis[j].rec_id);   rec_id+=item;
          sprintf(item, "\"%s\"", calc.analysis[j].rec_name); sprintf(pad, "%6s", item); rec_name+=pad;
          sprintf(item, "\"%s\"", calc.analysis[j].residue); sprintf(pad, "%6s", item);  residue+=pad;
          sprintf(item, "%5d ", calc.analysis[j].res_id);   res_id+=item;
          sprintf(item, "\"%s\"", calc.analysis[j].chain); sprintf(pad, "%6s", item);    chain+=pad;
        }
        fprintf(fp, "ANALYSIS: %s}\n", types.c_str());
        fprintf(fp, "ANALYSIS: %s}\n", lig_id.c_str());
        fprintf(fp, "ANALYSIS: %s}\n", ligname.c_str());
        fprintf(fp, "ANALYSIS: %s}\n", rec_id.c_str());
        fprintf(fp, "ANALYSIS: %s}\n", rec_name.c_str());
        fprintf(fp, "ANALYSIS: %s}\n", residue.c_str());
        fprintf(fp, "ANALYSIS: %s}\n", res_id.c_str());
        fprintf(fp, "ANALYSIS: %s}\n\n", chain.c_str());
      }
    }
    if(output_energy){
      fprintf(fp, "%s    Estimated Free Energy of Binding (bad) =", prefix);
      PRINT1000(fp, ((float) (calc.interE + calc.interflexE + calc.intraE + calc.intraflexE + calc.nts)));
      fprintf(fp, " kcal/mol  [=(1)+(2)+(3)]\n");
      fprintf(fp, "%s\n", prefix);
      fprintf(fp, "%s    (1) Final Intermolecular Energy     =", prefix);
      PRINT1000(fp, ((float) (calc.interE + calc.interflexE)));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s        vdW + Hbond + desolv Energy     =", prefix);
      PRINT1000(fp, ((float) (calc.interE - calc.interE_elec)));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s        Electrostatic Energy            =", prefix);
      PRINT1000(fp, ((float) calc.interE_elec));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s        Moving Ligand-Fixed Receptor    =", prefix);
      PRINT1000(fp, ((float) calc.interE));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s        Moving Ligand-Moving Receptor   =", prefix);
      PRINT1000(fp, ((float) calc.interflexE));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s    (2) Final Total Internal Energy     =", prefix);
      PRINT1000(fp, ((float) (calc.intraE + calc.intraflexE)));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s    (1)+(2) Enthalpy of Free Energy     =", prefix);
      PRINT1000(fp, ((float) (calc.interE + calc.interflexE + calc.intraE + calc.intraflexE)));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s    (3) -TS of Free Energy              =", prefix);
      PRINT1000(fp, ((float) calc.nts));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s    (4) Unbound System's Energy         =", prefix);
      PRINT1000(fp, ((float) (calc.intraE + calc.intraflexE + calc.ntsL)));
      fprintf(fp, " kcal/mol\n");
      fprintf(fp, "%s\n", prefix);
    }
    LiganddataFree(calc_lig);
    
    for(int i = 0; i < ligand->num_of_atoms; i++)
    {
      free(calc.atom_idxyzq[i]);
      free(calc.shakingGridE[i]);
    }
    free(calc.atom_idxyzq);
    free(calc.shakingGridE);
    free(calc.peratom_vdw);
    free(calc.peratom_elec);
  }
}

void generate_output(
                           Ligandresult *myresults,
                           int           num_of_runs,
                           IntraTables*  tables,
                           Liganddata*   ligand_ref,
                     const Liganddata*   ligand_xray,
                     const Dockpars*     mypars,
                     const Gridinfo*     mygrid,
                     const int*          argc,
                           char**        argv,
                     const double        docking_avg_runtime,
                           mpz_t         generations_used,
                           mpz_t         evals_performed,
                           double        exec_time,
                           double        idle_time
                    )
// The function performs ranked cluster analysis similar to that of AutoDock and creates a file with report_file_name name, the result
// will be written to it.
{
  int i, j, atom_cnt;
  int num_of_clusters = 0;
  int current_clust_center;
  double temp_rmsd;
  int result_clustered;
  int subrank;
  FILE* fp = stdout;
  FILE* fp_xml;
  char   *mpz_str;
  int    *cluster_sizes;
  double *sum_energy;
  double *best_energy;
  int    *best_energy_runid;
  float   lowestUnboundG;
  float   lowestUnboundG_H;
  float   lowestUnboundG_ntsL;
  
  if(num_of_runs > 1)
  {
    cluster_sizes     = (int *)    malloc(num_of_runs * sizeof(int));
    sum_energy        = (double *) malloc(num_of_runs * sizeof(double));
    best_energy       = (double *) malloc(num_of_runs * sizeof(double));
    best_energy_runid = (int *)    malloc(num_of_runs * sizeof(int));
  }
  else if(num_of_runs == 1)
  {
    cluster_sizes     = (int *)    malloc(mypars->topN4rescore * sizeof(int));
    sum_energy        = (double *) malloc(mypars->topN4rescore * sizeof(double));
    best_energy       = (double *) malloc(mypars->topN4rescore * sizeof(double));
    best_energy_runid = (int *)    malloc(mypars->topN4rescore * sizeof(int));
  }
  else
  {
    printf("ERROR: Bug? num_of_runs (%d) is not larger than zero!?\n", num_of_runs);
    exit(1);
  }
  char  *tempstr = NULL;
  size_t tempstrSize;
  char fileName4bestN[512];
  char **bestPdbqt;
  char **bestPdb;
  double cluster_tolerance = mypars->rmsd_tolerance;

  // // first of all, let's calculate the constant torsional free energy term
  // double torsional_energy = mypars->coeffs.AD4_coeff_tors * ligand_ref->true_ligand_rotbonds;

  int len = strlen(mypars->resname) + 4 + 1;

  int genotype_length_in_globmem = ligand_ref->num_of_rotbonds + 7;

  // GENERATING DLG FILE
  if(mypars->output_dlg){
    if(!mypars->dlg2stdout){
      char* report_file_name = (char*)malloc(len*sizeof(char));
      strcpy(report_file_name, mypars->resname);
      strcat(report_file_name, ".dlg");
      fp = fopen(report_file_name, "w");
      if(fp==NULL){
        printf("Error: Cannot create dlg output file %s: %s\n",report_file_name,strerror(errno));
        exit(7);
      }
      free(report_file_name);
    }

    // writing basic info
    write_basic_info_dlg(fp, ligand_ref, mypars, mygrid, argc, argv);

    if(!mypars->xml2dlg){
      fprintf(fp, "    COUNTER STATES\n");
      fprintf(fp, "    ________________________\n\n\n");
      mpz_str = mpz_get_str(NULL, 10, evals_performed);
      fprintf(fp, "Number of energy evaluations performed:    %s\n", mpz_str);
      free(mpz_str);
      mpz_str = mpz_get_str(NULL, 10, generations_used);
      fprintf(fp, "Number of generations used:                %s\n", mpz_str);
      free(mpz_str);
      fprintf(fp, "\n\n");
    }
    std::string pdbqt_template;
    std::vector<unsigned int> atom_data;
    std::string pdbqt_out;
    std::string pdb_out;
    std::vector<unsigned int> pdbqt_out_data;
    std::vector<unsigned int> pdb_out_data;
    char  *lineout = NULL;
    size_t lineoutSize;
    bool output_ref_calcs = mypars->reflig_en_required;
    if(mypars->given_xrayligandfile)
    {
      // writing xray ligand pdbqt file
      fprintf(fp, "    XRAY LIGAND PDBQT FILE:\n");
      fprintf(fp, "    ________________________\n\n\n");
      ligand_calc_output(fp, "XRAY-LIGAND-PDBQT: USER", tables, ligand_xray, mypars, mygrid, mypars->contact_analysis, output_ref_calcs);
      if(output_ref_calcs) output_ref_calcs=false;
      unsigned int line_count = 0;
      while (line_count < ligand_xray->ligand_line_count)
      {
        if(tempstr != NULL)
        {
          printf("ERROR: Bug? The tempstr should point to NULL.\n");
          exit(1);
        }
        tempstrSize = strlen(ligand_xray->fileData[line_count]) + 1;
        tempstr = (char *) malloc(tempstrSize);
        strncpy(tempstr, ligand_xray->fileData[line_count], tempstrSize);
        line_count++;
        fprintf(fp, "XRAY-LIGAND-PDBQT: %s", tempstr);
        free(tempstr);
        tempstr = NULL;
      }
      fprintf(fp, "\n\n");
    }
    // writing input pdbqt file
    fprintf(fp, "    INPUT LIGAND PDBQT FILE:\n    ________________________\n\n\n");
    ligand_calc_output(fp, "INPUT-LIGAND-PDBQT: USER", tables, ligand_ref, mypars, mygrid, mypars->contact_analysis, output_ref_calcs);
    unsigned int line_count = 0;
    while (line_count < ligand_ref->ligand_line_count)
    {
      if(tempstr != NULL)
      {
        printf("ERROR: Bug? The tempstr should point to NULL.\n");
        exit(1);
      }
      tempstrSize = strlen(ligand_ref->fileData[line_count]) + 1;
      tempstr = (char *) malloc(tempstrSize);
      if(tempstr == NULL)
      {
        printf("ERROR: Cannot allocate memory for tempstr.\n");
        exit(1);
      }
      strncpy(tempstr, ligand_ref->fileData[line_count], tempstrSize);
      line_count++;
      fprintf(fp, "INPUT-LIGAND-PDBQT: %s", tempstr);
      if ((strncmp("ATOM", tempstr, 4) == 0) || (strncmp("HETATM", tempstr, 6) == 0))
      {
        if(tempstrSize > 30)
        {
          tempstr[30] = '\0';
        }

        if(lineout != NULL)
        {
          printf("ERROR: Bug? The lineout should point to NULL.\n");
          exit(1);
        }
        lineoutSize = strlen(tempstr) + 9;
        lineout = (char *) malloc(lineoutSize);
        if(lineout == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout.\n");
          exit(1);
        }
        pdbqt_out += tempstr;
        pdbqt_out_data.push_back(pdbqt_out.length());
        pdb_out += tempstr;
        pdb_out_data.push_back(pdb_out.length());
        sprintf(lineout, "DOCKED: %s", tempstr);
        pdbqt_template += lineout;
        atom_data.push_back(pdbqt_template.length());
      }
      else
      {
        if (strncmp("ROOT", tempstr, 4) == 0)
        {
          pdbqt_template += "DOCKED: USER                              x       y       z     vdW  Elec       q    Type\n";
          pdbqt_template += "DOCKED: USER                           _______ _______ _______ _____ _____    ______ ____\n";
          pdbqt_out += "USER                              x       y       z     vdW  Elec       q    Type\n";
          pdbqt_out += "USER                           _______ _______ _______ _____ _____    ______ ____\n";
          pdb_out += "USER                              x       y       z     vdW  Elec       q    Type\n";
          pdb_out += "USER                           _______ _______ _______ _____ _____    ______ ____\n";
        }
        if(lineout != NULL)
        {
          printf("ERROR: Bug? The lineout should point to NULL.\n");
          exit(1);
        }
        lineoutSize = strlen(tempstr) + 9;
        lineout = (char *) malloc(lineoutSize);
        if(lineout == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout.\n");
          exit(1);
        }
        pdbqt_out += tempstr;
        if(strncmp("ROOT",      tempstr, 4) != 0 &&
           strncmp("ENDROOT",   tempstr, 7) != 0 &&
           strncmp("BRANCH",    tempstr, 6) != 0 &&
           strncmp("ENDBRANCH", tempstr, 9) != 0
          )
        {
          if(strncmp("TORSDOF", tempstr, 7) == 0)
          {
            pdb_out += "END\n";
          }
          else
          {
            pdb_out += tempstr;
          }
        }
        sprintf(lineout, "DOCKED: %s", tempstr);
        pdbqt_template += lineout;
      }
      free(tempstr);
      tempstr = NULL;
      free(lineout);
      lineout = NULL;
    }
    fprintf(fp, "\n\n");

    // writing input flexres pdbqt file if specified
    if (mypars->flexresfile!=NULL) {
      if ( strlen(mypars->flexresfile)>0 ) {
        fprintf(fp, "    INPUT FLEXRES PDBQT FILE:\n    ________________________\n\n\n");
        while (line_count < ligand_ref->fileDataLineNum)
        {
          if(tempstr != NULL)
          {
            printf("ERROR: Bug? The tempstr should point to NULL.\n");
            exit(1);
          }
          tempstrSize = strlen(ligand_ref->fileData[line_count]) + 1;
          tempstr = (char *) malloc(tempstrSize);
          if(tempstr == NULL)
          {
            printf("ERROR: Cannot allocate memory for tempstr.\n");
            exit(1);
          }
          strncpy(tempstr, ligand_ref->fileData[line_count], tempstrSize);
          line_count++;
          fprintf(fp, "INPUT-FLEXRES-PDBQT: %s", tempstr);
          if ((strncmp("ATOM", tempstr, 4) == 0) || (strncmp("HETATM", tempstr, 6) == 0))
          {
            tempstr[30] = '\0';
            if(lineout != NULL)
            {
              printf("ERROR: Bug? The lineout should point to NULL.\n");
              exit(1);
            }
            lineoutSize = strlen(tempstr) + 9;
            lineout = (char *) malloc(lineoutSize);
            if(lineout == NULL)
            {
              printf("ERROR: Cannot allocate memory for lineout.\n");
              exit(1);
            }
            pdbqt_out += tempstr;
            pdbqt_out_data.push_back(pdbqt_out.length());
            pdb_out += tempstr;
            pdb_out_data.push_back(pdb_out.length());
            sprintf(lineout, "DOCKED: %s", tempstr);
            pdbqt_template += lineout;
            atom_data.push_back(pdbqt_template.length());
          }
          else
          {
            if (strncmp("ROOT", tempstr, 4) == 0)
            {
              pdbqt_template += "DOCKED: USER                              x       y       z     vdW  Elec       q    Type\n";
              pdbqt_template += "DOCKED: USER                           _______ _______ _______ _____ _____    ______ ____\n";
              pdbqt_out += "USER                              x       y       z     vdW  Elec       q    Type\n";
              pdbqt_out += "USER                           _______ _______ _______ _____ _____    ______ ____\n";
              pdb_out += "USER                              x       y       z     vdW  Elec       q    Type\n";
              pdb_out += "USER                           _______ _______ _______ _____ _____    ______ ____\n";
            }
            if(lineout != NULL)
            {
              printf("ERROR: Bug? The lineout should point to NULL.\n");
              exit(1);
            }
            lineoutSize = strlen(tempstr) + 9;
            lineout = (char *) malloc(lineoutSize);
            if(lineout == NULL)
            {
              printf("ERROR: Cannot allocate memory for lineout.\n");
              exit(1);
            }
            pdbqt_out += tempstr;
            if(strncmp("ROOT",      tempstr, 4) != 0 &&
               strncmp("ENDROOT",   tempstr, 7) != 0 &&
               strncmp("BRANCH",    tempstr, 6) != 0 &&
               strncmp("ENDBRANCH", tempstr, 9) != 0
              )
            {
              if(strncmp("TORSDOF", tempstr, 7) == 0)
              {
                pdb_out += "END\n";
              }
              else
              {
                pdb_out += tempstr;
              }
            }
            sprintf(lineout, "DOCKED: %s", tempstr);
            pdbqt_template += lineout;
          }
          free(tempstr);
          tempstr = NULL;
          free(lineout);
          lineout = NULL;
        }
        fprintf(fp, "\n\n");
      }
    }

    // writing docked conformations
    std::string curr_model;
    std::string curr_pdbqt;
    std::string curr_pdb;
    if(num_of_runs > 1)
    {
      bestPdbqt = (char **) malloc(num_of_runs * sizeof(char*));
      if(bestPdbqt == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for bestPdbqt.\n");
        exit(1);
      }
      bestPdb = (char **) malloc(num_of_runs * sizeof(char*));
      if(bestPdb == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for bestPdb.\n");
        exit(1);
      }
      myresults[0].interE     *= mypars->enthalpyScaling;
      myresults[0].interflexE *= mypars->enthalpyScaling;
      myresults[0].intraE     *= mypars->enthalpyScaling;
      myresults[0].intraflexE *= mypars->enthalpyScaling;
      lowestUnboundG_H    = myresults[0].intraE + myresults[0].intraflexE;
      lowestUnboundG_ntsL = myresults[0].ntsL;
      lowestUnboundG      = lowestUnboundG_H + lowestUnboundG_ntsL;
      for(i = 1; i < num_of_runs; i++)
      {
        myresults[i].interE     *= mypars->enthalpyScaling;
        myresults[i].interflexE *= mypars->enthalpyScaling;
        myresults[i].intraE     *= mypars->enthalpyScaling;
        myresults[i].intraflexE *= mypars->enthalpyScaling;
        float UE = myresults[i].intraE + myresults[i].intraflexE + myresults[i].ntsL;
        if(UE < lowestUnboundG)
        {
          lowestUnboundG_H    = myresults[i].intraE + myresults[i].intraflexE;
          lowestUnboundG_ntsL = myresults[i].ntsL;
          lowestUnboundG      = UE;
        }
      }
      for (i=0; i<num_of_runs; i++)
      {
        fprintf(fp, "    FINAL DOCKED STATE:\n    ________________________\n\n\n");
      
        fprintf(fp, "Run:   %d / %lu\n", i+1, mypars->num_of_runs);
        fprintf(fp, "Time taken for this run:   %.3lfs\n\n", docking_avg_runtime);
        if(mypars->contact_analysis){
          if(myresults[i].analysis.size()>0){
            fprintf(fp, "ANALYSIS: COUNT %lu\n", myresults[i].analysis.size());
            std::string types    = "TYPE    {";
            std::string lig_id   = "LIGID   {";
            std::string ligname  = "LIGNAME {";
            std::string rec_id   = "RECID   {";
            std::string rec_name = "RECNAME {";
            std::string residue  = "RESIDUE {";
            std::string res_id   = "RESID   {";
            std::string chain    = "CHAIN   {";
            char item[8], pad[8];
            for(unsigned int j=0; j<myresults[i].analysis.size(); j++){
              if(j>0){
                types    += ",";
                lig_id   += ",";
                ligname  += ",";
                rec_id   += ",";
                rec_name += ",";
                residue  += ",";
                res_id   += ",";
                chain    += ",";
              }
              switch(myresults[i].analysis[j].type){
                case 0: types += "   \"R\"";
                        break;
                case 1: types += "   \"H\"";
                        break;
                default:
                case 2: types += "   \"V\"";
                        break;
              }
              sprintf(item, "%5d ", myresults[i].analysis[j].lig_id);   lig_id+=item;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].lig_name); sprintf(pad, "%6s", item); ligname+=pad;
              sprintf(item, "%5d ", myresults[i].analysis[j].rec_id);   rec_id+=item;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].rec_name); sprintf(pad, "%6s", item); rec_name+=pad;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].residue); sprintf(pad, "%6s", item);  residue+=pad;
              sprintf(item, "%5d ", myresults[i].analysis[j].res_id);   res_id+=item;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].chain); sprintf(pad, "%6s", item);    chain+=pad;
            }
            fprintf(fp, "ANALYSIS: %s}\n", types.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", lig_id.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", ligname.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", rec_id.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", rec_name.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", residue.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", res_id.c_str());
            fprintf(fp, "ANALYSIS: %s}\n\n", chain.c_str());
          }
        }
        
        fprintf(fp, "DOCKED: MODEL        %d\n", i+1);
        fprintf(fp, "DOCKED: USER    Run = %d\n", i+1);
        fprintf(fp, "DOCKED: USER\n");
        
        fprintf(fp, "DOCKED: USER    Estimated Free Energy of Binding    =");

        myresults[i].complexG_H   = myresults[i].interE + myresults[i].interflexE + myresults[i].intraE + myresults[i].intraflexE;
        myresults[i].complexG_nts = myresults[i].nts;
        myresults[i].complexG     = myresults[i].complexG_H + myresults[i].complexG_nts;
        PRINT1000(fp, (myresults[i].complexG - lowestUnboundG));
        fprintf(fp, "DOCKED: USER                                        = (1)+(2)+(3)-(4)\n");
        
        fprintf(fp, "DOCKED: USER\n");
        
        fprintf(fp, "DOCKED: USER    (1) Final Intermolecular Energy     =");
        PRINT1000(fp, ((float) (myresults[i].interE + myresults[i].interflexE)));
        
        fprintf(fp, "DOCKED: USER        vdW + Hbond + desolv Energy     =");
        PRINT1000(fp, ((float) (myresults[i].interE - myresults[i].interE_elec)));
        
        fprintf(fp, "DOCKED: USER        Electrostatic Energy            =");
        PRINT1000(fp, ((float) myresults[i].interE_elec));
        
        fprintf(fp, "DOCKED: USER        Moving Ligand-Fixed Receptor    =");
        PRINT1000(fp, ((float) myresults[i].interE));
        
        fprintf(fp, "DOCKED: USER        Moving Ligand-Moving Receptor   =");
        PRINT1000(fp, ((float) myresults[i].interflexE));
        
        fprintf(fp, "DOCKED: USER    (2) Final Total Internal Energy     =");
        PRINT1000(fp, ((float) (myresults[i].intraE + myresults[i].intraflexE)));
        
        fprintf(fp, "DOCKED: USER    (1)+(2)        Enthalpy of Complex  =");
        PRINT1000(fp, myresults[i].complexG_H);
        fprintf(fp, "DOCKED: USER    (3)                 -TS of Complex  =");
        PRINT1000(fp, myresults[i].complexG_nts);
        fprintf(fp, "DOCKED: USER    (1)+(2)+(3) Free Energy of Complex  =");
        PRINT1000(fp, myresults[i].complexG);
        
        fprintf(fp, "DOCKED: USER    (4) Lowest Energy of Unbound System =");
        PRINT1000(fp, lowestUnboundG);
        fprintf(fp, "DOCKED: USER             Enthalpy of Unbound System =");
        PRINT1000(fp, lowestUnboundG_H);
        fprintf(fp, "DOCKED: USER                  -TS of Unbound System =");
        PRINT1000(fp, lowestUnboundG_ntsL);
      
        fprintf(fp, "DOCKED: USER\n");
        if(mypars->xml2dlg || mypars->contact_analysis){
          fprintf(fp, "DOCKED: USER    NEWDPF about 0.0 0.0 0.0\n");
          fprintf(fp, "DOCKED: USER    NEWDPF tran0 %.6f %.6f %.6f\n", myresults[i].genotype[0]*mygrid->spacing, myresults[i].genotype[1]*mygrid->spacing, myresults[i].genotype[2]*mygrid->spacing);
          if(!mypars->xml2dlg){
            double phi = myresults[i].genotype[3]/180.0*PI;
            double theta = myresults[i].genotype[4]/180.0*PI;
            fprintf(fp, "DOCKED: USER    NEWDPF axisangle0 %.8f %.8f %.8f %.6f\n", sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), myresults[i].genotype[5]);
          } else fprintf(fp, "DOCKED: USER    NEWDPF axisangle0 %.8f %.8f %.8f %.6f\n", myresults[i].genotype[3], myresults[i].genotype[4], myresults[i].genotype[5], myresults[i].genotype[genotype_length_in_globmem-1]);
          fprintf(fp, "DOCKED: USER    NEWDPF dihe0");
          for(j=0; j<ligand_ref->num_of_rotbonds; j++)
            fprintf(fp, " %.6f", myresults[i].genotype[6+j]);
          fprintf(fp, "\n");
        }
        fprintf(fp, "DOCKED: USER\n");
        unsigned int lnr=1;
        if ( mypars->flexresfile!=NULL) {
          if ( strlen(mypars->flexresfile)>0 )
            lnr++;
        }
      
        curr_model = pdbqt_template;
        curr_pdbqt = pdbqt_out;
        curr_pdb   = pdb_out;
        // inserting text from the end means prior text positions won't shift
        // so there's less to keep track off ;-)
        if(lineout != NULL)
        {
          printf("ERROR: Bug? The lineout should point to NULL.\n");
          exit(1);
        }
        lineoutSize = 500;
        lineout       = (char *) malloc(lineoutSize);
        if(lineout == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout.\n");
          exit(1);
        }
        char * lineout_pdbqt;
        lineout_pdbqt = (char *) malloc(lineoutSize);
        if(lineout_pdbqt == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout_pdbqt.\n");
          exit(1);
        }
        char * lineout_pdb;
        lineout_pdb = (char *) malloc(lineoutSize);
        if(lineout_pdb == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout_pdb.\n");
          exit(1);
        }
        for(atom_cnt = ligand_ref->num_of_atoms - 1; atom_cnt >= 0; atom_cnt--)
        {
          char* line = lineout;
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][1]); // x
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][2]); // y
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][3]); // z
          line += sprintf(line, "%+6.2lf", copysign(fmin(fabs(myresults[i].peratom_vdw[atom_cnt]),99.99),myresults[i].peratom_vdw[atom_cnt])); // vdw
          line += sprintf(line, "%+6.2lf", copysign(fmin(fabs(myresults[i].peratom_elec[atom_cnt]),99.99),myresults[i].peratom_elec[atom_cnt])); // elec
          line += sprintf(line, "    %+6.3lf ", myresults[i].atom_idxyzq[atom_cnt][4]); // q
          line += sprintf(line, "%-2s\n", ligand_ref->atom_types[((int)myresults[i].atom_idxyzq[atom_cnt][0])]); // type
          curr_model.insert(atom_data[atom_cnt],lineout);
          line = lineout_pdbqt;
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][1]); // x
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][2]); // y
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][3]); // z
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "    %+6.3lf ", myresults[i].atom_idxyzq[atom_cnt][4]); // q
          line += sprintf(line, "%-6s\n", ligand_ref->atom_types[((int)myresults[i].atom_idxyzq[atom_cnt][0])]); // type
          curr_pdbqt.insert(pdbqt_out_data[atom_cnt], lineout_pdbqt);
          line = lineout_pdb;
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][1]); // x
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][2]); // y
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][3]); // z
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "    %+6.3lf ", myresults[i].atom_idxyzq[atom_cnt][4]); // q
          line += sprintf(line, "%-6s\n", ligand_ref->atom_types[((int)myresults[i].atom_idxyzq[atom_cnt][0])]); // type
          curr_pdb.insert(pdb_out_data[atom_cnt], lineout_pdb);
        }
        free(lineout);
        lineout = NULL;
        free(lineout_pdbqt);
        lineout_pdbqt = NULL;
        free(lineout_pdb);
        lineout_pdb = NULL;
        fprintf(fp, "%s", curr_model.c_str());
        bestPdbqt[i] = (char *) malloc(curr_pdbqt.size() + 1);
        if(bestPdbqt[i] == NULL)
        {
          printf("ERROR: Cannot allocate enough memory for bestPdbqt[i].\n");
          exit(1);
        }
        strncpy(bestPdbqt[i], curr_pdbqt.c_str(), curr_pdbqt.size() + 1);
        
        bestPdb[i] = (char *) malloc(curr_pdb.size() + 1);
        if(bestPdb[i] == NULL)
        {
          printf("ERROR: Cannot allocate enough memory for bestPdb[i].\n");
          exit(1);
        }
        strncpy(bestPdb[i], curr_pdb.c_str(), curr_pdb.size() + 1);
        fprintf(fp, "DOCKED: TER\n");
        fprintf(fp, "DOCKED: ENDMDL\n");
        fprintf(fp, "________________________________________________________________________________\n\n\n");
      }
    }
    else if(num_of_runs == 1)
    {
      bestPdbqt = (char **) malloc(mypars->topN4rescore * sizeof(char*));
      if(bestPdbqt == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for bestPdbqt.\n");
        exit(1);
      }
      bestPdb = (char **) malloc(mypars->topN4rescore * sizeof(char*));
      if(bestPdb == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for bestPdb.\n");
        exit(1);
      }
      
      myresults[0].interE     *= mypars->enthalpyScaling;
      myresults[0].interflexE *= mypars->enthalpyScaling;
      myresults[0].intraE     *= mypars->enthalpyScaling;
      myresults[0].intraflexE *= mypars->enthalpyScaling;
      lowestUnboundG_H    = myresults[0].intraE + myresults[0].intraflexE;
      lowestUnboundG_ntsL = myresults[0].ntsL;
      lowestUnboundG      = lowestUnboundG_H + lowestUnboundG_ntsL;
      for(i = 1; i < mypars->topN4rescore; i++)
      {
        myresults[i].interE     *= mypars->enthalpyScaling;
        myresults[i].interflexE *= mypars->enthalpyScaling;
        myresults[i].intraE     *= mypars->enthalpyScaling;
        myresults[i].intraflexE *= mypars->enthalpyScaling;
        float UE = myresults[i].intraE + myresults[i].intraflexE + myresults[i].ntsL;
        if(UE < lowestUnboundG)
        {
          lowestUnboundG_H    = myresults[i].intraE + myresults[i].intraflexE;
          lowestUnboundG_ntsL = myresults[i].ntsL;
          lowestUnboundG      = UE;
        }
      }
      
      for (i=0; i < mypars->topN4rescore; i++)
      {
        fprintf(fp, "    FINAL DOCKED STATE:\n    ________________________\n\n\n");
      
        fprintf(fp, "TopN:   %d / %lu\n", i+1, mypars->topN4rescore);
        fprintf(fp, "Time taken for this:   %.3lfs\n\n", docking_avg_runtime);
        if(mypars->contact_analysis){
          if(myresults[i].analysis.size()>0){
            fprintf(fp, "ANALYSIS: COUNT %lu\n", myresults[i].analysis.size());
            std::string types    = "TYPE    {";
            std::string lig_id   = "LIGID   {";
            std::string ligname  = "LIGNAME {";
            std::string rec_id   = "RECID   {";
            std::string rec_name = "RECNAME {";
            std::string residue  = "RESIDUE {";
            std::string res_id   = "RESID   {";
            std::string chain    = "CHAIN   {";
            char item[8], pad[8];
            for(unsigned int j=0; j<myresults[i].analysis.size(); j++){
              if(j>0){
                types    += ",";
                lig_id   += ",";
                ligname  += ",";
                rec_id   += ",";
                rec_name += ",";
                residue  += ",";
                res_id   += ",";
                chain    += ",";
              }
              switch(myresults[i].analysis[j].type){
                case 0: types += "   \"R\"";
                        break;
                case 1: types += "   \"H\"";
                        break;
                default:
                case 2: types += "   \"V\"";
                        break;
              }
              sprintf(item, "%5d ", myresults[i].analysis[j].lig_id);   lig_id+=item;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].lig_name); sprintf(pad, "%6s", item); ligname+=pad;
              sprintf(item, "%5d ", myresults[i].analysis[j].rec_id);   rec_id+=item;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].rec_name); sprintf(pad, "%6s", item); rec_name+=pad;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].residue); sprintf(pad, "%6s", item);  residue+=pad;
              sprintf(item, "%5d ", myresults[i].analysis[j].res_id);   res_id+=item;
              sprintf(item, "\"%s\"", myresults[i].analysis[j].chain); sprintf(pad, "%6s", item);    chain+=pad;
            }
            fprintf(fp, "ANALYSIS: %s}\n", types.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", lig_id.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", ligname.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", rec_id.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", rec_name.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", residue.c_str());
            fprintf(fp, "ANALYSIS: %s}\n", res_id.c_str());
            fprintf(fp, "ANALYSIS: %s}\n\n", chain.c_str());
          }
        }
      
        fprintf(fp, "DOCKED: MODEL        %d\n", i+1);
        fprintf(fp, "DOCKED: USER    Run = %d\n", i+1);
        fprintf(fp, "DOCKED: USER\n");
      
        fprintf(fp, "DOCKED: USER    Estimated Free Energy of Binding    =");
        
        myresults[i].complexG_H   = myresults[i].interE + myresults[i].interflexE + myresults[i].intraE + myresults[i].intraflexE;
        myresults[i].complexG_nts = myresults[i].nts;
        myresults[i].complexG     = myresults[i].complexG_H + myresults[i].complexG_nts;
        PRINT1000(fp, (myresults[i].complexG - lowestUnboundG));
        fprintf(fp, "DOCKED: USER                                        = (1)+(2)+(3)-(4)\n");
      
        fprintf(fp, "DOCKED: USER\n");
      
        fprintf(fp, "DOCKED: USER    (1) Final Intermolecular Energy     =");
        PRINT1000(fp, ((float) (myresults[i].interE + myresults[i].interflexE)));
      
        fprintf(fp, "DOCKED: USER        vdW + Hbond + desolv Energy     =");
        PRINT1000(fp, ((float) (myresults[i].interE - myresults[i].interE_elec)));
      
        fprintf(fp, "DOCKED: USER        Electrostatic Energy            =");
        PRINT1000(fp, ((float) myresults[i].interE_elec));
      
        fprintf(fp, "DOCKED: USER        Moving Ligand-Fixed Receptor    =");
        PRINT1000(fp, ((float) myresults[i].interE));
      
        fprintf(fp, "DOCKED: USER        Moving Ligand-Moving Receptor   =");
        PRINT1000(fp, ((float) myresults[i].interflexE));
      
        fprintf(fp, "DOCKED: USER    (2) Final Total Internal Energy     =");
        PRINT1000(fp, ((float) (myresults[i].intraE + myresults[i].intraflexE)));
      
        fprintf(fp, "DOCKED: USER    (1)+(2)        Enthalpy of Complex  =");
        PRINT1000(fp, myresults[i].complexG_H);
        fprintf(fp, "DOCKED: USER    (3)                 -TS of Complex  =");
        PRINT1000(fp, myresults[i].complexG_nts);
        fprintf(fp, "DOCKED: USER    (1)+(2)+(3) Free Energy of Complex  =");
        PRINT1000(fp, myresults[i].complexG);
      
        fprintf(fp, "DOCKED: USER    (4) Lowest Energy of Unbound System =");
        PRINT1000(fp, lowestUnboundG);
        fprintf(fp, "DOCKED: USER             Enthalpy of Unbound System =");
        PRINT1000(fp, lowestUnboundG_H);
        fprintf(fp, "DOCKED: USER                  -TS of Unbound System =");
        PRINT1000(fp, lowestUnboundG_ntsL);
      
        fprintf(fp, "DOCKED: USER\n");
        if(mypars->xml2dlg || mypars->contact_analysis){
          fprintf(fp, "DOCKED: USER    NEWDPF about 0.0 0.0 0.0\n");
          fprintf(fp, "DOCKED: USER    NEWDPF tran0 %.6f %.6f %.6f\n", myresults[i].genotype[0]*mygrid->spacing, myresults[i].genotype[1]*mygrid->spacing, myresults[i].genotype[2]*mygrid->spacing);
          if(!mypars->xml2dlg){
            double phi = myresults[i].genotype[3]/180.0*PI;
            double theta = myresults[i].genotype[4]/180.0*PI;
            fprintf(fp, "DOCKED: USER    NEWDPF axisangle0 %.8f %.8f %.8f %.6f\n", sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), myresults[i].genotype[5]);
          } else fprintf(fp, "DOCKED: USER    NEWDPF axisangle0 %.8f %.8f %.8f %.6f\n", myresults[i].genotype[3], myresults[i].genotype[4], myresults[i].genotype[5], myresults[i].genotype[genotype_length_in_globmem-1]);
          fprintf(fp, "DOCKED: USER    NEWDPF dihe0");
          for(j=0; j<ligand_ref->num_of_rotbonds; j++)
            fprintf(fp, " %.6f", myresults[i].genotype[6+j]);
          fprintf(fp, "\n");
        }
        fprintf(fp, "DOCKED: USER\n");
        unsigned int lnr=1;
        if ( mypars->flexresfile!=NULL) {
          if ( strlen(mypars->flexresfile)>0 )
            lnr++;
        }
      
        curr_model = pdbqt_template;
        curr_pdbqt = pdbqt_out;
        curr_pdb   = pdb_out;
        // inserting text from the end means prior text positions won't shift
        // so there's less to keep track off ;-)
        if(lineout != NULL)
        {
          printf("ERROR: Bug? The lineout should point to NULL.\n");
          exit(1);
        }
        lineoutSize = 500;
        lineout       = (char *) malloc(lineoutSize);
        if(lineout == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout.\n");
          exit(1);
        }
        char * lineout_pdbqt;
        lineout_pdbqt = (char *) malloc(lineoutSize);
        if(lineout_pdbqt == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout_pdbqt.\n");
          exit(1);
        }
        char * lineout_pdb;
        lineout_pdb = (char *) malloc(lineoutSize);
        if(lineout_pdb == NULL)
        {
          printf("ERROR: Cannot allocate memory for lineout_pdb.\n");
          exit(1);
        }
        for(atom_cnt = ligand_ref->num_of_atoms - 1; atom_cnt >= 0; atom_cnt--)
        {
          char* line = lineout;
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][1]); // x
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][2]); // y
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][3]); // z
          line += sprintf(line, "%+6.2lf", copysign(fmin(fabs(myresults[i].peratom_vdw[atom_cnt]),99.99),myresults[i].peratom_vdw[atom_cnt])); // vdw
          line += sprintf(line, "%+6.2lf", copysign(fmin(fabs(myresults[i].peratom_elec[atom_cnt]),99.99),myresults[i].peratom_elec[atom_cnt])); // elec
          line += sprintf(line, "    %+6.3lf ", myresults[i].atom_idxyzq[atom_cnt][4]); // q
          line += sprintf(line, "%-2s\n", ligand_ref->atom_types[((int)myresults[i].atom_idxyzq[atom_cnt][0])]); // type
          curr_model.insert(atom_data[atom_cnt],lineout);
          line = lineout_pdbqt;
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][1]); // x
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][2]); // y
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][3]); // z
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "    %+6.3lf ", myresults[i].atom_idxyzq[atom_cnt][4]); // q
          line += sprintf(line, "%-6s\n", ligand_ref->atom_types[((int)myresults[i].atom_idxyzq[atom_cnt][0])]); // type
          curr_pdbqt.insert(pdbqt_out_data[atom_cnt], lineout_pdbqt);
          line = lineout_pdb;
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][1]); // x
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][2]); // y
          line += sprintf(line, "%8.3lf", myresults[i].atom_idxyzq[atom_cnt][3]); // z
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "%6.2lf", 0.0);
          line += sprintf(line, "    %+6.3lf ", myresults[i].atom_idxyzq[atom_cnt][4]); // q
          line += sprintf(line, "%-6s\n", ligand_ref->atom_types[((int)myresults[i].atom_idxyzq[atom_cnt][0])]); // type
          curr_pdb.insert(pdb_out_data[atom_cnt], lineout_pdb);
        }
        free(lineout);
        lineout = NULL;
        free(lineout_pdbqt);
        lineout_pdbqt = NULL;
        free(lineout_pdb);
        lineout_pdb = NULL;
        fprintf(fp, "%s", curr_model.c_str());
        bestPdbqt[i] = (char *) malloc(curr_pdbqt.size() + 1);
        if(bestPdbqt[i] == NULL)
        {
          printf("ERROR: Cannot allocate enough memory for bestPdbqt[i].\n");
          exit(1);
        }
        strncpy(bestPdbqt[i], curr_pdbqt.c_str(), curr_pdbqt.size() + 1);
        
        bestPdb[i] = (char *) malloc(curr_pdb.size() + 1);
        if(bestPdb[i] == NULL)
        {
          printf("ERROR: Cannot allocate enough memory for bestPdb[i].\n");
          exit(1);
        }
        strncpy(bestPdb[i], curr_pdb.c_str(), curr_pdb.size() + 1);
        fprintf(fp, "DOCKED: TER\n");
        fprintf(fp, "DOCKED: ENDMDL\n");
        fprintf(fp, "________________________________________________________________________________\n\n\n");
      }
    }
    else
    {
      printf("ERROR: Bug? num_of_runs (%d) is not larger than zero!?\n", num_of_runs);
      exit(1);
    }
  }

  if(num_of_runs > 1)
  {
    // arranging results according to energy, myresults [energy_order[0]] will be the best one (with lowest energy)
    std::vector<int> energy_order(num_of_runs);
    std::vector<double> energies(num_of_runs);
    
    for (i=0; i<num_of_runs; i++){
      energy_order[i] = i;
      energies[i] = myresults[i].complexG - lowestUnboundG;
      myresults[i].clus_id = 0; // indicates that it hasn't been put into cluster yet (may as well do that here ...)
    }
    // sorting the indices instead of copying the results around will be faster
    for(i=0; i<num_of_runs-1; i++)
      for(j=0; j<num_of_runs-i-1; j++)
        if(energies[energy_order[j]]>energies[energy_order[j+1]]) // swap indices to percolate larger energies up
          std::swap(energy_order[j], energy_order[j+1]);
    // PERFORM CLUSTERING
    if(mypars->calc_clustering){
    
      // the best result is the center of the first cluster
      myresults[energy_order[0]].clus_id = 1;
      myresults[energy_order[0]].rmsd_from_cluscent = 0;
      num_of_clusters = 1;
    
      for (int w=1; w<num_of_runs; w++) // for each result
      {
        i=energy_order[w];
        current_clust_center = 0;
        result_clustered = 0;
    
        for (int u=0; u<w; u++) // results with lower id-s are clustered, look for cluster centers
        {
          j=energy_order[u];
          if (myresults[j].clus_id > current_clust_center) // it is the center of a new cluster
          {
            current_clust_center = myresults[j].clus_id;
            temp_rmsd = calc_rmsd(myresults[j].atom_idxyzq, myresults[i].atom_idxyzq, ligand_ref->true_ligand_atoms, mypars->handle_symmetry); // comparing current result with cluster center
            if (temp_rmsd <= cluster_tolerance) // in this case we put result i to cluster with center j
            {
              myresults[i].clus_id = current_clust_center;
              myresults[i].rmsd_from_cluscent = temp_rmsd;
              result_clustered = 1;
              break;
            }
          }
        }
    
        if (result_clustered != 1) // if no suitable cluster was found, this is the center of a new one
        {
          num_of_clusters++;
          myresults[i].clus_id = num_of_clusters; // new cluster id
          myresults[i].rmsd_from_cluscent = 0;
        }
      }
    
      for (i=1; i<=num_of_clusters; i++) // printing cluster info to file
      {
        subrank = 0;
        cluster_sizes [i-1] = 0;
        sum_energy [i-1] = 0;
        for (int u=0; u<num_of_runs; u++){
          j = energy_order[u];
          if (myresults [j].clus_id == i)
          {
            subrank++;
            cluster_sizes[i-1]++;
            sum_energy [i-1] += myresults[j].complexG - lowestUnboundG;
            myresults[j].clus_subrank = subrank;
            if (subrank == 1)
            {
              best_energy [i-1] = myresults[j].complexG - lowestUnboundG;
              best_energy_runid  [i-1] = myresults[j].run_number;
            }
          }
        }
      }
    
      if(mypars->output_dlg){
        // WRITING CLUSTER INFORMATION
        fprintf(fp, "    CLUSTERING HISTOGRAM\n    ____________________\n\n\n");
        fprintf(fp, "________________________________________________________________________________\n");
        fprintf(fp, "     |           |     |           |     |\n");
        fprintf(fp, "Clus | Lowest    | Run | Mean      | Num | Histogram\n");
        fprintf(fp, "-ter | Binding   |     | Binding   | in  |\n");
        fprintf(fp, "Rank | Energy    |     | Energy    | Clus|    5    10   15   20   25   30   35\n");
        fprintf(fp, "_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___\n");
    
        for (i=0; i<num_of_clusters; i++)
        {
          if(mypars->gen_best > i)
          {
            FILE* fp_best;
    
            sprintf(fileName4bestN, "%s_best_CR%d_run%d.pdbqt", mypars->resname, i + 1, best_energy_runid[i]);
            fp_best = fopen(fileName4bestN, "w");
            if (fp_best == NULL)
            {
              printf("ERROR: Can't open file %s.\n", fileName4bestN);
              exit(1);
            }
            fprintf(fp_best, "%s", bestPdbqt[best_energy_runid[i] - 1]);
            fclose(fp_best);
            
            sprintf(fileName4bestN, "%s_best_CR%d_run%d.pdb", mypars->resname, i + 1, best_energy_runid[i]);
            fp_best = fopen(fileName4bestN, "w");
            if (fp_best == NULL)
            {
              printf("ERROR: Can't open file %s.\n", fileName4bestN);
              exit(1);
            }
            fprintf(fp_best, "%s", bestPdb[best_energy_runid[i] - 1]);
            fclose(fp_best);
          }
          fprintf(fp, "%4d |", i+1);
    
          if (best_energy[i] > 999999.99)
            fprintf(fp, "%+10.2e", best_energy[i]);
          else
            fprintf(fp, "%+10.2f", best_energy[i]);
          fprintf(fp, " |%4d |", best_energy_runid[i]);
    
          if (sum_energy[i]/cluster_sizes[i] > 999999.99)
            fprintf(fp, "%+10.2e |", sum_energy[i]/cluster_sizes[i]);
          else
            fprintf(fp, "%+10.2f |", sum_energy[i]/cluster_sizes[i]);
    
          fprintf(fp, "%4d |", cluster_sizes [i]);
    
          for (j=0; j<cluster_sizes [i]; j++)
            fprintf(fp, "#");
    
          fprintf(fp, "\n");
        }
    
        fprintf(fp, "_____|___________|_____|___________|_____|______________________________________\n\n\n");
    
        // writing RMSD table
    
        fprintf(fp, "    RMSD TABLE\n");
        fprintf(fp, "    __________\n\n\n");
    
        fprintf(fp, "_______________________________________________________________________\n");
        fprintf(fp, "     |      |      |           |         |                 |\n");
        fprintf(fp, "Rank | Sub- | Run  | Binding   | Cluster | Reference       | Grep\n");
        fprintf(fp, "     | Rank |      | Energy    | RMSD    | RMSD            | Pattern\n");
        fprintf(fp, "_____|______|______|___________|_________|_________________|___________\n" );
    
        for (i=0; i<num_of_clusters; i++) // printing cluster info to file
        {
          for (int u=0; u<num_of_runs; u++){
            j = energy_order[u];
            if (myresults[j].clus_id == i+1) {
              if (myresults[j].complexG - lowestUnboundG > 999999.99)
                fprintf(fp, "%4d   %4d   %4d  %+10.2e  %8.2f  %8.2f           RANKING\n",
                             myresults[j].clus_id,
                                   myresults[j].clus_subrank,
                                         myresults[j].run_number,
                                              myresults[j].complexG - lowestUnboundG,
                                                             myresults[j].rmsd_from_cluscent,
                                                                 myresults[j].rmsd_from_ref);
              else
                fprintf(fp, "%4d   %4d   %4d  %+10.2f  %8.2f  %8.2f           RANKING\n",
                             myresults[j].clus_id,
                                   myresults[j].clus_subrank,
                                         myresults[j].run_number,
                                              myresults[j].complexG - lowestUnboundG,
                                                             myresults[j].rmsd_from_cluscent,
                                                                 myresults[j].rmsd_from_ref);
            }
          }
        }
      }
    }
    for(int i = 0; i < num_of_runs; i++)
    {
      free(bestPdbqt[i]);
      free(bestPdb[i]);
    }
    free(bestPdbqt);
    free(bestPdb);
    
    if(mypars->output_dlg){
      // Add execution and idle time information
      fprintf(fp, "\nRun time %.3f sec", exec_time);
      fprintf(fp, "\nIdle time %.3f sec\n", idle_time);
      if(!mypars->dlg2stdout){
        fclose(fp);
      }
    }
    
    // if xml has to be generated
    if (mypars->output_xml)
    {
      char* xml_file_name = (char*)malloc(len*sizeof(char));
      strcpy(xml_file_name, mypars->resname);
      strcat(xml_file_name, ".xml");
      fp_xml = fopen(xml_file_name, "w");
      if(fp==NULL){
        printf("Error: Cannot create xml output file %s: %s\n",xml_file_name,strerror(errno));
        exit(9);
      }
      fprintf(fp_xml, "<?xml version=\"1.0\" ?>\n");
      fprintf(fp_xml, "<omegadocker>\n");
      fprintf(fp_xml, "\t<version>%s</version>\n",VERSION);
      if((*argc)>1){
        fprintf(fp_xml, "\t<arguments>");
        for(i=1; i<(*argc); i++){
          fprintf(fp_xml, "%s%s", (i>1)?" ":"", argv[i]);
          if (argcmp("filelist", argv[i], 'B')){
            if(mypars->filelist_files>1){
              fprintf(fp_xml, " %s", mypars->ligandfile);
              i+=mypars->filelist_files; // skip ahead in case there are multiple entries here
            }
          }
          if (argcmp("xml2dlg", argv[i], 'X')){
            if(mypars->xml_files>1){
              fprintf(fp_xml, " %s", mypars->load_xml);
              i+=mypars->xml_files; // skip ahead in case there are multiple entries here
            }
          }
        }
        fprintf(fp_xml, "</arguments>\n");
      }
      if(mypars->dpffile)
        fprintf(fp_xml, "\t<dpf>%s</dpf>\n",mypars->dpffile);
      if(mypars->list_nr>1)
        fprintf(fp_xml, "\t<list_nr>%u</list_nr>\n",mypars->list_nr);
      fprintf(fp_xml, "\t<grid>%s</grid>\n", mypars->fldfile);
      fprintf(fp_xml, "\t<ligand>%s</ligand>\n", mypars->ligandfile);
      if(mypars->flexresfile)
        fprintf(fp_xml, "\t<flexres>%s</flexres>\n",mypars->flexresfile);
      fprintf(fp_xml, "\t<seed>");
      if(!mypars->seed[2]){
        if(!mypars->seed[1]){
          fprintf(fp_xml,"%d", mypars->seed[0]);
        } else fprintf(fp_xml,"%d %d", mypars->seed[0], mypars->seed[1]);
      } else fprintf(fp_xml,"%d %d %d", mypars->seed[0], mypars->seed[1], mypars->seed[2]);
      fprintf(fp_xml, "</seed>\n");
      fprintf(fp_xml, "\t<ls_method>%s</ls_method>\n",mypars->ls_method);
      fprintf(fp_xml, "\t<autostop>%s</autostop>\n",mypars->autostop ? "yes" : "no");
      fprintf(fp_xml, "\t<heuristics>%s</heuristics>\n",mypars->use_heuristics ? "yes" : "no");
      fprintf(fp_xml, "\t<run_requested>%lu</run_requested>\n",mypars->num_of_runs);
      fprintf(fp_xml, "\t<runs>\n");
      double phi, theta;
      for(int u=0; u<num_of_runs; u++){
        j = energy_order[u];
        fprintf(fp_xml, "\t\t<run id=\"%d\">\n",(myresults [j]).run_number);
        if(mypars->contact_analysis){
          if(myresults[j].analysis.size()>0){
            fprintf(fp_xml, "\t\t\t<contact_analysis count=\"%lu\">\n", myresults[j].analysis.size());
            std::string types;
            std::string lig_id;
            std::string ligname;
            std::string rec_id;
            std::string rec_name;
            std::string residue;
            std::string res_id;
            std::string chain;
            char item[8], pad[8];
            for(unsigned int i=0; i<myresults[j].analysis.size(); i++){
              if(i>0){
                types    += ",";
                lig_id   += ",";
                ligname  += ",";
                rec_id   += ",";
                rec_name += ",";
                residue  += ",";
                res_id   += ",";
                chain    += ",";
              }
              switch(myresults[j].analysis[i].type){
                case 0: types += "   \"R\"";
                        break;
                case 1: types += "   \"H\"";
                        break;
                default:
                case 2: types += "   \"V\"";
                        break;
              }
              sprintf(item, "%5d ", myresults[j].analysis[i].lig_id);   lig_id+=item;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].lig_name); sprintf(pad, "%6s", item); ligname+=pad;
              sprintf(item, "%5d ", myresults[j].analysis[i].rec_id);   rec_id+=item;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].rec_name); sprintf(pad, "%6s", item); rec_name+=pad;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].residue); sprintf(pad, "%6s", item);  residue+=pad;
              sprintf(item, "%5d ", myresults[j].analysis[i].res_id);   res_id+=item;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].chain); sprintf(pad, "%6s", item);    chain+=pad;
            }
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_types>  %s</contact_analysis_types>\n", types.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_ligid>  %s</contact_analysis_ligid>\n", lig_id.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_ligname>%s</contact_analsyis_ligname>\n", ligname.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_recid>  %s</contact_analysis_recid>\n", rec_id.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_recname>%s</contact_analysis_recname>\n", rec_name.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_residue>%s</contact_analysis_residue>\n", residue.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_resid>  %s</contact_analysis_resid>\n", res_id.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_chain>  %s</contact_analysis_chain>\n", chain.c_str());
            fprintf(fp_xml, "\t\t\t</contact_analysis>\n");
          }
        }
        fprintf(fp_xml, "\t\t\t<free_NRG_binding>%.2f</free_NRG_binding>\n", myresults[j].complexG - lowestUnboundG);
        fprintf(fp_xml, "\t\t\t<final_intermol_NRG> %.2f</final_intermol_NRG>\n", myresults[j].interE + myresults[j].interflexE);
        fprintf(fp_xml, "\t\t\t<internal_ligand_NRG>%.2f</internal_ligand_NRG>\n", myresults[j].intraE + myresults[j].intraflexE);
        fprintf(fp_xml, "\t\t\t<torsonial_free_NRG> %.2f</torsonial_free_NRG>\n", myresults[j].nts);
        fprintf(fp_xml, "\t\t\t<tran0>%.6f %.6f %.6f</tran0>\n", myresults[j].genotype[0]*mygrid->spacing, myresults[j].genotype[1]*mygrid->spacing, myresults[j].genotype[2]*mygrid->spacing);
        phi = myresults[j].genotype[3]/180.0*PI;
        theta = myresults[j].genotype[4]/180.0*PI;
        fprintf(fp_xml, "\t\t\t<axisangle0>%.8f %.8f %.8f %.6f</axisangle0>\n", sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), myresults[j].genotype[5]);
        fprintf(fp_xml, "\t\t\t<ndihe>%d</ndihe>\n", ligand_ref->num_of_rotbonds);
        fprintf(fp_xml, "\t\t\t<dihe0>");
        for(i=0; i<ligand_ref->num_of_rotbonds; i++)
          fprintf(fp_xml, "%s%.6f", (i>0)?" ":"", myresults[j].genotype[6+i]);
        fprintf(fp_xml, "\n\t\t\t</dihe0>\n");
        fprintf(fp_xml, "\t\t</run>\n");
      }
      fprintf(fp_xml, "\t</runs>\n");
      if(mypars->calc_clustering){
        fprintf(fp_xml, "\t<result>\n");
        fprintf(fp_xml, "\t\t<clustering_histogram>\n");
        for (i=0; i<num_of_clusters; i++)
        {
          fprintf(fp_xml, "\t\t\t<cluster cluster_rank=\"%d\" lowest_binding_energy=\"%.2lf\" run=\"%d\" mean_binding_energy=\"%.2lf\" num_in_clus=\"%d\" />\n",
            i+1, best_energy[i], best_energy_runid[i], sum_energy[i]/cluster_sizes[i], cluster_sizes [i]);
        }
        fprintf(fp_xml, "\t\t</clustering_histogram>\n");
    
        fprintf(fp_xml, "\t\t<rmsd_table>\n");
        for (i=0; i<num_of_clusters; i++)
        {
          for (int u=0; u<num_of_runs; u++){
            j = energy_order[u];
            if (myresults[j].clus_id == i+1)
            {
              fprintf(fp_xml, "\t\t\t<run rank=\"%d\" sub_rank=\"%d\" run=\"%d\" binding_energy=\"%.2lf\" cluster_rmsd=\"%.2lf\" reference_rmsd=\"%.2lf\" />\n",
                myresults[j].clus_id, myresults[j].clus_subrank, myresults[j].run_number, myresults[j].complexG - lowestUnboundG, myresults[j].rmsd_from_cluscent, myresults[j].rmsd_from_ref);
            }
          }
        }
        fprintf(fp_xml, "\t\t</rmsd_table>\n");
        fprintf(fp_xml, "\t</result>\n");
      }
      fprintf(fp_xml, "</omegadocker>\n");
      fclose(fp_xml);
      free(xml_file_name);
    }
  }
  else if(num_of_runs == 1)
  {
    // arranging results according to energy, myresults [energy_order[0]] will be the best one (with lowest energy)
    std::vector<int> energy_order(mypars->topN4rescore);
    std::vector<double> energies(mypars->topN4rescore);
    
    for(i = 0; i < mypars->topN4rescore; i++){
      energy_order[i] = i;
      energies[i] = myresults[i].complexG - lowestUnboundG;
      myresults[i].clus_id = 0; // indicates that it hasn't been put into cluster yet (may as well do that here ...)
    }
    // sorting the indices instead of copying the results around will be faster
    for(i=0; i<mypars->topN4rescore-1; i++)
    {
      for(j=0; j<mypars->topN4rescore-i-1; j++)
      {
        if(energies[energy_order[j]]>energies[energy_order[j+1]]) // swap indices to percolate larger energies up
        {
          std::swap(energy_order[j], energy_order[j+1]);
        }
      }
    }
    // PERFORM CLUSTERING
    if(mypars->calc_clustering){
      // the best result is the center of the first cluster
      myresults[energy_order[0]].clus_id = 1;
      myresults[energy_order[0]].rmsd_from_cluscent = 0;
      num_of_clusters = 1;
    
      for (int w=1; w<mypars->topN4rescore; w++) // for each result
      {
        i=energy_order[w];
        current_clust_center = 0;
        result_clustered = 0;
    
        for (int u=0; u<w; u++) // results with lower id-s are clustered, look for cluster centers
        {
          j=energy_order[u];
          if (myresults[j].clus_id > current_clust_center) // it is the center of a new cluster
          {
            current_clust_center = myresults[j].clus_id;
            temp_rmsd = calc_rmsd(myresults[j].atom_idxyzq, myresults[i].atom_idxyzq, ligand_ref->true_ligand_atoms, mypars->handle_symmetry); // comparing current result with cluster center
            if (temp_rmsd <= cluster_tolerance) // in this case we put result i to cluster with center j
            {
              myresults[i].clus_id = current_clust_center;
              myresults[i].rmsd_from_cluscent = temp_rmsd;
              result_clustered = 1;
              break;
            }
          }
        }
    
        if (result_clustered != 1) // if no suitable cluster was found, this is the center of a new one
        {
          num_of_clusters++;
          myresults[i].clus_id = num_of_clusters; // new cluster id
          myresults[i].rmsd_from_cluscent = 0;
        }
      }
    
      for (i=1; i<=num_of_clusters; i++) // printing cluster info to file
      {
        subrank = 0;
        cluster_sizes [i-1] = 0;
        sum_energy [i-1] = 0;
        for (int u=0; u<mypars->topN4rescore; u++){
          j = energy_order[u];
          if (myresults [j].clus_id == i)
          {
            subrank++;
            cluster_sizes[i-1]++;
            sum_energy [i-1] += myresults[j].complexG - lowestUnboundG;
            myresults[j].clus_subrank = subrank;
            if (subrank == 1)
            {
              best_energy [i-1] = myresults[j].complexG - lowestUnboundG;
              best_energy_runid  [i-1] = j + 1;
            }
          }
        }
      }
    
      if(mypars->output_dlg){
        // WRITING CLUSTER INFORMATION
        fprintf(fp, "    CLUSTERING HISTOGRAM\n    ____________________\n\n\n");
        fprintf(fp, "_________________________________________________________________________________\n");
        fprintf(fp, "     |           |      |           |     |\n");
        fprintf(fp, "Clus | Lowest    | Ori. | Mean      | Num | Histogram\n");
        fprintf(fp, "-ter | Binding   |      | Binding   | in  |\n");
        fprintf(fp, "Rank | Energy    | Rank | Energy    | Clus|    5    10   15   20   25   30   35\n");
        fprintf(fp, "_____|___________|______|___________|_____|____:____|____:____|____:____|____:___\n");
    
        for (i=0; i<num_of_clusters; i++)
        {
          if(mypars->gen_best > i)
          {
            FILE* fp_best;
    
            sprintf(fileName4bestN, "%s_best_CR%d_oriRank%d.pdbqt", mypars->resname, i + 1, best_energy_runid[i]);
            fp_best = fopen(fileName4bestN, "w");
            if (fp_best == NULL)
            {
              printf("ERROR: Can't open file %s.\n", fileName4bestN);
              exit(1);
            }
            fprintf(fp_best, "%s", bestPdbqt[best_energy_runid[i] - 1]);
            fclose(fp_best);
            
            sprintf(fileName4bestN, "%s_best_CR%d_oriRank%d.pdb", mypars->resname, i + 1, best_energy_runid[i]);
            fp_best = fopen(fileName4bestN, "w");
            if (fp_best == NULL)
            {
              printf("ERROR: Can't open file %s.\n", fileName4bestN);
              exit(1);
            }
            fprintf(fp_best, "%s", bestPdb[best_energy_runid[i] - 1]);
            fclose(fp_best);
          }
          fprintf(fp, "%4d |", i+1);
    
          if (best_energy[i] > 999999.99)
            fprintf(fp, "%+10.2e", best_energy[i]);
          else
            fprintf(fp, "%+10.2f", best_energy[i]);
          fprintf(fp, " |%5d |", best_energy_runid[i]);
    
          if (sum_energy[i]/cluster_sizes[i] > 999999.99)
            fprintf(fp, "%+10.2e |", sum_energy[i]/cluster_sizes[i]);
          else
            fprintf(fp, "%+10.2f |", sum_energy[i]/cluster_sizes[i]);
    
          fprintf(fp, "%4d |", cluster_sizes [i]);
    
          for (j=0; j<cluster_sizes [i]; j++)
            fprintf(fp, "#");
    
          fprintf(fp, "\n");
        }
    
        fprintf(fp, "_____|___________|______|___________|_____|______________________________________\n\n\n");
    
        // writing RMSD table
    
        fprintf(fp, "    RMSD TABLE\n");
        fprintf(fp, "    __________\n\n\n");
    
        fprintf(fp, "_______________________________________________________________________\n");
        fprintf(fp, "     |      |      |           |         |                 |\n");
        fprintf(fp, "Rank | Sub- | Ori. | Binding   | Cluster | Reference       | Grep\n");
        fprintf(fp, "     | Rank | Rank | Energy    | RMSD    | RMSD            | Pattern\n");
        fprintf(fp, "_____|______|______|___________|_________|_________________|___________\n" );
    
        for (i=0; i<num_of_clusters; i++) // printing cluster info to file
        {
          for (int u=0; u<mypars->topN4rescore; u++){
            j = energy_order[u];
            if (myresults[j].clus_id == i+1) {
              if (myresults[j].complexG - lowestUnboundG > 999999.99)
                fprintf(fp, "%4d   %4d   %4d  %+10.2e  %8.2f  %8.2f           RANKING\n",
                             myresults[j].clus_id,
                                   myresults[j].clus_subrank,
                                           j + 1,
                                              myresults[j].complexG - lowestUnboundG,
                                                             myresults[j].rmsd_from_cluscent,
                                                                 myresults[j].rmsd_from_ref);
              else
                fprintf(fp, "%4d   %4d   %4d  %+10.2f  %8.2f  %8.2f           RANKING\n",
                             myresults[j].clus_id,
                                   myresults[j].clus_subrank,
                                           j + 1,
                                              myresults[j].complexG - lowestUnboundG,
                                                             myresults[j].rmsd_from_cluscent,
                                                                 myresults[j].rmsd_from_ref);
            }
          }
        }
      }
    }
    for(int i = 0; i < mypars->topN4rescore; i++)
    {
      free(bestPdbqt[i]);
      free(bestPdb[i]);
    }
    free(bestPdbqt);
    free(bestPdb);
    
    if(mypars->output_dlg){
      // Add execution and idle time information
      fprintf(fp, "\nRun time %.3f sec", exec_time);
      fprintf(fp, "\nIdle time %.3f sec\n", idle_time);
      if(!mypars->dlg2stdout){
        fclose(fp);
      }
    }
    
    // if xml has to be generated
    if (mypars->output_xml)
    {
      char* xml_file_name = (char*)malloc(len*sizeof(char));
      strcpy(xml_file_name, mypars->resname);
      strcat(xml_file_name, ".xml");
      fp_xml = fopen(xml_file_name, "w");
      if(fp==NULL){
        printf("Error: Cannot create xml output file %s: %s\n",xml_file_name,strerror(errno));
        exit(9);
      }
      fprintf(fp_xml, "<?xml version=\"1.0\" ?>\n");
      fprintf(fp_xml, "<omegadocker>\n");
      fprintf(fp_xml, "\t<version>%s</version>\n",VERSION);
      if((*argc)>1){
        fprintf(fp_xml, "\t<arguments>");
        for(i=1; i<(*argc); i++){
          fprintf(fp_xml, "%s%s", (i>1)?" ":"", argv[i]);
          if (argcmp("filelist", argv[i], 'B')){
            if(mypars->filelist_files>1){
              fprintf(fp_xml, " %s", mypars->ligandfile);
              i+=mypars->filelist_files; // skip ahead in case there are multiple entries here
            }
          }
          if (argcmp("xml2dlg", argv[i], 'X')){
            if(mypars->xml_files>1){
              fprintf(fp_xml, " %s", mypars->load_xml);
              i+=mypars->xml_files; // skip ahead in case there are multiple entries here
            }
          }
        }
        fprintf(fp_xml, "</arguments>\n");
      }
      if(mypars->dpffile)
        fprintf(fp_xml, "\t<dpf>%s</dpf>\n",mypars->dpffile);
      if(mypars->list_nr>1)
        fprintf(fp_xml, "\t<list_nr>%u</list_nr>\n",mypars->list_nr);
      fprintf(fp_xml, "\t<grid>%s</grid>\n", mypars->fldfile);
      fprintf(fp_xml, "\t<ligand>%s</ligand>\n", mypars->ligandfile);
      if(mypars->flexresfile)
        fprintf(fp_xml, "\t<flexres>%s</flexres>\n",mypars->flexresfile);
      fprintf(fp_xml, "\t<seed>");
      if(!mypars->seed[2]){
        if(!mypars->seed[1]){
          fprintf(fp_xml,"%d", mypars->seed[0]);
        } else fprintf(fp_xml,"%d %d", mypars->seed[0], mypars->seed[1]);
      } else fprintf(fp_xml,"%d %d %d", mypars->seed[0], mypars->seed[1], mypars->seed[2]);
      fprintf(fp_xml, "</seed>\n");
      fprintf(fp_xml, "\t<ls_method>%s</ls_method>\n",mypars->ls_method);
      fprintf(fp_xml, "\t<autostop>%s</autostop>\n",mypars->autostop ? "yes" : "no");
      fprintf(fp_xml, "\t<heuristics>%s</heuristics>\n",mypars->use_heuristics ? "yes" : "no");
      fprintf(fp_xml, "\t<run_requested>%lu</run_requested>\n",mypars->num_of_runs);
      fprintf(fp_xml, "\t<runs>\n");
      double phi, theta;
      for(int u=0; u<mypars->topN4rescore; u++){
        j = energy_order[u];
        fprintf(fp_xml, "\t\t<run id=\"%d\">\n",(myresults [j]).run_number);
        if(mypars->contact_analysis){
          if(myresults[j].analysis.size()>0){
            fprintf(fp_xml, "\t\t\t<contact_analysis count=\"%lu\">\n", myresults[j].analysis.size());
            std::string types;
            std::string lig_id;
            std::string ligname;
            std::string rec_id;
            std::string rec_name;
            std::string residue;
            std::string res_id;
            std::string chain;
            char item[8], pad[8];
            for(unsigned int i=0; i<myresults[j].analysis.size(); i++){
              if(i>0){
                types    += ",";
                lig_id   += ",";
                ligname  += ",";
                rec_id   += ",";
                rec_name += ",";
                residue  += ",";
                res_id   += ",";
                chain    += ",";
              }
              switch(myresults[j].analysis[i].type){
                case 0: types += "   \"R\"";
                        break;
                case 1: types += "   \"H\"";
                        break;
                default:
                case 2: types += "   \"V\"";
                        break;
              }
              sprintf(item, "%5d ", myresults[j].analysis[i].lig_id);   lig_id+=item;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].lig_name); sprintf(pad, "%6s", item); ligname+=pad;
              sprintf(item, "%5d ", myresults[j].analysis[i].rec_id);   rec_id+=item;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].rec_name); sprintf(pad, "%6s", item); rec_name+=pad;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].residue); sprintf(pad, "%6s", item);  residue+=pad;
              sprintf(item, "%5d ", myresults[j].analysis[i].res_id);   res_id+=item;
              sprintf(item, "\"%s\"", myresults[j].analysis[i].chain); sprintf(pad, "%6s", item);    chain+=pad;
            }
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_types>  %s</contact_analysis_types>\n", types.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_ligid>  %s</contact_analysis_ligid>\n", lig_id.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_ligname>%s</contact_analsyis_ligname>\n", ligname.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_recid>  %s</contact_analysis_recid>\n", rec_id.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_recname>%s</contact_analysis_recname>\n", rec_name.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_residue>%s</contact_analysis_residue>\n", residue.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_resid>  %s</contact_analysis_resid>\n", res_id.c_str());
            fprintf(fp_xml, "\t\t\t\t<contact_analysis_chain>  %s</contact_analysis_chain>\n", chain.c_str());
            fprintf(fp_xml, "\t\t\t</contact_analysis>\n");
          }
        }
        fprintf(fp_xml, "\t\t\t<free_NRG_binding>   %.2f</free_NRG_binding>\n", myresults[j].complexG - lowestUnboundG);
        fprintf(fp_xml, "\t\t\t<final_intermol_NRG> %.2f</final_intermol_NRG>\n", myresults[j].interE + myresults[j].interflexE);
        fprintf(fp_xml, "\t\t\t<internal_ligand_NRG>%.2f</internal_ligand_NRG>\n", myresults[j].intraE + myresults[j].intraflexE);
        fprintf(fp_xml, "\t\t\t<torsonial_free_NRG> %.2f</torsonial_free_NRG>\n", myresults[j].nts);
        fprintf(fp_xml, "\t\t\t<tran0>%.6f %.6f %.6f</tran0>\n", myresults[j].genotype[0]*mygrid->spacing, myresults[j].genotype[1]*mygrid->spacing, myresults[j].genotype[2]*mygrid->spacing);
        phi = myresults[j].genotype[3]/180.0*PI;
        theta = myresults[j].genotype[4]/180.0*PI;
        fprintf(fp_xml, "\t\t\t<axisangle0>%.8f %.8f %.8f %.6f</axisangle0>\n", sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta), myresults[j].genotype[5]);
        fprintf(fp_xml, "\t\t\t<ndihe>%d</ndihe>\n", ligand_ref->num_of_rotbonds);
        fprintf(fp_xml, "\t\t\t<dihe0>");
        for(i=0; i<ligand_ref->num_of_rotbonds; i++)
          fprintf(fp_xml, "%s%.6f", (i>0)?" ":"", myresults[j].genotype[6+i]);
        fprintf(fp_xml, "\n\t\t\t</dihe0>\n");
        fprintf(fp_xml, "\t\t</run>\n");
      }
      fprintf(fp_xml, "\t</runs>\n");
      if(mypars->calc_clustering){
        fprintf(fp_xml, "\t<result>\n");
        fprintf(fp_xml, "\t\t<clustering_histogram>\n");
        for (i=0; i<num_of_clusters; i++)
        {
          fprintf(fp_xml, "\t\t\t<cluster cluster_rank=\"%d\" lowest_binding_energy=\"%.2lf\" run=\"%d\" mean_binding_energy=\"%.2lf\" num_in_clus=\"%d\" />\n",
            i+1, best_energy[i], best_energy_runid[i], sum_energy[i]/cluster_sizes[i], cluster_sizes [i]);
        }
        fprintf(fp_xml, "\t\t</clustering_histogram>\n");
    
        fprintf(fp_xml, "\t\t<rmsd_table>\n");
        for (i=0; i<num_of_clusters; i++)
        {
          for (int u=0; u<mypars->topN4rescore; u++){
            j = energy_order[u];
            if (myresults[j].clus_id == i+1)
            {
              fprintf(fp_xml, "\t\t\t<run rank=\"%d\" sub_rank=\"%d\" run=\"%d\" binding_energy=\"%.2lf\" cluster_rmsd=\"%.2lf\" reference_rmsd=\"%.2lf\" />\n",
                myresults[j].clus_id, myresults[j].clus_subrank, best_energy_runid[i], myresults[j].complexG - lowestUnboundG, myresults[j].rmsd_from_cluscent, myresults[j].rmsd_from_ref);
            }
          }
        }
        fprintf(fp_xml, "\t\t</rmsd_table>\n");
        fprintf(fp_xml, "\t</result>\n");
      }
      fprintf(fp_xml, "</omegadocker>\n");
      fclose(fp_xml);
      free(xml_file_name);
    }
  }
  else
  {
    printf("ERROR: Bug? num_of_runs (%d) is not larger than zero!?\n", num_of_runs);
    exit(1);
  }
  
  free(cluster_sizes);
  free(sum_energy);
  free(best_energy);
  free(best_energy_runid);
}

int qsortCompFloatInc(const void *a, const void *b)
{
  float fa = *(const float*) a;
  float fb = *(const float*) b;
  return (fa > fb) - (fa < fb);
}

void calcEntropy
(
  Ligandresult     *ligResult,
  const Gridinfo   *mygrid,
  const Dockpars   *mypars,
  const Liganddata *lig
)
{
  int shakingGridRadiusInNum; // this is not radius, but the half of (shakingGridEdgeLenInNum - 1)
  int shakingGridEdgeLenInNum;
  int shakingGridNum;
  int shakingGridCenterIndex;
  float centerE;
  double atomRandomness;
  double atomRandomnessLigOnly;
  double molRandomness;
  double molRandomnessLigOnly;
  
  shakingGridRadiusInNum  = (int) ceil(mypars->shakingRadius / mygrid->spacing);
  shakingGridEdgeLenInNum = shakingGridRadiusInNum * 2 + 1;
  shakingGridNum          = shakingGridEdgeLenInNum * shakingGridEdgeLenInNum * shakingGridEdgeLenInNum;
  shakingGridCenterIndex  = (shakingGridNum - 1) / 2;
  
  ligResult->nts  = 0.0;
  ligResult->ntsL = 0.0;
  
  for(int atomIndex = 0; atomIndex < lig->num_of_atoms; atomIndex++)
  {
    atomRandomness        = 0.0;
    atomRandomnessLigOnly = 0.0;
    molRandomness         = 1.0;
    molRandomnessLigOnly  = 1.0;
    centerE      = ligResult->shakingGridE[atomIndex][shakingGridCenterIndex];
    qsort(ligResult->shakingGridE[atomIndex], shakingGridNum, sizeof(float), qsortCompFloatInc);
    qsort(ligResult->shakingGridEL[atomIndex], shakingGridNum, sizeof(float), qsortCompFloatInc);
    centerE -= ligResult->shakingGridE[atomIndex][0];
    if(mypars->warningOn && centerE >= mypars->wellDepth)
    {
      printf("WARNING: The central potential energy (%f kcal/mol) of atom index %d is higher then the well depth (%f kcal/mol).\n", centerE, atomIndex, mypars->wellDepth);
    }
    
    for(int shakingGridIndex = shakingGridNum - 1; shakingGridIndex > -1; shakingGridIndex--)
    {
      ligResult->shakingGridE[atomIndex][shakingGridIndex]  -= ligResult->shakingGridE[atomIndex][0];
      ligResult->shakingGridEL[atomIndex][shakingGridIndex] -= ligResult->shakingGridEL[atomIndex][0];
    }
    
    for(int shakingGridIndex = 0; shakingGridIndex < shakingGridNum; shakingGridIndex++)
    {
      if(ligResult->shakingGridE[atomIndex][shakingGridIndex] < mypars->wellDepth)
      {
        atomRandomness += 1.0 - ((double) ligResult->shakingGridE[atomIndex][shakingGridIndex]) / mypars->wellDepth;
      }
      else
      {
        break;
      }
    }
    molRandomness *= atomRandomness;
    if(molRandomness > 1E290)
    {
      ligResult->nts += log(molRandomness);
      molRandomness = 1.0;
    }
    
    for(int shakingGridIndex = 0; shakingGridIndex < shakingGridNum; shakingGridIndex++)
    {
      if(ligResult->shakingGridEL[atomIndex][shakingGridIndex] < mypars->wellDepth)
      {
        atomRandomnessLigOnly += 1.0 - ((double) ligResult->shakingGridEL[atomIndex][shakingGridIndex]) / mypars->wellDepth;
      }
      else
      {
        break;
      }
    }
    molRandomnessLigOnly *= atomRandomnessLigOnly;
    if(molRandomnessLigOnly > 1E290)
    {
      ligResult->ntsL += log(molRandomnessLigOnly);
      molRandomnessLigOnly = 1.0;
    }
  }
  
  ligResult->nts  += log(molRandomness);
  ligResult->nts  *= mypars->entropyScaling;
  ligResult->ntsL += log(molRandomnessLigOnly);
  ligResult->ntsL *= mypars->entropyScaling;
}

void process_result(
                    const Gridinfo*        mygrid,
                    const Dockpars*        mypars,
                          Liganddata*      myligand_init,
                    const Liganddata*      myxrayligand,
                    const int*             argc,
                          char**           argv,
                          SimulationState& sim_state
                   )
{
  Ligandresult *cpu_result_ligands;
  if(mypars->num_of_runs > 1)
  {
    cpu_result_ligands = (Ligandresult *) malloc(mypars->num_of_runs * sizeof(Ligandresult));
    if(cpu_result_ligands == NULL)
    {
      printf("ERROR: Cannot allocat memory for cpu_result_ligands\n");
      exit(1);
    }
    for(long unsigned int i = 0; i < mypars->num_of_runs; i++)
    {
      cpu_result_ligands[i].genotype = NULL;
      cpu_result_ligands[i].atom_idxyzq = (double **) calloc(myligand_init->num_of_atoms, sizeof(double*));
      if(cpu_result_ligands[i].atom_idxyzq == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].atom_idxyzq.\n");
        exit(1);
      }
      cpu_result_ligands[i].shakingGridE = (float **) calloc(myligand_init->num_of_atoms, sizeof(float*));
      if(cpu_result_ligands[i].shakingGridE == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].shakingGridE.\n");
        exit(1);
      }
      cpu_result_ligands[i].shakingGridEL = (float **) calloc(myligand_init->num_of_atoms, sizeof(float*));
      if(cpu_result_ligands[i].shakingGridEL == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].shakingGridEL.\n");
        exit(1);
      }
      for(int j = 0; j < myligand_init->num_of_atoms; j++)
      {
        cpu_result_ligands[i].atom_idxyzq[j] = (double *) calloc(5, sizeof(double));
        if(cpu_result_ligands[i].atom_idxyzq[j] == NULL)
        {
          printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].atom_idxyzq[j].\n");
          exit(1);
        }
        
        cpu_result_ligands[i].shakingGridE[i] = NULL;
        cpu_result_ligands[i].shakingGridEL[i] = NULL;
      }
    
      cpu_result_ligands[i].peratom_vdw = (float *) calloc(myligand_init->num_of_atoms, sizeof(float));
      if(cpu_result_ligands[i].peratom_vdw == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].peratom_vdw.\n");
        exit(1);
      }
    
      cpu_result_ligands[i].peratom_elec = (float *) calloc(myligand_init->num_of_atoms, sizeof(float));
      if(cpu_result_ligands[i].peratom_elec == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].peratom_elec.\n");
        exit(1);
      }
    }
  }
  else if(mypars->num_of_runs == 1)
  {
    cpu_result_ligands = (Ligandresult *) malloc(mypars->topN4rescore * sizeof(Ligandresult));
    if(cpu_result_ligands == NULL)
    {
      printf("ERROR: Cannot allocat memory for cpu_result_ligands\n");
      exit(1);
    }
    
    for(long unsigned int i = 0; i < mypars->topN4rescore; i++)
    {
      cpu_result_ligands[i].genotype = NULL;
      cpu_result_ligands[i].atom_idxyzq = (double **) calloc(myligand_init->num_of_atoms, sizeof(double*));
      if(cpu_result_ligands[i].atom_idxyzq == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].atom_idxyzq.\n");
        exit(1);
      }
      cpu_result_ligands[i].shakingGridE = (float **) calloc(myligand_init->num_of_atoms, sizeof(float*));
      if(cpu_result_ligands[i].shakingGridE == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].shakingGridE.\n");
        exit(1);
      }
      cpu_result_ligands[i].shakingGridEL = (float **) calloc(myligand_init->num_of_atoms, sizeof(float*));
      if(cpu_result_ligands[i].shakingGridEL == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].shakingGridEL.\n");
        exit(1);
      }
      for(int j = 0; j < myligand_init->num_of_atoms; j++)
      {
        cpu_result_ligands[i].atom_idxyzq[j] = (double *) calloc(5, sizeof(double));
        if(cpu_result_ligands[i].atom_idxyzq[j] == NULL)
        {
          printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].atom_idxyzq[j].\n");
          exit(1);
        }
        
        cpu_result_ligands[i].shakingGridE[j] = NULL;
        cpu_result_ligands[i].shakingGridEL[j] = NULL;
      }
      
      cpu_result_ligands[i].peratom_vdw = (float *) calloc(myligand_init->num_of_atoms, sizeof(float));
      if(cpu_result_ligands[i].peratom_vdw == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].peratom_vdw.\n");
        exit(1);
      }
    
      cpu_result_ligands[i].peratom_elec = (float *) calloc(myligand_init->num_of_atoms, sizeof(float));
      if(cpu_result_ligands[i].peratom_elec == NULL)
      {
        printf("ERROR: Cannot allocate enough memory for cpu_result_ligands[i].peratom_elec.\n");
        exit(1);
      }
      
      cpu_result_ligands[i].interE_elec = 0.0;
    }
  }
  else
  {
    printf("ERROR: Bug? mypars->num_of_runs (%d) is not larger than zero!?\n", mypars->num_of_runs);
    exit(1);
  }
  
  // Fill in cpu_result_ligands
  float best_energy_of_all = 1000000000000.0;
  IntraTables *tables= new IntraTables(sim_state.myligand_reference, mypars->coeffs.scaled_AD4_coeff_elec, mypars->coeffs.AD4_coeff_desolv, mypars->qasp, mypars->nr_mod_atype_pairs, mypars->mod_atype_pairs);
  if(mypars->num_of_runs > 1)
  {
    for (unsigned long run_cnt=0; run_cnt < mypars->num_of_runs; run_cnt++)
    {
      make_resfiles(sim_state.cpu_populations[run_cnt],
                    sim_state.cpu_energies[run_cnt],
                    tables,
                    sim_state.myligand_reference,
                    myligand_init,
                    myxrayligand,
                    mypars,
                    sim_state.cpu_evals_of_runs[run_cnt],
                    sim_state.generation_cnt,
                    mygrid,
                    argc,
                    argv,
                    0,
                    run_cnt,
                    best_energy_of_all,
                    &(cpu_result_ligands [run_cnt]));
    }
  }
  else if(mypars->num_of_runs == 1)
  {
      make_resfiles(sim_state.cpu_populations[0],
                    sim_state.cpu_energies[0],
                    tables,
                    sim_state.myligand_reference,
                    myligand_init,
                    myxrayligand,
                    mypars,
                    sim_state.cpu_evals_of_runs[0],
                    sim_state.generation_cnt,
                    mygrid,
                    argc,
                    argv,
                    0,
                    0,
                    best_energy_of_all,
                    cpu_result_ligands);
  }
  else
  {
    printf("ERROR: Bug? mypars->num_of_runs (%d) is not larger than zero!?\n", mypars->num_of_runs);
    exit(1);
  }
  mpz_t Evaluations;
  mpz_init(Evaluations);
  mpz_fdiv_q_ui(Evaluations, sim_state.total_evals, mypars->num_of_runs);
  generate_output(cpu_result_ligands,
                  mypars->num_of_runs,
                  tables,
                  myligand_init,
                  myxrayligand,
                  mypars,
                  mygrid,
                  argc,
                  argv,
                  sim_state.sec_per_run,
                  sim_state.generation_cnt,
                  Evaluations,
                  sim_state.exec_time,
                  sim_state.idle_time);
  mpz_clear(Evaluations);
  delete tables;
  if(mypars->num_of_runs > 1)
  {
    for(long unsigned int i = 0; i < mypars->num_of_runs; i++)
    {
      for(int j = 0; j < myligand_init->num_of_atoms; j++)
      {
        free(cpu_result_ligands[i].atom_idxyzq[j]);
        free(cpu_result_ligands[i].shakingGridE[j]);
        free(cpu_result_ligands[i].shakingGridEL[j]);
      }
      free(cpu_result_ligands[i].atom_idxyzq);
      free(cpu_result_ligands[i].shakingGridE);
      free(cpu_result_ligands[i].shakingGridEL);
      free(cpu_result_ligands[i].peratom_vdw);
      free(cpu_result_ligands[i].peratom_elec);
    }
  }
  else if(mypars->num_of_runs == 1)
  {
    for(long unsigned int i = 0; i < mypars->topN4rescore; i++)
    {
      for(int j = 0; j < myligand_init->num_of_atoms; j++)
      {
        free(cpu_result_ligands[i].atom_idxyzq[j]);
        free(cpu_result_ligands[i].shakingGridE[j]);
        free(cpu_result_ligands[i].shakingGridEL[j]);
      }
      free(cpu_result_ligands[i].atom_idxyzq);
      free(cpu_result_ligands[i].shakingGridE);
      free(cpu_result_ligands[i].shakingGridEL);
      free(cpu_result_ligands[i].peratom_vdw);
      free(cpu_result_ligands[i].peratom_elec);
    }
  }
  else
  {
    printf("ERROR: Bug? mypars->num_of_runs (%d) is not larger than zero!?\n", mypars->num_of_runs);
    exit(1);
  }
  free(cpu_result_ligands);
}
