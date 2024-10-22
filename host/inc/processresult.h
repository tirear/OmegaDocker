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


#ifndef PROCESSRESULT_H_
#define PROCESSRESULT_H_

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "processligand.h"
#include "getparameters.h"
#include "simulation_state.hpp"

#include <gmp.h>

#define PRINT1000(file, x) fprintf(file,  ((fabs((x)) >= 0.0) && ((fabs(x)) <= 1000.)) ? "%+7.2f kcal/mol (%+7.2f kj/mol)\n" : "%+11.2e kcal/mol (%+11.2e kj/mol)\n", (x), (KCAL_TO_KJ * x));

typedef struct
{
  float*                    genotype; // a pointer here is sufficient and saves lots of memory copies
  double                  **atom_idxyzq;
  float                   **shakingGridE;  // shakingGridE[atomNum][shakingGridNum]
  float                   **shakingGridEL; // shakingGridE[atomNum][shakingGridNum] for ligand solo only
  double                    nts;  // -TS
  double                    ntsL; // -TS for ligand solo
  float                     complexG;
  float                     complexG_H;
  float                     complexG_nts;
  float                     interE;
  float                     interflexE;
  float                     interE_elec;
  float                     intraE;
  float                     intraflexE;
  float                    *peratom_vdw;
  float                    *peratom_elec;
  float                     rmsd_from_ref;
  float                     rmsd_from_cluscent;
  int                       clus_id;
  int                       clus_subrank;
  int                       run_number;
  std::vector<AnalysisData> analysis;
} Ligandresult;


void arrange_result(
                          float* final_population,
                          float* energies,
                    Liganddata *myligand_init,
                    const int    pop_size
                   );

void write_basic_info(
                            FILE*       fp,
                            Liganddata* ligand_ref,
                      const Dockpars*   mypars,
                      const Gridinfo*   mygrid,
                      const int*        argc,
                      char**            argv
                     );

void write_basic_info_dlg(
                                FILE*       fp,
                                Liganddata* ligand_ref,
                          const Dockpars*   mypars,
                          const Gridinfo*   mygrid,
                          const int*        argc,
                                char**      argv
                         );

void make_resfiles(
                         float*        final_population,
                         float*        energies,
                         IntraTables*  tables,
                         Liganddata*   ligand_ref,
                         Liganddata*   ligand_from_pdb,
                   const Liganddata*   ligand_xray,
                   const Dockpars*     mypars,
                         int           evals_performed,
                         int           generations_used,
                   const Gridinfo*     mygrid,
                   const int*          argc,
                         char**        argv,
                         int           debug,
                         int           run_cnt,
                         float&        best_energy_of_all,
                         Ligandresult* best_result
                  );

void ligand_calc_output(
                              FILE*         fp,
                        const char*         prefix,
                              IntraTables*  tables,
                        const Liganddata*   ligand,
                        const Dockpars*     mypars,
                        const Gridinfo*     mygrid,
                              bool          output_analysis,
                              bool          output_energy
                       );

void generate_output(
                           Ligandresult  myresults [],
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
                    );
int qsortCompFloatInc(const float *a, const float *b);

void calcEntropy
(
  Ligandresult     *ligResult,
  const Gridinfo   *mygrid,
  const Dockpars   *mypars,
  const Liganddata *lig
);

void process_result(
                    const Gridinfo*        mygrid,
                    const Dockpars*        mypars,
                          Liganddata*      myligand_init,
                    const Liganddata*      myxrayligand,
                    const int*             argc,
                          char**           argv,
                          SimulationState& sim_state
                   );

#endif /* PROCESSRESULT_H_ */
