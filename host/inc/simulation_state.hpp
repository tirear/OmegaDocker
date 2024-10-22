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


#ifndef SIMULATION_STATE_HPP
#define SIMULATION_STATE_HPP

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "profile.hpp"

#include "processgrid.h"
#include "miscellaneous.h"
#include "processligand.h"
#include "getparameters.h"
#include "calcenergy.h"
#include "processresult.h"

struct SimulationState{
  float             **cpu_populations;
  float             **cpu_energies;
  Liganddata        *myligand_reference;
  int               *cpu_evals_of_runs;
  mpz_t              generation_cnt;
  float              sec_per_run;
  mpz_t              total_evals;
  double             exec_time;
  double             idle_time;
};

#endif
