<center><img src="OmegaDocker_logo_01.png" width="300"></center>
OmegaDocker: an OligoMEr-supported and Gpu-Accelerated Docking program
======================================================================

# About

* OmegaDocker is a GPU-accelerated docking program developed for larger ligands, such as oligopeptides, oligonucleotides, and other oligomers.
* The first version of OmegaDocker is based on AutoDock-GPU v1.5.4.
* AutoDock-GPU is developed by the [Forli lab](https://forlilab.org/) at Scripps Research.

# CAUTION
This is still an alpha version.

# Citation

NOT YET

# Features

* The restrictions on the maximum number of atoms to 256 atoms and related parameters are removed.
* Cuda accelerated only ... for now.
* Modified searching algorithem
* Modified scoring function

# Compilation

The first step is to set environmental variables `GPU_INCLUDE_PATH` and `GPU_LIBRARY_PATH`,
as described here: https://github.com/ccsb-scripps/AutoDock-GPU/wiki/Guideline-for-users

```zsh
make DEVICE=<TYPE> NUMWI=<NWI>
```

| Parameters | Description                  | Values                                             |
|:----------:|:----------------------------:|:--------------------------------------------------:|
| `<TYPE>`   | Accelerator chosen           | `CUDA`                                             |
| `<NWI>`    | work-group/thread block size | `128`, `256`                                       |

# Usage

## Basic command
```zsh
./bin/autodock_<type>_<N>wi \
--ffile <protein>.maps.fld \
--lfile <ligand>.pdbqt \
--nrun <nruns>
```

| Mandatory options|   | Description   | Value                     |
|:----------------:|:-:|:-------------:|:-------------------------:|
|--ffile           |-M |Protein file   |&lt;protein&gt;.maps.fld   |
|--lfile           |-L |Ligand file    |&lt;ligand&gt;.pdbqt       |

Both options can alternatively be provided in the contents of the files specified with `--filelist (-B)` (see below for format) and `--import_dpf (-I)` (AD4 dpf file format).

## Example
```zsh
./bin/autodock_gpu_64wi \
--ffile ./input/1stp/derived/1stp_protein.maps.fld \
--lfile ./input/1stp/derived/1stp_ligand.pdbqt
```
By default the output log file is written in the current working folder. Examples of output logs can be found under [examples/output](examples/output/).

## Supported arguments
| Argument          |   | Description                                           | Default value    |
|:------------------|:--|:------------------------------------------------------|-----------------:|
|<tr><td colspan="4">**INPUT**</td></tr>
|--lfile            |-L | Ligand pdbqt file                                     | no default       |
|--lfile            |-L | Ligand pdbqt file                                     | no default       |
|--ffile            |-M | Grid map files descriptor fld file                    | no default       |
|--flexres          |-F | Flexible residue pdbqt file                           | no default       |
|--filelist         |-B | Batch file                                            | no default       |
|--import_dpf       |-I | Import AD4-type dpf input file (only partial support) | no default       |
|--xraylfile        |-R | reference ligand file for RMSD analysis               | ligand file      <tr><td colspan="4">**CONVERSION**</td></tr>
|--xml2dlg          |-X | One (or many) AD-GPU xml file(s) to convert to dlg(s) | no default       |
|--stateFileInName  |   | Name for the input docking state file.                | stateIn.bin      <tr><td colspan="4">**OUTPUT**</td></tr>
|--resnam           |-N | Name for docking output log                           | ligand basename  |
|--contact_analysis |-C | Perform distance-based analysis (description below)   | 0 (no)           |
|--xmloutput        |-x | Specify if xml output format is wanted                | 1 (yes)          |
|--dlgoutput        |-d | Control if dlg output is created                      | 1 (yes)          |
|--dlg2stdout       |-2 | Write dlg file output to stdout (if not OVERLAP=ON)   | 0 (no)           |
|--rlige            |   | Print reference ligand energies                       | 0 (no)           |
|--gfpop            |   | Output all poses from all populations of each LGA run | 0 (no)           |
|--npdb             |   | # pose pdbqt files from populations of each LGA run   | 0                |
|--gbest            |   | Number of output pdbqt files of the best pose of each run (up to the number of run) | 0                |
|--miniCoorOut      |   | The coordinates in the minimization stage will be     | 0 (no)           |
|                   |   | written to standard output if the value is not zero.  |                  |
|                   |   | If this value is not equal to 0, then the following   |                  |
|                   |   | is set: --psize 2 --nrun 1 --heurmax 10.              |                  |
|--stateFileOutPre  |   | Prefix for the output docking state file.             | stateOut         |
|--showTopRunNum    |   | The number of top runs of that the energy be showed.  | 0                |
|--showTopIdvNum    |   | The number of top individuals of that the energy be   | 0                |
|                   |   | showed.                                               |                  |
|--cautionOn        |   | Show the caution information.                         | 0 (no)           |
|--warningOn        |   | Show the warning information.                         | 0 (no)           |
|--clustering       |   | Output clustering analysis in dlg and/or xml file     | 1 (yes)          |
|--hsym             |   | Handle symmetry in RMSD calc.                         | 1 (yes)          |
|--rmstol           |   | RMSD clustering tolerance                             | 2 (Å)            <tr><td colspan="4">**SETUP**</td></tr>
|--devnum           |-D | OpenCL/Cuda device number (counting starts at 1)      | 1                |
|--loadxml          |-c | Load initial population from xml results file         | no default       |
|--seed             |-s | Random number seeds (up to three comma-sep. integers) | time, process id <tr><td colspan="4">**SEARCH**</td></tr>
|--heuristics       |-H | Ligand-based automatic search method and # evals      | 1 (yes)          |
|--heurmax          |-E | Asymptotic heuristics # evals limit (smooth limit)    | 12000000         |
|--autostop         |-A | Automatic stopping criterion based on convergence     | 1 (yes)          |
|--asfreq           |-a | AutoStop testing frequency (in # of generations)      | 5                |
|--nrun             |-n | # LGA runs                                            | 20               |
|--nev              |-e | # Score evaluations (max.) per LGA run                | 2500000          |
|--ngen             |-g | # Generations (max.) per LGA run                      | 42000            |
|--lsmet            |-l | Local-search method                                   | ad (ADADELTA)    |
|--lsit             |-i | # Local-search iterations (max.)                      | 300              |
|--psize            |-p | Population size                                       | 150              |
|--inputConfRatio   |   | The ratio of the population that the conformation     | 0.2              |
|                   |   | genes are initialized as file input.                  |                  |
|--mrat             |   | Mutation rate                                         | 2   (%)          |
|--crat             |   | Crossover rate                                        | 80  (%)          |
|--lsrat            |   | Local-search rate                                     | 100 (%)          |
|--trat             |   | Tournament (selection) rate                           | 60  (%)          |
|--dmov             |   | Maximum LGA movement delta                            | 6 (Å)            |
|--dang             |   | Maximum LGA angle delta                               | 90 (°)           |
|--rholb            |   | Solis-Wets lower bound of rho parameter               | 0.01             |
|--lsmov            |   | Solis-Wets movement delta                             | 2 (Å)            |
|--lsang            |   | Solis-Wets angle delta                                | 75 (°)           |
|--cslim            |   | Solis-Wets cons. success/failure limit to adjust rho  | 4                |
|--stopstd          |   | AutoStop energy standard deviation tolerance          | 0.15 (kcal/mol)  |
|--initswgens       |   | Initial # generations of Solis-Wets instead of -lsmet | 0 (no)           |
|--ligaStepSize     |   | The step size for ligand only minimization.           | 0.1              |
|--compStepSize     |   | The step size for complex minimization.               | 0.1              |
|--ligaRhoLimit     |   | The maximal rho value for ligand minimization. Larger | 100              |
|                   |   | rho, more steps for ligand minimization.              |                  |
|--compRhoLimit     |   | The maximal rho value for complex minimization.       | 1000             |
|                   |   | Larger rho, more steps for complex minimization.      |                  |
|--collRhoLimit     |   | The maximal rho value for collision-eliminating       | 5                |
|                   |   | minimization. Larger rho, more steps for              |                  |
|                   |   | collision-eliminating minimization.                   |                  |
|--collDistance     |   | The distance at which non-bonded atoms are considered | 1.3              |
|                   |   | to have collided (exclude H-H).                       |                  |
|--rebornRatio      |   | The ratio of population that is treated as bad        | 0.2              |
|                   |   | samples and will be regenerated.                      |                  |
|--miniOnly         |   | Do minimization only. If this value is not equal to   | 0 (no)           |
|                   |   | 0, then the following is set: --miniCoorOut 1.        |                  |
|--addonEng         |   | The amount of energy that needs to be added to set    | 400              |
|                   |   | the threshold of low energy.                          |                  |
|--lowEngNumRatio   |   | The maximal ratio of the low energy number in the     | 0.1              |
|                   |   | population.                                           |                  |
|--convergedThre    |   | The threshold of energy difference to define the      | 0.01             |
|                   |   | convergence of the energy used for comparison.        |                  |
|--ligandOnly       |   | Find the lowest energy conformation of ligand alone.  | 0 (no)           <tr><td colspan="4">**SCORING**</td></tr>
|--derivtype        |-T | Derivative atom types (e.g. C1,C2,C3=C/S4=S/H5=HD)    | no default       |
|--modpair          |-P | Modify vdW pair params (e.g. C1:S4,1.60,1.200,13,7)   | no default       |
|--ubmod            |-u | Unbound model: 0 (bound), 1 (extended), 2 (compact)   | 0 (same as bound)|
|--smooth           |   | Smoothing parameter for vdW interactions              | 0.5 (Å)          |
|--elecmindist      |   | Min. electrostatic potential distance (w/ dpf: 0.5 Å) | 0.01 (Å)         |
|--modqp            |   | Use modified QASP from VirtualDrug or AD4 original    | 0 (no, use AD4)  |
|--savetime         |   | Do not re-calculate final energy with CPU as the the  | 0 (no)           |
|                   |   | number of run is more than one. (0/1)                 |                  |
|--topN4rescore     |   | The number of top individuals of that the energy be   | 1                |
|                   |   | re-calculated as the number of run is one.            |                  |
|--enthalpyScaling  |   | The value for enthalpy scaling.                       | 0.7564           |
|--entropyScaling   |   | The value for entropy scaling.                        | 0.9              |
|--shakingRadius    |   | The radius of atom shaking for calculate entropy (in  | 0.75             |
|                   |   | Å).                                                   |                  |
|--wellDepth        |   | The depth of a energy landscape for calculating       | 0.35             |
|                   |   | entropy (in kcal/mol).                                |                  |

# Documentation

# Contributing
