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


#ifndef FILELIST_HPP
#define FILELIST_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "processgrid.h"

typedef struct _Dockpars Dockpars;

typedef struct _fld_files_data{
  std::string name;
  size_t grid_idx;
} fld_files_data;

class FileList{
  public:
    bool                        used;
    int                         nfiles;
    bool                        preload_maps;
    bool                        maps_are_loaded;
    char*                       filename;
    std::vector<std::string>    resnames;
    std::vector<fld_files_data> fld_files;
    std::vector<std::string>    ligand_files;
    std::vector<Dockpars>       mypars;
    std::vector<Gridinfo>       mygrids;
    std::vector<bool>           load_maps_gpu; // indicate which device needs to still load maps from cpu

    // Default to unused, with 1 file
    FileList() : used( false ), nfiles( 1 ), preload_maps( true ), maps_are_loaded( false ), filename( NULL ) {}
    ~FileList(){ if(filename) free(filename); }
};

#endif
