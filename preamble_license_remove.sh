#!/bin/bash
# Remove license preamble from all OmegaDocker source and header files
# if license preamble is present

# Reference
# sed '/PATTERNSTART/,/PATTERNEND/{0,/PATTERNEND/d}' file

# license preamble
LICENSE_PREAMBLE="./preamble_license"
if (grep -q "Copyright (C)" $LICENSE_PREAMBLE)
then
  echo "The content of preamble_license is as follows:"
  cat $LICENSE_PREAMBLE
  echo ""
  echo ""
else
  echo "ERROR: There is something wrong with $LICENSE_PREAMBLE?!"
  exit
fi

# kernel-header files
KRNL_HEADER_DIR="./common"
KRNL_HEADERS="$KRNL_HEADER_DIR/*.h"

# CU kernel-source2 files
CUKRNL_SOURCE_DIR="./cuda"
CUKRNL_SOURCE="$CUKRNL_SOURCE_DIR/*.cu $CUKRNL_SOURCE_DIR/*.h"

# CL kernel-source files
CLKRNL_SOURCE_DIR="./device"
CLKRNL_SOURCE="$CLKRNL_SOURCE_DIR/*.cl $CLKRNL_SOURCE_DIR/*.h"

# host-header files
HOST_HEADER_DIR="./host/inc"
HOST_HEADERS="$HOST_HEADER_DIR/*.h $HOST_HEADER_DIR/*.hpp $HOST_HEADER_DIR/*.Cuda $HOST_HEADER_DIR/*.OpenCL"

# host-source files
HOST_SOURCE_DIR="./host/src"
HOST_SOURCES="$HOST_SOURCE_DIR/*.cpp $HOST_SOURCE_DIR/*.Cuda $HOST_SOURCE_DIR/*.OpenCL"

# wrapcl-header files
WRAPCL_HEADER_DIR="./wrapcl/inc"
WRAPCL_HEADERS="$WRAPCL_HEADER_DIR/*.h"

# wrapcl-source files
WRAPCL_SOURCE_DIR="./wrapcl/src"
WRAPCL_SOURCES="$WRAPCL_SOURCE_DIR/*.cpp"

# full list of source files
OMEGADOCKER_SOURCE="$KRNL_HEADERS $CUKRNL_SOURCE $CLKRNL_SOURCE $HOST_HEADERS $HOST_SOURCES $WRAPCL_HEADERS $WRAPCL_SOURCES" # sans

# Print variables
#echo $KRNL_HEADERS
#echo $CUKRNL_SOURCE
#echo $CLKRNL_SOURCE
#echo $HOST_HEADERS
#echo $HOST_SOURCES
#echo $WRAPCL_HEADERS
#echo $WRAPCL_SOURCES
#echo $OMEGADOCKER_SOURCE

# Remove license-preamble
# Excluding sources that do not have it, and
# excluding the automatically-generated ./host/inc/stringify.h
rm -f host/inc/performdocking.h
rm -f host/src/performdocking.cpp
for f in $OMEGADOCKER_SOURCE
do
  if [ "$f" != "$HOST_HEADER_DIR/stringify.h" ]
  then
    if (grep -q "Copyright (C)" $f)
    then
      echo "Removing existing license-preamble from $f"
      sed '/\/\*/,/\*\//{0,/\*\//d}' "$f" > "$f.old"
      mv "$f.old" "$f"
      echo "Done!"
    else
      echo "License-preamble was not found in $f"
      echo "No license-preamble is removed."
    fi
    echo " "
  fi
done
