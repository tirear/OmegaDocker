#!/bin/bash
# Show license preamble from all OmegaDocker source and header files
# if license preamble is present

# Reference
# sed '/PATTERNSTART/,/PATTERNEND/{0,/PATTERNEND/d}' file

# license preamble
LICENSE_PREAMBLE="./preamble_license"

# kernel-header files
KRNL_HEADER_DIR="./common"
KRNL_HEADERS="$KRNL_HEADER_DIR/*.h"

# CU kernel-source files
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
OMEGADOCKER_SOURCE="$KRNL_HEADERS $CUKRNL_SOURCE $CLKRNL_SOURCE $HOST_HEADERS $HOST_SOURCES $WRAPCL_HEADERS $WRAPCL_SOURCES"

# Print variables
#echo $KRNL_HEADERS
#echo $CUKRNL_SOURCE
#echo $CLKRNL_SOURCE
#echo $HOST_HEADERS
#echo $HOST_SOURCES
#echo $WRAPCL_HEADERS
#echo $WRAPCL_SOURCES
#echo $OMEGADOCKER_SOURCE

# Show license-preamble
# Excluding sources that do not have it, and
# excluding the automatically-generated ./host/inc/stringify.h
rm -f host/inc/performdocking.h
rm -f host/src/performdocking.cpp
for f in $OMEGADOCKER_SOURCE
do
  if [ "$f" != "$HOST_HEADER_DIR/stringify.h" ]
  then
    if (grep -q "Lesser General Public License" $f)
    then
      echo "LesserG $f"
    elif (grep -q "General Public License" $f)
    then
      echo "G $f"
    else
      echo "N $f"
    fi
    echo " "
  fi
done
