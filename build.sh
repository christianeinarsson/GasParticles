#!/bin/bash
# Usage: build [-n y] [-N x]
#
# help_build() To print help
#
help_build()
{
  echo "Usage: build [-n y] [-N x]"
  echo "Options: These are optional argument"
  echo " -n initial number of particles in a node (default 500)"
  echo " -N maximum number of particles in a node (default 2*n)"
  exit 1
}

### Start main procedure
# Set default value for variable

n=500
N=`expr $n \* 2`

# if no argument
if [ $# -lt 1 ]; then
	help_build
fi
while getopts n:N: opt
do
	case "$opt" in
		n) n="$OPTARG";;
		N) N="$OPTARG";;
		\?) help_build;;
	esac
done

icc *.c -Nmpi -lVT -I$VT_ROOT/include -L$VT_LIB_DIR $VT_ADD_LIBS -DINIT_NO_PARTICLES=$n -DMAX_NO_PARTICLES=$N -o gaslaw_${n}_$N
