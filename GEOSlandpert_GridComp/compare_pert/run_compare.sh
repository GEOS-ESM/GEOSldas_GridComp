#!/bin/bash
# run_compare.sh -- build-free driver: timing sweep + visualization.
#
# Usage (from a run directory, with GEOS modules loaded):
#   run_compare.sh <bindir> [nlon] [nlat] [ndraw] ["ranks"] ["xcorrs"]
# Example:
#   run_compare.sh ../install/bin 2880 1440 5 "1 2 4 8 16" "0.2 0.5 1 2 5"
set -e
BIN=${1:?bindir}
NLON=${2:-2880}; NLAT=${3:-1440}; NDRAW=${4:-10}
RANKS=${5:-"1 2 4 8 16"}
XCORRS=${6:-"0.2 0.5 1 2 5"}
MPIRUN=${MPIRUN:-mpirun}

rm -f compare_pert_timing.csv
for xc in $XCORRS; do
  for np in $RANKS; do
    $MPIRUN -np $np "$BIN/compare_pert.x" $xc $NDRAW $NLON $NLAT
  done
  # keep one map/correlation figure per xcorr
  python3 "$BIN/compare_pert.py" "compare_pert_xcorr$xc"
done
python3 "$BIN/plot_timing.py" compare_pert_timing.csv compare_pert_timing
