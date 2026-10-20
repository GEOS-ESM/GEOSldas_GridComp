#!/bin/bash
# run_xcorr_stats.sh -- generate perturbations for several xcorr values and
# compute/plot the xcorr recovered from them (see XCORR_STATS.md).
#
# Usage (from an empty run directory):
#   module load ifort-stack
#   run_xcorr_stats.sh <bindir> [nlon] [nlat] [nranks] ["xcorrs"]
# Example:
#   run_xcorr_stats.sh ../install/bin 2880 1440 8 "0.2 0.4 0.8 1.0"
set -e
BIN=${1:?bindir}
NLON=${2:-2880}; NLAT=${3:-1440}; NP=${4:-8}
XCORRS=${5:-"0.2 0.4 0.8 1.0"}
MPIRUN=${MPIRUN:-mpirun}

for xc in $XCORRS; do
  mkdir -p x$xc
  ( cd x$xc && $MPIRUN -np $NP "$BIN/compare_pert.x" $xc 1 $NLON $NLAT )
done

export PYTHONPATH="$BIN:$PYTHONPATH"
python3 "$BIN/xcorr_stats.py"      $NLON $NLAT $XCORRS | tee xcorr_stats.txt
python3 "$BIN/pole_check.py"       $NLON $NLAT $XCORRS | tee pole_check.txt
python3 "$BIN/plot_xcorr_stats.py" $NLON $NLAT $XCORRS
