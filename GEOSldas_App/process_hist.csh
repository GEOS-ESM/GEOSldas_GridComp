#!/bin/csh -f

# process GEOSldas_HIST.rc (=$HISTRC) template
#
# - turn on/off HIST collections depending on tile space
#   - EASE:      turn on tile-space (1d) output
#   - otherwise: turn on gridded    (2d) output
# - turn on/off output variables depending on LSM_CHOICE, AEROSOL_DEPOSITION, and RUN_IRRIG
# - fill in source 'GridComp' info for variables depending on NENS

# process command line args

setenv HISTRC             $1         # file name of HIST rc template (GEOSldas_HIST.rc)
setenv OUTxd              $2         # "OUT1d" or "OUT2d" (to turn on/off collections)
setenv GRIDNAME         "'$3'"       # full name of grid associated with tile space
setenv LSM_CHOICE         $4
setenv AEROSOL_DEPOSITION $5
setenv RUN_IRRIG          $6
setenv NENS               $7

# -------------------------------------------------

echo $GRIDNAME

# uncomment 2d or 1d collections, depending on "OUT1d" (EASE tile space) or "OUT2d" (non-EASE tile space)

if($OUTxd == OUT1d) then
   sed -i 's|\#OUT1d|''|g' $HISTRC
else
   sed -i 's|\#OUT2d|''|g' $HISTRC
endif

# fill in name of grid associated with tile space

sed -i -e  s/\'GRIDNAME\'/$GRIDNAME/g $HISTRC

# set 'GridComp' based on LSM_CHOICE;
# turn on/off variables associated with CATCHCN, AEROSOL_DEPOSITION, RUN_IRRIG

if($LSM_CHOICE == 1) then
   set GridComp = CATCH
   sed -i '/^>>>HIST_CATCHCN<<</d'         $HISTRC
   sed -i '/^>>>HIST_CATCHCNCLM45<<</d'    $HISTRC
endif

if($LSM_CHOICE == 2) then
   set GridComp = CATCHCN
   sed -i '/^>>>HIST_CATCHCNCLM45<<</d'    $HISTRC
   sed -i 's/>>>HIST_CATCHCN<<</''/g'      $HISTRC
endif

if($LSM_CHOICE == 3) then
   set GridComp = CATCHCN
   sed -i 's/>>>HIST_CATCHCN<<</''/g'      $HISTRC
   sed -i 's/>>>HIST_CATCHCNCLM45<<</''/g' $HISTRC
endif

if($AEROSOL_DEPOSITION == 0) then
   sed -i '/^>>>HIST_AEROSOL<<</d'         $HISTRC
else
   sed -i 's/>>>HIST_AEROSOL<<</''/g'      $HISTRC
endif

if($RUN_IRRIG == 0) then
   sed -i '/^>>>HIST_IRRIG<<</d'           $HISTRC
else
   sed -i 's/>>>HIST_IRRIG<<</''/g'        $HISTRC
endif

# for ensemble simulations, set 'GridComp' to ENSAVG 

if($NENS > 1) then
   set GridComp = ENSAVG
   sed -i 's|VEGDYN|'VEGDYN_e0000'|g' $HISTRC
#   sed -i 's|TP1|'TSOIL1TILE'|g' $HISTRC
#   sed -i 's|TP2|'TSOIL2TILE'|g' $HISTRC
#   sed -i 's|TP3|'TSOIL3TILE'|g' $HISTRC
#   sed -i 's|TP4|'TSOIL4TILE'|g' $HISTRC
#   sed -i 's|TP5|'TSOIL5TILE'|g' $HISTRC
#   sed -i 's|TP6|'TSOIL6TILE'|g' $HISTRC
#   sed -i 's|DATAATM|'DATAATM0000'|g' $HISTRC
endif

# fill in source 'GridComp' information for output variables

sed -i 's|GridComp|'$GridComp'|g' $HISTRC
