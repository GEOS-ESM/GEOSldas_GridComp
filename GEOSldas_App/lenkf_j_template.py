job_directive = {"NCCS": '''#!/bin/csh -f

# GEOSldas job script ("lenkf" = Land Ensemble Kalman Filter)
#
# usage: lenkf.j [-debug]

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################

#SBATCH --output={MY_EXPDIR}/scratch/GEOSldas_log_txt
#SBATCH --error={MY_EXPDIR}/scratch/GEOSldas_err_txt
#SBATCH --account={MY_ACCOUNT}
#SBATCH --time={MY_WALLTIME}
#SBATCH --ntasks={MY_NTASKS_MODEL}
#SBATCH --job-name={MY_JOB}
#SBATCH --qos={MY_QOS}
#SBATCH --constraint={MY_CONSTRAINT}
'''
,
"NAS": '''#!/bin/csh -f

# GEOSldas job script ("lenkf" = Land Ensemble Kalman Filter)
#
# usage: lenkf.j [-debug]

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#PBS -l walltime={MY_WALLTIME}
#PBS -l select={MY_NODES}:ncpus=40:mpiprocs=40:model={MY_CONSTRAINT}
#PBS -N {MY_JOB}
#PBS -q {MY_QOS}
#PBS -W group_list={MY_ACCOUNT}
#PBS -o {MY_EXPDIR}/scratch/GEOSldas_log_txt
#PBS -e {MY_EXPDIR}/scratch/GEOSldas_err_txt
#PBS -j oe
'''
}

job_body = '''
#######################################################################
#    System Settings and Architecture Specific Environment Variables
#######################################################################
umask 022
limit stacksize unlimited
setenv ARCH `uname`

setenv EXPID      {MY_EXPID}
setenv EXPDOMAIN  {MY_EXPDOMAIN}
setenv EXPDIR     {MY_EXPDIR}
setenv ESMADIR    $EXPDIR/build/
setenv GEOSBIN    $ESMADIR/bin/
# need to unsetenv LD_LIBRARY_PATH for execution of LDAS within the coupled land-atm DAS
unsetenv LD_LIBRARY_PATH

set debug_flag = 0
if ( "$1" == "-debug" ) then
  set debug_flag = 1
endif
unset argv
setenv argv

source $GEOSBIN/g5_modules

setenv BASEBIN ${{BASEDIR}}/Linux/bin

setenv MPI_STACK {DETECTED_MPI_STACK}

if ( ${{MPI_STACK}} == "openmpi" ) then

    # OPENMPI flags
    # Turn off warning about TMPDIR on NFS
    setenv OMPI_MCA_shmem_mmap_enable_nfs_warning 0
    # pre-connect MPI procs on mpi_init
    setenv OMPI_MCA_mpi_preconnect_all 1
    setenv OMPI_MCA_coll_tuned_bcast_algorithm 7
    setenv OMPI_MCA_coll_tuned_scatter_algorithm 2
    setenv OMPI_MCA_coll_tuned_reduce_scatter_algorithm 3
    setenv OMPI_MCA_coll_tuned_allreduce_algorithm 3
    setenv OMPI_MCA_coll_tuned_allgather_algorithm 4
    setenv OMPI_MCA_coll_tuned_allgatherv_algorithm 3
    setenv OMPI_MCA_coll_tuned_gather_algorithm 1
    setenv OMPI_MCA_coll_tuned_barrier_algorithm 0
    # required for a tuned flag to be effective
    setenv OMPI_MCA_coll_tuned_use_dynamic_rules 1
    # disable file locks
    setenv OMPI_MCA_sharedfp "^lockedfile,individual"

else if ( ${{MPI_STACK}} == "intelmpi" ) then

    setenv I_MPI_FABRICS shm:ofi
    setenv I_MPI_OFI_PROVIDER psm3

endif # MPI_STACK

# By default, ensure 0-diff across processor architecture by limiting MKL's freedom to pick algorithms.
# As of June 2021, MKL_CBWR=AVX2 is fastest setting that works for both haswell and skylake at NCCS.
# Change to MKL_CBWR=AUTO for fastest execution at the expense of results becoming processor-dependent.
#setenv MKL_CBWR "COMPATIBLE"
#setenv MKL_CBWR "AUTO"
setenv MKL_CBWR "AVX2"

#setenv LD_LIBRARY_PATH ${{LD_LIBRARY_PATH}}:${{BASEDIR}}/${{ARCH}}/lib
# reversed sequence for LADAS_COUPLING (Sep 2020)  (needed when coupling with ADAS using different BASEDIR)
setenv LD_LIBRARY_PATH ${{BASEDIR}}/${{ARCH}}/lib:${{ESMADIR}}/lib:${{LD_LIBRARY_PATH}}

setenv RUN_CMD "$GEOSBIN/esma_mpirun -np "

#######################################################################
#             Experiment Specific Environment Variables
#######################################################################

setenv    HOMDIR         $EXPDIR/run/
setenv    SCRDIR         $EXPDIR/scratch
setenv    LANDMODEL      {MY_LANDMODEL}
setenv    MYNAME         `finger $USER | cut -d: -f3 | head -1`
setenv    POSTPROC_HIST  {MY_POSTPROC_HIST}

# LADAS_COUPLING : 0 -- stand-alone LDAS (no coupling to ADAS)
#                : 1 -- LDAS coupled to central (deterministic) component of ADAS
#                : 2 -- LDAS coupled to atmospheric ensemble component of ADAS

setenv    LADAS_COUPLING   {MY_LADAS_COUPLING}
setenv    ENSEMBLE_FORCING {MY_ENSEMBLE_FORCING}
setenv    ADAS_EXPDIR      {MY_ADAS_EXPDIR}

set NENS = `grep NUM_LDAS_ENSEMBLE:  $HOMDIR/LDAS.rc | cut -d':' -f2`
set END_DATE  = `grep     END_DATE:  $HOMDIR/CAP.rc | cut -d':' -f2`
set NUM_SGMT  = `grep     NUM_SGMT:  $HOMDIR/CAP.rc | cut -d':' -f2`

#######################################################################
#  if LADAS_COUPLING==2, compute ens avg of atmens forcing
#######################################################################

if ( $LADAS_COUPLING == 2 && $ENSEMBLE_FORCING == "NO" ) then
   cd $HOMDIR
   set force_in  = $ADAS_EXPDIR
   set force_out = `grep MET_PATH: $HOMDIR/LDAS.rc | cut -d ':' -f2`
   python  $GEOSBIN/average_ensemble_forcing.py $force_in $force_out $NENS
endif

/bin/rm -f $HOMDIR/lenkf_job_completed.txt

#######################################################################
#                   Set Experiment Run Parameters
#######################################################################

#######################################################################
#        Move to Scratch Directory and Copy .rc .nml .rst files
#######################################################################

if (! -e $SCRDIR            ) mkdir -p $SCRDIR
cd $SCRDIR
/bin/rm -rf *.*
/bin/cp     $HOMDIR/cap_restart .
/bin/cp -f  $HOMDIR/*.rc .
/bin/cp -f  $HOMDIR/*.nml .

set LSMCHOICE = `grep -n -m 1 "LSM_CHOICE" $HOMDIR/LDAS.rc | cut -d':' -f3`

#######################################################################
# if $LADAS_COUPLING == 1:  LDAS coupled to central ADAS simulation
#######################################################################

if ( $LADAS_COUPLING == 1 ) then

   if ( $ENSEMBLE_FORCING == "YES" ) then

      # create perturbed forcing from central simulation and atm ensemble

      # python should come with ESMA_env g5_modules
      #module load python/GEOSpyD/Ana2019.03_py3.7
      set forcgrid = `grep GEOSldas.GRIDNAME LDAS.rc | cut -d':' -f2 | awk '{{print $1}}'`
      setenv GRID $forcgrid
      $GEOSBIN/enpert_forc.csh
      cd $SCRDIR

   endif
endif

#######################################################################
#              Create HISTORY Collection Directories
#######################################################################

set collections = ''
foreach line ("`cat HISTORY.rc`")
   set firstword  = `echo $line | awk '{{print $1}}'`
   set firstchar  = `echo $firstword | cut -c1`
   set secondword = `echo $line | awk '{{print $2}}'`

   if ( $firstword == "::" ) goto done

   if ( $firstchar != "#" ) then
      set collection  = `echo $firstword | sed -e "s/'//g"`
      set collections = `echo $collections $collection`
      if ( $secondword == :: ) goto done
   endif

   if ( $firstword == COLLECTIONS: ) then
      set collections = `echo $secondword | sed -e "s/'//g"`
   endif
end

done:

@ n_c = 0
if ($POSTPROC_HIST > 0) then
   foreach ThisCol ($collections)
      set ref_t = `cat HISTORY.rc | grep ${{ThisCol}}.ref_time: | cut -d':' -f2 | cut -d',' -f1`
      if ( $ref_t != '000000' ) then
         echo ${{ThisCol}}.ref_time should be '000000'
         @ n_c = $n_c + 1
      endif
   end
endif
if ($n_c >= 1) then
   exit
endif

#######################################################################
#                      Domain Decomposition
#######################################################################
set npes_nx = `grep NX: LDAS.rc | cut -d':' -f2 `
set npes_ny = `grep NY: LDAS.rc | cut -d':' -f2 `
@ numprocs = $npes_nx * $npes_ny
if( -e IMS.rc ) then
   set oldtasks = `head -n 1 IMS.rc`
   if($numprocs != $oldtasks) then
      $GEOSBIN/preprocess_ldas.x optimize ../input/tile.data $numprocs nothing nothing nothing
   endif
endif

if( -e JMS.rc ) then
   set oldtasks = `head -n 1 JMS.rc | cut -c1-5`
   if($numprocs != $oldtasks) then
      $GEOSBIN/preprocess_ldas.x optimize ../input/tile.data $numprocs nothing nothing nothing
   endif
endif

set gridname = `grep GEOSldas.GRIDNAME LDAS.rc | cut -d':' -f2 | cut -d'-'  -f2 | awk '{{print $1}}'`
if ( "$gridname" == "CF" ) then
   set new_ny = `echo "NY:  "$numprocs`
   sed -i "/NY:/c\\\\$new_ny" LDAS.rc
else
   set new_nx = `echo "NX:  "$numprocs`
   sed -i "/NX:/c\\\\$new_nx" LDAS.rc
endif

#######################################################################
#         Create Strip Utility to Remove Multiple Blank Spaces
#######################################################################

set      FILE = strip
/bin/rm $FILE
cat << EOF > $FILE
#!/bin/ksh
/bin/mv \$1 \$1.tmp
touch   \$1
while read line
do
echo \$line >> \$1
done < \$1.tmp
exit
EOF
chmod +x $FILE

##################################################################
######
######         Perform multiple iterations of Model Run
######
##################################################################

@ counter    = 1
while ( $counter <= ${{NUM_SGMT}} )

   /bin/rm -f  EGRESS.ldas
   /bin/cp -f $HOMDIR/CAP.rc .
   ./strip            CAP.rc

   # Set Time Variables for Current_(c), Ending_(e), and Segment_(s) dates
   # ---------------------------------------------------------------------
   set nymdc = `cat cap_restart | cut -c1-8`
   set nhmsc = `cat cap_restart | cut -c10-15`
   set nymde = `cat CAP.rc | grep END_DATE:     | cut -d: -f2 | cut -c2-9`
   set nhmse = `cat CAP.rc | grep END_DATE:     | cut -d: -f2 | cut -c11-16`
   set nymds = `cat CAP.rc | grep JOB_SGMT:     | cut -d: -f2 | cut -c2-9`
   set nhmss = `cat CAP.rc | grep JOB_SGMT:     | cut -d: -f2 | cut -c11-16`

   # Compute Time Variables at the Finish_(f) of current segment
   # -----------------------------------------------------------
   set nyear   = `echo $nymds | cut -c1-4`
   set nmonth  = `echo $nymds | cut -c5-6`
   set nday    = `echo $nymds | cut -c7-8`
   set nhour   = `echo $nhmss | cut -c1-2`
   set nminute = `echo $nhmss | cut -c3-4`
   set nsec    = `echo $nhmss | cut -c5-6`
          @ dt = $nsec + 60 * $nminute + 3600 * $nhour + 86400 * $nday

   set nymdf = $nymdc
   set nhmsf = $nhmsc
   set date  = `$GEOSBIN/tick $nymdf $nhmsf $dt`
   set nymdf =  $date[1]
   set nhmsf =  $date[2]
   set year  = `echo $nymdf | cut -c1-4`
   set month = `echo $nymdf | cut -c5-6`
   set day   = `echo $nymdf | cut -c7-8`

        @  month = $month + $nmonth
   while( $month > 12 )
        @  month = $month - 12
        @  year  = $year  + 1
   end
        @  year  = $year  + $nyear
        @ nymdf  = $year * 10000 + $month * 100 + $day

   if( $nymdf >  $nymde )    set nymdf = $nymde
   if( $nymdf == $nymde )    then
       if( $nhmsf > $nhmse ) set nhmsf = $nhmse
   endif

   set yearc = `echo $nymdc | cut -c1-4`
   set yearf = `echo $nymdf | cut -c1-4`

   # Prescribed LAI/SAI for CATCHCN
   # -------------------------------

   set PRESCRIBE_DVG = `grep PRESCRIBE_DVG LDAS.rc | cut -d':' -f2`
   if( ${{PRESCRIBE_DVG}} == 3 ) then
       set FCSTDATE = `grep FCAST_BEGTIME  $HOMDIR/LDAS.rc | cut -d':' -f2`
       if( `echo $FCSTDATE | cut -d' ' -f1` == "" ) then
           set CAPRES = `cat cap_restart`
           set CAPRES1 = `echo $CAPRES | cut -d' ' -f1`
           set CAPRES2 = `echo $CAPRES | cut -d' ' -f2`
           set CAPRES = 'FCAST_BEGTIME: '`echo $CAPRES1``echo $CAPRES2`
           echo $CAPRES >> $HOMDIR/LDAS.rc
           /bin/cp -p $HOMDIR/LDAS.rc .
       endif
   endif

   if( ${{PRESCRIBE_DVG}} >= 1 ) then

       # Modify local CAP.rc Ending date if Finish time exceeds Current year boundary
       # ----------------------------------------------------------------------------

       if( $yearf > $yearc ) then
          @ yearf = $yearc + 1
          @ nymdf = $yearf * 10000 + 0101
           set oldstring = `cat CAP.rc | grep END_DATE:`
           set newstring = "END_DATE: $nymdf $nhmsf"
           /bin/mv CAP.rc CAP.tmp
           cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
       endif

       # Creaate VEGDATA FIle Links
       # --------------------------

       if( ${{PRESCRIBE_DVG}} == 1 ) set VEGYR = $yearc
       if( ${{PRESCRIBE_DVG}} >= 2 ) set VEGYR = CLIM

       set FILE = vegfile
       set   nz = 1
       /bin/rm CNLAI*
       /bin/rm CNSAI*

       while ( $nz <= 3 )
   	set   nv = 1
   	while ($nv <= 4 )
   	    /bin/ln -s ../VEGDATA/CNLAI${{nv}}${{nz}}_${{VEGYR}}.data CNLAI${{nv}}${{nz}}.data
   	    /bin/ln -s ../VEGDATA/CNSAI${{nv}}${{nz}}_${{VEGYR}}.data CNSAI${{nv}}${{nz}}.data
   	    echo "CNLAI${{nv}}${{nz}}_FILE:                       CNLAI${{nv}}${{nz}}.data" >> $FILE
   	    echo "CNSAI${{nv}}${{nz}}_FILE:                       CNSAI${{nv}}${{nz}}.data" >> $FILE
   	    @ nv++
           end
   	@ nz++
       end
       /bin/mv LDAS.rc LDAS.rc.tmp
       cat LDAS.rc.tmp $FILE >> LDAS.rc
       /bin/rm LDAS.rc.tmp $FILE
   endif

   # ----------------------------------------------------------------------------

   set bYEAR = `cat cap_restart | cut -c1-4`
   set bMON  = `cat cap_restart | cut -c5-6`
   set bDAY  = `cat cap_restart | cut -c7-8`
   set bHour = `cat cap_restart | cut -c10-11`
   set bMin  = `cat cap_restart | cut -c12-13`

   if($counter == 1) then
      set logYEAR = $bYEAR
      set logMON  = $bMON
      set logDAY  = $bDAY
      set logHour = $bHour
      set logMin  = $bMin
   endif

   set old_mwrtm_file  =  $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{bYEAR}}/M${{bMON}}/${{EXPID}}.ldas_mwRTMparam.${{bYEAR}}${{bMON}}${{bDAY}}_${{bHour}}${{bMin}}z.nc4
   set old_catch_param =  $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{bYEAR}}/M${{bMON}}/${{EXPID}}.ldas_catparam.${{bYEAR}}${{bMON}}${{bDAY}}_${{bHour}}${{bMin}}z.bin
   if ( -l "$old_mwrtm_file" ) then
      set old_mwrtm_file = `/usr/bin/readlink -f $old_mwrtm_file`
   endif
   if ( -l "$old_catch_param" ) then
      set old_catch_param = `/usr/bin/readlink -f $old_catch_param`
   endif


   /bin/cp LDAS.rc  $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{bYEAR}}/M${{bMON}}/${{EXPID}}.ldas_LDAS_rc.${{bYEAR}}${{bMON}}${{bDAY}}_${{bHour}}${{bMin}}z.txt
   /bin/cp CAP.rc  $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{bYEAR}}/M${{bMON}}/${{EXPID}}.ldas_CAP_rc.${{bYEAR}}${{bMON}}${{bDAY}}_${{bHour}}${{bMin}}z.txt

   # Run GEOSldas.x
   # --------------
   # clean up
   $GEOSBIN/RmShmKeys_sshmpi.csh

   # Debugging
   # ---------
   if ( $debug_flag == 1 ) then
      echo ""
      echo "------------------------------------------------------------------"
      echo ""
      echo "lenkf.j -debug:"
      echo ""
      echo "To start debugging, you must now go to the experiment's scratch directory."
      echo "From there, source g5_modules and launch your debugging tool with GEOSldas.x, e.g.,"
      echo ""
      echo "   cd        $SCRDIR
      echo "   source    $GEOSBIN/g5_modules                        [for bash or zsh: source g5_modules.[z]sh]"
      echo "   module    load tview                                 [at NCCS]
      echo "   totalview $GEOSBIN/GEOSldas.x"
      echo ""
      echo "Availability of tools depends on the computing system and may require"
      echo "loading modules.  For more information, check with your computing center."
      echo "See also GEOSldas Wiki at https://github.com/GEOS-ESM/GEOSldas/wiki"
      echo ""
      exit
   endif

   @ oserver_nodes = {MY_OSERVER_NODES}
   @ writers = {MY_WRITERS_NPES}

   if (! $?SLURM_NTASKS) then
     set total_npes = `wc -l $PBS_NODEFILE | awk '{{print $1}}'`
   else
     set total_npes = $SLURM_NTASKS
   endif

   if ($oserver_nodes == 0) then
      set oserver_options = ""
   else
      set oserver_options = "--oserver_type multigroup --nodes_output_server $oserver_nodes  --npes_backend_pernode $writers"
   endif

   $RUN_CMD $total_npes $GEOSBIN/GEOSldas.x --npes_model $numprocs $oserver_options

   if( -e EGRESS.ldas ) then
      set rc = 0
      echo GEOSldas Run Status: $rc
   else
      set rc = -1
      echo GEOSldas Run Status: $rc
      echo "ERROR: GEOSldas run FAILED, exit without post-processing"
      exit $rc
   endif


   #######################################################################
   #              Move Legacy LDASsa Files to ana/ens_avg Directory
   #######################################################################

   # must be done before moving HISTORY files

   set ObsFcses = `ls *.ldas_ObsFcstAna.*.bin`
   foreach obsfcs ( $ObsFcses )
      set ThisTime = `echo $obsfcs | rev | cut -d'.' -f2 | rev`
      set TY = `echo $ThisTime | cut -c1-4`
      set TM = `echo $ThisTime | cut -c5-6`
      set THISDIR = $EXPDIR/output/$EXPDOMAIN/ana/ens_avg/Y${{TY}}/M${{TM}}/
      if (! -e $THISDIR            ) mkdir -p $THISDIR
      /bin/mv $obsfcs ${{THISDIR}}$obsfcs
   end

   set smapL4s = `ls *.ldas_tile_inst_smapL4SMaup.*.bin`
   foreach smapl4 ( $smapL4s )
      set ThisTime = `echo $smapl4 | rev | cut -d'.' -f2 | rev`
      set TY = `echo $ThisTime | cut -c1-4`
      set TM = `echo $ThisTime | cut -c5-6`
      set THISDIR = $EXPDIR/output/$EXPDOMAIN/ana/ens_avg/Y${{TY}}/M${{TM}}/
      if (! -e $THISDIR            ) mkdir -p $THISDIR
      /bin/mv $smapl4 ${{THISDIR}}$smapl4
   end


   #######################################################################
   #              Move HISTORY Files to cat/ens Directory
   #######################################################################

   set outfiles = `ls $EXPID.*[bin,nc4]`
   set TILECOORD=`ls ../output/*/rc_out/*ldas_tilecoord.bin`

   # Move current files to /cat/ens
   # ------------------------------

   foreach ofile ( $outfiles )
      set ThisTime = `echo $ofile | rev | cut -d'.' -f2 | rev`
      set TY = `echo $ThisTime | cut -c1-4`
      set TM = `echo $ThisTime | cut -c5-6`
      if ($NENS == 1) then
         set THISDIR = $EXPDIR/output/$EXPDOMAIN/cat/ens0000/Y${{TY}}/M${{TM}}/
      else
         set THISDIR = $EXPDIR/output/$EXPDOMAIN/cat/ens_avg/Y${{TY}}/M${{TM}}/
      endif
      if (! -e $THISDIR            ) mkdir -p $THISDIR

      set file_ext = `echo $ofile | rev | cut -d'.' -f1 | rev`

      if($file_ext == nc4) then
         /bin/mv $ofile $THISDIR/.
      else
         set binfile   = `echo $ofile | rev | cut -d'.' -f2- | rev`
         set decr_file = `echo $ofile | rev | cut -d'.' -f3- | rev`.ctl
         ($GEOSBIN/tile_bin2nc4.x $binfile $decr_file $TILECOORD ; \\
         /bin/mv ${{binfile}}.nc4 $THISDIR/. ; \\
         /bin/rm ${{binfile}}.bin) &
      endif
   end
   wait

   #######################################################################
   #              Post-Process model diagnostic output
   #              (1) Concatenate sub-daily files to daily files
   #              (2) Write monthly means
   #######################################################################

   if ($POSTPROC_HIST > 0) then

     set PWD = `pwd`

     if ($NENS == 1) then
        set OUTDIR = $EXPDIR/output/$EXPDOMAIN/cat/ens0000/
     else
        set OUTDIR = $EXPDIR/output/$EXPDOMAIN/cat/ens_avg/
     endif

     set MONTHDIRS = `ls -d $OUTDIR/*/*`

     foreach THISMONTH ($MONTHDIRS)

       set MM = `echo $THISMONTH | rev | cut -d'/' -f1 | cut -c1-2 | rev`
       set YYYY = `echo $THISMONTH | rev | cut -d'/' -f2 | cut -c1-4 | rev`
       set NDAYS = `cal $MM $YYYY | awk 'NF {{DAYS = $NF}}; END {{print DAYS}}'`

       cd $THISMONTH

       foreach ThisCol ($collections)
          # if monthly exists, move on to the next collection
          if (-f $EXPID.${{ThisCol}}.monthly.$YYYY$MM.nc4) continue

          # create daily and remove the sub-daily
          # ------------------------------------------------------------------
          set day=1
          while ($day <= $NDAYS)
             if ( $day < 10  ) set DD=0${{day}}
             if ( $day >= 10 ) set DD=${{day}}
             @ day++
             set time_steps = `ls -1 $EXPID.$ThisCol.${{YYYY}}${{MM}}${{DD}}_* | rev | cut -d'.' -f2 | rev`
             set LEN_SUB = `echo $#time_steps`

             # no file or just one file? nothing to concatenate, move on to the next collection
             if ($LEN_SUB <= 1) continue

             # check if day is complete (get HISTORY time step from first two files)
             set hour1   = `echo $time_steps[1] | cut -c10-11`
             set min1    = `echo $time_steps[1] | cut -c12-13`
             set hour2   = `echo $time_steps[2] | cut -c10-11`
             set min2    = `echo $time_steps[2] | cut -c12-13`
             @ dt_hist   = ($hour2 - $hour1) * 60 + ($min2 - $min1)
             @ N_per_day = (24 * 60) / $dt_hist
             # not enough sub-daily files? move on to the next collection
             if($LEN_SUB < $N_per_day) continue

             set tstep2 = \\"`echo $time_steps | sed 's/\ /\\","/g'`\\"

# ----------------------------------------------------------------------------
#
# WARNING: The following block MUST begin in column 1!!!  Do NOT indent!!!

cat << EOF > timestamp.cdl
netcdf timestamp {{
dimensions:
time = UNLIMITED ; // (NT currently)
string_length = 14 ;
variables:
char time_stamp (time, string_length) ;

data:

time_stamp =
DATAVALUES;
}}
EOF

             sed -i -e "s/NT/$LEN_SUB/g" timestamp.cdl
             sed -i -e "s/DATAVALUES/$tstep2/g" timestamp.cdl
             $BASEBIN/ncgen -k4 -o timestamp.nc4 timestamp.cdl
             $BASEBIN/ncrcat -h $EXPID.$ThisCol.${{YYYY}}${{MM}}${{DD}}_* ${{EXPID}}.${{ThisCol}}.$YYYY$MM$DD.nc4
             $BASEBIN/ncks -4 -h -v time_stamp timestamp.nc4 -A ${{EXPID}}.${{ThisCol}}.$YYYY$MM$DD.nc4
             /bin/rm timestamp.cdl
             /bin/rm timestamp.nc4
             # rudimentary check for desired nc4 file;  if ok, delete sub-daily files
             if ( -f ${{EXPID}}.${{ThisCol}}.$YYYY$MM$DD.nc4 ) then
                if ( ! -z ${{EXPID}}.${{ThisCol}}.$YYYY$MM$DD.nc4 ) then
                   /bin/rm $EXPID.${{ThisCol}}.${{YYYY}}${{MM}}${{DD}}_*.nc4
                endif
             endif
          end # concatenate for each day

          # write monthly mean file and (optionally) remove daily files
          # ------------------------------------------------------------------

          # NOTE: Collections written with daily frequency ("tavg24" and "inst24") have not
          #       been concatenated into daily files.  There are two possibilities  for the
          #       time stamps of files to be averaged:
          #         *.YYYYMMDD.*       daily files from concatenation of sub-daily files
          #         *.YYYYMMDD_HHMM.*  daily (avg or inst) files written directly by HISTORY.rc

          set time_steps  = `ls -1 $EXPID.$ThisCol.${{YYYY}}${{MM}}??.* | rev | cut -d'.' -f2 | rev`
          set time_steps_ = `ls -1 $EXPID.$ThisCol.${{YYYY}}${{MM}}??_* | rev | cut -d'.' -f2 | cut -d'_' -f2 | rev`
          set LEN  = `echo $#time_steps`
          set LEN_ = `echo $#time_steps_`

          # check if month is complete
          if ($LEN != 0) then
            set dayl = `echo $time_steps[$LEN] | cut -c1-8`
            set day1 = `echo $time_steps[1] | cut -c1-8`
            @ NAVAIL = ($dayl - $day1) + 1
          else if( $LEN_ != 0 ) then
            set dayl = `echo $time_steps_[$LEN_] | cut -c1-8`
            set day1 = `echo $time_steps_[1] | cut -c1-8`
            @ NAVAIL = ($dayl - $day1) + 1
          else
            @ NAVAIL = 0
          endif

          # not enough days for monthly mean? move on to the next collection
          if($NAVAIL != $NDAYS) continue

          # create monthly-mean nc4 file
          $BASEBIN/ncra -h $EXPID.$ThisCol.${{YYYY}}${{MM}}*.nc4 ${{EXPID}}.${{ThisCol}}.monthly.$YYYY$MM.nc4

          if($POSTPROC_HIST == 2) then
             # rudimentary check for desired nc4 file;  if ok, delete daily files
             if ( -f ${{EXPID}}.${{ThisCol}}.monthly.$YYYY$MM.nc4 ) then
                if ( ! -z ${{EXPID}}.${{ThisCol}}.monthly.$YYYY$MM.nc4 ) then
                   /bin/rm $EXPID.${{ThisCol}}.${{YYYY}}${{MM}}*
                endif
             endif
             continue
          endif

       end # each collection
       cd $PWD
     end # each month
   endif # POSTPROC_HIST > 0

   #######################################################################
   #   Rename Final Checkpoints => Restarts for Next Segment and Archive
   #        Note: cap_restart contains the current NYMD and NHMS
   #######################################################################

   set eYEAR = `cat cap_restart | cut -c1-4`
   set eMON  = `cat cap_restart | cut -c5-6`
   set eDAY  = `cat cap_restart | cut -c7-8`
   set eHour = `cat cap_restart | cut -c10-11`
   set eMin  = `cat cap_restart | cut -c12-13`

   # Create rc_out/YYYY/MM
   # ---------------------

   set THISDIR = $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{eYEAR}}/M${{eMON}}/
   if (! -e $THISDIR  ) mkdir -p $THISDIR

   # Move mwrtm and cat_param

   set new_mwrtm_file  =  $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{eYEAR}}/M${{eMON}}/${{EXPID}}.ldas_mwRTMparam.${{eYEAR}}${{eMON}}${{eDAY}}_${{eHour}}${{eMin}}z.nc4
   set new_catch_param =  $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{eYEAR}}/M${{eMON}}/${{EXPID}}.ldas_catparam.${{eYEAR}}${{eMON}}${{eDAY}}_${{eHour}}${{eMin}}z.bin

   if (-f $old_mwrtm_file) then
     if ( -l "$new_mwrtm_file" ) then
        /bin/rm -f $new_mwrtm_file
     endif
     /bin/ln -rs $old_mwrtm_file $new_mwrtm_file
     /bin/rm ../input/restart/mwrtm_param_rst
     /bin/ln -rs $new_mwrtm_file ../input/restart/mwrtm_param_rst
   endif

   if (-f $old_catch_param) then
     if ( -l "$new_catch_param" ) then
        /bin/rm -f $new_catch_param
     endif
     /bin/ln -rs $old_catch_param $new_catch_param
   endif

   # Move Intermediate Checkpoints to RESTARTS directory
   # ---------------------------------------------------

   @ inens = {MY_FIRST_ENS_ID}
   @ enens = $inens + $NENS
   while ($inens < $enens)
       if ($inens <10) then
          set ENSDIR = `echo ens000${{inens}}`
       else if($inens<100) then
          set ENSDIR=`echo ens00${{inens}}`
       else if($inens < 1000) then
          set ENSDIR =`echo ens0${{inens}}`
       else
          set ENSDIR = `echo ens${{inens}}`
       endif
       set ENSID = `echo $ENSDIR | cut -c4-7`
       set ENSID = _e${{ENSID}}
       if ( $NENS == 1) set ENSID =''
       set THISDIR = $EXPDIR/output/$EXPDOMAIN/rs/$ENSDIR/Y${{eYEAR}}/M${{eMON}}/
       if (! -e $THISDIR            ) mkdir -p $THISDIR
   
       set rstfs = (${{LANDMODEL}} 'landice')
       foreach rstf ( $rstfs )
          if (-f ${{rstf}}${{ENSID}}_internal_checkpoint ) then
             set tmp_file = $EXPDIR/output/$EXPDOMAIN/rs/$ENSDIR/Y${{eYEAR}}/M${{eMON}}/${{EXPID}}.${{rstf}}_internal_rst.${{eYEAR}}${{eMON}}${{eDAY}}_${{eHour}}${{eMin}}
             /bin/mv ${{rstf}}${{ENSID}}_internal_checkpoint $tmp_file
             /bin/rm -f $EXPDIR/input/restart/${{rstf}}${{ENSID}}_internal_rst
             /bin/ln -rs  $tmp_file $EXPDIR/input/restart/${{rstf}}${{ENSID}}_internal_rst
          endif
       end

       set rstf = 'landassim_obspertrseed'
       if (-f ${{rstf}}${{ENSID}}_checkpoint ) then
         set tmp_file = $EXPDIR/output/$EXPDOMAIN/rs/$ENSDIR/Y${{eYEAR}}/M${{eMON}}/${{EXPID}}.${{rstf}}_rst.${{eYEAR}}${{eMON}}${{eDAY}}_${{eHour}}${{eMin}}
         /bin/mv ${{rstf}}${{ENSID}}_checkpoint $tmp_file
         /bin/rm -f $EXPDIR/input/restart/${{rstf}}${{ENSID}}_rst
         /bin/ln -rs  $tmp_file $EXPDIR/input/restart/${{rstf}}${{ENSID}}_rst
       endif

       set rstf = 'landpert'
       if (-f ${{rstf}}${{ENSID}}_internal_checkpoint ) then
          set tmp_file = $EXPDIR/output/$EXPDOMAIN/rs/$ENSDIR/Y${{eYEAR}}/M${{eMON}}/${{EXPID}}.${{rstf}}_internal_rst.${{eYEAR}}${{eMON}}${{eDAY}}_${{eHour}}${{eMin}}
	  # copy generic restart file to final location/name but remove lat/lon variables
	  #  (lat/lon variables are not correct when running in EASE-grid tile space)
          $BASEBIN/ncks -4 -O -C -x -v lat,lon ${{rstf}}${{ENSID}}_internal_checkpoint $tmp_file
          /bin/rm -f ${{rstf}}${{ENSID}}_internal_checkpoint
          set old_rst = `/usr/bin/readlink -f $EXPDIR/input/restart/${{rstf}}${{ENSID}}_internal_rst`
          /bin/rm -f $EXPDIR/input/restart/${{rstf}}${{ENSID}}_internal_rst
          /bin/ln -rs $tmp_file $EXPDIR/input/restart/${{rstf}}${{ENSID}}_internal_rst
          /usr/bin/gzip $old_rst &
       endif
   
   # move intermediate check point files to  output/$EXPDOMAIN/rs/$ENSDIR/Yyyyy/Mmm/ directories
   # -------------------------------------------------------------------------------------------
   
       set rstfiles1 = `ls ${{LANDMODEL}}${{ENSID}}_internal_checkpoint.*`
       set rstfiles2 = `ls landpert${{ENSID}}_internal_checkpoint.*`
       set rstfiles3 = `ls landassim_obspertrseed${{ENSID}}_checkpoint.*`
       set rstfiles4 = `ls landice${{ENSID}}_internal_checkpoint.*`
   
       foreach rfile ( $rstfiles1 $rstfiles4 ) 
          set ThisTime = `echo $rfile | rev | cut -d'.' -f2 | rev`
          set TY = `echo $ThisTime | cut -c1-4`
          set TM = `echo $ThisTime | cut -c5-6`
          set THISDIR = $EXPDIR/output/$EXPDOMAIN/rs/$ENSDIR/Y${{TY}}/M${{TM}}/
          if (! -e $THISDIR            ) mkdir -p $THISDIR
          /bin/mv $rfile ${{THISDIR}}${{EXPID}}.${{LANDMODEL}}_internal_rst.${{ThisTime}}.nc4
          /usr/bin/gzip ${{THISDIR}}${{EXPID}}.${{LANDMODEL}}_internal_rst.${{ThisTime}}.nc4 &
       end

       foreach rfile ( $rstfiles2 )
          set ThisTime = `echo $rfile | rev | cut -d'.' -f2 | rev`
          set TY = `echo $ThisTime | cut -c1-4`
          set TM = `echo $ThisTime | cut -c5-6`
          set THISDIR = $EXPDIR/output/$EXPDOMAIN/rs/$ENSDIR/Y${{TY}}/M${{TM}}/
          if (! -e $THISDIR            ) mkdir -p $THISDIR
             ($BASEBIN/ncks -4 -O -C -x -v lat,lon $rfile ${{THISDIR}}${{EXPID}}.landpert_internal_rst.${{ThisTime}}.nc4;\\
               /usr/bin/gzip ${{THISDIR}}${{EXPID}}.landpert_internal_rst.${{ThisTime}}.nc4; \\
               /bin/rm -f $rfile) &
       end

       foreach rfile ( $rstfiles3 )
          set ThisTime = `echo $rfile | rev | cut -d'.' -f2 | rev`
          set TY = `echo $ThisTime | cut -c1-4`
          set TM = `echo $ThisTime | cut -c5-6`
          set THISDIR = $EXPDIR/output/$EXPDOMAIN/rs/$ENSDIR/Y${{TY}}/M${{TM}}/
          if (! -e $THISDIR            ) mkdir -p $THISDIR
             /bin/mv $rfile ${{THISDIR}}${{EXPID}}.landassim_obspertrseed_rst.${{ThisTime}}.nc4
       end

       @ inens ++
   end  ## end of while ($inens < $NENS)
   wait
   #####################
   # update cap_restart
   # ##################

   set CO2_BEFORE = `grep CO2_YEAR: LDAS.rc | cut -d':' -f2`

   if ( $CO2_BEFORE >= 1 ) then

       # Update reference year for Carbon Tracker CO2
       ##############################################

       set CAP_BEFORE = `head -1 $HOMDIR/cap_restart | cut -c1-4`
       @ DY = $CAP_BEFORE - $CO2_BEFORE
       @ CO2_AFTER = `head -1 cap_restart | cut -c1-4` - $DY
       set CO2UPDATE = `echo "CO2_YEAR:  "$CO2_AFTER`
       sed -i "/CO2_YEAR:/c\\$CO2UPDATE" LDAS.rc
       /bin/rm -f $HOMDIR//LDAS.rc
       /bin/cp -p LDAS.rc $HOMDIR/LDAS.rc
   endif

   /bin/rm -f $HOMDIR/cap_restart
   /bin/cp cap_restart $HOMDIR/cap_restart

   #######################################################################
   #              Update Iteration Counter
   #######################################################################

   set enddate = `echo  $END_DATE | cut -c1-8`
   set endhour = `echo  $END_DATE | cut -c10-11`
   set capdate = `cat cap_restart | cut -c1-8`
   set caphour = `cat cap_restart | cut -c10-11`

   if ( $capdate < $enddate ) then
     @ counter = $counter + 1
   else if ( $capdate == $enddate && $caphour < $endhour ) then
     @ counter = $counter + 1
   else
     @ counter = ${{NUM_SGMT}} + 1
   endif

## End of the while ( $counter <= ${{NUM_SGMT}} ) loop ##
end

#######################################################################
#                 Set Next Log and Error Files
#######################################################################

set logfile = $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{logYEAR}}/M${{logMON}}/${{EXPID}}.ldas_log.${{logYEAR}}${{logMON}}${{logDAY}}_${{logHour}}${{logMin}}z.txt
set errfile = $EXPDIR/output/$EXPDOMAIN/rc_out/Y${{logYEAR}}/M${{logMON}}/${{EXPID}}.ldas_err.${{logYEAR}}${{logMON}}${{logDAY}}_${{logHour}}${{logMin}}z.txt

if (-f GEOSldas_log_txt) then
   /bin/cp GEOSldas_log_txt $logfile
   /bin/rm -f GEOSldas_log_txt
endif

if(-f GEOSldas_err_txt) then
  /bin/cp GEOSldas_err_txt $errfile
  /bin/rm -f GEOSldas_err_txt
endif

#######################################################################
#                 Re-Submit Job
#######################################################################

if ( $LADAS_COUPLING > 0 ) then
   if ( $rc == 0 ) then
      ##update CAP.rc END_DATE in $HOMDIR/
      set date  = `$GEOSBIN/tick $nymdf $nhmsf $dt`
      set nymdend =  $date[1]
      set nhmsend =  $date[2]
      cd  $HOMDIR 
      set oldstring = `cat CAP.rc | grep END_DATE:`
      set newstring = "END_DATE: $nymdend $nhmsend"
      /bin/mv CAP.rc CAP.tmp
      cat CAP.tmp | sed -e "s?$oldstring?$newstring?g" > CAP.rc
      /bin/rm -f CAP.tmp
      echo 'SUCCEEDED' > $HOMDIR/lenkf_job_completed.txt
   endif
else
   if ( $rc == 0 ) then
   cd   $HOMDIR
   #don't change below line(not even extra space)
   if($capdate<$enddate) {SBATCHQSUB} $HOMDIR/lenkf.j
  endif
endif
'''
