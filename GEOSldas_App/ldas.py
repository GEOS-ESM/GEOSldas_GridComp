#!/usr/bin/env python3

"""
ldas class is used by ldas_setup

The class is split into mixins (one module each):
  ldas_config.py  : LdasConfig  - read/verify inputs (steps called by __init__)
  ldas_staging.py : LdasStaging - createDirStructure, createLnRstBc
  ldas_rc.py      : LdasRc      - createRCFiles
  ldas.py         : ldas        - __init__ driver, createRunScripts, createBatchRun
This file is processed by cmake configure_file (@GEOS_SITE@, @MPI_STACK@); keep all
such placeholders in this file.
"""

import os
from lenkf_j_template             import job_directive, job_body
from setup_utils                  import (parseInputFile, printDictionary, require, ym_dir, ymd_hm)
from ldas_config                  import LdasConfig
from ldas_staging                 import LdasStaging
from ldas_rc                      import LdasRc

# Site-dependent settings for the job scripts (GEOS_SITE is set by cmake configure_file)
#   submit     : batch submit command
#   qos        : default quality-of-service
#   constraint : default node constraint
SITE_CONFIG = {
    'NCCS': {'submit': 'sbatch', 'qos': 'allnccs', 'constraint': '"[mil|cas]"'},
    'NAS' : {'submit': 'qsub',   'qos': 'normal',  'constraint': 'cas_ait'    },
}


class ldas(LdasConfig, LdasStaging, LdasRc):

    def __init__(self, cmdLineArgs):
        """
        Read and verify all inputs (command line, exeinp file, batinp file), derive
        settings, and locate input files.  Nothing is created on disk here.

        The order of the steps matters: ExeInputs is an OrderedDict and the order in
        which keys are added determines the order of entries in LDAS.rc.
        """
        # List of keywords, including:
        # - keywords that are excluded from LDAS.rc (i.e., only needed in pre- or post-processing)
        # - keywords whose values are re-assigned between the exeinp and LDAS.rc files
        # - keywords that are not in exeinp file but are in LDAS.rc file
        self.NoneLDASrcKeys=['EXP_ID', 'EXP_DOMAIN',
            'BEG_DATE', 'END_DATE','RESTART','RESTART_PATH',
            'RESTART_DOMAIN','RESTART_ID','BCS_PATH','TILING_FILE','GRN_FILE','LAI_FILE','LNFM_FILE','NIRDF_FILE',
            'VISDF_FILE','CATCH_DEF_FILE','NDVI_FILE', 'EASE_PFAF_TILE_FILE',
            'NML_INPUT_PATH','HISTRC_FILE','RST_FROM_GLOBAL','JOB_SGMT','NUM_SGMT','POSTPROC_HIST',
            'MINLON','MAXLON','MINLAT','MAXLAT','EXCLUDE_FILE','INCLUDE_FILE','MWRTM_PATH','GRIDNAME',
            'ADAS_EXPDIR', 'BCS_RESOLUTION', 'TILE_FILE_FORMAT' ]

        self.GEOS_SITE = "@GEOS_SITE@"
        require(self.GEOS_SITE in SITE_CONFIG,
                'GEOS_SITE [%s] is not supported; must be one of %s' % (self.GEOS_SITE, sorted(SITE_CONFIG)))
        self.site = SITE_CONFIG[self.GEOS_SITE]

        self._init_attributes(cmdLineArgs)
        self._check_ladas_cmdline(cmdLineArgs)

        # read exeinp file (empty dict for coupled land-atm DAS)
        self.ExeInputs = parseInputFile(cmdLineArgs['exeinpfile'], ladas_cpl = self.ladas_cpl )

        self._set_restart_option()
        if self.ladas_cpl > 0 :
            self._apply_ladas_coupling()

        if self.verbose:
            print ('\nInputs from exeinp file:\n')
            printDictionary(self.ExeInputs)

        self._set_model_options()
        self._verifyExeInputs()
        self._calculateJobSegments()
        self._resolve_bcs_files()
        self._set_exe_defaults()
        self._check_restart_files()
        self._set_mwrtm()
        self._read_resource_inputs(cmdLineArgs)
        self._set_paths()
        self._set_oserver()

    # -----------------------------------------------------------------------------------
    # job scripts
    # -----------------------------------------------------------------------------------

    def createBatchRun(self):
        """
        ldas_batchrun.j: submit all job segments at once, each depending on the previous one
        (the automatic re-submission at the end of lenkf.j is commented out while submitting).
        """
        submit = self.site['submit']
        lines  = ["#!/bin/bash -f", "",
                  "sed -i 's/if($capdate<$enddate) "+submit+"/#if($capdate<$enddate) "+submit+"/g' lenkf.j", ""]
        for iseg in range(self.nSegments):
            if self.GEOS_SITE == 'NAS':
                # PBS: qsub prints the job id only
                dep = '' if iseg == 0 else '-W depend=afterok:$jobid%d ' % (iseg-1)
                lines.append("jobid%d=$(qsub %slenkf.j)" % (iseg, dep))
            else:
                dep = '' if iseg == 0 else '--dependency=afterok:$jobid%d ' % (iseg-1)
                lines.append("jobid%d=$(echo $(sbatch %slenkf.j) | cut -d' ' -f 4)" % (iseg, dep))
            lines.append("echo $jobid%d" % iseg)
        lines += ["", "sed -i 's/#if($capdate<$enddate) "+submit+"/if($capdate<$enddate) "+submit+"/g' lenkf.j", "", ""]

        batchrun = self.rundir+'/ldas_batchrun.j'
        with open(batchrun, 'w') as fout:
            fout.write('\n'.join(lines))
        os.chmod(batchrun, 0o755)

    # -----------------------------------------------------------------------------------
    def createRunScripts(self):
        """ lenkf.j: job script (from lenkf_j_template.py) """
        my_qos     = self.RmInputs.get('qos',        self.site['qos'])
        my_job     = self.RmInputs.get('job_name',   self.ExeInputs['EXP_ID'])
        constraint = self.RmInputs.get('constraint', self.site['constraint'])

        start = self.begDates[0]
        expid = self.ExeInputs['EXP_ID']
        def _rc_out_txt(kind):
            # e.g., ../output/GLOBAL/rc_out/Y2020/M01/expid.ldas_log.20200101_0000z.txt (relative to rundir)
            return os.path.relpath('/'.join([self.outdir, self.ExeInputs['EXP_DOMAIN'], 'rc_out', ym_dir(start),
                                             '.'.join([expid, kind, ymd_hm(start)+'z', 'txt'])]), self.rundir)

        my_nodes=''
        if 'ntasks-per-node' in self.RmInputs:
           ntasks_per_node = int(self.RmInputs['ntasks-per-node'])
           ntasks          = int(self.RmInputs['ntasks_model'])
           require(ntasks%ntasks_per_node == 0, 'Please make ntasks_model a multiple of ntasks-per-node')
           nodes           = ntasks//ntasks_per_node
           my_nodes  = '#SBATCH --nodes=' + str(nodes) +' --ntasks-per-node=' + str(ntasks_per_node)
           if (ntasks_per_node > 46):
              require(constraint != 'cas', "Make sure constraint is compatible with ntasks-per-node")

        lenkf_str= (job_directive[self.GEOS_SITE]+job_body).format(
                  SBATCHQSUB         = self.site['submit'],
                  MY_ACCOUNT         = self.RmInputs['account'],
                  MY_WALLTIME        = self.RmInputs['walltime'],
                  MY_NTASKS_MODEL    = str(self.RmInputs['ntasks_model']),
                  MY_NODES           = my_nodes,
                  MY_CONSTRAINT      = constraint,
                  MY_OSERVER_NODES   = str(self.RmInputs['oserver_nodes']),
                  MY_WRITERS_NPES    = str(self.RmInputs['writers-per-node']),
                  MY_QOS             = my_qos,
                  MY_JOB             = my_job,
                  MY_EXPID           = expid,
                  MY_EXPDOMAIN       = self.ExeInputs['EXP_DOMAIN'],
                  MY_LOGFILE         = _rc_out_txt('ldas_log'),
                  MY_ERRFILE         = _rc_out_txt('ldas_err'),
                  MY_LANDMODEL       = self.catch,
                  MY_POSTPROC_HIST   = str(self.ExeInputs['POSTPROC_HIST']),
                  MY_FIRST_ENS_ID    = str(self.first_ens_id),
                  MY_LADAS_COUPLING  = str(self.ladas_cpl),
                  MY_ENSEMBLE_FORCING= self.ExeInputs.get('ENSEMBLE_FORCING', 'NO').upper(),
                  MY_ADAS_EXPDIR     = self.adas_expdir,
                  MY_EXPDIR          = self.expdir,
                  DETECTED_MPI_STACK = "@MPI_STACK@",
                  )

        lenkf = self.rundir+'/lenkf.j'
        with open(lenkf,'wt') as fout :
           fout.write(lenkf_str)
        os.chmod(lenkf, 0o755)

        print ('\nExperiment directory: %s' % self.expdir)
        print ()
