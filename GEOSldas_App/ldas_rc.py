#!/usr/bin/env python3

"""
ldas_rc.py: nml and rc files in the run directory (mixin class for ldas.ldas)
"""

import os
import glob
import shutil
import tempfile
from collections                  import OrderedDict
from remap_utils                  import STRETCH_GRID
from setup_utils                  import (parseInputFile, require, run_cmd, replace_in_file,
                                          copy_with_suffix)


class LdasRc:
    """ mixin: createRCFiles for class ldas (see ldas.py) """

    # -----------------------------------------------------------------------------------
    # createRCFiles: nml and rc files in run/
    # -----------------------------------------------------------------------------------

    def createRCFiles(self):
        """
        Write resource files to the run directory:
          - copies of exeinp/batinp files (for the record)
          - LDASsa_DEFAULT/SPECIAL nml files
          - IMS.rc/JMS.rc (optimized domain decomposition)
          - HISTORY.rc, ExtData.rc, CAP.rc, LDAS.rc, cap_restart
        """
        for mydir in [self.blddirLn, self.rundir]:
            require(os.path.isdir(mydir), 'dir [%s] does not exist!' % mydir)

        etcdir = self.blddirLn + '/etc'
        for f in ['CAP', 'ExtData', 'HIST', 'HISTdet', 'HISTens', 'LDAS']:
            require(os.path.isfile(etcdir+'/GEOSldas_'+f+'.rc'), 'Missing rc template [%s]' % (etcdir+'/GEOSldas_'+f+'.rc'))

        if self.ladas_cpl == 0:
            # copy ldas_setup exeinp and batinp input files to rundir (for the record)
            copy_with_suffix(self.exeinpfile, self.rundir)
            copy_with_suffix(self.batinpfile, self.rundir)

        self._copy_nml(etcdir)
        optinxny = self._optimize_layout()
        self._write_history_rc(etcdir)
        # just copy an empty ExtData.rc
        shutil.copy2(etcdir+'/GEOSldas_ExtData.rc', self.rundir+'/ExtData.rc')
        self._write_cap_rc(etcdir)
        self._write_ldas_rc(etcdir, optinxny)

        with open(self.rundir+'/'+'cap_restart','w') as fout:
            fout.write(self.begDates[0].strftime('%Y%m%d %H%M%S'))

    # -----------------------------------------------------------------------------------
    def _copy_nml(self, etcdir):
        """ DEFAULT nml files from etc/; SPECIAL nml files from NML_INPUT_PATH (or etc/ if coupled) """
        for nmlfile in glob.glob(etcdir+'/LDASsa_DEFAULT_inputs_*.nml'):
            shutil.copy2(nmlfile, self.rundir+'/'+os.path.basename(nmlfile))

        special_nml=[]
        if self.ladas_cpl > 0:
            special_nml= glob.glob(etcdir+'/LDASsa_SPECIAL_inputs_*.nml')
        elif 'NML_INPUT_PATH' in self.ExeInputs :
            special_nml = glob.glob(self.ExeInputs['NML_INPUT_PATH']+'/LDASsa_SPECIAL_inputs_*.nml')
        for nmlfile in special_nml:
            shutil.copy2(nmlfile, self.rundir+'/'+os.path.basename(nmlfile))

        if self.ladas_cpl > 0:
            # edit resolution info in ensupd nml file (exact replacement, like 'sed -i s/../../g')
            nml_ = self.rundir+'/LDASsa_SPECIAL_inputs_ensupd.nml'
            with open(nml_) as f:
                txt_ = f.read()
            with open(nml_, 'w') as f:
                f.write(txt_.replace('<CFnnnn>', self.agcm_res))

    # -----------------------------------------------------------------------------------
    def _optimize_layout(self):
        """
        Run preprocess_ldas.x optimize: writes IMS.rc/JMS.rc to rundir and returns the
        NX/NY layout (dict).  For cube-sphere tile space, ntasks_model may be adjusted.
        """
        tmp_ = tempfile.NamedTemporaryFile(delete=False)
        tmp_.close()
        try:
            print  ("Optimizing... decomposition of processes.... \n")
            print ("IMS.rc or JMS.rc would be generated on " + self.rundir)
            run_cmd(self.bindir + '/preprocess_ldas.x optimize '+ self.inpdir+'/tile.data '
                    + str(self.RmInputs['ntasks_model']) + ' ' + tmp_.name + ' ' + self.rundir + ' '
                    + '_'.join(self.tile_types))
            optinxny = parseInputFile(tmp_.name)
        finally:
            os.remove(tmp_.name)
        if (int(optinxny['NX']) == 1):
            if int(optinxny['NY']) != int(self.RmInputs['ntasks_model']):
                self.RmInputs['ntasks_model']=optinxny['NY']
                print ('adjust ntasks_model %d for cubed-sphere grid' % int(self.RmInputs['ntasks_model']))
        return optinxny

    # -----------------------------------------------------------------------------------
    def _write_history_rc(self, etcdir):
        """
        HISTORY.rc:
          - coupled land-atm DAS: GEOSldas_HISTdet.rc (ladas_cpl=1) or GEOSldas_HISTens.rc (ladas_cpl=2)
          - otherwise: user's HISTRC_FILE (if given), or GEOSldas_HIST.rc processed by process_hist.csh
        """
        tmprcfile = self.rundir+'/HISTORY.rc'

        if self.ladas_cpl > 0:
            shutil.copy2(etcdir + ('/GEOSldas_HISTdet.rc' if self.ladas_cpl == 1 else '/GEOSldas_HISTens.rc'), tmprcfile)
            replace_in_file(tmprcfile, 'GEOSldas_expid', self.ExeInputs['EXP_ID'])
            replace_in_file(tmprcfile, 'GRIDNAME',       self.ExeInputs['GRIDNAME'])
            return

        user_histrc = ''
        if 'HISTRC_FILE' in self.ExeInputs :
            _tmpfile = self.ExeInputs['HISTRC_FILE'].replace("'",'').replace('"','')
            if os.path.isfile(_tmpfile) :
                user_histrc = self.ExeInputs['HISTRC_FILE']
            else :
                require(not _tmpfile.strip(), "HISTRC_FILE: %s is NOT a file. " %_tmpfile)

        if user_histrc :
            shutil.copy2(user_histrc, tmprcfile)
        else :
            shutil.copy2(etcdir+'/GEOSldas_HIST.rc', tmprcfile)
            OUTxd = 'OUT1d' if 'EASE' in self.ExeInputs['GRIDNAME'] else 'OUT2d'
            run_cmd(' '.join([self.bindir +'/process_hist.csh', tmprcfile, OUTxd, self.ExeInputs['GRIDNAME'],
                              str(self.ExeInputs['LSM_CHOICE']), str(self.ExeInputs['AEROSOL_DEPOSITION']),
                              str(self.ExeInputs['RUN_IRRIG']), str(self.nens), str(self.do_issm)]))
        replace_in_file(tmprcfile, 'GEOSldas_expid', self.ExeInputs['EXP_ID'])

    # -----------------------------------------------------------------------------------
    def _write_cap_rc(self, etcdir):
        tmprcfile = self.rundir+'/CAP.rc'
        shutil.copy2(etcdir+'/GEOSldas_CAP.rc', tmprcfile)
        replace_in_file(tmprcfile, 'JOB_SGMT:', self.job_sgmt[0])
        replace_in_file(tmprcfile, 'NUM_SGMT:', 'NUM_SGMT: %d'% int(self.ExeInputs['NUM_SGMT']))
        replace_in_file(tmprcfile, 'BEG_DATE:', self.begDates[ 0].strftime('BEG_DATE: %Y%m%d %H%M%S'))
        replace_in_file(tmprcfile, 'END_DATE:', self.endDates[-1].strftime('END_DATE: %Y%m%d %H%M%S'))

    # -----------------------------------------------------------------------------------
    def _write_ldas_rc(self, etcdir, optinxny):
        """
        LDAS.rc = optimized layout (NX, NY)
                + defaults from GEOS_SurfaceGridComp.rc, overwritten by GEOSldas_LDAS.rc,
                  overwritten by ExeInputs (except NoneLDASrcKeys), plus derived entries
                + OUT_PATH, EXP_ID, TILING_FILE
        """
        surf_defaults = parseInputFile(etcdir+'/GEOS_SurfaceGridComp.rc', ladas_cpl=self.ladas_cpl)
        ldas_defaults = parseInputFile(etcdir+'/GEOSldas_LDAS.rc',        ladas_cpl=self.ladas_cpl)
        ldasrcInp = self._ldas_rc_dict(surf_defaults, ldas_defaults)

        with open(self.rundir+'/LDAS.rc','w') as fout:
            for key,val in list(optinxny.items()) + list(ldasrcInp.items()):
                fout.write((key+":").ljust(36)+str(val)+'\n')
            fout.write("OUT_PATH:".ljust(36)+self.out_path+'\n')
            fout.write("EXP_ID:".ljust(36)+self.ExeInputs['EXP_ID']+'\n')
            fout.write("TILING_FILE:".ljust(36)+"../input/tile.data\n")

    # -----------------------------------------------------------------------------------
    def _ldas_rc_dict(self, surf_defaults, ldas_defaults):
        """ assemble LDAS.rc entries (no file I/O); the order of entries is preserved """
        ldasrcInp = OrderedDict()
        ldasrcInp.update(surf_defaults)              # land defaults
        ldasrcInp.update(ldas_defaults)              # ldas defaults, may overwrite land defaults
        for key,val in self.ExeInputs.items():       # exeinp, may overwrite defaults
            if key not in self.NoneLDASrcKeys:
                ldasrcInp[key]= val

        # '%s' template for ensemble member id in file names
        tmpl_ = '%s' if self.nens > 1 else ''

        if self.perturb == 1:
           ldasrcInp['PERTURBATIONS'] ='1'

        # restart components: (RC key prefix, file name prefix)
        rst = []
        if self.with_land :
           # BC files
           for key in ['GREEN','LAI','LNFM','NDVI','NIRDF','VISDF']:
              ldasrcInp[key+'_FILE'] = '../input/'+key.lower()+'.data'
           if self.ExeInputs['EASE_PFAF_TILE_FILE'] !='':
              ldasrcInp['EASE_PFAF_TILE_FILE'] = '../input/EASE_pfaf_tile_file.nc4'

           if('catchcn' in self.catch):
              ldasrcInp['CO2_MonthlyMean_DiurnalCycle_FILE']= '../input/CO2_MonthlyMean_DiurnalCycle.nc4'
           else:
              # remove catchcn-specific entries that do not apply to catch model
              for key in ['DTCN', 'ATM_CO2', 'CO2', 'CO2_YEAR', 'PRESCRIBE_DVG']:
                 ldasrcInp.pop(key,None)

           catch_ = self.catch.upper()
           # avoid duplicates
           for key in [catch_+'_INTERNAL_RESTART_TYPE', catch_+'_INTERNAL_CHECKPOINT_TYPE', 'VEGDYN_INTERNAL_RESTART_TYPE']:
              ldasrcInp.pop(key, None)

           rst += [(catch_, self.catch), ('VEGDYN', 'vegdyn')]

        if self.with_lake:
           rst.append(('LAKE', 'lake'))
        if self.with_landice:
           rst.append(('LANDICE', 'landice'))
           if self.with_issm:
              rst.append(('ISSM', 'issm'))
        if self.run_route > 0:
           rst.append(('ROUTE', 'route'))

        if self.has_mwrtm :
           ldasrcInp['LANDASSIM_INTERNAL_RESTART_FILE'] = '../input/restart/mwrtm_param_rst'
        if self.has_vegopacity :
           ldasrcInp['VEGOPACITY_FILE'] = '../input/vegopacity.data'
        if self.nens > 1 :
           ldasrcInp['ENS_ID_WIDTH'] = str(self.ens_id_width)
        if self.has_landassim_seed and self.assim :
           ldasrcInp['LANDASSIM_OBSPERTRSEED_RESTART_FILE'] = '../input/restart/landassim_obspertrseed'+tmpl_+'_rst'
        if self.assim:
           ldasrcInp['LANDASSIM_OBSPERTRSEED_CHECKPOINT_FILE'] = 'landassim_obspertrseed'+tmpl_+'_checkpoint'

        for key, val in rst:
            ldasrcInp[key+'_INTERNAL_RESTART_FILE'] = '../input/restart/'+val+tmpl_+'_internal_rst'
            if key != 'VEGDYN':                       # no checkpoint for vegdyn
               ldasrcInp[key+'_INTERNAL_CHECKPOINT_FILE'] = val+tmpl_+'_internal_checkpoint'

        # LANDPERT restart file
        if (self.perturb == 1):
            ldasrcInp['LANDPERT_INTERNAL_RESTART_FILE'] = '../input/restart/landpert'+tmpl_+'_internal_rst'
            # for lat/lon and EASE tile space, specify LANDPERT checkpoint file here (via MAPL);
            # for cube-sphere tile space, Landpert GC will set up LANDPERT checkpoint file
            if ('-CF' not in self.ExeInputs['GRIDNAME']):
                ldasrcInp['LANDPERT_INTERNAL_CHECKPOINT_FILE'] = 'landpert'+tmpl_+'_internal_checkpoint'

        # add items for stretched grid
        if '-SG' in self.ExeInputs['BCS_RESOLUTION']:
           pos_ = self.ExeInputs['BCS_RESOLUTION'].find('-SG')
           SG   = self.ExeInputs['BCS_RESOLUTION'][pos_+1:pos_+6]    # get ID of stretched grid (e.g., SG002)
           ldasrcInp['STRETCH_FACTOR'] = STRETCH_GRID[SG][0]
           ldasrcInp['TARGET_LAT']     = STRETCH_GRID[SG][1]
           ldasrcInp['TARGET_LON']     = STRETCH_GRID[SG][2]

        return ldasrcInp
