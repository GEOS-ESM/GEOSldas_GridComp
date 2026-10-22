#!/usr/bin/env python3

"""
ldas_config.py: reading and verifying ldas_setup inputs (mixin class for ldas.ldas)

LdasConfig provides the steps called by ldas.__init__: command line args, exeinp
file (or coupled land-atm DAS settings), model options, BCS files, restart checks,
resource manager inputs, and experiment paths.  Nothing is created on disk here.
"""

import os
import glob
from datetime                     import datetime
from datetime                     import timedelta
from collections                  import OrderedDict
from dateutil.relativedelta       import relativedelta
from setup_utils                  import (parseInputFile, printDictionary, hours_to_hhmmss,
                                          getExeKeys, getResourceKeys, read_gridname,
                                          LdasSetupError, require, glob_first,
                                          bcs_version_from_path, ym_dir, ymd_hm)

# -----------------------------------------------------------------------------------
# Hard-wired (NCCS/discover) paths used for the coupled land-atm DAS (ladas_cpl > 0)
# and in user messages.
# -----------------------------------------------------------------------------------
LADAS_BCS_PATH   = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles"
LADAS_MWRTM_PATH = "/discover/nobackup/projects/gmao/smap/LDAS_inputs_for_LADAS/RTM_params/RTMParam_SMAP_L4SM_v006/"
ROUTE_RST_DIR    = "/discover/nobackup/projects/gmao/bcs_shared/restarts/surface/route/"

# LSM_CHOICE -> name of land model (Catchment) gridded component / restart prefix
LSM_CHOICE_TO_CATCH = {1: 'catch', 2: 'catchcnclm40', 3: 'catchcnclm45', 4: 'catchcnclm51'}


class LdasConfig:
    """ mixin: input processing for class ldas (see ldas.py) """

    # -----------------------------------------------------------------------------------
    def _init_attributes(self, cmdLineArgs):
        """ copy command line args and initialize all instance variables """
        self.exeinpfile         = cmdLineArgs['exeinpfile']
        self.batinpfile         = cmdLineArgs['batinpfile']
        exphome_                = cmdLineArgs['exphome'].rstrip('/')
        require(os.path.isdir(exphome_), 'exphome [%s] does not exist' % exphome_)
        self.exphome            = os.path.abspath(exphome_)
        self.verbose            = cmdLineArgs['verbose']

        # command line args for coupled land-atm DAS (see "help" strings in parseCmdLine() for details)
        self.ladas_cpl          = cmdLineArgs['ladas_cpl']
        self.nymdb              = cmdLineArgs['nymdb']
        self.nhmsb              = cmdLineArgs['nhmsb']
        self.agcm_res           = cmdLineArgs['agcm_res']
        self.bcs_version        = cmdLineArgs['bcs_version']
        self.rstloc             = cmdLineArgs['rstloc']
        self.varwindow          = cmdLineArgs['varwindow']
        self.nens               = cmdLineArgs['nens']

        # obsolete command line args
        self.runmodel           = cmdLineArgs['runmodel']
        if self.runmodel :
          print('\n The option "--runmodel" is out of date, not necessary anymore. \n')

        self.daysperjob         = cmdLineArgs['daysperjob']
        self.monthsperjob       = cmdLineArgs['monthsperjob']

        self.ExeInputs          = OrderedDict()
        self.RmInputs           = OrderedDict()
        self.rundir             = None
        self.blddir             = None
        self.blddirLn           = None
        self.outdir             = None
        self.out_path           = None
        self.inpdir             = None
        self.exefyl             = None
        self.isZoomIn           = False
        self.isZoomInLand       = False
        self.isZoomInLake       = False
        self.isZoomInLandice    = False
        self.catch              = ''
        self.has_mwrtm          = False
        self.has_vegopacity     = False
        self.assim              = False
        self.has_landassim_seed = False
        self.has_geos_pert      = False
        self.nSegments          = 1
        self.perturb            = 0
        self.first_ens_id       = 0
        self.in_rstdir          = None
        self.in_tilefile        = None                 # default string
        self.ens_id_width       = 6                    # _eXXXX
        self.bcs_dir_land       = ''
        self.bcs_dir_geom       = ''
        self.bcs_dir_landshared = ''
        self.bcs_dir_landiceshared = ''
        self.tile_types         = ''
        self.with_land          = False
        self.with_lake          = False
        self.with_landice       = False
        self.run_route          = 0
        self.with_issm          = False
        self.adas_expdir        = ''
        self.do_issm            = 0
        self.restart            = None                 # RESTART option as string: '1', '2', '3', or 'M'

    # -----------------------------------------------------------------------------------
    def _check_ladas_cmdline(self, cmdLineArgs):
        """
        Set self.ladas_cpl.  If the exeinp file does not exist, this must be a coupled
        land-atm DAS setup and all coupled-mode command line args are required.
        If the exeinp file exists, ladas_cpl=0 (--ladas_cpl is ignored; intended).
        """
        if not os.path.exists(cmdLineArgs['exeinpfile']):
           # make sure all necessary command line arguments were supplied
           require(self.ladas_cpl is not None,
                   "exeinpfile [%s] does not exist (allowed only for coupled land-atm DAS with --ladas_cpl)" % cmdLineArgs['exeinpfile'])
           self.ladas_cpl = int(self.ladas_cpl)
           require(self.ladas_cpl > 0, "If not ladas coupling, exeinpfile must be provided.")
           for name_ in ['nymdb', 'nhmsb', 'agcm_res', 'bcs_version', 'rstloc', 'varwindow', 'nens']:
              require(getattr(self, name_) is not None,
                      "Must have command line arg --%s for coupled land-atm DAS." % name_)
        else:
           self.ladas_cpl = 0

    # -----------------------------------------------------------------------------------
    def _set_restart_option(self):
        """ RESTART option (default 1); for RESTART=M no restart exp is needed """
        if 'RESTART' not in self.ExeInputs :
            self.ExeInputs['RESTART'] = "1"

        RESTART_str  = str(self.ExeInputs['RESTART'])
        require(RESTART_str in ["1", "2", "3", "M"], "RESTART must be 1, 2, 3, or M")
        self.restart = RESTART_str

        if RESTART_str == 'M':
            # use archived MERRA-2 restarts, no need for RESTART_ID/DOMAIN/PATH
            self.ExeInputs['RESTART_ID']     = 'None'
            self.ExeInputs['RESTART_DOMAIN'] = 'None'
            self.ExeInputs['RESTART_PATH']   = 'None'

    # -----------------------------------------------------------------------------------
    def _apply_ladas_coupling(self):
        """
        Coupled land-atm DAS (ladas_cpl > 0): set/overwrite ExeInputs from the command
        line args and hard-wired values.
        """
        self.ExeInputs['BEG_DATE']                 = f"{self.nymdb} {self.nhmsb}"
        rstloc_                                    = self.rstloc.rstrip('/')                  # remove trailing '/'
        require(os.path.isdir(rstloc_), 'rstloc [%s] is not a directory' % rstloc_)
        self.rstloc                                = os.path.abspath(rstloc_)
        self.ExeInputs['RESTART_PATH']             = os.path.dirname( self.rstloc)
        self.ExeInputs['RESTART_ID']               = os.path.basename(self.rstloc)
        self.adas_expdir                           = os.path.dirname( self.exphome)
        self.ExeInputs['ADAS_EXPDIR']              = self.adas_expdir
        self.adas_expid                            = os.path.basename(self.adas_expdir)
        self.ExeInputs['MET_TAG']                  = self.adas_expid + '__bkg'

        if   self.ladas_cpl == 1 :
            # ldas coupled with determistic component of ADAS
            self.ExeInputs['EXP_ID']               = self.adas_expid  + '_LDAS'
            self.ExeInputs['MET_PATH']             = self.adas_expdir + '/recycle/holdpredout'
            self.ExeInputs['ENSEMBLE_FORCING']     = 'NO'
        elif self.ladas_cpl == 2 :
            # ldas coupled with ensemble component of ADAS
            self.ExeInputs['EXP_ID']               = self.adas_expid  + '_LDAS4ens'
            self.ExeInputs['MET_PATH']             = self.adas_expdir + '/atmens/mem'
            self.ExeInputs['ENSEMBLE_FORCING']     = 'YES'
        else :
            raise LdasSetupError("Unknown value of ladas_cpl: %s" % self.ladas_cpl)

        self.ExeInputs['NUM_LDAS_ENSEMBLE']        = self.nens                                # fvsetup finds Nens by counting restart files
        self.first_ens_id                          = 1                                        # match ADAS convention
        self.ExeInputs['FIRST_ENS_ID']             = self.first_ens_id

        self.agcm_res                              = 'CF' + self.agcm_res                     # change format to "CFnnnn"
        self.ExeInputs['EXP_DOMAIN']               = self.agcm_res +'x6C_GLOBAL'

        # when coupled to ADAS, "BCS_PATH" EXCLUDE bcs version info
        # hard-wired BCS_PATH for now
        self.ExeInputs['BCS_PATH']                 = LADAS_BCS_PATH
        self.ExeInputs['BCS_PATH']                 = self.ExeInputs['BCS_PATH'].rstrip('/') + '/' + self.bcs_version
        if self.bcs_version == "Icarus-NLv3" :
            self.ExeInputs['BCS_PATH']             = self.ExeInputs['BCS_PATH']  + '_new_layout'
        self.ExeInputs['BCS_RESOLUTION']           = self.agcm_res +'x6C_'       + self.agcm_res +'x6C'
        self.ExeInputs['RESTART_DOMAIN']           = self.agcm_res +'x6C_GLOBAL'

        # the following are not in default ExeInputs list; hardwire for now
        self.ExeInputs['MWRTM_PATH']               = LADAS_MWRTM_PATH
        self.ExeInputs['LAND_ASSIM']               = "YES"
        self.ExeInputs['MET_HINTERP']              = 0
        self.landassim_dt                          = 10800                                    # seconds
        # make sure ADAS analysis window [minutes] is multiple of LANDASSIM_DT [seconds]
        if int(self.varwindow) % (self.landassim_dt/60) == 0 :
            self.ExeInputs['LANDASSIM_DT']         = self.landassim_dt
        else :
            raise LdasSetupError("LANDASSIM_DT is inconsistent with ADAS analysis window.")
        self.ExeInputs['LANDASSIM_T0']             = "013000"                                 # HHMMSS
        jsgmt1                                     = "00000000"
        jsgmt2                                     = hours_to_hhmmss(int(self.varwindow)/60)  # convert minutes to HHMMSS
        self.ExeInputs['JOB_SGMT']                 = f"{jsgmt1} {jsgmt2}"
        self.ExeInputs['NUM_SGMT']                 = 1
        self.ExeInputs['FORCE_DTSTEP']             = 3600

        # determine END_DATE = BEG_DATE + TIME_STEP_OF_ADAS_CYCLE
        _beg_date                                  = datetime.strptime( self.ExeInputs['BEG_DATE'], "%Y%m%d %H%M%S")
        _end_date                                  = _beg_date + timedelta(hours=int(self.varwindow)/60)
        self.ExeInputs['END_DATE']                 = _end_date.strftime("%Y%m%d %H%M%S")

    # -----------------------------------------------------------------------------------
    def _set_model_options(self):
        """ land model, tile types, ISSM, routing, ensemble members/ids """
        RESTART_str = self.restart

        if 'LSM_CHOICE' not in self.ExeInputs:
            self.ExeInputs['LSM_CHOICE'] = 1
        _lsm_choice_int = int(self.ExeInputs['LSM_CHOICE'])
        require(_lsm_choice_int in LSM_CHOICE_TO_CATCH,
                "Unknown LSM_CHOICE=%d; must be 1 (Catchment) or 2 (Catchment-CN4.0)" % _lsm_choice_int)
        self.catch = LSM_CHOICE_TO_CATCH[_lsm_choice_int]

        self.tile_types = self.ExeInputs.get('TILE_TYPES',"100").split()
        if "100" in self.tile_types :
            self.with_land = True
            require(int(self.ExeInputs['LSM_CHOICE']) <= 2, "LSM_CHOICE=3 (Catchment-CN4.5) is no longer supported. Please set LSM_CHOICE to 1 (Catchment) or 2 (Catchment-CN4.0)")

        if "19" in self.tile_types :
            self.with_lake = True

        if "20" in self.tile_types :
            self.with_landice = True

        # DO_ISSM is optional in exeinp (default 0, as in GEOS_SurfaceGridComp.rc);
        # self.do_issm is used only by ldas_setup and not written to LDAS.rc
        self.do_issm = int(self.ExeInputs.get('DO_ISSM', 0))
        if self.with_landice and self.do_issm == 1:
            self.with_issm = True

        self.run_route = int(self.ExeInputs.get('RUN_ROUTE',0))

        # fail early (before any remapping) for unsupported option combinations
        require(not (self.run_route > 0 and RESTART_str == 'M'),
                "RUN_ROUTE does not support MERRA-2 option (RESTART=M)")
        require(not (self.with_issm and RESTART_str == 'M'),
                "RESTART=M does not support DO_ISSM=1. Please provide a compatible ISSM restart and use RESTART=2.")

        self.nens = int(self.ExeInputs['NUM_LDAS_ENSEMBLE'])             # fails if value of Nens is not an integer
        self.first_ens_id = int(self.ExeInputs.get('FIRST_ENS_ID',0))
        self.perturb = int(self.ExeInputs.get('PERTURBATIONS',0))
        if self.nens > 1:
           self.perturb = 1
        self.ensdirs  = ['ens%04d'%iens for iens in range(self.first_ens_id, self.nens  + self.first_ens_id)]
        # if self.ens_id_width = 4, _width = '_e%04d'
        _width = '_e%0{}d'.format(self.ens_id_width-2)
        # self.ensids will be a list of [_e0000, _e0001, ...]
        self.ensids   = [ _width%iens   for iens in range(self.first_ens_id, self.nens  + self.first_ens_id)]
        if (self.nens == 1) :
            self.ensdirs_avg = self.ensdirs
            self.ensids=['']
        else :
            self.ensdirs_avg = self.ensdirs + ['ens_avg']

    # -----------------------------------------------------------------------------------
    def _resolve_bcs_files(self):
        """
        BCS directories, catchment.def, tile file, land BC files, GRIDNAME, EASE Pfafstetter tile file.
        For RESTART=1, tile and BC files are taken from the rc_out dir of the restart exp.
        """
        RESTART_str = self.restart

        # assemble bcs sub-directories
        self.bcs_dir_land       = self.ExeInputs['BCS_PATH']+ '/land/'       + self.ExeInputs['BCS_RESOLUTION']+'/'
        self.bcs_dir_geom       = self.ExeInputs['BCS_PATH']+ '/geometry/'   + self.ExeInputs['BCS_RESOLUTION']+'/'
        self.bcs_dir_landshared = self.ExeInputs['BCS_PATH']+ '/land/shared/'
        self.bcs_dir_landiceshared = self.ExeInputs['BCS_PATH']+ '/landice/shared/'


        # make sure MET_PATH and RESTART_PATH have trailing '/'
        if self.ExeInputs['MET_PATH'][-1] != '/':
           self.ExeInputs['MET_PATH'] = self.ExeInputs['MET_PATH']+'/'
        if self.ExeInputs['RESTART_PATH'][-1] != '/':
           self.ExeInputs['RESTART_PATH'] = self.ExeInputs['RESTART_PATH']+'/'

        # make sure catchment and vegdyn restart files exist (at least one of each)
        if 'CATCH_DEF_FILE' not in self.ExeInputs :
           self.ExeInputs['CATCH_DEF_FILE']= self.bcs_dir_land + 'clsm/catchment.def'
        if self.with_land :
           require(os.path.isfile(self.ExeInputs['CATCH_DEF_FILE']), "CATCH_DEF_FILE [%s] does not exist" % self.ExeInputs['CATCH_DEF_FILE'])

        # assigning BC files
        self.ExeInputs['LNFM_FILE'] = ''
        self.ExeInputs['EASE_PFAF_TILE_FILE'] = ''
        tile_file_format = self.ExeInputs.get('TILE_FILE_FORMAT', 'DEFAULT')
        domain_      = ''  # tile geometry domain
        land_domain_ = ''  # Land BC domain
        inpdir_  = self.bcs_dir_land
        inpgeom_ = self.bcs_dir_geom


        # find restart files and tile file in restart dir (if necessary)
        inp_ = self.ExeInputs['RESTART_PATH']+'/'.join([self.ExeInputs['RESTART_ID'],'output',
                    self.ExeInputs['RESTART_DOMAIN'], 'rc_out/'])
        if RESTART_str == '1' :
           inpdir_  = inp_
           inpgeom_ = inp_

           # verify that BCS_VERSION and BCS_RESOLUTION of restart experiment match those of new experiment
           BCS_txt  = glob.glob(inp_ + 'BCS_info.txt')
           if len(BCS_txt)== 0:
              print("Warning: BCS_info.txt not found for restart experiment. User is responsible for ensuring consistency of restart and experiment BCS.")
           elif len(BCS_txt)== 1:
              BCS_tmp   = parseInputFile(BCS_txt[0])
              bcs_version_ = bcs_version_from_path(self.ExeInputs['BCS_PATH'])
              require(bcs_version_ == BCS_tmp['BCS_VERSION'], "BCS version (=trailing dir of BCS_PATH) does not match version from restart dir ("  + bcs_version_ + ")")
              require(self.ExeInputs['BCS_RESOLUTION'] == BCS_tmp['BCS_RESOLUTION'], "BCS_RESOLUTION does not match resolution from restart dir (" + BCS_tmp['BCS_RESOLUTION'] + ")")

           # Any reduced tile type produces a .domain tile file
           txt_tile = glob.glob(inp_ + '*.domain')
           if len(txt_tile) > 0:
              domain_  = '.domain'
           # Land BCs have .domain only when Land itself was reduced
           land_bc = glob.glob(inp_ + 'green_clim_*.data.domain')
           
           if len(land_bc) > 0:
               land_domain_ = '.domain'              

        elif RESTART_str == '2':
           txt_tile = glob.glob(inp_ + '*.domain')
           require(len(txt_tile) == 0, "For RESTART=2, must restart from and run on global domain")
           in_tilefiles_ = glob.glob(inp_+'MAPL_*.til')
           if len(in_tilefiles_) == 0 :
              in_tilefiles_ = glob.glob(inp_+'/*.nc4')
           if len(in_tilefiles_) == 0 :
              in_tilefiles_ = glob.glob(inp_+'/*.til')
           require(len(in_tilefiles_) > 0, "For RESTART=2, cannot find tile file (MAPL_*.til, *.nc4, *.til) in [%s]" % inp_)
           # for EASE tile space, pick "standard" tile file (one tile per grid cell) and not "-Pfafstetter" tile file
           if 'EASEv' in in_tilefiles_[0]:
              in_tilefiles_ = [ item_ for item_ in in_tilefiles_ if '-Pfafstetter' not in item_ ]
           self.in_tilefile =os.path.realpath(in_tilefiles_[0])

        inpdir_  = os.path.realpath(inpdir_)+'/'
        inpgeom_ = os.path.realpath(inpgeom_)+'/'

        txt_tile = glob.glob(inpgeom_ + '*.til' + domain_)
        for f in txt_tile:
           if 'MAPL_' in os.path.basename(f):
              txt_tile = [f]
              break
        nc4_tile = glob.glob(inpgeom_ + '*.nc4' + domain_)
        if len(nc4_tile) > 0:
           if 'EASEv' in nc4_tile[0]:
             nc4_tile = [ item_ for item_ in nc4_tile if '-Pfafstetter' not in item_ ]
        require(len(txt_tile + nc4_tile) > 0, "Cannot find tile file (*.til, *.nc4) in [%s]" % inpgeom_)
        if tile_file_format.upper() == 'TXT'     :
           require(len(txt_tile) > 0, "TILE_FILE_FORMAT=TXT but no *.til file in [%s]" % inpgeom_)
           self.ExeInputs['TILING_FILE'] =  txt_tile[0]
        if tile_file_format.upper() == 'DEFAULT' : self.ExeInputs['TILING_FILE'] = (txt_tile+nc4_tile)[-1]

        self.ExeInputs['GRN_FILE']    = glob_first(inpdir_ + 'green_clim_*.data'+land_domain_ , 'GRN_FILE')
        self.ExeInputs['LAI_FILE']    = glob_first(inpdir_ + 'lai_clim_*.data'  +land_domain_ , 'LAI_FILE')
        tmp_ = glob.glob(inpdir_ + 'lnfm_clim_*.data'+land_domain_)
        if (len(tmp_) ==1) :
           self.ExeInputs['LNFM_FILE'] = tmp_[0]

        self.ExeInputs['NDVI_FILE']   = glob_first(inpdir_ + 'ndvi_clim_*.data'+land_domain_ , 'NDVI_FILE')
        self.ExeInputs['NIRDF_FILE']  = glob_first(inpdir_ + 'nirdf_*.dat'     +land_domain_ , 'NIRDF_FILE')
        self.ExeInputs['VISDF_FILE']  = glob_first(inpdir_ + 'visdf_*.dat'     +land_domain_ , 'VISDF_FILE')

        # assigning Gridname (read_gridname handles *.til, *.nc4, and trailing '.domain')
        if 'GRIDNAME' not in self.ExeInputs :
           self.ExeInputs['GRIDNAME'] = read_gridname(self.ExeInputs['TILING_FILE'])

        # to run routing on standard EASE tile space, need EASE_PFAF_TILE_FILE
        if (self.run_route > 0 and 'EASE' in self.ExeInputs['GRIDNAME']):
           tmp_ =  glob.glob(inpgeom_ + '*Pfafstetter.nc4' + land_domain_ )
           if (len(tmp_) > 0) :
              self.ExeInputs['EASE_PFAF_TILE_FILE'] = tmp_[0]

    # -----------------------------------------------------------------------------------
    def _set_exe_defaults(self):
        if 'POSTPROC_HIST' not in self.ExeInputs:
            self.ExeInputs['POSTPROC_HIST'] = 0

        if 'RUN_IRRIG' not in self.ExeInputs:
            self.ExeInputs['RUN_IRRIG'] = 0

        if 'AEROSOL_DEPOSITION' not in self.ExeInputs:
            self.ExeInputs['AEROSOL_DEPOSITION'] = 0

    # -----------------------------------------------------------------------------------
    def _check_restart_files(self):
        """ make sure the restart files to be used exist (RESTART=1,2,3) """
        RESTART_str = self.restart

        if RESTART_str not in ['1','2','3']:
           return
        rstid = self.ExeInputs['RESTART_ID']
        stamp = ymd_hm(self.begDates[0])
        self.in_rstdir = self.ExeInputs['RESTART_PATH']+'/'.join([rstid,'output',
                                                                self.ExeInputs['RESTART_DOMAIN'],'rs',self.ensdirs[0],
                                                                ym_dir(self.begDates[0])])
        def _rst_file(comp):
           return self.in_rstdir+'/'+rstid+'.'+comp+'_internal_rst.'+stamp

        required = []
        if self.with_land     : required.append(self.catch)
        if self.with_lake     : required.append('lake')
        if self.with_landice  : required.append('landice')
        if self.run_route > 0 : required.append('route')
        if self.with_issm     : required.append('issm')
        for comp in required:
           require(os.path.isfile(_rst_file(comp)), comp+'_internal_rst file [%s] does not exist!' % _rst_file(comp))

        # use landpert restarts if available (RESTART=1 only)
        if self.with_land and RESTART_str == '1' and os.path.isfile(_rst_file('landpert')):
           self.has_geos_pert = True

    # -----------------------------------------------------------------------------------
    def _set_mwrtm(self):
        """ land assimilation flag; mwRTM parameter and vegopacity files """
        # DEAL WITH mwRTM input from exec
        self.assim = True if self.ExeInputs.get('LAND_ASSIM', 'NO').upper() == 'YES' and self.with_land else False
        # verify mwrtm file
        if 'MWRTM_PATH' in self.ExeInputs and self.with_land :
           self.ExeInputs['MWRTM_PATH'] = self.ExeInputs['MWRTM_PATH']+'/'+ self.ExeInputs['BCS_RESOLUTION']+'/'
           mwrtm_param_file_ = self.ExeInputs['MWRTM_PATH']+'mwRTM_param.nc4'
           vegopacity_file_  = self.ExeInputs['MWRTM_PATH']+'vegopacity.bin'
           if os.path.isfile(mwrtm_param_file_) :
              self.has_mwrtm = True
              self.mwrtm_file = mwrtm_param_file_
           else :
              raise LdasSetupError('MWRTM_PATH/BCS_RESOLUTION [%s] must contain mwRTM_param.nc4' % self.ExeInputs['MWRTM_PATH'])
           if os.path.isfile(vegopacity_file_) :
              self.has_vegopacity = True
              self.ExeInputs['VEGOPACITY_FILE'] = vegopacity_file_

    # -----------------------------------------------------------------------------------
    def _read_resource_inputs(self, cmdLineArgs):
        """ resource manager (SLURM/PBS) inputs: batinp file, or defaults for coupled land-atm DAS """
        if self.ladas_cpl == 0 :
           self.RmInputs  = parseInputFile(cmdLineArgs['batinpfile'])
        else :
           self.RmInputs['account']         = cmdLineArgs['account']
           self.RmInputs['walltime']        = "01:00:00"
           self.RmInputs['ntasks_model']    = 120

        self._verifyResourceInputs()

        # print rm inputs
        if self.verbose:
            print ('\n\nInputs for resource manager:')
            printDictionary(self.RmInputs)
            print ('\n\n')

    # -----------------------------------------------------------------------------------
    def _set_paths(self):
        """ installation (bindir, blddir, executable) and experiment directories """
        self.bindir     = os.path.dirname(os.path.realpath(__file__))
        self.blddir     = self.bindir.rsplit('/',1)[0]
        exefyl          = '/bin/GEOSldas.x'
        tmp_execfyl     = self.blddir   + exefyl
        require(os.path.isfile(tmp_execfyl), 'Executable [%s] does not exist!' % tmp_execfyl)
        self.expdir     = self.exphome  + '/' + self.ExeInputs['EXP_ID']
        self.rundir     = self.expdir   + '/run'
        self.inpdir     = self.expdir   + '/input'
        self.outdir     = self.expdir   + '/output'
        self.scratchdir = self.expdir   + '/scratch'
        self.blddirLn   = self.expdir   + '/build'
        self.out_path   = self.outdir   + '/'+self.ExeInputs['EXP_DOMAIN']
        self.rc_out     = self.outdir   + '/'+self.ExeInputs['EXP_DOMAIN']+'/rc_out/'
        self.rstdir     = self.outdir   + '/'+self.ExeInputs['EXP_DOMAIN']+'/rs/'
        self.exefyl     = self.blddirLn + exefyl

    # -----------------------------------------------------------------------------------
    def _set_oserver(self):
        """ output server defaults """
        # default is set to 0 ( no output server)
        if 'oserver_nodes' not in self.RmInputs :
           self.RmInputs['oserver_nodes'] = 0

        if (int(self.RmInputs['oserver_nodes']) >=1) :
           self.ExeInputs['WRITE_RESTART_BY_OSERVER'] = "YES"
           # set default for now
           if 'writers-per-node' not in self.RmInputs:
              self.RmInputs['writers-per-node'] = 5
        else:
           self.RmInputs['writers-per-node'] = 0

    # -----------------------------------------------------------------------------------
    def _verifyExeInputs(self):
       ExeInputs = self.ExeInputs
       #) verify keys
       option = '1'

       rqdExeInpKeys = getExeKeys(option)
       for key in rqdExeInpKeys:
          require(key in ExeInputs, '"%s" is required in the inputs (from exeinpfile or command line)' % (key))

       _mydir = self.exphome + '/' + self.ExeInputs['EXP_ID']
       require(not os.path.isdir(_mydir), 'Dir [%s] already exists!' % _mydir)
       _mydir = None

       # nens is an integer and =1 for model run
       require(self.nens>0, 'NUM_LDAS_ENSEMBLE [%d] <= 0' % self.nens)
    # -----------------------------------------------------------------------------------
    def _verifyResourceInputs(self):
       #-----
       # verify resource input keys are correct
       #-----
       ResourceInputs  = self.RmInputs
       rqdRmInpKeys    = getResourceKeys('required')
       optSlurmInpKeys = getResourceKeys('optional')
       allKeys = rqdRmInpKeys + optSlurmInpKeys
       for key in rqdRmInpKeys:
         require(key in ResourceInputs, '"%s" is required in the inputs (from batinpfile or command line)' % (key))

       for key in ResourceInputs:
         require(key in allKeys, '"%s" in batinp file is not recognized' % key)

    # -----------------------------------------------------------------------------------
    def _calculateJobSegments(self):
        """
        begDates/endDates: start/end of each job segment; nSegments; job_sgmt (CAP.rc line).
        NOTE: one job segment = NUM_SGMT * JOB_SGMT, where JOB_SGMT is converted with
              1 month = 30 days and 1 year = 365 days (intended; see REFACTOR_PLAN.md B7).
        """
        self.begDates = [datetime.strptime(self.ExeInputs['BEG_DATE'], '%Y%m%d %H%M%S')]
        self.endDates = [datetime.strptime(self.ExeInputs['END_DATE'], '%Y%m%d %H%M%S')]

        require(self.endDates[0]>self.begDates[0], 'END_DATE <= BEG_DATE')

        self.job_sgmt = []
        if 'JOB_SGMT' in self.ExeInputs:
            self.job_sgmt.append("JOB_SGMT: "+self.ExeInputs['JOB_SGMT'])
        else:
            _datediff = relativedelta(self.endDates[0],self.begDates[0])
            self.ExeInputs['JOB_SGMT'] = "%04d%02d%02d %02d%02d%02d" %(_datediff.years,
                                               _datediff.months,
                                               _datediff.days,
                                               _datediff.hours,
                                               _datediff.minutes,
                                               _datediff.seconds)
            self.job_sgmt.append("JOB_SGMT: "+self.ExeInputs['JOB_SGMT'])

        if 'NUM_SGMT' not in self.ExeInputs:
           self.ExeInputs['NUM_SGMT'] = 1

        _years  = int(self.ExeInputs['JOB_SGMT'][ 0: 4])
        _months = int(self.ExeInputs['JOB_SGMT'][ 4: 6])
        _days   = int(self.ExeInputs['JOB_SGMT'][ 6: 8])
        require(self.ExeInputs['JOB_SGMT'][8] == ' ' and self.ExeInputs['JOB_SGMT'][9] != ' ', "JOB_SGMT format is not right (expected 'yyyymmdd hhmmss')")
        _hours  = int(self.ExeInputs['JOB_SGMT'][ 9:11])
        _mins   = int(self.ExeInputs['JOB_SGMT'][11:13])
        _seconds= int(self.ExeInputs['JOB_SGMT'][13:15])


        _difftime =timedelta(days = _years*365+_months*30+_days,hours = _hours,minutes=_mins,seconds=_seconds)
        _difftime = int(self.ExeInputs['NUM_SGMT'])*_difftime
        _d = self.begDates[0]
        _endDate = self.endDates[0]
        _d = _d + _difftime
        while _d < _endDate :
            self.nSegments +=1
            self.begDates.append(_d)
            self.endDates.insert(-1,_d)
            _d = _d+ _difftime
