#!/usr/bin/env python3

import os
import glob
import shutil
import subprocess as sp
import tempfile
from datetime                     import datetime
from datetime                     import timedelta
from collections                  import OrderedDict
from dateutil.relativedelta       import relativedelta
from remap_utils                  import STRETCH_GRID, init_merra2, get_config_from_answers
from remap_other_restarts         import other_restarts
from remap_catchANDcn             import catchANDcn
from lenkf_j_template             import job_directive, job_body
from setup_utils                  import (parseInputFile, printDictionary, hours_to_hhmmss,
                                          getExeKeys, getResourceKeys, read_gridname,
                                          LdasSetupError, require, run_cmd, glob_first,
                                          bcs_version_from_path, replace_in_file, ym_dir, ymd_hm,
                                          copy_with_suffix)

"""
ldas class is used by ldas_setup
"""

# -----------------------------------------------------------------------------------
# Hard-wired (NCCS/discover) paths used for the coupled land-atm DAS (ladas_cpl > 0)
# and in user messages.
# -----------------------------------------------------------------------------------
LADAS_BCS_PATH   = "/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles"
LADAS_MWRTM_PATH = "/discover/nobackup/projects/gmao/smap/LDAS_inputs_for_LADAS/RTM_params/RTMParam_SMAP_L4SM_v006/"
ROUTE_RST_DIR    = "/discover/nobackup/projects/gmao/bcs_shared/restarts/surface/route/"

# LSM_CHOICE -> name of land model (Catchment) gridded component / restart prefix
LSM_CHOICE_TO_CATCH = {1: 'catch', 2: 'catchcnclm40', 3: 'catchcnclm45', 4: 'catchcnclm51'}

class ldas:

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
    def _write_domain_def(self):
        """ write simulation domain namelist (&domain_inputs) to a temporary file for preprocess_ldas.x """
        # default is global
        _domain_dic=OrderedDict()
        _domain_dic['MINLON']=-180.
        _domain_dic['MAXLON']= 180.
        _domain_dic['MINLAT']= -90.
        _domain_dic['MAXLAT']=  90.
        _domain_dic['EXCLUDE_FILE']= "''"
        _domain_dic['INCLUDE_FILE']= "''"

        for key,val in _domain_dic.items() :
            if key in self.ExeInputs :
               _domain_dic[key]= self.ExeInputs[key]

        self.domain_def = tempfile.NamedTemporaryFile(mode='w', delete=False)
        self.domain_def.write('&domain_inputs\n')
        for key,val in _domain_dic.items() :
            keyn=(key+"  = ").ljust(16)
            valn = str(val)
            if '_FILE' in key:
               self.domain_def.write(keyn+ "'"+valn+"'"+'\n')
            else :
               self.domain_def.write(keyn+ valn +'\n')
        self.domain_def.write('/\n')
        self.domain_def.close()

    # -----------------------------------------------------------------------------------
    def _check_restart_files(self):
        """ make sure the restart files to be used exist (RESTART=1,2,3) """
        RESTART_str = self.restart

        if RESTART_str in ['1','2','3']:
           y4m2='Y%4d/M%02d' % (self.begDates[0].year, self.begDates[0].month)
           y4m2d2_h2m2='%4d%02d%02d_%02d%02d' % (self.begDates[0].year, self.begDates[0].month,
                                                  self.begDates[0].day,self.begDates[0].hour,self.begDates[0].minute)
           self.in_rstdir = self.ExeInputs['RESTART_PATH']+'/'.join([self.ExeInputs['RESTART_ID'],'output',
                                                                   self.ExeInputs['RESTART_DOMAIN'],'rs',self.ensdirs[0],y4m2])
           if self.with_land:
              tmpFile=self.ExeInputs['RESTART_ID']+'.'+self.catch+'_internal_rst.'+y4m2d2_h2m2
              catchRstFile=self.in_rstdir+'/'+tmpFile
              require(os.path.isfile(catchRstFile), self.catch+'_internal_rst file [%s] does not exist!' %(catchRstFile))
              if RESTART_str == '1':
                 tmpFile=self.ExeInputs['RESTART_ID']+'.landpert_internal_rst.'+y4m2d2_h2m2
                 landpertRstFile=self.in_rstdir+'/'+tmpFile
                 if ( os.path.isfile(landpertRstFile)) :
                    self.has_geos_pert = True

           if self.with_lake:
             tmpFile=self.ExeInputs['RESTART_ID']+'.lake_internal_rst.'+y4m2d2_h2m2
             lakeRstFile=self.in_rstdir+'/'+tmpFile
             require(os.path.isfile(lakeRstFile), 'lake_internal_rst file [%s] does not exist!' %(lakeRstFile))

           if self.with_landice:
             tmpFile=self.ExeInputs['RESTART_ID']+'.landice_internal_rst.'+y4m2d2_h2m2
             landiceRstFile=self.in_rstdir+'/'+tmpFile
             require(os.path.isfile(landiceRstFile), 'landice_internal_rst file [%s] does not exist!' %(landiceRstFile))

           if self.run_route > 0:
              tmpFile=self.ExeInputs['RESTART_ID']+'.route_internal_rst.'+y4m2d2_h2m2
              routeRstFile=self.in_rstdir+'/'+tmpFile
              require(os.path.isfile(routeRstFile), 'route_internal_rst file [%s] does not exist!' %(routeRstFile))

           if self.with_issm:
              tmpFile=self.ExeInputs['RESTART_ID']+'.issm_internal_rst.'+y4m2d2_h2m2
              issmRstFile=self.in_rstdir+'/'+tmpFile
              require(os.path.isfile(issmRstFile), 'issm_internal_rst file [%s] does not exist!' %(issmRstFile))

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

    # ---------------------
    # calculate JobSegments
    # ---------------------
    def _calculateJobSegments(self):
        ## convert date-time strings to datetime object
        ## start/end_time are converted to lists
        ## ensure end>start

        self.begDates=[]
        self.endDates=[]
        self.begDates.append(
                datetime.strptime(
                    self.ExeInputs['BEG_DATE'],
                    '%Y%m%d %H%M%S'
                    )
                )
        self.endDates.append(
                datetime.strptime(
                    self.ExeInputs['END_DATE'],
                    '%Y%m%d %H%M%S'
                    )
                )

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

    def createDirStructure(self):
        """
        Create required dir structure
        """

        status = False

        # run/inp/wrk dirs
        os.makedirs(self.exphome+'/'+self.ExeInputs['EXP_ID'], exist_ok=True)
        os.makedirs(self.rundir, exist_ok=True)
        os.makedirs(self.inpdir, exist_ok=True)
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.scratchdir, exist_ok=True)

        #-start-shorthand-function-
        def _getDirName(outtyp, ensdir, yyyymm):
            return '/'.join([
                    self.outdir,
                    self.ExeInputs['EXP_DOMAIN'],
                    outtyp, # ana/diag/rs/rc_out
                    ensdir,
                    yyyymm
                    ])
        #-end-shorthand-function-

        # met forcing dir
        myMetDir = self.inpdir + '/met_forcing'
        os.makedirs(myMetDir, exist_ok=True)

        # ensxxxx directories
        nSegments = self.nSegments
        for iseg in range(nSegments):
            _start = self.begDates[iseg]
            _end = self.endDates[iseg]

            # Yyyyy/Mmm between StartDateTime and EndDateTime
            newDate = _start
            y4m2_list = [('Y%4d/M%02d' % (newDate.year, newDate.month))]
            while newDate<_end:
                newDate += relativedelta(months=1)
                y4m2_list.append('Y%4d/M%02d' % (newDate.year, newDate.month))

            # ExpDomain/ana/, /diag/ directories
            for ensdir in self.ensdirs_avg:
                for y4m2 in y4m2_list:
                    os.makedirs(_getDirName('ana',  ensdir, y4m2), exist_ok=True)
                    os.makedirs(_getDirName('diag', ensdir, y4m2), exist_ok=True)
            # ExpDomain/rs/ directories
            for ensdir in self.ensdirs:
                for y4m2 in y4m2_list:
                    os.makedirs(_getDirName('rs', ensdir, y4m2), exist_ok=True)

            # ExpDomain/rc_out/ - only for _start
            os.makedirs(_getDirName('rc_out', '', y4m2_list[0]), exist_ok=True)

        # create restart dir
        os.makedirs(self.inpdir + '/restart', exist_ok=True)

        # output dir: create link from ./cat to ./diag for backward compatibility
        os.symlink(_getDirName('', '', 'diag'), _getDirName('', '', 'cat'))

        status = True
        return status


    # -----------------------------------------------------------------------------------
    # createLnRstBc: link/copy build dir, met forcing, tile file, BCs, and restarts
    # -----------------------------------------------------------------------------------

    def createLnRstBc(self) :
        """
        Create links to build dir, met forcing, tile file, BCs, and restarts in input/;
        copy (and zoom in, if necessary) tile/BC files to rc_out/ and restarts to rs/.
        """
        self._write_domain_def()
        tmp_f2g_file = tempfile.NamedTemporaryFile(delete=False)
        tmp_f2g_file.close()
        try:
            os.symlink(self.blddir, self.blddirLn)

            self._link_met_forcing()

            dzsf      = self.ExeInputs.get('SURFLAY',   '50.0')
            # These are dummy values for *cold* restart:
            wemin_in  = self.ExeInputs.get('WEMIN_IN',  '13')    #  WEmin input/output for scale_catch(cn)
            wemin_out = self.ExeInputs.get('WEMIN_OUT', '13')

            tile = self._stage_tile_file(tmp_f2g_file.name, dzsf)

            myTile=self.inpdir+'/tile.data'
            os.symlink(tile,myTile)

            if self.with_land:
               self._stage_land_bcs(tmp_f2g_file.name)

            if self.with_issm:
               self._link_issm_bcs()

            # create and link restart
            print ("Creating and linking restart...")
            myRstDir  = self.inpdir + '/restart/'
            mk_outdir = self.exphome+'/'+self.ExeInputs['EXP_ID']+'/mk_restarts/'

            self._copy_landassim_seeds(myRstDir)

            # (RESTART=M with DO_ISSM=1 or RUN_ROUTE>0 is rejected in __init__)
            if (self.restart in ['2', 'M'] and (self.with_land or self.with_lake or self.with_landice)):
               self._remap_restarts(mk_outdir, dzsf, wemin_in, wemin_out)

            self._stage_restarts(myRstDir, mk_outdir, tmp_f2g_file.name)

            # catch_param restart file (created by preprocess_ldas.x c_f2g)
            catch_param_file = self.rc_out+'/'+ ym_dir(self.begDates[0])+'/'+self.ExeInputs['EXP_ID']+'.ldas_catparam.'+ymd_hm(self.begDates[0])+'z.bin'
            if self.with_land:
               require(os.path.isfile(catch_param_file), "need catch_param file %s" % catch_param_file)

            if self.has_mwrtm :
               self._stage_mwrtm(myRstDir, tmp_f2g_file.name)

            self._write_bcs_info()

            # update 'restart_path' to use relative path from outdir
            print ("Updating restart path...")
            self.ExeInputs['RESTART_PATH'] = myRstDir
        finally:
            for f in [self.domain_def.name, tmp_f2g_file.name]:
                if os.path.exists(f):
                    os.remove(f)
        return True

    # -----------------------------------------------------------------------------------
    def _zoom_or_copy(self, zoom_op, src, dst, f2g_file):
        """ zoom_op = preprocess_ldas.x sub-command (e.g., 'zoomin_catchrst') or None (plain copy) """
        if zoom_op:
            run_cmd(self.bindir + '/preprocess_ldas.x ' + zoom_op + ' ' + src + ' ' + dst + ' ' + f2g_file)
        else:
            shutil.copy(src, dst)

    # -----------------------------------------------------------------------------------
    def _link_met_forcing(self):
        """ link met forcing dir(s) into input/met_forcing; MET_PATH becomes relative to rundir """
        _nens = self.nens
        self.ensemble_forcing = True if self.ExeInputs.get('ENSEMBLE_FORCING', 'NO').upper() == 'YES' else False

        myMetPath =''
        for _i in range(self.first_ens_id, _nens + self.first_ens_id) :
           str_ens = ''
           if ( _nens != 1 and self.ensemble_forcing):
              str_ens = '%03d'%(_i)
           metpath = self.ExeInputs['MET_PATH'].rstrip('/')+str_ens
           myMetDir = self.inpdir + '/met_forcing'
           myMetPath = myMetDir + '/' + metpath.split('/')[-1]
           os.symlink(metpath, myMetPath)
           if ( not self.ensemble_forcing):
              break
        if ( _nens !=1 and self.ensemble_forcing) :
           # replace last three character with '%s"
           self.ExeInputs['MET_PATH'] = os.path.relpath(myMetPath, self.rundir)[:-3]+'%s'
        else:
           self.ExeInputs['MET_PATH'] = os.path.relpath(myMetPath, self.rundir)

    # -----------------------------------------------------------------------------------
    def _stage_tile_file(self, f2g_file, dzsf):
        """
        Copy tile file to rc_out (correcting old EASE format), run c_f2g to determine the
        simulation domain (sets isZoomIn*), and create the domain tile file if needed.
        Returns path of tile file to be used.
        """
        tile       = self.ExeInputs['TILING_FILE']
        short_tile = os.path.basename(tile)
        newtile    = self.rc_out+'/'+short_tile
        shutil.copy(tile, newtile)
        tile       = newtile

        print ('\nCorrect the tile file if it is an old EASE tile format... \n')
        EASEtile=self.rc_out+'/MAPL_'+short_tile
        run_cmd(self.bindir + '/preprocess_ldas.x correctease  '+ tile + ' '+ EASEtile)
        if os.path.isfile(EASEtile) :
            tile=EASEtile

        cmd = self.bindir +'/preprocess_ldas.x c_f2g ' + tile + ' ' + self.domain_def.name + ' '+ self.out_path + ' ' \
            + self.ExeInputs['CATCH_DEF_FILE'] + ' ' + self.ExeInputs['EXP_ID'] + ' ' \
            + self.begDates[0].strftime('%Y%m%d%H%M') + ' '+ dzsf + ' ' + f2g_file + ' ' + '_'.join(self.tile_types)
        print ('Creating f2g file if necessary: '+ f2g_file +'....\n')
        run_cmd(cmd)

        self._read_zoom_info(f2g_file)

        # if running routing or ISSM, make sure domain is global
        if self.run_route > 0 and self.isZoomInLand:
            raise LdasSetupError(f"Must have global domain to run routing model, RUN_ROUTE={self.run_route}")
        if self.with_issm and self.isZoomInLandice:
            raise LdasSetupError("Must have global domain to run ISSM (DO_ISSM: 1)")

        if self.isZoomIn:
            newZoominTile = tile+'.domain'
            print ("\nCreating local tile file :"+ newZoominTile)
            print ("\nAdding 1000 to type of tiles to be excluded from domain...\n")
            self._zoom_or_copy('zoomin_tile', tile, newZoominTile, f2g_file)
            tile = newZoominTile
        return tile

    # -----------------------------------------------------------------------------------
    def _read_zoom_info(self, f2g_file):
        """
        The f2g file is empty if the simulation domain equals the input domain.  Otherwise
        its header lists, for each tile type, the number of tiles in the reduced (simulation)
        and full domains.  Sets isZoomIn, isZoomInLand, isZoomInLake, isZoomInLandice.
        """
        if os.path.getsize(f2g_file) == 0:
            return
        with open(f2g_file, 'r') as f:
            n_types      = int(f.readline().split()[0])
            tile_types   = [int(x) for x in f.readline().split()]
            n_tiles_r    = [int(x) for x in f.readline().split()]
            n_tiles_f    = [int(x) for x in f.readline().split()]

        require(len(tile_types) == n_types and len(n_tiles_r) == n_types and len(n_tiles_f) == n_types,
                "Inconsistent f2g file [%s]" % f2g_file)

        zoom_by_type = { tile_type: (nr != nf) for tile_type, nr, nf in zip(tile_types, n_tiles_r, n_tiles_f) }

        self.isZoomInLand    = zoom_by_type.get(100, False)
        self.isZoomInLake    = zoom_by_type.get(19,  False)
        self.isZoomInLandice = zoom_by_type.get(20,  False)
        self.isZoomIn        = any(zoom_by_type.values())

    # -----------------------------------------------------------------------------------
    def _stage_land_bcs(self, f2g_file):
        """ copy land BC files to rc_out (zoom in if needed) and link them into input/ """
        # (link name, ExeInputs key); optional files only if present
        bcs = [('green', 'GRN_FILE'), ('lai', 'LAI_FILE'), ('ndvi', 'NDVI_FILE'),
               ('nirdf', 'NIRDF_FILE'), ('visdf', 'VISDF_FILE')]
        if self.ExeInputs['LNFM_FILE'] != '':
           bcs.append(('lnfm', 'LNFM_FILE'))
        if self.has_vegopacity:
           bcs.append(('vegopacity', 'VEGOPACITY_FILE'))
        if self.ExeInputs['EASE_PFAF_TILE_FILE'] != '':
           bcs.append(('EASE_pfaf_tile_file.nc4', 'EASE_PFAF_TILE_FILE'))

        if self.isZoomInLand:
           print ("Creating the boundary files for the simulation domain...\n")
        staged = []
        for bcln, key in bcs:
           bcf = self.rc_out+'/'+os.path.basename(self.ExeInputs[key])
           shutil.copy(self.ExeInputs[key], self.rc_out+'/')
           if self.isZoomInLand:
              self._zoom_or_copy('zoomin_bc', bcf, bcf+'.domain', f2g_file)
              bcf = bcf+'.domain'
           staged.append((bcln, bcf))

        print ("linking bcs...")
        for bcln, bc in staged:
           myBC = self.inpdir+'/'+bcln if '.nc4' in bcln else self.inpdir+'/'+bcln+'.data'
           os.symlink(bc,myBC)

        if ("catchcn" in self.catch):
           os.symlink(self.bcs_dir_landshared + 'CO2_MonthlyMean_DiurnalCycle.nc4', \
                      self.inpdir+'/CO2_MonthlyMean_DiurnalCycle.nc4')

    # -----------------------------------------------------------------------------------
    def _link_issm_bcs(self):
        for pattern in ('*.bin', '*.toolkits', '*.nc*'):
           for issmbc in glob.glob(self.bcs_dir_landiceshared + pattern):
              os.symlink(issmbc, self.inpdir+'/'+os.path.basename(issmbc))

    # -----------------------------------------------------------------------------------
    def _rst_src_dir(self):
        """ rs/ dir of the restart experiment (with trailing '/') """
        return self.ExeInputs['RESTART_PATH'] + self.ExeInputs['RESTART_ID'] + '/output/'+ self.ExeInputs['RESTART_DOMAIN'] + '/rs/'

    # -----------------------------------------------------------------------------------
    def _copy_landassim_seeds(self, myRstDir):
        """ copy the landassim obs perturbation seed restarts if they exist (and LAND_ASSIM) """
        start   = self.begDates[0]
        rstid   = self.ExeInputs['RESTART_ID']
        exp_id  = self.ExeInputs['EXP_ID']
        rstpath = self._rst_src_dir()
        for iens in range(self.nens) :
           _ensdir = self.ensdirs[iens]
           _ensid  = self.ensids[iens]
           landassim_seeds = rstpath + _ensdir + '/' + ym_dir(start)+'/' + rstid + '.landassim_obspertrseed_rst.'+ymd_hm(start)
           if os.path.isfile(landassim_seeds) and self.assim :
               _seeds = self.rstdir + _ensdir + '/' + ym_dir(start)+'/' + exp_id + '.landassim_obspertrseed_rst.'+ymd_hm(start)
               shutil.copy(landassim_seeds, _seeds)
               os.symlink(_seeds, myRstDir+ '/landassim_obspertrseed'+ _ensid +'_rst')
               self.has_landassim_seed = True

    # -----------------------------------------------------------------------------------
    def _remap_restarts(self, mk_outdir, dzsf, wemin_in, wemin_out):
        """ RESTART=2 or M: remap restarts into mk_outdir (uses remap_restarts package) """
        RESTART_str = self.restart
        YYYYMMDDHH  = self.begDates[0].strftime('%Y%m%d%H')
        bc_base     = os.path.dirname(self.ExeInputs['BCS_PATH'].rstrip('/'))
        bc_version  = bcs_version_from_path(self.ExeInputs['BCS_PATH'])
        answers={}
        if RESTART_str == 'M' :
           answers['input:shared:MERRA-2']   = True
           answers['input:shared:yyyymmddhh']= YYYYMMDDHH
           answers['output:shared:out_dir']  = mk_outdir
           init_merra2(answers)
        config = get_config_from_answers(answers, config_tpl = True)

        config['slurm_pbs']['account'] = self.RmInputs['account']
        config['slurm_pbs']['qos'] = 'debug'

        if RESTART_str == '2' :
           config['input']['surface']['catch_tilefile']  = self.in_tilefile
           config['input']['shared']['expid']            = self.ExeInputs['RESTART_ID']
           config['input']['shared']['yyyymmddhh']       = YYYYMMDDHH
           config['input']['shared']['rst_dir']          = self.in_rstdir
           config['input']['surface']['wemin']           = wemin_in
           config['input']['surface']['catch_model']     = self.catch

        config['output']['shared']['out_dir']         = mk_outdir
        config['output']['surface']['catch_remap']    = True
        config['output']['surface']['catch_tilefile'] = self.ExeInputs['TILING_FILE']
        config['output']['shared']['bc_base']         = bc_base
        config['output']['shared']['bc_version']      = bc_version
        config['output']['surface']['EASE_grid']      = self.ExeInputs['BCS_RESOLUTION']

        config['output']['shared']['expid']           = self.ExeInputs['EXP_ID']
        config['output']['surface']['surflay']        = dzsf
        config['output']['surface']['wemin']          = wemin_out

        if self.with_land:
          catch_obj = catchANDcn(config_obj = config)
          catch_obj.remap()
        if self.with_lake or self.with_landice:
            config['output']['surface']['remap_water'] = True
            config['input']['surface']['zoom'] = '2'
            other_rst_obj = other_restarts(config_obj = config)
            other_rst_obj.remap()

    # -----------------------------------------------------------------------------------
    def _restart_components(self, mk_outdir):
        """
        Table of restart components to be staged (in this order) for each ensemble member.

        Keys:
          name      : component name; link is input/restart/<name><ensid>_internal_rst
          src       : function(ensdir) -> source file (restart exp, mk_restarts, or BCS)
          local     : function(ensdir) -> local copy in output/.../rs/
          zoom      : preprocess_ldas.x zoom sub-command, or None (plain copy)
          first_only: only for first ensemble member (single-member components)
          fallback0 : if source is missing, re-use local file of first member
          required0 : source must exist for first member (else error)
          always    : copy without checking that the source exists
          link      : create the link
        If the source does not exist (and no fallback0), the link points to the source path.
        """
        R      = self.restart
        start  = self.begDates[0]
        y4m2   = ym_dir(start)
        stamp  = ymd_hm(start)
        rstid  = self.ExeInputs['RESTART_ID']
        exp_id = self.ExeInputs['EXP_ID']
        rstpath= self._rst_src_dir()
        YYYYMMDD = start.strftime('%Y%m%d')

        def src_rs(comp):
            return lambda ensdir: rstpath+ensdir+'/'+y4m2+'/'+rstid+'.'+comp+'_internal_rst.'+stamp
        def local_rs(comp):
            return lambda ensdir: self.rstdir+ensdir+'/'+y4m2+'/'+exp_id+'.'+comp+'_internal_rst.'+stamp
        def mk(comp, what):
            return lambda ensdir: glob_first(mk_outdir+'*'+comp+'_internal_rst.'+YYYYMMDD+'*', what)
        def vegdyn_bcs(ensdir):
            return glob_first(self.bcs_dir_land + 'vegdyn_*.dat', 'vegdyn BCS file')
        def nofile(ensdir):
            return ''

        def comp(name, src, local=None, zoom=None, first_only=False, fallback0=False,
                 required0=False, always=False, link=True):
            return dict(name=name, src=src, local=local or local_rs(name), zoom=zoom,
                        first_only=first_only, fallback0=fallback0, required0=required0,
                        always=always, link=link)

        C = []
        if self.with_land:
            C.append(comp(self.catch,
                          mk(self.catch, 'remapped '+self.catch+' restart in mk_restarts') if R in ['2','M'] else src_rs(self.catch),
                          zoom='zoomin_catchrst' if self.isZoomInLand else None, fallback0=True))
            C.append(comp('vegdyn',
                          vegdyn_bcs if R in ['2','M','3'] else (lambda ensdir: rstpath+ensdir+'/'+rstid+'.vegdyn_internal_rst'),
                          local=lambda ensdir: self.rstdir+ensdir+'/'+exp_id+'.vegdyn_internal_rst',
                          zoom='zoomin_vegrst' if self.isZoomInLand else None, fallback0=True))
        if self.with_lake:
            C.append(comp('lake', mk('lake', 'lake restart in mk_restarts') if R in ['2','M'] else src_rs('lake'),
                          zoom='zoomin_lakerst' if self.isZoomInLake else None, first_only=True))
        if self.with_landice:
            C.append(comp('landice', mk('landice', 'landice restart in mk_restarts') if R in ['2','M'] else src_rs('landice'),
                          zoom='zoomin_landicerst' if self.isZoomInLandice else None, first_only=True))
        if self.with_issm:
            issm_src = {'1': src_rs('issm'), '3': src_rs('issm'), 'M': nofile,
                        '2': mk('issm', 'ISSM restart in mk_restarts (RESTART=2 with DO_ISSM=1)')}[R]
            C.append(comp('issm', issm_src, first_only=True, link=(R in ['1','2','3'])))
        if self.run_route > 0:
            C.append(comp('route', mk('route', 'route restart in mk_restarts') if R == '2' else src_rs('route'),
                          fallback0=True, required0=True))
        if self.has_geos_pert and self.perturb == 1:
            C.append(comp('landpert', src_rs('landpert'), always=True))
        return C

    # -----------------------------------------------------------------------------------
    def _stage_restarts(self, myRstDir, mk_outdir, f2g_file):
        """ copy (or zoom in) restarts to output/.../rs/ and link them into input/restart/ """
        print ("restart: " + self.restart)
        member0 = {}                                   # name -> local file of first member
        for C in self._restart_components(mk_outdir):
            for iens in range(self.nens) :
                if C['first_only'] and iens > 0:
                    break
                ensdir = self.ensdirs[iens]
                src    = C['src'](ensdir)
                if C['required0'] and iens == 0 and not os.path.isfile(src):
                    raise LdasSetupError("Please copy a " + C['name'] + " restart (e.g., from " + ROUTE_RST_DIR + ") to " + src)
                if C['always'] or os.path.isfile(src):
                    target = C['local'](ensdir)
                    self._zoom_or_copy(C['zoom'], src, target, f2g_file)
                    if iens == 0:
                        member0[C['name']] = target
                elif C['fallback0']:
                    target = member0.get(C['name'], '')
                else:
                    target = src
                if C['link']:
                    link = myRstDir+'/'+ C['name'] + self.ensids[iens] +'_internal_rst'
                    print ("link " + C['name'] + " restart: " + link + " -> " + target)
                    os.symlink(target, link)

    # -----------------------------------------------------------------------------------
    def _stage_mwrtm(self, myRstDir, f2g_file):
        start = self.begDates[0]
        mwRTMLocal = self.rc_out+'/'+ ym_dir(start)+'/'+self.ExeInputs['EXP_ID']+'.ldas_mwRTMparam.'+ymd_hm(start)+'z.nc4'
        if self.isZoomInLand:
           print ("Creating the local mwRTM restart file... \n")
        self._zoom_or_copy('zoomin_mwrtmrst' if self.isZoomInLand else None, self.mwrtm_file, mwRTMLocal, f2g_file)
        os.symlink(mwRTMLocal, myRstDir+'/mwrtm_param_rst')

    # -----------------------------------------------------------------------------------
    def _write_bcs_info(self):
        """ BCS_info.txt facilitates BCS consistency check when exp is later used as restart for another exp """
        with open(self.rc_out+'/BCS_info.txt','wt') as fout :
           fout.write("BCS_VERSION: "    + bcs_version_from_path(self.ExeInputs['BCS_PATH']) + '\n')
           fout.write("BCS_RESOLUTION: " + self.ExeInputs['BCS_RESOLUTION'])

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
        return True

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

    # -----------------------------------------------------------------------------------

    def createBatchRun(self):
        """
        """

        status = False

        os.chdir(self.rundir)
        fout =open(self.rundir+'/ldas_batchrun.j','w')
        fout.write("#!/bin/bash -f\n")
        SBATCHQSUB = 'sbatch'
        if self.GEOS_SITE == 'NAS':
          SBATCHQSUB = 'qsub'
        fout.write("\nsed -i 's/if($capdate<$enddate) "+SBATCHQSUB+"/#if($capdate<$enddate) "+SBATCHQSUB+"/g' lenkf.j\n\n")
        nSegments = self.nSegments
        for iseg in range(nSegments):
            if self.GEOS_SITE == 'NAS':
                # PBS: qsub prints the job id only
                if iseg ==0 :
                    fout.write("jobid%d=$(qsub lenkf.j)\n"%(iseg))
                else :
                    fout.write("jobid%d=$(qsub -W depend=afterok:$jobid%d lenkf.j)\n"%(iseg,iseg-1))
                fout.write("echo $jobid%d\n"%iseg )
                continue
            if iseg ==0 :
                fout.write("jobid%d=$(echo $(sbatch lenkf.j) | cut -d' ' -f 4)\n"%(iseg))
                fout.write("echo $jobid%d\n"%iseg )
            else :
                #fout.write("jobid%d=$(echo $(sbatch --dependency=afterany:$jobid%d --output=%s --error=%s lenkf.j) | cut -d' ' -f 4)\n"%(iseg,iseg-1,_logfile, _errfile))
                fout.write("jobid%d=$(echo $(sbatch --dependency=afterok:$jobid%d lenkf.j) | cut -d' ' -f 4)\n"%(iseg,iseg-1))
                fout.write("echo $jobid%d\n"%iseg )
        fout.write("\nsed -i 's/#if($capdate<$enddate) "+SBATCHQSUB+"/if($capdate<$enddate) "+SBATCHQSUB+"/g' lenkf.j\n\n")
        fout.close()

        sp.call(['chmod', '755', self.rundir+'/ldas_batchrun.j'])
        status = True
        return status

    # -----------------------------------------------------------------------------------

    def createRunScripts(self):
        """
        """

        status = False

        os.chdir(self.rundir)

        my_qos='allnccs'
        if self.GEOS_SITE == 'NAS': my_qos = 'normal'
        if 'qos' in self.RmInputs :
            my_qos = self.RmInputs['qos']

        my_job=self.ExeInputs['EXP_ID']
        if 'job_name' in self.RmInputs :
            my_job = self.RmInputs['job_name']

        start = self.begDates[0]
        expid = self.ExeInputs['EXP_ID']
        myDateTime = '%04d%02d%02d_%02d%02dz' % \
                (start.year, start.month, start.day,start.hour,start.minute)
        my_logfile = os.path.relpath(
                '/'.join([
                        self.outdir,
                        self.ExeInputs['EXP_DOMAIN'],
                        'rc_out',
                        'Y%04d' % start.year,
                        'M%02d' % start.month,
                        '.'.join([expid, 'ldas_log', myDateTime, 'txt']),
                        ]),
                self.rundir)
        my_errfile = os.path.relpath(
                '/'.join([
                        self.outdir,
                        self.ExeInputs['EXP_DOMAIN'],
                        'rc_out',
                        'Y%04d' % start.year,
                        'M%02d' % start.month,
                        '.'.join([expid, 'ldas_err', myDateTime, 'txt']),
                        ]),
                self.rundir)

        constraint = '"[mil|cas]"'
        if self.GEOS_SITE == "NAS" :
            constraint = 'cas_ait'

        if 'constraint' in self.RmInputs:
            constraint = self.RmInputs['constraint']

        my_nodes=''
        if 'ntasks-per-node' in self.RmInputs:
           ntasks_per_node = int(self.RmInputs['ntasks-per-node'])
           ntasks          = int(self.RmInputs['ntasks_model'])
           require(ntasks%ntasks_per_node == 0, 'Please make ntasks_model a multiple of ntasks-per-node')
           nodes           = ntasks//ntasks_per_node
           my_nodes  = '#SBATCH --nodes=' + str(nodes) +' --ntasks-per-node=' + str(ntasks_per_node)
           if (ntasks_per_node > 46):
              require(constraint != 'cas', "Make sure constraint is compatible with ntasks-per-node")

        SBATCHQSUB = 'sbatch'
        if self.GEOS_SITE == 'NAS':
            SBATCHQSUB = 'qsub'

        DETECTED_MPI_STACK = "@MPI_STACK@"

        job_head = job_directive[self.GEOS_SITE]
        lenkf_str= (job_head+job_body).format(
                  SBATCHQSUB         = SBATCHQSUB,
                  MY_ACCOUNT         = self.RmInputs['account'],
                  MY_WALLTIME        = self.RmInputs['walltime'],
                  MY_NTASKS_MODEL    = str(self.RmInputs['ntasks_model']),
                  MY_NODES           = my_nodes,
                  MY_CONSTRAINT      = constraint,
                  MY_OSERVER_NODES   = str(self.RmInputs['oserver_nodes']),
                  MY_WRITERS_NPES    = str(self.RmInputs['writers-per-node']),
                  MY_QOS             = my_qos,
                  MY_JOB             = my_job,
                  MY_EXPID           = self.ExeInputs['EXP_ID'],
                  MY_EXPDOMAIN       = self.ExeInputs['EXP_DOMAIN'],
                  MY_LOGFILE         = my_logfile,
                  MY_ERRFILE         = my_errfile,
                  MY_LANDMODEL       = self.catch,
                  MY_POSTPROC_HIST   = str(self.ExeInputs['POSTPROC_HIST']),
                  MY_FIRST_ENS_ID    = str(self.first_ens_id),
                  MY_LADAS_COUPLING  = str(self.ladas_cpl),
                  MY_ENSEMBLE_FORCING= self.ExeInputs.get('ENSEMBLE_FORCING', 'NO').upper(),
                  MY_ADAS_EXPDIR     = self.adas_expdir,
                  MY_EXPDIR          = self.expdir,
                  DETECTED_MPI_STACK = DETECTED_MPI_STACK,
                  )

        with open('lenkf.j','wt') as fout :
           fout.write(lenkf_str)
        sp.call(['chmod', '755', 'lenkf.j'])

        print ('\nExperiment directory: %s' % self.expdir)
        print ()
        status = True
        return status
