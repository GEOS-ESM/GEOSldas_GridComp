#!/usr/bin/env python3

"""
ldas_staging.py: experiment directory structure, links and copies of met forcing,
tile file, BCs, and restarts (mixin class for ldas.ldas)
"""

import os
import glob
import shutil
import tempfile
from collections                  import OrderedDict
from dateutil.relativedelta       import relativedelta
from remap_utils                  import init_merra2, get_config_from_answers
from remap_other_restarts         import other_restarts
from remap_catchANDcn             import catchANDcn
from setup_utils                  import (LdasSetupError, require, run_cmd, glob_first,
                                          bcs_version_from_path, ym_dir, ymd_hm)
from ldas_config                  import ROUTE_RST_DIR


class LdasStaging:
    """ mixin: createDirStructure and createLnRstBc for class ldas (see ldas.py) """

    # -----------------------------------------------------------------------------------
    def createDirStructure(self):
        """
        Create required dir structure
        """
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
            y4m2_list = [ym_dir(newDate)]
            while newDate<_end:
                newDate += relativedelta(months=1)
                y4m2_list.append(ym_dir(newDate))

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
