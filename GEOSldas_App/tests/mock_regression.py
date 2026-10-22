#!/usr/bin/env python3
"""
Local end-to-end regression harness for ldas_setup (no Discover / no build needed).

It builds two fake installations ("ref" from a git commit, "new" from the working
tree), a fake BCS / restart / met-forcing data tree, and runs `ldas_setup setup`
for a set of cases with both installations.  The resulting experiment directories
are compared with compare_exp.py.

Faked pieces:
  - preprocess_ldas.x : python stub (c_f2g, zoomin_*, optimize, correctease)
  - GEOSldas.x        : empty file
  - remap_utils / remap_catchANDcn / remap_other_restarts : python stubs that
                        write remapped restarts into mk_restarts/
  - BCS / restart / MWRTM files : small text files
Real pieces: ldas_setup, ldas.py, setup_utils.py, lenkf_j_template.py,
             process_hist.csh (needs csh), GEOSldas_*.rc, nml files.

usage (from GEOSldas_App):
   python3 tests/mock_regression.py [--ref COMMIT] [--keep] [--case NAME ...] [-v]

Default reference commit: REF_COMMIT below.
Exit status: 0 if every case matches its expectation, 1 otherwise.

NOT covered: coupled land-atm DAS (ladas_cpl>0; uses hard-wired /discover paths),
NAS site, real Fortran preprocessing (zoom, optimize), real restart remapping.
"""

import os
import sys
import shutil
import subprocess as sp
import tempfile
import textwrap

HERE    = os.path.dirname(os.path.abspath(__file__))
APPDIR  = os.path.dirname(HERE)                       # GEOSldas_App
REPO    = os.path.dirname(APPDIR)                     # @GEOSldas_GridComp
SURF_RC = os.path.join(REPO, '..', '@GEOSgcm_GridComp', 'GEOSagcm_GridComp',
                       'GEOSphysics_GridComp', 'GEOSsurface_GridComp', 'Shared',
                       'GEOS_SurfaceGridComp.rc')

REF_COMMIT = '2ec60a8'            # end of phase 5; later refactoring must be identical to it
# (history: ed23aa6 = end of phase 2 bug fixes)
# (use --ref 0f6155f to compare against the original pre-refactor code)

sys.path.insert(0, HERE)
from compare_exp import compare  # noqa: E402

APP_FILES = ['ldas_setup', 'ldas.py', 'setup_utils.py', 'lenkf_j_template.py', 'process_hist.csh']

# ---------------------------------------------------------------------------------
# stubs
# ---------------------------------------------------------------------------------

PREPROCESS_STUB = r'''#!/usr/bin/env python3
import sys, os, shutil
op = sys.argv[1]
a  = sys.argv[2:]
if op == 'correctease':
    pass
elif op == 'c_f2g':
    tile, domain_def, out_path, catchdef, expid, ymdhm, dzsf, f2g, types = a
    txt = open(domain_def).read()
    vals = {}
    for line in txt.splitlines():
        if '=' in line:
            k, v = line.split('=', 1)
            vals[k.strip()] = v.strip()
    zoom = (float(vals['MINLON']), float(vals['MAXLON']), float(vals['MINLAT']),
            float(vals['MAXLAT'])) != (-180., 180., -90., 90.)
    if zoom:
        tts = [int(t) for t in types.split('_')]
        with open(f2g, 'w') as f:
            f.write('%d\n' % len(tts))
            f.write(' '.join(str(t) for t in tts) + '\n')
            f.write(' '.join('5' for t in tts) + '\n')
            f.write(' '.join('10' if t == 100 else '5' for t in tts) + '\n')
    y, m = ymdhm[0:4], ymdhm[4:6]
    d = os.path.join(out_path, 'rc_out', 'Y'+y, 'M'+m)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, expid+'.ldas_catparam.'+ymdhm[0:8]+'_'+ymdhm[8:12]+'z.bin'), 'w').write('catparam\n')
elif op.startswith('zoomin_'):
    src, dst = a[0], a[1]
    with open(dst, 'w') as f:
        f.write(open(src).read() + op + '\n')
elif op == 'optimize':
    tile, ntasks, outfile, rundir = a[0], a[1], a[2], a[3]
    with open(outfile, 'w') as f:
        f.write('NX: 1\nNY: %s\n' % ntasks)
    open(os.path.join(rundir, 'IMS.rc'), 'w').write('ims\n')
else:
    sys.exit('unknown op ' + op)
'''

REMAP_UTILS_STUB = '''
STRETCH_GRID = {'SG001': [2.5, 39.5, -98.35], 'SG002': [3.0, 39.5, -98.35]}
def init_merra2(answers):
    answers['input:shared:rst_dir'] = 'MERRA2'
def get_config_from_answers(answers, config_tpl=True):
    return {'slurm_pbs': {}, 'input':  {'surface': {}, 'shared': dict(answers)},
            'output': {'surface': {}, 'shared': {}}}
'''

REMAP_CATCH_STUB = '''
import os
class catchANDcn:
    def __init__(self, config_obj):
        self.c = config_obj
    def remap(self):
        o   = self.c['output']['shared']
        inp = self.c['input']
        ymdh = inp['shared'].get('yyyymmddhh') or inp['shared'].get('input:shared:yyyymmddhh')
        model = inp['surface'].get('catch_model', 'catch')
        os.makedirs(o['out_dir'], exist_ok=True)
        fn = os.path.join(o['out_dir'], o['expid'] + '.' + model + '_internal_rst.' + ymdh[0:8] + '_' + ymdh[8:10] + 'z.nc4')
        open(fn, 'w').write('remapped catch %s\\n' % sorted(self.c['output']['surface'].items()))
'''

REMAP_OTHER_STUB = '''
import os
class other_restarts:
    def __init__(self, config_obj):
        self.c = config_obj
    def remap(self):
        o   = self.c['output']['shared']
        inp = self.c['input']['shared']
        ymdh = inp.get('yyyymmddhh') or inp.get('input:shared:yyyymmddhh')
        os.makedirs(o['out_dir'], exist_ok=True)
        for t in ['lake', 'landice', 'issm', 'route']:
            fn = os.path.join(o['out_dir'], o['expid'] + '.' + t + '_internal_rst.' + ymdh[0:8] + '_' + ymdh[8:10] + 'z.nc4')
            open(fn, 'w').write('remapped ' + t + '\\n')
'''

# ---------------------------------------------------------------------------------
# fake installation
# ---------------------------------------------------------------------------------

def git_show(commit, relpath):
    return sp.check_output(['git', 'show', '%s:GEOSldas_App/%s' % (commit, relpath)], cwd=REPO)


def list_py_modules(commit):
    """all GEOSldas_App/*.py files in the commit (or working tree if commit is None)"""
    if commit is None:
        return [f for f in os.listdir(APPDIR) if f.endswith('.py')]
    out = sp.check_output(['git', 'ls-tree', '--name-only', commit, 'GEOSldas_App/'], cwd=REPO).decode()
    return [os.path.basename(l) for l in out.split() if l.endswith('.py')]


DISCOVER_BCS   = b"/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles"
DISCOVER_MWRTM = b"/discover/nobackup/projects/gmao/smap/LDAS_inputs_for_LADAS/RTM_params/RTMParam_SMAP_L4SM_v006/"


def make_install(dest, commit, site='NCCS', data_root=None):
    """commit=None -> working tree.  Hard-wired /discover paths (coupled mode) are
    redirected to the fake data tree under data_root."""
    bindir = os.path.join(dest, 'bin')
    etcdir = os.path.join(dest, 'etc')
    os.makedirs(bindir); os.makedirs(etcdir)
    files = sorted(set(APP_FILES) | set(list_py_modules(commit)))
    for f in files:
        if commit is None:
            data = open(os.path.join(APPDIR, f), 'rb').read()
        else:
            data = git_show(commit, f)
        # emulate cmake configure_file(@ONLY) for all files (only ldas.py has placeholders)
        data = data.replace(b'@GEOS_SITE@', site.encode()).replace(b'@MPI_STACK@', b'openmpi')
        if data_root:
            data = data.replace(DISCOVER_BCS,   os.path.join(data_root, 'bcs').encode())
            data = data.replace(DISCOVER_MWRTM, os.path.join(data_root, 'mwrtm').encode() + b'/')
        with open(os.path.join(bindir, f), 'wb') as fo:
            fo.write(data)
        os.chmod(os.path.join(bindir, f), 0o755)
    for name, txt in [('preprocess_ldas.x', PREPROCESS_STUB),
                      ('remap_utils.py', REMAP_UTILS_STUB),
                      ('remap_catchANDcn.py', REMAP_CATCH_STUB),
                      ('remap_other_restarts.py', REMAP_OTHER_STUB),
                      ('GEOSldas.x', '')]:
        with open(os.path.join(bindir, name), 'w') as fo:
            fo.write(txt)
        os.chmod(os.path.join(bindir, name), 0o755)
    # etc: rc templates and nml (always from the working tree; not part of the refactor)
    for f in os.listdir(APPDIR):
        if (f.startswith('GEOSldas_') and f.endswith('.rc')) or f.startswith('LDASsa_DEFAULT_inputs_'):
            shutil.copy(os.path.join(APPDIR, f), etcdir)
    for f in os.listdir(os.path.join(APPDIR, 'LandAtmDAS_nml')):
        shutil.copy(os.path.join(APPDIR, 'LandAtmDAS_nml', f), etcdir)
    shutil.copy(SURF_RC, etcdir)
    return bindir

# ---------------------------------------------------------------------------------
# fake data
# ---------------------------------------------------------------------------------

BCS_VERSION = 'v12'
BCS_RES     = 'EASEv2_M36'
GRIDNAME    = 'EASEv2_M36'
BEG         = ('2020', '01', '01', '00', '00')


def w(path, text='x\n'):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(text)


def make_data(root):
    bcs = os.path.join(root, 'bcs', BCS_VERSION)
    land = os.path.join(bcs, 'land', BCS_RES)
    geom = os.path.join(bcs, 'geometry', BCS_RES)
    tile_txt = 'header1\nheader2\n%s\nrest\n' % GRIDNAME
    w(os.path.join(geom, BCS_RES + '.til'), tile_txt)
    w(os.path.join(land, 'clsm', 'catchment.def'), 'catchdef\n')
    for f in ['green_clim_x.data', 'lai_clim_x.data', 'ndvi_clim_x.data', 'nirdf_x.dat',
              'visdf_x.dat', 'lnfm_clim_x.data', 'vegdyn_x.dat']:
        w(os.path.join(land, f), f + '\n')
    w(os.path.join(bcs, 'land', 'shared', 'CO2_MonthlyMean_DiurnalCycle.nc4'), 'co2\n')
    for f in ['a.bin', 'b.toolkits', 'c.nc4']:
        w(os.path.join(bcs, 'landice', 'shared', f), f + '\n')

    # restart experiment "rst_exp" with ens0000, ens0001
    y, m, d, h, mi = BEG
    rpath = os.path.join(root, 'restarts')
    rid, rdom = 'rst_exp', 'GLOBAL'
    rc_out = os.path.join(rpath, rid, 'output', rdom, 'rc_out')
    w(os.path.join(rc_out, BCS_RES + '.til'), tile_txt)
    for f in ['green_clim_x.data', 'lai_clim_x.data', 'ndvi_clim_x.data', 'nirdf_x.dat',
              'visdf_x.dat', 'lnfm_clim_x.data']:
        w(os.path.join(rc_out, f), 'rst ' + f + '\n')
    w(os.path.join(rc_out, 'BCS_info.txt'), 'BCS_VERSION: %s\nBCS_RESOLUTION: %s' % (BCS_VERSION, BCS_RES))
    stamp = '%s%s%s_%s%s' % (y, m, d, h, mi)
    for ens in ['ens0000', 'ens0001']:
        rs = os.path.join(rpath, rid, 'output', rdom, 'rs', ens)
        w(os.path.join(rs, rid + '.vegdyn_internal_rst'), 'vegdyn ' + ens + '\n')
        for comp in ['catch', 'catchcnclm40', 'lake', 'landice', 'issm', 'route', 'landpert',
                     'landassim_obspertrseed']:
            suffix = '_rst.' if comp == 'landassim_obspertrseed' else '_internal_rst.'
            w(os.path.join(rs, 'Y' + y, 'M' + m, rid + '.' + comp + suffix + stamp), comp + ' ' + ens + '\n')

    # coupled land-atm DAS (ladas_cpl): CF0012 BCS, LDAS restarts at rstloc with ens0001, ens0002
    cres = 'CF0012x6C_CF0012x6C'
    cdom = 'CF0012x6C_GLOBAL'
    ctile = 'header1\nheader2\nPE12x72-CF\nrest\n'
    w(os.path.join(bcs, 'geometry', cres, cres + '.til'), ctile)
    w(os.path.join(bcs, 'land', cres, 'clsm', 'catchment.def'), 'catchdef\n')
    w(os.path.join(bcs, 'land', cres, 'vegdyn_x.dat'), 'vegdyn\n')
    crid = 'ldas_rst'
    crc = os.path.join(root, 'adas_rst', crid, 'output', cdom, 'rc_out')
    w(os.path.join(crc, cres + '.til'), ctile)
    for f in ['green_clim_x.data', 'lai_clim_x.data', 'ndvi_clim_x.data', 'nirdf_x.dat', 'visdf_x.dat']:
        w(os.path.join(crc, f), 'cpl ' + f + '\n')
    for ens in ['ens0001', 'ens0002']:
        rs = os.path.join(root, 'adas_rst', crid, 'output', cdom, 'rs', ens)
        w(os.path.join(rs, crid + '.vegdyn_internal_rst'), 'vegdyn ' + ens + '\n')
        w(os.path.join(rs, 'Y2020', 'M01', crid + '.catch_internal_rst.20200101_2100'), 'catch ' + ens + '\n')
    w(os.path.join(root, 'mwrtm', cres, 'mwRTM_param.nc4'), 'mwrtm cpl\n')

    w(os.path.join(root, 'mwrtm', BCS_RES, 'mwRTM_param.nc4'), 'mwrtm\n')
    w(os.path.join(root, 'mwrtm', BCS_RES, 'vegopacity.bin'), 'vegop\n')
    os.makedirs(os.path.join(root, 'met', 'MERRA2'), exist_ok=True)
    for i in range(3):
        os.makedirs(os.path.join(root, 'met', 'ensmet%03d' % i), exist_ok=True)
    w(os.path.join(root, 'nml', 'LDASsa_SPECIAL_inputs_ensupd.nml'), '&ens_upd_inputs\n/\n')
    return dict(bcs=bcs, rpath=rpath, rid=rid, rdom=rdom, root=root,
                rstloc=os.path.join(root, 'adas_rst', crid))

# ---------------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------------

def base_exeinp(D):
    return {
        'EXP_ID': 'testexp', 'EXP_DOMAIN': 'GLOBAL', 'NUM_LDAS_ENSEMBLE': '1',
        'BEG_DATE': '20200101 000000', 'END_DATE': '20200201 000000',
        'RESTART': '1', 'RESTART_ID': D['rid'], 'RESTART_PATH': D['rpath'],
        'RESTART_DOMAIN': D['rdom'],
        'MET_TAG': 'M2COR_cross', 'MET_PATH': os.path.join(D['root'], 'met', 'MERRA2'),
        'FORCE_DTSTEP': '3600',
        'BCS_PATH': D['bcs'], 'BCS_RESOLUTION': BCS_RES,
        'DO_ISSM': '0',
    }


def cases(D):
    C = {}
    def add(name, upd=(), drop=(), batinp=None, expect='same', site='NCCS', cpl=None):
        e = base_exeinp(D)
        e.update(dict(upd))
        for k in drop:
            e.pop(k, None)
        b = {'account': 's1234', 'walltime': '01:00:00', 'ntasks_model': '126'}
        b.update(batinp or {})
        C[name] = dict(exeinp=e, batinp=b, expect=expect, site=site, cpl=cpl)

    add('r1_global')
    add('r1_zoom', {'MINLAT': '20.', 'MAXLAT': '50.', 'MINLON': '-120.', 'MAXLON': '-60.'})
    add('r1_ens_assim', {'NUM_LDAS_ENSEMBLE': '2', 'LAND_ASSIM': 'YES', 'PERTURBATIONS': '1',
                         'MWRTM_PATH': os.path.join(D['root'], 'mwrtm'),
                         'ENSEMBLE_FORCING': 'YES',
                         'MET_PATH': os.path.join(D['root'], 'met', 'ensmet')})
    add('r1_zoom_ens_assim', {'NUM_LDAS_ENSEMBLE': '2', 'LAND_ASSIM': 'YES',
                              'MWRTM_PATH': os.path.join(D['root'], 'mwrtm'),
                              'MINLAT': '20.', 'MAXLAT': '50.'})
    add('r3_catchcn', {'RESTART': '3', 'LSM_CHOICE': '2'})
    add('r1_lake_landice_route', {'TILE_TYPES': '100 19 20', 'RUN_ROUTE': '1'})
    add('r1_lake_landice_zoom', {'TILE_TYPES': '100 19 20', 'MINLAT': '20.'})
    add('r1_issm', {'TILE_TYPES': '100 20', 'DO_ISSM': '1'})
    add('r2_remap', {'RESTART': '2', 'TILE_TYPES': '100 19 20'})
    add('rM_merra2', {'RESTART': 'M'})
    add('multi_sgmt', {'END_DATE': '20200501 000000', 'JOB_SGMT': '00000100 000000', 'NUM_SGMT': '1'})
    add('job_sgmt_days', {'END_DATE': '20200301 000000', 'JOB_SGMT': '00000015 000000', 'NUM_SGMT': '2'})
    add('nml_hist_custom', {'NML_INPUT_PATH': os.path.join(D['root'], 'nml'),
                            'POSTPROC_HIST': '1', 'RUN_IRRIG': '1', 'AEROSOL_DEPOSITION': '1'},
        batinp={'qos': 'debug', 'job_name': 'myjob', 'ntasks-per-node': '42', 'constraint': 'mil',
                'oserver_nodes': '1', 'writers-per-node': '3'})
    add('tile_txt', {'TILE_FILE_FORMAT': 'TXT'})
    # restart staging edge cases (phase 4)
    add('r3_ens3_fallback', {'RESTART': '3', 'NUM_LDAS_ENSEMBLE': '3'})          # no ens0002 -> member-0 files
    add('r1_ens2_all_types', {'NUM_LDAS_ENSEMBLE': '2', 'TILE_TYPES': '100 19 20', 'DO_ISSM': '1',
                              'RUN_ROUTE': '1'})
    add('r2_ens2_all_types', {'RESTART': '2', 'NUM_LDAS_ENSEMBLE': '2', 'TILE_TYPES': '100 19 20',
                              'DO_ISSM': '1', 'RUN_ROUTE': '1'})
    add('r1_lake_only_zoom', {'TILE_TYPES': '19', 'MINLAT': '20.'})
    add('rM_lake', {'RESTART': 'M', 'TILE_TYPES': '100 19'})
    # coupled land-atm DAS: exeinp/batinp are dummy (non-existing) files
    for lc, nens in [('1', '1'), ('2', '2')]:
        add('ladas_cpl' + lc, cpl=['--ladas_cpl', lc, '--account', 's1234', '--nymdb', '20200101',
                                   '--nhmsb', '210000', '--agcm_res', '0012', '--bcs_version', BCS_VERSION,
                                   '--rstloc', D['rstloc'], '--varwindow', '360', '--nens', nens])
    add('B2_mwrtm_missing', {'MWRTM_PATH': os.path.join(D['root'], 'nomwrtm')}, expect='new_fails')
    # expect: 'same'      -> ref and new must succeed and be identical
    #         'diff'      -> ref and new must succeed; differences are printed (intended change)
    #         'ref_fails' -> ref setup fails (bug), new must succeed
    #         'new_fails' -> new setup must fail (input error detected); ref result ignored
    # Phase 2 bug fixes:
    add('B1_no_do_issm', {'TILE_TYPES': '100 20'}, drop=['DO_ISSM'])
    add('B6_nas_multi_sgmt', {'END_DATE': '20200501 000000', 'JOB_SGMT': '00000100 000000'},
        site='NAS')
    add('nas_single', site='NAS')
    add('B3_route_merra2', {'RESTART': 'M', 'RUN_ROUTE': '1'}, expect='new_fails')
    add('B5_bad_lsm_choice', {'LSM_CHOICE': '7'}, expect='new_fails')
    add('unknown_site', site='OTHER', expect='new_fails')            # phase 6: early GEOS_SITE check
    return C


def write_kv(path, d):
    with open(path, 'w') as f:
        for k, v in d.items():
            f.write('%s: %s\n' % (k, v))


def _last_error(exphome):
    if not os.path.exists(os.path.join(exphome, 'setup.log')):
        exphome = os.path.join(exphome, 'ldas')
    log = open(os.path.join(exphome, 'setup.log')).read().strip().splitlines()
    errs = [l for l in log if 'ERROR' in l or 'Error' in l]
    return (errs or log or [''])[-1]


def run_case(bindir, exphome, case, verbose, quiet=False):
    os.makedirs(exphome)
    exeinp = os.path.join(exphome, 'exeinp.txt')
    batinp = os.path.join(exphome, 'batinp.txt')
    if case['cpl']:
        # ldas exphome is a sub-dir of the ADAS exp dir; exeinp/batinp are dummy
        exphome = os.path.join(exphome, 'ldas')
        os.makedirs(exphome)
        cmd = [sys.executable, os.path.join(bindir, 'ldas_setup'), 'setup', exphome,
               'dummy_exeinp', 'dummy_batinp'] + case['cpl']
    else:
        write_kv(exeinp, case['exeinp'])
        write_kv(batinp, case['batinp'])
        cmd = [sys.executable, os.path.join(bindir, 'ldas_setup'), 'setup', exphome, exeinp, batinp]
    p = sp.run(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, cwd=exphome)
    log = p.stdout.decode(errors='replace')
    with open(os.path.join(exphome, 'setup.log'), 'w') as f:
        f.write(log)
    if verbose or (p.returncode != 0 and not quiet):
        print(textwrap.indent(log[-3000:], '      | '))
    return p.returncode


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', default=REF_COMMIT)
    ap.add_argument('--keep', action='store_true', help='keep work dir')
    ap.add_argument('--case', nargs='*')
    ap.add_argument('-v', '--verbose', action='store_true')
    o = ap.parse_args(argv)

    work = tempfile.mkdtemp(prefix='ldas_regr_')
    ok = True
    try:
        installs = {}
        def inst(which, site):
            key = (which, site)
            if key not in installs:
                installs[key] = make_install(os.path.join(work, 'inst_%s_%s' % key),
                                             o.ref if which == 'ref' else None, site, D['root'])
            return installs[key]
        D = make_data(os.path.join(work, 'data'))
        C = cases(D)
        def expdir(which, name):
            if C[name]['cpl']:
                return os.path.join(work, which, name, 'ldas', name + '_LDAS' + ('4ens' if C[name]['cpl'][1] == '2' else ''))
            return os.path.join(work, which, name, 'testexp')
        names = o.case or list(C)
        if not o.case:
            # 'ldas_setup sample' output must be identical
            for opt in ['--exeinp', '--batinp']:
                outs = [sp.run([sys.executable, os.path.join(inst(w, 'NCCS'), 'ldas_setup'), 'sample', opt],
                               stdout=sp.PIPE, stderr=sp.STDOUT).stdout for w in ['ref', 'new']]
                same = outs[0] == outs[1] and len(outs[0]) > 1000 if opt == '--exeinp' else outs[0] == outs[1]
                print('== %-24s %s' % ('sample ' + opt, 'IDENTICAL' if same else 'DIFFERS'))
                ok = ok and same
        for name in names:
            case = C[name]
            expect = case['expect']
            print('== %-24s' % name, end=' ', flush=True)
            quiet = expect in ('ref_fails', 'new_fails') and not o.verbose
            rc_ref = run_case(inst('ref', case['site']), os.path.join(work, 'ref', name), case, o.verbose, quiet)
            rc_new = run_case(inst('new', case['site']), os.path.join(work, 'new', name), case, o.verbose, quiet)
            if expect == 'new_fails':
                good = rc_new != 0
                print('new failed as expected (ref rc=%d)' % rc_ref if good else 'ERROR: new setup did not fail')
                if good:
                    print(textwrap.indent(_last_error(os.path.join(work, 'new', name)), '      '))
                ok = ok and good
                continue
            if expect == 'ref_fails':
                good = rc_ref != 0 and rc_new == 0
                print('ref failed, new OK (as expected)' if good else
                      'UNEXPECTED: ref rc=%d new rc=%d' % (rc_ref, rc_new))
                ok = ok and good
                continue
            if rc_ref != 0 or rc_new != 0:
                print('setup failed: ref rc=%d new rc=%d' % (rc_ref, rc_new))
                ok = False
                continue
            import io
            buf = io.StringIO()
            nd = compare(expdir('ref', name), expdir('new', name),
                         same=[(os.path.dirname(inst('ref', case['site'])),
                                os.path.dirname(inst('new', case['site']))),
                               (os.path.join(work, 'ref', name), os.path.join(work, 'new', name))],
                         out=buf)
            if nd == 0:
                print('IDENTICAL' + ('  (WARNING: a difference was expected)' if expect == 'diff' else ''))
            else:
                print('%d DIFFERENCE(S)%s' % (nd, '  (expected)' if expect == 'diff' else ''))
                print(textwrap.indent(buf.getvalue(), '      '))
                ok = ok and expect == 'diff'
    finally:
        if o.keep:
            print('work dir kept: ' + work)
        else:
            shutil.rmtree(work, ignore_errors=True)
    print('\nRESULT: ' + ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
