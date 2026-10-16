#!/usr/bin/env python3

# Compute observation scaling parameters (climatological mean and stdv of obs and
# model forecasts) from the ObsFcstAna files of a GEOSldas experiment.
#
#   --method latlon : soil moisture obs (e.g., ASCAT, CYGNSS) on a global lat/lon grid,
#                     nc4 output for GEOSldas scale_obs_sfmc_zscore
#   --method tb     : Tb obs (e.g., SMAP, SMOS) in tile space, one orbit per run,
#                     binary output for GEOSldas scale_obs_Tb_zscore
#
# Output goes to <exp-path>/<exp-run>/output/<domain>/stats/<out-dir>/.
#
# Examples:
#
#   python run_obs_scaling_params.py --method latlon --exp-path /path/to/exps --exp-run OL1 \
#       --domain CF0360x6C_GLOBAL --start 2019-02 --end 2020-08 \
#       --species ASCAT_HSAF_META_SM,ASCAT_HSAF_METB_SM,ASCAT_HSAF_METC_SM \
#       --prefix CF0360_OL1_HSAF_zscore_stats_
#
#   python run_obs_scaling_params.py --method tb --exp-path /path/to/exps --exp-run OL1 \
#       --domain SMAP_EASEv2_M36_GLOBAL --start 2015-04 --end 2021-03 \
#       --descr-contains SMAP_L1C --orbit D --prefix OL1_SMAP_zscore_stats_
#
# amfox - Oct 2026

import argparse
import csv
import os
import sys;        sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../shared/python'))

import numpy as np

from datetime import datetime

from read_GEOSldas           import read_obs_param
from obs_scaling_utils      import get_year_bounds
from get_obs_scaling_latlon import get_obs_scaling_latlon
from get_obs_scaling_tb     import get_obs_scaling_tb

def parse_args(argv=None):

    p = argparse.ArgumentParser(description='Compute observation scaling parameters from GEOSldas ObsFcstAna files.')

    p.add_argument('--method',      required=True, choices=['latlon', 'tb'])
    p.add_argument('--exp-path',    required=True, help='directory that contains the experiment')
    p.add_argument('--exp-run',     required=True, help='experiment ID')
    p.add_argument('--domain',      required=True, help='e.g. SMAP_EASEv2_M36_GLOBAL')
    p.add_argument('--start',       required=True, help='first month of the period, YYYY-MM')
    p.add_argument('--end',         required=True, help='last month of the period, YYYY-MM')
    p.add_argument('--prefix',      required=True, help='output file name prefix')
    p.add_argument('--out-dir',     default='z_score_clim', help='output subdirectory of output/<domain>/stats/')
    p.add_argument('--run-months',  default='',
                   help='comma-separated months to process (default: 1-12 plus the months needed to wrap the window)')
    p.add_argument('--window-days', type=int, default=75, help='length of moving window [days]')
    p.add_argument('--ndata-min',   type=int, default=20, help='minimum number of obs in window')
    p.add_argument('--dt-assim',    type=int, default=10800, help='ObsFcstAna file interval [s]')
    p.add_argument('--t0-assim',    type=int, default=0,     help='time of first ObsFcstAna file of the day [s]')
    p.add_argument('--obsparam',    default='',
                   help='ldas_obsparam file (default: rc_out/Y<start>/M<start>/<exp-run>.ldas_obsparam.<start>01_0000z.txt)')
    p.add_argument('--each-doy',    action='store_true', help='also write one file per DOY')

    g = p.add_argument_group('latlon method')
    g.add_argument('--species',              default='', help='comma-separated obs descr from ldas_obsparam')
    g.add_argument('--grid-resolution',      type=float, default=0.25, help='lat/lon grid spacing [deg]')
    g.add_argument('--no-combine-species',   dest='combine_species', action='store_false',
                   help='one set of statistics per species (default: all species combined)')
    g.add_argument('--each-pentad',          action='store_true', help='also write one file per pentad')
    g.add_argument('--no-all-pentads',       dest='all_pentads', action='store_false',
                   help='do not write the file with all pentads')
    g.add_argument('--exclude-tiles-csv',    default='',
                   help="CSV file with a 'tilenum' column and a column (see --exclude-tiles-column) that is "
                        "true for tiles whose obs are not used, e.g. to apply a QC that postdates the ObsFcstAna files")
    g.add_argument('--exclude-tiles-column', default='reject')

    g = p.add_argument_group('tb method')
    g.add_argument('--descr-contains',  default='SMAP_L1C', help='select Tb species whose descr contains this string')
    g.add_argument('--orbit',           choices=['A', 'D'], default='D')
    g.add_argument('--angles',          default='40', help='comma-separated incidence angles [deg]')
    g.add_argument('--no-convert-grid', dest='convert_grid', action='store_false',
                   help='use all model tiles (default: EASEv2 M36 administering tiles)')

    args = p.parse_args(argv)

    # check and convert the arguments

    for name in ('start', 'end'):
        try:
            setattr(args, name, datetime.strptime(getattr(args, name), '%Y-%m'))
        except ValueError:
            p.error(f'--{name} must be YYYY-MM, not {getattr(args, name)}')
    if args.start > args.end:
        p.error('--start must not be after --end')

    try:
        args.run_months = [int(x) for x in split_list(args.run_months)]
    except ValueError:
        p.error(f'--run-months must be comma-separated integers, not {args.run_months}')
    if any(m < 1 or m > 12 for m in args.run_months):
        p.error(f'--run-months must be between 1 and 12: {args.run_months}')

    args.species = split_list(args.species)
    if args.method == 'latlon' and not args.species:
        p.error('--species is required for --method latlon')

    try:
        args.angles = [float(x) for x in split_list(args.angles)]
    except ValueError:
        p.error(f'--angles must be comma-separated numbers, not {args.angles}')
    if args.method == 'tb' and not args.angles:
        p.error('--angles is required for --method tb')
    if len(set(args.angles)) != len(args.angles):
        p.error(f'duplicate --angles: {args.angles}')

    return args

def split_list(value):

    # comma-separated list -> list of stripped, non-empty strings

    return [x.strip() for x in value.split(',') if x.strip()]

def main(argv=None):

    args = parse_args(argv)

    start_year, start_month = args.start.year, args.start.month
    end_year,   end_month   = args.end.year,   args.end.month

    if args.run_months:
        run_months = args.run_months
    else:
        run_months = list(range(1, 13)) + list(range(1, int(np.ceil(args.window_days / 30)) + 1))

    first_year, last_year = get_year_bounds(run_months, start_year, start_month, end_year, end_month)

    obsparam_file = args.obsparam or (
        f'{args.exp_path}/{args.exp_run}/output/{args.domain}/rc_out/Y{start_year:04d}/M{start_month:02d}/'
        f'{args.exp_run}.ldas_obsparam.{start_year:04d}{start_month:02d}01_0000z.txt')

    obs_param = read_obs_param(obsparam_file)

    print('Observation scaling configuration:')
    for key, value in vars(args).items():
        if isinstance(value, datetime):
            value = value.strftime('%Y-%m')
        print(f'  {key:22s}: {value}')

    if args.method == 'latlon':

        # one species ID per name, in the order of the names (the output files of the
        # separate-species statistics are named with species_names)

        species_names = args.species
        species       = []
        for name in species_names:
            ids = [int(p['species']) for p in obs_param if p['descr'] == name]
            if len(ids) != 1:
                raise ValueError(f'expected one obsparam entry for species {name} in {obsparam_file}, found {len(ids)}')
            species += ids
        if len(set(species)) != len(species):
            raise ValueError(f'duplicate species in --species: {species_names}')
        print(f'  species ids             : {species}')

        exclude_tiles = None
        if args.exclude_tiles_csv:
            exclude_tiles = read_exclude_tiles(args.exclude_tiles_csv, args.exclude_tiles_column)

        get_obs_scaling_latlon(args.exp_path, args.exp_run, args.domain,
                               species           = species,
                               species_names     = species_names,
                               combine_species   = args.combine_species,
                               run_months        = run_months,
                               first_year        = first_year,
                               last_year         = last_year,
                               dt_assim          = args.dt_assim,
                               t0_assim          = args.t0_assim,
                               resol             = args.grid_resolution,
                               w_days            = args.window_days,
                               ndata_min         = args.ndata_min,
                               prefix            = args.prefix,
                               out_dir           = args.out_dir,
                               print_each_doy    = args.each_doy,
                               print_each_pentad = args.each_pentad,
                               print_all_pentads = args.all_pentads,
                               exclude_tiles     = exclude_tiles)

    else:

        get_obs_scaling_tb(args.exp_path, args.exp_run, args.domain,
                           obs_param      = obs_param,
                           descr_contains = args.descr_contains,
                           orbit          = 1 if args.orbit == 'A' else 2,
                           angles         = args.angles,
                           run_months     = run_months,
                           first_year     = first_year,
                           last_year      = last_year,
                           dt_assim       = args.dt_assim,
                           t0_assim       = args.t0_assim,
                           w_days         = args.window_days,
                           ndata_min      = args.ndata_min,
                           prefix         = args.prefix,
                           out_dir        = args.out_dir,
                           convert_grid   = args.convert_grid,
                           print_each_doy = args.each_doy)

# ----------------------------------------------------------------------------
#
# read CSV file with tilenum and exclude flag; returns the tile numbers (1-based) whose
# obs are not used, or None if no tile is flagged

def read_exclude_tiles(fname, column):

    tilenum = []
    reject  = []
    with open(fname, newline='') as fid:
        reader = csv.DictReader(fid)
        if column not in (reader.fieldnames or []):
            raise ValueError(f'column {column} not found in {fname} (columns: {reader.fieldnames})')
        for row in reader:
            tilenum.append(int(row['tilenum']))
            reject.append(str(row[column]).strip().lower() in ('true', '1', 'yes'))

    tilenum = np.array(tilenum, dtype=int)
    reject  = np.array(reject,  dtype=bool)

    print(f'  excluding obs at {reject.sum()} tiles listed in {fname}')

    if not np.any(reject):
        return None

    return np.unique(tilenum[reject])

if __name__ == '__main__':
    main()

# ====================== EOF =========================================================
