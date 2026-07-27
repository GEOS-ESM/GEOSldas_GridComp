
# Self cross-masked ObsFcstAna sums: matches two ASCAT product species (e.g. legacy
# vs H121, same platform) *within a single experiment's own files*, instead of
# cross-masking one experiment's obs against a second experiment's (as
# postproc_ObsFcstAna.compute_monthly_sums does for OL-vs-DA obs substitution).
#
# For each cycle and each declared species pair (species_a, species_b), a tile only
# contributes to species_a's sums if species_b also has a valid obs on that same
# tile in that same cycle, and vice versa. All other var_list fields (obs_fcst,
# obs_obsvar, obs_fcstvar, obs_ana, obs_anavar) are taken from the *same* file/run
# for both species -- there is no cross-experiment substitution here, so O_mean and
# F_mean for the two paired species are always drawn from a common, jointly-matched
# sample and can be differenced directly (e.g. H121 ASCAT minus legacy ASCAT).
#
# Output sums/stats files are written in the exact same (tile, species) shape as
# postproc_ObsFcstAna, via the same write_sums_nc4/write_stats_nc4 helpers, so the
# inherited calc_temporal_stats_from_sums/calc_spatial_stats_from_sums methods work
# unmodified.
#
# amfox - July 2026

import os
import numpy as np

from datetime import timedelta
from dateutil.relativedelta import relativedelta

from read_GEOSldas import read_ObsFcstAna, read_ObsFcstAna_nc4
from postproc_ObsFcstAna import postproc_ObsFcstAna


class postproc_ObsFcstAna_SelfXMask(postproc_ObsFcstAna):

    def __init__(self, exp, species_pairs, start_time, end_time, sum_path='./'):
        # species_pairs: dict of {species_a: species_b} (obsparam 'species' values,
        # not array indices); each pair is cross-masked against the other within
        # this single experiment's own ObsFcstAna files.
        super().__init__([exp], start_time, end_time, sum_path=sum_path)
        self.species_pairs = species_pairs

    # ------------------------------------------------------------------------------
    #
    # Override: cross-mask matched species pairs within one experiment's own files,
    # instead of cross-masking against a second experiment's (scaled) obs.

    def compute_monthly_sums(self, date_time):

        expdir = self.expdir_list[0]
        expid = self.expid_list[0]
        tc = self.tilecoord
        obs_param = self.obsparam_list[0]
        var_list = self.var_list
        da_dt = self.da_dt

        n_tile = tc['N_tile']
        n_spec = len(obs_param)
        species_to_ispec = {int(obs_param[i]['species']): i for i in range(n_spec)}

        date_time = date_time.replace(hour=int(self.da_t0), minute=int(np.mod(self.da_t0, 1) * 60))
        month_string = date_time.strftime('%Y%m')
        stop_time = date_time + relativedelta(months=1)

        data_sum = {}
        data2_sum = {}
        N_data = np.zeros((n_tile, n_spec))
        oxf_sum = np.zeros((n_tile, n_spec))
        oxa_sum = np.zeros((n_tile, n_spec))
        fxa_sum = np.zeros((n_tile, n_spec))
        for var in var_list:
            data_sum[var] = np.zeros((n_tile, n_spec))
            data2_sum[var] = np.zeros((n_tile, n_spec))

        n_files_read = 0

        while date_time < stop_time:

            fname_base = (expdir + expid + '/output/' + self.domain + '/ana/ens_avg/Y' +
                          date_time.strftime('%Y') + '/M' + date_time.strftime('%m') + '/' +
                          expid + '.ens_avg.ldas_ObsFcstAna.' + date_time.strftime('%Y%m%d_%H%M') + 'z')
            fname_nc4 = fname_base + '.nc4'
            fname_bin = fname_base + '.bin'

            OFA = None
            if os.path.isfile(fname_nc4):
                OFA = read_ObsFcstAna_nc4(fname_nc4)
                n_files_read += 1
            elif os.path.isfile(fname_bin):
                OFA = read_ObsFcstAna(fname_bin)
                n_files_read += 1

            if OFA is not None and len(OFA['obs_tilenum']) > 0:

                data_tile = {var: np.full((n_tile, n_spec), np.nan) for var in var_list}

                for ispec in np.arange(n_spec):
                    this_species = int(obs_param[ispec]['species'])
                    if obs_param[ispec]['assim'] == 'T':
                        is_valid = np.logical_and(OFA['obs_species'] == this_species, OFA['obs_assim'] == 1)
                    else:
                        is_valid = OFA['obs_species'] == this_species

                    tile_idx = OFA['obs_tilenum'][is_valid] - 1
                    for var in var_list:
                        data_tile[var][tile_idx, ispec] = OFA[var][is_valid]

                # cross-mask each declared species pair against each other, within
                # this single file/cycle (tile-for-tile; no cross-experiment obs
                # substitution)
                has_obs = ~np.isnan(data_tile['obs_obs'])
                for sp_a, sp_b in self.species_pairs.items():
                    ia = species_to_ispec[sp_a]
                    ib = species_to_ispec[sp_b]
                    cross_valid = has_obs[:, ia] & has_obs[:, ib]
                    drop_a = has_obs[:, ia] & ~cross_valid
                    drop_b = has_obs[:, ib] & ~cross_valid
                    for var in var_list:
                        data_tile[var][drop_a, ia] = np.nan
                        data_tile[var][drop_b, ib] = np.nan

                is_valid = ~np.isnan(data_tile['obs_obs'])
                N_data[is_valid] += 1
                oxf_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_fcst'][is_valid]
                oxa_sum[is_valid] += data_tile['obs_obs'][is_valid] * data_tile['obs_ana'][is_valid]
                fxa_sum[is_valid] += data_tile['obs_fcst'][is_valid] * data_tile['obs_ana'][is_valid]

                for var in var_list:
                    data_sum[var][is_valid] += data_tile[var][is_valid]
                    data2_sum[var][is_valid] += data_tile[var][is_valid] ** 2

            date_time = date_time + timedelta(seconds=da_dt)

        if n_files_read == 0:
            raise FileNotFoundError('No ObsFcstAna .nc4 or .bin files found for ' + month_string)

        return N_data, data_sum, data2_sum, oxf_sum, oxa_sum, fxa_sum

# ============== EOF ====================================================================
