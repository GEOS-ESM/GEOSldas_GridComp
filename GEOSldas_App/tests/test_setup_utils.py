#!/usr/bin/env python3
"""
Unit tests for setup_utils.py (pure helpers used by ldas_setup / ldas.py).

Run from GEOSldas_App:
    python3 -m unittest discover -s tests -v
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import setup_utils as su


def _write(dirname, name, text):
    path = os.path.join(dirname, name)
    with open(path, 'w') as f:
        f.write(text)
    return path


class TestHelpers(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    # require / glob_first ------------------------------------------------
    def test_require(self):
        su.require(True, 'never raised')
        with self.assertRaises(su.LdasSetupError) as cm:
            su.require(False, 'boom')
        self.assertEqual(str(cm.exception), 'boom')

    def test_glob_first(self):
        _write(self.d, 'green_clim_a.data', '')
        self.assertTrue(su.glob_first(self.d + '/green_clim_*.data', 'green').endswith('green_clim_a.data'))
        with self.assertRaises(su.LdasSetupError):
            su.glob_first(self.d + '/lai_clim_*.data', 'lai')

    # bcs_version_from_path -----------------------------------------------
    def test_bcs_version_from_path(self):
        self.assertEqual(su.bcs_version_from_path('/a/b/v12'),    'v12')
        self.assertEqual(su.bcs_version_from_path('/a/b/v12///'), 'v12')

    # hours_to_hhmmss -------------------------------------------------------
    def test_hours_to_hhmmss(self):
        self.assertEqual(su.hours_to_hhmmss(6),    '060000')
        self.assertEqual(su.hours_to_hhmmss(1.5),  '013000')
        self.assertEqual(su.hours_to_hhmmss(0),    '000000')
        self.assertEqual(su.hours_to_hhmmss(360/60), '060000')

    # replace_in_file -------------------------------------------------------
    def test_replace_in_file(self):
        f = _write(self.d, 'x.rc', 'EXP: GEOSldas_expid   \nother\n')
        su.replace_in_file(f, 'GEOSldas_expid', 'myexp')
        with open(f) as fin:
            self.assertEqual(fin.read(), 'EXP: myexp\nother\n')

    # read_gridname -------------------------------------------------------
    def test_read_gridname_til(self):
        f = _write(self.d, 'a.til', 'line1\nline2\n  SMAP-EASEv2-M36  \nline4\n')
        self.assertEqual(su.read_gridname(f), 'EASEv2_M36')

    def test_read_gridname_til_domain(self):
        f = _write(self.d, 'a.til.domain', 'l1\nl2\nPE720x4320-CF\n')
        self.assertEqual(su.read_gridname(f), 'PE720x4320-CF')

    # date helpers ----------------------------------------------------------
    def test_date_helpers(self):
        from datetime import datetime
        dt = datetime(2020, 3, 5, 21, 30)
        self.assertEqual(su.ym_dir(dt), 'Y2020/M03')
        self.assertEqual(su.ymd_hm(dt), '20200305_2130')

    # copy_with_suffix ------------------------------------------------------
    def test_copy_with_suffix(self):
        src = _write(self.d, 'exeinp', 'a\n')
        dst = os.path.join(self.d, 'run'); os.makedirs(dst)
        self.assertTrue(su.copy_with_suffix(src, dst).endswith('/run/exeinp'))
        self.assertTrue(su.copy_with_suffix(src, dst).endswith('/run/exeinp.1'))
        self.assertTrue(su.copy_with_suffix(src, dst).endswith('/run/exeinp.2'))

    # run_cmd ---------------------------------------------------------------
    def test_run_cmd(self):
        self.assertEqual(su.run_cmd('true'), 0)
        self.assertEqual(su.run_cmd('false'), 1)            # warning only
        with self.assertRaises(su.LdasSetupError):
            su.run_cmd('false', check=True)

    # resource keys -----------------------------------------------------
    def test_resource_keys(self):
        self.assertIn('account', su.getResourceKeys('required'))
        self.assertIn('qos',     su.getResourceKeys('optional'))
        with self.assertRaises(su.LdasSetupError):
            su.getResourceKeys('bogus')
        with self.assertRaises(su.LdasSetupError):
            su.getExeKeys('9')


class TestParseInputFile(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_basic(self):
        f = _write(self.d, 'exeinp', (
            '# comment\n'
            '\n'
            'EXP_ID:  myexp   # trailing comment\n'
            'BEG_DATE: 20200101 000000\n'
            'MET_PATH: /a:b/c\n'          # only first ':' splits
            'EMPTY:\n'                    # warning, skipped
            ))
        d = su.parseInputFile(f)
        self.assertEqual(list(d.keys()), ['EXP_ID', 'BEG_DATE', 'MET_PATH'])
        self.assertEqual(d['EXP_ID'],   'myexp')
        self.assertEqual(d['BEG_DATE'], '20200101 000000')
        self.assertEqual(d['MET_PATH'], '/a:b/c')

    def test_rc_defaults_prefix(self):
        f = _write(self.d, 'x.rc', (
            '#GEOSldas=>NUM_LDAS_ENSEMBLE: 1\n'
            '#GEOSagcm=>NUM_LDAS_ENSEMBLE: 32\n'
            ))
        self.assertEqual(su.parseInputFile(f)['NUM_LDAS_ENSEMBLE'], '1')
        self.assertEqual(su.parseInputFile(f, ladas_cpl=1)['NUM_LDAS_ENSEMBLE'], '32')

    def test_duplicate_key(self):
        f = _write(self.d, 'dup', 'A: 1\nA: 2\n')
        with self.assertRaises(su.LdasSetupError):
            su.parseInputFile(f)

    def test_bad_line(self):
        f = _write(self.d, 'bad', 'no colon here\n')
        with self.assertRaises(su.LdasSetupError):
            su.parseInputFile(f)

    def test_missing_file(self):
        missing = os.path.join(self.d, 'nope')
        self.assertEqual(su.parseInputFile(missing, ladas_cpl=1), {})
        with self.assertRaises(su.LdasSetupError):
            su.parseInputFile(missing)


if __name__ == '__main__':
    unittest.main()
