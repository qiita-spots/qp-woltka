# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main, TestCase
from os import environ
from os.path import join
from tempfile import TemporaryDirectory
import gzip
import io
import pandas as pd

from qp_woltka.util import (
    get_dbs, generate_woltka_dflt_params, mux, demux, search_by_filename,
    merge_ranges, coverage_percentage)


class UtilTests(TestCase):
    def setUp(self):
        self.db_path = environ["QC_WOLTKA_DB_DP"]

    def test_get_dbs(self):
        db_path = self.db_path
        obs = get_dbs(db_path)
        exp = {'wol': join(db_path, 'wol/WoLmin'),
               'rep82': join(db_path, 'rep82/5min')}

        self.assertDictEqual(obs, exp)

    def test_generate_woltka_dflt_params(self):
        obs = generate_woltka_dflt_params()
        exp = {'wol': {'Database': join(self.db_path, 'wol/WoLmin')},
               'rep82': {'Database': join(self.db_path, 'rep82/5min')}}

        self.assertDictEqual(obs, exp)

    def test_mux(self):
        f1 = b"@foo\nATGC\n+\nIIII\n"
        f2 = b"@bar\nAAAA\n+\nIIII\n"
        exp = b"@foo@@@foofile\nATGC\n+\nIIII\n@bar@@@barfile\nAAAA\n+\nIIII\n"
        with TemporaryDirectory() as d:
            f1fp = d + '/foofile.fastq'
            f2fp = d + '/barfile.fastq'
            ofp = d + '/output'
            with gzip.open(f1fp, 'wb') as fp:
                fp.write(f1)
            with gzip.open(f2fp, 'wb') as fp:
                fp.write(f2)
            with open(ofp, 'wb') as output:
                mux([f1fp, f2fp], output)
            with open(ofp, 'rb') as result:
                obs = result.read()

        self.assertEqual(obs, exp)

    def test_demux(self):
        prep = pd.DataFrame([['sample_foo', 'foofile'],
                             ['sample_bar', 'barfile']],
                            columns=['sample_name', 'run_prefix'])
        input_ = io.BytesIO(b"foo@@@foofile_R1\tATGC\t+\tIIII\nbar@@@"
                            b"barfile_R2\tAAAA\t+\tIIII\n")
        expfoo = b"foo\tATGC\t+\tIIII\n"
        expbar = b"bar\tAAAA\t+\tIIII\n"
        with TemporaryDirectory() as d:
            demux(input_, d.encode('ascii'), prep)
            foo = open(d + '/sample_foo.sam', 'rb').read()
            bar = open(d + '/sample_bar.sam', 'rb').read()

        self.assertEqual(foo, expfoo)
        self.assertEqual(bar, expbar)

    def test_search_by_filename(self):
        lookup = {'foo_bar': 'baz',
                  'foo': 'bar'}
        tests = [('foo_bar_thing', 'baz'),
                 ('foo_stuff_blah', 'bar'),
                 ('foo.stuff.blah', 'bar'),
                 ('foo_bar', 'baz')]
        for test, exp in tests:
            obs = search_by_filename(test, lookup)
            self.assertEqual(obs, exp)

        with self.assertRaises(KeyError):
            search_by_filename('does_not_exist', lookup)

    def test_merge_ranges(self):
        files = [
            'qp_woltka/support_files/coverages/coverage_1.cov',
            'qp_woltka/support_files/coverages/coverage_2.cov',
        ]
        exp = ['GXXX\t100\t400\t500\t600', 'GYYY\t100\t400\t500\t1000']
        self.assertEqual(merge_ranges(files), exp)

        files = [
            'qp_woltka/support_files/coverages/coverage_1.cov',
            'qp_woltka/support_files/coverages/coverage_2.cov',
            'qp_woltka/support_files/coverages/coverage_3.cov',
        ]
        exp = ['GXXX\t100\t400\t500\t600', 'GYYY\t100\t400\t500\t1000',
               'GZZZ\t200\t400\t500\t1500']
        self.assertEqual(merge_ranges(files), exp)

    def test_coverage_percentage(self):
        # testing with 3 cov
        files = [
            'qp_woltka/support_files/coverages/coverage_1.cov',
            'qp_woltka/support_files/coverages/coverage_2.cov',
            'qp_woltka/support_files/coverages/coverage_3.cov',
        ]
        length_map = 'qp_woltka/support_files/coverages/length_map.map'
        exp = ['GXXX\t67.00', 'GYYY\t72.91', 'GZZZ\t60.10']
        self.assertEqual(coverage_percentage(files, length_map), exp)

        # testing with the archive.cov | note same exp values
        self.assertEqual(coverage_percentage(
            ['qp_woltka/support_files/coverages/artifact.cov'], length_map),
            exp)

        # testing errors
        length_map = 'qp_woltka/support_files/coverages/length_map_bad.map'
        with self.assertRaises(ValueError):
            coverage_percentage(files, length_map)


if __name__ == '__main__':
    main()
