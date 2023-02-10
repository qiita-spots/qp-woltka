# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo
from os import remove, environ
from os.path import exists, isdir, join, dirname, basename
from shutil import rmtree, copyfile
from tempfile import mkdtemp, TemporaryDirectory
from json import dumps
import gzip
import io
import pandas as pd

from qp_woltka import plugin
from qp_woltka.util import (
    get_dbs, generate_woltka_dflt_params, mux, demux, search_by_filename)
from qp_woltka.woltka import woltka_to_array, woltka


class WoltkaTests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        out_dir = mkdtemp()
        self.maxDiff = None
        self.out_dir = out_dir
        self.db_path = environ["QC_WOLTKA_DB_DP"]
        self.params = {
            'Database': join(self.db_path, 'rep82/5min'),
        }
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)
        self.environment = environ["ENVIRONMENT"]

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

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

    def _helper_woltka_bowtie(self, prep_info_dict, database=None):
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'S22205_S104_L001_R1_001.fastq.gz')
        fp1_2 = join(in_dir, 'S22205_S104_L001_R2_001.fastq.gz')
        fp2_1 = join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz')
        fp2_2 = join(in_dir, 'S22282_S102_L001_R2_001.fastq.gz')
        source_dir = 'qp_woltka/support_files'
        copyfile(f'{source_dir}/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
        copyfile(f'{source_dir}/S22205_S104_L001_R2_001.fastq.gz', fp1_2)
        copyfile(f'{source_dir}/S22282_S102_L001_R1_001.fastq.gz', fp2_1)
        copyfile(f'{source_dir}/S22282_S102_L001_R2_001.fastq.gz', fp2_2)

        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs'),
                (fp2_1, 'raw_forward_seqs'),
                (fp2_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test Woltka artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid

        if database is not None:
            self.params['Database'] = database

        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-woltka', '2022.09', 'Woltka v0.1.4']),
                'status': 'running',
                'parameters': dumps(self.params)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        return pid, aid, job_id

    def test_woltka_to_array_rep82(self):
        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        files, prep = self.qclient.artifact_and_preparation_files(aid)

        url = 'this-is-my-url'
        database = self.params['Database']
        main_fp, merge_fp = woltka_to_array(
            files, out_dir, database, prep, url, job_id)

        self.assertEqual(join(out_dir, f'{job_id}.slurm'), main_fp)
        self.assertEqual(join(out_dir, f'{job_id}.merge.slurm'), merge_fp)

        with open(main_fp) as f:
            main = f.readlines()
        with open(merge_fp) as f:
            merge = f.readlines()

        prep_file = join(out_dir, 'prep_info.tsv')
        exp_main = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 40:00:00\n',
            '#SBATCH --mem 90g\n',
            f'#SBATCH --output {out_dir}/{job_id}_%a.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_%a.err\n',
            '#SBATCH --array 1-1%8\n',
            f'cd {out_dir}\n',
            f'prep_full_path={prep_file}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            f'dbbase={dirname(database)}\n',
            f'dbname={basename(database)}\n',
            f'output={out_dir}\n',
            'files=`cat sample_details_${SLURM_ARRAY_TASK_ID}.txt`\n',
            'mux ${files} | bowtie2 -p 8 -x '
            f'{database} -q - --seed 42 '
            '--very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal | demux ${output} '
            f'{prep_file} | sort | uniq > '
            'sample_processing_${SLURM_ARRAY_TASK_ID}.log\n',
            '# for each one of our input files, form woltka commands, \n',
            '# and farm off to gnu parallel\n',
            'for f in `cat sample_processing_${SLURM_ARRAY_TASK_ID}.log`\n',
            'do\n',
            '  echo "woltka classify -i ${f} -o ${f}.woltka-taxa --no-demux '
            f'--lineage {database}.tax --rank free,none"\n',
            'done | parallel -j 8\n',
            'for f in `cat sample_processing_${SLURM_ARRAY_TASK_ID}.log`\n',
            'do\n',
            '  # compress it\n',
            '  echo "xz -1 -T1 ${f}"\n',
            'done | parallel -j 8\n',
            'date\n']
        self.assertEqual(main, exp_main)

        exp_merge = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name merge-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 2\n',
            '#SBATCH --time 30:00:00\n',
            '#SBATCH --mem 140g\n',
            f'#SBATCH --output {out_dir}/merge-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/merge-{job_id}.err\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            'set -e\n',
            "sruns=`grep 'overall alignment rate' *.err | wc -l`\n",
            "sjobs=`ls sample_details_* | wc -l`\n",
            'if [[ ! -f "errors.log" && $sruns -eq $sjobs ]]; then\n',
            f'woltka_merge --base {out_dir}  --name '
            'free --glob "*.woltka-taxa/free.biom" &\n',
            f'woltka_merge --base {out_dir}  --name '
            'none --glob "*.woltka-taxa/none.biom" --rename &\n',
            'wait\n',
            '\n',
            f'cd {out_dir}; tar -cvf alignment.tar *.sam.xz\n',
            'fi\n',
            f'finish_woltka {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(merge, exp_merge)

        # now let's test that if finished correctly
        sdir = 'qp_woltka/support_files/'
        copyfile(f'{sdir}/none.biom', f'{out_dir}/none.biom')
        copyfile(f'{sdir}/free.biom', f'{out_dir}/free.biom')
        copyfile(f'{sdir}/alignment.tar', f'{out_dir}/alignment.tar')

        success, ainfo, msg = woltka(
            self.qclient, job_id, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        exp = [
            ArtifactInfo('Alignment Profile', 'BIOM',
                         [(f'{out_dir}/free.biom', 'biom'),
                          (f'{out_dir}/alignment.tar', 'log')]),
            ArtifactInfo('Per genome Predictions', 'BIOM',
                         [(f'{out_dir}/none.biom', 'biom')])]

        self.assertCountEqual(ainfo, exp)

    def test_woltka_to_array_wol(self):
        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104_L001_R'},
            'SKD8.640184': {'run_prefix': 'S22282_S102_L001_R'}}
        database = join(self.db_path, 'wol/WoLmin')

        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict, database)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        files, prep = self.qclient.artifact_and_preparation_files(aid)

        url = 'this-is-my-url'
        main_fp, merge_fp = woltka_to_array(
             files, out_dir, database, prep, url, job_id)

        self.assertEqual(join(out_dir, f'{job_id}.slurm'), main_fp)
        self.assertEqual(join(out_dir, f'{job_id}.merge.slurm'), merge_fp)

        with open(main_fp) as f:
            main = f.readlines()
        with open(merge_fp) as f:
            merge = f.readlines()

        prep_file = join(out_dir, 'prep_info.tsv')
        exp_main = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 40:00:00\n',
            '#SBATCH --mem 90g\n',
            f'#SBATCH --output {out_dir}/{job_id}_%a.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_%a.err\n',
            '#SBATCH --array 1-1%8\n',
            f'cd {out_dir}\n',
            f'prep_full_path={prep_file}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            f'dbbase={dirname(database)}\n',
            f'dbname={basename(database)}\n',
            f'output={out_dir}\n',
            'files=`cat sample_details_${SLURM_ARRAY_TASK_ID}.txt`\n',
            'mux ${files} | bowtie2 -p 8 -x '
            f'{database} -q - --seed 42 '
            '--very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal | demux ${output} '
            f'{prep_file} | sort | uniq > '
            'sample_processing_${SLURM_ARRAY_TASK_ID}.log\n',
            '# for each one of our input files, form woltka commands, \n',
            '# and farm off to gnu parallel\n',
            'for f in `cat sample_processing_${SLURM_ARRAY_TASK_ID}.log`\n',
            'do\n',
            '  echo "woltka classify -i ${f} -o ${f}.woltka-taxa --no-demux '
            f'--lineage {database}.tax --rank free,none"\n',
            f'  echo "woltka classify -i ${{f}} -c {database}.coords '
            '-o ${f}.woltka-per-gene --no-demux"\n',
            'done | parallel -j 8\n',
            'for f in `cat sample_processing_${SLURM_ARRAY_TASK_ID}.log`\n',
            'do\n',
            '  # compress it\n',
            '  echo "xz -1 -T1 ${f}"\n',
            'done | parallel -j 8\n',
            'date\n']
        self.assertEqual(main, exp_main)

        exp_merge = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name merge-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 3\n',
            '#SBATCH --time 30:00:00\n',
            '#SBATCH --mem 140g\n',
            f'#SBATCH --output {out_dir}/merge-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/merge-{job_id}.err\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            'set -e\n',
            "sruns=`grep 'overall alignment rate' *.err | wc -l`\n",
            "sjobs=`ls sample_details_* | wc -l`\n",
            'if [[ ! -f "errors.log" && $sruns -eq $sjobs ]]; then\n',
            f'woltka_merge --base {out_dir}  --name '
            'free --glob "*.woltka-taxa/free.biom" &\n',
            f'woltka_merge --base {out_dir}  --name '
            'none --glob "*.woltka-taxa/none.biom" &\n',
            f'woltka_merge --base {out_dir}  --name '
            'per-gene --glob "*.woltka-per-gene" --rename &\n',
            'wait\n',
            '\n',
            f'cd {out_dir}; tar -cvf alignment.tar *.sam.xz\n',
            'fi\n',
            f'finish_woltka {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(merge, exp_merge)

        # now let's test that if finished correctly
        sdir = 'qp_woltka/support_files/'
        copyfile(f'{sdir}/none.biom', f'{out_dir}/none.biom')
        copyfile(f'{sdir}/per-gene.biom', f'{out_dir}/per-gene.biom')
        copyfile(f'{sdir}/free.biom', f'{out_dir}/free.biom')
        copyfile(f'{sdir}/alignment.tar', f'{out_dir}/alignment.tar')
        success, ainfo, msg = woltka(
            self.qclient, job_id, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        exp = [
            ArtifactInfo('Alignment Profile', 'BIOM',
                         [(f'{out_dir}/free.biom', 'biom'),
                          (f'{out_dir}/alignment.tar', 'log')]),
            ArtifactInfo('Per genome Predictions', 'BIOM',
                         [(f'{out_dir}/none.biom', 'biom')]),
            ArtifactInfo('Per gene Predictions', 'BIOM',
                         [(f'{out_dir}/per-gene.biom', 'biom')])]

        self.assertCountEqual(ainfo, exp)

    def test_woltka_to_array_error(self):
        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
        database = join(self.db_path, 'wol/WoLmin')
        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict, database)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = woltka(
            self.qclient, job_id, self.params, out_dir)

        exp_msg = '\n'.join([
            'Missing files from the "Alignment Profile"; please contact '
            'qiita.help@gmail.com for more information',
            'Table none/per-genome was not created, please contact '
            'qiita.help@gmail.com for more information',
            'Table per-gene was not created, please contact '
            'qiita.help@gmail.com for more information'
        ])
        self.assertEqual(exp_msg, msg)
        self.assertFalse(success)

    def test_creation_error_no_run_prefix(self):
        # run prefix doesn't exist
        prep_info_dict = {
            'SKB8.640193': {'run_prefix_bla': 'S22205_S104'},
            'SKD8.640184': {'run_prefix_bla': 'S22282_S102'}}
        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        with self.assertRaises(KeyError) as error:
            self.qclient.artifact_and_preparation_files(aid)
        self.assertEqual(str(error.exception), "'run_prefix'")

    def test_creation_error_no_unique_run_prefix(self):
        # run prefix doesn't exist
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22205_S104'}}
        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # retriving info of the prep/artifact just created
        with self.assertRaises(ValueError) as error:
            files, prep = self.qclient.artifact_and_preparation_files(aid)
        self.assertEqual(str(error.exception), "Multiple run prefixes match "
                         "this fwd read: S22205_S104_L001_R1_001.fastq.gz")

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


if __name__ == '__main__':
    main()
