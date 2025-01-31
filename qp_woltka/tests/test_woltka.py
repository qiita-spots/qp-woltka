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
from os import remove, environ, mkdir
from os.path import exists, isdir, join, dirname, basename
from shutil import rmtree, copyfile
from tempfile import mkdtemp
from json import dumps

from qp_woltka import plugin
from qp_woltka.woltka import (
    woltka_to_array, woltka, woltka_syndna_to_array, woltka_syndna,
    calculate_cell_counts, calculate_rna_copy_counts)


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
        fp_summary = join(in_dir, 'summary.html')
        source_dir = 'qp_woltka/support_files'
        copyfile(f'{source_dir}/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
        copyfile(f'{source_dir}/S22205_S104_L001_R2_001.fastq.gz', fp1_2)
        copyfile(f'{source_dir}/S22282_S102_L001_R1_001.fastq.gz', fp2_1)
        copyfile(f'{source_dir}/S22282_S102_L001_R2_001.fastq.gz', fp2_2)
        copyfile(f'{source_dir}/summary.html', fp_summary)

        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs'),
                (fp2_1, 'raw_forward_seqs'),
                (fp2_2, 'raw_reverse_seqs'),
                (fp_summary, 'html_summary')]),
            'type': "per_sample_FASTQ",
            'name': "Test Woltka artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid

        if database is not None:
            self.params['Database'] = database

        data = {'user': 'demo@microbio.me',
                'command': dumps(
                    ['qp-woltka', '2024.09',
                     'Woltka v0.1.7, paired-end']),
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
        html_summary = self.qclient.get_artifact_html_summary(aid)
        files['html_summary'] = html_summary

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
            '#SBATCH --mail-type=FAIL,TIME_LIMIT_80,INVALID_DEPEND\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 80:00:00\n',
            '#SBATCH --mem 100g\n',
            f'#SBATCH --output {out_dir}/{job_id}_%a.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_%a.err\n',
            '#SBATCH --array 0-0%12\n',
            '#SBATCH --constraint="amd"\n',
            f'cd {out_dir}\n',
            f'prep_full_path={prep_file}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            f'dbbase={dirname(database)}\n',
            f'dbname={basename(database)}\n',
            f'output={out_dir}\n',
            'bt2_cores=6\n',
            f'mxdx mux --file-map {out_dir}/files_list.tsv --batch '
            '${SLURM_ARRAY_TASK_ID} --batch-size 50000000 '
            '--paired-handling interleave | '
            'bowtie2 -p ${bt2_cores} -x '
            f'{database} --interleaved - --seed 42 --very-sensitive -k 16 '
            '--np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal --no-exact-upfront '
            "--no-1mm-upfront | cut -f1-9 | sed \'s/$/\t*\t*/' | mxdx demux "
            f'--file-map {out_dir}/files_list.tsv '
            '--batch ${SLURM_ARRAY_TASK_ID} --batch-size 50000000 '
            f'--output-base {out_dir}/alignments --extension sam.xz\n',
            'date\n']
        self.assertEqual(main, exp_main)

        exp_merge = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name merge-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 12\n',
            '#SBATCH --time 25:00:00\n',
            '#SBATCH --mem 80g\n',
            f'#SBATCH --output {out_dir}/merge-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/merge-{job_id}.err\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            'set -e\n',
            "sruns=`grep 'overall alignment rate' *.err | wc -l`\n",
            'if [[ ! -f "errors.log" && $sruns -eq "1" ]]; then\n',
            f'woltka_merge mxdx --base {out_dir}\n',
            f'mkdir -p {out_dir}/bioms\n',
            f'for f in `ls {out_dir}/alignments/*.sam.xz`; do bname=`basename '
            '${f/.sam.xz/}`; '
            f'mkdir -p {out_dir}/bioms/'
            '${bname}; echo woltka classify -i $f -o '
            f'{out_dir}/bioms/'
            '${bname}/none.biom --no-demux --lineage '
            f'{database}.tax --rank none --outcov {out_dir}/coverages/; '
            'done | parallel -j 12\n',
            'wait\n',
            f'woltka_merge biom --base {out_dir}\n',
            f'cd {out_dir};\n',
            '\n',
            'cd alignments; tar -cvf ../alignment.tar *.sam.xz; cd ..;\n',
            'fi\n',
            f'finish_woltka {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(merge, exp_merge)

        # now let's test that if finished correctly
        sdir = 'qp_woltka/support_files/'
        mkdir(f'{out_dir}/woltka')
        copyfile(f'{sdir}/none.biom', f'{out_dir}/none.biom')
        copyfile(f'{sdir}/alignment.tar', f'{out_dir}/alignment.tar')
        copyfile(f'{sdir}/coverages.tgz', f'{out_dir}/coverages.tgz')

        success, ainfo, msg = woltka(
            self.qclient, job_id, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        exp = [
            ArtifactInfo('Per genome Predictions', 'BIOM',
                         [(f'{out_dir}/none.biom', 'biom'),
                          (f'{out_dir}/alignment.tar', 'log'),
                          (f'{out_dir}/none/coverages.tgz', 'plain_text')])]
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
        html_summary = self.qclient.get_artifact_html_summary(aid)
        files['html_summary'] = html_summary

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
            '#SBATCH --mail-type=FAIL,TIME_LIMIT_80,INVALID_DEPEND\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 80:00:00\n',
            '#SBATCH --mem 100g\n',
            f'#SBATCH --output {out_dir}/{job_id}_%a.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_%a.err\n',
            '#SBATCH --array 0-0%12\n',
            '#SBATCH --constraint="amd"\n',
            f'cd {out_dir}\n',
            f'prep_full_path={prep_file}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            f'dbbase={dirname(database)}\n',
            f'dbname={basename(database)}\n',
            f'output={out_dir}\n',
            'bt2_cores=6\n',
            f'mxdx mux --file-map {out_dir}/files_list.tsv --batch '
            '${SLURM_ARRAY_TASK_ID} --batch-size 50000000 '
            '--paired-handling interleave | '
            'bowtie2 -p ${bt2_cores} -x '
            f'{database} --interleaved - --seed 42 --very-sensitive -k 16 '
            '--np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal --no-exact-upfront '
            "--no-1mm-upfront | cut -f1-9 | sed \'s/$/\t*\t*/' | mxdx demux "
            f'--file-map {out_dir}/files_list.tsv '
            '--batch ${SLURM_ARRAY_TASK_ID} --batch-size 50000000 '
            f'--output-base {out_dir}/alignments --extension sam.xz\n',
            'date\n']
        self.assertEqual(main, exp_main)

        exp_merge = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name merge-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 12\n',
            '#SBATCH --time 25:00:00\n',
            '#SBATCH --mem 80g\n',
            f'#SBATCH --output {out_dir}/merge-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/merge-{job_id}.err\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            'set -e\n',
            "sruns=`grep 'overall alignment rate' *.err | wc -l`\n",
            'if [[ ! -f "errors.log" && $sruns -eq "1" ]]; then\n',
            f'woltka_merge mxdx --base {out_dir}\n',
            f'mkdir -p {out_dir}/bioms\n',
            f'for f in `ls {out_dir}/alignments/*.sam.xz`; do bname=`basename '
            '${f/.sam.xz/}`; '
            f'mkdir -p {out_dir}/bioms/'
            '${bname}; echo woltka classify -i $f -o '
            f'{out_dir}/bioms/'
            '${bname}/none.biom --no-demux --lineage '
            f'{database}.tax --rank none --outcov {out_dir}/coverages/; '
            'done | parallel -j 12\n',
            'wait\n',
            f'for f in `ls {out_dir}/alignments/*.sam.xz`; do bname=`basename '
            '${f/.sam.xz/}`; mkdir -p '
            f'{out_dir}/bioms/'
            '${bname}; echo woltka classify -i $f -o '
            f'{out_dir}/bioms/'
            '${bname}/per-gene.biom --no-demux -c '
            f'{database}.coords; done | parallel -j 12\n',
            'wait\n',
            f'woltka_merge biom --base {out_dir}\n',
            f'cd {out_dir};\n',
            '\n',
            'cd alignments; tar -cvf ../alignment.tar *.sam.xz; cd ..;\n',
            'fi\n',
            f'finish_woltka {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(merge, exp_merge)

        # now let's test that if finished correctly
        sdir = 'qp_woltka/support_files/'
        mkdir(f'{out_dir}/woltka')
        copyfile(f'{sdir}/none.biom', f'{out_dir}/none.biom')
        copyfile(f'{sdir}/per-gene.biom', f'{out_dir}/per-gene.biom')
        copyfile(f'{sdir}/alignment.tar', f'{out_dir}/alignment.tar')
        copyfile(f'{sdir}/coverages.tgz', f'{out_dir}/coverages.tgz')
        success, ainfo, msg = woltka(
            self.qclient, job_id, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        exp = [
            ArtifactInfo('Per genome Predictions', 'BIOM',
                         [(f'{out_dir}/none.biom', 'biom'),
                          (f'{out_dir}/alignment.tar', 'log'),
                          (f'{out_dir}/none/coverages.tgz', 'plain_text')]),
            ArtifactInfo('Per gene Predictions', 'BIOM',
                         [(f'{out_dir}/per-gene.biom', 'biom'),
                          (f'{out_dir}/per_gene/coverages.tgz',
                           'plain_text')])]

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

    def test_woltka_syndna_to_array(self):
        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {
                'run_prefix': 'S22205_S104', 'syndna_pool_number': 1,
                'raw_reads_r1r2': 10000, 'mass_syndna_input_ng': 120},
            'SKD8.640184': {
                'run_prefix': 'S22282_S102', 'syndna_pool_number': 1,
                'raw_reads_r1r2': 10002, 'mass_syndna_input_ng': 120}}
        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        files, prep = self.qclient.artifact_and_preparation_files(aid)

        url = 'this-is-my-url'
        database = self.params['Database']
        main_fp, finish_fp = woltka_syndna_to_array(
            files, out_dir, database, prep, url, job_id)

        self.assertEqual(join(out_dir, f'{job_id}.slurm'), main_fp)
        self.assertEqual(join(out_dir, f'{job_id}.finish.slurm'), finish_fp)

        with open(main_fp) as f:
            main = f.readlines()
        with open(finish_fp) as f:
            finish = f.readlines()

        exp_main = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 12:00:00\n',
            '#SBATCH --mem 190g\n',
            f'#SBATCH --output {out_dir}/{job_id}_%a.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_%a.err\n',
            '#SBATCH --array 1-1%12\n',
            f'cd {out_dir}\n',
            'mkdir -p reads/uneven sams\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            f'dbbase={dirname(database)}\n',
            f'dbname={basename(database)}\n',
            f'output={out_dir}\n',
            'while read -r sn f;\n',
            '  do\n',
            '    fn=`basename $f`; \n',
            f'    bowtie2 -p 8 -x {database.rsplit("/", 1)[0]}/pUC57 '
            '-q ${f} --seed 42 --very-sensitive -k 16 --np 1 --mp "1,1" '
            '--rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" '
            '--no-head --no-unal --no-exact-upfront --no-1mm-upfront '
            "| sam_filter -i 0.98 -r 0.90 | awk '{print $1}' > "
            '$PWD/reads/uneven/${fn/.gz/}.seqID.txt; seqkit grep -v '
            '-f $PWD/reads/uneven/${fn/.gz/}.seqID.txt ${f} > '
            '$PWD/reads/uneven/no-plasmid-${fn/.gz/}\n',
            '    bowtie2 -p 8 -x {database} -q '
            '$PWD/reads/uneven/no-plasmid-${fn/.gz/} '
            '-S $PWD/sams/${sn}.sam --seed 42 --very-sensitive -k 16 '
            '--np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min '
            'L,0,-0.05" --no-head --no-unal --no-exact-upfront '
            '--no-1mm-upfront --un $PWD/reads/uneven/${fn/.gz/}\n',
            '  done < sample_details_${SLURM_ARRAY_TASK_ID}.txt\n',
            'date']
        self.assertEqual(main, exp_main)

        exp_finish = [
            '#!/bin/bash\n',
            '#SBATCH -p qiita\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name finish-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 15\n',
            '#SBATCH --mem 1g\n',
            f'#SBATCH --output {out_dir}/finish-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/finish-{job_id}.err\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            "sruns=`grep 'overall alignment rate' *.err | wc -l`\n",
            'sjobs=`ls sams/*.sam | wc -l`\n',
            'if [[ $sruns -eq $((2*sjobs)) ]]; then\n',
            '  mkdir -p sams/final\n',
            '  while read -r fwd rev; do echo "fastq_pair -t 50000000 '
            'reads/uneven/${fwd} reads/uneven/${rev}; mv '
            'reads/uneven/${fwd}.paired.fq reads/${fwd}; '
            'mv reads/uneven/${rev}.paired.fq reads/${rev}; '
            'gzip reads/${fwd} reads/${rev}"; done < '
            'finish_sample_details.txt | parallel -j 8\n',
            '  for f in `ls sams/fwd_*`;\n',
            '    do\n',
            '      fn=`basename $f`;\n',
            '      lines=`head $f | wc -l`\n',
            '      if [[ "$lines" != "0" ]]; then echo cp ${f} '
            'sams/final/${fn:4} ; fi ;\n',
            '  done | parallel -j 8\n',
            '  for f in `ls sams/rev_*`;\n',
            '    do\n',
            '      fn=`basename $f`;\n',
            '      lines=`head $f | wc -l`\n',
            '      if [[ "$lines" != "0" ]]; then echo "cat ${f} >> '
            'sams/final/${fn:4}" ; fi ;\n',
            '  done | parallel -j 8\n',
            '  woltka classify -i sams/final/ -o syndna.biom --no-demux\n',
            '  for f in `ls sams/final/*.sam`;\n',
            '    do\n',
            '      echo "cat ${f} | cut -f1-9 | sed \'s/$/\t*\t*/\' | '
            'xz -1 -T1 > ${f}.xz; rm ${f}"\n',
            '  done | parallel -j 8\n',
            '  cd sams/final/; tar -cvf alignment.tar *.sam.xz; cd ../../;\n',
            'fi\n',
            f'finish_woltka {url} {job_id} {out_dir}\n',
            'set -e\n',
            'date']
        self.assertEqual(finish, exp_finish)

        # now let's test that if finished correctly
        sdir = 'qp_woltka/support_files/'
        copyfile(f'{sdir}/syndna.biom', f'{out_dir}/syndna.biom')
        # we just need the file to exists so is fine to use the biom as tar
        mkdir(f'{out_dir}/sams')
        mkdir(f'{out_dir}/sams/final')
        copyfile(f'{sdir}/none.biom', f'{out_dir}/sams/final/alignment.tar')
        reads_fp = f'{out_dir}/reads/'
        mkdir(reads_fp)
        reads = []
        for i in range(2):
            for j in range(2):
                ft = 'raw_reverse_seqs' if j else 'raw_forward_seqs'
                iname = files[i][j]['filepath']
                jname = f'{reads_fp}{basename(iname)}'
                copyfile(iname, jname)
                reads.append((jname, ft))

        params = self.params.copy()
        params['min_sample_counts'] = 1
        success, ainfo, msg = woltka_syndna(
            self.qclient, job_id, params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        exp = [
            ArtifactInfo('SynDNA hits', 'BIOM',
                         [(f'{out_dir}/syndna.biom', 'biom'),
                          (f'{out_dir}/sams/final/alignment.tar', 'log'),
                          (f'{out_dir}/lin_regress_by_sample_id.yaml', 'log'),
                          (f'{out_dir}/fit_syndna_models_log.txt', 'log')]),
            ArtifactInfo('reads without SynDNA', 'per_sample_FASTQ', reads)]
        self.assertCountEqual(ainfo, exp)

    def test_calculate_cell_counts(self):
        params = {'SynDNA hits': 5, 'Woltka per-genome': 6,
                  'min_coverage': 1, 'read_length': 150,
                  'min_rsquared': 0.8}
        job_id = 'my-job-id'
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # this should fail cause we don't have valid data
        success, ainfo, msg = calculate_cell_counts(
            self.qclient, job_id, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, "No logs found, are you sure you selected the "
                         "correct artifact for 'SynDNA hits'?")

        # this should fail too because but now we are getting deeper into
        # the validation
        prep_info_dict = {
            'SKB8.640193': {
                'run_prefix': 'S22205_S104', 'syndna_pool_number': 1,
                'raw_reads_r1r2': 10000, 'mass_syndna_input_ng': 120,
                'extracted_gdna_concentration_ng_ul': 10,
                'vol_extracted_elution_ul': 5},
            'SKD8.640184': {
                'run_prefix': 'S22282_S102', 'syndna_pool_number': 1,
                'raw_reads_r1r2': 10002, 'mass_syndna_input_ng': 120,
                'extracted_gdna_concentration_ng_ul': 11,
                'vol_extracted_elution_ul': 6}}
        data = {'prep_info': dumps(prep_info_dict),
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        sdir = 'qp_woltka/support_files/'
        fp_to_cp = [
            'fit_syndna_models_log.txt', 'lin_regress_by_sample_id.yaml',
            'syndna.biom']
        for fp in fp_to_cp:
            copyfile(f'{sdir}/{fp}', f'{out_dir}/{fp}')
        data = {
            'filepaths': dumps([
                (f'{out_dir}/syndna.biom', 'biom'),
                (f'{out_dir}/lin_regress_by_sample_id.yaml', 'log'),
                (f'{out_dir}/fit_syndna_models_log.txt', 'log')]),
            'type': "BIOM",
            'name': "SynDNA Hits - Test",
            'prep': pid}
        params['SynDNA hits'] = self.qclient.post(
            '/apitest/artifact/', data=data)['artifact']

        success, ainfo, msg = calculate_cell_counts(
            self.qclient, job_id, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, "The selected 'Woltka per-genome' artifact "
                         "doesn't look like one, did you select the correct "
                         "file?")

        # Finally, adding a full test is close to impossible - too many steps.

    def test_calculate_rna_copy_counts(self):
        params = {'Woltka per-gene': 6}
        job_id = 'my-job-id'
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # this should fail cause we don't have valid data
        success, ainfo, msg = calculate_rna_copy_counts(
            self.qclient, job_id, params, out_dir)
        self.assertFalse(success)
        self.assertEqual(msg, "The selected 'Woltka per-gene' artifact "
                         "doesn't look like one, did you select the "
                         "correct file?")

        # Finally, adding a full test is close to impossible - too many steps.


if __name__ == '__main__':
    main()
