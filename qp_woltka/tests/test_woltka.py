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
from os.path import exists, isdir, join, dirname
from shutil import rmtree, copyfile
from tempfile import mkdtemp
from json import dumps
from biom import load_table

from qp_woltka import plugin
from qp_woltka.util import get_dbs, generate_woltka_dflt_params
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
                'command': dumps(['qp-woltka', '2020.11', 'Woltka v0.1.1']),
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

        # retriving info of the prep/artifact just created
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)
        directory = {dirname(ffs) for _, fs in artifact_info['files'].items()
                     for ffs in fs}
        directory = directory.pop()
        prep_info = artifact_info['prep_information']
        prep_info = self.qclient.get(
            '/qiita_db/prep_template/%s/' % prep_info[0])
        prep_file = prep_info['prep-file']

        url = 'this-is-my-url'
        database = self.params['Database']
        main_qsub_fp, merge_qsub_fp = woltka_to_array(
            directory, out_dir, database, prep_file, url, job_id)

        self.assertEqual(join(out_dir, f'{job_id}.qsub'), main_qsub_fp)
        self.assertEqual(join(out_dir, f'{job_id}.merge.qsub'), merge_qsub_fp)

        with open(main_qsub_fp) as f:
            main_qsub = f.readlines()
        with open(merge_qsub_fp) as f:
            merge_qsub = f.readlines()

        exp_main_qsub = [
            '#!/bin/bash\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 30:00:00\n',
            '#SBATCH --mem 90g\n',
            f'#SBATCH --output {out_dir}/{job_id}_'
            '${SLURM_ARRAY_TASK_ID}.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_'
            '${SLURM_ARRAY_TASK_ID}.err\n',
            '#SBATCH --array 1-2%8\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            'offset=${SLURM_ARRAY_TASK_ID}\n',
            'step=$(( $offset - 0 ))\n',
            'if [[ $step -gt 2 ]]; then exit 0; fi\n',
            f'args0=$(head -n $step {out_dir}/{job_id}.array-details'
            ' | tail -n 1)\n',
            "infile0=$(echo -e $args0 | awk '{ print $1 }')\n",
            "outfile0=$(echo -e $args0 | awk '{ print $2 }')\n",
            'set -e\n',
            'cat $infile0*.fastq.gz > $outfile0.fastq.gz; bowtie2 -p 8 -x '
            f'{database} -q $outfile0.fastq.gz -S $outfile0.sam --seed 42 '
            '--very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal; woltka classify -i '
            '$outfile0.sam -o $outfile0.woltka-taxa --no-demux --lineage '
            f'{database}.tax --rank phylum,genus,species,free,none; xz -9 -T8 '
            '-c $outfile0.sam > $outfile0.xz\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main_qsub, exp_main_qsub)

        exp_merge_qsub = [
            '#!/bin/bash\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name merge-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 5\n',
            '#SBATCH --time 10:00:00\n',
            '#SBATCH --mem 48g\n',
            f'#SBATCH --output {out_dir}/merge-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/merge-{job_id}.err\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            'set -e\n',
            "PROCESS=1; COUNTER=0; for f in `awk '{print $NF}' "
            f'{out_dir}/*.array-details`; do let COUNTER=COUNTER+1; '
            "if [ ! -f ${f}*/species.biom ]; then if ! grep -xq "
            "'0.00% overall alignment rate' *_${COUNTER}.err-${COUNTER}; "
            "then PROCESS=0; fi; fi; done\n",
            "if [ 1 -eq $PROCESS ]; then \n",
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'phylum --glob "*.woltka-taxa/phylum.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'genus --glob "*.woltka-taxa/genus.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'species --glob "*.woltka-taxa/species.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'free --glob "*.woltka-taxa/free.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'none --glob "*.woltka-taxa/none.biom" --rename &\n',
            'wait\n',
            f'cd {out_dir}; tar -cvf alignment.tar *.sam.xz\n',
            'fi\n',
            f'finish_woltka {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(merge_qsub, exp_merge_qsub)

        # now let's test that if finished correctly
        sdir = 'qp_woltka/support_files/'
        copyfile(f'{sdir}/genus.biom', f'{out_dir}/genus.biom')
        copyfile(f'{sdir}/none.biom', f'{out_dir}/none.biom')
        copyfile(f'{sdir}/species.biom', f'{out_dir}/species.biom')
        copyfile(f'{sdir}/phylum.biom', f'{out_dir}/phylum.biom')
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
            ArtifactInfo('Taxonomic Predictions - phylum', 'BIOM',
                         [(f'{out_dir}/phylum.biom', 'biom')]),
            ArtifactInfo('Taxonomic Predictions - genus', 'BIOM',
                         [(f'{out_dir}/genus.biom', 'biom')]),
            ArtifactInfo('Taxonomic Predictions - species', 'BIOM',
                         [(f'{out_dir}/species.biom', 'biom')]),
            ArtifactInfo('Per genome Predictions', 'BIOM',
                         [(f'{out_dir}/none.biom', 'biom')])]

        self.assertCountEqual(ainfo, exp)

    def test_woltka_to_array_wol(self):
        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104_L001_R1'},
            'SKD8.640184': {'run_prefix': 'S22282_S102_L001_R1'}}
        database = join(self.db_path, 'wol/WoLmin')

        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict, database)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # retriving info of the prep/artifact just created
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)
        directory = {dirname(ffs) for _, fs in artifact_info['files'].items()
                     for ffs in fs}
        directory = directory.pop()
        prep_info = artifact_info['prep_information']
        prep_info = self.qclient.get(
            '/qiita_db/prep_template/%s/' % prep_info[0])
        prep_file = prep_info['prep-file']

        url = 'this-is-my-url'
        main_qsub_fp, merge_qsub_fp = woltka_to_array(
            directory, out_dir, database, prep_file, url, job_id)

        self.assertEqual(join(out_dir, f'{job_id}.qsub'), main_qsub_fp)
        self.assertEqual(join(out_dir, f'{job_id}.merge.qsub'), merge_qsub_fp)

        with open(main_qsub_fp) as f:
            main_qsub = f.readlines()
        with open(merge_qsub_fp) as f:
            merge_qsub = f.readlines()

        exp_main_qsub = [
            '#!/bin/bash\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name {job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 8\n',
            '#SBATCH --time 30:00:00\n',
            '#SBATCH --mem 90g\n',
            f'#SBATCH --output {out_dir}/{job_id}_'
            '${SLURM_ARRAY_TASK_ID}.log\n',
            f'#SBATCH --error {out_dir}/{job_id}_'
            '${SLURM_ARRAY_TASK_ID}.err\n',
            '#SBATCH --array 1-2%8\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}\n',
            'offset=${SLURM_ARRAY_TASK_ID}\n',
            'step=$(( $offset - 0 ))\n',
            'if [[ $step -gt 2 ]]; then exit 0; fi\n',
            f'args0=$(head -n $step {out_dir}/{job_id}.array-details'
            ' | tail -n 1)\n',
            "infile0=$(echo -e $args0 | awk '{ print $1 }')\n",
            "outfile0=$(echo -e $args0 | awk '{ print $2 }')\n",
            'set -e\n',
            'cat $infile0*.fastq.gz > $outfile0.fastq.gz; bowtie2 -p 8 -x '
            f'{database} -q $outfile0.fastq.gz -S $outfile0.sam --seed 42 '
            '--very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal; woltka classify '
            '-i $outfile0.sam -o $outfile0.woltka-taxa --no-demux '
            f'--lineage {database}.tax --rank phylum,genus,species,free,none; '
            f'woltka classify -i $outfile0.sam -c {database}.coords '
            '-o $outfile0.woltka-per-gene --no-demux; xz -9 -T8 -c '
            '$outfile0.sam > $outfile0.xz\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main_qsub, exp_main_qsub)

        exp_merge_qsub = [
            '#!/bin/bash\n',
            '#SBATCH --mail-user "qiita.help@gmail.com"\n',
            f'#SBATCH --job-name merge-{job_id}\n',
            '#SBATCH -N 1\n',
            '#SBATCH -n 6\n',
            '#SBATCH --time 10:00:00\n',
            '#SBATCH --mem 48g\n',
            f'#SBATCH --output {out_dir}/merge-{job_id}.log\n',
            f'#SBATCH --error {out_dir}/merge-{job_id}.err\n',
            f'cd {out_dir}\n',
            f'{self.environment}\n',
            'date\n',
            'hostname\n',
            'echo $SLURM_JOBID\n',
            'set -e\n',
            "PROCESS=1; COUNTER=0; for f in `awk '{print $NF}' "
            f'{out_dir}/*.array-details`; do let COUNTER=COUNTER+1; '
            "if [ ! -f ${f}*/species.biom ]; then if ! grep -xq "
            "'0.00% overall alignment rate' *_${COUNTER}.err-${COUNTER}; "
            "then PROCESS=0; fi; fi; done\n",
            "if [ 1 -eq $PROCESS ]; then \n",
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'phylum --glob "*.woltka-taxa/phylum.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'genus --glob "*.woltka-taxa/genus.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'species --glob "*.woltka-taxa/species.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'free --glob "*.woltka-taxa/free.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'none --glob "*.woltka-taxa/none.biom" &\n',
            f'woltka_merge --prep {prep_file} --base {out_dir}  --name '
            'per-gene --glob "*.woltka-per-gene" --rename &\n',
            'wait\n',
            f'cd {out_dir}; tar -cvf alignment.tar *.sam.xz\n',
            'fi\n',
            f'finish_woltka {url} {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(merge_qsub, exp_merge_qsub)

        # now let's test that if finished correctly
        sdir = 'qp_woltka/support_files/'
        copyfile(f'{sdir}/genus.biom', f'{out_dir}/genus.biom')
        copyfile(f'{sdir}/none.biom', f'{out_dir}/none.biom')
        copyfile(f'{sdir}/per-gene.biom', f'{out_dir}/per-gene.biom')
        copyfile(f'{sdir}/species.biom', f'{out_dir}/species.biom')
        copyfile(f'{sdir}/phylum.biom', f'{out_dir}/phylum.biom')
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
            ArtifactInfo('Taxonomic Predictions - phylum', 'BIOM',
                         [(f'{out_dir}/phylum.biom', 'biom')]),
            ArtifactInfo('Taxonomic Predictions - genus', 'BIOM',
                         [(f'{out_dir}/genus.biom', 'biom')]),
            ArtifactInfo('Taxonomic Predictions - species', 'BIOM',
                         [(f'{out_dir}/species.biom', 'biom')]),
            ArtifactInfo('Per genome Predictions', 'BIOM',
                         [(f'{out_dir}/none.biom', 'biom')]),
            ArtifactInfo('Per gene Predictions', 'BIOM',
                         [(f'{out_dir}/per-gene.biom', 'biom')])]

        self.assertCountEqual(ainfo, exp)

        # check that the produced table have feature taxonomy
        bt = load_table(f'{out_dir}/phylum.biom')
        self.assertCountEqual(
            bt.metadata_to_dataframe('observation').columns,
            ['taxonomy_0', 'taxonomy_1'])

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
            'Table phylum was not created, please contact '
            'qiita.help@gmail.com for more information',
            'Table genus was not created, please contact qiita.help@gmail.com '
            'for more information',
            'Table species was not created, please contact '
            'qiita.help@gmail.com for more information',
            'Table none/per-genome was not created, please contact '
            'qiita.help@gmail.com for more information',
            'Table per-gene was not created, please contact '
            'qiita.help@gmail.com for more information'])
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

        # retriving info of the prep/artifact just created
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)
        directory = {dirname(ffs) for _, fs in artifact_info['files'].items()
                     for ffs in fs}
        directory = directory.pop()
        prep_info = artifact_info['prep_information']
        prep_info = self.qclient.get(
            '/qiita_db/prep_template/%s/' % prep_info[0])
        prep_file = prep_info['prep-file']

        url = 'this-is-my-url'
        with self.assertRaises(ValueError) as error:
            woltka_to_array(directory, out_dir, self.params['Database'],
                            prep_file, url, job_id)
        self.assertEqual(str(error.exception), "Prep information is missing "
                         "the required run_prefix column")

    def test_creation_error_no_unique_run_prefix(self):
        # run prefix doesn't exist
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22205_S104'}}
        pid, aid, job_id = self._helper_woltka_bowtie(prep_info_dict)

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # retriving info of the prep/artifact just created
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)
        directory = {dirname(ffs) for _, fs in artifact_info['files'].items()
                     for ffs in fs}
        directory = directory.pop()
        prep_info = artifact_info['prep_information']
        prep_info = self.qclient.get(
            '/qiita_db/prep_template/%s/' % prep_info[0])
        prep_file = prep_info['prep-file']

        url = 'this-is-my-url'
        with self.assertRaises(ValueError) as error:
            woltka_to_array(directory, out_dir, self.params['Database'],
                            prep_file, url, job_id)
        self.assertEqual(str(error.exception), "The run_prefix values are "
                         "not unique for each sample")


if __name__ == '__main__':
    main()
