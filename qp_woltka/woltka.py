# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import re
from os import environ, mkdir
from os.path import join, basename, exists, dirname
from glob import glob
from shutil import copy2
from math import ceil
from biom import load_table
from biom.util import biom_open
import pandas as pd
from tarfile import open as topen
from pysyndna import fit_linear_regression_models_for_qiita
from pysyndna import calc_ogu_cell_counts_per_g_of_sample_for_qiita
from pysyndna import calc_copies_of_ogu_orf_ssrna_per_g_sample_for_qiita
from pysyndna import OGU_ID_KEY, OGU_PERCENT_COVERAGE_KEY

from qp_woltka.util import search_by_filename

from qiita_client import ArtifactInfo
from qiita_client.util import system_call

# resources per job
PPN = 8
MAX_RUNNING = 12
TASKS_IN_SCRIPT = 10

MEMORY = '100g'
LARGE_MEMORY = '180g'
MERGE_MEMORY = '80g'
SYNDNA_MEMORY = '190g'
# setting so an iSeq run, generates 2 jobs
BATCHSIZE = 50000000

WALLTIME = '80:00:00'
MERGE_WALLTIME = '25:00:00'
SYNDNA_WALLTIME = '12:00:00'


def _process_database_files(database_fp):
    files = glob(f'{database_fp}*')

    database_files = {
        'taxonomy': None,
        'gene_coordinates': None,
        'kegg': {
            'orf-to-ko.map.xz': None,
            'ko-to-ec.map': None,
            'ko-to-reaction.map': None,
            'reaction-to-module.map': None,
            'module-to-pathway.map': None
        }
    }
    database_files['taxonomy'] = [f for f in files if f.endswith('.tax')][0]
    gene_coordinates = [f for f in files if f.endswith('.coords')]
    # not all databases have their coordinates fp
    if gene_coordinates:
        database_files['gene_coordinates'] = gene_coordinates[0]
        # if there are gene_coordinates, there might be function translations
        dname = dirname(database_fp)
        if f'{dname}/function' in glob(f'{dname}/*'):
            if f'{dname}/function/kegg' in glob(f'{dname}/function/*'):
                dt = f'{dname}/function/kegg'
                files = glob(f'{dt}/*')
                for k in database_files['kegg'].keys():
                    if f'{dt}/{k}' in files:
                        database_files['kegg'][k] = f'{dt}/{k}'
    lmap_fp = f'{dirname(database_fp)}/genomes/length.map'
    if exists(lmap_fp):
        database_files['length.map'] = lmap_fp

    return database_files


def woltka_to_array(files, output, database_bowtie2, prep, url, name):
    """Creates files for submission of per sample bowtie2 and woltka
    """
    environment = environ["ENVIRONMENT"]

    # processing html_summary
    html_summary = files.pop('html_summary')
    try:
        df = pd.read_html(html_summary)[0]
    except ValueError:
        txt = ('The summary table could not parsed; please send an email '
               'to qiita.help@ucsd.edu')
        raise ValueError(txt)

    dname = dirname(html_summary)
    rev_exists = 'raw_reverse_seqs' in df.file_type.unique()

    fwd = dict(df[df.file_type == 'raw_forward_seqs'].apply(
        lambda x: (x.filename.rsplit('_R1')[0],
                   (x.filename, x.reads)), axis=1).values)
    if rev_exists:
        rev = dict(df[df.file_type == 'raw_reverse_seqs'].apply(
            lambda x: (x.filename.rsplit('_R2')[0],
                       (x.filename, x.reads)), axis=1).values)
    # let's check that there is some overlap and if not try something different
    if rev_exists and not set(fwd) & set(rev):
        fwd = dict(df[df.file_type == 'raw_forward_seqs'].apply(
            lambda x: (x.filename.rsplit('.R1.')[0],
                       (x.filename, x.reads)), axis=1).values)
        rev = dict(df[df.file_type == 'raw_reverse_seqs'].apply(
            lambda x: (x.filename.rsplit('.R2.')[0],
                       (x.filename, x.reads)), axis=1).values)
        if not set(fwd) & set(rev):
            raise ValueError('There is no overlap between fwd/rev reads, if '
                             'you think that not correct please send an email '
                             'to qiita.help@gmail.com')
    if rev_exists:
        failed_reads = []
        lines = ['filename_1\tfilename_2\trecord_count']
        for k, (fn, reads) in fwd.items():
            rfn, rreads = rev.pop(k)
            if int(rreads) != int(reads):
                failed_reads.append(f'{basename(fn)} {basename(rfn)}')
            lines.append(f'{dname}/{fn}\t{dname}/{rfn}\t{reads}')
        if failed_reads:
            failed_reads = '\n'.join(failed_reads)
            raise ValueError(
                'Some of the fwd/rev do not have the same number of reads; '
                'are you using an artifact created with a newer command?\n\n'
                f'Failed files:\n {failed_reads}')
    else:
        lines = ['filename_1\trecord_count']
        for k, (fn, reads) in fwd.items():
            lines.append(f'{dname}/{fn}\t{reads}')

    files_list_fp = f'{output}/files_list.tsv'
    with open(files_list_fp, 'w') as fp:
        fp.write('\n'.join(lines))

    cmd = (f'mxdx get-max-batch-number --file-map {files_list_fp} '
           f'--batch-size {BATCHSIZE}')
    n_files, stderr, return_value = system_call(cmd)
    if return_value != 0 or stderr:
        raise ValueError('`mxdx get-max-batch-number` failed '
                         f'{return_value}: {stderr}')
    # just making sure that n_files is an int
    n_files = int(n_files)

    db_files = _process_database_files(database_bowtie2)
    db_folder = dirname(database_bowtie2)
    db_name = basename(database_bowtie2)

    woltka_merge = f'woltka_merge mxdx --base {output}'
    extra_commands = ''
    if 'length.map' in db_files:
        woltka_merge += f' --length_map {db_files["length.map"]}'
        extra_commands = (
            'find ${PWD}/coverages -iname "*.cov" > ${PWD}/cov_files.txt\n'
            'micov consolidate --paths ${PWD}/cov_files.txt --lengths '
            f'{db_files["length.map"]} --output '
            '${PWD}/coverages.tgz\n')

    woltka_cmds = [
        # creating the output folder
        f'mkdir -p {output}/bioms',
        # executing the parallel classify
        f'for f in `ls {output}/alignments/*.sam.xz`; '
        'do bname=`basename ${f/.sam.xz/}`; '
        f'mkdir -p {output}/bioms/' '${bname}; echo woltka classify -i $f '
        f'-o {output}/bioms/' '${bname}/none.biom --no-demux --lineage '
        f'{db_files["taxonomy"]} --rank none --outcov {output}/coverages/; '
        'done | parallel -j 12',
        'wait']

    if db_files['gene_coordinates']:
        woltka_cmds.append(
            f'for f in `ls {output}/alignments/*.sam.xz`; '
            'do bname=`basename ${f/.sam.xz/}`; '
            f'mkdir -p {output}/bioms/' '${bname}; echo '
            'woltka classify -i $f '
            f'-o {output}/'
            'bioms/${bname}/per-gene.biom --no-demux -c '
            f'{db_files["gene_coordinates"]}; done | parallel -j 12')
        woltka_cmds.append('wait')
        woltka_cmds.append(f'woltka_merge biom --base {output}')

        wcdm = 'woltka collapse -i '
        dbfk = db_files['kegg']
        if dbfk["orf-to-ko.map.xz"] is not None:
            woltka_cmds.append(f'{wcdm} per-gene.biom -m '
                               f'{dbfk["orf-to-ko.map.xz"]} -o ko.biom')
        if dbfk["ko-to-ec.map"] is not None:
            woltka_cmds.append(f'{wcdm} ko.biom -m {dbfk["ko-to-ec.map"]} '
                               '-o ec.biom')
        if dbfk["ko-to-reaction.map"] is not None and \
                dbfk["reaction-to-module.map"] is not None and \
                dbfk["module-to-pathway.map"] is not None:
            woltka_cmds.append(
                f'{wcdm} ko.biom -m {dbfk["ko-to-reaction.map"]} '
                '-o reaction.biom')
            woltka_cmds.append(
                f'{wcdm} reaction.biom -m '
                f'{dbfk["reaction-to-module.map"]} -o module.biom')
            woltka_cmds.append(
                f'{wcdm} module.biom -m '
                f'{dbfk["module-to-pathway.map"]} -o pathway.biom')
    else:
        woltka_cmds.append(f'woltka_merge biom --base {output}')

    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             f'#SBATCH --job-name merge-{name}',
             '#SBATCH -N 1',
             '#SBATCH -n 12',
             f'#SBATCH --time {MERGE_WALLTIME}',
             f'#SBATCH --mem {MERGE_MEMORY}',
             f'#SBATCH --output {output}/merge-{name}.log',
             f'#SBATCH --error {output}/merge-{name}.err',
             f'cd {output}',
             f'{environment}',
             'date',  # start time
             'hostname',  # executing system
             'echo $SLURM_JOBID',
             'set -e',
             # if the error file doesn't exi and the number of alignment
             # reports is equal to the numbers of jobs that started : process
             # the bioms
             "sruns=`grep 'overall alignment rate' *.err | wc -l`",
             f'if [[ ! -f "errors.log" && $sruns -eq "{n_files + 1}" ]]; then',
             woltka_merge,
             '\n'.join(woltka_cmds),
             f'cd {output};',
             extra_commands,
             'cd alignments; tar -cvf ../alignment.tar *.sam.xz; cd ..;\n'
             'fi',
             f'finish_woltka {url} {name} {output}\n'
             "date"]  # end time
    # write out the merge script
    merge_fp = join(output, f'{name}.merge.slurm')
    with open(merge_fp, 'w') as out:
        out.write('\n'.join(lines))
        out.write('\n')

    # Bowtie2 command structure based on
    # https://github.com/BenLangmead/bowtie2/issues/311
    preparation_information = join(output, 'prep_info.tsv')
    prep.set_index('sample_name').to_csv(preparation_information, sep='\t')
    if rev_exists:
        bowtie2 = (
            f'mxdx mux --file-map {files_list_fp} --batch '
            '${SLURM_ARRAY_TASK_ID} '
            f'--batch-size {BATCHSIZE} --paired-handling interleave | '
            'bowtie2 -p ${bt2_cores} '
            f'-x {database_bowtie2} --interleaved - --seed 42 '
            '--very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal --no-exact-upfront '
            "--no-1mm-upfront | cut -f1-9 | sed 's/$/\t*\t*/' | "
            f'mxdx demux --file-map {files_list_fp} '
            '--batch ${SLURM_ARRAY_TASK_ID} '
            f'--batch-size {BATCHSIZE} --output-base {output}/alignments '
            '--extension sam.xz')
    else:
        bowtie2 = (
            f'mxdx mux --file-map {files_list_fp} --batch '
            '${SLURM_ARRAY_TASK_ID} '
            f'--batch-size {BATCHSIZE} | '
            'bowtie2 -p ${bt2_cores} '
            f'-x {database_bowtie2} -q - --seed 42 '
            '--very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" '
            '--score-min "L,0,-0.05" --no-head --no-unal --no-exact-upfront '
            "--no-1mm-upfront | cut -f1-9 | sed 's/$/\t*\t*/' | "
            f'mxdx demux --file-map {files_list_fp} '
            '--batch ${SLURM_ARRAY_TASK_ID} '
            f'--batch-size {BATCHSIZE} --output-base {output}/alignments '
            '--extension sam.xz')

    memory = MEMORY
    if 'RS2' in database_bowtie2:
        memory = LARGE_MEMORY

    # all the setup pieces
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             '#SBATCH --mail-type=FAIL,TIME_LIMIT_80,INVALID_DEPEND',
             f'#SBATCH --job-name {name}',
             '#SBATCH -N 1',
             f'#SBATCH -n {PPN}',
             f'#SBATCH --time {WALLTIME}',
             f'#SBATCH --mem {memory}',
             f'#SBATCH --output {output}/{name}_%a.log',
             f'#SBATCH --error {output}/{name}_%a.err',
             f'#SBATCH --array 0-{n_files}%{MAX_RUNNING}',
             '#SBATCH --constraint="amd"',
             f'cd {output}',
             f'prep_full_path={preparation_information}',
             f'{environment}',
             'date',  # start time
             'hostname',  # executing system
             'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}',
             f'dbbase={db_folder}',
             f'dbname={db_name}',
             f'output={output}',
             f'bt2_cores={PPN - 2}',
             bowtie2]

    lines.append('date')  # end time

    # write out the script
    main_fp = join(output, f'{name}.slurm')
    with open(main_fp, 'w') as job:
        job.write('\n'.join(lines))
        job.write('\n')

    return main_fp, merge_fp


def woltka(qclient, job_id, parameters, out_dir):
    """Run Woltka with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    db_files = _process_database_files(parameters['Database'])

    def _coverage_copy(dest):
        fp_coverages = join(out_dir, 'coverages.tgz')
        mkdir(dest)
        dest = join(dest, 'coverages.tgz')
        copy2(fp_coverages, dest)

        return dest

    errors = []
    ainfo = []

    fp_biom = f'{out_dir}/none.biom'
    fp_alng = f'{out_dir}/alignment.tar'
    if exists(fp_biom) and exists(fp_alng):
        ainfo.append(ArtifactInfo('Per genome Predictions', 'BIOM', [
            (fp_biom, 'biom'), (fp_alng, 'log'),
            (_coverage_copy(f'{out_dir}/none/'), 'plain_text')]))
    else:
        errors.append('Table none/per-genome was not created, please contact '
                      'qiita.help@gmail.com for more information')

    if db_files['gene_coordinates'] is not None:
        fp_biom = f'{out_dir}/per-gene.biom'
        if exists(fp_biom):
            ainfo.append(ArtifactInfo('Per gene Predictions', 'BIOM', [
                (fp_biom, 'biom'),
                (_coverage_copy(f'{out_dir}/per_gene/'), 'plain_text')]))
        else:
            errors.append('Table per-gene was not created, please contact '
                          'qiita.help@gmail.com for more information')

        dbfk = db_files['kegg']

        if dbfk["orf-to-ko.map.xz"] is not None:
            fp_biom = f'{out_dir}/ko.biom'
            if exists(fp_biom):
                ainfo.append(ArtifactInfo('KEGG Ontology (KO)', 'BIOM', [
                    (fp_biom, 'biom'),
                    (_coverage_copy(f'{out_dir}/ko/'), 'plain_text')]))
            else:
                errors.append('Table KEGG Ontology was not created, please '
                              'contact qiita.help@gmail.com for more '
                              'information')

        if dbfk["ko-to-ec.map"] is not None:
            fp_biom = f'{out_dir}/ec.biom'
            if exists(fp_biom):
                ainfo.append(ArtifactInfo('KEGG Enzyme (EC)', 'BIOM', [
                    (fp_biom, 'biom'),
                    (_coverage_copy(f'{out_dir}/ec/'), 'plain_text')]))
            else:
                errors.append('Table KEGG Enzyme was not created, please '
                              'contact qiita.help@gmail.com for more '
                              'information')

        if dbfk["ko-to-reaction.map"] is not None and \
                dbfk["reaction-to-module.map"] is not None and \
                dbfk["module-to-pathway.map"] is not None:
            fp_biom = f'{out_dir}/pathway.biom'
            if exists(fp_biom):
                ainfo.append(
                    ArtifactInfo('KEGG Pathway', 'BIOM', [
                        (fp_biom, 'biom'),
                        (_coverage_copy(f'{out_dir}/pathway/'),
                         'plain_text')]))
            else:
                errors.append('Table KEGG Pathway was not created, please '
                              'contact qiita.help@gmail.com for more '
                              'information')

    if errors:
        return False, ainfo, '\n'.join(errors)
    else:

        return True, ainfo, ""


def woltka_syndna_to_array(files, output, database_bowtie2, prep, url, name):
    """Creates files for submission of per sample bowtie2 and woltka_syndna
    """
    environment = environ["ENVIRONMENT"]
    db_folder = dirname(database_bowtie2)
    db_name = basename(database_bowtie2)

    # storing the prep so we can use later
    preparation_information = join(output, 'prep_info.tsv')
    prep.set_index('sample_name').to_csv(preparation_information, sep='\t')

    lookup = prep.set_index('run_prefix')['sample_name'].to_dict()
    n_files = 1
    for i, (k, (f, r)) in enumerate(files.items()):
        if i >= n_files*TASKS_IN_SCRIPT:
            n_files += 1

        sname = search_by_filename(basename(f['filepath']), lookup)
        line = f'fwd_{sname} {f["filepath"]}\n'
        fline = f'{basename(f["filepath"])[:-3]}\n'
        fastq_pair_cmd = ('  while read -r fwd; do echo "'
                          'mv reads/uneven/${fwd} reads/${fwd}; '
                          'gzip reads/${fwd} reads/${rev}";')
        if r is not None:
            line += f'rev_{sname} {r["filepath"]}\n'
            fline = (f'{basename(f["filepath"])[:-3]}\t'
                     f'{basename(r["filepath"])[:-3]}\n')
            fastq_pair_cmd = (
                '  while read -r fwd rev; do echo "fastq_pair -t 50000000 '
                'reads/uneven/${fwd} reads/uneven/${rev}; '
                'mv reads/uneven/${fwd}.paired.fq reads/${fwd}; '
                'mv reads/uneven/${rev}.paired.fq reads/${rev}; '
                'gzip reads/${fwd} reads/${rev}";')

        with open(join(output, 'finish_sample_details.txt'), 'a+') as fh:
            fh.write(fline)

        with open(join(output, f'sample_details_{n_files}.txt'), 'a+') as fh:
            fh.write(line)

    # the plasmid database should live in the same location than
    # synDNA_metagenomic
    plamid_db = join(database_bowtie2.rsplit('/', 1)[0], 'pUC57')
    bowtie2_plasmids = f'bowtie2 -p {PPN} -x {plamid_db} ' + \
        '-q ${f} ' +\
        '--seed 42 --very-sensitive -k 16 --np 1 --mp "1,1" ' + \
        '--rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" ' + \
        '--no-head --no-unal --no-exact-upfront --no-1mm-upfront' + \
        " | sam_filter -i 0.98 -r 0.90 | awk '{print $1}' > " + \
        '$PWD/reads/uneven/${fn/.gz/}.seqID.txt; seqkit grep -v -f ' + \
        '$PWD/reads/uneven/${fn/.gz/}.seqID.txt ${f} > ' + \
        '$PWD/reads/uneven/no-plasmid-${fn/.gz/}'

    bowtie2_inserts = f'bowtie2 -p {PPN} -x {database_bowtie2} ' + \
        '-q $PWD/reads/uneven/no-plasmid-${fn/.gz/} -S $PWD/sams/${sn}.sam ' +\
        '--seed 42 --very-sensitive -k 16 --np 1 --mp "1,1" ' + \
        '--rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" ' + \
        '--no-head --no-unal --no-exact-upfront --no-1mm-upfront ' + \
        '--un $PWD/reads/uneven/${fn/.gz/}'

    # all the setup pieces
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             f'#SBATCH --job-name {name}',
             '#SBATCH -N 1',
             f'#SBATCH -n {PPN}',
             f'#SBATCH --time {SYNDNA_WALLTIME}',
             f'#SBATCH --mem {SYNDNA_MEMORY}',
             f'#SBATCH --output {output}/{name}_%a.log',
             f'#SBATCH --error {output}/{name}_%a.err',
             f'#SBATCH --array 1-{n_files}%{MAX_RUNNING}',
             f'cd {output}',
             'mkdir -p reads/uneven sams',
             f'{environment}',
             'date',  # start time
             'hostname',  # executing system
             'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}',
             f'dbbase={db_folder}',
             f'dbname={db_name}',
             f'output={output}',
             'while read -r sn f;',
             '  do',
             '    fn=`basename $f`; ',
             f'    {bowtie2_plasmids}',
             f'    {bowtie2_inserts}',
             '  done < sample_details_${SLURM_ARRAY_TASK_ID}.txt',
             'date']

    # write out the script
    main_fp = join(output, f'{name}.slurm')
    with open(main_fp, 'w') as job:
        job.write('\n'.join(lines))

    memory = ceil(n_files * .8)
    time = ceil(n_files * 15)
    # creating finish job
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             f'#SBATCH --job-name finish-{name}',
             '#SBATCH -N 1',
             f'#SBATCH -n {PPN}',
             f'#SBATCH --time {time}',  # this is in minutes
             f'#SBATCH --mem {memory}g',
             f'#SBATCH --output {output}/finish-{name}.log',
             f'#SBATCH --error {output}/finish-{name}.err',
             f'cd {output}',
             f'{environment}',
             'date',  # start time
             'hostname',  # executing system
             'echo $SLURM_JOBID',
             "sruns=`grep 'overall alignment rate' *.err | wc -l`",
             'sjobs=`ls sams/*.sam | wc -l`',
             'if [[ $sruns -eq $((2*sjobs)) ]]; then',
             '  mkdir -p sams/final',
             f'{fastq_pair_cmd} done < '
             f'finish_sample_details.txt | parallel -j {PPN}',
             '  for f in `ls sams/fwd_*`;',
             '    do',
             '      fn=`basename $f`;',
             '      lines=`head $f | wc -l`',
             '      if [[ "$lines" != "0" ]]; then echo cp ${f} '
             'sams/final/${fn:4} ; fi ;',
             f'  done | parallel -j {PPN}',
             '  for f in `ls sams/rev_*`;',
             '    do',
             '      fn=`basename $f`;',
             '      lines=`head $f | wc -l`',
             '      if [[ "$lines" != "0" ]]; then echo "cat ${f} >> '
             'sams/final/${fn:4}" ; fi ;',
             f'  done | parallel -j {PPN}',
             '  woltka classify -i sams/final/ -o syndna.biom --no-demux',
             '  for f in `ls sams/final/*.sam`;',
             '    do',
             '      echo "cat ${f} | cut -f1-9 | sed \'s/$/\t*\t*/\' | '
             'xz -1 -T1 > ${f}.xz; rm ${f}"',
             f'  done | parallel -j {PPN}',
             '  cd sams/final/; tar -cvf alignment.tar *.sam.xz; cd ../../;',
             'fi',
             f'finish_woltka {url} {name} {output}',
             'set -e',
             "date"]

    finish_fp = join(output, f'{name}.finish.slurm')
    with open(finish_fp, 'w') as job:
        job.write('\n'.join(lines))

    return main_fp, finish_fp


def woltka_syndna(qclient, job_id, parameters, out_dir):
    """Run Woltka againts the SynDNA default

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to wolka syndna
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    errors = []
    ainfo = []
    fp_biom = f'{out_dir}/syndna.biom'
    fp_alng = f'{out_dir}/sams/final/alignment.tar'
    if exists(fp_biom) and exists(fp_alng):
        # if we got to this point a preparation file should exist in
        # the output folder
        prep = pd.read_csv(
            f'{out_dir}/prep_info.tsv', index_col=None, sep='\t')
        output = fit_linear_regression_models_for_qiita(
            prep, load_table(fp_biom), int(parameters['min_sample_counts']))
        # saving results to disk
        lin_regress_results_fp = f'{out_dir}/lin_regress_by_sample_id.yaml'
        fit_syndna_models_log_fp = f'{out_dir}/fit_syndna_models_log.txt'
        with open(lin_regress_results_fp, 'w') as fp:
            fp.write(output['lin_regress_by_sample_id'])
        with open(fit_syndna_models_log_fp, 'w') as fp:
            fp.write(output['fit_syndna_models_log'])
        ainfo = [ArtifactInfo('SynDNA hits', 'BIOM', [
            (fp_biom, 'biom'), (fp_alng, 'log'),
            (lin_regress_results_fp, 'log'),
            (fit_syndna_models_log_fp, 'log')])]
    else:
        ainfo = []
        errors.append('Missing files from the "SynDNA hits"; please '
                      'contact qiita.help@gmail.com for more information')

    fp_seqs = f'{out_dir}/reads'
    reads = []
    regex_fwd = re.compile(r'(.*)_(S\d{1,4})_(L\d{1,3})_([RI]1).*')
    regex_rev = re.compile(r'(.*)_(S\d{1,4})_(L\d{1,3})_([RI]2).*')
    for f in glob(f'{fp_seqs}/*.fastq.gz'):
        fwd = regex_fwd.match(f)
        rev = regex_rev.match(f)
        if fwd == rev:
            errors.append(f"Is this fwd or rev read? {basename(f)}; please "
                          'contact qiita.help@gmail.com for more information')
            # resetting ainfo
            ainfo = []
        elif fwd is not None:
            reads.append((f, 'raw_forward_seqs'))
        else:
            reads.append((f, 'raw_reverse_seqs'))

    if not errors:
        ainfo.append(
            ArtifactInfo('reads without SynDNA', 'per_sample_FASTQ', reads))

    if errors:
        return False, ainfo, '\n'.join(errors)
    else:

        return True, ainfo, ""


def calculate_cell_counts(qclient, job_id, parameters, out_dir):
    """Run calc_ogu_cell_counts_per_g_of_sample_for_qiita

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to wolka syndna
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """
    error = ''
    # let's get the syndna_id and prep in a single go
    syndna_id = parameters['SynDNA hits']
    syndna_files, prep = qclient.artifact_and_preparation_files(syndna_id)
    if 'log' not in syndna_files.keys():
        error = ("No logs found, are you sure you selected the correct "
                 "artifact for 'SynDNA hits'?")
    else:

        lin_regress_by_sample_id_fp = [f for f in syndna_files['log']
                                       if 'lin_regress_by_sample_id' in f]
        if not lin_regress_by_sample_id_fp:
            error = ("No 'lin_regress_by_sample_id' log found, are you sure "
                     " you selected the correct artifact for 'SynDNA hits'?")
        else:
            lin_regress_by_sample_id_fp = lin_regress_by_sample_id_fp[0]

            # for per_genome_id let's do it separately so we can also ge the
            # sample information
            per_genome_id = parameters['Woltka per-genome']
            ainfo = qclient.get("/qiita_db/artifacts/%s/" % per_genome_id)
            aparams = ainfo['processing_parameters']
            ogu_fp = ainfo['files']['biom'][0]['filepath']

            if 'Database' not in aparams or not ogu_fp.endswith('none.biom'):
                error = ("The selected 'Woltka per-genome' artifact doesn't "
                         "look like one, did you select the correct file?")
            elif 'plain_text' not in ainfo['files']:
                error = ("'Woltka per-genome' artifact is missing "
                         'coverage information')
            else:
                coverages_fp = ainfo['files']['plain_text'][0]['filepath']
                with topen(coverages_fp, 'r:gz') as tgz:
                    member = tgz.getmember('coverage_percentage.txt')
                    coverages = tgz.extractfile(member)
                    coverages_df = pd.read_csv(
                        coverages, sep='\t', header=None)
                # this is the micov generated format; the legacy format only
                # has 2 columns
                if coverages_df.shape[1] == 4:
                    coverages_df = coverages_df.drop(columns=[1, 2])[1:]
                coverages_df.columns = [OGU_ID_KEY, OGU_PERCENT_COVERAGE_KEY]
                coverages_df[OGU_PERCENT_COVERAGE_KEY] = coverages_df[
                    OGU_PERCENT_COVERAGE_KEY].astype(float)

                ogu_counts_per_sample = load_table(ogu_fp)

                db_files = _process_database_files(aparams['Database'])
                ogu_lengths_fp = db_files["length.map"]

                sample_info = qclient.get(
                    '/qiita_db/prep_template/%s/data/?sample_information=true'
                    % ainfo['prep_information'][0])
                sample_info = pd.DataFrame.from_dict(
                    sample_info['data'], orient='index')
                sample_info.reset_index(names='sample_name', inplace=True)

    if error:
        return False, None, error

    try:
        output = calc_ogu_cell_counts_per_g_of_sample_for_qiita(
            sample_info, prep, lin_regress_by_sample_id_fp,
            ogu_counts_per_sample, coverages_df, ogu_lengths_fp,
            min_coverage=float(parameters['min_coverage']),
            min_rsquared=float(parameters['min_rsquared']))
    except Exception as e:
        return False, None, str(e)

    log_fp = f'{out_dir}/cell_counts.log'
    with open(log_fp, 'w') as f:
        f.write(output['calc_cell_counts_log'])
    biom_fp = f'{out_dir}/cell_counts.biom'
    with biom_open(biom_fp, 'w') as f:
        output['cell_count_biom'].to_hdf5(f, f"Cell Counts - {job_id}")
    ainfo = [
        ArtifactInfo(
            'Cell counts', 'BIOM', [(biom_fp, 'biom'), (log_fp, 'log')])]

    return True, ainfo, ""


def calculate_rna_copy_counts(qclient, job_id, parameters, out_dir):
    """Run calc_copies_of_ogu_orf_ssrna_per_g_sample_for_qiita

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to wolka syndna
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """

    per_gene_id = parameters['Woltka per-gene']
    ainfo = qclient.get("/qiita_db/artifacts/%s/" % per_gene_id)
    aparams = ainfo['processing_parameters']
    pg_fp = ainfo['files']['biom'][0]['filepath']

    if 'Database' not in aparams or not pg_fp.endswith('per-gene.biom'):
        error = ("The selected 'Woltka per-gene' artifact doesn't "
                 "look like one, did you select the correct file?")
        return False, None, error

    pergene = load_table(pg_fp)
    db_files = _process_database_files(aparams['Database'])
    ogu_orf_coords_fp = db_files["gene_coordinates"]

    _, prep_info = qclient.artifact_and_preparation_files(per_gene_id)

    sample_info = qclient.get(
        '/qiita_db/prep_template/%s/data/?sample_information=true'
        % ainfo['prep_information'][0])
    sample_info = pd.DataFrame.from_dict(
        sample_info['data'], orient='index')
    sample_info.reset_index(names='sample_name', inplace=True)

    try:
        output, log_msgs = calc_copies_of_ogu_orf_ssrna_per_g_sample_for_qiita(
            sample_info, prep_info, pergene, ogu_orf_coords_fp)
    except Exception as e:
        return False, None, str(e)

    log_fp = f'{out_dir}/rna_copy_counts.log'
    with open(log_fp, 'w') as f:
        f.write(''.join(log_msgs))
    biom_fp = f'{out_dir}/rna_copy_counts.biom'
    with biom_open(biom_fp, 'w') as f:
        output.to_hdf5(f, f"RNA copy counts - {job_id}")
    ainfo = [
        ArtifactInfo(
            'RNA copy counts', 'BIOM', [(biom_fp, 'biom'), (log_fp, 'log')])]

    return True, ainfo, ""
