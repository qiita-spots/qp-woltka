# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import re
from os import environ, mkdir, listdir
from os.path import join, basename, exists, dirname
from glob import glob
from shutil import copy2
from math import ceil
from biom import load_table
import pandas as pd
from src.fit_syndna_models import fit_linear_regression_models_for_qiita

from qp_woltka.util import search_by_filename

from qiita_client import ArtifactInfo

# resources per job
PPN = 8
MAX_RUNNING = 8
TASKS_IN_SCRIPT = 10

MEMORY = '90g'
LARGE_MEMORY = '150g'
MERGE_MEMORY = '140g'
SYNDNA_MEMORY = '190g'

WALLTIME = '40:00:00'
MERGE_WALLTIME = '30:00:00'
SYNDNA_WALLTIME = '8:00:00'


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

    db_files = _process_database_files(database_bowtie2)
    db_folder = dirname(database_bowtie2)
    db_name = basename(database_bowtie2)

    n_files = 1
    for i, (k, (f, r)) in enumerate(files.items()):
        if i >= n_files*TASKS_IN_SCRIPT:
            n_files += 1
        with open(join(output, f'sample_details_{n_files}.txt'), 'a+') as fh:
            fh.write(f'{f["filepath"]}\n')
            if r is not None:
                fh.write(f'{r["filepath"]}\n')

    ranks = ["free", "none"]
    # now, let's establish the merge script.
    merges = []
    merge_inv = f'woltka_merge --base {output} '
    fcmds = []
    for r in ranks:
        cmd = [merge_inv, f'--name {r}', f'--glob "*.woltka-taxa/{r}.biom"']
        if r == 'free' and 'length.map' in db_files:
            cmd.append(f'--length_map {db_files["length.map"]}')
        cmd.append('&')
        merges.append(" ".join(cmd))
    if db_files['gene_coordinates'] is not None:
        merges.append(" ".join([merge_inv, '--name per-gene',
                                '--glob "*.woltka-per-gene"',
                                '--rename &']))  # run all at once
        wcdm = 'woltka tools collapse -i '
        dbfk = db_files['kegg']
        if dbfk["orf-to-ko.map.xz"] is not None:
            fcmds.append(f'{wcdm} per-gene.biom -m {dbfk["orf-to-ko.map.xz"]} '
                         '-o ko.biom')
        if dbfk["ko-to-ec.map"] is not None:
            fcmds.append(f'{wcdm} ko.biom -m {dbfk["ko-to-ec.map"]} '
                         '-o ec.biom')
        if dbfk["ko-to-reaction.map"] is not None and \
                dbfk["reaction-to-module.map"] is not None and \
                dbfk["module-to-pathway.map"] is not None:
            fcmds.append(f'{wcdm} ko.biom -m {dbfk["ko-to-reaction.map"]} '
                         '-o reaction.biom')
            fcmds.append(f'{wcdm} reaction.biom -m '
                         f'{dbfk["reaction-to-module.map"]} -o module.biom')
            fcmds.append(f'{wcdm} module.biom -m '
                         f'{dbfk["module-to-pathway.map"]} -o pathway.biom')
    else:
        # for "simplicity" we will inject the `--rename` flag to the last
        # merge command (between all the parameters and the last &)
        m = merges[-1].split(' ')
        merges[-1] = " ".join(m[:-1] + ['--rename'] + [m[-1]])

    # The merge for a HiSeq 2000 lane was 40 seconds and ~150MB of memory.
    # But, let's over request just in case (and this is a very small request
    # relative to the rest of the work).
    n_merges = len(merges)
    assert n_merges < 32  # 32 merges would be crazy...

    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             f'#SBATCH --job-name merge-{name}',
             '#SBATCH -N 1',
             f'#SBATCH -n {n_merges}',
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
             'sjobs=`ls sample_details_* | wc -l`',
             'if [[ ! -f "errors.log" && $sruns -eq $sjobs ]]; then',
             '\n'.join(merges),
             "wait",
             '\n'.join(fcmds),
             f'cd {output}; tar -cvf alignment.tar *.sam.xz; '
             'tar zcvf coverages.tgz coverage_percentage.txt artifact.cov '
             'coverages\n'
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
    bowtie2 = 'mux ${files} | ' + \
              f'bowtie2 -p {PPN} -x {database_bowtie2} ' + \
              '-q - --seed 42 ' + \
              '--very-sensitive -k 16 --np 1 --mp "1,1" ' + \
              '--rdg "0,1" --rfg "0,1" --score-min ' + \
              '"L,0,-0.05" --no-head --no-unal' + \
              " | cut -f1-9 | sed 's/$/\t*\t*/'" + \
              ' | demux ${output} ' + preparation_information + \
              ' | sort | uniq > sample_processing_${SLURM_ARRAY_TASK_ID}.log'
    woltka = 'woltka classify -i ${f} ' + \
             '-o ${f}.woltka-taxa ' + \
             '--no-demux ' + \
             f'--lineage {db_files["taxonomy"]} ' + \
             f'--rank {",".join(ranks)} --outcov coverages/'

    memory = MEMORY
    if 'RS210' in database_bowtie2:
        memory = LARGE_MEMORY

    # all the setup pieces
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             f'#SBATCH --job-name {name}',
             '#SBATCH -N 1',
             f'#SBATCH -n {PPN}',
             f'#SBATCH --time {WALLTIME}',
             f'#SBATCH --mem {memory}',
             f'#SBATCH --output {output}/{name}_%a.log',
             f'#SBATCH --error {output}/{name}_%a.err',
             f'#SBATCH --array 1-{n_files}%{MAX_RUNNING}',
             f'cd {output}',
             f'prep_full_path={preparation_information}',
             f'{environment}',
             'date',  # start time
             'hostname',  # executing system
             'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}',
             f'dbbase={db_folder}',
             f'dbname={db_name}',
             f'output={output}',
             'files=`cat sample_details_${SLURM_ARRAY_TASK_ID}.txt`',
             bowtie2,
             '# for each one of our input files, form woltka commands, ',
             '# and farm off to gnu parallel',
             'for f in `cat sample_processing_${SLURM_ARRAY_TASK_ID}.log`',
             'do',
             f'  echo "{woltka}"']

    if db_files['gene_coordinates'] is not None:
        lines.append('  echo "woltka classify -i ${f} '
                     f'-c {db_files["gene_coordinates"]} '
                     '-o ${f}.woltka-per-gene --no-demux"')
    lines.append('done | parallel -j 8')

    # finally, compress each one of our sam files
    lines.extend([
        'for f in `cat sample_processing_${SLURM_ARRAY_TASK_ID}.log`',
        'do',
        '  # compress it',
        '  echo "xz -1 -T1 ${f}"',
        'done | parallel -j 8'])

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
    fp_biom = f'{out_dir}/free.biom'
    fp_alng = f'{out_dir}/alignment.tar'
    if exists(fp_biom) and exists(fp_alng):
        ainfo = [ArtifactInfo('Alignment Profile', 'BIOM', [
            (fp_biom, 'biom'), (fp_alng, 'log'),
            (_coverage_copy(f'{out_dir}/alignment/'), 'plain_text')])]
    else:
        ainfo = []
        errors.append('Missing files from the "Alignment Profile"; please '
                      'contact qiita.help@gmail.com for more information')

    fp_biom = f'{out_dir}/none.biom'
    if exists(fp_biom):
        ainfo.append(ArtifactInfo('Per genome Predictions', 'BIOM', [
            (fp_biom, 'biom'),
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
        if r is not None:
            line += f'rev_{sname} {r["filepath"]}\n'

        with open(join(output, f'sample_details_{n_files}.txt'), 'a+') as fh:
            fh.write(line)

    # Bowtie2 command structure based on
    # https://github.com/BenLangmead/bowtie2/issues/311
    bowtie2 = f'bowtie2 -p {PPN} -x {database_bowtie2} ' + \
              '-q ${f} -S $PWD/sams/${sn}.sam ' +\
              '--seed 42 --very-sensitive -k 16 --np 1 --mp "1,1" ' + \
              '--rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" ' + \
              '--no-head --no-unal --un-gz $PWD/reads/${fn}'

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
             'mkdir -p reads sams',
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
             f'    {bowtie2}',
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
             'if [[ $sruns -eq $sjobs ]]; then',
             '  mkdir -p sams/final',
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
            f'{out_dir}/prep_info.tsv', index_col='sample_name', sep='\t')
        output = fit_linear_regression_models_for_qiita(
            prep, load_table(fp_biom), parameters['min_sample_counts'])
        # saving results to disk
        lin_regress_results_fp = f'{out_dir}/lin_regress_results.json'
        fit_syndna_models_log_fp = f'{out_dir}/fit_syndna_models_log.txt'
        with open(lin_regress_results_fp, 'w') as fp:
            fp.write(output['LIN_REGRESS_RESULT_KEY'])
        with open(fit_syndna_models_log_fp, 'w') as fp:
            fp.write(output['FIT_SYNDNA_MODELS_LOG_KEY'])
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
    for f in listdir(fp_seqs):
        fwd = regex_fwd.match(f)
        rev = regex_rev.match(f)
        if fwd == rev:
            errors.append(f"Is this fwd or rev read? {basename(f)}; please "
                          'contact qiita.help@gmail.com for more information')
            # resetting ainfo
            ainfo = []
        elif fwd is not None:
            reads.append((f'{fp_seqs}/{f}', 'raw_forward_seqs'))
        else:
            reads.append((f'{fp_seqs}/{f}', 'raw_reverse_seqs'))

    if not errors:
        ainfo.append(
            ArtifactInfo('reads without SynDNA', 'per_sample_FASTQ', reads))

    if errors:
        return False, ainfo, '\n'.join(errors)
    else:

        return True, ainfo, ""


def calculate_cell_counts(qclient, job_id, parameters, out_dir):
    """
    """
    raise ValueError('Not implemented')
