#!/usr/bin/env python

import click
import glob as glob_
import re
import os
from math import ceil
import subprocess


@click.group()
def cli():
    pass


# https://stackoverflow.com/a/40195800
_common_options = [
    click.option('--directory', type=click.Path(exists=True),
                 help='The path to obtain files from', required=False),
    click.option('--files', multiple=True, type=str, required=False,
                 help='Specific entries to use as input files, assumes '
                 'absolute paths, incompatible with --directory'),
    click.option('--output', type=click.Path(exists=False),
                 help='The path to write output to', required=True),
    click.option('--output-extension', type=str,
                 help='An extension to add to output files', required=True),
    click.option('--glob', type=str,
                 help='A file matching pattern (e.g., "*.fastq.gz"). '
                      'Defaults to "*"',
                 required=False, default='*'),
    click.option('--max-running', type=int,
                 help='The maximum number of jobs to run at any given time',
                 required=True),
    click.option('--ppn', type=int,
                 help='The number of processors to use per task',
                 required=True),
    click.option('--walltime', type=str,
                 help='The maximum walltime for a given task',
                 required=True),
    click.option('--memory', type=str, help='The amount of memory per task '
                                            '(e.g., 8g)',
                 required=True),
    click.option('--environment', type=str,
                 help='A conda environment to activate', required=True),
    click.option('--queue', type=str, help='The submission queue',
                 required=False),
    click.option('--name', type=str, help='The submission name', required=True)
]


def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options


@cli.command()
@click.option('--submit', is_flag=True, default=False,
              help="Actually submit the jobs")
@click.option('--preparation-information', type=click.Path(exists=True),
              required=True, help="The preparation file which describes the "
                                  "relationship between a run_prefix and the "
                                  "sample_name")
@click.option('--database-bowtie2', type=str, required=True,
              help="The database to use")
@click.option('--database-gene-coordinates', type=str,
              required=False,
              help="The gene coordinates relative to the bowtie2 database")
@click.option('--database-taxonomy', type=str, required=True,
              help="The taxonomy relative to the bowtie2 database")
@add_options(_common_options)
@click.pass_context
def woltka_from_per_sample(ctx, submit, preparation_information,
                           database_bowtie2, database_gene_coordinates,
                           database_taxonomy, **kwargs):
    """Compute Woltka genome mapping per sample"""
    if kwargs.get('directory') is None:
        raise ValueError("woltka-from-per-sample requires --directory")

    name = kwargs['name']
    ppn = kwargs['ppn']
    output = kwargs['output']
    queue = kwargs['queue']
    environment = kwargs['environment']
    directory = kwargs.pop('directory')

    # determine per-sample prefixes. force *.fastq.gz here for the glob
    # as that's what the per-sample shotgun data are...
    splitter = re.compile(r'[._]R[12][._]')
    basenames = set()
    for f in glob_.glob(os.path.join(directory, '*.fastq.gz')):
        sample = os.path.basename(f)
        hit = splitter.search(sample)
        if hit is None:
            raise ValueError('%s appears malformed' % sample)
        else:
            prefix = sample[:hit.span()[0]]
            basenames.add(os.path.join(directory, prefix))
    kwargs['files'] = list(basenames)

    # woltka assumes R1 and R2 are combined even though it doesn't use the
    # paired end data, so let's concatenate based off the prefixes first.
    # note: 'cat' is safe with gzip'd data, see:
    # https://stackoverflow.com/a/8005155
    concat = 'cat {infile}*.fastq.gz > {outfile}.fastq.gz'

    # Bowtie2 command structure based on shogun settings
    # https://github.com/knights-lab/SHOGUN/blob/ff1aabe772469d6a1c2c83cf146140b5341df83c/shogun/wrappers/bowtie2_wrapper.py#L21-L37
    # And as described in:
    # https://github.com/BenLangmead/bowtie2/issues/311
    bowtie2 = f'bowtie2 -p {ppn} -x {database_bowtie2} ' + \
              '-q {outfile}.fastq.gz -S {outfile}.sam --seed 42 ' + \
              '--very-sensitive -k 16 --np 1 --mp "1,1" ' + \
              '--rdg "0,1" --rfg "0,1" --score-min ' + \
              '"L,0,-0.05" --no-head --no-unal'
    xz = f'xz -9 -T{ppn} -c ' + '{outfile}.sam > {outfile}.xz'

    # Not performing demux as this is per sample, so no need
    ranks = ["phylum", "genus", "species", "free", "none"]
    woltka = 'woltka classify -i {outfile}.sam ' + \
             '-o {outfile}.woltka-taxa ' + \
             '--no-demux ' + \
             f'--lineage {database_taxonomy} ' + \
             f'--rank {",".join(ranks)}'

    # compute per-gene results
    if database_gene_coordinates is not None:
        woltka_per_gene = 'woltka classify -i {outfile}.sam ' + \
                          f'-c {database_gene_coordinates} ' + \
                          '-o {outfile}.woltka-per-gene ' + \
                          '--no-demux'
        cmd_fmt = f'{concat}; {bowtie2}; {woltka}; {woltka_per_gene}; {xz}'
    else:
        cmd_fmt = f'{concat}; {bowtie2}; {woltka}; {xz}'

    # first we'll use the concatenation command, then run bowtie2,
    # finally we'll run woltka
    kwargs['command_format'] = cmd_fmt

    # now, let's establish the merge script.
    merges = []
    merge_inv = f'woltka_merge --prep {preparation_information} ' + \
                f'--base {output} '
    for r in ranks:
        merges.append(" ".join([merge_inv,
                                f'--name {r}',
                                f'--glob "*.woltka-taxa/{r}.biom"',
                                '&']))  # run all at once
    if database_gene_coordinates is not None:
        merges.append(" ".join([merge_inv, '--name per-gene',
                                '--glob "*.woltka-per-gene"',
                                '--rename &']))  # run all at once
    else:
        m = merges[-1].split(' ')
        merges[-1] = " ".join(m[:-1] + ['--rename'] + [m[-1]])
    merges.append(f'cd {output}; tar -cvf alignment.tar *.sam.xz &')

    # The merge for a HiSeq 2000 lane was 40 seconds and ~150MB of memory.
    # But, let's over request just in case (and this is a very small request
    # relative to the rest of the work).
    n_merges = len(merges)
    assert n_merges < 32  # 32 merges would be crazy...

    lines = ['#!/bin/bash',
             '#PBS -M qiita.help@gmail.com',
             f'#PBS -N merge-{name}',
             f'#PBS -l nodes=1:ppn={n_merges}',
             '#PBS -l walltime=4:00:00',
             '#PBS -l mem=48g',
             f'#PBS -o {output}/merge-{name}.log',
             f'#PBS -e {output}/merge-{name}.err',
             f'#PBS -q {queue}' if queue is not None else '',
             f'cd {output}',
             f'source activate {environment}',
             'date',  # start time
             'hostname',  # executing system
             'set -e',
             '\n'.join(merges),
             "wait",
             "date"]  # end time

    # construct the job array
    ctx.invoke(to_array, **kwargs)

    # write out the merge script
    merge_name = f'{name}.merge.qsub'
    with open(os.path.join(output, merge_name), 'w') as out:
        out.write('\n'.join(lines))
        out.write('\n')

    if submit:
        bowtie = subprocess.run(['qsub', os.path.join(output, f'{name}.qsub')],
                                stdout=subprocess.PIPE)
        job = bowtie.stdout.decode('utf8')
        subprocess.run(['qsub', f'-W depend=afterokarray:{job}',
                        os.path.join(output, merge_name)])


@cli.command()
@click.option('--command-format', type=str,
              help='A structured command format string. Arguments are '
                   'specified using f-string format. Only {infile} and '
                   '{outfile} are allowed', required=True)
@add_options(_common_options)
def to_array(directory, output, glob, max_running, ppn, walltime, queue,
             command_format, memory, name, output_extension, environment,
             files):
    if directory is None and len(files) == 0:
        raise ValueError('Must have a --directory or --files')
    elif directory is not None and len(files) > 0:
        raise ValueError('Cannot have both --directory and --files')
    elif directory is not None:
        files = [f for f in glob_.glob(os.path.join(directory, glob))]

    # sanity checking
    assert len(files) > 0
    assert ppn > 0
    assert re.match(r'\d+:\d\d:\d\d', walltime) is not None
    assert '{infile}' in command_format
    assert '{outfile}' in command_format

    # 1024 -> maximum number of job array IDs
    max_jobs = 1024
    n_files = len(files)

    # if we have too many files, break pack them into the individual
    # jobs. So, if we had 2000 files, we would create 1000 jobs of
    # each which process 2 files. If we had 1500, we would create
    # 750 jobs each processing 2 files.
    if n_files > max_jobs:
        per_job = int(ceil(n_files / max_jobs))
        n_jobs = int(ceil(n_files / per_job))
    else:
        per_job = 1
        n_jobs = n_files

    os.makedirs(output)

    # the details describe each input file, and its output file
    details_name = os.path.join(output, f'{name}.array-details')
    with open(details_name, 'w') as details:
        for f in files:
            basename = os.path.basename(f)
            details.write(f'{f}\t{output}/{basename}.{output_extension}\n')

    # all the setup pieces
    lines = ['#!/bin/bash',
             '#PBS -M qiita.help@gmail.com',
             f'#PBS -N {name}',
             f'#PBS -l nodes=1:ppn={ppn}',
             f'#PBS -l walltime={walltime}',
             f'#PBS -l mem={memory}',
             f'#PBS -o {output}/{name}' + '_${PBS_ARRAYID}.log',
             f'#PBS -e {output}/{name}' + '_${PBS_ARRAYID}.err',
             f'#PBS -q {queue}' if queue is not None else '',
             f'#PBS -t 1-{n_jobs}%{max_running}',
             f'cd {output}',
             f'source activate {environment}',
             'date',  # start time
             'hostname',  # executing system
             'offset=${PBS_ARRAYID}']

    # if we have more than one file per job, we need to adjust our offset
    # position accordingly. If we had three files per job, then the first
    # job processes lines 1, 2, and 3 of the details. The second job
    # processes lines 4, 5, 6. Note that the PBS_ARRAYID is 1-based not
    # 0-based.
    if per_job > 1:
        lines.append(f"offset=$(( $offset * {per_job} ))")

    # reversed due to the substraction with the offset, so that we process
    # lines N, N+1, N+2, etc in the details file.
    for i in reversed(range(per_job)):
        lines.append(f'step=$(( $offset - {i} ))')

        # do not let the last job in the array overstep
        lines.append(f'if [[ $step -gt {n_files} ]]; then exit 0; fi')

        # if we're okay, get the next set of arguments
        lines.append(f'args{i}=$(head -n $step {details_name} | tail -n 1)')

        # f-string is broken up so the awk program is not interpreted
        # as a component of the f-string
        lines.append(
            f"infile{i}=$(echo -e $args{i}" + " | awk '{ print $1 }')")
        lines.append(
            f"outfile{i}=$(echo -e $args{i}" + " | awk '{ print $2 }')")

        # wrap the command calls in "fail on error", and then disable it. The
        # reason to enable (-e) and disable (+e) is some of the other shell
        # scripting may *correctly* produce a nonzero exit status, like the
        # calls to the [ program.
        lines.append('set -e')
        cmd_args = {'infile': f"$infile{i}", 'outfile': f"$outfile{i}"}
        lines.append(command_format.format(**cmd_args))
        lines.append('set +e')
    lines.append('date')  # end time

    # write out the script
    with open(os.path.join(output, f'{name}.qsub'), 'w') as job:
        job.write('\n'.join(lines))
        job.write('\n')


if __name__ == '__main__':
    cli()
