# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os import environ, listdir
from os.path import join, isdir, expanduser
from configparser import ConfigParser

import gzip
from signal import signal, SIGPIPE, SIG_DFL

from qiita_client import QiitaClient


plugin_details = {'name': 'qp-woltka',
                  'version': '2022.09',
                  'description': 'Woltka'}


def get_dbs(db_folder):
    dbs = dict()
    # Loop through the databases and create a dict of them
    for folder in listdir(db_folder):
        folder_path = join(db_folder, folder)
        if isdir(folder_path):
            files = listdir(folder_path)
            # the bowtie2 db format is name.#.bt2 so getting just the name
            db_fp = [f for f in files if (f.endswith('.bt2') or
                     f.endswith('.bt2l')) and 'rev' not in f][0].rsplit(
                        '.', 2)[0]
            dbs[folder] = join(folder_path, db_fp)

    return (dbs)


def generate_woltka_dflt_params():
    dflt_param_set = {}
    db_parent_path = environ["QC_WOLTKA_DB_DP"]
    # Get a the databases available and the database name
    dbs = get_dbs(db_parent_path)
    # Create dict with command options per database
    for name, fp in dbs.items():
        dflt_param_set[name] = {'Database': fp}

    return (dflt_param_set)


def client_connect(url):
    name = plugin_details['name']
    version = plugin_details['version']

    config = ConfigParser()
    conf_dir = environ.get(
        'QIITA_PLUGINS_DIR', join(expanduser('~'), '.qiita_plugins'))
    conf_fp = join(conf_dir, f'{name}_{version}.conf')

    with open(conf_fp, 'U') as conf_file:
        config.readfp(conf_file)
    qclient = QiitaClient(url, config.get('oauth2', 'CLIENT_ID'),
                          config.get('oauth2', 'CLIENT_SECRET'),
                          server_cert=config.get('oauth2', 'SERVER_CERT'))

    return qclient


def mux(files, output):
    # https://linuxpip.org/broken-pipe-python-error/
    # Ignore SIG_PIPE and don't throw exceptions on it
    # http://docs.python.org/library/signal.html
    signal(SIGPIPE, SIG_DFL)

    delimiter = b'@@@'
    newline = b'\n'

    errors = []
    # the name used here is the filename, it is not read orientation agnostic,
    # should it be?
    for f in files:
        name = f.split('/')[-1]
        name = name.split('.fastq')[0].encode('ascii')

        try:
            fp = gzip.open(f)
            id_ = iter(fp)
            seq = iter(fp)
            dumb = iter(fp)
            qual = iter(fp)
            for i, s, d, q in zip(id_, seq, dumb, qual):
                base_i = i.strip().split(b' ', 1)[0]
                new_i = base_i + delimiter + name + newline

                output.write(new_i)
                output.write(s)
                output.write(d)
                output.write(q)
        except Exception as e:
            errors.append(f'{f}\t{e}')

    if errors:
        with open('errors.log', 'w') as error_log:
            error_log.write('\n'.join(errors))


def search_by_filename(fname, lookup):
    if fname in lookup:
        return lookup[fname]

    original = fname
    while '_' in fname:
        fname = fname.rsplit('_', 1)[0]
        if fname in lookup:
            return lookup[fname]

    fname = original
    while '.' in fname:
        fname = fname.rsplit('.', 1)[0]
        if fname in lookup:
            return lookup[fname]

    raise KeyError("Cannot determine run_prefix for %s" % original)


def demux(input_, output, prep):
    lookup = prep.set_index('run_prefix')['sample_name'].to_dict()
    delimiter = b'@@@'
    tab = b'\t'
    mode = 'ab'  # IMPORTANT: we are opening in append not write
    ext = b'.sam'
    sep = b'/'

    # read each record
    # parse the filename out
    # if the file is not open, open it
    # we are assuming the records are coming in grouped by file
    # however the method will work if records are for whatever
    # reason interleaved
    current_fname = None
    current_fp = None

    for line in input_:
        id_, remainder = line.split(tab, 1)
        id_, fname = id_.rsplit(delimiter, 1)
        fname = fname.strip()

        if fname != current_fname:
            if current_fp is not None:
                current_fp.close()

            sample = search_by_filename(
                fname.decode('utf8'), lookup).encode('ascii')
            fullname = output + sep + sample + ext
            current_fp = open(fullname, mode)
            current_fname = fname
            print(fullname.decode('ascii'))

        current_fp.write(id_)
        current_fp.write(tab)
        current_fp.write(remainder)
