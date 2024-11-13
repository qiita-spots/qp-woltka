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

from qiita_client import QiitaClient


plugin_details = {'name': 'qp-woltka',
                  'version': '2024.09',
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
                          config.get('oauth2', 'SERVER_CERT'))

    return qclient


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

    for rp in lookup:
        if original.startswith(rp):
            return lookup[rp]

    raise KeyError("Cannot determine run_prefix for %s" % original)
