# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .woltka import woltka
from qp_woltka.util import generate_woltka_dflt_params, get_dbs, plugin_details
from os import environ

# Initialize the plugin
plugin = QiitaPlugin(**plugin_details)

db_list = list(get_dbs(environ["QC_WOLTKA_DB_DP"]).values())
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # database
    'Database': ["choice: [%s]" % ','.join([f'"{db}"' for db in db_list]),
                 # making the first option default and rm quotes
                 db_list[0]],
    }
outputs = {
    'Alignment Profile': 'BIOM',
    'Taxonomic Predictions - phylum': 'BIOM',
    'Taxonomic Predictions - genus': 'BIOM',
    'Taxonomic Predictions - species': 'BIOM',
    'Per genome Predictions': 'BIOM',
    'Per gene Predictions': 'BIOM',
    }
dflt_param_set = generate_woltka_dflt_params()

woltka_cmd = QiitaCommand(
    'Woltka v0.1.1', "Functional and Taxonomic Predictions", woltka,
    req_params, opt_params, outputs, dflt_param_set)

plugin.register_command(woltka_cmd)
