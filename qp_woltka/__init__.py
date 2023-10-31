# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .woltka import woltka, syndna_woltka
from qp_woltka.util import generate_woltka_dflt_params, get_dbs, plugin_details
from os import environ

# Initialize the plugin
plugin = QiitaPlugin(**plugin_details)

# Main woltka command
db_list = list(get_dbs(environ["QC_WOLTKA_DB_DP"]).values())
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    # database
    'Database': ["choice: [%s]" % ','.join([f'"{db}"' for db in db_list]),
                 # making the first option default and rm quotes
                 db_list[0]],
    }
outputs = {
    # taxonomic
    'Alignment Profile': 'BIOM',
    'Per genome Predictions': 'BIOM',
    'Per gene Predictions': 'BIOM',
    # functional
    'KEGG Ontology (KO)': 'BIOM',
    'KEGG Enzyme (EC)': 'BIOM',
    'KEGG Pathway': 'BIOM',
    }
dflt_param_set = generate_woltka_dflt_params()
woltka_cmd = QiitaCommand(
    'Woltka v0.1.4', "Functional and Taxonomic Predictions", woltka,
    req_params, opt_params, outputs, dflt_param_set)
plugin.register_command(woltka_cmd)

# abs quant
db_path = environ["QC_WOLTKA_SYNDNA_DB_DP"]
req_params = {
    'input': ('artifact', ['per_sample_FASTQ'])
}
opt_params = {
    'Database': [f"choice: [{db_path}]", db_path]
}
outputs = {
    'SynDNA hits': 'BIOM',
    'reads without SynDNA': 'per_sample_FASTQ',
}
dflt_param_set = {
    'SynDNA': {'Database': db_path},
}
syndna_cmd = QiitaCommand(
    'SynDNA Woltka v0.1', "Process SynDNA reads using woltka", syndna_woltka,
    req_params, opt_params, outputs, dflt_param_set)
plugin.register_command(syndna_cmd)
