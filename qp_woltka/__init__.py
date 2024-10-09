# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .woltka import (woltka, woltka_syndna, calculate_cell_counts,
                     calculate_rna_copy_counts)
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
    'Per genome Predictions': 'BIOM',
    'Per gene Predictions': 'BIOM',
    # functional
    'KEGG Ontology (KO)': 'BIOM',
    'KEGG Enzyme (EC)': 'BIOM',
    'KEGG Pathway': 'BIOM',
    }
dflt_param_set = generate_woltka_dflt_params()
woltka_cmd = QiitaCommand(
    'Woltka v0.1.7, paired-end',
    "Functional and Taxonomic Predictions", woltka,
    req_params, opt_params, outputs, dflt_param_set)
plugin.register_command(woltka_cmd)

# split synDNA in reads
db_path = environ["QC_WOLTKA_SYNDNA_DB_DP"]
req_params = {
    'input': ('artifact', ['per_sample_FASTQ'])
}
opt_params = {
    'Database': [f'choice: ["{db_path}"]', db_path],
    'min_sample_counts': ('integer', '1')
}
outputs = {
    'SynDNA hits': 'BIOM',
    'reads without SynDNA': 'per_sample_FASTQ',
}
dflt_param_set = {
    'SynDNA': {'Database': db_path, 'min_sample_counts': 1},
}
syndna_cmd = QiitaCommand(
    'Remove SynDNA reads', "Remove SynDNA reads using bowtie2 & woltka, and "
    "generate even mates", woltka_syndna,
    req_params, opt_params, outputs, dflt_param_set)
plugin.register_command(syndna_cmd)

# WGS cell counts
req_params = {
    'SynDNA hits': ('artifact', ['BIOM']),
    'Woltka per-genome': ('artifact', ['BIOM'])
}
opt_params = {
    'min_coverage': ('integer', '1'),
    'read_length': ('integer', '150'),
    'min_rsquared': ('float', '0.8'),
}
outputs = {
    'Cell counts': 'BIOM'
}
dflt_param_set = {
    '150bp @ min_coverage:1 R^2:0.8': {
        'min_coverage': 1,
        'read_length': 150,
        'min_rsquared': 0.8
    }
}
calculate_cell_counts_cmd = QiitaCommand(
    'Calculate Cell Counts', "Calculate cell counts per-genome",
    calculate_cell_counts, req_params, opt_params, outputs, dflt_param_set)
plugin.register_command(calculate_cell_counts_cmd)


# MTX calculate RNA copy counts
req_params = {
    'Woltka per-gene': ('artifact', ['BIOM'])
}
opt_params = {}
outputs = {
    'RNA copy counts': 'BIOM'
}
dflt_param_set = {'default': {}}
calculate_rna_copy_counts_cmd = QiitaCommand(
    'Calculate RNA Copy Counts', "Calculate RNA copy counts per-gene",
    calculate_rna_copy_counts, req_params, opt_params, outputs, dflt_param_set)
plugin.register_command(calculate_rna_copy_counts_cmd)
