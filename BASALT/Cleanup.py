#!/usr/bin/env python

"""Conservative cleanup helpers for BASALT workspaces."""

import fnmatch
import glob
import os
import shutil


def cleanup_enabled():
    value = str(os.environ.get('BASALT_CLEANUP_ENABLED', '1')).strip().lower()
    return value not in ('0', 'false', 'no', 'off', '')


def _remove_path(path):
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)


def _is_preserved_pipeline_state(path):
    """Return whether a path is required for downstream work or resuming."""
    name = os.path.basename(os.path.normpath(path))
    preserved_patterns = (
        '*_genomes',
        '*_BestBinsSet',
        '*_BestBinset',
        '*_comparison_files',
        '*_checkm',
        '*_checkm2',
        'BestBinset*',
        'Final_binset*',
        'Final_bestbinset*',
        'Data_feeded*',
        'Coverage_matrix_*',
        '*.depth.txt',
        'Connections_*',
        'Combat_*',
        'condense_connections_*.txt',
        '*checkpoint*.txt',
        'BASALT_command.txt',
        'Basalt_log.txt',
        'BASALT_log.txt',
        '*_list.txt',
        'PE_r1_*',
        'PE_r2_*',
        '*quality_report*.tsv',
        '*_quality_report.tsv',
        '*_contigs_summary.txt',
        'Bins_change_ID_*.txt',
        'Bins_total_connections_*.txt',
        'Genome_*.txt',
        '*_deep_retrieval*',
        '*_MP_1',
        '*_MP_2',
        '*_gf_lr*',
        '*_long_read',
        '*_sr_bins_seq',
        '*_lr_bins_seq',
        'Polish_*',
        'Remained_seq*',
        'Remapping*',
        'Re-mapped_depth.txt',
        'Total_contigs_after_OLC_reassembly.fa',
        'Hybrid_re-assembly_status.txt',
        'Reassembled_bins_comparison.txt',
    )
    return any(fnmatch.fnmatch(name, pattern) for pattern in preserved_patterns)


def _remove_patterns(base_dir, patterns, preserve_pipeline_state=False):
    for pattern in patterns:
        for path in glob.glob(os.path.join(base_dir, pattern)):
            if preserve_pipeline_state and _is_preserved_pipeline_state(path):
                continue
            _remove_path(path)


def cleanup_checkm2_output(checkm_dir):
    """Remove bulky CheckM2 scratch folders only after its report exists."""
    if not cleanup_enabled() or not os.path.isdir(checkm_dir):
        return
    if not os.path.isfile(os.path.join(checkm_dir, 'quality_report.tsv')):
        return
    _remove_patterns(
        checkm_dir,
        ['diamond_output', 'protein_files', 'genes', 'tmp', 'checkm2_res'],
    )


def cleanup(assembly_list):
    """Remove disposable scratch while retaining all resume-critical state."""
    del assembly_list  # Kept in the public API for compatibility.
    if not cleanup_enabled():
        return
    scratch_patterns = [
        '*.sam', '*.bam', '*.bai', '*.bt2', '*.bt2l', '*.njs', '*.ndb',
        '*.nto', '*.ntf', '*.not', '*.nos', '*.nhr', '*.nin', '*.nsq',
        '*.nog', '*.nsd', '*.nsi', '*.nhd', '*.nhi', 'temp.orfs.*',
        'temp_db.txt', '*.tmp', '*_kmer', 'bin_coverage',
        'Bin_coverage_after_contamination_removal', 'bin_comparison_folder',
        'bin_extract-eleminated-selected_contig', 'Bins_blast_output',
        'split_blast_output', 'SPAdes_corrected_reads', 'Concoct_*',
    ]
    _remove_patterns(os.getcwd(), scratch_patterns, preserve_pipeline_state=True)


if __name__ == '__main__':
    cleanup([])
