#!/usr/bin/env python

"""
Helper utilities for cleaning up intermediate BASALT files.
"""

import glob
import fnmatch
import os
import shutil


def cleanup_enabled():
    """
    Return whether cleanup behavior is enabled for this BASALT run.
    """
    value = str(os.environ.get('BASALT_CLEANUP_ENABLED', '1')).strip().lower()
    return value not in ('0', 'false', 'no', 'off')


def _remove_path(path):
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path, ignore_errors=True)
    elif os.path.exists(path):
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def _is_preserved_pipeline_state(path):
    """
    Return True for artifacts that a later BASALT stage, resume, or recovery
    may need even if a broad scratch glob happens to match them.
    """
    name = os.path.basename(os.path.normpath(path))
    if not name:
        return False

    preserved_patterns = [
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
        '*_assembly.depth.txt',
        'Connections_*',
        'Connections_total_dict.txt',
        'Combat_*',
        'condense_connections_*.txt',
        'Basalt_checkpoint.txt',
        'Autobinner_checkpoint.txt',
        'De-rep_checkpoint.txt',
        'BASALT_command.txt',
        'Basalt_log.txt',
        '*_list.txt',
        'PE_r1_*',
        'PE_r2_*',
        '*quality_report*.tsv',
        '*_quality_report.tsv',
        '*_contigs_summary.txt',
        'prebinned_genomes_output_for_dataframe_*.txt',
        'Bins_change_ID_*.txt',
        'Bins_total_connections_*.txt',
        'Genome_*.txt',
        'Deep_retrieved_bins*',
        'coverage_deep_refined_bins*',
        'TNFs_deep_refined_bins*',
        'S6_coverage_filtration_matrix*',
        'S6_TNF_filtration_matrix*',
        'Merged_seqs_*',
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
    ]
    return any(fnmatch.fnmatch(name, pattern) for pattern in preserved_patterns)


def _remove_patterns(base_dir, patterns, preserve_pipeline_state=False):
    for pattern in patterns:
        for path in glob.glob(os.path.join(base_dir, pattern)):
            if preserve_pipeline_state and _is_preserved_pipeline_state(path):
                continue
            _remove_path(path)


def cleanup_binner_workspace(bin_dir, assembly_files=None, depth_files=None,
                             extra_patterns=None):
    """
    Remove large temporary files from a single binner output directory.

    Parameters
    ----------
    bin_dir : str
        Binner output directory.
    assembly_files : list[str] | None
        Copied assembly files to delete from the binner directory.
    depth_files : list[str] | None
        Copied depth / coverage files to delete from the binner directory.
    extra_patterns : list[str] | None
        Additional glob patterns to remove from the binner directory.
    """
    if not os.path.isdir(bin_dir):
        return
    if not cleanup_enabled():
        return

    patterns = [
        'Coverage_list*.txt',
        '*.sam',
        '*.bam',
        '*.bt2',
        '*.njs',
        '*.ndb',
        '*.nto',
        '*.ntf',
        '*.not',
        '*.nos',
        '*.seed',
        '*.err',
    ]

    if assembly_files:
        patterns.extend(assembly_files)
    if depth_files:
        patterns.extend(depth_files)
    if extra_patterns:
        patterns.extend(extra_patterns)

    _remove_patterns(bin_dir, patterns)


def cleanup_semibin_workspace(bin_dir):
    """
    Remove SemiBin scratch tables after bins and QC outputs were created.
    """
    cleanup_binner_workspace(
        bin_dir,
        extra_patterns=[
            'data.csv',
            'data_split.csv',
            '*_data_cov.csv',
            '*_data_split_cov.csv',
            'SemiBinRun.log',
            'output_bins',
            'pre_reclustering_bins',
            'recluster_bins',
        ],
    )


def cleanup_checkm2_output(checkm_dir):
    """
    Remove bulky CheckM2 intermediate folders once summary outputs exist.
    """
    if not os.path.isdir(checkm_dir):
        return
    if not cleanup_enabled():
        return

    quality_report = os.path.join(checkm_dir, 'quality_report.tsv')
    if not os.path.exists(quality_report):
        return

    _remove_patterns(
        checkm_dir,
        [
            'diamond_output',
            'protein_files',
            'genes',
            'tmp',
            'checkm2_res',
        ],
    )


def cleanup_autobinner_assembly_workspace(work_dir, group, assembly_name):
    """
    Remove large top-level autobinner intermediates for one finished assembly.

    This is intentionally conservative: it keeps the modified assembly FASTA,
    the merged PE-connection summary and the main depth file because later
    BASALT stages still consume those, but removes the large mapping outputs
    and temporary indices once a per-assembly autobinning run is complete.
    """
    if not os.path.isdir(work_dir):
        return
    if not cleanup_enabled():
        return

    assembly_prefix = f'{group}_{assembly_name}'
    depth_prefix = f'{group}_assembly.depth'
    patterns = [
        f'{group}_DNA-*.sam',
        f'{group}_DNA-*.bam',
        f'{group}_lr*.sam',
        f'{group}_lr*.bam',
        f'{group}_LR-*.sam',
        f'{group}_LR-*.bam',
        f'{assembly_prefix}*.bt2',
        f'{assembly_prefix}.mmi',
        f'condensed.cytoscape.connections_{group}_DNA-*.tab',
        f'Coverage_list_{group}_{assembly_name}.txt',
        f'Concoct_{assembly_prefix}',
        f'Concoct_{depth_prefix}.txt',
        f'{depth_prefix}_1.txt',
        f'{depth_prefix}_2.txt',
    ]
    _remove_patterns(work_dir, patterns)


def cleanup_redundant_short_read_inputs(original_reads, modified_reads):
    """
    Drop original short-read files once PE-tracking-ready copies exist.
    """
    if not cleanup_enabled():
        return
    for original_read, modified_read in zip(original_reads, modified_reads):
        if os.path.exists(modified_read) and os.path.exists(original_read):
            _remove_path(original_read)


def cleanup(assembly_list):
    """
    Remove only disposable scratch files from a completed BASALT run.

    This cleanup is intentionally conservative.  Downstream BASALT stages and
    checkpoint resumes still need bin FASTAs, best-bin folders, coverage
    matrices, comparison files, retrieval/reassembly folders, checkpoints,
    logs, and PE-tracking read copies. Earlier cleanup code archived and
    removed those files, which made later steps unrecoverable without
    rebuilding bins. This function saves space by removing regenerable
    alignment/index/temp outputs while leaving all pipeline state needed for a
    complete run or resume in place.

    Parameters
    ----------
    assembly_list : list
        List of assemblies. Currently unused but kept for API compatibility.
    """
    if not cleanup_enabled():
        return

    scratch_patterns = [
        '*.sam',
        '*.bam',
        '*.bai',
        '*.bt2',
        '*.bt2l',
        '*.njs',
        '*.ndb',
        '*.nto',
        '*.ntf',
        '*.not',
        '*.nos',
        '*.nhr',
        '*.nin',
        '*.nsq',
        '*.nog',
        '*.nsd',
        '*.nsi',
        '*.nhd',
        '*.nhi',
        'temp.orfs.*',
        'temp_db.txt',
        '*.tmp',
    ]
    scratch_dirs = [
        '*_kmer',
        'bin_coverage',
        'Bin_coverage_after_contamination_removal',
        'bin_comparison_folder',
        'bin_extract-eleminated-selected_contig',
        'Bins_blast_output',
        'split_blast_output',
        'SPAdes_corrected_reads',
        'Concoct_*',
    ]

    _remove_patterns(
        os.getcwd(),
        scratch_patterns + scratch_dirs,
        preserve_pipeline_state=True,
    )

    # Keep these required state families uncompressed and in-place:
    # *_genomes, *_BestBinsSet, BestBinset*, *_comparison_files, *_checkm,
    # Coverage_matrix_*, *.depth.txt, Connections_*, Combat_*, checkpoints,
    # retrieval/reassembly folders, list files, logs, PE_r1_*, PE_r2_*,
    # and numbered assembly FASTAs.

if __name__ == '__main__':
    assembly_list=['8_medium_S001_SPAdes_scaffolds.fasta','10_medium_cat_SPAdes_scaffolds.fasta']
    cleanup(assembly_list)
