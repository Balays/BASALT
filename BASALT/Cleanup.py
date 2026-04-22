#!/usr/bin/env python

"""
Helper utilities for cleaning up intermediate BASALT files.
"""

import glob
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


def _remove_patterns(base_dir, patterns):
    for pattern in patterns:
        for path in glob.glob(os.path.join(base_dir, pattern)):
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
    Remove intermediate index/coverage files and archive key matrices.

    Parameters
    ----------
    assembly_list : list
        List of assemblies. Currently unused but kept for API compatibility.
    """
    if not cleanup_enabled():
        return
    os.system('rm *.njs *.ndb *.nto *.ntf *.not *.nos')
    os.mkdir('Coverage_depth_connection_SimilarBin_files_backup')
    os.system(
        'mv *.depth.txt Coverage_matrix_* Combat_* condense_connections_* '
        'Connections_* Similar_bins.txt Coverage_depth_connection_SimilarBin_files_backup'
    )
    os.system(
        'tar -zcvf Coverage_depth_connection_SimilarBin_files_backup.tar.gz '
        'Coverage_depth_connection_SimilarBin_files_backup'
    )
    os.system('rm -rf Coverage_depth_connection_SimilarBin_files_backup')
    os.system('rm -rf *_kmer bin_coverage Bin_coverage_after_contamination_removal bin_comparison_folder bin_extract-eleminated-selected_contig Bins_blast_output')
    os.system('tar -zcvf Group_comparison_files.tar.gz *_comparison_files')
    os.system('tar -zcvf Group_Bestbinset.tar.gz *_BestBinsSet')
    # os.system('tar -zcvf Group_checkm.tar.gz *_checkm')
    os.system('tar -zcvf Group_genomes.tar.gz *_genomes')
    os.system('rm -rf *_sr_bins_seq')
    os.system('tar -zcvf Binsets_backup.tar.gz BestBinse*') ###
    os.system('rm -rf *_comparison_files *_checkm *_genomes *_BestBinset *_BestBinsSet BestBinse* Deep_retrieved_bins coverage_deep_refined_bins S6_coverage_filtration_matrix S6_TNF_filtration_matrix split_blast_output TNFs_deep_refined_bins')
    os.system('rm *_checkpoint.txt')
    os.system('rm -rf Merged_seqs_*')
    os.system('rm -rf *.bt2 Outlier_in_threshold* Summary_threshold* Refined_total_bins_contigs.fa Total_bins.fa') 
    for i in range(1,20):
        os.system('rm -rf *_deep_retrieval_'+str(i))
    os.system('rm -rf *_MP_1 *_MP_2 *_gf_lr_polished *_gf_lr *_gf_lr_mod *_gf_lr_checkm *_long_read')
    os.system('rm Bin_reads_summary.txt Depth_total.txt Basalt_log.txt Assembly_mo_list.txt Assembly_MoDict.txt *_gf_lr_blasted.txt Bestbinset_list.txt Bin_extract_contigs_after_coverage_filtration.txt Bin_lw.txt Bin_record_error.txt')
    os.system('rm Bins_folder.txt BLAST_output_error.txt Concoct_* condensed.cytoscape* cytoscape.*')
    os.system('rm Hybrid_re-assembly_status.txt Mapping_log_* OLC_merged_error_blast_results.txt Potential_contaminted_seq_vari.txt Reassembled_bins_comparison.txt Rejudge_clean.txt')
    os.system('rm Remained_seq* Remapped_depth_test.txt Re-mapped_depth.txt Remapping.fasta TNFs_exceptional_contigs.txt Total_contigs_after_OLC_reassembly.fa')
    os.system('rm PE_r1_* PE_r2_*')
    for i in range(1, len(assembly_list)+1):
        os.system('rm '+str(i)+'_'+str(assembly_list[i-1]))

if __name__ == '__main__':
    assembly_list=['8_medium_S001_SPAdes_scaffolds.fasta','10_medium_cat_SPAdes_scaffolds.fasta']
    cleanup(assembly_list)
