#!/usr/bin/env python

"""Checked wrappers around external quality and depth tools."""

import os
import shlex
import shutil
import subprocess

from Cleanup import cleanup_checkm2_output
from bin_utils import strip_fasta_suffix


CHECKM2_QUALITY_REPORT_HEADER = (
    'Name\tCompleteness\tContamination\tCompleteness_Model_Used\t'
    'Translation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\t'
    'Genome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\t'
    'Max_Contig_Length\tAdditional_Notes\n'
)


def _checkm2_report_ids(report):
    if not os.path.isfile(report) or os.path.getsize(report) == 0:
        return set()
    identifiers = set()
    with open(report, 'r', encoding='utf-8', errors='replace') as handle:
        for line_number, line in enumerate(handle):
            if line_number == 0:
                continue
            name = line.rstrip('\n').split('\t', 1)[0].strip()
            if name:
                identifiers.add(strip_fasta_suffix(name))
    return identifiers


def run_checkm2_predict(binset, extension, output_folder, threads, runner=None):
    """Run CheckM2 once, reuse complete output, and reject partial output."""
    binset = os.fspath(binset)
    output_folder = os.fspath(output_folder)
    report = os.path.join(output_folder, 'quality_report.tsv')
    suffix = '.'+str(extension).lstrip('.').lower()
    bin_files = []
    if os.path.isdir(binset):
        bin_files = [
            name for name in os.listdir(binset)
            if name.lower().endswith(suffix)
            and os.path.isfile(os.path.join(binset, name))
            and os.path.getsize(os.path.join(binset, name)) > 0
        ]
    expected_ids = {strip_fasta_suffix(name) for name in bin_files}
    if expected_ids and expected_ids.issubset(_checkm2_report_ids(report)):
        print('Reusing complete CheckM2 quality report: '+report)
        cleanup_checkm2_output(output_folder)
        return report
    if not bin_files:
        if os.path.isdir(output_folder):
            shutil.rmtree(output_folder)
        os.makedirs(output_folder)
        with open(report, 'w', encoding='utf-8') as handle:
            handle.write(CHECKM2_QUALITY_REPORT_HEADER)
        print('No '+suffix+' bins found in '+binset+'; wrote an empty CheckM2 report.')
        return report

    if os.path.isdir(output_folder):
        print('Removing incomplete CheckM2 output: '+output_folder)
        shutil.rmtree(output_folder)

    command = [
        'checkm2', 'predict', '-t', str(max(1, int(threads))),
        '-i', binset, '-x', str(extension).lstrip('.'), '-o', output_folder,
    ]
    if runner is None:
        runner = subprocess.run
    runner(command, check=True)
    if not os.path.isfile(report) or os.path.getsize(report) == 0:
        raise RuntimeError('CheckM2 did not create a usable report: '+report)
    missing_ids = sorted(expected_ids.difference(_checkm2_report_ids(report)))
    if missing_ids:
        raise RuntimeError(
            'CheckM2 report is incomplete; {} of {} bin IDs are missing, including {}'.format(
                len(missing_ids), len(expected_ids), ', '.join(missing_ids[:5])
            )
        )
    cleanup_checkm2_output(output_folder)
    return report


def run_depth_summarizer(output_depth, bam_files, attempts=2, runner=None):
    """Run MetaBAT's depth summarizer and require a header plus data row."""
    output_depth = os.fspath(output_depth)
    if isinstance(bam_files, str):
        bam_files = shlex.split(bam_files)
    command = [
        'jgi_summarize_bam_contig_depths', '--outputDepth', output_depth,
    ] + [os.fspath(path) for path in bam_files]
    if runner is None:
        runner = subprocess.run

    last_returncode = None
    for attempt in range(1, max(1, int(attempts))+1):
        if os.path.exists(output_depth):
            os.remove(output_depth)
        result = runner(command, check=False)
        last_returncode = getattr(result, 'returncode', 0)
        if last_returncode == 0 and os.path.isfile(output_depth):
            with open(output_depth, 'r', encoding='utf-8', errors='replace') as handle:
                if sum(1 for _ in handle) >= 2:
                    return output_depth
        if attempt < attempts:
            print('Depth generation failed; retrying '+output_depth)
    raise RuntimeError(
        'jgi_summarize_bam_contig_depths failed for {} (exit status {})'.format(
            output_depth, last_returncode
        )
    )
