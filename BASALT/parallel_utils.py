#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Helpers for parallelizing independent BASALT file-preparation tasks.
"""

import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor


def _detect_archive_type(path):
    parts = str(path).split(".")
    if len(parts) < 2:
        raise ValueError("Input format error! Please check the input file.")
    if parts[-1] in ("fq", "fastq", "fa", "fna", "fasta"):
        return 1
    if parts[-1] == "zip":
        return 2
    if parts[-1] == "gz":
        if len(parts) >= 3 and parts[-2] == "tar":
            return 3
        return 4
    raise ValueError("Input format error! Please check the input file.")


def _run_shell(command):
    completed = subprocess.run(command, shell=True)
    if completed.returncode != 0:
        raise RuntimeError("Command failed: " + command)


def _expand_one(input_path, archive_type, exists_probe, expanded_path):
    if os.path.exists(exists_probe) or os.path.exists(expanded_path):
        return
    quoted_input = shlex.quote(str(input_path))
    quoted_output = shlex.quote(str(expanded_path))
    if archive_type == 2:
        command = "unzip " + quoted_input
    elif archive_type == 3:
        command = "tar -zxf " + quoted_input
    elif archive_type == 4:
        command = "gunzip -c " + quoted_input + " > " + quoted_output
    else:
        return
    _run_shell(command)


def _run_parallel(tasks, max_workers):
    if not tasks:
        return
    workers = max(1, min(int(max_workers), len(tasks)))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_expand_one, *task) for task in tasks]
        for future in futures:
            future.result()


def _normalize_result_path(expanded_path, default_extension):
    if default_extension == ".fq":
        if ".fq" not in expanded_path and ".fastq" not in expanded_path:
            return expanded_path + ".fq"
    elif default_extension == ".fa":
        if ".fa" not in expanded_path and ".fna" not in expanded_path and ".fasta" not in expanded_path:
            return expanded_path + ".fa"
    return expanded_path


def prepare_paired_datasets(datasets, pwd, max_workers):
    if not datasets:
        return datasets
    first_key = next(iter(datasets.keys()))
    archive_type = _detect_archive_type(datasets[first_key][0])
    if archive_type == 1:
        return datasets

    prepared = {}
    tasks = []
    for item in datasets.keys():
        prepared[item] = []
        for read_index, prefix in ((0, "PE_r1_"), (1, "PE_r2_")):
            source_path = str(datasets[item][read_index])
            if archive_type == 2:
                expanded_path = source_path.split(".zip")[0]
            elif archive_type == 3:
                expanded_path = source_path.split(".tar.gz")[0]
            else:
                expanded_path = source_path.split(".gz")[0]
            prepared[item].append(_normalize_result_path(expanded_path, ".fq"))
            tasks.append(
                (
                    source_path,
                    archive_type,
                    os.path.join(pwd, prefix + expanded_path),
                    expanded_path,
                )
            )

    _run_parallel(tasks, max_workers)
    return prepared


def prepare_sequence_files(sequence_paths, pwd, max_workers, default_extension):
    if not sequence_paths:
        return sequence_paths
    archive_type = _detect_archive_type(sequence_paths[0])
    if archive_type == 1:
        return sequence_paths

    prepared = []
    tasks = []
    for item in sequence_paths:
        source_path = str(item)
        if archive_type == 2:
            expanded_path = source_path.split(".zip")[0]
        elif archive_type == 3:
            expanded_path = source_path.split(".tar.gz")[0]
        else:
            expanded_path = source_path.split(".gz")[0]
        prepared.append(_normalize_result_path(expanded_path, default_extension))
        tasks.append(
            (
                source_path,
                archive_type,
                os.path.join(pwd, expanded_path),
                expanded_path,
            )
        )

    _run_parallel(tasks, max_workers)
    return prepared
