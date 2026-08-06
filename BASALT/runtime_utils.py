#!/usr/bin/env python

"""Shared runtime helpers for large and restartable BASALT runs."""

import gzip
import os
import shutil
import subprocess
from collections import OrderedDict


FALSE_VALUES = {"", "0", "false", "no", "off"}


class DiskSpaceError(RuntimeError):
    """Raised before launching more work when free disk space is too low."""


def env_flag(name, default=False):
    """Return a boolean environment option using common true/false spellings."""
    fallback = "1" if default else "0"
    return str(os.environ.get(name, fallback)).strip().lower() not in FALSE_VALUES


def build_bowtie2_index(reference, threads, threaded=True, runner=None):
    """Build a Bowtie2 index and fail immediately if the command fails."""
    minimum_gb = float(os.environ.get("BASALT_MIN_FREE_GB", "0"))
    ensure_min_free_gb(reference, minimum_gb, context="Bowtie2 index build")
    command = ["bowtie2-build"]
    if threaded:
        command.extend(["--threads", str(max(1, int(threads)))])
    command.extend([str(reference), str(reference)])
    if runner is None:
        runner = subprocess.run
    runner(command, check=True)
    return command


def gzip_level(default=1):
    """Read the optional gzip level while keeping it in gzip's valid range."""
    try:
        value = int(os.environ.get("BASALT_GZIP_LEVEL", str(default)))
    except ValueError:
        value = default
    return max(1, min(9, value))


def atomic_gzip(path, enabled=True, compresslevel=None, chunk_size=1024 * 1024):
    """Atomically gzip a completed file and remove the original on success."""
    path = os.fspath(path)
    if not enabled or path.endswith(".gz"):
        return path

    gz_path = path + ".gz"
    if not os.path.exists(path):
        return gz_path if os.path.exists(gz_path) else path
    if os.path.getsize(path) == 0:
        return path

    level = gzip_level() if compresslevel is None else compresslevel
    tmp_path = gz_path + ".tmp"
    try:
        with open(path, "rb") as source, gzip.open(
            tmp_path, "wb", compresslevel=max(1, min(9, int(level)))
        ) as destination:
            shutil.copyfileobj(source, destination, length=chunk_size)
        os.replace(tmp_path, gz_path)
        os.remove(path)
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise
    return gz_path


def resolve_text_path(path):
    """Resolve a text intermediate stored either plain or as ``.gz``."""
    path = os.fspath(path)
    if os.path.exists(path):
        return path
    gz_path = path + ".gz"
    return gz_path if os.path.exists(gz_path) else path


def open_text_auto(path, mode="rt", **kwargs):
    """Open plain or gzipped text using the same call site."""
    resolved = resolve_text_path(path)
    if resolved.endswith(".gz"):
        return gzip.open(resolved, mode, **kwargs)
    return open(resolved, mode, **kwargs)


def free_space_gb(path):
    """Return free filesystem capacity in decimal gigabytes."""
    probe = os.path.abspath(os.fspath(path))
    while not os.path.exists(probe):
        parent = os.path.dirname(probe)
        if parent == probe:
            break
        probe = parent
    return shutil.disk_usage(probe).free / 1_000_000_000


def ensure_min_free_gb(path, minimum_gb, context="operation"):
    """Abort before an operation when the target filesystem is nearly full."""
    minimum_gb = float(minimum_gb)
    if minimum_gb <= 0:
        return free_space_gb(path)
    available = free_space_gb(path)
    if available <= minimum_gb:
        raise DiskSpaceError(
            "Refusing to start {}: {:.2f} GB free at {}, threshold {:.2f} GB".format(
                context, available, os.path.abspath(os.fspath(path)), minimum_gb
            )
        )
    return available


def load_id_set(path):
    """Load the first whitespace-delimited token from a comment-aware ID file."""
    if not path or not os.path.isfile(path):
        return set()
    identifiers = set()
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                identifiers.add(stripped.split()[0])
    return identifiers


def hybrid_skip_reason(output_path, item_id=None, skip_ids=None):
    """Return why a hybrid assembly may be skipped, or ``None`` to run it."""
    if os.path.isfile(output_path) and os.path.getsize(output_path) > 0:
        return "completed output exists"
    if item_id is not None and skip_ids and str(item_id) in skip_ids:
        return "listed in hybrid skip file"
    return None


def bounded_parallel_resources(total_threads, total_ram, item_count,
                               target_threads=24, target_ram=48):
    """Choose a conservative worker count and divide resources between jobs."""
    item_count = max(0, int(item_count))
    if item_count == 0:
        return 0, 0, 0
    total_threads = max(1, int(total_threads))
    total_ram = max(1, int(total_ram))
    jobs_by_threads = max(1, total_threads // max(1, int(target_threads)))
    jobs_by_ram = max(1, total_ram // max(1, int(target_ram)))
    jobs = min(item_count, jobs_by_threads, jobs_by_ram)
    return jobs, max(1, total_threads // jobs), max(1, total_ram // jobs)


def split_paired_sam_by_bin(sam_file, dataset_index, output_dir='.',
                            max_open_files=64):
    """Split query-oriented paired SAM records into per-bin FASTQ files.

    Secondary, supplementary, and unmapped alignments are ignored. Pending
    records are removed as soon as both mates assigned to the same bin arrive.
    """
    output_dir = os.path.abspath(os.fspath(output_dir))
    os.makedirs(output_dir, exist_ok=True)
    pending = {}
    pair_counts = {}
    writers = OrderedDict()

    def writer_for(bin_id, mate):
        key = (bin_id, mate)
        handle = writers.pop(key, None)
        if handle is None:
            if len(writers) >= max(1, int(max_open_files)):
                _, stale = writers.popitem(last=False)
                stale.close()
            path = os.path.join(output_dir, '{}_seq_R{}.fq'.format(bin_id, mate))
            handle = open(path, 'a', encoding='utf-8')
        writers[key] = handle
        return handle

    try:
        with open(sam_file, 'r', encoding='utf-8', errors='replace') as sam:
            for line_number, line in enumerate(sam, 1):
                if not line or line.startswith('@'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 11 or fields[2] == '*':
                    continue
                try:
                    flag = int(fields[1])
                except ValueError:
                    continue
                if flag & (0x4 | 0x100 | 0x800):
                    continue
                if flag & 0x40:
                    mate = '1'
                elif flag & 0x80:
                    mate = '2'
                else:
                    mate = fields[0].rsplit('_', 1)[-1]
                    if mate not in {'1', '2'}:
                        continue
                bin_id = fields[2].split('_', 1)[0]
                read_base = fields[0]
                if read_base.endswith(('_1', '_2', '/1', '/2')):
                    read_base = read_base[:-2]
                read_key = (bin_id, read_base)
                bucket = pending.setdefault(read_key, {})
                bucket.setdefault(mate, (read_base, fields[9], fields[10]))
                if '1' in bucket and '2' in bucket:
                    for read_mate in ('1', '2'):
                        read_id, sequence, quality = bucket[read_mate]
                        header = '@{}_{} {}'.format(dataset_index, read_id, read_mate)
                        writer_for(bin_id, read_mate).write(
                            '{}\n{}\n+\n{}\n'.format(header, sequence, quality)
                        )
                    pair_counts[bin_id] = pair_counts.get(bin_id, 0) + 1
                    del pending[read_key]
                if line_number % 1000000 == 0:
                    print('Parsed', line_number, 'SAM lines')
    finally:
        for handle in writers.values():
            handle.close()
    return pair_counts, len(pending)
