#!/usr/bin/env python

"""Bin identifier and FASTA resolution helpers shared by BASALT stages."""

import os
import re


FASTA_SUFFIXES = ('.fasta', '.fna', '.fas', '.fsa', '.fa')


def is_fasta_file(filename):
    return str(filename).lower().endswith(FASTA_SUFFIXES)


def strip_fasta_suffix(filename):
    """Strip only a real FASTA suffix, preserving numeric bin suffixes."""
    value = os.path.basename(os.fspath(filename))
    lower = value.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return value[:-len(suffix)]
    return value


def _numeric_suffix(identifier):
    base = strip_fasta_suffix(identifier)
    token = base.rsplit('_genomes.', 1)[-1] if '_genomes.' in base else base
    if re.fullmatch(r'\d+', token):
        return int(token)
    return None


def is_legacy_zero_bin(identifier):
    """Identify the exact historical ``_genomes.0`` placeholder bin."""
    base = strip_fasta_suffix(identifier)
    return base.endswith('_genomes.0') and '_semibin_genomes.0' not in base


def _match_score(candidate, requested):
    candidate_base = strip_fasta_suffix(candidate)
    requested_base = strip_fasta_suffix(requested)
    if candidate_base == requested_base:
        return 3
    candidate_number = _numeric_suffix(candidate_base)
    requested_number = _numeric_suffix(requested_base)
    if candidate_number is not None and candidate_number == requested_number:
        candidate_qualified = '_genomes.' in candidate_base
        requested_qualified = '_genomes.' in requested_base
        if candidate_qualified == requested_qualified:
            return 2
        return 1
    return 0


def resolve_bin_filename(folder, bin_id):
    """Resolve qualified, padded, or bare numeric IDs to one FASTA filename."""
    folder = os.fspath(folder)
    if not os.path.isdir(folder):
        return None
    requested_name = os.path.basename(os.fspath(bin_id))
    exact = os.path.join(folder, requested_name)
    if os.path.isfile(exact) and is_fasta_file(exact):
        return requested_name

    matches = []
    for filename in sorted(os.listdir(folder)):
        path = os.path.join(folder, filename)
        if not os.path.isfile(path) or not is_fasta_file(filename):
            continue
        score = _match_score(filename, requested_name)
        if score:
            matches.append((score, filename))
    if not matches:
        return None
    highest = max(score for score, _ in matches)
    best = [filename for score, filename in matches if score == highest]
    return best[0] if len(best) == 1 else None


def resolve_bin_path(folder, bin_id):
    filename = resolve_bin_filename(folder, bin_id)
    return os.path.join(os.fspath(folder), filename) if filename else None


def resolve_quality_key(quality, bin_id):
    """Resolve a FASTA/bin identifier against CheckM or CheckM2 table keys."""
    requested = os.fspath(bin_id)
    if requested in quality:
        return requested
    matches = []
    for key in quality:
        score = _match_score(key, requested)
        if score:
            matches.append((score, key))
    if not matches:
        return None
    highest = max(score for score, _ in matches)
    best = [key for score, key in matches if score == highest]
    return best[0] if len(best) == 1 else None
