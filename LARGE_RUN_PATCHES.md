# BASALT large-run patches

This branch starts from EMBL-PKU/BASALT commit `5f51ba3` and contains a
reviewed subset of changes for large, restartable CheckM2-based runs. It does
not contain dataset-specific paths, identifiers, manifests, recovery scripts,
or checkpoint edits.

## Scope

- Build Bowtie2 indexes with the requested thread count and checked exit status.
- Stream paired SAM parsing once instead of retaining and rereading the full
  mapping three times.
- Bound concurrent SPAdes/IDBA jobs and divide threads and RAM between workers.
- Skip completed hybrid reassemblies and optionally skip IDs listed in a
  user-supplied text file.
- Validate CheckM2 reports against all non-empty input bins before reuse,
  cleanup, or downstream filtering.
- Retry depth summarization and reject header-only depth output.
- Preserve numeric and zero-padded bin IDs and resolve common FASTA suffixes.
- Keep cleanup conservative: disposable mapping/index scratch is removed, but
  checkpoints, quality reports, binsets, matrices, logs, and reassembly state
  remain available for resume and downstream analysis.
- Optionally gzip completed per-bin long-read intermediates atomically.
- Refuse new large operations when free disk space is below a configurable
  threshold.

## New command-line options

```text
--min-free-gb N
--cleanup / --no-cleanup
--gzip-long-read-intermediates
--gzip-level 1..9
--skip-semibin
--hybrid-skip-file FILE
--initial-drep-start-iteration N
```

`--hybrid-skip-file` accepts one bin ID per line; text after the first
whitespace-delimited field is ignored. `--initial-drep-start-iteration`
requires the corresponding `Iteration_N_genomes`, `Contigs_iteration_N.fa`,
and `Coverage_matrix_for_binning_iteration_N.txt` artifacts.

## Deliberately excluded

- Dataset-specific paths, sample names, skip lists, and resource settings.
- One-off rescue, salvage, and live-workspace patch scripts.
- Debug snapshots and checkpoint rewrites.
- Broad refactors unrelated to correctness, runtime, storage, or restart safety.
- Behavioral changes to the legacy CheckM quality-control backend.

## Validation

The focused unit suite covers atomic compression, disk-aware threaded indexing,
conservative cleanup, CheckM2 completeness checks, depth retry, padded bin-ID
resolution, bounded worker resources, hybrid skip behavior, and one-pass paired
SAM splitting.

```bash
python -m unittest discover -s tests -v
python -m compileall -q BASALT
```

Full end-to-end validation still requires the external BASALT toolchain and
representative sequencing input data.
