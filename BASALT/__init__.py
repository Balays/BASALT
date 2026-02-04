"""
BASALT: A modular metagenomic binning and refinement pipeline.

This package exposes the high-level CLI entry point in ``BASALT.py``
and a number of internal modules implementing individual steps of
the pipeline (autobinning, dereplication, outlier removal, contig
retrieval, reassembly, and model utilities).

Typical users should invoke BASALT via the command-line script
``BASALT.py`` rather than importing these internals directly. The
package namespace is kept minimal on purpose.
"""

__all__ = []
