import os
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "BASALT"))

from runtime_utils import (  # noqa: E402
    atomic_gzip,
    bounded_parallel_resources,
    build_bowtie2_index,
    hybrid_skip_reason,
    load_id_set,
    open_text_auto,
    resolve_text_path,
    split_paired_sam_by_bin,
)
from Cleanup import cleanup, cleanup_checkm2_output  # noqa: E402
from bin_utils import (  # noqa: E402
    is_legacy_zero_bin,
    resolve_bin_filename,
    resolve_quality_key,
    strip_fasta_suffix,
)
from qc_utils import run_checkm2_predict, run_depth_summarizer  # noqa: E402


class RuntimeUtilsTests(unittest.TestCase):
    def test_bowtie2_index_uses_requested_threads(self):
        calls = []

        def fake_runner(command, check):
            calls.append((command, check))

        command = build_bowtie2_index("reference.fa", 12, runner=fake_runner)
        self.assertEqual(
            command,
            ["bowtie2-build", "--threads", "12", "reference.fa", "reference.fa"],
        )
        self.assertEqual(calls, [(command, True)])

    def test_atomic_gzip_is_resumable(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "reads.fq"
            source.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
            gz_path = atomic_gzip(source, compresslevel=1)
            self.assertFalse(source.exists())
            self.assertTrue(Path(gz_path).exists())
            self.assertEqual(resolve_text_path(source), gz_path)
            with open_text_auto(source, encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "@r1\nACGT\n+\n!!!!\n")
            self.assertEqual(atomic_gzip(source), gz_path)

    def test_failed_gzip_does_not_remove_source(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "reads.fq"
            source.write_text("data", encoding="utf-8")
            with self.assertRaises(ValueError):
                atomic_gzip(source, compresslevel="invalid")
            self.assertTrue(source.exists())
            self.assertFalse(Path(str(source) + ".gz.tmp").exists())

    def test_load_id_set_ignores_comments_and_extra_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "skip.txt"
            path.write_text("# reason\nbin1\nbin2\tfailed\n\n", encoding="utf-8")
            self.assertEqual(load_id_set(path), {"bin1", "bin2"})

    def test_hybrid_skip_reason_prefers_completed_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "bin1.fa"
            output.write_text(">contig\nACGT\n", encoding="utf-8")
            self.assertEqual(
                hybrid_skip_reason(output, "bin1", {"bin1"}),
                "completed output exists",
            )
            os.remove(output)
            self.assertEqual(
                hybrid_skip_reason(output, "bin1", {"bin1"}),
                "listed in hybrid skip file",
            )

    def test_cleanup_preserves_state_and_removes_scratch(self):
        with tempfile.TemporaryDirectory() as directory:
            previous = os.getcwd()
            try:
                os.chdir(directory)
                Path("mapping.sam").write_text("scratch", encoding="utf-8")
                Path("Basalt_checkpoint.txt").write_text("1st done", encoding="utf-8")
                Path("BestBinset_outlier_refined").mkdir()
                cleanup([])
                self.assertFalse(Path("mapping.sam").exists())
                self.assertTrue(Path("Basalt_checkpoint.txt").exists())
                self.assertTrue(Path("BestBinset_outlier_refined").exists())
            finally:
                os.chdir(previous)

    def test_checkm2_cleanup_requires_completed_report(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "diamond_output").mkdir()
            cleanup_checkm2_output(output)
            self.assertTrue((output / "diamond_output").exists())
            (output / "quality_report.tsv").write_text("Name\n", encoding="utf-8")
            cleanup_checkm2_output(output)
            self.assertFalse((output / "diamond_output").exists())

    def test_numeric_bin_suffix_is_not_treated_as_extension(self):
        self.assertEqual(strip_fasta_suffix("sample_genomes.001.fa"), "sample_genomes.001")
        self.assertFalse(is_legacy_zero_bin("sample_genomes.001"))
        self.assertTrue(is_legacy_zero_bin("sample_genomes.0.fa"))

    def test_padded_and_bare_bin_filenames_resolve(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "001.fa").write_text(">c\nACGT\n", encoding="utf-8")
            self.assertEqual(
                resolve_bin_filename(directory, "sample_genomes.1"),
                "001.fa",
            )

    def test_quality_keys_resolve_across_fasta_suffixes(self):
        quality = {"sample_genomes.001": {"Completeness": 90}}
        self.assertEqual(
            resolve_quality_key(quality, "sample_genomes.001.fa"),
            "sample_genomes.001",
        )

    def test_checkm2_empty_binset_gets_valid_empty_report(self):
        with tempfile.TemporaryDirectory() as directory:
            binset = Path(directory) / "bins"
            output = Path(directory) / "checkm2"
            binset.mkdir()
            report = run_checkm2_predict(binset, "fa", output, 4)
            self.assertTrue(Path(report).is_file())
            self.assertEqual(len(Path(report).read_text(encoding="utf-8").splitlines()), 1)

    def test_checkm2_replaces_partial_report(self):
        with tempfile.TemporaryDirectory() as directory:
            binset = Path(directory) / "bins"
            output = Path(directory) / "checkm2"
            binset.mkdir()
            output.mkdir()
            Path(binset, "bin1.fa").write_text(">c1\nACGT\n", encoding="utf-8")
            Path(binset, "bin2.fa").write_text(">c2\nTGCA\n", encoding="utf-8")
            Path(output, "quality_report.tsv").write_text(
                "Name\tCompleteness\n" "bin1\t90\n", encoding="utf-8"
            )
            calls = []

            def fake_runner(command, check):
                calls.append((command, check))
                output.mkdir(exist_ok=True)
                Path(output, "quality_report.tsv").write_text(
                    "Name\tCompleteness\n" "bin1\t90\n" "bin2\t80\n",
                    encoding="utf-8",
                )

            report = run_checkm2_predict(binset, "fa", output, 4, runner=fake_runner)
            self.assertEqual(len(calls), 1)
            self.assertIn("bin2", Path(report).read_text(encoding="utf-8"))

    def test_depth_summarizer_retries_until_output_is_usable(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "depth.txt"
            attempts = []

            class Result:
                returncode = 0

            def fake_runner(command, check):
                attempts.append((command, check))
                if len(attempts) == 1:
                    output.write_text("header\n", encoding="utf-8")
                else:
                    output.write_text("header\ncontig\t1\n", encoding="utf-8")
                return Result()

            self.assertEqual(
                run_depth_summarizer(output, ["reads.bam"], runner=fake_runner),
                str(output),
            )
            self.assertEqual(len(attempts), 2)

    def test_parallel_resources_are_bounded(self):
        self.assertEqual(bounded_parallel_resources(96, 256, 100), (4, 24, 64))
        self.assertEqual(bounded_parallel_resources(4, 32, 0), (0, 0, 0))

    def test_sam_is_split_in_one_pass_by_bin_and_mate(self):
        with tempfile.TemporaryDirectory() as directory:
            sam = Path(directory) / "reads.sam"
            sam.write_text(
                "@HD\tVN:1.6\n"
                "readA_1\t99\tbin1_contig\t1\t42\t4M\t=\t1\t0\tACGT\t!!!!\n"
                "readA_2\t147\tbin1_contig\t1\t42\t4M\t=\t1\t0\tTGCA\t####\n"
                "readB_1\t99\tbin2_contig\t1\t42\t4M\t=\t1\t0\tAAAA\t!!!!\n",
                encoding="utf-8",
            )
            counts, unmatched = split_paired_sam_by_bin(sam, 3, directory)
            self.assertEqual(counts, {"bin1": 1})
            self.assertEqual(unmatched, 1)
            r1 = Path(directory, "bin1_seq_R1.fq").read_text()
            r2 = Path(directory, "bin1_seq_R2.fq").read_text()
            self.assertIn("@3_readA 1", r1)
            self.assertIn("@3_readA 2", r2)
            self.assertIn("ACGT", r1)
            self.assertIn("TGCA", r2)

    def test_sam_pairing_preserves_legitimate_query_name_suffixes(self):
        with tempfile.TemporaryDirectory() as directory:
            sam = Path(directory) / "reads.sam"
            sam.write_text(
                "sample_read_A\t99\tbin1_contig\t1\t42\t4M\t=\t1\t0\tACGT\t!!!!\n"
                "sample_read_A\t147\tbin1_contig\t1\t42\t4M\t=\t1\t0\tTGCA\t####\n"
                "sample_read_B\t99\tbin1_contig\t1\t42\t4M\t=\t1\t0\tAAAA\t!!!!\n"
                "sample_read_B\t147\tbin1_contig\t1\t42\t4M\t=\t1\t0\tTTTT\t####\n",
                encoding="utf-8",
            )
            counts, unmatched = split_paired_sam_by_bin(sam, 2, directory)
            self.assertEqual(counts, {"bin1": 2})
            self.assertEqual(unmatched, 0)
            r1 = Path(directory, "bin1_seq_R1.fq").read_text()
            self.assertIn("@2_sample_read_A 1", r1)
            self.assertIn("@2_sample_read_B 1", r1)


if __name__ == "__main__":
    unittest.main()
