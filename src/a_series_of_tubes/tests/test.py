import unittest
import io
import sys
from unittest.mock import patch, MagicMock
from pathlib import Path

from ..utils.logger import setup_logger
from ..utils.progressbar import ProgressBar
from ..core.genomemanager import GenomeManager
from ..utils.run_script import run_script
from ..utils.cli import parse_args
from ..config import REFERENCE, FILES


class TestLogger(unittest.TestCase):
    def setUp(self):
        self.logger = setup_logger("test_logger", level="DEBUG")
        self.stream = io.StringIO()
        self.handler = logging.StreamHandler(self.stream)
        self.logger.addHandler(self.handler)

    def test_logger_levels(self):
        self.logger.debug("Debug message")
        self.logger.info("Info message")
        self.logger.warning("Warning message")
        self.logger.error("Error message")
        self.logger.critical("Critical message")

        log_output = self.stream.getvalue()
        self.assertIn("Debug message", log_output)
        self.assertIn("Info message", log_output)
        self.assertIn("Warning message", log_output)
        self.assertIn("Error message", log_output)
        self.assertIn("Critical message", log_output)


class TestProgressBar(unittest.TestCase):
    def test_progress_bar_update(self):
        with patch("sys.stdout", new=io.StringIO()) as fake_out:
            bar = ProgressBar(
                total=100, prefix="Progress:", suffix="Complete", length=20
            )
            bar.update(50)
            output = fake_out.getvalue()
            self.assertIn("Progress:", output)
            self.assertIn("50.0%", output)
            self.assertIn("Complete", output)

    def test_progress_bar_finish(self):
        with patch("sys.stdout", new=io.StringIO()) as fake_out:
            bar = ProgressBar(
                total=100, prefix="Progress:", suffix="Complete", length=20
            )
            bar.finish()
            output = fake_out.getvalue()
            self.assertIn("Progress:", output)
            self.assertIn("100.0%", output)
            self.assertIn("Complete", output)


class TestGenomeManager(unittest.TestCase):
    def setUp(self):
        self.genome_manager = GenomeManager()

    def test_resolve_url(self):
        for species in REFERENCE.keys():
            for file in FILES.keys():
                url = self.genome_manager._resolve_url(species, file)
                self.assertIsInstance(url, str)
                self.assertTrue(url.startswith("ftp://"))

    def test_resolve_url_invalid_species(self):
        with self.assertRaises(ValueError):
            self.genome_manager._resolve_url("invalid_species", "fasta")

    def test_resolve_url_invalid_file(self):
        with self.assertRaises(ValueError):
            self.genome_manager._resolve_url("human", "invalid_file")


class TestRunScript(unittest.TestCase):
    @patch("subprocess.run")
    def test_run_script_success(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stdout="Test output")
        result = run_script("test_script.sh")
        self.assertEqual(result, 0)

    @patch("subprocess.run")
    def test_run_script_failure(self, mock_run):
        mock_run.side_effect = subprocess.CalledProcessError(
            1, "test_script.sh", stderr="Test error"
        )
        with self.assertRaises(subprocess.CalledProcessError):
            run_script("test_script.sh")

    def test_run_script_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            run_script("non_existent_script.sh")


class TestCLI(unittest.TestCase):
    def test_parse_args_download(self):
        args = parse_args(["download", "--species", "human", "fasta"])
        self.assertEqual(args.command, "download")
        self.assertEqual(args.species, "human")
        self.assertEqual(args.files, ["fasta"])

    def test_parse_args_align(self):
        args = parse_args(
            [
                "align",
                "bwa",
                "-i",
                "/input/dir",
                "-o",
                "/output/dir",
                "-r",
                "/reference/file",
            ]
        )
        self.assertEqual(args.command, "align")
        self.assertEqual(args.aligner, "bwa")
        self.assertEqual(str(args.input_directory), "/input/dir")
        self.assertEqual(str(args.output_directory), "/output/dir")
        self.assertEqual(str(args.reference), "/reference/file")

    def test_parse_args_verbosity(self):
        args = parse_args(["-v", "download", "--species", "human", "fasta"])
        self.assertEqual(args.log_level, logging.INFO)

        args = parse_args(["-vv", "download", "--species", "human", "fasta"])
        self.assertEqual(args.log_level, logging.DEBUG)


if __name__ == "__main__":
    unittest.main()
