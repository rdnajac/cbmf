import unittest
import io
import sys
import subprocess
from unittest.mock import patch, MagicMock
import logging
from pathlib import Path

from a_series_of_tubes.utils.logger import setup_logger
from a_series_of_tubes.utils.progressbar import ProgressBar
from a_series_of_tubes.genomemanager import GenomeManager
from a_series_of_tubes.run_script import run_script
from a_series_of_tubes.config import REFERENCE, FILES


class TestLogger(unittest.TestCase):
    def setUp(self):
        self.logger = setup_logger("test_logger", level=logging.DEBUG)
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
    @patch("a_series_of_tubes.run_script.subprocess.run")
    def test_run_script_success(self, mock_run):
        mock_run.return_value = MagicMock(returncode=0, stdout="Test output")
        result = run_script("test_script.sh")
        self.assertEqual(result, 0)

    def test_run_script_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            run_script("non_existent_script.sh")


if __name__ == "__main__":
    unittest.main()
