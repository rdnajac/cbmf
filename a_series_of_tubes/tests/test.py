import unittest
import io
import sys
from unittest.mock import patch, MagicMock
from pathlib import Path

# Import the modules to test
from cbmf.utils import ColorPrinter, ProgressBar
from cbmf.core import GenomeManager
from cbmf.config import REFERENCE, FILES


class TestColorPrinter(unittest.TestCase):
    def setUp(self):
        self.held_output = io.StringIO()
        sys.stderr = self.held_output

    def tearDown(self):
        sys.stderr = sys.__stderr__

    def test_print_status(self):
        ColorPrinter.print_status("success", "Test message")
        self.assertIn("[Success]: Test message", self.held_output.getvalue())

    def test_success(self):
        ColorPrinter.success("Success message")
        self.assertIn("[Success]: Success message", self.held_output.getvalue())

    def test_info(self):
        ColorPrinter.info("Info message")
        self.assertIn("[Info]: Info message", self.held_output.getvalue())

    def test_warning(self):
        ColorPrinter.warning("Warning message")
        self.assertIn("[Warning]: Warning message", self.held_output.getvalue())

    def test_error(self):
        ColorPrinter.error("Error message")
        self.assertIn("[Error]: Error message", self.held_output.getvalue())


class TestProgressBar(unittest.TestCase):
    def test_update(self):
        with patch("sys.stdout", new=io.StringIO()) as fake_out:
            bar = ProgressBar(
                total=100, prefix="Progress:", suffix="Complete", length=20
            )
            bar.update(50)
            output = fake_out.getvalue()
            self.assertIn("Progress:", output)
            self.assertIn("50.0%", output)
            self.assertIn("Complete", output)

    def test_finish(self):
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

    @patch("cbmf.core.genomemanager.urllib.request.urlopen")
    @patch("cbmf.core.genomemanager.open")
    def test_fetch_file(self, mock_open, mock_urlopen):
        mock_response = MagicMock()
        mock_response.info.return_value.get.return_value = "1000"
        mock_response.read.side_effect = [b"data"] * 10 + [b""]
        mock_urlopen.return_value.__enter__.return_value = mock_response

        self.genome_manager._fetch_file("human", "fasta", show_progress=False)

        mock_urlopen.assert_called_once()
        mock_open.assert_called_once()


if __name__ == "__main__":
    unittest.main()
