import os
import subprocess
import urllib.request
import threading
from pathlib import Path
from typing import Union, List
from ..utils.logger import logger
from ..utils import ProgressBar
from ..config import GENOMES_MIRROR, REFERENCE, FILES


class GenomeManager:
    def __init__(self):
        self.lock = threading.Lock()

    def _resolve_url(self, species: str, file: str) -> str:
        valid_species = ["mouse", "human"]
        if species not in valid_species:
            raise ValueError(f"Invalid species: {species}")
        if file not in FILES:
            raise ValueError(f"Invalid file type: {file}")

        species_data = REFERENCE[species]
        file_suffix = FILES[file]
        return f"{GENOMES_MIRROR}/{species_data['uri']}/seqs_for_alignment_pipelines.ucsc_ids/{species_data['ref']}_full_analysis_set.{file_suffix}"

    def _fetch_file(self, species: str, file: str, show_progress: bool) -> None:
        """Download a file for a given species."""

        file_url = self._resolve_url(species, file)
        file_path = Path(f"~/genomes/{species}/{FILES[file]}").expanduser()

        try:
            with urllib.request.urlopen(file_url) as response:
                total_size = int(response.info().get("Content-Length", -1))
                downloaded_size = 0

                with open(file_path, "wb") as out_file:
                    progress_bar = (
                        ProgressBar(
                            total_size, prefix=f"{file:<15} ({total_size})", length=30
                        )
                        if show_progress and total_size > 0
                        else None
                    )

                    while True:
                        buffer = response.read(8192)
                        if not buffer:
                            break
                        downloaded_size += len(buffer)
                        out_file.write(buffer)

                        if progress_bar:
                            progress_bar.update(downloaded_size)

            with self.lock:
                print(f"Downloaded {file} for {species}")

        except urllib.error.URLError as e:
            with self.lock:
                logger.error(f"Failed to download {file} for {species}: {str(e)}")
        except IOError as e:
            with self.lock:
                logger.error(f"Failed to write {file} for {species}: {str(e)}")

    def download(self, species: str, files: Union[List[str], str]) -> None:
        """Download files for a given species."""

        # Validate species
        if species not in REFERENCE:
            raise ValueError(f"Invalid species: {species}")

        # Determine files to download
        if isinstance(files, str):
            if files == "ALL":
                files = list(FILES.keys())
            else:
                raise ValueError(
                    "Invalid string input for files. Use 'ALL' to download all files or provide a list of files."
                )

        if not isinstance(files, list):
            raise ValueError(
                "Invalid input for files parameter. Must be a list or 'ALL'."
            )

        # Validate each file type
        invalid_files = [file for file in files if file not in FILES]
        if invalid_files:
            raise ValueError(f"Invalid file types: {', '.join(invalid_files)}")

        # Download files
        if len(files) == 1:
            self._fetch_file(species, files[0], show_progress=True)
        else:
            threads = []
            for file in files:
                # if the file exists, skip it
                if os.path.exists(f"~/genomes/{species}/{FILES[file]}"):
                    logger.warning(f"Skipping {file} as it already exists.")
                    continue

                thread = threading.Thread(
                    target=self._fetch_file, args=(species, file, False)
                )
                threads.append(thread)
                thread.start()

            # Wait for all threads to complete
            for thread in threads:
                thread.join()

        # Decompress files
        self.decompress(f"~/genomes/{species}")

    @staticmethod
    def decompress(path: str):
        """Decompress files in the given path."""

        files = os.listdir(path)

        for file in files:
            file_path = os.path.join(path, file)
            if file.endswith("index.tar.gz"):
                try:
                    subprocess.run(["tar", "-xvzf", file_path, "-C", path], check=True)
                except subprocess.CalledProcessError as e:
                    logger.error(f"Failed to decompress {file}: {str(e)}")
            elif file.endswith(".gz"):
                try:
                    subprocess.run(["gunzip", "-v", file_path], check=True)
                except subprocess.CalledProcessError as e:
                    logger.error(f"Failed to decompress {file}: {str(e)}")
            else:
                logger.warning(f"Skipping {file} as it is not a compressed file.")
