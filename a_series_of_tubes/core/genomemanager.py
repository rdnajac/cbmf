import os
import subprocess
import urllib.request
import threading
from pathlib import Path
from typing import Literal, Union, List
from ..utils import ColorPrinter as pr
from ..utils import ProgressBar
from ..config import GENOMES_MIRROR, REFERENCE, FILES

Species = Literal[tuple(REFERENCE.keys())]
Files = Literal[tuple(FILES.keys())]


class GenomeManager:
    def __init__(self):
        self.lock = threading.Lock()

    def _resolve_url(self, species: Species, file: Files) -> str:
        if species not in REFERENCE:
            raise ValueError(f"Invalid species: {species}")
        if file not in FILES:
            raise ValueError(f"Invalid file type: {file}")

        species_data = REFERENCE[species]
        file_suffix = FILES[file]
        return f"{GENOMES_MIRROR}/{species_data['uri']}/seqs_for_alignment_pipelines.ucsc_ids/{species_data['ref']}_full_analysis_set.{file_suffix}"

    def _fetch_file(self, species: Species, file: Files, show_progress: bool) -> None:
        try:
            file_url = self._resolve_url(species, file)
        except ValueError as e:
            pr.error(str(e))
            return

        file_path = Path(f"~/cbmf/genomes/{species}/{FILES[file]}").expanduser()

        pr.info(f"Downloading {species} {file} from {file_url}...")
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
                pr.success(f"Downloaded {species} {file} to {file_path}")

        except urllib.error.URLError as e:
            with self.lock:
                pr.error(f"Error downloading {file} for {species}: {str(e)}")
        except IOError as e:
            with self.lock:
                pr.error(f"Error writing {file} for {species}: {str(e)}")

    def download(self, species: Species, files: Union[List[Files], str]) -> None:
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

        if len(files) == 1:
            self._fetch_file(species, files[0], show_progress=True)
        else:
            threads = []
            for file in files:
                thread = threading.Thread(
                    target=self._fetch_file, args=(species, file, False)
                )
                threads.append(thread)
                thread.start()

            for thread in threads:
                thread.join()

        pr.success(f"All specified downloads completed for {species}")

    @staticmethod
    def decompress(path: str):
        files = os.listdir(path)

        for file in files:
            file_path = os.path.join(path, file)
            if file.endswith("index.tar.gz"):
                pr.info(f"Decompressing {file} with tar...")
                try:
                    subprocess.run(["tar", "-xvzf", file_path, "-C", path], check=True)
                    pr.success(f"Successfully extracted {file}")
                except subprocess.CalledProcessError as e:
                    pr.error(f"Failed to extract {file}: {e}")
            elif file.endswith(".gz"):
                pr.info(f"Decompressing {file} with gunzip...")
                try:
                    subprocess.run(["gunzip", "-v", file_path], check=True)
                    pr.success(f"Successfully decompressed {file}")
                except subprocess.CalledProcessError as e:
                    pr.error(f"Failed to decompress {file}: {e}")
            else:
                pr.info(f"Skipping {file}.")
