#!/bin/bash
bcl2fastq --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls --auto-set-to-zero-barcode-mismatches --find-adapters-with-sliding-window --adapter-stringency 0.9 --mask-short-adapter-reads 35 --minimum-trimmed-read-length 35 -R run5 --sample-sheet ./run5/SampleSheet.csv -o ./run5_fastq
