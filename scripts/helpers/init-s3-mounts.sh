#!/bin/bash
# lab-aaf-ngs-data-archive
# lab-aaf-scratch
# lab-aaf-server-backup
# lab-tp-rstudio-scratch
# labaf-missionbio-tapestri

sudo mkdir -p /mnt/S3/{lab-aaf-ngs-data-archive,lab-aaf-scratch,lab-aaf-server-backup,lab-tp-rstudio-scratch,labaf-missionbio-tapestri}

sudo mount-s3 --allow-other lab-aaf-ngs-data-archive /mnt/S3/lab-aaf-ngs-data-archive
sudo mount-s3 --allow-other lab-aaf-scratch /mnt/S3/lab-aaf-scratch
sudo mount-s3 --allow-other lab-aaf-server-backup /mnt/S3/lab-aaf-server-backup
sudo mount-s3 --allow-other lab-tp-rstudio-scratch /mnt/S3/lab-tp-rstudio-scratch
sudo mount-s3 --allow-other labaf-missionbio-tapestri /mnt/S3/labaf-missionbio-tapestri
