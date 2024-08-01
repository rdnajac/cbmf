#!/bin/bash
#
## Setup for rnaseq pipeline

scp my-ec2:~/rnaseq/agx/conts.tsv ./agx_conts.tsv
scp my-ec2:~/rnaseq/agx/counts.tsv.summary ./agx_featureCounts.summary
scp my-ec2:~/rnaseq/ra/ra_counts.tsv ./ra_counts.tsv
scp my-ec2:~/rnaseq/ra/ra_counts.tsv.summary ./ra_featureCounts.summary
scp my-ec2:~/rnaseq/clone/counts.tsv ./clone_counts.tsv
scp my-ec2:~/rnaseq/clone/counts.tsv.summary ./clone_featureCounts.summary

scp my-ec2:~/rnaseq/ra/fin_dmso.tsv ./fin_dmso.tsv
scp my-ec2:~/rnaseq/ra/fin_dmso.tsv.summary ./fin_dmso.summary
scp my-ec2:~/rnaseq/ra/oza_dmso.tsv ./oza_dmso.tsv
scp my-ec2:~/rnaseq/ra/oza_dmso.tsv.summary ./oza_dmso.summary
scp my-ec2:~/rnaseq/ra/pon_dmso.tsv ./pon_dmso.tsv
scp my-ec2:~/rnaseq/ra/pon_dmso.tsv.summary ./pon_dmso.summary

fin_dmso.tsv
fin_dmso.tsv.summary
oza_dmso.tsv
oza_dmso.tsv.summary
pon_dmso.tsv
pon_dmso.tsv.summary
ra_counts.tsv
ra_counts.tsv.summary

