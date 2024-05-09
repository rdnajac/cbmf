#!/bin/bash

# aws cli list snapshots
aws ec2 describe-snapshots --owner-ids self --query 'Snapshots[*].{ID:SnapshotId,VolumeID:VolumeId,StartTime:StartTime,State:State,Progress:Progress,VolumeSize:VolumeSize,Description:Description}' --output table

# run that again but turn it into a tsv file:
aws ec2 describe-snapshots --owner-ids self --query 'Snapshots[*].{ID:SnapshotId,VolumeID:VolumeId,StartTime:StartTime,State:State,Progress:Progress,VolumeSize:VolumeSize,Description:Description}' --output text > snapshots.tsv


