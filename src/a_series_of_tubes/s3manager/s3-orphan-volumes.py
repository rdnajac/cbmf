#!/usr/bin/env python3
import boto3

def get_attached_volume_ids(ec2):
    attached_volumes = set()
    reservations = ec2.describe_instances()['Reservations']
    for reservation in reservations:
        for instance in reservation['Instances']:
            for mapping in instance.get('BlockDeviceMappings', []):
                ebs = mapping.get('Ebs')
                if ebs and 'VolumeId' in ebs:
                    attached_volumes.add(ebs['VolumeId'])
    return attached_volumes

def get_ami_volume_ids(ec2):
    ami_volumes = set()
    images = ec2.describe_images(Owners=['self'])['Images']
    for image in images:
        for mapping in image.get('BlockDeviceMappings', []):
            ebs = mapping.get('Ebs')
            if ebs and 'SnapshotId' in ebs:
                # AMIs use snapshots, not volumes, but you might want to correlate snapshots to volumes if needed
                pass
    return ami_volumes  # Not directly used here, unless snapshot correlation is desired

def get_all_volumes(ec2):
    all_volumes = set()
    paginator = ec2.get_paginator('describe_volumes')
    for page in paginator.paginate():
        for volume in page['Volumes']:
            all_volumes.add(volume['VolumeId'])
    return all_volumes

def main():
    session = boto3.Session()
    ec2 = session.client('ec2')

    print("Fetching attached volume IDs from instances...")
    attached_volumes = get_attached_volume_ids(ec2)

    print("Fetching all volume IDs in the region...")
    all_volumes = get_all_volumes(ec2)

    orphaned_volumes = all_volumes - attached_volumes

    if orphaned_volumes:
        print("Orphaned volumes:")
        for vol in sorted(orphaned_volumes):
            print(f"  {vol}")
    else:
        print("No orphaned volumes found.")

if __name__ == "__main__":
    main()
