#!/usr/bin/env python3

import argparse
from a_series_of_tubes.ec2manager.client import EC2Helper


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        nargs="?",
        choices=["start", "stop", "ssh", "update"],
        help="Action to perform",
    )
    parser.add_argument("instance_id", nargs="?", help="Optional instance id")
    args = parser.parse_args()

    ec2 = EC2Helper()

    instance = (
        ec2.find_by_id(args.instance_id) if args.instance_id else ec2.pick_instance()
    )

    action = args.action if args.action else ec2.pick_action(instance)

    if action == "start":
        ec2.start(instance)
    elif action == "stop":
        ec2.stop(instance)
    elif action == "ssh":
        ec2.ssh(instance)
    elif action == "update":
        ec2.update_ssh_config(instance)


if __name__ == "__main__":
    main()
