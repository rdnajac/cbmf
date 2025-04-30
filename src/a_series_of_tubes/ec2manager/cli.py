#!/usr/bin/env python3

import argparse
from a_series_of_tubes.ec2manager.ec2_helper import EC2Helper

dispatch = {
    "start": EC2Helper.start,
    "stop": EC2Helper.stop,
    "ssh": EC2Helper.ssh,
    "update": EC2Helper.update_ssh_config,
}


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Manage EC2 instances",
    )
    parser.add_argument(
        "action",
        nargs="?",
        choices=dispatch.keys(),
        help="Action to perform",
    )
    parser.add_argument("instance_id", nargs="?", help="Optional instance id")
    args = parser.parse_args(argv)

    if argv in (["-h"], ["--help"]):
        parser.print_help()
        return 0

    ec2 = EC2Helper()
    instance = (
        ec2.find_by_id(args.instance_id) if args.instance_id else ec2.pick_instance()
    )
    action = args.action if args.action else ec2.pick_action(instance)

    dispatch[action](ec2, instance)


if __name__ == "__main__":
    main()
