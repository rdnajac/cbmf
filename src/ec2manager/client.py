import os
import boto3
from pathlib import Path

from .instance import EC2Instance
from mypyfzf import pick

class EC2Helper:
    """Encapsulation of AWS EC2 operations"""

    def __init__(self):
        self.ssh_config = os.path.expanduser("~/.ssh/config")
        self.ec2 = boto3.client("ec2")
        self.instances = self._load_instances()

    def _load_instances(self):
        instances = []
        response = self.ec2.describe_instances()
        for reservation in response["Reservations"]:
            for inst in reservation["Instances"]:
                instance_id = inst["InstanceId"]
                instance_type = inst["InstanceType"]
                state = inst["State"]["Name"]
                dns = inst.get("PublicDnsName", "")
                name = ""
                for tag in inst.get("Tags", []):
                    if tag["Key"] == "Name":
                        name = tag["Value"]
                        break
                instance = EC2Instance(
                    instance_id, name, instance_type, state, dns, helper=self
                )
                instances.append(instance)
        return instances

    def pick_instance(self):
        result = pick(self.instances, header="Select an instance")
        if not result:
            raise SystemExit("Nothing picked")
        return self.find_by_id(result.split()[1])

    def pick_action(self, instance=None):
        actions = ["start", "stop", "ssh", "update", "cancel"]
        header = (
            f"Select action for {instance.name} ({instance.id})" if instance else None
        )
        result = pick(actions, header=header)
        if result == "cancel" or not result:
            raise SystemExit(f"No action selected for instance:\n{instance}\n")
        return result

    def find_by_id(self, id):
        for instance in self.instances:
            if instance.id == id:
                return instance
        raise SystemExit(f"Instance not found: {id}")

    def update_ssh_config(self, instance):
        if instance.state != "running" or not instance.dns.startswith("ec2-"):
            return

        host = instance.name.replace(" ", "_")
        Path(self.ssh_config).touch(exist_ok=True)

        with open(self.ssh_config, "r") as f:
            lines = f.readlines()

        updated_lines = []
        in_host_block = False

        for line in lines:
            stripped = line.strip()
            if stripped == f"Host {host}":
                in_host_block = True
                updated_lines.append(line)
                continue
            if in_host_block:
                if stripped.startswith("Host "):
                    in_host_block = False
                elif stripped.startswith("HostName"):
                    line = f"  HostName {instance.dns}\n"
                    in_host_block = False
            updated_lines.append(line)

        with open(self.ssh_config, "w") as f:
            f.writelines(updated_lines)

        print(f"Updated Hostname in SSH config for {host} to {instance.dns}")

    def start(self, instance):
        print(f"Starting {instance.id}...")
        self.ec2.start_instances(InstanceIds=[instance.id])
        self.ec2.get_waiter("instance_running").wait(InstanceIds=[instance.id])
        print(f"Started: {instance.name}")

        instance.refresh()
        self.update_ssh_config(instance)

    def stop(self, instance):
        print(f"Stopping {instance.id}...")
        self.ec2.stop_instances(InstanceIds=[instance.id])
        self.ec2.get_waiter("instance_stopped").wait(InstanceIds=[instance.id])
        print(f"Stopped: {instance.name}")

    def ssh(self, instance):
        host = instance.name.replace(" ", "_")
        print(f"Connecting to: {host}\n")
        os.execvp("ssh", ["ssh", "-t", host, 'tmux new-session -A -s main'])
