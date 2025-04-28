class EC2Instance:
    """Minimal class to represent an EC2 instance"""

    def __init__(self, instance_id, name, instance_type, state, dns, helper=None):
        self.id = instance_id
        self.name = name
        self.type = instance_type
        self.state = state
        self.dns = dns
        self.helper = helper  # EC2Helper assigned

    def __repr__(self):
        return f"[{self.id}, ({self.name}), type={self.type}, state={self.state}, dns={self.dns}]"

    def __str__(self):
        return (
            f"{self.name:<25} {self.id:<20} {self.type:<15} {self.state:<10} {self.dns}"
        )

    def refresh(self):
        """Reload this instance's state and DNS from AWS"""
        if not self.helper:
            raise RuntimeError("No EC2Helper attached to instance for refresh")

        refreshed = self.helper.ec2.describe_instances(InstanceIds=[self.id])
        data = refreshed["Reservations"][0]["Instances"][0]
        self.state = data["State"]["Name"]
        self.dns = data.get("PublicDnsName", "")
