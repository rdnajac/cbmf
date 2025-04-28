import subprocess


def pick(items, header=None):
    """Launch fzf with given items and return the selected one."""
    input_data = "\n".join(str(item) for item in items)

    cmd = ["fzf"]
    if header:
        cmd.extend(["--header", header])

    try:
        result = subprocess.run(
            cmd, input=input_data, capture_output=True, text=True, check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError:
        return None
