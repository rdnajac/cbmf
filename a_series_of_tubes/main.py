import sys
from pathlib import Path
from .core import GenomeManager
from .config import SCRIPTS_DIR
from .utils import parse_args, ColorPrinter as pr, run_script

def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]  # Skip the script name

    args = parse_args(argv)
    pr.info(f"Received command: {args.command}")

    try:
        if args.command == "download":
            GenomeManager.download(args.species, args.files)
        elif args.command == "test":
            # Implement test logic here
            pr.info("Running tests...")
            # You might want to implement a test runner or use a testing framework
        else:
            script_path = SCRIPTS_DIR / f"{args.command}.sh"
            return run_script(script_path)
    except ValueError as e:
        pr.error(f"Error: {str(e)}")
        return 1
    except FileNotFoundError as e:
        pr.error(f"Script not found: {str(e)}")
        return 1
    except Exception as e:
        pr.error(f"An unexpected error occurred: {str(e)}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
