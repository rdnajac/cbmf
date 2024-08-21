import sys
import unittest
from pathlib import Path
from .core import align_reads, GenomeManager
from .utils import parse_args, run_script
from .utils.logger import logger, setup_logger
from .config import SCRIPTS_DIR


def run_tests():
    loader = unittest.TestLoader()
    start_dir = Path(__file__).parent / "tests"
    suite = loader.discover(start_dir, pattern="test_*.py")

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return 0 if result.wasSuccessful() else 1


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]  # Skip the script name

    args = parse_args(argv)
    setup_logger(level=args.log_level)
    logger.info(f"Received command: {args.command}")

    try:
        if args.command == "download":
            GenomeManager.download(args.species, args.files)
        # TODO
        # elif args.command == "align":
        #     align_reads(args)
        elif args.command == "test":
            return run_tests()
        # else:
        # script_path = SCRIPTS_DIR / f"{args.command}.sh"
        # return run_script(script_path)
    except ValueError as e:
        logger.error(f"Error: {str(e)}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"Script not found: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"An unexpected error occurred: {str(e)}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
