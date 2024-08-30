import sys
import unittest
from pathlib import Path
from .core import align_reads, GenomeManager  # Adjust import paths as necessary
from .utils import parse_args, setup_logger
from .utils.logger import logger


def test() -> int:
    """Run the test suite."""
    loader = unittest.TestLoader()
    start_dir = Path(__file__).parent / "tests"
    suite = loader.discover(start_dir, pattern="test_*.py")
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return 0 if result.wasSuccessful() else 1


def download(args) -> None:
    """Handle the download command."""
    manager = GenomeManager()
    if args.files == ['ALL']:
        files = 'ALL'
    else:
        files = args.files

    manager.download(args.species, files)


def align(args) -> None:
    """Handle the align command."""
    align_reads(args)


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]  # Skip the script name

    args = parse_args(argv)
    setup_logger(level=args.log_level)
    logger.info(f"Received command: {args.command}")

    try:
        if args.command == "download":
            download(args)
        elif args.command == "align":
            align(args)
        elif args.command == "test":
            return test()
    except ValueError as e:
        logger.error(f"Error: {str(e)}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"File not found: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"An unexpected error occurred: {str(e)}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
