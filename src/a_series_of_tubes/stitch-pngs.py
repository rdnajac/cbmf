#!/usr/bin/env python3
"""
Image Stitcher - A script to stitch multiple images together horizontally.
This is a Python implementation of the shell script that uses ImageMagick's montage command.
"""

import argparse
import subprocess
import sys
import shlex


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Stitch multiple images together horizontally using ImageMagick's montage command.",
        epilog="Example: python image_stitcher.py -o output.jpg input1.jpg input2.jpg input3.jpg",
    )

    # Add arguments
    parser.add_argument(
        "-o", "--output", required=True, help="Output filename (required)"
    )

    parser.add_argument(
        "input_files",
        nargs="+",
        help="Input image files to stitch together (at least one required)",
    )

    # Optional arguments that extend functionality
    parser.add_argument(
        "-t",
        "--tile",
        default="x1",
        help="Tiling specification (default: 'x1' for horizontal layout)",
    )

    parser.add_argument(
        "-g",
        "--geometry",
        default="+0+0",
        help="Geometry specification for montage (default: '+0+0' for no gaps)",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print verbose output"
    )

    # Parse arguments
    args = parser.parse_args()

    if args.verbose:
        print(f"Output file: {args.output}")
        print(f"Input files: {' '.join(args.input_files)}")
        print(f"Tile configuration: {args.tile}")
        print(f"Geometry: {args.geometry}")

    # Construct and run the montage command
    try:
        cmd = ["montage"]
        cmd.extend(args.input_files)
        cmd.extend(["-tile", args.tile])
        cmd.extend(["-geometry", args.geometry])
        cmd.append(args.output)

        if args.verbose:
            print(f"Running command: {shlex.join(cmd)}")

        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        if args.verbose and result.stdout:
            print(result.stdout)

        print(f"Successfully created {args.output}")
        return 0

    except subprocess.CalledProcessError as e:
        print(f"Error executing montage command: {e}", file=sys.stderr)
        if e.stderr:
            print(e.stderr, file=sys.stderr)
        return 1

    except FileNotFoundError:
        print(
            "Error: ImageMagick's 'montage' command not found. Please ensure ImageMagick is installed.",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
