#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Class containing functions to print color msgs to the console.
import sys


class ANSIcolor:
    class Foreground:
        """ANSI escape codes for foreground colors."""

        BLACK = "\033[30m"
        RED = "\033[31m"
        GREEN = "\033[32m"
        YELLOW = "\033[33m"
        BLUE = "\033[34m"
        MAGENTA = "\033[35m"
        CYAN = "\033[36m"
        WHITE = "\033[37m"

    class Background:
        """ANSI escape codes for background colors."""

        BLACK = "\033[40m"
        RED = "\033[41m"
        GREEN = "\033[42m"
        YELLOW = "\033[43m"
        BLUE = "\033[44m"
        MAGENTA = "\033[45m"
        CYAN = "\033[46m"
        WHITE = "\033[47m"


# Aliases
fg = ANSIcolor.Foreground
bg = ANSIcolor.Background


class ColorPrinter:
    """Class containing functions to print color msgs to the console."""

    RESET = "\033[0m"

    @staticmethod
    def pr_msg(msg, fg="", bg=""):
        """Prints a msg with given foreground and background colors."""
        print(f"{fg}{bg}{msg}{ColorPrinter.RESET}", file=sys.stderr)

    @staticmethod
    def success(msg):
        """Prints a success msg in green."""
        ColorPrinter.pr_msg("Success: " + msg, fg.GREEN)

    @staticmethod
    def info(msg):
        """Prints an info msg in blue."""
        ColorPrinter.pr_msg("Info: " + msg, fg.BLUE)

    @staticmethod
    def warn(msg):
        """Prints a warning msg in yellow."""
        ColorPrinter.pr_msg("Warning: " + msg, fg.YELLOW)

    @staticmethod
    def error(msg):
        """Prints an error msg in red."""
        ColorPrinter.pr_msg("Error: " + msg, fg.RED)

