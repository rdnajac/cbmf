# utils/colorprinter.py
import sys
from enum import Enum
from typing import Dict, Optional


class Color(Enum):
    BLACK = 0
    RED = 1
    GREEN = 2
    YELLOW = 3
    BLUE = 4
    MAGENTA = 5
    CYAN = 6
    WHITE = 7


class ColorCode:
    FOREGROUND_BASE = 30
    BACKGROUND_BASE = 40
    RESET = "\033[0m"

    @staticmethod
    def fg(color: Color) -> str:
        return f"\033[{ColorCode.FOREGROUND_BASE + color.value}m"

    @staticmethod
    def bg(color: Color) -> str:
        return f"\033[{ColorCode.BACKGROUND_BASE + color.value}m"


class ColorPrinter:
    """Class for printing colored messages to the console."""

    _color_map: Dict[str, Color] = {
        "success": Color.BLUE,
        "info": Color.GREEN,
        "warning": Color.YELLOW,
        "error": Color.RED,
    }

    @staticmethod
    def print_color(
        msg: str,
        fg: Optional[Color] = None,
        bg: Optional[Color] = None,
        file=sys.stderr,
    ):
        """Prints a message with given foreground and background colors."""
        fg_code = ColorCode.fg(fg) if fg else ""
        bg_code = ColorCode.bg(bg) if bg else ""
        print(f"{fg_code}{bg_code}{msg}{ColorCode.RESET}", file=file)

    @classmethod
    def print_status(cls, status: str, msg: str):
        """Prints a status message with predefined color."""
        color = cls._color_map.get(status.lower())
        if color:
            cls.print_color(f"{status.capitalize()}: {msg}", fg=color, file=sys.stderr)
        else:
            print(f"{status.capitalize()}: {msg}", file=sys.stderr)

    @classmethod
    def success(cls, msg: str):
        """Prints a success message in green."""
        cls.print_status("success", msg)

    @classmethod
    def info(cls, msg: str):
        """Prints an info message in blue."""
        cls.print_status("info", msg)

    @classmethod
    def warning(cls, msg: str):
        """Prints a warning message in yellow."""
        cls.print_status("warning", msg)

    @classmethod
    def error(cls, msg: str):
        """Prints an error message in red."""
        cls.print_status("error", msg)

