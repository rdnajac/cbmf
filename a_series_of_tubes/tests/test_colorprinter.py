# tests/test_colorprinter.py
from ..utils import ColorPrinter as pr


def smoke_test():
    for status in ["success", "info", "warning", "error"]:
        getattr(pr, status)(f"This is a {status} message")
